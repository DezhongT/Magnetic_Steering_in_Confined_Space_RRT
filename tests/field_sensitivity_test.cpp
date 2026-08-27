#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <option-file>\n";
		return 1;
	}

	setInput input;
	if (input.LoadOptions(argv[1]) != 0)
	{
		return 1;
	}
	world simulation(input);
	simulation.setRodStepper();
	const EquilibriumResult baseResult = simulation.solveStaticEquilibrium();
	if (!baseResult.success)
	{
		std::cerr << "Base magnetic equilibrium failed.\n";
		return 1;
	}

	const RodState baseState = simulation.captureRodState();
	const Vector3d baseField = simulation.getAppliedField();
	const MatrixXd analyticSensitivity =
		simulation.computeConfigurationFieldSensitivity();
	constexpr double step = 1.0e-5;
	MatrixXd finiteDifference = MatrixXd::Zero(
		baseState.configuration.size(), 3);

	for (int component = 0; component < 3; ++component)
	{
		Vector3d plusField = baseField;
		Vector3d minusField = baseField;
		plusField[component] += step;
		minusField[component] -= step;

		simulation.restoreRodState(baseState);
		simulation.setAppliedField(plusField);
		const EquilibriumResult plusResult = simulation.solveStaticEquilibrium();
		simulation.restoreRodState(baseState);
		simulation.setAppliedField(minusField);
		const EquilibriumResult minusResult = simulation.solveStaticEquilibrium();
		if (!plusResult.success || !minusResult.success)
		{
			std::cerr << "Perturbed magnetic equilibrium failed for component "
					  << component << ".\n";
			return 1;
		}
		finiteDifference.col(component) =
			(plusResult.state.configuration - minusResult.state.configuration) /
			(2.0 * step);
	}

	simulation.restoreRodState(baseState);
	simulation.setAppliedField(baseField);
	const double scale = std::max(
		1.0, std::max(analyticSensitivity.norm(), finiteDifference.norm()));
	const double relativeError =
		(analyticSensitivity - finiteDifference).norm() / scale;
	const int tipIndex = 4 * (simulation.numPoints() - 1);
	const Matrix3d tipSensitivity =
		analyticSensitivity.block<3, 3>(tipIndex, 0);

	std::cout << "Field sensitivity relative error=" << relativeError
			  << ", tip sensitivity=\n" << tipSensitivity << '\n';
	if (!std::isfinite(relativeError) || relativeError > 1.0e-4)
	{
		return 1;
	}
	return 0;
}
