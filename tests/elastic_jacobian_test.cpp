#include "setInput.h"
#include "world.h"

#include <array>
#include <cmath>
#include <iostream>
#include <random>

namespace
{
VectorXd evaluatePerturbation(
	world &simulation,
	const RodState &baseState,
	const VectorXd &direction,
	double step)
{
	simulation.restoreRodState(baseState);
	simulation.applyFreeDofIncrement(step * direction);
	return simulation.evaluateStaticSystem().residual;
}
}

int main()
{
	setInput input;
	world simulation(input);
	simulation.setRodStepper();

	// Move away from the exactly straight configuration so bending, twisting,
	// and frame-transport derivatives are exercised at a generic state.
	VectorXd baseIncrement = VectorXd::Zero(simulation.numFreeDofs());
	for (int i = 0; i < baseIncrement.size(); ++i)
	{
		baseIncrement[i] = 1.0e-5 * std::sin(0.173 * static_cast<double>(i + 1));
	}
	simulation.applyFreeDofIncrement(baseIncrement);
	simulation.evaluateStaticSystem();
	const RodState baseState = simulation.captureRodState();
	const StaticEvaluation analyticEvaluation = simulation.evaluateStaticSystem();

	std::mt19937 generator(20260825u);
	std::normal_distribution<double> normal(0.0, 1.0);
	constexpr std::array<double, 4> steps = {1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
	constexpr int directionCount = 3;
	double worstBestRelativeError = 0.0;

	for (int sample = 0; sample < directionCount; ++sample)
	{
		VectorXd direction(simulation.numFreeDofs());
		for (int i = 0; i < direction.size(); ++i)
		{
			direction[i] = normal(generator);
		}
		direction.normalize();

		const VectorXd analyticProduct = analyticEvaluation.multiplyJacobian(direction);
		double bestRelativeError = std::numeric_limits<double>::infinity();
		double previousError = std::numeric_limits<double>::infinity();
		bool observedDecrease = false;

		for (const double step : steps)
		{
			const VectorXd residualPlus =
				evaluatePerturbation(simulation, baseState, direction, step);
			const VectorXd residualMinus =
				evaluatePerturbation(simulation, baseState, direction, -step);
			const VectorXd finiteDifference =
				(residualPlus - residualMinus) / (2.0 * step);

			const double scale = std::max(
				1.0, std::max(analyticProduct.norm(), finiteDifference.norm()));
			const double relativeError =
				(finiteDifference - analyticProduct).norm() / scale;
			bestRelativeError = std::min(bestRelativeError, relativeError);
			observedDecrease = observedDecrease || relativeError < previousError;
			previousError = relativeError;

			std::cout << "direction=" << sample << ", epsilon=" << step
					  << ", relative_error=" << relativeError << '\n';
		}

		if (!observedDecrease || !std::isfinite(bestRelativeError))
		{
			std::cerr << "Finite-difference error did not decrease for direction "
					  << sample << ".\n";
			return 1;
		}
		worstBestRelativeError = std::max(worstBestRelativeError, bestRelativeError);
	}

	simulation.restoreRodState(baseState);
	if (worstBestRelativeError > 1.0e-4)
	{
		std::cerr << "Elastic Jacobian directional error is too large: "
				  << worstBestRelativeError << '\n';
		return 1;
	}

	std::cout << "Worst best elastic Jacobian relative error: "
			  << worstBestRelativeError << '\n';
	return 0;
}
