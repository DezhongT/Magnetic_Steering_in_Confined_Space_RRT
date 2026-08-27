#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <option-file>\n";
		return 1;
	}

	setInput staticInput;
	setInput dynamicInput;
	if (staticInput.LoadOptions(argv[1]) != 0 ||
		dynamicInput.LoadOptions(argv[1]) != 0)
	{
		return 1;
	}

	world staticSimulation(staticInput);
	staticSimulation.setRodStepper();
	const EquilibriumResult staticResult = staticSimulation.solveStaticEquilibrium();
	if (!staticResult.success)
	{
		std::cerr << "Gravity-only static solve failed with residual "
				  << staticResult.finalResidualNorm
				  << ", iterations=" << staticResult.newtonIterations
				  << ", backtracks=" << staticResult.lineSearchBacktracks
				  << ", minimum_step=" << staticResult.minimumAcceptedStepLength
				  << ", line_search_failed=" << staticResult.lineSearchFailed
				  << ".\n";
		return 1;
	}

	world relaxedSimulation(dynamicInput);
	relaxedSimulation.setRodStepper();
	constexpr double residualTolerance = 1.0e-6;
	constexpr double velocityTolerance = 1.0e-5;
	double relaxedResidual = relaxedSimulation.getStaticResidualNorm();
	bool relaxed = false;
	int steps = 0;

	while (relaxedSimulation.simulationRunning() > 0)
	{
		relaxedSimulation.updateTimeStep();
		++steps;
		relaxedResidual = relaxedSimulation.getStaticResidualNorm();
		if (relaxedResidual <= residualTolerance &&
			relaxedSimulation.getVelocityNorm() <= velocityTolerance)
		{
			relaxed = true;
			break;
		}
	}

	if (!relaxed)
	{
		std::cerr << "Dynamic relaxation did not reach the requested tolerances. "
				  << "residual=" << relaxedResidual
				  << ", velocity=" << relaxedSimulation.getVelocityNorm()
				  << ", steps=" << steps << '\n';
		return 1;
	}

	const RodState relaxedState = relaxedSimulation.captureRodState();
	const double configurationDifference =
		(relaxedState.configuration - staticResult.state.configuration).norm();
	constexpr double configurationTolerance = 1.0e-4;
	if (configurationDifference > configurationTolerance)
	{
		std::cerr << "Static and dynamically relaxed equilibria differ by "
				  << configurationDifference << ".\n";
		return 1;
	}

	std::cout << "Relaxation comparison: steps=" << steps
			  << ", time=" << relaxedSimulation.getCurrentTime()
			  << ", static_residual=" << staticResult.finalResidualNorm
			  << ", relaxed_residual=" << relaxedResidual
			  << ", velocity=" << relaxedSimulation.getVelocityNorm()
			  << ", configuration_difference=" << configurationDifference
			  << '\n';
	return 0;
}
