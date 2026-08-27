#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

namespace
{
constexpr double stateTolerance = 1.0e-8;

bool checkSuccessfulSolve()
{
	setInput input;
	world simulation(input);
	simulation.setRodStepper();

	const RodState equilibrium = simulation.captureRodState();
	const double initialTime = simulation.getCurrentTime();
	VectorXd increment = VectorXd::Zero(simulation.numFreeDofs());
	increment[increment.size() - 1] = 1.0e-3;
	simulation.applyFreeDofIncrement(increment);

	const EquilibriumResult result = simulation.solveStaticEquilibrium();
	if (!result.success)
	{
		std::cerr << "Static Newton failed with LAPACK info="
				  << result.linearSolverInfo << ", residual="
				  << result.finalResidualNorm << '\n';
		return false;
	}

	if (result.newtonIterations <= 0 || result.initialResidualNorm <= 1.0e-8)
	{
		std::cerr << "Static Newton did not exercise a perturbed configuration.\n";
		return false;
	}

	if (!std::isfinite(result.finalResidualNorm) || result.finalResidualNorm > 1.0e-7)
	{
		std::cerr << "Static Newton residual is too large: "
				  << result.finalResidualNorm << '\n';
		return false;
	}

	if (std::abs(simulation.getCurrentTime() - initialTime) > 0.0)
	{
		std::cerr << "Static Newton advanced physical time.\n";
		return false;
	}

	if ((result.state.configuration.head(7) - equilibrium.configuration.head(7)).norm() != 0.0)
	{
		std::cerr << "Static Newton moved a clamped degree of freedom.\n";
		return false;
	}

	if ((result.state.configuration - equilibrium.configuration).norm() > stateTolerance)
	{
		std::cerr << "Unloaded rod did not return to its straight equilibrium.\n";
		return false;
	}

	if (result.state.velocity.norm() != 0.0)
	{
		std::cerr << "Committed static equilibrium has nonzero velocity.\n";
		return false;
	}

	std::cout << "Static Newton: initial residual=" << result.initialResidualNorm
			  << ", final residual=" << result.finalResidualNorm
			  << ", iterations=" << result.newtonIterations << '\n';

	return true;
}

bool checkFailureRollback(const char *failureOptionFile)
{
	setInput input;
	if (input.LoadOptions(failureOptionFile) != 0)
	{
		return false;
	}

	world simulation(input);
	simulation.setRodStepper();
	VectorXd increment = VectorXd::Zero(simulation.numFreeDofs());
	increment[increment.size() - 1] = 1.0e-3;
	simulation.applyFreeDofIncrement(increment);
	const RodState inputState = simulation.captureRodState();
	const double initialTime = simulation.getCurrentTime();

	const EquilibriumResult result = simulation.solveStaticEquilibrium();
	if (result.success)
	{
		std::cerr << "A zero-iteration static solve unexpectedly succeeded.\n";
		return false;
	}

	const RodState restored = simulation.captureRodState();
	if ((restored.configuration - inputState.configuration).norm() != 0.0 ||
		(restored.previousConfiguration - inputState.previousConfiguration).norm() != 0.0 ||
		(restored.velocity - inputState.velocity).norm() != 0.0)
	{
		std::cerr << "Failed static solve did not restore its input state.\n";
		return false;
	}

	if (simulation.getCurrentTime() != initialTime)
	{
		std::cerr << "Failed static solve advanced physical time.\n";
		return false;
	}

	return true;
}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <failure-option-file>\n";
		return 1;
	}

	if (!checkSuccessfulSolve() || !checkFailureRollback(argv[1]))
	{
		return 1;
	}

	std::cout << "Static Newton convergence and rollback checks passed.\n";
	return 0;
}
