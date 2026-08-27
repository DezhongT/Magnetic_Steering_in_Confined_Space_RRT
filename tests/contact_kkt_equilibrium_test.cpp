#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

namespace
{
setInput makeInput(double gravity, int maxIterations = 200)
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = maxIterations;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, gravity);
	return input;
}

bool checkRetainedContact()
{
	setInput input = makeInput(0.22);
	world simulation(input);
	simulation.setRodStepper();
	const ContactKktEquilibriumResult result =
		simulation.solvePlanarContactKktEquilibrium();
	if (!result.success || result.rolledBack || !result.tipSafe ||
		result.activeContacts.empty() ||
		(result.multipliers.array() <= 0.0).any() ||
		result.stationarityNorm > 1.0e-6 ||
		result.activeConstraintNorm > 1.0e-8 ||
		result.complementarityNorm > 1.0e-10 ||
		result.minimumBodyGap < -1.0e-8)
	{
		std::cerr << "Retained-contact nonlinear KKT solve failed: success="
				  << result.success << ", rollback=" << result.rolledBack
				  << ", contacts=" << result.activeContacts.size()
				  << ", stationarity=" << result.stationarityNorm
				  << ", constraint=" << result.activeConstraintNorm
				  << ", complementarity=" << result.complementarityNorm
				  << ", minimum_gap=" << result.minimumBodyGap
				  << ", minimum_multiplier=" << result.minimumMultiplier
				  << ", tip_safe=" << result.tipSafe
				  << ", iterations=" << result.nonlinearIterations
				  << ", backtracks=" << result.lineSearchBacktracks << ".\n";
		return false;
	}
	std::cout << "Retained-contact KKT: iterations="
			  << result.nonlinearIterations
			  << ", active_updates=" << result.activeSetUpdates
			  << ", backtracks=" << result.lineSearchBacktracks
			  << ", stationarity=" << result.stationarityNorm
			  << ", constraint=" << result.activeConstraintNorm
			  << ", multiplier=" << result.minimumMultiplier << '\n';
	return true;
}

bool checkReleasedContact()
{
	setInput input = makeInput(0.2);
	world simulation(input);
	simulation.setRodStepper();
	const ContactKktEquilibriumResult result =
		simulation.solvePlanarContactKktEquilibrium();
	if (!result.success || result.rolledBack || !result.tipSafe ||
		!result.activeContacts.empty() || result.activeSetUpdates <= 0 ||
		result.stationarityNorm > 1.0e-6 || result.minimumBodyGap <= 0.0)
	{
		std::cerr << "Negative-multiplier contact release failed: success="
				  << result.success << ", rollback=" << result.rolledBack
				  << ", contacts=" << result.activeContacts.size()
				  << ", active_updates=" << result.activeSetUpdates
				  << ", stationarity=" << result.stationarityNorm
				  << ", minimum_gap=" << result.minimumBodyGap << ".\n";
		return false;
	}
	std::cout << "Released-contact KKT: iterations="
			  << result.nonlinearIterations
			  << ", active_updates=" << result.activeSetUpdates
			  << ", minimum_gap=" << result.minimumBodyGap << '\n';
	return true;
}

bool checkFailureRollback()
{
	setInput input = makeInput(0.22, 200);
	input.GetScalarOpt("tipSafeDistance") = 0.02;
	world simulation(input);
	simulation.setRodStepper();
	const RodState start = simulation.captureRodState();
	const ContactKktEquilibriumResult result =
		simulation.solvePlanarContactKktEquilibrium();
	const RodState end = simulation.captureRodState();
	if (result.success || !result.rolledBack || result.tipSafe ||
		(end.configuration - start.configuration).norm() != 0.0 ||
		(end.previousConfiguration - start.previousConfiguration).norm() != 0.0)
	{
		std::cerr << "Unsafe-tip KKT rejection did not roll back exactly: tip_safe="
				  << result.tipSafe << ", iterations="
				  << result.nonlinearIterations << ".\n";
		return false;
	}
	return true;
}
}

int main()
{
	if (!checkRetainedContact() || !checkReleasedContact() ||
		!checkFailureRollback())
	{
		return 1;
	}
	return 0;
}
