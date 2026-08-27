#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

int main()
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 200;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetScalarOpt("tipSafeDistance") = 0.0;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.2);

	world simulation(input);
	simulation.setRodStepper();
	const EquilibriumResult result = simulation.solveStaticEquilibrium();
	if (!result.success)
	{
		std::cerr << "Loaded planar contact solve failed: residual="
				  << result.finalResidualNorm
				  << ", iterations=" << result.newtonIterations
				  << ", backtracks=" << result.lineSearchBacktracks
				  << ", infeasible_trials=" << result.infeasibleTrialRejections
				  << ", minimum_step=" << result.minimumAcceptedStepLength
				  << ", line_search_failed=" << result.lineSearchFailed << ".\n";
		return 1;
	}

	// Refresh and inspect the candidate/tip diagnostics at the accepted state.
	const double residual = simulation.getStaticResidualNorm();
	const ContactDetectionResult &contact = simulation.getLastContactDetection();
	const bool barrierActive = !contact.bodyCandidates.empty() &&
		contact.minimumBodyGap < input.GetScalarOpt("dBar");
	if (!std::isfinite(residual) || residual > 1.0e-6 ||
		!barrierActive || !(contact.minimumBodyGap > 0.0) || !contact.tip.safe)
	{
		std::cerr << "Loaded planar equilibrium did not satisfy acceptance checks: "
				  << "residual=" << residual
				  << ", candidates=" << contact.bodyCandidates.size()
				  << ", minimum_gap=" << contact.minimumBodyGap
				  << ", tip_clearance=" << contact.tip.clearance
				  << ", tip_safe=" << contact.tip.safe << ".\n";
		return 1;
	}
	if (result.lineSearchBacktracks <= 0 ||
		!(result.minimumAcceptedStepLength < 1.0))
	{
		std::cerr << "Loaded contact solve did not exercise backtracking.\n";
		return 1;
	}

	// A deliberately over-large load makes the first full Newton trial cross
	// the slab boundary. With only one Newton iteration available, the solve
	// must reject that trial, remain feasible during backtracking, then roll
	// back exactly when the iteration budget is exhausted.
	setInput rejectionInput;
	rejectionInput.GetStringOpt("contactModel") = "planar_barrier";
	rejectionInput.GetIntOpt("numVertices") = 20;
	rejectionInput.GetIntOpt("maxIter") = 1;
	rejectionInput.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 2.0);
	world rejectionSimulation(rejectionInput);
	rejectionSimulation.setRodStepper();
	const RodState rejectionStart = rejectionSimulation.captureRodState();
	const EquilibriumResult rejectionResult =
		rejectionSimulation.solveStaticEquilibrium();
	const RodState rejectionEnd = rejectionSimulation.captureRodState();
	if (rejectionResult.success ||
		rejectionResult.infeasibleTrialRejections <= 0 ||
		(rejectionEnd.configuration - rejectionStart.configuration).norm() != 0.0)
	{
		std::cerr << "Infeasible contact-trial rejection or rollback failed: "
				  << "success=" << rejectionResult.success
				  << ", infeasible_trials="
				  << rejectionResult.infeasibleTrialRejections
				  << ", rollback_error="
				  << (rejectionEnd.configuration -
					  rejectionStart.configuration).norm() << ".\n";
		return 1;
	}

	std::cout << "Loaded planar equilibrium: residual=" << residual
			  << ", iterations=" << result.newtonIterations
			  << ", backtracks=" << result.lineSearchBacktracks
			  << ", infeasible_trials=" << result.infeasibleTrialRejections
			  << ", minimum_step=" << result.minimumAcceptedStepLength
			  << ", minimum_gap=" << contact.minimumBodyGap
			  << ", tip_clearance=" << contact.tip.clearance << '\n';
	return 0;
}
