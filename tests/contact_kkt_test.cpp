#include "contact/contactKktSystem.h"
#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool checkManufacturedSystem()
{
	ContactKktLinearization linearization;
	linearization.hessian.resize(2, 2);
	linearization.hessian << 4.0, 1.0, 1.0, 3.0;
	linearization.residual.resize(2);
	linearization.residual << 1.0, -2.0;
	linearization.constraintJacobian.resize(1, 2);
	linearization.constraintJacobian << 1.0, 0.0;
	linearization.gaps = VectorXd::Constant(1, 0.1);
	linearization.multipliers = VectorXd::Constant(1, 0.5);

	const ContactKktSystem system = buildContactKktSystem(linearization);
	if ((system.matrix - system.matrix.transpose()).norm() > 1.0e-14 ||
		system.matrix.rows() != 3 || system.primalDofs != 2 ||
		system.contactCount != 1)
	{
		std::cerr << "Manufactured KKT block assembly failed.\n";
		return false;
	}

	const ContactKktSolution solution = solveContactKktSystem(linearization);
	const VectorXd primalResidual =
		linearization.residual -
		linearization.constraintJacobian.transpose() * linearization.multipliers +
		linearization.hessian * solution.configurationIncrement -
		linearization.constraintJacobian.transpose() * solution.multiplierIncrement;
	const VectorXd constraintResidual =
		linearization.gaps +
		linearization.constraintJacobian * solution.configurationIncrement;
	if (!solution.success || solution.linearResidualNorm > 1.0e-12 ||
		primalResidual.norm() > 1.0e-12 || constraintResidual.norm() > 1.0e-12 ||
		(solution.updatedMultipliers -
		 (linearization.multipliers + solution.multiplierIncrement)).norm() != 0.0)
	{
		std::cerr << "Manufactured KKT solve or multiplier sign failed.\n";
		return false;
	}

	bool rejectedDimensions = false;
	try
	{
		ContactKktLinearization invalid = linearization;
		invalid.gaps.resize(2);
		(void)buildContactKktSystem(invalid);
	}
	catch (const std::invalid_argument &)
	{
		rejectedDimensions = true;
	}
	if (!rejectedDimensions)
	{
		std::cerr << "KKT builder accepted inconsistent dimensions.\n";
		return false;
	}
	return true;
}

bool checkActiveSetUpdate()
{
	ContactCandidate retained;
	retained.rodVertex = 3;
	retained.boundaryId = 1;
	retained.gap = 0.0;
	ContactCandidate released = retained;
	released.rodVertex = 4;
	ContactCandidate added = retained;
	added.rodVertex = 5;
	added.gap = -2.0e-5;
	ContactCandidate separated = retained;
	separated.rodVertex = 6;
	separated.gap = 2.0e-5;
	VectorXd multipliers(2);
	multipliers << 0.3, -0.2;

	const ContactActiveSetUpdate update = updateContactActiveSet(
		{retained, released}, multipliers, {added, separated}, 1.0e-6, 1.0e-6);
	if (!update.changed() || update.added != 1 || update.released != 1 ||
		update.contacts.size() != 2 || update.contacts[0].rodVertex != 3 ||
		update.contacts[1].rodVertex != 5 || update.multipliers[0] != 0.3 ||
		update.multipliers[1] != 0.0)
	{
		std::cerr << "Contact insertion/release rule failed.\n";
		return false;
	}
	return true;
}

bool checkBarrierSeededSystem()
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 200;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.22);

	world simulation(input);
	simulation.setRodStepper();
	const EquilibriumResult barrierResult = simulation.solveStaticEquilibrium();
	if (!barrierResult.success)
	{
		std::cerr << "Could not construct the barrier-equilibrium KKT seed.\n";
		return false;
	}

	const PlanarContactKktSeed seed = simulation.buildPlanarContactKktSeed();
	const ContactDetectionResult &detection = simulation.getLastContactDetection();
	const ContactKktLinearization &linearization = seed.linearization;
	const int contacts = static_cast<int>(seed.contacts.size());
	if (contacts <= 0 || linearization.constraintJacobian.rows() != contacts ||
		(linearization.multipliers.array() <= 0.0).any() || !detection.tip.safe)
	{
		std::cerr << "Barrier seed has no active contact or nonpositive multipliers.\n";
		return false;
	}

	const VectorXd seededStationarity =
		linearization.residual -
		linearization.constraintJacobian.transpose() * linearization.multipliers;
	const double stationarityScale = std::max(1.0, linearization.residual.norm());
	if (seededStationarity.norm() / stationarityScale > 1.0e-8)
	{
		std::cerr << "Barrier force did not map to the KKT multiplier convention: "
				  << seededStationarity.norm() / stationarityScale << ".\n";
		return false;
	}

	const ContactKktSystem system = buildContactKktSystem(linearization);
	const double symmetryScale = std::max(1.0, system.matrix.norm());
	if ((system.matrix - system.matrix.transpose()).norm() / symmetryScale > 1.0e-10)
	{
		std::cerr << "Seeded planar KKT matrix is not symmetric.\n";
		return false;
	}

	const ContactKktSolution solution = solveContactKktSystem(linearization);
	const VectorXd predictedGaps = linearization.gaps +
		linearization.constraintJacobian * solution.configurationIncrement;
	const VectorXd predictedStationarity = seededStationarity +
		linearization.hessian * solution.configurationIncrement -
		linearization.constraintJacobian.transpose() * solution.multiplierIncrement;
	if (!solution.success || predictedGaps.norm() > 1.0e-9 ||
		predictedStationarity.norm() > 1.0e-8 ||
		(solution.updatedMultipliers.array() <= 0.0).any())
	{
		std::cerr << "Barrier-seeded KKT correction failed: linear_residual="
				  << solution.linearResidualNorm
				  << ", gap_residual=" << predictedGaps.norm()
				  << ", stationarity=" << predictedStationarity.norm()
				  << ", minimum_multiplier="
				  << solution.updatedMultipliers.minCoeff() << ".\n";
		return false;
	}

	std::cout << "Barrier-seeded KKT: contacts=" << contacts
			  << ", seed_stationarity=" << seededStationarity.norm()
			  << ", linear_residual=" << solution.linearResidualNorm
			  << ", predicted_gap=" << predictedGaps.norm()
			  << ", minimum_multiplier="
			  << solution.updatedMultipliers.minCoeff() << '\n';
	return true;
}
}

int main()
{
	if (!checkManufacturedSystem() || !checkActiveSetUpdate() ||
		!checkBarrierSeededSystem())
	{
		return 1;
	}
	return 0;
}
