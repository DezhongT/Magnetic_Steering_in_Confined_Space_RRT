#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
setInput makeInput(const Actuation &actuation)
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetStringOpt("magneticModel") = "axial_tip";
	input.GetStringOpt("insertionModel") = "proximal_guide";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 200;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetScalarOpt("tipDipoleMoment") = 1.0e-3;
	input.GetScalarOpt("insertionCoordinate") = actuation.xi;
	input.GetScalarOpt("insertionStiffness") = 1.0e3;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.22);
	input.GetVecOpt("baVector") = actuation.field;
	return input;
}

double positionDifference(const RodState &left, const RodState &right)
{
	double squared = 0.0;
	for (int dof = 0; dof < left.configuration.size(); ++dof)
	{
		if (dof % 4 != 3)
		{
			squared += std::pow(
				left.configuration[dof] - right.configuration[dof], 2);
		}
	}
	return std::sqrt(squared);
}

bool sameContactIds(
	const std::vector<ContactCandidate> &left,
	const std::vector<ContactCandidate> &right)
{
	if (left.size() != right.size())
	{
		return false;
	}
	for (int i = 0; i < static_cast<int>(left.size()); ++i)
	{
		if (left[i].rodVertex != right[i].rodVertex ||
			left[i].boundaryId != right[i].boundaryId)
		{
			return false;
		}
	}
	return true;
}
}

int main()
{
	Actuation start;
	start.xi = 1.0e-3;
	start.field = Vector3d(0.0, 0.1, 0.0);
	setInput input = makeInput(start);
	world simulation(input);
	simulation.setRodStepper();
	const ContactKktEquilibriumResult base =
		simulation.solvePlanarContactKktEquilibrium();
	if (!base.success)
	{
		std::cerr << "Insertion-actuated base equilibrium failed.\n";
		return 1;
	}
	const ContactStabilityAnalysis stability =
		simulation.analyzePlanarContactStability(base);
	const ActuationEquilibriumSensitivity analytic =
		simulation.computePlanarContactActuationSensitivity(base);
	if (!stability.valid || !analytic.success)
	{
		std::cerr << "Insertion stability or sensitivity failed.\n";
		return 1;
	}

	constexpr double xiStep = 1.0e-6;
	Actuation plusActuation = start;
	Actuation minusActuation = start;
	plusActuation.xi += xiStep;
	minusActuation.xi -= xiStep;
	simulation.setActuation(plusActuation);
	const ContactKktEquilibriumResult plus =
		simulation.correctPlanarContactKktEquilibrium(base);
	simulation.setActuation(minusActuation);
	const ContactKktEquilibriumResult minus =
		simulation.correctPlanarContactKktEquilibrium(base);
	if (!plus.success || !minus.success ||
		!sameContactIds(base.activeContacts, plus.activeContacts) ||
		!sameContactIds(base.activeContacts, minus.activeContacts))
	{
		std::cerr << "Insertion finite-difference equilibria failed.\n";
		return 1;
	}
	const VectorXd configurationDifference =
		(plus.state.configuration - minus.state.configuration) / (2.0 * xiStep);
	const VectorXd multiplierDifference =
		(plus.multipliers - minus.multipliers) / (2.0 * xiStep);
	const double configurationError =
		(configurationDifference - analytic.configurationDerivative.col(0)).norm() /
		std::max(1.0, configurationDifference.norm());
	const double multiplierError =
		(multiplierDifference - analytic.multiplierDerivative.col(0)).norm() /
		std::max(1.0, multiplierDifference.norm());
	if (configurationError > 1.0e-4 || multiplierError > 1.0e-4)
	{
		std::cerr << "Insertion equilibrium sensitivity failed: q_error="
				  << configurationError << ", lambda_error="
				  << multiplierError << ".\n";
		return 1;
	}

	simulation.setActuation(start);
	simulation.restoreRodState(base.state);
	const PlannerState plannerState =
		simulation.capturePlannerState(base, stability);
	Actuation changed = start;
	changed.xi += 5.0e-4;
	changed.field += Vector3d(0.01, -0.02, 0.03);
	simulation.setActuation(changed);
	VectorXd perturbation = VectorXd::Constant(simulation.numFreeDofs(), 1.0e-7);
	simulation.applyFreeDofIncrement(perturbation);
	simulation.restorePlannerState(plannerState);
	if ((simulation.captureRodState().configuration -
		 plannerState.rodState.configuration).norm() != 0.0 ||
		std::abs(simulation.getActuation().xi - start.xi) > 0.0 ||
		(simulation.getActuation().field - start.field).norm() != 0.0)
	{
		std::cerr << "PlannerState actuation/rod round trip failed.\n";
		return 1;
	}

	Actuation target = start;
	target.xi += 2.0e-4;
	target.field += Vector3d(0.01, 0.01, 0.005);
	FieldContinuationOptions options;
	options.initialStepFraction = 0.25;
	options.maximumStepFraction = 0.4;
	options.easyCorrectorIterations = 10;
	const ActuationContinuationResult continuation =
		simulation.continuePlanarContactActuation(plannerState, target, options);
	if (!continuation.success || continuation.rolledBack ||
		continuation.points.size() < 3)
	{
		std::cerr << "Combined actuation continuation failed: success="
				  << continuation.success << ", rolled_back="
				  << continuation.rolledBack << ", points="
				  << continuation.points.size() << ", attempts="
				  << continuation.attemptedSteps << ", rejected="
				  << continuation.rejectedSteps << ", minimum_step="
				  << continuation.minimumAttemptedStepFraction << ".\n";
		return 1;
	}
	for (int point = 1; point < static_cast<int>(continuation.points.size()); ++point)
	{
		if (continuation.points[point].actuation.xi <
			continuation.points[point - 1].actuation.xi)
		{
			std::cerr << "Insertion continuation retracted unexpectedly.\n";
			return 1;
		}
	}

	setInput directInput = makeInput(target);
	world directWorld(directInput);
	directWorld.setRodStepper();
	const ContactKktEquilibriumResult direct =
		directWorld.solvePlanarContactKktEquilibrium();
	const ContactKktEquilibriumResult &continued =
		continuation.points.back().equilibrium;
	if (!direct.success ||
		!sameContactIds(direct.activeContacts, continued.activeContacts) ||
		positionDifference(direct.state, continued.state) > 1.0e-6 ||
		(direct.multipliers - continued.multipliers).norm() > 1.0e-6)
	{
		std::cerr << "Combined continuation endpoint failed direct comparison.\n";
		return 1;
	}

	std::cout << "Actuation sensitivity/continuation: q_error="
			  << configurationError << ", lambda_error=" << multiplierError
			  << ", points=" << continuation.points.size()
			  << ", endpoint_position_error="
			  << positionDifference(direct.state, continued.state) << '\n';
	return 0;
}
