#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

namespace
{
setInput makeInput(const Vector3d &field)
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetStringOpt("magneticModel") = "axial_tip";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 200;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetScalarOpt("tipDipoleMoment") = 1.0e-3;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.22);
	input.GetVecOpt("baVector") = field;
	return input;
}

bool sameContactIds(
	const std::vector<ContactCandidate> &left,
	const std::vector<ContactCandidate> &right)
{
	if (left.size() != right.size())
	{
		return false;
	}
	for (int contact = 0; contact < static_cast<int>(left.size()); ++contact)
	{
		if (left[contact].rodVertex != right[contact].rodVertex ||
			left[contact].boundaryId != right[contact].boundaryId)
		{
			return false;
		}
	}
	return true;
}

double positionDifference(const RodState &left, const RodState &right)
{
	double squaredDifference = 0.0;
	for (int dof = 0; dof < left.configuration.size(); ++dof)
	{
		if (dof % 4 != 3)
		{
			squaredDifference += std::pow(
				left.configuration[dof] - right.configuration[dof], 2);
		}
	}
	return std::sqrt(squaredDifference);
}

bool checkFixedActiveContinuation()
{
	const Vector3d startField(0.0, 0.1, 0.0);
	const Vector3d targetField(0.02, 0.12, 0.01);
	setInput continuationInput = makeInput(startField);
	world continuationWorld(continuationInput);
	continuationWorld.setRodStepper();
	FieldContinuationOptions options;
	options.initialStepFraction = 0.2;
	options.maximumStepFraction = 0.35;
	options.easyCorrectorIterations = 10;
	const FieldContinuationResult continuation =
		continuationWorld.continuePlanarContactField(targetField, options);
	if (!continuation.success || continuation.rolledBack ||
		continuation.points.size() < 3 ||
		(continuation.points.back().field - targetField).norm() > 1.0e-14)
	{
		std::cerr << "Fixed-active field continuation failed: success="
				  << continuation.success << ", rollback="
				  << continuation.rolledBack << ", points="
				  << continuation.points.size() << ", attempts="
				  << continuation.attemptedSteps << ", rejected="
				  << continuation.rejectedSteps << ".\n";
		return false;
	}
	for (int point = 1; point < static_cast<int>(continuation.points.size()); ++point)
	{
		const FieldContinuationPoint &current = continuation.points[point];
		if (current.pathFraction <= continuation.points[point - 1].pathFraction ||
			current.contactSetChanged || !current.stability.stable ||
			!sameContactIds(
				continuation.points[0].equilibrium.activeContacts,
				current.equilibrium.activeContacts))
		{
			std::cerr << "Fixed-active continuation point diagnostics failed.\n";
			return false;
		}
	}

	setInput directInput = makeInput(targetField);
	world directWorld(directInput);
	directWorld.setRodStepper();
	const ContactKktEquilibriumResult direct =
		directWorld.solvePlanarContactKktEquilibrium();
	const ContactKktEquilibriumResult &continued =
		continuation.points.back().equilibrium;
	if (!direct.success ||
		!sameContactIds(direct.activeContacts, continued.activeContacts) ||
		positionDifference(direct.state, continued.state) > 1.0e-7 ||
		(direct.multipliers - continued.multipliers).norm() > 1.0e-7)
	{
		std::cerr << "Continuation endpoint disagrees with direct solve: q_error="
				  << positionDifference(direct.state, continued.state)
				  << ", lambda_error="
				  << (direct.multipliers - continued.multipliers).norm() << ".\n";
		return false;
	}

	std::cout << "Fixed-active continuation: points="
			  << continuation.points.size()
			  << ", attempts=" << continuation.attemptedSteps
			  << ", rejected=" << continuation.rejectedSteps
			  << ", endpoint_q_error="
			  << positionDifference(direct.state, continued.state) << '\n';
	return true;
}

bool checkStabilityRollback()
{
	const Vector3d startField(0.0, 0.1, 0.0);
	setInput input = makeInput(startField);
	world simulation(input);
	simulation.setRodStepper();
	const RodState startState = simulation.captureRodState();
	FieldContinuationOptions options;
	options.stabilityTolerance = 1.0e6;
	const FieldContinuationResult result = simulation.continuePlanarContactField(
		Vector3d(0.0, 0.11, 0.0), options);
	const RodState restored = simulation.captureRodState();
	if (result.success || !result.rolledBack || !result.points.empty() ||
		(restored.configuration - startState.configuration).norm() != 0.0 ||
		(simulation.getAppliedField() - startField).norm() != 0.0)
	{
		std::cerr << "Unstable continuation did not roll back exactly.\n";
		return false;
	}
	return true;
}

bool checkContactReleaseEvent()
{
	const Vector3d startField(0.0, 0.1, 0.0);
	const Vector3d targetField(0.0, 0.1, -5.0);
	setInput input = makeInput(startField);
	world simulation(input);
	simulation.setRodStepper();
	FieldContinuationOptions options;
	options.initialStepFraction = 0.15;
	options.maximumStepFraction = 0.25;
	options.easyCorrectorIterations = 10;
	const FieldContinuationResult continuation =
		simulation.continuePlanarContactField(targetField, options);
	int releasedContacts = 0;
	int eventPoints = 0;
	for (const FieldContinuationPoint &point : continuation.points)
	{
		releasedContacts += point.contactsReleased;
		eventPoints += point.contactSetChanged ? 1 : 0;
	}
	if (!continuation.success || continuation.rolledBack ||
		continuation.points.front().equilibrium.activeContacts.empty() ||
		!continuation.points.back().equilibrium.activeContacts.empty() ||
		releasedContacts <= 0 || eventPoints <= 0)
	{
		std::cerr << "Contact-release continuation event failed: success="
				  << continuation.success << ", points="
				  << continuation.points.size() << ", events=" << eventPoints
				  << ", released=" << releasedContacts << ".\n";
		return false;
	}

	setInput directInput = makeInput(targetField);
	world directWorld(directInput);
	directWorld.setRodStepper();
	const ContactKktEquilibriumResult direct =
		directWorld.solvePlanarContactKktEquilibrium();
	const ContactKktEquilibriumResult &continued =
		continuation.points.back().equilibrium;
	if (!direct.success || positionDifference(direct.state, continued.state) > 1.0e-6)
	{
		std::cerr << "Event-continuation endpoint disagrees with direct solve.\n";
		return false;
	}
	std::cout << "Contact-event continuation: points="
			  << continuation.points.size() << ", events=" << eventPoints
			  << ", released=" << releasedContacts
			  << ", rejected=" << continuation.rejectedSteps << '\n';
	return true;
}
}

int main()
{
	if (!checkFixedActiveContinuation() || !checkContactReleaseEvent() ||
		!checkStabilityRollback())
	{
		return 1;
	}
	return 0;
}
