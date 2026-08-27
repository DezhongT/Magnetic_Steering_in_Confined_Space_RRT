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
	input.GetStringOpt("contactModel") = "spherical_shell_barrier";
	input.GetStringOpt("magneticModel") = "axial_tip";
	input.GetStringOpt("insertionModel") = "proximal_guide";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 300;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 100.0;
	input.GetScalarOpt("tipDipoleMoment") = 1.0e-3;
	input.GetScalarOpt("insertionCoordinate") = actuation.xi;
	input.GetScalarOpt("insertionStiffness") = 1.0e3;
	input.GetVecOpt("insertionAxis") = Vector3d::UnitX();
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.5);
	input.GetVecOpt("baVector") = actuation.field;
	input.GetVecOpt("shellCenter") = Vector3d(0.0, 0.0, 1.0);
	input.GetScalarOpt("shellRadius") = 1.0;
	input.GetScalarOpt("shellMinusThickness") = 0.1;
	input.GetScalarOpt("shellPlusThickness") = 1.0;
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

bool checkEquilibriumAndSensitivity(
	world &simulation,
	const Actuation &baseActuation,
	ContactKktEquilibriumResult &base,
	ContactStabilityAnalysis &stability)
{
	base = simulation.solvePlanarContactKktEquilibrium();
	if (!base.success || base.rolledBack || base.activeContacts.size() != 1 ||
		base.activeContacts[0].boundaryId != 0 ||
		base.activeContacts[0].gapHessian.norm() <= 0.0 ||
		base.multipliers.size() != 1 || base.multipliers[0] <= 0.0 ||
		base.minimumBodyGap < -1.0e-8 || base.tipClearance <= 0.0 ||
		!base.tipSafe)
	{
		std::cerr << "Curved KKT equilibrium failed: success=" << base.success
			<< ", contacts=" << base.activeContacts.size()
			<< ", gap=" << base.minimumBodyGap
			<< ", tip=" << base.tipClearance
			<< ", multiplier=" << base.minimumMultiplier << ".\n";
		return false;
	}

	stability = simulation.analyzePlanarContactStability(base);
	const ContactEquilibriumSensitivity analytic =
		simulation.computePlanarContactFieldSensitivity(base);
	if (!stability.valid || !stability.stable ||
		stability.minimumEigenvalue <= 0.0 ||
		stability.nullSpaceResidualNorm > 1.0e-10 ||
		stability.hessianSymmetryError > 1.0e-10 || !analytic.success)
	{
		std::cerr << "Curved stability or sensitivity solve failed.\n";
		return false;
	}

	constexpr double step = 1.0e-4;
	MatrixXd configurationDifference = MatrixXd::Zero(
		base.state.configuration.size(), 3);
	MatrixXd multiplierDifference = MatrixXd::Zero(base.multipliers.size(), 3);
	for (int component = 0; component < 3; ++component)
	{
		Vector3d plusField = baseActuation.field;
		Vector3d minusField = baseActuation.field;
		plusField[component] += step;
		minusField[component] -= step;
		simulation.setAppliedField(plusField);
		const ContactKktEquilibriumResult plus =
			simulation.correctPlanarContactKktEquilibrium(base);
		simulation.setAppliedField(minusField);
		const ContactKktEquilibriumResult minus =
			simulation.correctPlanarContactKktEquilibrium(base);
		if (!plus.success || !minus.success ||
			!sameContactIds(base.activeContacts, plus.activeContacts) ||
			!sameContactIds(base.activeContacts, minus.activeContacts))
		{
			std::cerr << "Curved finite-difference correction failed.\n";
			return false;
		}
		configurationDifference.col(component) =
			(plus.state.configuration - minus.state.configuration) / (2.0 * step);
		multiplierDifference.col(component) =
			(plus.multipliers - minus.multipliers) / (2.0 * step);
	}
	simulation.setAppliedField(baseActuation.field);
	const double configurationError =
		(analytic.configurationDerivative - configurationDifference).norm() /
		std::max({1.0, analytic.configurationDerivative.norm(),
			configurationDifference.norm()});
	const double multiplierError =
		(analytic.multiplierDerivative - multiplierDifference).norm() /
		std::max({1.0, analytic.multiplierDerivative.norm(),
			multiplierDifference.norm()});
	if (configurationError > 1.0e-4 || multiplierError > 1.0e-4)
	{
		std::cerr << "Curved KKT sensitivity finite difference failed: q_error="
			<< configurationError << ", lambda_error=" << multiplierError << ".\n";
		return false;
	}
	std::cout << "Curved KKT: vertex=" << base.activeContacts[0].rodVertex
		<< ", multiplier=" << base.multipliers[0]
		<< ", stability=" << stability.minimumEigenvalue
		<< ", q_error=" << configurationError
		<< ", lambda_error=" << multiplierError << '\n';
	return true;
}

bool checkContinuationAndRollback(
	world &simulation,
	const ContactKktEquilibriumResult &base,
	const ContactStabilityAnalysis &stability,
	const Actuation &target)
{
	const PlannerState start = simulation.capturePlannerState(base, stability);
	FieldContinuationOptions options;
	options.initialStepFraction = 0.2;
	options.maximumStepFraction = 0.35;
	options.easyCorrectorIterations = 10;
	const ActuationContinuationResult continuation =
		simulation.continuePlanarContactActuation(start, target, options);
	if (!continuation.success || continuation.rolledBack ||
		continuation.points.size() < 3 ||
		!sameContactIds(
			base.activeContacts,
			continuation.points.back().equilibrium.activeContacts) ||
		(continuation.points.back().actuation.field - target.field).norm() > 1.0e-14 ||
		std::abs(continuation.points.back().actuation.xi - target.xi) > 1.0e-14)
	{
		std::cerr << "Curved combined-actuation continuation failed.\n";
		return false;
	}

	setInput directInput = makeInput(target);
	world direct(directInput);
	direct.setRodStepper();
	const ContactKktEquilibriumResult directResult =
		direct.solvePlanarContactKktEquilibrium();
	const ContactKktEquilibriumResult &continued =
		continuation.points.back().equilibrium;
	if (!directResult.success ||
		!sameContactIds(directResult.activeContacts, continued.activeContacts) ||
		positionDifference(directResult.state, continued.state) > 1.0e-6 ||
		(directResult.multipliers - continued.multipliers).norm() > 1.0e-6)
	{
		std::cerr << "Curved continuation endpoint disagrees with direct solve.\n";
		return false;
	}

	simulation.restorePlannerState(start);
	FieldContinuationOptions rejectingOptions = options;
	rejectingOptions.stabilityTolerance = 1.0e6;
	const RodState before = simulation.captureRodState();
	const ActuationContinuationResult rejected =
		simulation.continuePlanarContactActuation(start, target, rejectingOptions);
	const RodState after = simulation.captureRodState();
	if (rejected.success || !rejected.rolledBack || !rejected.points.empty() ||
		(after.configuration - before.configuration).norm() != 0.0 ||
		std::abs(simulation.getActuation().xi - start.actuation.xi) > 0.0 ||
		(simulation.getAppliedField() - start.actuation.field).norm() != 0.0)
	{
		std::cerr << "Curved continuation rollback was not exact.\n";
		return false;
	}
	std::cout << "Curved continuation: points=" << continuation.points.size()
		<< ", attempts=" << continuation.attemptedSteps
		<< ", rejected=" << continuation.rejectedSteps << '\n';
	return true;
}
}

int main()
{
	Actuation baseActuation;
	baseActuation.xi = 1.0e-3;
	baseActuation.field = Vector3d(0.0, 0.1, 0.0);
	setInput input = makeInput(baseActuation);
	world simulation(input);
	simulation.setRodStepper();
	ContactKktEquilibriumResult base;
	ContactStabilityAnalysis stability;
	if (!checkEquilibriumAndSensitivity(
			simulation, baseActuation, base, stability))
	{
		return 1;
	}
	Actuation target;
	target.xi = 1.2e-3;
	target.field = Vector3d(0.01, 0.11, 0.005);
	if (!checkContinuationAndRollback(simulation, base, stability, target))
	{
		return 1;
	}
	return 0;
}
