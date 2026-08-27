#include "planning/mechanicsSession.h"

#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
void validateConfig(const MechanicsConfig &config)
{
	if (!std::isfinite(config.rodLength) || config.rodLength <= 0.0 ||
		!std::isfinite(config.rodRadius) || config.rodRadius <= 0.0 ||
		config.numVertices < 4 ||
		!std::isfinite(config.youngModulus) || config.youngModulus <= 0.0 ||
		!std::isfinite(config.poissonRatio) || config.poissonRatio <= -1.0 ||
		config.maximumNewtonIterations <= 0 || !config.gravity.allFinite() ||
		(config.domainType != "planar_slab" &&
		 config.domainType != "spherical_shell" &&
		 config.domainType != "spherical_obstacle" &&
		 config.domainType != "double_spherical_obstacle") ||
		!std::isfinite(config.planeHalfThickness) ||
		config.planeHalfThickness <= config.rodRadius ||
		!config.shellCenter.allFinite() || !std::isfinite(config.shellRadius) ||
		!std::isfinite(config.shellMinusThickness) ||
		!std::isfinite(config.shellPlusThickness) || config.shellRadius <= 0.0 ||
		config.shellMinusThickness < 0.0 || config.shellPlusThickness < 0.0 ||
		config.shellRadius - config.shellMinusThickness <= config.rodRadius ||
		!config.cavityCenter.allFinite() || !std::isfinite(config.cavityRadius) ||
		!config.obstacleCenter.allFinite() ||
		!std::isfinite(config.obstacleRadius) || config.cavityRadius <= 0.0 ||
		config.obstacleRadius <= 0.0 ||
		(config.obstacleCenter - config.cavityCenter).norm() +
			config.obstacleRadius >= config.cavityRadius ||
		!config.secondObstacleCenter.allFinite() ||
		!std::isfinite(config.secondObstacleRadius) ||
		config.secondObstacleRadius <= 0.0 ||
		(config.secondObstacleCenter - config.cavityCenter).norm() +
			config.secondObstacleRadius >= config.cavityRadius ||
		!std::isfinite(config.barrierDistance) || config.barrierDistance <= 0.0 ||
		!std::isfinite(config.barrierStiffness) || config.barrierStiffness <= 0.0 ||
		!std::isfinite(config.tipSafeDistance) || config.tipSafeDistance < 0.0 ||
		!std::isfinite(config.tipDipoleMoment) ||
		!std::isfinite(config.insertionStiffness) ||
		config.insertionStiffness <= 0.0 || !config.insertionAxis.allFinite() ||
		config.insertionAxis.norm() == 0.0 ||
		!std::isfinite(config.initialActuation.xi) ||
		config.initialActuation.xi < 0.0 ||
		!config.initialActuation.field.allFinite())
	{
		throw std::invalid_argument("Invalid mechanics-session configuration");
	}
}

setInput makeInput(const MechanicsConfig &config)
{
	setInput input;
	input.GetScalarOpt("RodLength") = config.rodLength;
	input.GetScalarOpt("rodRadius") = config.rodRadius;
	input.GetIntOpt("numVertices") = config.numVertices;
	input.GetScalarOpt("youngM") = config.youngModulus;
	input.GetScalarOpt("Poisson") = config.poissonRatio;
	input.GetIntOpt("maxIter") = config.maximumNewtonIterations;
	input.GetVecOpt("gVector") = config.gravity;
	input.GetScalarOpt("thickness") = config.planeHalfThickness;
	input.GetVecOpt("shellCenter") = config.shellCenter;
	input.GetScalarOpt("shellRadius") = config.shellRadius;
	input.GetScalarOpt("shellMinusThickness") = config.shellMinusThickness;
	input.GetScalarOpt("shellPlusThickness") = config.shellPlusThickness;
	input.GetVecOpt("cavityCenter") = config.cavityCenter;
	input.GetScalarOpt("cavityRadius") = config.cavityRadius;
	input.GetVecOpt("obstacleCenter") = config.obstacleCenter;
	input.GetScalarOpt("obstacleRadius") = config.obstacleRadius;
	input.GetVecOpt("secondObstacleCenter") = config.secondObstacleCenter;
	input.GetScalarOpt("secondObstacleRadius") = config.secondObstacleRadius;
	input.GetScalarOpt("dBar") = config.barrierDistance;
	input.GetScalarOpt("stiffness") = config.barrierStiffness;
	input.GetScalarOpt("tipSafeDistance") = config.tipSafeDistance;
	if (config.domainType == "planar_slab")
	{
		input.GetStringOpt("contactModel") = "planar_barrier";
	}
	else if (config.domainType == "spherical_shell")
	{
		input.GetStringOpt("contactModel") = "spherical_shell_barrier";
	}
	else if (config.domainType == "spherical_obstacle")
	{
		input.GetStringOpt("contactModel") = "spherical_obstacle_barrier";
	}
	else
	{
		input.GetStringOpt("contactModel") =
			"double_spherical_obstacle_barrier";
	}
	input.GetStringOpt("magneticModel") = "axial_tip";
	input.GetScalarOpt("tipDipoleMoment") = config.tipDipoleMoment;
	input.GetStringOpt("insertionModel") = "proximal_guide";
	input.GetScalarOpt("insertionCoordinate") = config.initialActuation.xi;
	input.GetScalarOpt("insertionStiffness") = config.insertionStiffness;
	input.GetVecOpt("insertionAxis") = config.insertionAxis.normalized();
	input.GetVecOpt("baVector") = config.initialActuation.field;
	return input;
}

ContactKktEquilibriumResult equilibriumFromPlannerState(
	const PlannerState &state)
{
	if (state.multipliers.size() !=
		static_cast<int>(state.activeContacts.size()))
	{
		throw std::invalid_argument(
			"PlannerState contact and multiplier counts must agree");
	}
	ContactKktEquilibriumResult equilibrium;
	equilibrium.success = true;
	equilibrium.state = state.rodState;
	equilibrium.activeContacts = state.activeContacts;
	equilibrium.multipliers = state.multipliers;
	equilibrium.stationarityNorm = state.stationarityNorm;
	equilibrium.activeConstraintNorm = state.activeConstraintNorm;
	equilibrium.complementarityNorm = state.complementarityNorm;
	equilibrium.minimumBodyGap = state.minimumBodyGap;
	equilibrium.tipClearance = state.tipClearance;
	equilibrium.tipSafe = state.tipSafe;
	return equilibrium;
}
}

MechanicsSession::MechanicsSession(const MechanicsConfig &config)
{
	validateConfig(config);
	setInput input = makeInput(config);
	simulation = std::make_unique<world>(input);
	simulation->setRodStepper();
}

MechanicsSession::~MechanicsSession() = default;

PlannerState MechanicsSession::solveInitialState()
{
	const ContactKktEquilibriumResult equilibrium =
		simulation->solvePlanarContactKktEquilibrium();
	if (!equilibrium.success)
	{
		throw std::runtime_error("Initial planar KKT equilibrium solve failed");
	}
	const ContactStabilityAnalysis stability =
		simulation->analyzePlanarContactStability(equilibrium);
	if (!stability.valid)
	{
		throw std::runtime_error("Initial equilibrium stability analysis failed");
	}
	return simulation->capturePlannerState(equilibrium, stability);
}

LocalSteeringResult MechanicsSession::evaluateLocalSteering(
	const PlannerState &state)
{
	if (state.rodState.configuration.size() < 3)
	{
		throw std::invalid_argument("PlannerState configuration has no distal tip");
	}
	simulation->restorePlannerState(state);
	const ContactKktEquilibriumResult equilibrium =
		equilibriumFromPlannerState(state);
	const ActuationEquilibriumSensitivity sensitivity =
		simulation->computePlanarContactActuationSensitivity(equilibrium);
	if (!sensitivity.success)
	{
		throw std::runtime_error("Local actuation sensitivity solve failed");
	}
	const int tipDof = state.rodState.configuration.size() - 3;
	LocalSteeringResult result;
	result.tipPosition = state.rodState.configuration.segment<3>(tipDof);
	result.tipActuationJacobian =
		sensitivity.configurationDerivative.block<3, 4>(tipDof, 0);
	result.linearResidualNorm = sensitivity.linearResidualNorm;
	return result;
}

ContinuationEdgeResult MechanicsSession::attemptContinuation(
	const PlannerState &startState,
	const Actuation &targetActuation,
	const FieldContinuationOptions &options)
{
	const ActuationContinuationResult continuation =
		simulation->continuePlanarContactActuation(
			startState, targetActuation, options);

	ContinuationEdgeResult result;
	result.success = continuation.success;
	result.rolledBack = continuation.rolledBack;
	result.attemptedSteps = continuation.attemptedSteps;
	result.rejectedSteps = continuation.rejectedSteps;
	result.storedPoints = static_cast<int>(continuation.points.size());
	result.state = startState;
	for (const ActuationContinuationPoint &point : continuation.points)
	{
		result.contactsAdded += point.contactsAdded;
		result.contactsReleased += point.contactsReleased;
		result.reachedPathFraction = std::max(
			result.reachedPathFraction, point.pathFraction);
		if (point.stability.valid)
		{
			result.minimumStabilityMargin = std::min(
				result.minimumStabilityMargin,
				point.stability.minimumEigenvalue);
		}
	}
	result.failureReason = continuation.failureReason;
	if (continuation.success && !continuation.points.empty())
	{
		const ActuationContinuationPoint &finalPoint = continuation.points.back();
		result.state = simulation->capturePlannerState(
			finalPoint.equilibrium, finalPoint.stability);
	}
	return result;
}
