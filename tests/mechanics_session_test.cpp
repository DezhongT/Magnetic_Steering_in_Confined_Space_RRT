#include "planning/mechanicsSession.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
MechanicsConfig makeConfig(const Actuation *actuation = nullptr)
{
	MechanicsConfig config;
	config.gravity = Vector3d(0.0, 0.0, 0.22);
	config.initialActuation.xi = 1.0e-3;
	config.initialActuation.field = Vector3d(0.0, 0.1, 0.0);
	if (actuation != nullptr)
	{
		config.initialActuation = *actuation;
	}
	return config;
}

bool sameState(const PlannerState &left, const PlannerState &right)
{
	return (left.rodState.configuration - right.rodState.configuration).norm() == 0.0 &&
		(left.rodState.previousConfiguration -
		 right.rodState.previousConfiguration).norm() == 0.0 &&
		std::abs(left.actuation.xi - right.actuation.xi) == 0.0 &&
		(left.actuation.field - right.actuation.field).norm() == 0.0 &&
		(left.multipliers - right.multipliers).norm() == 0.0;
}
}

int main()
{
	MechanicsSession session(makeConfig());
	const PlannerState start = session.solveInitialState();
	if (!start.tipSafe || !std::isfinite(start.stabilityMargin) ||
		start.stabilityMargin <= 0.0 || !std::isfinite(start.tipClearance) ||
		start.tipClearance <= 0.0 || !std::isfinite(start.minimumBodyGap))
	{
		std::cerr << "MechanicsSession returned an invalid initial state.\n";
		return 1;
	}
	const LocalSteeringResult steering = session.evaluateLocalSteering(start);
	if ((steering.tipPosition - start.rodState.configuration.tail<3>()).norm() != 0.0 ||
		steering.tipActuationJacobian.rows() != 3 ||
		steering.tipActuationJacobian.cols() != 4 ||
		steering.linearResidualNorm > 1.0e-10)
	{
		std::cerr << "MechanicsSession local steering query failed.\n";
		return 1;
	}

	constexpr double derivativeStep = 1.0e-6;
	Matrix<double, 3, 4> finiteDifference;
	for (int parameter = 0; parameter < 4; ++parameter)
	{
		Actuation plus = start.actuation;
		Actuation minus = start.actuation;
		if (parameter == 0)
		{
			plus.xi += derivativeStep;
			minus.xi -= derivativeStep;
		}
		else
		{
			plus.field[parameter - 1] += derivativeStep;
			minus.field[parameter - 1] -= derivativeStep;
		}
		MechanicsSession plusSession(makeConfig(&plus));
		MechanicsSession minusSession(makeConfig(&minus));
		const PlannerState plusState = plusSession.solveInitialState();
		const PlannerState minusState = minusSession.solveInitialState();
		finiteDifference.col(parameter) =
			(plusState.rodState.configuration.tail<3>() -
			 minusState.rodState.configuration.tail<3>()) /
			(2.0 * derivativeStep);
	}
	const double steeringError =
		(finiteDifference - steering.tipActuationJacobian).norm() /
		std::max(1.0, finiteDifference.norm());
	if (steeringError > 1.0e-4)
	{
		std::cerr << "Tip-actuation Jacobian finite difference failed: error="
				  << steeringError << ".\n";
		return 1;
	}

	Actuation target = start.actuation;
	target.xi += 2.0e-4;
	target.field += Vector3d(0.01, 0.01, 0.005);
	FieldContinuationOptions options;
	options.initialStepFraction = 0.25;
	options.maximumStepFraction = 0.4;
	options.easyCorrectorIterations = 10;
	const ContinuationEdgeResult accepted =
		session.attemptContinuation(start, target, options);
	if (!accepted.success || accepted.rolledBack || accepted.storedPoints < 3 ||
		accepted.state.actuation.xi < start.actuation.xi ||
		std::abs(accepted.state.actuation.xi - target.xi) > 0.0 ||
		(accepted.state.actuation.field - target.field).norm() > 0.0 ||
		accepted.minimumStabilityMargin <= 0.0)
	{
		std::cerr << "MechanicsSession accepted-edge attempt failed.\n";
		return 1;
	}

	FieldContinuationOptions rejectingOptions = options;
	rejectingOptions.stabilityTolerance = start.stabilityMargin + 1.0;
	const ContinuationEdgeResult rejected =
		session.attemptContinuation(start, target, rejectingOptions);
	if (rejected.success || !rejected.rolledBack || !sameState(rejected.state, start))
	{
		std::cerr << "MechanicsSession rejected-edge rollback failed.\n";
		return 1;
	}
	const ContinuationEdgeResult recovered =
		session.attemptContinuation(start, target, options);
	if (!recovered.success)
	{
		std::cerr << "MechanicsSession could not continue after rollback.\n";
		return 1;
	}

	std::cout << "MechanicsSession: points=" << accepted.storedPoints
			  << ", attempts=" << accepted.attemptedSteps
			  << ", minimum_stability=" << accepted.minimumStabilityMargin
			  << ", tip_jacobian_error=" << steeringError
			  << ", rollback_verified=1\n";
	return 0;
}
