#include "contact/barrierPotential.h"
#include "contact/planarBarrierContactForce.h"
#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
double relativeError(double left, double right)
{
	return std::abs(left - right) /
		std::max(1.0, std::max(std::abs(left), std::abs(right)));
}

bool checkScalarLaw()
{
	constexpr double activation = 0.02;
	constexpr double stiffness = 37.0;
	constexpr double step = 1.0e-7;
	for (const double gap : {0.002, 0.007, 0.015})
	{
		const BarrierPotentialEvaluation center =
			evaluateBarrierPotential(gap, activation, stiffness);
		const BarrierPotentialEvaluation plus =
			evaluateBarrierPotential(gap + step, activation, stiffness);
		const BarrierPotentialEvaluation minus =
			evaluateBarrierPotential(gap - step, activation, stiffness);
		const double firstDifference = (plus.energy - minus.energy) / (2.0 * step);
		const double secondDifference =
			(plus.firstDerivative - minus.firstDerivative) / (2.0 * step);
		if (!center.valid || !center.active ||
			relativeError(firstDifference, center.firstDerivative) > 2.0e-8 ||
			relativeError(secondDifference, center.secondDerivative) > 2.0e-8)
		{
			std::cerr << "Barrier scalar derivative check failed at gap="
					  << gap << ".\n";
			return false;
		}
	}

	const BarrierPotentialEvaluation cutoff =
		evaluateBarrierPotential(activation, activation, stiffness);
	const BarrierPotentialEvaluation outside =
		evaluateBarrierPotential(0.03, activation, stiffness);
	const BarrierPotentialEvaluation infeasible =
		evaluateBarrierPotential(0.0, activation, stiffness);
	if (!cutoff.valid || cutoff.active || cutoff.energy != 0.0 ||
		!outside.valid || outside.active || infeasible.valid)
	{
		std::cerr << "Barrier cutoff or feasibility convention failed.\n";
		return false;
	}
	return true;
}

PlanarBarrierContactEvaluation evaluateAtPosition(
	const Vector3d &position,
	const Vector3d &normal,
	double baseGap,
	const Vector3d &basePosition,
	double activation,
	double stiffness)
{
	ContactCandidate candidate;
	candidate.gap = baseGap + normal.dot(position - basePosition);
	candidate.normal = normal;
	return evaluatePlanarBarrierContact(candidate, activation, stiffness);
}

bool checkPlanarAssembly()
{
	constexpr double activation = 0.02;
	constexpr double stiffness = 37.0;
	constexpr double baseGap = 0.007;
	constexpr double step = 1.0e-7;
	const Vector3d normal = Vector3d(1.0, -2.0, 3.0).normalized();
	const Vector3d basePosition(0.3, -0.1, 0.8);
	const Vector3d direction = Vector3d(-0.4, 0.7, 0.2).normalized();
	const PlanarBarrierContactEvaluation center = evaluateAtPosition(
		basePosition, normal, baseGap, basePosition, activation, stiffness);
	const PlanarBarrierContactEvaluation plus = evaluateAtPosition(
		basePosition + step * direction, normal, baseGap, basePosition,
		activation, stiffness);
	const PlanarBarrierContactEvaluation minus = evaluateAtPosition(
		basePosition - step * direction, normal, baseGap, basePosition,
		activation, stiffness);
	const double energyDerivative = (plus.energy - minus.energy) / (2.0 * step);
	const Vector3d residualDerivative =
		(plus.residual - minus.residual) / (2.0 * step);
	const double energyScale = std::max(
		1.0, std::max(std::abs(energyDerivative),
						 std::abs(center.residual.dot(direction))));
	const double jacobianScale = std::max(
		1.0, std::max(residualDerivative.norm(),
						 (center.jacobian * direction).norm()));
	if (!center.valid || !center.active ||
		std::abs(energyDerivative - center.residual.dot(direction)) /
			energyScale > 2.0e-8 ||
		(residualDerivative - center.jacobian * direction).norm() /
			jacobianScale > 2.0e-8)
	{
		std::cerr << "Planar barrier residual/Jacobian finite difference failed.\n";
		return false;
	}
	return true;
}

bool checkWorldSelection()
{
	setInput noContactInput;
	noContactInput.GetStringOpt("contactModel") = "none";
	world noContact(noContactInput);
	noContact.setRodStepper();

	setInput barrierInput;
	barrierInput.GetStringOpt("contactModel") = "planar_barrier";
	world barrierContact(barrierInput);
	barrierContact.setRodStepper();

	RodState noContactState = noContact.captureRodState();
	RodState barrierState = barrierContact.captureRodState();
	constexpr int vertex = 2;
	constexpr double zPosition = 0.088;
	noContactState.configuration[4 * vertex + 2] = zPosition;
	barrierState.configuration[4 * vertex + 2] = zPosition;
	noContact.restoreRodState(noContactState);
	barrierContact.restoreRodState(barrierState);

	const StaticEvaluation withoutBarrier = noContact.evaluateStaticSystem();
	const StaticEvaluation withBarrier = barrierContact.evaluateStaticSystem();
	const VectorXd residualDifference =
		withBarrier.residual - withoutBarrier.residual;
	const BarrierPotentialEvaluation expected = evaluateBarrierPotential(
		0.1 - zPosition - 0.01, 0.005, 1.0e4);
	if (!expected.active ||
		relativeError(residualDifference.norm(),
			std::abs(expected.firstDerivative)) > 1.0e-10)
	{
		std::cerr << "Selectable world contact residual assembly failed.\n";
		return false;
	}

	Eigen::Index activeDof = 0;
	residualDifference.cwiseAbs().maxCoeff(&activeDof);
	VectorXd direction = VectorXd::Zero(residualDifference.size());
	direction[activeDof] = 1.0;
	const VectorXd jacobianDifference =
		withBarrier.multiplyJacobian(direction) -
		withoutBarrier.multiplyJacobian(direction);
	if (relativeError(jacobianDifference.norm(), expected.secondDerivative) >
		1.0e-10)
	{
		std::cerr << "Selectable world contact Jacobian assembly failed.\n";
		return false;
	}
	return true;
}
}

int main()
{
	if (!checkScalarLaw() || !checkPlanarAssembly() || !checkWorldSelection())
	{
		return 1;
	}
	std::cout << "Barrier law and planar contact assembly checks passed.\n";
	return 0;
}
