#include "contact/planarBarrierContactForce.h"

#include <cmath>
#include <stdexcept>
#include <vector>

PlanarBarrierContactEvaluation evaluatePlanarBarrierContact(
	const ContactCandidate &candidate,
	double activationDistance,
	double stiffness)
{
	if (!candidate.normal.allFinite() ||
		!candidate.gapHessian.allFinite() ||
		std::abs(candidate.normal.norm() - 1.0) > 1.0e-10)
	{
		throw std::invalid_argument(
			"Barrier contact requires finite gap derivatives and a unit normal");
	}

	const BarrierPotentialEvaluation barrier = evaluateBarrierPotential(
		candidate.gap, activationDistance, stiffness);
	PlanarBarrierContactEvaluation result;
	result.valid = barrier.valid;
	result.active = barrier.active;
	result.energy = barrier.energy;
	if (barrier.active)
	{
		result.residual = barrier.firstDerivative * candidate.normal;
		result.jacobian = barrier.secondDerivative *
			(candidate.normal * candidate.normal.transpose()) +
			barrier.firstDerivative * candidate.gapHessian;
	}
	return result;
}

PlanarBarrierContactForce::PlanarBarrierContactForce(
	elasticRod &m_rod,
	timeStepper &m_stepper,
	const ConfinedDomain &domain,
	double rodRadius,
	double m_activationDistance,
	double m_stiffness,
	double tipSafeDistance)
	: rod(&m_rod),
	  stepper(&m_stepper),
	  detector(domain, rodRadius, m_activationDistance, tipSafeDistance),
	  activationDistance(m_activationDistance),
	  stiffness(m_stiffness)
{
	// Validate the scalar-law parameters at construction time.
	(void)evaluateBarrierPotential(
		m_activationDistance, m_activationDistance, m_stiffness);
}

void PlanarBarrierContactForce::computeFc()
{
	detection = detectCurrent();

	for (const ContactCandidate &candidate : detection.bodyCandidates)
	{
		const PlanarBarrierContactEvaluation contact =
			evaluatePlanarBarrierContact(
				candidate, activationDistance, stiffness);
		if (!contact.valid)
		{
			throw InfeasibleContactError(
				"Planar barrier contact encountered a nonpositive rod-wall gap");
		}
		if (!contact.active)
		{
			continue;
		}

		for (int row = 0; row < 3; ++row)
		{
			const int rowDof = 4 * candidate.rodVertex + row;
			stepper->addForce(rowDof, contact.residual[row]);
			for (int column = 0; column < 3; ++column)
			{
				const int columnDof = 4 * candidate.rodVertex + column;
				stepper->addJacobian(
					rowDof, columnDof, contact.jacobian(row, column));
			}
		}
	}
}

ContactDetectionResult PlanarBarrierContactForce::detectCurrent() const
{
	std::vector<Vector3d> vertices;
	vertices.reserve(rod->nv);
	for (int vertex = 0; vertex < rod->nv; ++vertex)
	{
		vertices.push_back(rod->getVertex(vertex));
	}
	return detector.detect(vertices);
}

const ContactDetectionResult &PlanarBarrierContactForce::lastDetection() const
{
	return detection;
}
