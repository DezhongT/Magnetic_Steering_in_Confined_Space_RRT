#include "geometry/sphericalShellDomain.h"
#include "contact/planarBarrierContactForce.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool checkDerivatives(
	const SphericalShellDomain &domain,
	const Vector3d &point,
	int boundaryId)
{
	constexpr double gradientStep = 1.0e-6;
	constexpr double hessianStep = 2.0e-6;
	const DomainQuery query = domain.queryBoundary(point, boundaryId);
	Vector3d finiteDifferenceGradient;
	Matrix3d finiteDifferenceHessian;
	for (int component = 0; component < 3; ++component)
	{
		Vector3d offset = Vector3d::Zero();
		offset[component] = gradientStep;
		finiteDifferenceGradient[component] =
			(domain.queryBoundary(point + offset, boundaryId).clearance -
			 domain.queryBoundary(point - offset, boundaryId).clearance) /
			(2.0 * gradientStep);

		offset[component] = hessianStep;
		finiteDifferenceHessian.col(component) =
			(domain.queryBoundary(point + offset, boundaryId).clearanceNormal -
			 domain.queryBoundary(point - offset, boundaryId).clearanceNormal) /
			(2.0 * hessianStep);
	}
	return (finiteDifferenceGradient - query.clearanceNormal).norm() < 1.0e-9 &&
		(finiteDifferenceHessian - query.clearanceHessian).norm() < 1.0e-9 &&
		(query.clearanceHessian - query.clearanceHessian.transpose()).norm() <
			1.0e-14;
}

ContactCandidate candidateAt(
	const SphericalShellDomain &domain,
	const Vector3d &point,
	int boundaryId)
{
	const DomainQuery query = domain.queryBoundary(point, boundaryId);
	ContactCandidate candidate;
	candidate.gap = query.clearance;
	candidate.normal = query.clearanceNormal;
	candidate.gapHessian = query.clearanceHessian;
	return candidate;
}

bool checkCurvedBarrierJacobian(
	const SphericalShellDomain &domain,
	const Vector3d &point)
{
	constexpr double step = 1.0e-7;
	constexpr double activationDistance = 5.0e-2;
	constexpr double stiffness = 10.0;
	const Vector3d direction = Vector3d(0.3, -0.4, 0.2).normalized();
	const PlanarBarrierContactEvaluation evaluation =
		evaluatePlanarBarrierContact(
			candidateAt(domain, point, 1), activationDistance, stiffness);
	const Vector3d finiteDifference =
		(evaluatePlanarBarrierContact(
			 candidateAt(domain, point + step * direction, 1),
			 activationDistance, stiffness).residual -
		 evaluatePlanarBarrierContact(
			 candidateAt(domain, point - step * direction, 1),
			 activationDistance, stiffness).residual) /
		(2.0 * step);
	return evaluation.active &&
		(finiteDifference - evaluation.jacobian * direction).norm() /
			std::max(1.0, finiteDifference.norm()) < 1.0e-8;
}
}

int main()
{
	const Vector3d center(0.2, -0.3, 0.4);
	const SphericalShellDomain domain(center, 2.0, 0.25, 0.4);
	const Vector3d radial = Vector3d(1.0, 2.0, -1.0).normalized();
	const Vector3d innerPoint = center + 1.9 * radial;
	const Vector3d outerPoint = center + 2.2 * radial;
	const DomainQuery inner = domain.query(innerPoint);
	const DomainQuery outer = domain.query(outerPoint);
	if (std::abs(domain.innerRadius() - 1.75) > 1.0e-14 ||
		std::abs(domain.outerRadius() - 2.4) > 1.0e-14 ||
		inner.boundaryId != 0 || std::abs(inner.clearance - 0.15) > 1.0e-14 ||
		(inner.clearanceNormal - radial).norm() > 1.0e-14 ||
		(inner.closestPoint - (center + 1.75 * radial)).norm() > 1.0e-14 ||
		outer.boundaryId != 1 || std::abs(outer.clearance - 0.2) > 1.0e-14 ||
		(outer.clearanceNormal + radial).norm() > 1.0e-14 ||
		(outer.closestPoint - (center + 2.4 * radial)).norm() > 1.0e-14 ||
		!inner.admissible || !outer.admissible)
	{
		std::cerr << "Spherical-shell clearance or boundary selection failed.\n";
		return 1;
	}

	const Vector3d derivativePoint = center + Vector3d(1.1, -0.7, 1.3);
	if (!checkDerivatives(domain, derivativePoint, 0) ||
		!checkDerivatives(domain, derivativePoint, 1))
	{
		std::cerr << "Spherical-shell gradient/Hessian finite difference failed.\n";
		return 1;
	}
	if (!checkCurvedBarrierJacobian(domain, center + 2.37 * radial))
	{
		std::cerr << "Curved barrier Jacobian finite difference failed.\n";
		return 1;
	}
	const DomainQuery forcedInner = domain.queryBoundary(outerPoint, 0);
	const DomainQuery forcedOuter = domain.queryBoundary(innerPoint, 1);
	if (forcedInner.boundaryId != 0 || forcedOuter.boundaryId != 1 ||
		(forcedInner.clearanceHessian +
		 forcedOuter.clearanceHessian * (1.9 / 2.2)).norm() > 1.0e-12)
	{
		std::cerr << "Spherical-shell stable boundary query failed.\n";
		return 1;
	}

	bool rejectedCenter = false;
	bool rejectedThickness = false;
	try
	{
		(void)domain.query(center);
	}
	catch (const std::domain_error &)
	{
		rejectedCenter = true;
	}
	try
	{
		SphericalShellDomain invalid(center, 1.0, 1.0, 0.2);
		(void)invalid;
	}
	catch (const std::invalid_argument &)
	{
		rejectedThickness = true;
	}
	if (!rejectedCenter || !rejectedThickness)
	{
		std::cerr << "Spherical-shell singular/invalid input was accepted.\n";
		return 1;
	}

	std::cout << "Spherical-shell clearance, gradient, and Hessian checks passed.\n";
	return 0;
}
