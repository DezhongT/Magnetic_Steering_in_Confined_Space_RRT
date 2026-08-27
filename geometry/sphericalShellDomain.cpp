#include "geometry/sphericalShellDomain.h"

#include <cmath>
#include <stdexcept>

SphericalShellDomain::SphericalShellDomain(
	const Vector3d &m_center,
	double m_referenceRadius,
	double m_minusThickness,
	double m_plusThickness)
	: center(m_center),
	  referenceRadius(m_referenceRadius),
	  minusThickness(m_minusThickness),
	  plusThickness(m_plusThickness)
{
	if (!center.allFinite() || !std::isfinite(referenceRadius) ||
		!std::isfinite(minusThickness) || !std::isfinite(plusThickness) ||
		!(referenceRadius > 0.0) || minusThickness < 0.0 ||
		plusThickness < 0.0 || !(innerRadius() > 0.0))
	{
		throw std::invalid_argument(
			"Spherical shell requires finite parameters and a positive inner radius");
	}
}

double SphericalShellDomain::innerRadius() const
{
	return referenceRadius - minusThickness;
}

double SphericalShellDomain::outerRadius() const
{
	return referenceRadius + plusThickness;
}

DomainQuery SphericalShellDomain::query(const Vector3d &point) const
{
	const double radius = (point - center).norm();
	const double innerClearance = radius - innerRadius();
	const double outerClearance = outerRadius() - radius;
	return queryBoundary(point, innerClearance <= outerClearance ? 0 : 1);
}

DomainQuery SphericalShellDomain::queryBoundary(
	const Vector3d &point,
	int boundaryId) const
{
	if (boundaryId != 0 && boundaryId != 1)
	{
		throw std::invalid_argument(
			"Spherical shell boundary ID must be 0 or 1");
	}
	if (!point.allFinite())
	{
		throw std::invalid_argument(
			"Spherical shell query point must be finite");
	}
	const Vector3d offset = point - center;
	const double radius = offset.norm();
	if (!(radius > 0.0))
	{
		throw std::domain_error(
			"Spherical-shell clearance derivatives are undefined at the center");
	}

	const Vector3d radial = offset / radius;
	const Matrix3d radialHessian =
		(Matrix3d::Identity() - radial * radial.transpose()) / radius;
	DomainQuery result;
	result.boundaryId = boundaryId;
	if (boundaryId == 0)
	{
		result.clearance = radius - innerRadius();
		result.closestPoint = center + innerRadius() * radial;
		result.clearanceNormal = radial;
		result.clearanceHessian = radialHessian;
	}
	else
	{
		result.clearance = outerRadius() - radius;
		result.closestPoint = center + outerRadius() * radial;
		result.clearanceNormal = -radial;
		result.clearanceHessian = -radialHessian;
	}
	result.admissible = result.clearance >= 0.0;
	return result;
}
