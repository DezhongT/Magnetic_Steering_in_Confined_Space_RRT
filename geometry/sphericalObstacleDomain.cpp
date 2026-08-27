#include "geometry/sphericalObstacleDomain.h"

#include <cmath>
#include <stdexcept>

SphericalObstacleDomain::SphericalObstacleDomain(
	const Vector3d &m_outerCenter,
	double m_outerRadius,
	const Vector3d &m_obstacleCenter,
	double m_obstacleRadius)
	: outerCenter(m_outerCenter),
	  outerRadius(m_outerRadius),
	  obstacleCenter(m_obstacleCenter),
	  obstacleRadius(m_obstacleRadius)
{
	if (!outerCenter.allFinite() || !obstacleCenter.allFinite() ||
		!std::isfinite(outerRadius) || !std::isfinite(obstacleRadius) ||
		!(outerRadius > 0.0) || !(obstacleRadius > 0.0) ||
		(obstacleCenter - outerCenter).norm() + obstacleRadius >= outerRadius)
	{
		throw std::invalid_argument(
			"Spherical obstacle must be finite, positive, and strictly inside its cavity");
	}
}

DomainQuery SphericalObstacleDomain::query(const Vector3d &point) const
{
	if (!point.allFinite())
	{
		throw std::invalid_argument("Spherical-obstacle query point must be finite");
	}
	const double outerClearance = outerRadius - (point - outerCenter).norm();
	const double obstacleClearance =
		(point - obstacleCenter).norm() - obstacleRadius;
	return queryBoundary(point, outerClearance <= obstacleClearance ? 0 : 1);
}

DomainQuery SphericalObstacleDomain::queryBoundary(
	const Vector3d &point,
	int boundaryId) const
{
	if (boundaryId != 0 && boundaryId != 1)
	{
		throw std::invalid_argument(
			"Spherical-obstacle boundary ID must be 0 or 1");
	}
	if (!point.allFinite())
	{
		throw std::invalid_argument("Spherical-obstacle query point must be finite");
	}
	const Vector3d center = boundaryId == 0 ? outerCenter : obstacleCenter;
	const double boundaryRadius = boundaryId == 0 ? outerRadius : obstacleRadius;
	const Vector3d offset = point - center;
	const double radius = offset.norm();
	if (!(radius > 0.0))
	{
		throw std::domain_error(
			"Spherical-boundary clearance derivatives are undefined at its center");
	}
	const Vector3d radial = offset / radius;
	const Matrix3d radialHessian =
		(Matrix3d::Identity() - radial * radial.transpose()) / radius;
	DomainQuery result;
	result.boundaryId = boundaryId;
	result.closestPoint = center + boundaryRadius * radial;
	if (boundaryId == 0)
	{
		result.clearance = outerRadius - radius;
		result.clearanceNormal = -radial;
		result.clearanceHessian = -radialHessian;
	}
	else
	{
		result.clearance = radius - obstacleRadius;
		result.clearanceNormal = radial;
		result.clearanceHessian = radialHessian;
	}
	result.admissible = result.clearance >= 0.0;
	return result;
}
