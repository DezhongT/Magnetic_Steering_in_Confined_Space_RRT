#include "geometry/doubleSphericalObstacleDomain.h"

#include <cmath>
#include <stdexcept>

DoubleSphericalObstacleDomain::DoubleSphericalObstacleDomain(
	const Vector3d &m_outerCenter,
	double m_outerRadius,
	const Vector3d &m_firstObstacleCenter,
	double m_firstObstacleRadius,
	const Vector3d &m_secondObstacleCenter,
	double m_secondObstacleRadius)
	: outerCenter(m_outerCenter),
	  outerRadius(m_outerRadius),
	  firstObstacleCenter(m_firstObstacleCenter),
	  firstObstacleRadius(m_firstObstacleRadius),
	  secondObstacleCenter(m_secondObstacleCenter),
	  secondObstacleRadius(m_secondObstacleRadius)
{
	const bool firstInside =
		(firstObstacleCenter - outerCenter).norm() + firstObstacleRadius <
		outerRadius;
	const bool secondInside =
		(secondObstacleCenter - outerCenter).norm() + secondObstacleRadius <
		outerRadius;
	if (!outerCenter.allFinite() || !firstObstacleCenter.allFinite() ||
		!secondObstacleCenter.allFinite() || !std::isfinite(outerRadius) ||
		!std::isfinite(firstObstacleRadius) ||
		!std::isfinite(secondObstacleRadius) || !(outerRadius > 0.0) ||
		!(firstObstacleRadius > 0.0) || !(secondObstacleRadius > 0.0) ||
		!firstInside || !secondInside)
	{
		throw std::invalid_argument(
			"Double spherical obstacles must be finite, positive, and strictly inside their cavity");
	}
}

DomainQuery DoubleSphericalObstacleDomain::query(const Vector3d &point) const
{
	if (!point.allFinite())
	{
		throw std::invalid_argument(
			"Double-spherical-obstacle query point must be finite");
	}
	const double clearances[] = {
		outerRadius - (point - outerCenter).norm(),
		(point - firstObstacleCenter).norm() - firstObstacleRadius,
		(point - secondObstacleCenter).norm() - secondObstacleRadius};
	int nearestBoundary = 0;
	for (int boundary = 1; boundary < 3; ++boundary)
	{
		if (clearances[boundary] < clearances[nearestBoundary])
		{
			nearestBoundary = boundary;
		}
	}
	return queryBoundary(point, nearestBoundary);
}

DomainQuery DoubleSphericalObstacleDomain::queryBoundary(
	const Vector3d &point,
	int boundaryId) const
{
	if (boundaryId < 0 || boundaryId > 2)
	{
		throw std::invalid_argument(
			"Double-spherical-obstacle boundary ID must be 0, 1, or 2");
	}
	if (!point.allFinite())
	{
		throw std::invalid_argument(
			"Double-spherical-obstacle query point must be finite");
	}
	const Vector3d center = boundaryId == 0 ? outerCenter :
		(boundaryId == 1 ? firstObstacleCenter : secondObstacleCenter);
	const double boundaryRadius = boundaryId == 0 ? outerRadius :
		(boundaryId == 1 ? firstObstacleRadius : secondObstacleRadius);
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
		result.clearance = radius - boundaryRadius;
		result.clearanceNormal = radial;
		result.clearanceHessian = radialHessian;
	}
	result.admissible = result.clearance >= 0.0;
	return result;
}
