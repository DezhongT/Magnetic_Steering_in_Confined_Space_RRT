#include "geometry/planarSlabDomain.h"

#include <stdexcept>

PlanarSlabDomain::PlanarSlabDomain(
	const Vector3d &m_planePoint,
	const Vector3d &orientedNormal,
	double m_minusThickness,
	double m_plusThickness)
	: planePoint(m_planePoint),
	  minusThickness(m_minusThickness),
	  plusThickness(m_plusThickness)
{
	if (!(orientedNormal.norm() > 0.0))
	{
		throw std::invalid_argument("Planar slab normal must be nonzero");
	}
	if (minusThickness < 0.0 || plusThickness < 0.0)
	{
		throw std::invalid_argument("Planar slab thicknesses must be nonnegative");
	}
	normal = orientedNormal.normalized();
}

DomainQuery PlanarSlabDomain::query(const Vector3d &point) const
{
	const double signedOffset = (point - planePoint).dot(normal);
	const double lowerClearance = signedOffset + minusThickness;
	const double upperClearance = plusThickness - signedOffset;
	if (lowerClearance <= upperClearance)
	{
		return queryBoundary(point, 0);
	}
	return queryBoundary(point, 1);
}

DomainQuery PlanarSlabDomain::queryBoundary(
	const Vector3d &point,
	int boundaryId) const
{
	if (boundaryId != 0 && boundaryId != 1)
	{
		throw std::invalid_argument("Planar slab boundary ID must be 0 or 1");
	}
	const double signedOffset = (point - planePoint).dot(normal);
	DomainQuery result;
	result.boundaryId = boundaryId;
	if (boundaryId == 0)
	{
		result.clearance = signedOffset + minusThickness;
		result.closestPoint = planePoint - minusThickness * normal +
			((point - planePoint) - signedOffset * normal);
		result.clearanceNormal = normal;
	}
	else
	{
		result.clearance = plusThickness - signedOffset;
		result.closestPoint = planePoint + plusThickness * normal +
			((point - planePoint) - signedOffset * normal);
		result.clearanceNormal = -normal;
	}
	result.admissible = result.clearance >= 0.0;
	return result;
}
