#ifndef PLANARSLABDOMAIN_H
#define PLANARSLABDOMAIN_H

#include "geometry/confinedDomain.h"

class PlanarSlabDomain : public ConfinedDomain
{
public:
	PlanarSlabDomain(
		const Vector3d &planePoint,
		const Vector3d &orientedNormal,
		double minusThickness,
		double plusThickness);

	DomainQuery query(const Vector3d &point) const override;
	DomainQuery queryBoundary(
		const Vector3d &point, int boundaryId) const override;

private:
	Vector3d planePoint;
	Vector3d normal;
	double minusThickness;
	double plusThickness;
};

#endif
