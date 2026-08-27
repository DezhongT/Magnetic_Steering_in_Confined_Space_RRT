#ifndef SPHERICALSHELLDOMAIN_H
#define SPHERICALSHELLDOMAIN_H

#include "geometry/confinedDomain.h"

class SphericalShellDomain : public ConfinedDomain
{
public:
	SphericalShellDomain(
		const Vector3d &center,
		double referenceRadius,
		double minusThickness,
		double plusThickness);

	DomainQuery query(const Vector3d &point) const override;
	DomainQuery queryBoundary(
		const Vector3d &point, int boundaryId) const override;

	double innerRadius() const;
	double outerRadius() const;

private:
	Vector3d center;
	double referenceRadius;
	double minusThickness;
	double plusThickness;
};

#endif
