#ifndef SPHERICALOBSTACLEDOMAIN_H
#define SPHERICALOBSTACLEDOMAIN_H

#include "geometry/confinedDomain.h"

// Admissible region: inside an outer sphere and outside an eccentric sphere.
class SphericalObstacleDomain : public ConfinedDomain
{
public:
	SphericalObstacleDomain(
		const Vector3d &outerCenter,
		double outerRadius,
		const Vector3d &obstacleCenter,
		double obstacleRadius);

	DomainQuery query(const Vector3d &point) const override;
	DomainQuery queryBoundary(
		const Vector3d &point, int boundaryId) const override;

private:
	Vector3d outerCenter;
	double outerRadius;
	Vector3d obstacleCenter;
	double obstacleRadius;
};

#endif
