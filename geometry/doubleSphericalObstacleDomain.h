#ifndef DOUBLESPHERICALOBSTACLEDOMAIN_H
#define DOUBLESPHERICALOBSTACLEDOMAIN_H

#include "geometry/confinedDomain.h"

// Admissible region: inside an outer sphere and outside two excluded spheres.
class DoubleSphericalObstacleDomain : public ConfinedDomain
{
public:
	DoubleSphericalObstacleDomain(
		const Vector3d &outerCenter,
		double outerRadius,
		const Vector3d &firstObstacleCenter,
		double firstObstacleRadius,
		const Vector3d &secondObstacleCenter,
		double secondObstacleRadius);

	DomainQuery query(const Vector3d &point) const override;
	DomainQuery queryBoundary(
		const Vector3d &point, int boundaryId) const override;

private:
	Vector3d outerCenter;
	double outerRadius;
	Vector3d firstObstacleCenter;
	double firstObstacleRadius;
	Vector3d secondObstacleCenter;
	double secondObstacleRadius;
};

#endif
