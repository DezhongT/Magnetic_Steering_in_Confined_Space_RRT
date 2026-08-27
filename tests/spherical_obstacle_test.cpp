#include "geometry/sphericalObstacleDomain.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool checkDerivatives(
	const SphericalObstacleDomain &domain,
	const Vector3d &point,
	int boundaryId)
{
	constexpr double step = 1.0e-6;
	const DomainQuery query = domain.queryBoundary(point, boundaryId);
	Vector3d gradient;
	Matrix3d hessian;
	for (int component = 0; component < 3; ++component)
	{
		Vector3d offset = Vector3d::Zero();
		offset[component] = step;
		const DomainQuery plus = domain.queryBoundary(point + offset, boundaryId);
		const DomainQuery minus = domain.queryBoundary(point - offset, boundaryId);
		gradient[component] = (plus.clearance - minus.clearance) / (2.0 * step);
		hessian.col(component) =
			(plus.clearanceNormal - minus.clearanceNormal) / (2.0 * step);
	}
	return (gradient - query.clearanceNormal).norm() < 1.0e-9 &&
		(hessian - query.clearanceHessian).norm() < 1.0e-9;
}
}

int main()
{
	const Vector3d outerCenter(0.1, -0.2, 0.3);
	const Vector3d obstacleCenter(0.4, -0.1, 0.2);
	const SphericalObstacleDomain domain(
		outerCenter, 2.0, obstacleCenter, 0.25);
	const Vector3d obstaclePoint = obstacleCenter + Vector3d(0.3, 0.0, 0.0);
	const Vector3d outerPoint = outerCenter + Vector3d(1.9, 0.0, 0.0);
	const DomainQuery nearObstacle = domain.query(obstaclePoint);
	const DomainQuery nearOuter = domain.query(outerPoint);
	if (nearObstacle.boundaryId != 1 ||
		std::abs(nearObstacle.clearance - 0.05) > 1.0e-14 ||
		(nearObstacle.clearanceNormal - Vector3d::UnitX()).norm() > 1.0e-14 ||
		nearOuter.boundaryId != 0 ||
		std::abs(nearOuter.clearance - 0.1) > 1.0e-14 ||
		(nearOuter.clearanceNormal + Vector3d::UnitX()).norm() > 1.0e-14 ||
		!nearObstacle.admissible || !nearOuter.admissible)
	{
		std::cerr << "Spherical-obstacle clearance selection failed.\n";
		return 1;
	}
	const Vector3d derivativePoint(0.9, 0.6, -0.4);
	if (!checkDerivatives(domain, derivativePoint, 0) ||
		!checkDerivatives(domain, derivativePoint, 1))
	{
		std::cerr << "Spherical-obstacle derivative finite difference failed.\n";
		return 1;
	}

	bool rejectedOutside = false;
	try
	{
		SphericalObstacleDomain invalid(
			Vector3d::Zero(), 1.0, Vector3d(0.9, 0.0, 0.0), 0.2);
		(void)invalid;
	}
	catch (const std::invalid_argument &)
	{
		rejectedOutside = true;
	}
	if (!rejectedOutside)
	{
		std::cerr << "Spherical obstacle outside its cavity was accepted.\n";
		return 1;
	}
	std::cout << "Spherical-obstacle clearance, gradient, and Hessian passed.\n";
	return 0;
}
