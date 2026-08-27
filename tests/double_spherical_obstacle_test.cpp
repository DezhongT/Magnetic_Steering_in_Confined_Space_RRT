#include "geometry/doubleSphericalObstacleDomain.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool checkDerivatives(
	const DoubleSphericalObstacleDomain &domain,
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
		const DomainQuery plus = domain.queryBoundary(
			point + offset, boundaryId);
		const DomainQuery minus = domain.queryBoundary(
			point - offset, boundaryId);
		gradient[component] =
			(plus.clearance - minus.clearance) / (2.0 * step);
		hessian.col(component) =
			(plus.clearanceNormal - minus.clearanceNormal) / (2.0 * step);
	}
	return (gradient - query.clearanceNormal).norm() < 1.0e-9 &&
		(hessian - query.clearanceHessian).norm() < 1.0e-9;
}
}

int main()
{
	const DoubleSphericalObstacleDomain domain(
		Vector3d::Zero(), 2.0,
		Vector3d(0.4, 0.2, 0.0), 0.1,
		Vector3d(0.4, -0.2, 0.0), 0.15);
	const DomainQuery first = domain.query(Vector3d(0.4, 0.31, 0.0));
	const DomainQuery second = domain.query(Vector3d(0.4, -0.36, 0.0));
	const DomainQuery outer = domain.query(Vector3d(1.9, 0.0, 0.0));
	if (first.boundaryId != 1 || std::abs(first.clearance - 0.01) > 1.0e-14 ||
		second.boundaryId != 2 || std::abs(second.clearance - 0.01) > 1.0e-14 ||
		outer.boundaryId != 0 || std::abs(outer.clearance - 0.1) > 1.0e-14 ||
		!checkDerivatives(domain, Vector3d(0.8, 0.5, -0.3), 0) ||
		!checkDerivatives(domain, Vector3d(0.8, 0.5, -0.3), 1) ||
		!checkDerivatives(domain, Vector3d(0.8, -0.5, 0.3), 2))
	{
		std::cerr << "Double-spherical-obstacle query validation failed.\n";
		return 1;
	}

	bool rejected = false;
	try
	{
		DoubleSphericalObstacleDomain invalid(
			Vector3d::Zero(), 1.0,
			Vector3d(0.2, 0.0, 0.0), 0.1,
			Vector3d(0.95, 0.0, 0.0), 0.1);
		(void)invalid;
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	if (!rejected)
	{
		std::cerr << "Obstacle outside the cavity was accepted.\n";
		return 1;
	}
	std::cout << "Double-spherical-obstacle clearance derivatives passed.\n";
	return 0;
}
