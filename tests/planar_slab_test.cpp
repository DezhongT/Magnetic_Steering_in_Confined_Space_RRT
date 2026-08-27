#include "geometry/planarSlabDomain.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
bool checkDirectionalDerivative(
	const PlanarSlabDomain &domain,
	const Vector3d &point,
	const Vector3d &direction)
{
	constexpr double step = 1.0e-7;
	const DomainQuery query = domain.query(point);
	const double finiteDifference =
		(domain.query(point + step * direction).clearance -
		 domain.query(point - step * direction).clearance) /
		(2.0 * step);
	return std::abs(finiteDifference - query.clearanceNormal.dot(direction)) < 1.0e-9;
}
}

int main()
{
	const PlanarSlabDomain domain(
		Vector3d::Zero(), Vector3d(0.0, 0.0, 2.0), 0.2, 0.3);
	const DomainQuery center = domain.query(Vector3d(0.7, -0.4, 0.0));
	const DomainQuery nearUpper = domain.query(Vector3d(0.1, 0.2, 0.25));
	const DomainQuery below = domain.query(Vector3d(0.0, 0.0, -0.3));
	const DomainQuery upperBoundary = domain.query(Vector3d(0.0, 0.0, 0.3));
	const DomainQuery forcedUpper =
		domain.queryBoundary(Vector3d(0.7, -0.4, 0.0), 1);
	if (std::abs(center.clearance - 0.2) > 1.0e-12 || center.boundaryId != 0 ||
		(center.clearanceNormal - Vector3d::UnitZ()).norm() > 1.0e-12 ||
		(center.closestPoint - Vector3d(0.7, -0.4, -0.2)).norm() > 1.0e-12 ||
		std::abs(nearUpper.clearance - 0.05) > 1.0e-12 || nearUpper.boundaryId != 1 ||
		(nearUpper.clearanceNormal + Vector3d::UnitZ()).norm() > 1.0e-12 ||
		below.admissible || std::abs(below.clearance + 0.1) > 1.0e-12 ||
		!upperBoundary.admissible || std::abs(upperBoundary.clearance) > 1.0e-12 ||
		std::abs(forcedUpper.clearance - 0.3) > 1.0e-12 ||
		(forcedUpper.clearanceNormal + Vector3d::UnitZ()).norm() > 1.0e-12)
	{
		std::cerr << "Planar slab clearance/sign convention failed.\n";
		return 1;
	}

	const Vector3d direction = Vector3d(0.3, -0.2, 0.7).normalized();
	if (!checkDirectionalDerivative(domain, Vector3d(0.1, 0.2, -0.1), direction) ||
		!checkDirectionalDerivative(domain, Vector3d(-0.2, 0.1, 0.2), direction))
	{
		std::cerr << "Planar slab clearance normal failed finite differences.\n";
		return 1;
	}

	const DomainQuery bodyQuery = applyClearanceMargin(center, 0.05);
	if (std::abs(bodyQuery.clearance - 0.15) > 1.0e-12 || !bodyQuery.admissible)
	{
		std::cerr << "Body-radius clearance margin failed.\n";
		return 1;
	}

	bool rejectedInvalidNormal = false;
	try
	{
		PlanarSlabDomain invalid(Vector3d::Zero(), Vector3d::Zero(), 0.1, 0.1);
		(void)invalid;
	}
	catch (const std::invalid_argument &)
	{
		rejectedInvalidNormal = true;
	}
	if (!rejectedInvalidNormal)
	{
		std::cerr << "Planar slab accepted a zero normal.\n";
		return 1;
	}

	std::cout << "Planar slab clearance and normal checks passed.\n";
	return 0;
}
