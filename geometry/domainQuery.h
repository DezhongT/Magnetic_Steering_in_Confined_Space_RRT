#ifndef DOMAINQUERY_H
#define DOMAINQUERY_H

#include "eigenIncludes.h"

#include <limits>

struct DomainQuery
{
	double clearance = -std::numeric_limits<double>::infinity();
	Vector3d closestPoint = Vector3d::Zero();
	// Gradient of clearance: this points toward increasing admissibility.
	Vector3d clearanceNormal = Vector3d::Zero();
	// Hessian of clearance with respect to the query point.
	Matrix3d clearanceHessian = Matrix3d::Zero();
	int boundaryId = -1;
	bool admissible = false;
};

inline DomainQuery applyClearanceMargin(DomainQuery query, double margin)
{
	query.clearance -= margin;
	query.admissible = query.clearance >= 0.0;
	return query;
}

#endif
