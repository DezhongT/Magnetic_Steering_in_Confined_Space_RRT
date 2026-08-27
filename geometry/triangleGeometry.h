#ifndef TRIANGLEGEOMETRY_H
#define TRIANGLEGEOMETRY_H

#include "eigenIncludes.h"

struct TriangleProjection
{
	Vector3d closestPoint = Vector3d::Zero();
	Vector3d barycentric = Vector3d::Zero();
	double squaredDistance = 0.0;
	bool degenerate = false;
};

TriangleProjection projectPointToTriangle(
	const Vector3d &point,
	const Vector3d &a,
	const Vector3d &b,
	const Vector3d &c);

#endif
