#ifndef SURFACEPROJECTION_H
#define SURFACEPROJECTION_H

#include "eigenIncludes.h"

#include <limits>

struct SurfaceProjection
{
	double distance = std::numeric_limits<double>::infinity();
	double signedOffset = std::numeric_limits<double>::quiet_NaN();
	Vector3d closestPoint = Vector3d::Zero();
	Vector3d surfaceNormal = Vector3d::Zero();
	int primitiveId = -1;
};

#endif
