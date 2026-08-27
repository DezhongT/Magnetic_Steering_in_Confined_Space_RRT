#ifndef SURFACE_H
#define SURFACE_H

#include "geometry/surfaceProjection.h"

class Surface
{
public:
	virtual ~Surface() = default;
	virtual SurfaceProjection project(const Vector3d &point) const = 0;
};

#endif
