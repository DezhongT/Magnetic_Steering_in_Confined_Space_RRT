#ifndef TRIANGLEMESHSURFACE_H
#define TRIANGLEMESHSURFACE_H

#include "geometry/surface.h"

#include <string>
#include <vector>

class TriangleMeshSurface : public Surface
{
public:
	TriangleMeshSurface(
		std::vector<Vector3d> vertices,
		std::vector<Vector3i> triangles);

	static TriangleMeshSurface fromFiles(
		const std::string &vertexFile,
		const std::string &triangleFile);

	SurfaceProjection project(const Vector3d &point) const override;
	int vertexCount() const;
	int triangleCount() const;

private:
	std::vector<Vector3d> vertices;
	std::vector<Vector3i> triangles;
};

#endif
