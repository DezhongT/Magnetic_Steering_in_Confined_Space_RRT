#include "geometry/triangleGeometry.h"
#include "geometry/triangleMeshSurface.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
bool near(const Vector3d &left, const Vector3d &right, double tolerance = 1.0e-12)
{
	return (left - right).norm() <= tolerance;
}
}

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " <vertex-file> <triangle-file>\n";
		return 1;
	}

	const Vector3d a(0.0, 0.0, 0.0);
	const Vector3d b(1.0, 0.0, 0.0);
	const Vector3d c(0.0, 1.0, 0.0);
	const TriangleProjection interior =
		projectPointToTriangle(Vector3d(0.2, 0.3, 1.0), a, b, c);
	const TriangleProjection edge =
		projectPointToTriangle(Vector3d(0.4, -0.2, 0.5), a, b, c);
	const TriangleProjection vertex =
		projectPointToTriangle(Vector3d(1.2, -0.1, 0.2), a, b, c);
	if (!near(interior.closestPoint, Vector3d(0.2, 0.3, 0.0)) ||
		!near(interior.barycentric, Vector3d(0.5, 0.2, 0.3)) ||
		!near(edge.closestPoint, Vector3d(0.4, 0.0, 0.0)) ||
		!near(vertex.closestPoint, b))
	{
		std::cerr << "Point-to-triangle region projection failed.\n";
		return 1;
	}

	const TriangleProjection degenerate = projectPointToTriangle(
		Vector3d(0.6, 1.0, 0.0), a, b, Vector3d(2.0, 0.0, 0.0));
	if (!degenerate.degenerate ||
		!near(degenerate.closestPoint, Vector3d(0.6, 0.0, 0.0)))
	{
		std::cerr << "Degenerate triangle fallback failed.\n";
		return 1;
	}

	std::vector<Vector3d> vertices = {
		Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0),
		Vector3d(1.0, 1.0, 0.0), Vector3d(0.0, 1.0, 0.0)};
	std::vector<Vector3i> triangles = {
		Vector3i(0, 1, 2), Vector3i(0, 2, 3)};
	const TriangleMeshSurface surface(vertices, triangles);
	const SurfaceProjection above = surface.project(Vector3d(0.25, 0.25, 0.4));
	const SurfaceProjection below = surface.project(Vector3d(0.75, 0.25, -0.3));
	const SurfaceProjection outside = surface.project(Vector3d(1.5, 0.5, 0.2));
	if (!near(above.closestPoint, Vector3d(0.25, 0.25, 0.0)) ||
		std::abs(above.distance - 0.4) > 1.0e-12 || above.signedOffset <= 0.0 ||
		std::abs(below.distance - 0.3) > 1.0e-12 || below.signedOffset >= 0.0 ||
		!near(outside.closestPoint, Vector3d(1.0, 0.5, 0.0)) ||
		std::abs(outside.distance - std::sqrt(0.29)) > 1.0e-12)
	{
		std::cerr << "Triangle mesh projection or orientation failed.\n";
		return 1;
	}

	bool rejectedDegenerateMesh = false;
	try
	{
		TriangleMeshSurface invalid(
			{a, b, Vector3d(2.0, 0.0, 0.0)}, {Vector3i(0, 1, 2)});
		(void)invalid;
	}
	catch (const std::invalid_argument &)
	{
		rejectedDegenerateMesh = true;
	}
	if (!rejectedDegenerateMesh)
	{
		std::cerr << "TriangleMeshSurface accepted a degenerate primitive.\n";
		return 1;
	}

	const TriangleMeshSurface fileSurface =
		TriangleMeshSurface::fromFiles(argv[1], argv[2]);
	if (fileSurface.vertexCount() != 177 || fileSurface.triangleCount() != 312 ||
		!std::isfinite(fileSurface.project(Vector3d::Zero()).distance))
	{
		std::cerr << "Repository mesh loading/projection failed.\n";
		return 1;
	}

	std::cout << "Triangle geometry and mesh projection checks passed.\n";
	return 0;
}
