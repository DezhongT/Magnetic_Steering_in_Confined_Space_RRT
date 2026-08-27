#include "geometry/triangleMeshSurface.h"

#include "geometry/triangleGeometry.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <utility>

TriangleMeshSurface::TriangleMeshSurface(
	std::vector<Vector3d> m_vertices,
	std::vector<Vector3i> m_triangles)
	: vertices(std::move(m_vertices)), triangles(std::move(m_triangles))
{
	if (vertices.empty() || triangles.empty())
	{
		throw std::invalid_argument("Triangle mesh must contain vertices and triangles");
	}
	for (const Vector3i &triangle : triangles)
	{
		for (int corner = 0; corner < 3; ++corner)
		{
			if (triangle[corner] < 0 || triangle[corner] >= static_cast<int>(vertices.size()))
			{
				throw std::out_of_range("Triangle mesh index is out of range");
			}
		}
		const Vector3d edge01 = vertices[triangle[1]] - vertices[triangle[0]];
		const Vector3d edge02 = vertices[triangle[2]] - vertices[triangle[0]];
		const Vector3d edge12 = vertices[triangle[2]] - vertices[triangle[1]];
		const double maximumSquaredEdge = std::max(
			{edge01.squaredNorm(), edge02.squaredNorm(), edge12.squaredNorm()});
		const double squaredAreaVector = edge01.cross(edge02).squaredNorm();
		if (maximumSquaredEdge == 0.0 ||
			squaredAreaVector <= 1.0e-24 * maximumSquaredEdge * maximumSquaredEdge)
		{
			throw std::invalid_argument("Triangle mesh contains a degenerate triangle");
		}
	}
}

TriangleMeshSurface TriangleMeshSurface::fromFiles(
	const std::string &vertexFile,
	const std::string &triangleFile)
{
	std::ifstream vertexInput(vertexFile);
	if (!vertexInput)
	{
		throw std::runtime_error("Could not open mesh vertex file: " + vertexFile);
	}
	std::vector<Vector3d> loadedVertices;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	while (vertexInput >> x >> y >> z)
	{
		loadedVertices.emplace_back(x, y, z);
	}
	if (!vertexInput.eof())
	{
		throw std::runtime_error("Invalid data in mesh vertex file: " + vertexFile);
	}

	std::ifstream triangleInput(triangleFile);
	if (!triangleInput)
	{
		throw std::runtime_error("Could not open mesh triangle file: " + triangleFile);
	}
	std::vector<Vector3i> loadedTriangles;
	int a = 0;
	int b = 0;
	int c = 0;
	while (triangleInput >> a >> b >> c)
	{
		loadedTriangles.emplace_back(a, b, c);
	}
	if (!triangleInput.eof())
	{
		throw std::runtime_error("Invalid data in mesh triangle file: " + triangleFile);
	}
	return TriangleMeshSurface(
		std::move(loadedVertices), std::move(loadedTriangles));
}

SurfaceProjection TriangleMeshSurface::project(const Vector3d &point) const
{
	SurfaceProjection best;
	double bestSquaredDistance = std::numeric_limits<double>::infinity();
	for (int triangleId = 0; triangleId < static_cast<int>(triangles.size()); ++triangleId)
	{
		const Vector3i &indices = triangles[triangleId];
		const Vector3d &a = vertices[indices[0]];
		const Vector3d &b = vertices[indices[1]];
		const Vector3d &c = vertices[indices[2]];
		const TriangleProjection projection =
			projectPointToTriangle(point, a, b, c);
		if (projection.squaredDistance < bestSquaredDistance)
		{
			bestSquaredDistance = projection.squaredDistance;
			best.closestPoint = projection.closestPoint;
			best.surfaceNormal = (b - a).cross(c - a).normalized();
			best.primitiveId = triangleId;
		}
	}
	best.distance = std::sqrt(bestSquaredDistance);
	best.signedOffset = (point - best.closestPoint).dot(best.surfaceNormal);
	return best;
}

int TriangleMeshSurface::vertexCount() const
{
	return static_cast<int>(vertices.size());
}

int TriangleMeshSurface::triangleCount() const
{
	return static_cast<int>(triangles.size());
}
