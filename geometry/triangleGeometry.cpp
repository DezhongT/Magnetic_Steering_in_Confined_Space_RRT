#include "geometry/triangleGeometry.h"

#include <algorithm>
#include <array>

namespace
{
struct SegmentProjection
{
	Vector3d point;
	double parameter;
	double squaredDistance;
};

SegmentProjection projectPointToSegment(
	const Vector3d &point,
	const Vector3d &start,
	const Vector3d &end)
{
	const Vector3d edge = end - start;
	const double squaredLength = edge.squaredNorm();
	double parameter = 0.0;
	if (squaredLength > 0.0)
	{
		parameter = std::clamp((point - start).dot(edge) / squaredLength, 0.0, 1.0);
	}
	const Vector3d closest = start + parameter * edge;
	return {closest, parameter, (point - closest).squaredNorm()};
}
}

TriangleProjection projectPointToTriangle(
	const Vector3d &point,
	const Vector3d &a,
	const Vector3d &b,
	const Vector3d &c)
{
	const Vector3d ab = b - a;
	const Vector3d ac = c - a;
	const Vector3d bc = c - b;
	const double maximumSquaredEdge =
		std::max({ab.squaredNorm(), ac.squaredNorm(), bc.squaredNorm()});
	const double squaredAreaVector = ab.cross(ac).squaredNorm();

	if (maximumSquaredEdge == 0.0 ||
		squaredAreaVector <= 1.0e-24 * maximumSquaredEdge * maximumSquaredEdge)
	{
		const SegmentProjection abProjection = projectPointToSegment(point, a, b);
		const SegmentProjection bcProjection = projectPointToSegment(point, b, c);
		const SegmentProjection caProjection = projectPointToSegment(point, c, a);
		TriangleProjection result;
		result.degenerate = true;
		if (abProjection.squaredDistance <= bcProjection.squaredDistance &&
			abProjection.squaredDistance <= caProjection.squaredDistance)
		{
			result.closestPoint = abProjection.point;
			result.barycentric << 1.0 - abProjection.parameter, abProjection.parameter, 0.0;
			result.squaredDistance = abProjection.squaredDistance;
		}
		else if (bcProjection.squaredDistance <= caProjection.squaredDistance)
		{
			result.closestPoint = bcProjection.point;
			result.barycentric << 0.0, 1.0 - bcProjection.parameter, bcProjection.parameter;
			result.squaredDistance = bcProjection.squaredDistance;
		}
		else
		{
			result.closestPoint = caProjection.point;
			result.barycentric << caProjection.parameter, 0.0, 1.0 - caProjection.parameter;
			result.squaredDistance = caProjection.squaredDistance;
		}
		return result;
	}

	const Vector3d ap = point - a;
	const double d1 = ab.dot(ap);
	const double d2 = ac.dot(ap);
	if (d1 <= 0.0 && d2 <= 0.0)
	{
		return {a, Vector3d(1.0, 0.0, 0.0), (point - a).squaredNorm(), false};
	}

	const Vector3d bp = point - b;
	const double d3 = ab.dot(bp);
	const double d4 = ac.dot(bp);
	if (d3 >= 0.0 && d4 <= d3)
	{
		return {b, Vector3d(0.0, 1.0, 0.0), (point - b).squaredNorm(), false};
	}

	const double vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
	{
		const double v = d1 / (d1 - d3);
		const Vector3d closest = a + v * ab;
		return {closest, Vector3d(1.0 - v, v, 0.0),
			(point - closest).squaredNorm(), false};
	}

	const Vector3d cp = point - c;
	const double d5 = ab.dot(cp);
	const double d6 = ac.dot(cp);
	if (d6 >= 0.0 && d5 <= d6)
	{
		return {c, Vector3d(0.0, 0.0, 1.0), (point - c).squaredNorm(), false};
	}

	const double vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
	{
		const double w = d2 / (d2 - d6);
		const Vector3d closest = a + w * ac;
		return {closest, Vector3d(1.0 - w, 0.0, w),
			(point - closest).squaredNorm(), false};
	}

	const double va = d3 * d6 - d5 * d4;
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
	{
		const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		const Vector3d closest = b + w * bc;
		return {closest, Vector3d(0.0, 1.0 - w, w),
			(point - closest).squaredNorm(), false};
	}

	const double denominator = 1.0 / (va + vb + vc);
	const double v = vb * denominator;
	const double w = vc * denominator;
	const Vector3d closest = a + ab * v + ac * w;
	return {closest, Vector3d(1.0 - v - w, v, w),
		(point - closest).squaredNorm(), false};
}
