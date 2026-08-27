#ifndef CONTACTCANDIDATE_H
#define CONTACTCANDIDATE_H

#include "eigenIncludes.h"

#include <limits>
#include <vector>

struct ContactCandidate
{
	int rodVertex = -1;
	int boundaryId = -1;
	double centerlineClearance = std::numeric_limits<double>::infinity();
	double gap = std::numeric_limits<double>::infinity();
	Vector3d rodCenter = Vector3d::Zero();
	Vector3d closestBoundaryPoint = Vector3d::Zero();
	// Gradient of gap with respect to the rod-center position.
	Vector3d normal = Vector3d::Zero();
	// Hessian of gap with respect to the rod-center position.
	Matrix3d gapHessian = Matrix3d::Zero();
};

struct TipClearanceResult
{
	double clearance = -std::numeric_limits<double>::infinity();
	double requiredClearance = 0.0;
	Vector3d closestBoundaryPoint = Vector3d::Zero();
	Vector3d normal = Vector3d::Zero();
	int boundaryId = -1;
	bool safe = false;
};

struct ContactDetectionResult
{
	std::vector<ContactCandidate> bodyCandidates;
	TipClearanceResult tip;
	double minimumBodyGap = std::numeric_limits<double>::infinity();
};

#endif
