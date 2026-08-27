#ifndef RODSTATE_H
#define RODSTATE_H

#include "eigenIncludes.h"

// A value snapshot of all mutable DER configuration and frame-history fields.
// The history fields are required when restoring a branch for continuation:
// positions and twist alone do not define the time-parallel-transport state.
struct RodState
{
	VectorXd configuration;
	VectorXd previousConfiguration;
	VectorXd velocity;

	MatrixXd referenceDirector1;
	MatrixXd referenceDirector2;
	MatrixXd previousReferenceDirector1;
	MatrixXd previousReferenceDirector2;

	MatrixXd materialDirector1;
	MatrixXd materialDirector2;
	MatrixXd previousMaterialDirector1;
	MatrixXd previousMaterialDirector2;

	MatrixXd tangent;
	MatrixXd previousTangent;
	VectorXd referenceTwist;
	VectorXd previousReferenceTwist;

	VectorXd edgeLength;
	MatrixXd curvatureBinormal;
	MatrixXd curvature;
};

#endif
