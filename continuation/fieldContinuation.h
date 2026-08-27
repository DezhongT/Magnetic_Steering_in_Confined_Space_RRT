#ifndef FIELDCONTINUATION_H
#define FIELDCONTINUATION_H

#include "contact/contactConstraintAnalysis.h"
#include "contact/contactKktEquilibriumResult.h"
#include "eigenIncludes.h"

#include <limits>
#include <vector>

struct FieldContinuationOptions
{
	double initialStepFraction = 0.25;
	double minimumStepFraction = 1.0e-3;
	double maximumStepFraction = 0.5;
	double stepReduction = 0.5;
	double stepGrowth = 1.5;
	double stabilityTolerance = 0.0;
	int easyCorrectorIterations = 5;
	int maximumAttempts = 200;
};

struct FieldContinuationPoint
{
	Vector3d field = Vector3d::Zero();
	ContactKktEquilibriumResult equilibrium;
	ContactStabilityAnalysis stability;
	double pathFraction = 0.0;
	double acceptedStepFraction = 0.0;
	double configurationPredictorError = 0.0;
	double multiplierPredictorError = 0.0;
	int contactsAdded = 0;
	int contactsReleased = 0;
	bool contactSetChanged = false;
};

struct FieldContinuationResult
{
	bool success = false;
	bool rolledBack = false;
	Vector3d startField = Vector3d::Zero();
	Vector3d targetField = Vector3d::Zero();
	int attemptedSteps = 0;
	int rejectedSteps = 0;
	double minimumAttemptedStepFraction =
		std::numeric_limits<double>::infinity();
	std::vector<FieldContinuationPoint> points;
};

#endif
