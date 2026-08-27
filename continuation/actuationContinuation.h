#ifndef ACTUATIONCONTINUATION_H
#define ACTUATIONCONTINUATION_H

#include "actuation.h"
#include "continuation/fieldContinuation.h"

#include <limits>
#include <string>
#include <vector>

struct ActuationContinuationPoint
{
	Actuation actuation;
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

struct ActuationContinuationResult
{
	bool success = false;
	bool rolledBack = false;
	Actuation startActuation;
	Actuation targetActuation;
	int attemptedSteps = 0;
	int rejectedSteps = 0;
	double minimumAttemptedStepFraction =
		std::numeric_limits<double>::infinity();
	std::string failureReason;
	std::vector<ActuationContinuationPoint> points;
};

#endif
