#ifndef PLANNERSTATE_H
#define PLANNERSTATE_H

#include "actuation.h"
#include "contact/contactCandidate.h"
#include "rodState.h"

#include <limits>
#include <vector>

struct PlannerState
{
	RodState rodState;
	Actuation actuation;
	std::vector<ContactCandidate> activeContacts;
	VectorXd multipliers;
	double stabilityMargin = std::numeric_limits<double>::quiet_NaN();
	double stationarityNorm = std::numeric_limits<double>::infinity();
	double activeConstraintNorm = std::numeric_limits<double>::infinity();
	double complementarityNorm = std::numeric_limits<double>::infinity();
	double minimumBodyGap = -std::numeric_limits<double>::infinity();
	double tipClearance = -std::numeric_limits<double>::infinity();
	bool tipSafe = false;
};

#endif
