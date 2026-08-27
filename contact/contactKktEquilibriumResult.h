#ifndef CONTACTKKTEQUILIBRIUMRESULT_H
#define CONTACTKKTEQUILIBRIUMRESULT_H

#include "contact/contactCandidate.h"
#include "rodState.h"

#include <limits>
#include <vector>

struct ContactKktEquilibriumResult
{
	bool success = false;
	bool rolledBack = false;
	bool tipSafe = false;
	RodState state;
	std::vector<ContactCandidate> activeContacts;
	VectorXd multipliers;
	int nonlinearIterations = 0;
	int activeSetUpdates = 0;
	int lineSearchBacktracks = 0;
	double stationarityNorm = std::numeric_limits<double>::infinity();
	double activeConstraintNorm = std::numeric_limits<double>::infinity();
	double complementarityNorm = std::numeric_limits<double>::infinity();
	double minimumBodyGap = -std::numeric_limits<double>::infinity();
	double tipClearance = -std::numeric_limits<double>::infinity();
	double minimumMultiplier = std::numeric_limits<double>::infinity();
};

#endif
