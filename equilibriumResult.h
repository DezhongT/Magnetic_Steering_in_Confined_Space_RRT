#ifndef EQUILIBRIUMRESULT_H
#define EQUILIBRIUMRESULT_H

#include "rodState.h"

#include <limits>

struct EquilibriumResult
{
	bool success = false;
	RodState state;
	double initialResidualNorm = std::numeric_limits<double>::quiet_NaN();
	double finalResidualNorm = std::numeric_limits<double>::quiet_NaN();
	int newtonIterations = 0;
	int linearSolverInfo = 0;
	int lineSearchBacktracks = 0;
	int infeasibleTrialRejections = 0;
	double minimumAcceptedStepLength = 1.0;
	bool lineSearchFailed = false;
};

#endif
