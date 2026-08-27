#ifndef BARRIERPOTENTIAL_H
#define BARRIERPOTENTIAL_H

#include <stdexcept>

class InfeasibleContactError : public std::runtime_error
{
public:
	explicit InfeasibleContactError(const char *message)
		: std::runtime_error(message)
	{
	}
};

struct BarrierPotentialEvaluation
{
	bool valid = false;
	bool active = false;
	double energy = 0.0;
	double firstDerivative = 0.0;
	double secondDerivative = 0.0;
};

// IPC-style C2 barrier on the open feasible interval gap > 0.
// The energy and its first two derivatives vanish at activationDistance.
BarrierPotentialEvaluation evaluateBarrierPotential(
	double gap,
	double activationDistance,
	double stiffness);

#endif
