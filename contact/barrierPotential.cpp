#include "contact/barrierPotential.h"

#include <cmath>
#include <stdexcept>

BarrierPotentialEvaluation evaluateBarrierPotential(
	double gap,
	double activationDistance,
	double stiffness)
{
	if (!std::isfinite(activationDistance) || activationDistance <= 0.0)
	{
		throw std::invalid_argument(
			"Barrier activation distance must be positive and finite");
	}
	if (!std::isfinite(stiffness) || stiffness < 0.0)
	{
		throw std::invalid_argument(
			"Barrier stiffness must be nonnegative and finite");
	}
	if (!std::isfinite(gap))
	{
		throw std::invalid_argument("Barrier gap must be finite");
	}

	BarrierPotentialEvaluation result;
	result.valid = gap > 0.0;
	if (!result.valid || gap >= activationDistance || stiffness == 0.0)
	{
		return result;
	}

	result.active = true;
	const double offset = gap - activationDistance;
	const double logRatio = std::log(gap / activationDistance);
	result.energy = -stiffness * offset * offset * logRatio;
	result.firstDerivative = -stiffness *
		(2.0 * offset * logRatio + offset * offset / gap);
	result.secondDerivative = -stiffness *
		(2.0 * logRatio + 4.0 * offset / gap -
		 offset * offset / (gap * gap));
	return result;
}
