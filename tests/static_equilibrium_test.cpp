#include "setInput.h"
#include "world.h"

#include <cmath>
#include <iostream>

int main()
{
	setInput input;
	world simulation(input);
	simulation.setRodStepper();

	const double residualNorm = simulation.getStaticResidualNorm();
	constexpr double tolerance = 1.0e-10;

	if (!std::isfinite(residualNorm))
	{
		std::cerr << "Static residual is not finite.\n";
		return 1;
	}

	if (residualNorm > tolerance)
	{
		std::cerr << "Expected the unloaded straight rod to be in equilibrium, "
				  << "but the residual norm is " << residualNorm << ".\n";
		return 1;
	}

	std::cout << "Zero-load static residual norm: " << residualNorm << '\n';
	return 0;
}
