#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
double stateDifference(const RodState &a, const RodState &b)
{
	double difference = 0.0;
	auto include = [&difference](const auto &left, const auto &right) {
		difference = std::max(difference, (left - right).norm());
	};

	include(a.configuration, b.configuration);
	include(a.previousConfiguration, b.previousConfiguration);
	include(a.velocity, b.velocity);
	include(a.referenceDirector1, b.referenceDirector1);
	include(a.referenceDirector2, b.referenceDirector2);
	include(a.previousReferenceDirector1, b.previousReferenceDirector1);
	include(a.previousReferenceDirector2, b.previousReferenceDirector2);
	include(a.materialDirector1, b.materialDirector1);
	include(a.materialDirector2, b.materialDirector2);
	include(a.previousMaterialDirector1, b.previousMaterialDirector1);
	include(a.previousMaterialDirector2, b.previousMaterialDirector2);
	include(a.tangent, b.tangent);
	include(a.previousTangent, b.previousTangent);
	include(a.referenceTwist, b.referenceTwist);
	include(a.previousReferenceTwist, b.previousReferenceTwist);
	include(a.edgeLength, b.edgeLength);
	include(a.curvatureBinormal, b.curvatureBinormal);
	include(a.curvature, b.curvature);
	return difference;
}

bool checkInitialRoundTrip()
{
	setInput input;
	world simulation(input);
	simulation.setRodStepper();

	const double initialResidual = simulation.getStaticResidualNorm();
	const RodState initialState = simulation.captureRodState();

	VectorXd increment = VectorXd::Zero(simulation.numFreeDofs());
	increment[increment.size() - 1] = 1.0e-3;
	simulation.applyFreeDofIncrement(increment);
	const double perturbedResidual = simulation.getStaticResidualNorm();

	if (!std::isfinite(perturbedResidual) || perturbedResidual <= initialResidual + 1.0e-8)
	{
		std::cerr << "A free-end perturbation did not increase the static residual: "
				  << perturbedResidual << '\n';
		return false;
	}

	simulation.restoreRodState(initialState);
	const RodState restoredState = simulation.captureRodState();
	const double restoredResidual = simulation.getStaticResidualNorm();

	if (stateDifference(initialState, restoredState) != 0.0)
	{
		std::cerr << "Initial RodState was not restored exactly.\n";
		return false;
	}

	if (std::abs(restoredResidual - initialResidual) > 1.0e-12)
	{
		std::cerr << "Static residual changed after state restoration.\n";
		return false;
	}

	return true;
}

bool checkPostStepRoundTrip(const char *optionFile)
{
	setInput input;
	if (input.LoadOptions(optionFile) != 0)
	{
		return false;
	}

	world simulation(input);
	simulation.setRodStepper();
	simulation.updateTimeStep();
	simulation.updateTimeStep();

	const RodState checkpoint = simulation.captureRodState();
	simulation.updateTimeStep();

	if (stateDifference(checkpoint, simulation.captureRodState()) == 0.0)
	{
		std::cerr << "Loaded dynamic state did not advance after the checkpoint.\n";
		return false;
	}

	simulation.restoreRodState(checkpoint);
	if (stateDifference(checkpoint, simulation.captureRodState()) != 0.0)
	{
		std::cerr << "Post-step RodState was not restored exactly.\n";
		return false;
	}

	return true;
}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <option-file>\n";
		return 1;
	}

	if (!checkInitialRoundTrip() || !checkPostStepRoundTrip(argv[1]))
	{
		return 1;
	}

	std::cout << "RodState capture/restore checks passed.\n";
	return 0;
}
