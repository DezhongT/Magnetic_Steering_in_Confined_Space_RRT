#ifndef CONTACTKKTSYSTEM_H
#define CONTACTKKTSYSTEM_H

#include "contact/contactCandidate.h"
#include "eigenIncludes.h"

#include <limits>
#include <vector>

// Sign convention for a positive compressive multiplier:
//     stationarity = residual - G^T lambda = 0,
//     gap >= 0, lambda >= 0, gap_i lambda_i = 0.
// The symmetric KKT system solves for dmu = -dlambda.
struct ContactKktLinearization
{
	MatrixXd hessian;
	VectorXd residual;
	MatrixXd constraintJacobian;
	VectorXd gaps;
	VectorXd multipliers;
};

struct ContactKktSystem
{
	MatrixXd matrix;
	VectorXd rightHandSide;
	int primalDofs = 0;
	int contactCount = 0;
};

struct ContactKktSolution
{
	bool success = false;
	VectorXd configurationIncrement;
	VectorXd multiplierIncrement;
	VectorXd updatedMultipliers;
	double linearResidualNorm = std::numeric_limits<double>::infinity();
};

struct PlanarContactKktSeed
{
	ContactKktLinearization linearization;
	std::vector<ContactCandidate> contacts;
};

struct ContactActiveSetUpdate
{
	std::vector<ContactCandidate> contacts;
	VectorXd multipliers;
	int added = 0;
	int released = 0;

	bool changed() const
	{
		return added > 0 || released > 0;
	}
};

ContactKktSystem buildContactKktSystem(
	const ContactKktLinearization &linearization);

ContactKktSolution solveContactKktSystem(
	const ContactKktLinearization &linearization);

ContactActiveSetUpdate updateContactActiveSet(
	const std::vector<ContactCandidate> &activeContacts,
	const VectorXd &multipliers,
	const std::vector<ContactCandidate> &detectedCandidates,
	double gapTolerance,
	double multiplierTolerance);

#endif
