#ifndef CONTACTCONSTRAINTANALYSIS_H
#define CONTACTCONSTRAINTANALYSIS_H

#include "eigenIncludes.h"
#include "contact/contactCandidate.h"

#include <limits>
#include <vector>

struct ContactStabilityAnalysis
{
	bool valid = false;
	bool stable = false;
	int constraintRank = 0;
	MatrixXd nullSpaceBasis;
	MatrixXd reducedHessian;
	VectorXd eigenvalues;
	MatrixXd eigenvectors;
	double minimumEigenvalue = std::numeric_limits<double>::infinity();
	double nullSpaceResidualNorm = std::numeric_limits<double>::infinity();
	double hessianSymmetryError = std::numeric_limits<double>::infinity();
};

struct ConstrainedSensitivityResult
{
	bool success = false;
	MatrixXd configurationDerivative;
	MatrixXd multiplierDerivative;
	double linearResidualNorm = std::numeric_limits<double>::infinity();
};

struct ContactEquilibriumSensitivity
{
	bool success = false;
	std::vector<ContactCandidate> activeContacts;
	MatrixXd configurationDerivative;
	MatrixXd multiplierDerivative;
	double linearResidualNorm = std::numeric_limits<double>::infinity();
};

struct ActuationEquilibriumSensitivity
{
	bool success = false;
	std::vector<ContactCandidate> activeContacts;
	// Columns are ordered as (xi, Bx, By, Bz).
	MatrixXd configurationDerivative;
	MatrixXd multiplierDerivative;
	double linearResidualNorm = std::numeric_limits<double>::infinity();
};

ContactStabilityAnalysis analyzeContactStability(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian,
	double eigenvalueTolerance = 0.0);

ConstrainedSensitivityResult solveConstrainedSensitivity(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian,
	const MatrixXd &residualParameterDerivative);

#endif
