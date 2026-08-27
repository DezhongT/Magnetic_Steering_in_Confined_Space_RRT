#include "contact/contactConstraintAnalysis.h"

#include <cmath>
#include <iostream>

int main()
{
	MatrixXd hessian = MatrixXd::Zero(3, 3);
	hessian.diagonal() << 2.0, -1.0, 4.0;
	MatrixXd constrainNegative(1, 3);
	constrainNegative << 0.0, 1.0, 0.0;
	const ContactStabilityAnalysis stable = analyzeContactStability(
		hessian, constrainNegative, 1.0e-12);
	if (!stable.valid || !stable.stable || stable.constraintRank != 1 ||
		stable.nullSpaceBasis.cols() != 2 ||
		stable.nullSpaceResidualNorm > 1.0e-12 ||
		std::abs(stable.minimumEigenvalue - 2.0) > 1.0e-12)
	{
		std::cerr << "Manufactured reduced-Hessian stability check failed.\n";
		return 1;
	}

	MatrixXd constrainPositive(1, 3);
	constrainPositive << 1.0, 0.0, 0.0;
	const ContactStabilityAnalysis unstable = analyzeContactStability(
		hessian, constrainPositive, 1.0e-12);
	if (!unstable.valid || unstable.stable ||
		std::abs(unstable.minimumEigenvalue + 1.0) > 1.0e-12)
	{
		std::cerr << "Manufactured unstable mode was not detected.\n";
		return 1;
	}

	MatrixXd sensitivityHessian = MatrixXd::Zero(2, 2);
	sensitivityHessian.diagonal() << 2.0, 3.0;
	MatrixXd constraint(1, 2);
	constraint << 1.0, 0.0;
	MatrixXd residualDerivative(2, 2);
	residualDerivative << 5.0, -2.0, 6.0, 3.0;
	const ConstrainedSensitivityResult sensitivity =
		solveConstrainedSensitivity(
			sensitivityHessian, constraint, residualDerivative);
	MatrixXd expectedConfiguration(2, 2);
	expectedConfiguration << 0.0, 0.0, -2.0, -1.0;
	MatrixXd expectedMultiplier(1, 2);
	expectedMultiplier << 5.0, -2.0;
	if (!sensitivity.success || sensitivity.linearResidualNorm > 1.0e-12 ||
		(sensitivity.configurationDerivative - expectedConfiguration).norm() > 1.0e-12 ||
		(sensitivity.multiplierDerivative - expectedMultiplier).norm() > 1.0e-12 ||
		(constraint * sensitivity.configurationDerivative).norm() > 1.0e-12)
	{
		std::cerr << "Manufactured constrained sensitivity check failed.\n";
		return 1;
	}

	std::cout << "Constraint analysis: stable_minimum="
			  << stable.minimumEigenvalue
			  << ", unstable_minimum=" << unstable.minimumEigenvalue
			  << ", sensitivity_residual="
			  << sensitivity.linearResidualNorm << '\n';
	return 0;
}
