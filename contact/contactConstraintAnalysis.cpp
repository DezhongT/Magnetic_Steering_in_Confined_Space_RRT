#include "contact/contactConstraintAnalysis.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
void validateMechanicsMatrices(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian)
{
	if (hessian.rows() <= 0 || hessian.rows() != hessian.cols() ||
		constraintJacobian.cols() != hessian.cols())
	{
		throw std::invalid_argument(
			"Constrained mechanics matrix dimensions are inconsistent");
	}
	if (!hessian.allFinite() || !constraintJacobian.allFinite())
	{
		throw std::invalid_argument(
			"Constrained mechanics matrices contain nonfinite values");
	}
}

MatrixXd symmetricKktMatrix(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian)
{
	const int primalDofs = static_cast<int>(hessian.rows());
	const int contacts = static_cast<int>(constraintJacobian.rows());
	MatrixXd matrix = MatrixXd::Zero(
		primalDofs + contacts, primalDofs + contacts);
	matrix.topLeftCorner(primalDofs, primalDofs) = hessian;
	if (contacts > 0)
	{
		matrix.topRightCorner(primalDofs, contacts) =
			constraintJacobian.transpose();
		matrix.bottomLeftCorner(contacts, primalDofs) = constraintJacobian;
	}
	return matrix;
}
}

ContactStabilityAnalysis analyzeContactStability(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian,
	double eigenvalueTolerance)
{
	validateMechanicsMatrices(hessian, constraintJacobian);
	if (!std::isfinite(eigenvalueTolerance) || eigenvalueTolerance < 0.0)
	{
		throw std::invalid_argument(
			"Stability eigenvalue tolerance must be finite and nonnegative");
	}

	ContactStabilityAnalysis result;
	const int primalDofs = static_cast<int>(hessian.rows());
	if (constraintJacobian.rows() == 0)
	{
		result.constraintRank = 0;
		result.nullSpaceBasis = MatrixXd::Identity(primalDofs, primalDofs);
	}
	else
	{
		Eigen::JacobiSVD<MatrixXd> svd(
			constraintJacobian, Eigen::ComputeFullV);
		const double largestSingularValue = svd.singularValues().size() == 0 ?
			0.0 : svd.singularValues()[0];
		const double rankTolerance = std::max(
			constraintJacobian.rows(), constraintJacobian.cols()) *
			std::numeric_limits<double>::epsilon() * largestSingularValue;
		result.constraintRank = 0;
		for (int index = 0; index < svd.singularValues().size(); ++index)
		{
			if (svd.singularValues()[index] > rankTolerance)
			{
				++result.constraintRank;
			}
		}
		result.nullSpaceBasis = svd.matrixV().rightCols(
			primalDofs - result.constraintRank);
	}

	result.nullSpaceResidualNorm =
		(constraintJacobian * result.nullSpaceBasis).norm();
	const double hessianScale = std::max(1.0, hessian.norm());
	result.hessianSymmetryError =
		(hessian - hessian.transpose()).norm() / hessianScale;
	const MatrixXd symmetricHessian = 0.5 * (hessian + hessian.transpose());
	result.reducedHessian = result.nullSpaceBasis.transpose() *
		symmetricHessian * result.nullSpaceBasis;

	if (result.reducedHessian.rows() == 0)
	{
		result.eigenvalues = VectorXd(0);
		result.eigenvectors = MatrixXd(0, 0);
		result.minimumEigenvalue = std::numeric_limits<double>::infinity();
		result.stable = true;
		result.valid = true;
		return result;
	}

	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(result.reducedHessian);
	if (eigenSolver.info() != Eigen::Success)
	{
		return result;
	}
	result.eigenvalues = eigenSolver.eigenvalues();
	result.eigenvectors = eigenSolver.eigenvectors();
	result.minimumEigenvalue = result.eigenvalues.minCoeff();
	result.stable = result.minimumEigenvalue > eigenvalueTolerance;
	result.valid = std::isfinite(result.minimumEigenvalue) &&
		std::isfinite(result.nullSpaceResidualNorm) &&
		std::isfinite(result.hessianSymmetryError);
	return result;
}

ConstrainedSensitivityResult solveConstrainedSensitivity(
	const MatrixXd &hessian,
	const MatrixXd &constraintJacobian,
	const MatrixXd &residualParameterDerivative)
{
	validateMechanicsMatrices(hessian, constraintJacobian);
	if (residualParameterDerivative.rows() != hessian.rows() ||
		!residualParameterDerivative.allFinite())
	{
		throw std::invalid_argument(
			"Residual parameter derivative dimensions are inconsistent");
	}
	const int primalDofs = static_cast<int>(hessian.rows());
	const int contacts = static_cast<int>(constraintJacobian.rows());
	const int parameters = static_cast<int>(residualParameterDerivative.cols());
	const MatrixXd matrix = symmetricKktMatrix(hessian, constraintJacobian);
	MatrixXd rightHandSide = MatrixXd::Zero(
		primalDofs + contacts, parameters);
	rightHandSide.topRows(primalDofs) = -residualParameterDerivative;

	ConstrainedSensitivityResult result;
	Eigen::FullPivLU<MatrixXd> factorization(matrix);
	if (!factorization.isInvertible())
	{
		return result;
	}
	const MatrixXd solution = factorization.solve(rightHandSide);
	if (!solution.allFinite())
	{
		return result;
	}
	result.configurationDerivative = solution.topRows(primalDofs);
	result.multiplierDerivative = -solution.bottomRows(contacts);
	result.linearResidualNorm =
		(matrix * solution - rightHandSide).norm();
	result.success = std::isfinite(result.linearResidualNorm) &&
		result.linearResidualNorm <= 1.0e-9 *
			std::max(1.0, rightHandSide.norm());
	return result;
}
