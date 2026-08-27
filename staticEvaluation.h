#ifndef STATICEVALUATION_H
#define STATICEVALUATION_H

#include "eigenIncludes.h"

#include <algorithm>
#include <stdexcept>

// An owned snapshot of a static residual and its LAPACK general-banded
// Jacobian. LAPACK dgbsv overwrites both inputs, so evaluations must remain
// separate from solver workspaces.
struct StaticEvaluation
{
	VectorXd residual;
	MatrixXd bandedJacobian;
	int lowerBandwidth = 0;
	int upperBandwidth = 0;

	double residualNorm() const
	{
		return residual.norm();
	}

	MatrixXd denseJacobian() const
	{
		const int systemSize = static_cast<int>(residual.size());
		const int expectedRows = 2 * lowerBandwidth + upperBandwidth + 1;
		if (bandedJacobian.rows() != expectedRows ||
			bandedJacobian.cols() != systemSize)
		{
			throw std::invalid_argument(
				"StaticEvaluation dimensions are inconsistent with its band structure");
		}

		MatrixXd dense = MatrixXd::Zero(systemSize, systemSize);
		for (int column = 0; column < systemSize; ++column)
		{
			const int firstRow = std::max(0, column - upperBandwidth);
			const int lastRow = std::min(systemSize - 1, column + lowerBandwidth);
			for (int row = firstRow; row <= lastRow; ++row)
			{
				const int bandRow = lowerBandwidth + upperBandwidth + row - column;
				dense(row, column) = bandedJacobian(bandRow, column);
			}
		}
		return dense;
	}

	VectorXd multiplyJacobian(const VectorXd &direction) const
	{
		const int systemSize = static_cast<int>(residual.size());
		const int expectedRows = 2 * lowerBandwidth + upperBandwidth + 1;
		if (direction.size() != systemSize ||
			bandedJacobian.rows() != expectedRows ||
			bandedJacobian.cols() != systemSize)
		{
			throw std::invalid_argument(
				"StaticEvaluation dimensions are inconsistent with its band structure");
		}

		VectorXd product = VectorXd::Zero(systemSize);
		for (int column = 0; column < systemSize; ++column)
		{
			const int firstRow = std::max(0, column - upperBandwidth);
			const int lastRow = std::min(systemSize - 1, column + lowerBandwidth);
			for (int row = firstRow; row <= lastRow; ++row)
			{
				const int bandRow = lowerBandwidth + upperBandwidth + row - column;
				product[row] += bandedJacobian(bandRow, column) * direction[column];
			}
		}
		return product;
	}
};

#endif
