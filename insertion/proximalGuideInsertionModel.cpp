#include "insertion/proximalGuideInsertionModel.h"

#include <cmath>
#include <stdexcept>

ProximalGuideInsertionModel::ProximalGuideInsertionModel(
	elasticRod &m_rod,
	timeStepper &m_stepper,
	const Vector3d &m_axis,
	double m_stiffness,
	double initialCoordinate)
	: rod(&m_rod),
	  stepper(&m_stepper),
	  stiffness(m_stiffness),
	  xi(initialCoordinate),
	  drivenVertex(2)
{
	if (rod->nv <= drivenVertex)
	{
		throw std::invalid_argument(
			"Proximal-guide insertion requires at least three rod vertices");
	}
	if (!m_axis.allFinite() || !(m_axis.norm() > 0.0))
	{
		throw std::invalid_argument("Insertion axis must be finite and nonzero");
	}
	if (!std::isfinite(stiffness) || stiffness <= 0.0 ||
		!std::isfinite(xi) || xi < 0.0)
	{
		throw std::invalid_argument(
			"Insertion stiffness must be positive and coordinate nonnegative");
	}
	axis = m_axis.normalized();
	zeroCoordinateTarget = rod->getVertex(drivenVertex);
}

ProximalGuideEvaluation ProximalGuideInsertionModel::evaluate(
	const Vector3d &position,
	const Vector3d &zeroCoordinateTarget,
	const Vector3d &axis,
	double stiffness,
	double xi)
{
	if (!position.allFinite() || !zeroCoordinateTarget.allFinite() ||
		!axis.allFinite() || std::abs(axis.norm() - 1.0) > 1.0e-10 ||
		!std::isfinite(stiffness) || stiffness <= 0.0 ||
		!std::isfinite(xi) || xi < 0.0)
	{
		throw std::invalid_argument("Invalid proximal-guide insertion evaluation");
	}
	const double axialError =
		axis.dot(position - zeroCoordinateTarget) - xi;
	ProximalGuideEvaluation result;
	result.energy = 0.5 * stiffness * axialError * axialError;
	result.residual = stiffness * axialError * axis;
	result.hessian = stiffness * axis * axis.transpose();
	result.residualDerivativeCoordinate = -stiffness * axis;
	return result;
}

void ProximalGuideInsertionModel::setCoordinate(double coordinate)
{
	if (!std::isfinite(coordinate) || coordinate < 0.0)
	{
		throw std::invalid_argument("Insertion coordinate must be finite and nonnegative");
	}
	xi = coordinate;
}

double ProximalGuideInsertionModel::coordinate() const
{
	return xi;
}

Vector3d ProximalGuideInsertionModel::guideTarget() const
{
	return zeroCoordinateTarget + xi * axis;
}

std::string ProximalGuideInsertionModel::modelName() const
{
	return "proximal_guide";
}

void ProximalGuideInsertionModel::assemble()
{
	const ProximalGuideEvaluation evaluation = evaluate(
		rod->getVertex(drivenVertex), zeroCoordinateTarget, axis, stiffness, xi);
	for (int row = 0; row < 3; ++row)
	{
		const int rowDof = 4 * drivenVertex + row;
		stepper->addForce(rowDof, evaluation.residual[row]);
		for (int column = 0; column < 3; ++column)
		{
			stepper->addJacobian(
				rowDof, 4 * drivenVertex + column,
				evaluation.hessian(row, column));
		}
	}
}

VectorXd ProximalGuideInsertionModel::freeResidualDerivativeCoordinate() const
{
	const ProximalGuideEvaluation evaluation = evaluate(
		rod->getVertex(drivenVertex), zeroCoordinateTarget, axis, stiffness, xi);
	VectorXd derivative = VectorXd::Zero(rod->uncons);
	for (int component = 0; component < 3; ++component)
	{
		const int fullDof = 4 * drivenVertex + component;
		if (rod->getIfConstrained(fullDof) == 0)
		{
			derivative[rod->fullToUnconsMap[fullDof]] =
				evaluation.residualDerivativeCoordinate[component];
		}
	}
	return derivative;
}
