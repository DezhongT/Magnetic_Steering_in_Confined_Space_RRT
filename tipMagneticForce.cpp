#include "tipMagneticForce.h"

#include <stdexcept>

tipMagneticForce::tipMagneticForce(
	elasticRod &m_rod,
	timeStepper &m_stepper,
	const Vector3d &m_field,
	double m_dipoleMagnitude)
	: rod(&m_rod),
	  stepper(&m_stepper),
	  field(m_field),
	  dipoleMagnitude(m_dipoleMagnitude)
{
}

void tipMagneticForce::setField(const Vector3d &newField)
{
	field = newField;
}

MatrixXd tipMagneticForce::computeFreeDofResidualDerivativeField() const
{
	const int edgeIndex = rod->ne - 1;
	const AxialTipMagneticEvaluation evaluation = evaluate(
		rod->getVertex(edgeIndex), rod->getVertex(edgeIndex + 1),
		field, dipoleMagnitude);
	MatrixXd derivative = MatrixXd::Zero(rod->uncons, 3);
	for (int local = 0; local < 7; ++local)
	{
		const int fullIndex = 4 * edgeIndex + local;
		if (rod->getIfConstrained(fullIndex) == 0)
		{
			derivative.row(rod->fullToUnconsMap[fullIndex]) =
				evaluation.residualDerivativeField.row(local);
		}
	}
	return derivative;
}

AxialTipMagneticEvaluation tipMagneticForce::evaluate(
	const Vector3d &firstVertex,
	const Vector3d &secondVertex,
	const Vector3d &appliedField,
	double momentMagnitude)
{
	const Vector3d edge = secondVertex - firstVertex;
	const double length = edge.norm();
	if (!(length > 0.0))
	{
		throw std::invalid_argument("Tip magnetic edge must have positive length");
	}

	const Vector3d tangent = edge / length;
	const Matrix3d projection = Matrix3d::Identity() - tangent * tangent.transpose();
	const double tangentField = tangent.dot(appliedField);
	const Vector3d transverseField = appliedField - tangentField * tangent;

	AxialTipMagneticEvaluation evaluation;
	evaluation.energy = -momentMagnitude * tangentField;
	const Vector3d edgeGradient =
		-momentMagnitude / length * transverseField;
	evaluation.residual.segment<3>(0) = -edgeGradient;
	evaluation.residual.segment<3>(4) = edgeGradient;

	const Matrix3d edgeHessian = momentMagnitude / (length * length) *
		(transverseField * tangent.transpose() +
		 tangent * transverseField.transpose() +
		 tangentField * projection);
	evaluation.hessian.block<3, 3>(0, 0) = edgeHessian;
	evaluation.hessian.block<3, 3>(0, 4) = -edgeHessian;
	evaluation.hessian.block<3, 3>(4, 0) = -edgeHessian;
	evaluation.hessian.block<3, 3>(4, 4) = edgeHessian;

	const Matrix3d edgeGradientField = -momentMagnitude / length * projection;
	evaluation.residualDerivativeField.block<3, 3>(0, 0) = -edgeGradientField;
	evaluation.residualDerivativeField.block<3, 3>(4, 0) = edgeGradientField;
	return evaluation;
}

void tipMagneticForce::computeFm()
{
	const int edgeIndex = rod->ne - 1;
	const AxialTipMagneticEvaluation evaluation = evaluate(
		rod->getVertex(edgeIndex), rod->getVertex(edgeIndex + 1),
		field, dipoleMagnitude);
	for (int local = 0; local < 7; ++local)
	{
		stepper->addForce(4 * edgeIndex + local, evaluation.residual[local]);
	}
}

void tipMagneticForce::computeJm()
{
	const int edgeIndex = rod->ne - 1;
	const AxialTipMagneticEvaluation evaluation = evaluate(
		rod->getVertex(edgeIndex), rod->getVertex(edgeIndex + 1),
		field, dipoleMagnitude);
	for (int row = 0; row < 7; ++row)
	{
		for (int column = 0; column < 7; ++column)
		{
			stepper->addJacobian(
				4 * edgeIndex + row,
				4 * edgeIndex + column,
				evaluation.hessian(row, column));
		}
	}
}
