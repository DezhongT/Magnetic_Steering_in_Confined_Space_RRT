#ifndef PROXIMALGUIDEINSERTIONMODEL_H
#define PROXIMALGUIDEINSERTIONMODEL_H

#include "elasticRod.h"
#include "insertion/insertionModel.h"
#include "timeStepper.h"

struct ProximalGuideEvaluation
{
	double energy = 0.0;
	Vector3d residual = Vector3d::Zero();
	Matrix3d hessian = Matrix3d::Zero();
	Vector3d residualDerivativeCoordinate = Vector3d::Zero();
};

class ProximalGuideInsertionModel : public InsertionModel
{
public:
	ProximalGuideInsertionModel(
		elasticRod &rod,
		timeStepper &stepper,
		const Vector3d &axis,
		double stiffness,
		double initialCoordinate);

	void setCoordinate(double xi) override;
	double coordinate() const override;
	void assemble() override;
	VectorXd freeResidualDerivativeCoordinate() const override;
	Vector3d guideTarget() const override;
	std::string modelName() const override;

	static ProximalGuideEvaluation evaluate(
		const Vector3d &position,
		const Vector3d &zeroCoordinateTarget,
		const Vector3d &axis,
		double stiffness,
		double xi);

private:
	elasticRod *rod;
	timeStepper *stepper;
	Vector3d axis;
	Vector3d zeroCoordinateTarget;
	double stiffness;
	double xi;
	int drivenVertex;
};

#endif
