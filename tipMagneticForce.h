#ifndef TIPMAGNETICFORCE_H
#define TIPMAGNETICFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

struct AxialTipMagneticEvaluation
{
	double energy = 0.0;
	Matrix<double, 7, 1> residual = Matrix<double, 7, 1>::Zero();
	Matrix<double, 7, 7> hessian = Matrix<double, 7, 7>::Zero();
	Matrix<double, 7, 3> residualDerivativeField = Matrix<double, 7, 3>::Zero();
};

// Permanent distal-tip dipole aligned with the final rod edge. This is the
// first planning-model sanity case. A general material-frame dipole can replace
// it after symbolic code generation is available.
class tipMagneticForce
{
public:
	tipMagneticForce(
		elasticRod &rod,
		timeStepper &stepper,
		const Vector3d &field,
		double dipoleMagnitude);

	void computeFm();
	void computeJm();
	void setField(const Vector3d &field);
	MatrixXd computeFreeDofResidualDerivativeField() const;

	static AxialTipMagneticEvaluation evaluate(
		const Vector3d &firstVertex,
		const Vector3d &secondVertex,
		const Vector3d &field,
		double dipoleMagnitude);

private:
	elasticRod *rod;
	timeStepper *stepper;
	Vector3d field;
	double dipoleMagnitude;
};

#endif
