#include "setInput.h"
#include "tipMagneticForce.h"
#include "world.h"

#include <cmath>
#include <iostream>

namespace
{
bool checkLocalDerivativesAndInvariants()
{
	Vector3d x1(-0.1, 0.02, 0.03);
	Vector3d x2(0.3, -0.04, 0.2);
	Vector3d field(0.4, -0.7, 1.2);
	constexpr double moment = 0.7;
	constexpr double step = 1.0e-6;
	const AxialTipMagneticEvaluation evaluation =
		tipMagneticForce::evaluate(x1, x2, field, moment);

	Matrix<double, 7, 1> finiteDifferenceGradient;
	for (int dof = 0; dof < 7; ++dof)
	{
		Vector3d plusX1 = x1;
		Vector3d plusX2 = x2;
		Vector3d minusX1 = x1;
		Vector3d minusX2 = x2;
		if (dof < 3)
		{
			plusX1[dof] += step;
			minusX1[dof] -= step;
		}
		else if (dof > 3)
		{
			plusX2[dof - 4] += step;
			minusX2[dof - 4] -= step;
		}
		const double plusEnergy =
			tipMagneticForce::evaluate(plusX1, plusX2, field, moment).energy;
		const double minusEnergy =
			tipMagneticForce::evaluate(minusX1, minusX2, field, moment).energy;
		finiteDifferenceGradient[dof] = (plusEnergy - minusEnergy) / (2.0 * step);
	}

	Matrix<double, 7, 7> finiteDifferenceHessian;
	for (int dof = 0; dof < 7; ++dof)
	{
		Vector3d plusX1 = x1;
		Vector3d plusX2 = x2;
		Vector3d minusX1 = x1;
		Vector3d minusX2 = x2;
		if (dof < 3)
		{
			plusX1[dof] += step;
			minusX1[dof] -= step;
		}
		else if (dof > 3)
		{
			plusX2[dof - 4] += step;
			minusX2[dof - 4] -= step;
		}
		finiteDifferenceHessian.col(dof) =
			(tipMagneticForce::evaluate(plusX1, plusX2, field, moment).residual -
			 tipMagneticForce::evaluate(minusX1, minusX2, field, moment).residual) /
			(2.0 * step);
	}

	Matrix<double, 7, 3> finiteDifferenceFieldDerivative;
	for (int component = 0; component < 3; ++component)
	{
		Vector3d plusField = field;
		Vector3d minusField = field;
		plusField[component] += step;
		minusField[component] -= step;
		finiteDifferenceFieldDerivative.col(component) =
			(tipMagneticForce::evaluate(x1, x2, plusField, moment).residual -
			 tipMagneticForce::evaluate(x1, x2, minusField, moment).residual) /
			(2.0 * step);
	}

	const double gradientError =
		(finiteDifferenceGradient - evaluation.residual).norm();
	const double hessianError =
		(finiteDifferenceHessian - evaluation.hessian).norm();
	const double fieldDerivativeError =
		(finiteDifferenceFieldDerivative - evaluation.residualDerivativeField).norm();
	const double symmetryError =
		(evaluation.hessian - evaluation.hessian.transpose()).norm();

	const Vector3d physicalForce1 = -evaluation.residual.segment<3>(0);
	const Vector3d physicalForce2 = -evaluation.residual.segment<3>(4);
	const Vector3d midpoint = 0.5 * (x1 + x2);
	const Vector3d torque =
		(x1 - midpoint).cross(physicalForce1) +
		(x2 - midpoint).cross(physicalForce2);
	const Vector3d expectedTorque = moment * ((x2 - x1).normalized()).cross(field);

	std::cout << "Tip magnetic derivative errors: gradient=" << gradientError
			  << ", hessian=" << hessianError
			  << ", field=" << fieldDerivativeError << '\n';

	return gradientError < 1.0e-9 && hessianError < 1.0e-8 &&
		fieldDerivativeError < 1.0e-9 && symmetryError < 1.0e-12 &&
		(physicalForce1 + physicalForce2).norm() < 1.0e-12 &&
		(torque - expectedTorque).norm() < 1.0e-12;
}

bool checkMagneticBending(const char *optionFile)
{
	setInput input;
	if (input.LoadOptions(optionFile) != 0)
	{
		return false;
	}
	world simulation(input);
	simulation.setRodStepper();
	const RodState initial = simulation.captureRodState();
	const EquilibriumResult result = simulation.solveStaticEquilibrium();
	if (!result.success)
	{
		std::cerr << "Axial tip magnetic static solve failed: residual="
				  << result.finalResidualNorm << '\n';
		return false;
	}
	const int tipPositionIndex = 4 * (simulation.numPoints() - 1);
	const Vector3d displacement =
		result.state.configuration.segment<3>(tipPositionIndex) -
		initial.configuration.segment<3>(tipPositionIndex);
	std::cout << "Tip magnetic equilibrium: iterations=" << result.newtonIterations
			  << ", residual=" << result.finalResidualNorm
			  << ", tip_displacement=" << displacement.transpose() << '\n';
	return displacement.norm() > 1.0e-4 && result.finalResidualNorm < 1.0e-7;
}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <option-file>\n";
		return 1;
	}
	if (!checkLocalDerivativesAndInvariants() || !checkMagneticBending(argv[1]))
	{
		return 1;
	}
	return 0;
}
