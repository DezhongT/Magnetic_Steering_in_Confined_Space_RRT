#include "setInput.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
setInput makeInput(const Vector3d &field)
{
	setInput input;
	input.GetStringOpt("contactModel") = "planar_barrier";
	input.GetStringOpt("magneticModel") = "axial_tip";
	input.GetIntOpt("numVertices") = 20;
	input.GetIntOpt("maxIter") = 200;
	input.GetScalarOpt("dBar") = 0.015;
	input.GetScalarOpt("stiffness") = 1.0e2;
	input.GetScalarOpt("tipDipoleMoment") = 1.0e-3;
	input.GetVecOpt("gVector") = Vector3d(0.0, 0.0, 0.22);
	input.GetVecOpt("baVector") = field;
	return input;
}

bool sameContactIds(
	const std::vector<ContactCandidate> &left,
	const std::vector<ContactCandidate> &right)
{
	if (left.size() != right.size())
	{
		return false;
	}
	for (int index = 0; index < static_cast<int>(left.size()); ++index)
	{
		if (left[index].rodVertex != right[index].rodVertex ||
			left[index].boundaryId != right[index].boundaryId)
		{
			return false;
		}
	}
	return true;
}
}

int main()
{
	const Vector3d baseField(0.0, 0.1, 0.0);
	setInput baseInput = makeInput(baseField);
	world baseSimulation(baseInput);
	baseSimulation.setRodStepper();
	const ContactKktEquilibriumResult base =
		baseSimulation.solvePlanarContactKktEquilibrium();
	if (!base.success || base.activeContacts.empty())
	{
		std::cerr << "Base contact equilibrium for sensitivity failed.\n";
		return 1;
	}

	const ContactStabilityAnalysis stability =
		baseSimulation.analyzePlanarContactStability(base);
	const ContactEquilibriumSensitivity analytic =
		baseSimulation.computePlanarContactFieldSensitivity(base);
	if (!stability.valid || !stability.stable ||
		stability.nullSpaceResidualNorm > 1.0e-10 ||
		stability.hessianSymmetryError > 1.0e-10 || !analytic.success ||
		!sameContactIds(base.activeContacts, analytic.activeContacts))
	{
		std::cerr << "Contact stability or analytic sensitivity failed: stable="
				  << stability.stable << ", minimum_eigenvalue="
				  << stability.minimumEigenvalue << ", null_residual="
				  << stability.nullSpaceResidualNorm << ", symmetry="
				  << stability.hessianSymmetryError << ", sensitivity_success="
				  << analytic.success << ".\n";
		return 1;
	}

	constexpr double step = 1.0e-4;
	MatrixXd configurationDifference = MatrixXd::Zero(
		base.state.configuration.size(), 3);
	MatrixXd multiplierDifference = MatrixXd::Zero(base.multipliers.size(), 3);
	for (int component = 0; component < 3; ++component)
	{
		Vector3d plusField = baseField;
		Vector3d minusField = baseField;
		plusField[component] += step;
		minusField[component] -= step;
		baseSimulation.setAppliedField(plusField);
		const ContactKktEquilibriumResult plus =
			baseSimulation.correctPlanarContactKktEquilibrium(base);
		baseSimulation.setAppliedField(minusField);
		const ContactKktEquilibriumResult minus =
			baseSimulation.correctPlanarContactKktEquilibrium(base);
		if (!plus.success || !minus.success ||
			!sameContactIds(base.activeContacts, plus.activeContacts) ||
			!sameContactIds(base.activeContacts, minus.activeContacts))
		{
			std::cerr << "Perturbed field changed the active set or failed at component "
					  << component << ".\n";
			return 1;
		}
		configurationDifference.col(component) =
			(plus.state.configuration - minus.state.configuration) / (2.0 * step);
		multiplierDifference.col(component) =
			(plus.multipliers - minus.multipliers) / (2.0 * step);
	}

	const double configurationScale = std::max(
		1.0, std::max(analytic.configurationDerivative.norm(),
					  configurationDifference.norm()));
	const double multiplierScale = std::max(
		1.0, std::max(analytic.multiplierDerivative.norm(),
					  multiplierDifference.norm()));
	const double configurationError =
		(analytic.configurationDerivative - configurationDifference).norm() /
		configurationScale;
	const double multiplierError =
		(analytic.multiplierDerivative - multiplierDifference).norm() /
		multiplierScale;
	double positionDifference = 0.0;
	double thetaDifference = 0.0;
	for (int dof = 0; dof < analytic.configurationDerivative.rows(); ++dof)
	{
		if (dof % 4 == 3)
		{
			thetaDifference +=
				(analytic.configurationDerivative.row(dof) -
				 configurationDifference.row(dof)).squaredNorm();
		}
		else
		{
			positionDifference +=
				(analytic.configurationDerivative.row(dof) -
				 configurationDifference.row(dof)).squaredNorm();
		}
	}
	if (!std::isfinite(configurationError) || !std::isfinite(multiplierError) ||
		configurationError > 1.0e-4 || multiplierError > 1.0e-4)
	{
		std::cerr << "Contact field sensitivity finite difference failed: q_error="
				  << configurationError << ", lambda_error="
				  << multiplierError << ", analytic_q_norm="
				  << analytic.configurationDerivative.norm() << ", fd_q_norm="
				  << configurationDifference.norm() << ", difference_by_component="
				  << (analytic.configurationDerivative - configurationDifference)
					 .colwise().norm() << ", analytic_by_component="
				  << analytic.configurationDerivative.colwise().norm()
				  << ", fd_by_component="
				  << configurationDifference.colwise().norm()
				  << ", position_difference=" << std::sqrt(positionDifference)
				  << ", theta_difference=" << std::sqrt(thetaDifference) << ".\n";
		return 1;
	}

	std::cout << "Contact stability/sensitivity: minimum_eigenvalue="
			  << stability.minimumEigenvalue
			  << ", q_error=" << configurationError
			  << ", lambda_error=" << multiplierError
			  << ", linear_residual=" << analytic.linearResidualNorm << '\n';
	return 0;
}
