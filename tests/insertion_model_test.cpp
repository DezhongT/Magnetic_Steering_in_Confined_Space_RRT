#include "insertion/proximalGuideInsertionModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int main()
{
	const Vector3d position(0.3, -0.2, 0.4);
	const Vector3d target(-0.1, 0.2, 0.05);
	const Vector3d axis = Vector3d(1.0, 2.0, -1.0).normalized();
	constexpr double stiffness = 31.0;
	constexpr double xi = 0.07;
	constexpr double step = 1.0e-7;
	const ProximalGuideEvaluation center = ProximalGuideInsertionModel::evaluate(
		position, target, axis, stiffness, xi);
	const Vector3d direction = Vector3d(-0.4, 0.3, 0.8).normalized();
	const ProximalGuideEvaluation positionPlus =
		ProximalGuideInsertionModel::evaluate(
			position + step * direction, target, axis, stiffness, xi);
	const ProximalGuideEvaluation positionMinus =
		ProximalGuideInsertionModel::evaluate(
			position - step * direction, target, axis, stiffness, xi);
	const ProximalGuideEvaluation xiPlus = ProximalGuideInsertionModel::evaluate(
		position, target, axis, stiffness, xi + step);
	const ProximalGuideEvaluation xiMinus = ProximalGuideInsertionModel::evaluate(
		position, target, axis, stiffness, xi - step);

	const double energyDerivative =
		(positionPlus.energy - positionMinus.energy) / (2.0 * step);
	const Vector3d residualPositionDerivative =
		(positionPlus.residual - positionMinus.residual) / (2.0 * step);
	const Vector3d residualXiDerivative =
		(xiPlus.residual - xiMinus.residual) / (2.0 * step);
	const double scale = std::max(1.0, center.residual.norm());
	if (std::abs(energyDerivative - center.residual.dot(direction)) / scale > 1.0e-8 ||
		(residualPositionDerivative - center.hessian * direction).norm() / scale >
			1.0e-8 ||
		(residualXiDerivative - center.residualDerivativeCoordinate).norm() /
			std::max(1.0, center.residualDerivativeCoordinate.norm()) > 1.0e-8)
	{
		std::cerr << "Proximal-guide insertion derivatives failed finite differences.\n";
		return 1;
	}
	std::cout << "Proximal-guide insertion derivative checks passed.\n";
	return 0;
}
