#include "planning/mechanicsSession.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <vector>

namespace py = pybind11;

namespace
{
Vector3d vector3FromPython(const py::sequence &values)
{
	if (py::len(values) != 3)
	{
		throw py::value_error("Expected a three-component vector");
	}
	return Vector3d(
		py::cast<double>(values[0]),
		py::cast<double>(values[1]),
		py::cast<double>(values[2]));
}

std::array<double, 3> vector3ToPython(const Vector3d &value)
{
	return {value[0], value[1], value[2]};
}

std::vector<double> vectorToPython(const VectorXd &value)
{
	return std::vector<double>(value.data(), value.data() + value.size());
}

std::array<std::array<double, 4>, 3> matrix34ToPython(
	const Matrix<double, 3, 4> &value)
{
	std::array<std::array<double, 4>, 3> result{};
	for (int row = 0; row < 3; ++row)
	{
		for (int column = 0; column < 4; ++column)
		{
			result[row][column] = value(row, column);
		}
	}
	return result;
}
}

PYBIND11_MODULE(simder, module)
{
	module.doc() = "High-level mechanics interface for simDER planning";

	py::class_<Actuation>(module, "Actuation")
		.def(py::init<>())
		.def(py::init([](double xi, const py::sequence &field) {
			Actuation actuation;
			actuation.xi = xi;
			actuation.field = vector3FromPython(field);
			return actuation;
		}), py::arg("xi"), py::arg("field"))
		.def_readwrite("xi", &Actuation::xi)
		.def_property("field",
			[](const Actuation &actuation) {
				return vector3ToPython(actuation.field);
			},
			[](Actuation &actuation, const py::sequence &field) {
				actuation.field = vector3FromPython(field);
			});

	py::class_<MechanicsConfig>(module, "MechanicsConfig")
		.def(py::init<>())
		.def_readwrite("rod_length", &MechanicsConfig::rodLength)
		.def_readwrite("rod_radius", &MechanicsConfig::rodRadius)
		.def_readwrite("num_vertices", &MechanicsConfig::numVertices)
		.def_readwrite("young_modulus", &MechanicsConfig::youngModulus)
		.def_readwrite("poisson_ratio", &MechanicsConfig::poissonRatio)
		.def_readwrite("maximum_newton_iterations",
			&MechanicsConfig::maximumNewtonIterations)
		.def_property("gravity",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.gravity);
			},
			[](MechanicsConfig &config, const py::sequence &gravity) {
				config.gravity = vector3FromPython(gravity);
			})
		.def_readwrite("domain_type", &MechanicsConfig::domainType)
		.def_readwrite("plane_half_thickness",
			&MechanicsConfig::planeHalfThickness)
		.def_property("shell_center",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.shellCenter);
			},
			[](MechanicsConfig &config, const py::sequence &center) {
				config.shellCenter = vector3FromPython(center);
			})
		.def_readwrite("shell_radius", &MechanicsConfig::shellRadius)
		.def_readwrite("shell_minus_thickness",
			&MechanicsConfig::shellMinusThickness)
		.def_readwrite("shell_plus_thickness",
			&MechanicsConfig::shellPlusThickness)
		.def_property("cavity_center",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.cavityCenter);
			},
			[](MechanicsConfig &config, const py::sequence &center) {
				config.cavityCenter = vector3FromPython(center);
			})
		.def_readwrite("cavity_radius", &MechanicsConfig::cavityRadius)
		.def_property("obstacle_center",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.obstacleCenter);
			},
			[](MechanicsConfig &config, const py::sequence &center) {
				config.obstacleCenter = vector3FromPython(center);
			})
		.def_readwrite("obstacle_radius", &MechanicsConfig::obstacleRadius)
		.def_property("second_obstacle_center",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.secondObstacleCenter);
			},
			[](MechanicsConfig &config, const py::sequence &center) {
				config.secondObstacleCenter = vector3FromPython(center);
			})
		.def_readwrite("second_obstacle_radius",
			&MechanicsConfig::secondObstacleRadius)
		.def_readwrite("barrier_distance", &MechanicsConfig::barrierDistance)
		.def_readwrite("barrier_stiffness", &MechanicsConfig::barrierStiffness)
		.def_readwrite("tip_safe_distance", &MechanicsConfig::tipSafeDistance)
		.def_readwrite("tip_dipole_moment", &MechanicsConfig::tipDipoleMoment)
		.def_readwrite("insertion_stiffness",
			&MechanicsConfig::insertionStiffness)
		.def_property("insertion_axis",
			[](const MechanicsConfig &config) {
				return vector3ToPython(config.insertionAxis);
			},
			[](MechanicsConfig &config, const py::sequence &axis) {
				config.insertionAxis = vector3FromPython(axis);
			})
		.def_readwrite("initial_actuation", &MechanicsConfig::initialActuation);

	py::class_<FieldContinuationOptions>(module, "ContinuationOptions")
		.def(py::init<>())
		.def_readwrite("initial_step_fraction",
			&FieldContinuationOptions::initialStepFraction)
		.def_readwrite("minimum_step_fraction",
			&FieldContinuationOptions::minimumStepFraction)
		.def_readwrite("maximum_step_fraction",
			&FieldContinuationOptions::maximumStepFraction)
		.def_readwrite("step_reduction", &FieldContinuationOptions::stepReduction)
		.def_readwrite("step_growth", &FieldContinuationOptions::stepGrowth)
		.def_readwrite("stability_tolerance",
			&FieldContinuationOptions::stabilityTolerance)
		.def_readwrite("easy_corrector_iterations",
			&FieldContinuationOptions::easyCorrectorIterations)
		.def_readwrite("maximum_attempts", &FieldContinuationOptions::maximumAttempts);

	py::class_<PlannerState>(module, "PlannerState")
		.def_property_readonly("configuration", [](const PlannerState &state) {
			return vectorToPython(state.rodState.configuration);
		})
		.def_property_readonly("actuation", [](const PlannerState &state) {
			return state.actuation;
		})
		.def_property_readonly("tip_position", [](const PlannerState &state) {
			if (state.rodState.configuration.size() < 3)
			{
				throw py::value_error("PlannerState configuration has no distal tip");
			}
			return vector3ToPython(state.rodState.configuration.tail<3>());
		})
		.def_property_readonly("active_contact_ids", [](const PlannerState &state) {
			std::vector<std::array<int, 2>> ids;
			ids.reserve(state.activeContacts.size());
			for (const ContactCandidate &candidate : state.activeContacts)
			{
				ids.push_back({candidate.rodVertex, candidate.boundaryId});
			}
			return ids;
		})
		.def_property_readonly("multipliers", [](const PlannerState &state) {
			return vectorToPython(state.multipliers);
		})
		.def_readonly("stability_margin", &PlannerState::stabilityMargin)
		.def_readonly("stationarity_norm", &PlannerState::stationarityNorm)
		.def_readonly("active_constraint_norm", &PlannerState::activeConstraintNorm)
		.def_readonly("complementarity_norm", &PlannerState::complementarityNorm)
		.def_readonly("minimum_body_gap", &PlannerState::minimumBodyGap)
		.def_readonly("tip_clearance", &PlannerState::tipClearance)
		.def_readonly("tip_safe", &PlannerState::tipSafe);

	py::class_<LocalSteeringResult>(module, "LocalSteeringResult")
		.def_property_readonly("tip_position", [](const LocalSteeringResult &result) {
			return vector3ToPython(result.tipPosition);
		})
		.def_property_readonly("tip_actuation_jacobian",
			[](const LocalSteeringResult &result) {
				return matrix34ToPython(result.tipActuationJacobian);
			})
		.def_readonly("linear_residual_norm",
			&LocalSteeringResult::linearResidualNorm);

	py::class_<ContinuationEdgeResult>(module, "ContinuationResult")
		.def_readonly("success", &ContinuationEdgeResult::success)
		.def_readonly("rolled_back", &ContinuationEdgeResult::rolledBack)
		.def_readonly("state", &ContinuationEdgeResult::state)
		.def_readonly("attempted_steps", &ContinuationEdgeResult::attemptedSteps)
		.def_readonly("rejected_steps", &ContinuationEdgeResult::rejectedSteps)
		.def_readonly("stored_points", &ContinuationEdgeResult::storedPoints)
		.def_readonly("contacts_added", &ContinuationEdgeResult::contactsAdded)
		.def_readonly("contacts_released", &ContinuationEdgeResult::contactsReleased)
		.def_readonly("reached_path_fraction",
			&ContinuationEdgeResult::reachedPathFraction)
		.def_readonly("minimum_stability_margin",
			&ContinuationEdgeResult::minimumStabilityMargin)
		.def_readonly("failure_reason", &ContinuationEdgeResult::failureReason);

	py::class_<MechanicsSession>(module, "MechanicsSession")
		.def(py::init<const MechanicsConfig &>(),
			py::arg("config") = MechanicsConfig())
		.def("solve_initial_state", &MechanicsSession::solveInitialState)
		.def("evaluate_local_steering",
			&MechanicsSession::evaluateLocalSteering,
			py::arg("state"))
		.def("attempt_continuation", &MechanicsSession::attemptContinuation,
			py::arg("start_state"), py::arg("target_actuation"),
			py::arg("options") = FieldContinuationOptions());
}
