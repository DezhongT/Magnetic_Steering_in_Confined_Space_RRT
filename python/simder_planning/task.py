from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Dict, Sequence, Tuple

import simder

from .rrt import RRTConfig


TASK_SCHEMA_VERSION = "simder-planning-task-v1"
Vector3 = Tuple[float, float, float]


def _triple(values: Sequence[float]) -> Vector3:
    if len(values) != 3:
        raise ValueError("Expected a three-component vector")
    return tuple(float(value) for value in values)


@dataclass(frozen=True)
class MechanicsTaskConfig:
    rod_length: float = 1.0
    rod_radius: float = 1.0e-2
    num_vertices: int = 20
    young_modulus: float = 1.0e7
    poisson_ratio: float = 0.5
    maximum_newton_iterations: int = 200
    gravity: Vector3 = (0.0, 0.0, 0.0)
    domain_type: str = "planar_slab"
    plane_half_thickness: float = 1.0e-1
    shell_center: Vector3 = (0.0, 0.0, 0.0)
    shell_radius: float = 1.0
    shell_minus_thickness: float = 0.2
    shell_plus_thickness: float = 0.2
    cavity_center: Vector3 = (0.0, 0.0, 0.0)
    cavity_radius: float = 2.0
    obstacle_center: Vector3 = (0.0, 0.5, 0.0)
    obstacle_radius: float = 0.1
    second_obstacle_center: Vector3 = (0.0, -0.5, 0.0)
    second_obstacle_radius: float = 0.1
    barrier_distance: float = 1.5e-2
    barrier_stiffness: float = 1.0e2
    tip_safe_distance: float = 0.0
    tip_dipole_moment: float = 1.0e-3
    insertion_stiffness: float = 1.0e3
    insertion_axis: Vector3 = (1.0, 0.0, 0.0)
    initial_insertion: float = 1.0e-3
    initial_field: Vector3 = (0.0, 0.05, 0.0)

    def validated(self) -> "MechanicsTaskConfig":
        if self.domain_type not in (
            "planar_slab", "spherical_shell", "spherical_obstacle",
            "double_spherical_obstacle"
        ):
            raise ValueError("Unsupported confined-domain type")
        if self.plane_half_thickness <= self.rod_radius:
            raise ValueError("Planar half thickness must exceed rod radius")
        if (
            self.shell_radius <= 0.0
            or self.shell_minus_thickness < 0.0
            or self.shell_plus_thickness < 0.0
            or self.shell_radius - self.shell_minus_thickness <= self.rod_radius
        ):
            raise ValueError("Invalid spherical-shell dimensions")
        _triple(self.shell_center)
        cavity_center = _triple(self.cavity_center)
        obstacle_center = _triple(self.obstacle_center)
        second_obstacle_center = _triple(self.second_obstacle_center)
        if (
            self.cavity_radius <= 0.0
            or self.obstacle_radius <= 0.0
            or math.dist(cavity_center, obstacle_center) + self.obstacle_radius
            >= self.cavity_radius
            or self.second_obstacle_radius <= 0.0
            or math.dist(cavity_center, second_obstacle_center)
            + self.second_obstacle_radius >= self.cavity_radius
        ):
            raise ValueError("Invalid spherical-obstacle dimensions")
        return self

    def to_binding(self) -> object:
        self.validated()
        config = simder.MechanicsConfig()
        config.rod_length = self.rod_length
        config.rod_radius = self.rod_radius
        config.num_vertices = self.num_vertices
        config.young_modulus = self.young_modulus
        config.poisson_ratio = self.poisson_ratio
        config.maximum_newton_iterations = self.maximum_newton_iterations
        config.gravity = self.gravity
        config.domain_type = self.domain_type
        config.plane_half_thickness = self.plane_half_thickness
        config.shell_center = self.shell_center
        config.shell_radius = self.shell_radius
        config.shell_minus_thickness = self.shell_minus_thickness
        config.shell_plus_thickness = self.shell_plus_thickness
        config.cavity_center = self.cavity_center
        config.cavity_radius = self.cavity_radius
        config.obstacle_center = self.obstacle_center
        config.obstacle_radius = self.obstacle_radius
        config.second_obstacle_center = self.second_obstacle_center
        config.second_obstacle_radius = self.second_obstacle_radius
        config.barrier_distance = self.barrier_distance
        config.barrier_stiffness = self.barrier_stiffness
        config.tip_safe_distance = self.tip_safe_distance
        config.tip_dipole_moment = self.tip_dipole_moment
        config.insertion_stiffness = self.insertion_stiffness
        config.insertion_axis = self.insertion_axis
        config.initial_actuation = simder.Actuation(
            self.initial_insertion, self.initial_field
        )
        return config

    @classmethod
    def from_dict(cls, values: Dict[str, Any]) -> "MechanicsTaskConfig":
        converted = dict(values)
        for name in (
            "gravity",
            "shell_center",
            "cavity_center",
            "obstacle_center",
            "second_obstacle_center",
            "insertion_axis",
            "initial_field",
        ):
            if name in converted:
                converted[name] = _triple(converted[name])
        return cls(**converted)


@dataclass(frozen=True)
class ContinuationTaskConfig:
    initial_step_fraction: float = 0.25
    minimum_step_fraction: float = 1.0e-3
    maximum_step_fraction: float = 0.5
    step_reduction: float = 0.5
    step_growth: float = 1.5
    stability_tolerance: float = 0.0
    easy_corrector_iterations: int = 10
    maximum_attempts: int = 200

    def to_binding(self) -> object:
        options = simder.ContinuationOptions()
        options.initial_step_fraction = self.initial_step_fraction
        options.minimum_step_fraction = self.minimum_step_fraction
        options.maximum_step_fraction = self.maximum_step_fraction
        options.step_reduction = self.step_reduction
        options.step_growth = self.step_growth
        options.stability_tolerance = self.stability_tolerance
        options.easy_corrector_iterations = self.easy_corrector_iterations
        options.maximum_attempts = self.maximum_attempts
        return options


@dataclass(frozen=True)
class PlanningTask:
    name: str
    description: str
    mechanics: MechanicsTaskConfig
    planner: RRTConfig
    continuation: ContinuationTaskConfig
    target_tip: Vector3
    seeds: Tuple[int, ...]
    required_contact_event: str = "none"
    schema_version: str = field(default=TASK_SCHEMA_VERSION, init=False)

    def validated(self) -> "PlanningTask":
        if not self.name:
            raise ValueError("Task name must not be empty")
        if not self.seeds or any(seed < 0 for seed in self.seeds):
            raise ValueError("Task requires at least one nonnegative seed")
        if len(set(self.seeds)) != len(self.seeds):
            raise ValueError("Task seeds must be unique")
        if self.required_contact_event not in (
            "none", "insertion", "release", "any", "sustained"
        ):
            raise ValueError("Unsupported required contact event")
        _triple(self.target_tip)
        self.mechanics.validated()
        self.planner.validated()
        return self

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "name": self.name,
            "description": self.description,
            "mechanics": asdict(self.mechanics),
            "planner": asdict(self.planner),
            "continuation": asdict(self.continuation),
            "target_tip": list(self.target_tip),
            "seeds": list(self.seeds),
            "required_contact_event": self.required_contact_event,
        }

    @classmethod
    def from_dict(cls, values: Dict[str, Any]) -> "PlanningTask":
        if values.get("schema_version") != TASK_SCHEMA_VERSION:
            raise ValueError("Unsupported planning-task schema version")
        planner_values = dict(values["planner"])
        planner_values["workspace_min"] = _triple(planner_values["workspace_min"])
        planner_values["workspace_max"] = _triple(planner_values["workspace_max"])
        return cls(
            name=str(values["name"]),
            description=str(values["description"]),
            mechanics=MechanicsTaskConfig.from_dict(values["mechanics"]),
            planner=RRTConfig(**planner_values),
            continuation=ContinuationTaskConfig(**values["continuation"]),
            target_tip=_triple(values["target_tip"]),
            seeds=tuple(int(seed) for seed in values["seeds"]),
            required_contact_event=str(values.get("required_contact_event", "none")),
        ).validated()

    def save(self, path: Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path: Path) -> "PlanningTask":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def contact_free_task_a(seeds: Sequence[int] = (0, 1, 2)) -> PlanningTask:
    return PlanningTask(
        name="task_a_contact_free",
        description=(
            "Contact-free planar sanity benchmark for steering, terminal "
            "connection, and deterministic mechanics cost."
        ),
        mechanics=MechanicsTaskConfig(),
        planner=RRTConfig(
            workspace_min=(0.48, -0.01, -0.01),
            workspace_max=(0.52, 0.01, 0.01),
            goal_bias=0.5,
            local_tip_step=1.0e-3,
            goal_neighborhood_radius=2.0e-4,
            terminal_tolerance=1.0e-6,
            maximum_terminal_iterations=10,
            maximum_iterations=100,
            maximum_insertion_step=2.0e-4,
            maximum_field_step=2.0e-2,
            maximum_field_magnitude=0.25,
        ),
        continuation=ContinuationTaskConfig(),
        target_tip=(
            0.5000318800053356,
            0.0011974095086801746,
            0.0003592229347563601,
        ),
        seeds=tuple(int(seed) for seed in seeds),
    ).validated()


def planar_contact_release_task(
    seeds: Sequence[int] = (0, 1, 2),
) -> PlanningTask:
    return PlanningTask(
        name="task_planar_contact_release",
        description=(
            "Planar KKT benchmark starting in upper-wall body contact and "
            "ending on the contact-free manifold after a contact release."
        ),
        mechanics=MechanicsTaskConfig(
            gravity=(0.0, 0.0, 0.22),
            initial_field=(0.0, 0.1, 0.0),
        ),
        planner=RRTConfig(
            workspace_min=(0.48, -0.01, 0.09),
            workspace_max=(0.51, 0.01, 0.099),
            goal_bias=0.5,
            local_tip_step=1.0e-3,
            goal_neighborhood_radius=3.0e-4,
            terminal_tolerance=1.0e-5,
            maximum_terminal_iterations=10,
            maximum_iterations=100,
            maximum_insertion_step=2.0e-4,
            maximum_field_step=5.0e-2,
            maximum_field_magnitude=0.5,
        ),
        continuation=ContinuationTaskConfig(),
        target_tip=(
            0.4946049542481765,
            0.0005961682283530806,
            0.09580065071523321,
        ),
        seeds=tuple(int(seed) for seed in seeds),
        required_contact_event="release",
    ).validated()


def spherical_shell_task_b(
    seeds: Sequence[int] = (0, 1, 2),
) -> PlanningTask:
    return PlanningTask(
        name="task_b_spherical_shell",
        description=(
            "Curved-contact benchmark following the inner surface of an "
            "analytic spherical shell with one sustained body contact."
        ),
        mechanics=MechanicsTaskConfig(
            gravity=(0.0, 0.0, 0.5),
            domain_type="spherical_shell",
            shell_center=(0.0, 0.0, 1.0),
            shell_radius=1.0,
            shell_minus_thickness=0.1,
            shell_plus_thickness=1.0,
            initial_insertion=1.0e-3,
            initial_field=(0.0, 0.1, 0.0),
        ),
        planner=RRTConfig(
            workspace_min=(0.4785, -0.001, 0.1835),
            workspace_max=(0.481, 0.003, 0.186),
            goal_bias=0.5,
            local_tip_step=2.0e-4,
            goal_neighborhood_radius=1.0e-4,
            terminal_tolerance=2.0e-6,
            maximum_terminal_iterations=12,
            maximum_iterations=100,
            maximum_insertion_step=5.0e-4,
            maximum_field_step=5.0e-2,
            maximum_field_magnitude=0.5,
        ),
        continuation=ContinuationTaskConfig(),
        target_tip=(
            0.47972378621565775,
            0.0011594835657040834,
            0.18460557852902632,
        ),
        seeds=tuple(int(seed) for seed in seeds),
        required_contact_event="sustained",
    ).validated()


def spherical_obstacle_dead_end_task(
    seeds: Sequence[int] = (0, 1, 2),
) -> PlanningTask:
    return PlanningTask(
        name="task_d_spherical_obstacle_dead_end",
        description=(
            "Tip-obstacle dead end: direct goal seeking violates the forbidden "
            "tip clearance, while a successful route first moves around an "
            "eccentric spherical obstacle."
        ),
        mechanics=MechanicsTaskConfig(
            domain_type="spherical_obstacle",
            cavity_center=(0.0, 0.0, 0.0),
            cavity_radius=2.0,
            obstacle_center=(
                0.5000485204990591,
                0.0008915111473809573,
                0.00029650602658044003,
            ),
            obstacle_radius=2.5e-4,
            tip_safe_distance=1.0e-4,
            initial_insertion=1.0e-3,
            initial_field=(0.0, 0.05, 0.0),
        ),
        planner=RRTConfig(
            workspace_min=(0.4995, -0.0005, -0.002),
            workspace_max=(0.5005, 0.0025, 0.002),
            goal_bias=0.2,
            local_tip_step=2.0e-4,
            goal_neighborhood_radius=1.0e-4,
            terminal_tolerance=2.0e-6,
            maximum_terminal_iterations=15,
            maximum_iterations=300,
            maximum_insertion_step=5.0e-4,
            maximum_field_step=5.0e-2,
            maximum_field_magnitude=0.5,
        ),
        continuation=ContinuationTaskConfig(
            initial_step_fraction=0.15,
            minimum_step_fraction=1.0e-4,
            maximum_step_fraction=0.3,
            maximum_attempts=300,
        ),
        target_tip=(
            0.5000806258711333,
            0.0014825224746173426,
            0.0005930089898471434,
        ),
        seeds=tuple(int(seed) for seed in seeds),
    ).validated()


def double_obstacle_contact_branching_task(
    seeds: Sequence[int] = (0, 1, 2),
) -> PlanningTask:
    return PlanningTask(
        name="task_e_double_obstacle_contact_branching",
        description=(
            "Body-contact-history branch: a direct upward actuation becomes "
            "pinned between two obstacles, while positive- and negative-side "
            "routes insert and release different stable contacts before "
            "reaching the same contact-free goal."
        ),
        mechanics=MechanicsTaskConfig(
            domain_type="double_spherical_obstacle",
            cavity_center=(0.0, 0.0, 0.0),
            cavity_radius=3.0,
            obstacle_center=(0.23684, 0.0108, 0.003),
            obstacle_radius=1.0e-3,
            second_obstacle_center=(0.23684, -0.0108, 0.003),
            second_obstacle_radius=1.0e-3,
            barrier_distance=1.0e-5,
            tip_dipole_moment=1.0e-2,
            initial_insertion=1.0e-3,
            initial_field=(0.0, 0.0, 0.0),
        ),
        planner=RRTConfig(
            workspace_min=(0.495, -0.055, -0.005),
            workspace_max=(0.505, 0.055, 0.022),
            goal_bias=0.2,
            local_tip_step=2.0e-3,
            goal_neighborhood_radius=5.0e-4,
            terminal_tolerance=2.0e-6,
            maximum_terminal_iterations=20,
            maximum_iterations=1000,
            maximum_insertion_step=5.0e-4,
            maximum_field_step=1.0e-1,
            maximum_field_magnitude=1.0,
        ),
        continuation=ContinuationTaskConfig(
            initial_step_fraction=0.05,
            minimum_step_fraction=1.0e-7,
            maximum_step_fraction=0.3,
            maximum_attempts=1000,
        ),
        target_tip=(
            0.49986371556769144,
            0.0,
            0.01779815499963202,
        ),
        seeds=tuple(int(seed) for seed in seeds),
        required_contact_event="any",
    ).validated()
