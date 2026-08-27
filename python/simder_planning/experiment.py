from dataclasses import asdict, dataclass, replace
from collections import Counter
import json
import math
from pathlib import Path
import statistics
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import simder

from .rrt import LocalGoalSeekingPlanner, RRTNode, TipSpaceRRT
from .task import PlanningTask


EXPERIMENT_SCHEMA_VERSION = "simder-planning-experiment-v1"
COMPARISON_SCHEMA_VERSION = "simder-planning-comparison-v1"
PLANNING_METHODS = (
    "mechanics_informed_rrt",
    "random_actuation_rrt",
    "local_goal_seeking",
)


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _distance(left: Sequence[float], right: Sequence[float]) -> float:
    return _norm([a - b for a, b in zip(left, right)])


@dataclass(frozen=True)
class LoggedState:
    node_index: int
    parent_index: Optional[int]
    tip_position: Tuple[float, float, float]
    insertion: float
    field: Tuple[float, float, float]
    active_contact_ids: Tuple[Tuple[int, int], ...]
    contact_history_boundaries: Tuple[int, ...]
    minimum_contact_multiplier: Optional[float]
    stability_margin: float
    minimum_body_gap: float
    tip_clearance: float
    tip_safe: bool


@dataclass(frozen=True)
class RunRecord:
    method: str
    seed: int
    success: bool
    termination_reason: str
    continuation_step_budget: Optional[int]
    continuation_step_budget_exhausted: bool
    final_tip_error: float
    rrt_iterations: int
    attempted_edges: int
    accepted_edges: int
    accepted_node_count: int
    failed_continuations: int
    tip_safety_failures: int
    continuation_attempted_steps: int
    continuation_rejected_steps: int
    mean_continuation_steps_per_accepted_edge: float
    terminal_attempts: int
    terminal_failures: int
    terminal_steps: int
    terminal_failure_reasons: Tuple[str, ...]
    minimum_terminal_projected_residual: Optional[float]
    final_terminal_projected_residual: Optional[float]
    final_terminal_projected_residual_ratio: Optional[float]
    terminal_insertion_constraint_activations: int
    nearest_neighbor_queries: int
    nearest_neighbor_distance_evaluations: int
    shortlist_candidates_evaluated: int
    shortlist_non_nearest_selections: int
    shortlist_triggered_queries: int
    tree_contact_history_counts: Tuple[Tuple[str, int], ...]
    minimum_goal_distance_by_contact_history: Tuple[Tuple[str, float], ...]
    contact_transitions: int
    path_contact_transitions: int
    peak_field: float
    total_insertion: float
    minimum_stability_margin: float
    minimum_body_gap: float
    minimum_tip_clearance: float
    minimum_contact_multiplier: Optional[float]
    planning_time_seconds: float
    path: Tuple[LoggedState, ...]
    tree: Optional[Tuple[LoggedState, ...]] = None

    def deterministic_dict(self) -> Dict[str, Any]:
        values = asdict(self)
        values.pop("planning_time_seconds")
        return values


@dataclass(frozen=True)
class ExperimentRecord:
    task: Dict[str, Any]
    method: str
    runs: Tuple[RunRecord, ...]
    schema_version: str = EXPERIMENT_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "task": self.task,
            "method": self.method,
            "runs": [asdict(run) for run in self.runs],
        }

    def save(self, path: Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


@dataclass(frozen=True)
class MethodAggregate:
    method: str
    run_count: int
    success_rate: float
    median_final_tip_error: float
    median_planning_time_seconds: float
    median_continuation_attempted_steps: float

    def deterministic_dict(self) -> Dict[str, Any]:
        values = asdict(self)
        values.pop("median_planning_time_seconds")
        return values


@dataclass(frozen=True)
class ComparisonRecord:
    task: Dict[str, Any]
    experiments: Tuple[ExperimentRecord, ...]
    aggregates: Tuple[MethodAggregate, ...]
    schema_version: str = COMPARISON_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "task": self.task,
            "experiments": [experiment.to_dict() for experiment in self.experiments],
            "aggregates": [asdict(aggregate) for aggregate in self.aggregates],
        }

    def deterministic_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "task": self.task,
            "experiments": [
                {
                    "method": experiment.method,
                    "runs": [run.deterministic_dict() for run in experiment.runs],
                }
                for experiment in self.experiments
            ],
            "aggregates": [
                aggregate.deterministic_dict() for aggregate in self.aggregates
            ],
        }

    def save(self, path: Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def _logged_state(index: int, node: RRTNode) -> LoggedState:
    return LoggedState(
        node_index=index,
        parent_index=node.parent_index,
        tip_position=tuple(float(value) for value in node.tip_position),
        insertion=float(node.state.actuation.xi),
        field=tuple(float(value) for value in node.state.actuation.field),
        active_contact_ids=tuple(
            (int(contact[0]), int(contact[1]))
            for contact in node.state.active_contact_ids
        ),
        contact_history_boundaries=tuple(int(value) for value in node.contact_history),
        minimum_contact_multiplier=(
            min(float(value) for value in node.state.multipliers)
            if node.state.multipliers else None
        ),
        stability_margin=float(node.state.stability_margin),
        minimum_body_gap=float(node.state.minimum_body_gap),
        tip_clearance=float(node.state.tip_clearance),
        tip_safe=bool(node.state.tip_safe),
    )


def _path_indices(nodes: Sequence[RRTNode], final_index: int) -> List[int]:
    path = []
    index: Optional[int] = final_index
    while index is not None:
        path.append(index)
        index = nodes[index].parent_index
    path.reverse()
    return path


def _contact_event_satisfied(path: Sequence[LoggedState], required: str) -> bool:
    if required == "none":
        return True
    if required == "sustained":
        return bool(path) and all(state.active_contact_ids for state in path)
    inserted = False
    released = False
    for previous, current in zip(path, path[1:]):
        previous_ids = set(previous.active_contact_ids)
        current_ids = set(current.active_contact_ids)
        inserted = inserted or bool(current_ids - previous_ids)
        released = released or bool(previous_ids - current_ids)
    if required == "insertion":
        return inserted
    if required == "release":
        return released
    return inserted or released


def _terminal_failure_reason(diagnostic: object) -> str:
    if diagnostic.success:
        return ""
    if not diagnostic.steps:
        return "no_steps"
    step = diagnostic.steps[-1]
    if step.continuation_failed:
        return f"continuation:{step.failure_reason or 'unknown'}"
    if step.jacobian_rank < 3:
        return "rank_deficient"
    if step.accepted:
        return "terminal_iteration_limit"
    return "no_progress"


def _history_label(boundaries: Sequence[int]) -> str:
    relevant = set(boundaries) & {1, 2}
    if relevant == {1}:
        return "boundary_1_only"
    if relevant == {2}:
        return "boundary_2_only"
    if relevant:
        return "mixed"
    return "none"


def run_task(
    task: PlanningTask,
    include_full_tree: bool = False,
    method: str = "mechanics_informed_rrt",
) -> ExperimentRecord:
    task.validated()
    if method not in PLANNING_METHODS:
        raise ValueError("Unsupported planning method")
    runs = []
    for seed in task.seeds:
        mechanics = simder.MechanicsSession(task.mechanics.to_binding())
        initial = mechanics.solve_initial_state()
        planner_config = replace(task.planner, random_seed=seed)
        if method == "local_goal_seeking":
            planner = LocalGoalSeekingPlanner(
                mechanics,
                planner_config,
                task.continuation.to_binding(),
            )
        else:
            planner = TipSpaceRRT(
                mechanics,
                planner_config,
                task.continuation.to_binding(),
                proposal_policy=(
                    "mechanics_informed"
                    if method == "mechanics_informed_rrt"
                    else "random_actuation"
                ),
            )
        start_time = time.perf_counter()
        result = planner.plan(initial, task.target_tip)
        elapsed = time.perf_counter() - start_time

        final_index = result.goal_index
        if final_index is None:
            final_index = min(
                range(len(result.nodes)),
                key=lambda index: _distance(
                    result.nodes[index].tip_position, task.target_tip
                ),
            )
        path_indices = _path_indices(result.nodes, final_index)
        path = tuple(_logged_state(index, result.nodes[index]) for index in path_indices)
        tree = (
            tuple(_logged_state(index, node) for index, node in enumerate(result.nodes))
            if include_full_tree
            else None
        )
        terminal_steps = [
            step
            for diagnostic in result.terminal_diagnostics
            for step in diagnostic.steps
        ]
        accepted_edges = sum(
            diagnostic.accepted for diagnostic in result.edge_diagnostics
        ) + sum(step.accepted for step in terminal_steps)
        attempted_edges = len(result.edge_diagnostics) + len(terminal_steps)
        failed_continuations = sum(
            diagnostic.continuation_failed
            for diagnostic in result.edge_diagnostics
        ) + sum(step.continuation_failed for step in terminal_steps)
        tip_safety_failures = sum(
            diagnostic.failure_reason == "tip_safety"
            for diagnostic in result.edge_diagnostics
        ) + sum(
            step.failure_reason == "tip_safety" for step in terminal_steps
        )
        continuation_attempted_steps = sum(
            diagnostic.attempted_steps
            for diagnostic in result.edge_diagnostics
        ) + sum(step.attempted_steps for step in terminal_steps)
        continuation_rejected_steps = sum(
            diagnostic.rejected_steps
            for diagnostic in result.edge_diagnostics
        ) + sum(step.rejected_steps for step in terminal_steps)
        if continuation_attempted_steps != result.continuation_attempted_steps:
            raise RuntimeError("Planner and experiment continuation work disagree")
        contact_transitions = sum(
            diagnostic.contacts_added + diagnostic.contacts_released
            for diagnostic in result.edge_diagnostics
        ) + sum(
            step.contacts_added + step.contacts_released
            for step in terminal_steps
        )
        peak_field = max(_norm(state.field) for state in path)
        minimum_stability = min(state.stability_margin for state in path)
        path_contact_transitions = sum(
            path[index].active_contact_ids != path[index - 1].active_contact_ids
            for index in range(1, len(path))
        )
        contact_multipliers = [
            state.minimum_contact_multiplier
            for state in path
            if state.minimum_contact_multiplier is not None
        ]
        tree_history_counts = Counter(
            _history_label(node.contact_history) for node in result.nodes
        )
        minimum_history_distance: Dict[str, float] = {}
        for node in result.nodes:
            label = _history_label(node.contact_history)
            distance = _distance(node.tip_position, task.target_tip)
            minimum_history_distance[label] = min(
                distance, minimum_history_distance.get(label, math.inf)
            )
        runs.append(RunRecord(
            method=method,
            seed=seed,
            success=(
                result.success
                and _contact_event_satisfied(path, task.required_contact_event)
            ),
            termination_reason=result.termination_reason,
            continuation_step_budget=planner_config.maximum_continuation_steps,
            continuation_step_budget_exhausted=(
                result.termination_reason == "mechanics_work_budget"
            ),
            final_tip_error=_distance(
                result.nodes[final_index].tip_position, task.target_tip
            ),
            rrt_iterations=result.iterations,
            attempted_edges=attempted_edges,
            accepted_edges=accepted_edges,
            accepted_node_count=len(result.nodes),
            failed_continuations=failed_continuations,
            tip_safety_failures=tip_safety_failures,
            continuation_attempted_steps=continuation_attempted_steps,
            continuation_rejected_steps=continuation_rejected_steps,
            mean_continuation_steps_per_accepted_edge=(
                continuation_attempted_steps / accepted_edges
                if accepted_edges > 0 else 0.0
            ),
            terminal_attempts=len(result.terminal_diagnostics),
            terminal_failures=sum(
                not diagnostic.success
                for diagnostic in result.terminal_diagnostics
            ),
            terminal_steps=len(terminal_steps),
            terminal_failure_reasons=tuple(
                _terminal_failure_reason(diagnostic)
                for diagnostic in result.terminal_diagnostics
                if not diagnostic.success
            ),
            minimum_terminal_projected_residual=(
                min(
                    step.projected_linear_residual_norm
                    for step in terminal_steps
                    if math.isfinite(step.projected_linear_residual_norm)
                )
                if any(
                    math.isfinite(step.projected_linear_residual_norm)
                    for step in terminal_steps
                )
                else None
            ),
            final_terminal_projected_residual=(
                terminal_steps[-1].projected_linear_residual_norm
                if terminal_steps
                and math.isfinite(
                    terminal_steps[-1].projected_linear_residual_norm
                )
                else None
            ),
            final_terminal_projected_residual_ratio=(
                terminal_steps[-1].projected_linear_residual_ratio
                if terminal_steps
                and math.isfinite(
                    terminal_steps[-1].projected_linear_residual_ratio
                )
                else None
            ),
            terminal_insertion_constraint_activations=sum(
                step.insertion_constraint_active for step in terminal_steps
            ),
            nearest_neighbor_queries=result.nearest_neighbor_queries,
            nearest_neighbor_distance_evaluations=(
                result.nearest_neighbor_distance_evaluations
            ),
            shortlist_candidates_evaluated=(
                result.shortlist_candidates_evaluated
            ),
            shortlist_non_nearest_selections=(
                result.shortlist_non_nearest_selections
            ),
            shortlist_triggered_queries=result.shortlist_triggered_queries,
            tree_contact_history_counts=tuple(sorted(tree_history_counts.items())),
            minimum_goal_distance_by_contact_history=tuple(
                sorted(minimum_history_distance.items())
            ),
            contact_transitions=contact_transitions,
            path_contact_transitions=path_contact_transitions,
            peak_field=peak_field,
            total_insertion=path[-1].insertion - path[0].insertion,
            minimum_stability_margin=minimum_stability,
            minimum_body_gap=min(state.minimum_body_gap for state in path),
            minimum_tip_clearance=min(state.tip_clearance for state in path),
            minimum_contact_multiplier=(
                min(contact_multipliers) if contact_multipliers else None
            ),
            planning_time_seconds=elapsed,
            path=path,
            tree=tree,
        ))
    return ExperimentRecord(task=task.to_dict(), method=method, runs=tuple(runs))


def run_comparison(
    task: PlanningTask,
    methods: Sequence[str] = PLANNING_METHODS,
    include_full_tree: bool = False,
) -> ComparisonRecord:
    if not methods or len(set(methods)) != len(methods):
        raise ValueError("Comparison methods must be nonempty and unique")
    experiments = tuple(
        run_task(task, include_full_tree=include_full_tree, method=method)
        for method in methods
    )
    aggregates = []
    for experiment in experiments:
        runs = experiment.runs
        aggregates.append(MethodAggregate(
            method=experiment.method,
            run_count=len(runs),
            success_rate=sum(run.success for run in runs) / len(runs),
            median_final_tip_error=statistics.median(
                run.final_tip_error for run in runs
            ),
            median_planning_time_seconds=statistics.median(
                run.planning_time_seconds for run in runs
            ),
            median_continuation_attempted_steps=statistics.median(
                run.continuation_attempted_steps for run in runs
            ),
        ))
    return ComparisonRecord(
        task=task.to_dict(),
        experiments=experiments,
        aggregates=tuple(aggregates),
    )
