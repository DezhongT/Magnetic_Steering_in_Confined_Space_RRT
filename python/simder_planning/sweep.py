from dataclasses import asdict, dataclass, replace
from collections import Counter
import json
import math
from pathlib import Path
import statistics
from typing import Any, Dict, Optional, Sequence, Tuple

from .experiment import ComparisonRecord, PLANNING_METHODS, run_comparison
from .task import (
    PlanningTask,
    double_obstacle_contact_branching_task,
    spherical_obstacle_dead_end_task,
)


SWEEP_SCHEMA_VERSION = "simder-planning-robustness-sweep-v1"


@dataclass(frozen=True)
class SweepCase:
    name: str
    parameter: str
    value: Optional[float]
    task: PlanningTask


@dataclass(frozen=True)
class SweepMethodSummary:
    method: str
    run_count: int
    success_count: int
    success_rate: float
    success_rate_ci95_low: float
    success_rate_ci95_high: float
    median_final_tip_error: float
    median_successful_tip_error: Optional[float]
    median_continuation_attempted_steps: float
    median_failed_continuations: float
    median_tip_safety_failures: float
    budget_exhaustion_count: int
    iteration_exhaustion_count: int
    termination_reason_counts: Tuple[Tuple[str, int], ...]
    terminal_failure_reason_counts: Tuple[Tuple[str, int], ...]
    successful_contact_history_counts: Tuple[Tuple[str, int], ...]


@dataclass(frozen=True)
class SweepCaseRecord:
    name: str
    parameter: str
    value: Optional[float]
    comparison: ComparisonRecord
    summaries: Tuple[SweepMethodSummary, ...]


@dataclass(frozen=True)
class RobustnessSweepRecord:
    cases: Tuple[SweepCaseRecord, ...]
    schema_version: str = SWEEP_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "cases": [
                {
                    "name": case.name,
                    "parameter": case.parameter,
                    "value": case.value,
                    "comparison": case.comparison.to_dict(),
                    "summaries": [asdict(summary) for summary in case.summaries],
                }
                for case in self.cases
            ],
        }

    def deterministic_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "cases": [
                {
                    "name": case.name,
                    "parameter": case.parameter,
                    "value": case.value,
                    "comparison": case.comparison.deterministic_dict(),
                    "summaries": [asdict(summary) for summary in case.summaries],
                }
                for case in self.cases
            ],
        }

    def save(self, path: Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def wilson_score_interval(
    success_count: int,
    run_count: int,
    z_score: float = 1.959963984540054,
) -> Tuple[float, float]:
    if run_count <= 0:
        raise ValueError("Wilson interval requires at least one run")
    if success_count < 0 or success_count > run_count:
        raise ValueError("Success count must be between zero and run count")
    probability = success_count / run_count
    z_squared = z_score * z_score
    denominator = 1.0 + z_squared / run_count
    center = (probability + z_squared / (2.0 * run_count)) / denominator
    half_width = z_score * math.sqrt(
        probability * (1.0 - probability) / run_count
        + z_squared / (4.0 * run_count * run_count)
    ) / denominator
    return center - half_width, center + half_width


def dead_end_sweep_cases(
    seeds: Sequence[int] = tuple(range(10)),
    profile: str = "full",
) -> Tuple[SweepCase, ...]:
    if profile not in ("full", "regression"):
        raise ValueError("Sweep profile must be 'full' or 'regression'")
    baseline = spherical_obstacle_dead_end_task(seeds)
    center = baseline.mechanics.obstacle_center

    definitions = (
        ("baseline", "baseline", None, baseline.mechanics, baseline.planner),
        ("obstacle_radius_low", "obstacle_radius", 2.0e-4,
         replace(baseline.mechanics, obstacle_radius=2.0e-4), baseline.planner),
        ("obstacle_radius_high", "obstacle_radius", 3.0e-4,
         replace(baseline.mechanics, obstacle_radius=3.0e-4), baseline.planner),
        ("obstacle_offset_z_low", "obstacle_center_z", center[2] - 1.5e-4,
         replace(baseline.mechanics,
                 obstacle_center=(center[0], center[1], center[2] - 1.5e-4)),
         baseline.planner),
        ("obstacle_offset_z_high", "obstacle_center_z", center[2] + 1.5e-4,
         replace(baseline.mechanics,
                 obstacle_center=(center[0], center[1], center[2] + 1.5e-4)),
         baseline.planner),
        ("tip_safety_low", "tip_safe_distance", 5.0e-5,
         replace(baseline.mechanics, tip_safe_distance=5.0e-5), baseline.planner),
        ("tip_safety_high", "tip_safe_distance", 1.5e-4,
         replace(baseline.mechanics, tip_safe_distance=1.5e-4), baseline.planner),
        ("field_authority_low", "maximum_field_magnitude", 0.3,
         baseline.mechanics, replace(baseline.planner, maximum_field_magnitude=0.3)),
        ("field_authority_high", "maximum_field_magnitude", 0.7,
         baseline.mechanics, replace(baseline.planner, maximum_field_magnitude=0.7)),
        ("goal_bias_low", "goal_bias", 0.05,
         baseline.mechanics, replace(baseline.planner, goal_bias=0.05)),
        ("goal_bias_high", "goal_bias", 0.5,
         baseline.mechanics, replace(baseline.planner, goal_bias=0.5)),
    )

    if profile == "regression":
        selected = {
            "baseline",
            "obstacle_radius_high",
            "obstacle_offset_z_high",
            "tip_safety_high",
            "field_authority_low",
            "goal_bias_low",
        }
        definitions = tuple(value for value in definitions if value[0] in selected)

    return tuple(
        SweepCase(
            name=name,
            parameter=parameter,
            value=value,
            task=replace(
                baseline,
                name=f"{baseline.name}_{name}",
                description=f"{baseline.description} Robustness case: {name}.",
                mechanics=mechanics,
                planner=planner,
            ).validated(),
        )
        for name, parameter, value, mechanics, planner in definitions
    )


def contact_branching_budget_cases(
    seeds: Sequence[int] = tuple(range(10)),
    profile: str = "full",
    use_contact_history_shortlist: bool = False,
    contact_history_shortlist_start_iteration: Optional[int] = None,
    contact_history_shortlist_period: Optional[int] = None,
    contact_history_shortlist_stagnation_iterations: Optional[int] = None,
    contact_history_shortlist_stagnation_progress: float = 1.0e-4,
    contact_history_shortlist_burst_iterations: int = 1,
    use_dual_frontier: bool = False,
    dual_frontier_auxiliary_period: int = 4,
) -> Tuple[SweepCase, ...]:
    if profile not in ("full", "regression"):
        raise ValueError("Sweep profile must be 'full' or 'regression'")
    baseline = double_obstacle_contact_branching_task(seeds)
    if (
        use_contact_history_shortlist
        or contact_history_shortlist_start_iteration is not None
        or contact_history_shortlist_period is not None
        or contact_history_shortlist_stagnation_iterations is not None
        or use_dual_frontier
    ):
        baseline = replace(
            baseline,
            description=(
                f"{baseline.description} Uses the contact-history-stratified "
                "mechanics shortlist."
            ),
            planner=replace(
                baseline.planner,
                use_contact_history_shortlist=use_contact_history_shortlist,
                contact_history_shortlist_start_iteration=(
                    contact_history_shortlist_start_iteration
                ),
                contact_history_shortlist_period=(
                    contact_history_shortlist_period
                ),
                contact_history_shortlist_stagnation_iterations=(
                    contact_history_shortlist_stagnation_iterations
                ),
                contact_history_shortlist_stagnation_progress=(
                    contact_history_shortlist_stagnation_progress
                ),
                contact_history_shortlist_burst_iterations=(
                    contact_history_shortlist_burst_iterations
                ),
                use_dual_frontier=use_dual_frontier,
                dual_frontier_auxiliary_period=dual_frontier_auxiliary_period,
            ),
        )
    definitions = (
        ("iterations_500", "maximum_iterations", 500,
         replace(baseline.planner, maximum_iterations=500,
                 maximum_continuation_steps=15000)),
        ("iterations_1000", "maximum_iterations", 1000,
         replace(baseline.planner, maximum_iterations=1000,
                 maximum_continuation_steps=15000)),
        ("work_3000", "maximum_continuation_steps", 3000,
         replace(baseline.planner, maximum_iterations=2000,
                 maximum_continuation_steps=3000)),
        ("work_7000", "maximum_continuation_steps", 7000,
         replace(baseline.planner, maximum_iterations=2000,
                 maximum_continuation_steps=7000)),
        ("work_15000", "maximum_continuation_steps", 15000,
         replace(baseline.planner, maximum_iterations=2000,
                 maximum_continuation_steps=15000)),
    )
    if profile == "regression":
        definitions = definitions[2:3]
    return tuple(
        SweepCase(
            name=name,
            parameter=parameter,
            value=float(value),
            task=replace(
                baseline,
                name=f"{baseline.name}_{name}",
                description=f"{baseline.description} Budget case: {name}.",
                planner=planner,
            ).validated(),
        )
        for name, parameter, value, planner in definitions
    )


def _contact_history(run: object) -> str:
    boundaries = {
        boundary
        for state in run.path
        for boundary in state.contact_history_boundaries
        if boundary in (1, 2)
    }
    if boundaries == {1}:
        return "boundary_1_only"
    if boundaries == {2}:
        return "boundary_2_only"
    if boundaries:
        return "mixed"
    return "none"


def _summaries(comparison: ComparisonRecord) -> Tuple[SweepMethodSummary, ...]:
    summaries = []
    for experiment in comparison.experiments:
        runs = experiment.runs
        success_count = sum(run.success for run in runs)
        ci_low, ci_high = wilson_score_interval(success_count, len(runs))
        successful_errors = [run.final_tip_error for run in runs if run.success]
        termination_reasons = Counter(run.termination_reason for run in runs)
        terminal_failure_reasons = Counter(
            reason for run in runs for reason in run.terminal_failure_reasons
        )
        successful_histories = Counter(
            _contact_history(run) for run in runs if run.success
        )
        summaries.append(SweepMethodSummary(
            method=experiment.method,
            run_count=len(runs),
            success_count=success_count,
            success_rate=success_count / len(runs),
            success_rate_ci95_low=ci_low,
            success_rate_ci95_high=ci_high,
            median_final_tip_error=statistics.median(
                run.final_tip_error for run in runs
            ),
            median_successful_tip_error=(
                statistics.median(successful_errors) if successful_errors else None
            ),
            median_continuation_attempted_steps=statistics.median(
                run.continuation_attempted_steps for run in runs
            ),
            median_failed_continuations=statistics.median(
                run.failed_continuations for run in runs
            ),
            median_tip_safety_failures=statistics.median(
                run.tip_safety_failures for run in runs
            ),
            budget_exhaustion_count=sum(
                run.continuation_step_budget_exhausted for run in runs
            ),
            iteration_exhaustion_count=sum(
                run.termination_reason == "iteration_budget" for run in runs
            ),
            termination_reason_counts=tuple(sorted(termination_reasons.items())),
            terminal_failure_reason_counts=tuple(
                sorted(terminal_failure_reasons.items())
            ),
            successful_contact_history_counts=tuple(
                sorted(successful_histories.items())
            ),
        ))
    return tuple(summaries)


def run_dead_end_sweep(
    seeds: Sequence[int] = tuple(range(10)),
    profile: str = "full",
    methods: Sequence[str] = PLANNING_METHODS,
) -> RobustnessSweepRecord:
    return run_sweep(dead_end_sweep_cases(seeds, profile=profile), methods)


def run_sweep(
    cases: Sequence[SweepCase],
    methods: Sequence[str] = PLANNING_METHODS,
) -> RobustnessSweepRecord:
    records = []
    for case in cases:
        comparison = run_comparison(case.task, methods=methods)
        records.append(SweepCaseRecord(
            name=case.name,
            parameter=case.parameter,
            value=case.value,
            comparison=comparison,
            summaries=_summaries(comparison),
        ))
    return RobustnessSweepRecord(cases=tuple(records))


def run_contact_branching_budget_sweep(
    seeds: Sequence[int] = tuple(range(10)),
    profile: str = "full",
    methods: Sequence[str] = ("mechanics_informed_rrt",),
    use_contact_history_shortlist: bool = False,
    contact_history_shortlist_start_iteration: Optional[int] = None,
    contact_history_shortlist_period: Optional[int] = None,
    contact_history_shortlist_stagnation_iterations: Optional[int] = None,
    contact_history_shortlist_stagnation_progress: float = 1.0e-4,
    contact_history_shortlist_burst_iterations: int = 1,
    use_dual_frontier: bool = False,
    dual_frontier_auxiliary_period: int = 4,
) -> RobustnessSweepRecord:
    return run_sweep(
        contact_branching_budget_cases(
            seeds,
            profile=profile,
            use_contact_history_shortlist=use_contact_history_shortlist,
            contact_history_shortlist_start_iteration=(
                contact_history_shortlist_start_iteration
            ),
            contact_history_shortlist_period=contact_history_shortlist_period,
            contact_history_shortlist_stagnation_iterations=(
                contact_history_shortlist_stagnation_iterations
            ),
            contact_history_shortlist_stagnation_progress=(
                contact_history_shortlist_stagnation_progress
            ),
            contact_history_shortlist_burst_iterations=(
                contact_history_shortlist_burst_iterations
            ),
            use_dual_frontier=use_dual_frontier,
            dual_frontier_auxiliary_period=dual_frontier_auxiliary_period,
        ),
        methods,
    )
