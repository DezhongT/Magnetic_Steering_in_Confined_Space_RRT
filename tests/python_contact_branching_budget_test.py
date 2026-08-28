import json
from pathlib import Path
import tempfile

from simder_planning import (
    SWEEP_SCHEMA_VERSION,
    contact_branching_budget_cases,
    run_contact_branching_budget_sweep,
)


cases = contact_branching_budget_cases((0,), profile="full")
assert [case.name for case in cases] == [
    "iterations_500",
    "iterations_1000",
    "work_3000",
    "work_7000",
    "work_15000",
]
assert cases[0].task.planner.maximum_iterations == 500
assert cases[2].task.planner.maximum_continuation_steps == 3000
assert not cases[0].task.planner.use_contact_history_shortlist
shortlist_cases = contact_branching_budget_cases(
    (0,), profile="regression", use_contact_history_shortlist=True
)
assert shortlist_cases[0].task.planner.use_contact_history_shortlist
periodic_cases = contact_branching_budget_cases(
    (0,),
    profile="regression",
    contact_history_shortlist_period=10,
    contact_history_shortlist_burst_iterations=5,
)
assert periodic_cases[0].task.planner.contact_history_shortlist_period == 10
assert periodic_cases[0].task.planner.contact_history_shortlist_burst_iterations == 5
fallback_cases = contact_branching_budget_cases(
    (0,), profile="regression", contact_history_shortlist_start_iteration=500
)
assert fallback_cases[0].task.planner.contact_history_shortlist_start_iteration == 500
stagnation_cases = contact_branching_budget_cases(
    (0,),
    profile="regression",
    contact_history_shortlist_stagnation_iterations=50,
)
assert (
    stagnation_cases[0].task.planner.contact_history_shortlist_stagnation_iterations
    == 50
)
dual_cases = contact_branching_budget_cases(
    (0,),
    profile="regression",
    use_dual_frontier=True,
    dual_frontier_auxiliary_period=4,
)
assert dual_cases[0].task.planner.use_dual_frontier
assert dual_cases[0].task.planner.dual_frontier_auxiliary_period == 4

record = run_contact_branching_budget_sweep(
    seeds=(0,),
    profile="regression",
    methods=("mechanics_informed_rrt",),
)
assert len(record.cases) == 1
run = record.cases[0].comparison.experiments[0].runs[0]
summary = record.cases[0].summaries[0]
assert run.continuation_step_budget == 3000
assert run.termination_reason in ("goal_reached", "mechanics_work_budget")
assert summary.run_count == 1
assert dict(summary.termination_reason_counts)[run.termination_reason] == 1
assert summary.budget_exhaustion_count == int(
    run.termination_reason == "mechanics_work_budget"
)
assert sum(dict(run.tree_contact_history_counts).values()) == run.accepted_node_count
assert run.nearest_neighbor_queries > 0
assert run.nearest_neighbor_distance_evaluations > 0
assert run.shortlist_candidates_evaluated == 0
assert run.shortlist_non_nearest_selections == 0
assert run.shortlist_triggered_queries == 0
assert run.backbone_expansion_queries == run.nearest_neighbor_queries
assert run.auxiliary_expansion_queries == 0
assert run.backbone_node_count == run.accepted_node_count
assert run.auxiliary_node_count == 0
assert run.local_steering_evaluations > 0

dual_record = run_contact_branching_budget_sweep(
    seeds=(0,),
    profile="regression",
    methods=("mechanics_informed_rrt",),
    use_dual_frontier=True,
    dual_frontier_auxiliary_period=2,
)
dual_run = dual_record.cases[0].comparison.experiments[0].runs[0]
assert dual_run.backbone_expansion_queries > 0
assert dual_run.auxiliary_expansion_queries > 0
assert dual_run.shortlist_triggered_queries == dual_run.auxiliary_expansion_queries
assert dual_run.backbone_node_count + dual_run.auxiliary_node_count == (
    dual_run.accepted_node_count
)
assert dual_run.auxiliary_node_count > 0
assert dual_run.local_steering_evaluations > 0
assert {state.frontier for state in dual_run.path} <= {"backbone", "auxiliary"}

with tempfile.TemporaryDirectory() as temporary_directory:
    output = Path(temporary_directory) / "contact_budget.json"
    record.save(output)
    stored = json.loads(output.read_text(encoding="utf-8"))
    assert stored["schema_version"] == SWEEP_SCHEMA_VERSION
    assert stored["cases"][0]["summaries"][0]["run_count"] == 1

print(
    "Contact budget sweep: cases=5, soft_cap=3000, "
    f"termination={run.termination_reason}"
)
