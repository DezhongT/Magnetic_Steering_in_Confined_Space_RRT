import json
from pathlib import Path
import tempfile

from simder_planning import (
    COMPARISON_SCHEMA_VERSION,
    PLANNING_METHODS,
    contact_free_task_a,
    run_comparison,
)


task = contact_free_task_a((0, 1, 2))
first = run_comparison(task)
assert tuple(experiment.method for experiment in first.experiments) == PLANNING_METHODS
assert all(experiment.task == task.to_dict() for experiment in first.experiments)
assert all(
    run.method == experiment.method
    for experiment in first.experiments
    for run in experiment.runs
)

by_method = {
    experiment.method: experiment
    for experiment in first.experiments
}
assert all(run.success for run in by_method["mechanics_informed_rrt"].runs)
assert not any(run.success for run in by_method["random_actuation_rrt"].runs)
assert all(run.success for run in by_method["local_goal_seeking"].runs)
assert all(
    run.final_tip_error <= task.planner.terminal_tolerance
    for method in ("mechanics_informed_rrt", "local_goal_seeking")
    for run in by_method[method].runs
)

for experiment in first.experiments:
    for run in experiment.runs:
        assert run.total_insertion >= 0.0
        assert run.peak_field <= task.planner.maximum_field_magnitude
        assert run.minimum_stability_margin > 0.0
        assert run.contact_transitions == 0
        assert all(
            state.insertion >= run.path[index - 1].insertion
            for index, state in enumerate(run.path)
            if index > 0
        )

local_paths = [
    [state.tip_position for state in run.path]
    for run in by_method["local_goal_seeking"].runs
]
assert local_paths[0] == local_paths[1] == local_paths[2]

replay = run_comparison(task)
assert replay.deterministic_dict() == first.deterministic_dict()

with tempfile.TemporaryDirectory() as temporary_directory:
    output = Path(temporary_directory) / "comparison.json"
    first.save(output)
    stored = json.loads(output.read_text(encoding="utf-8"))
    assert stored["schema_version"] == COMPARISON_SCHEMA_VERSION
    assert len(stored["experiments"]) == 3
    assert len(stored["aggregates"]) == 3

aggregates = {aggregate.method: aggregate for aggregate in first.aggregates}
assert aggregates["mechanics_informed_rrt"].success_rate == 1.0
assert aggregates["random_actuation_rrt"].success_rate == 0.0
assert aggregates["local_goal_seeking"].success_rate == 1.0

print(
    "Task A comparison: mechanics=3/3, random=0/3, local=3/3, "
    "replay_exact=1"
)
