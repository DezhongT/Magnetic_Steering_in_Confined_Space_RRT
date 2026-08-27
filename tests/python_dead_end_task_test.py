import json
import math
from pathlib import Path
import tempfile

import simder
from simder_planning import (
    COMPARISON_SCHEMA_VERSION,
    PlanningTask,
    run_comparison,
    spherical_obstacle_dead_end_task,
)


task = spherical_obstacle_dead_end_task((0, 1, 2))
assert task.mechanics.domain_type == "spherical_obstacle"
session = simder.MechanicsSession(task.mechanics.to_binding())
start = session.solve_initial_state()
start_configuration = list(start.configuration)

# The direct actuation path is a real dead end: it intersects the forbidden
# tip-clearance region and rolls back after making only partial progress.
target_actuation = simder.Actuation(5.0e-3, [0.1, 0.25, 0.1])
direct = session.attempt_continuation(
    start, target_actuation, task.continuation.to_binding()
)
assert not direct.success
assert direct.rolled_back
assert direct.failure_reason == "tip_safety"
assert 0.0 < direct.reached_path_fraction < 0.5
assert direct.state.configuration == start_configuration

# A deterministic route first moves below the obstacle, advances laterally,
# and only then approaches the target from the safe side.
route = (
    simder.Actuation(1.0e-3, [0.0, 0.05, -0.2]),
    simder.Actuation(5.0e-3, [0.1, 0.25, -0.2]),
    target_actuation,
)
current = start
for waypoint in route:
    edge = session.attempt_continuation(
        current, waypoint, task.continuation.to_binding()
    )
    assert edge.success
    assert edge.state.tip_safe
    assert edge.state.tip_clearance >= task.mechanics.tip_safe_distance
    current = edge.state
assert math.dist(current.tip_position, task.target_tip) < 1.0e-14

with tempfile.TemporaryDirectory() as temporary_directory:
    directory = Path(temporary_directory)
    task_path = directory / "dead_end_task.json"
    comparison_path = directory / "dead_end_comparison.json"
    task.save(task_path)
    loaded = PlanningTask.load(task_path)
    assert loaded.to_dict() == task.to_dict()

    first = run_comparison(loaded)
    first.save(comparison_path)
    stored = json.loads(comparison_path.read_text(encoding="utf-8"))
    assert stored["schema_version"] == COMPARISON_SCHEMA_VERSION
    by_method = {
        experiment.method: experiment
        for experiment in first.experiments
    }
    assert sum(run.success for run in by_method["mechanics_informed_rrt"].runs) == 3
    assert sum(run.success for run in by_method["random_actuation_rrt"].runs) == 3
    assert sum(run.success for run in by_method["local_goal_seeking"].runs) == 0
    assert all(
        run.tip_safety_failures >= 1
        for run in by_method["local_goal_seeking"].runs
    )

    obstacle_center = task.mechanics.obstacle_center
    required_tip_distance = (
        task.mechanics.obstacle_radius + task.mechanics.tip_safe_distance
    )
    for method in ("mechanics_informed_rrt", "random_actuation_rrt"):
        for run in by_method[method].runs:
            assert run.failed_continuations > 0
            assert run.minimum_body_gap > 0.0
            assert run.minimum_tip_clearance >= task.mechanics.tip_safe_distance
            assert run.minimum_stability_margin > 0.0
            assert all(not state.active_contact_ids for state in run.path)
            assert all(state.tip_safe for state in run.path)
            assert all(
                math.dist(state.tip_position, obstacle_center)
                >= required_tip_distance - 1.0e-12
                for state in run.path
            )

    aggregate = {value.method: value for value in first.aggregates}
    assert (
        aggregate["mechanics_informed_rrt"].median_continuation_attempted_steps
        < aggregate["random_actuation_rrt"].median_continuation_attempted_steps
    )

    replay = run_comparison(PlanningTask.load(task_path))
    assert replay.deterministic_dict() == first.deterministic_dict()

print(
    "Dead-end task: direct_blocked=1, reference_route=1, mechanics=3/3, "
    "random=3/3, local=0/3, replay_exact=1"
)
