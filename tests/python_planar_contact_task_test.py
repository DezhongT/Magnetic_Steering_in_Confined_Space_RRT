import json
import math
from pathlib import Path
import tempfile

import simder
from simder_planning import (
    COMPARISON_SCHEMA_VERSION,
    PlanningTask,
    planar_contact_release_task,
    run_comparison,
)


task = planar_contact_release_task((0, 1, 2))
assert task.required_contact_event == "release"

# Validate the calibrated reference path independently of the planners.
reference_session = simder.MechanicsSession(task.mechanics.to_binding())
reference_start = reference_session.solve_initial_state()
assert reference_start.active_contact_ids
assert reference_start.multipliers
assert min(reference_start.multipliers) >= 0.0
assert math.isfinite(reference_start.minimum_body_gap)
assert reference_start.tip_clearance > 0.0
reference = reference_session.attempt_continuation(
    reference_start,
    simder.Actuation(reference_start.actuation.xi, [0.0, 0.1, -0.25]),
    task.continuation.to_binding(),
)
assert reference.success
assert reference.contacts_released >= 1
assert not reference.state.active_contact_ids
assert math.dist(reference.state.tip_position, task.target_tip) < 1.0e-14

with tempfile.TemporaryDirectory() as temporary_directory:
    directory = Path(temporary_directory)
    task_path = directory / "planar_contact_task.json"
    comparison_path = directory / "planar_contact_comparison.json"
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
    assert sum(run.success for run in by_method["local_goal_seeking"].runs) == 3

    for experiment in first.experiments:
        assert experiment.task == task.to_dict()
        for run in experiment.runs:
            assert run.path[0].active_contact_ids
            assert not run.path[-1].active_contact_ids
            assert run.path_contact_transitions >= 1
            assert run.minimum_contact_multiplier is not None
            assert run.minimum_contact_multiplier >= -1.0e-8
            assert run.minimum_body_gap >= -1.0e-8
            assert run.minimum_tip_clearance >= 0.0
            assert run.minimum_stability_margin > 0.0
            assert all(state.tip_safe for state in run.path)
            assert all(
                state.minimum_contact_multiplier is None
                or state.minimum_contact_multiplier >= -1.0e-8
                for state in run.path
            )

    replay = run_comparison(PlanningTask.load(task_path))
    assert replay.deterministic_dict() == first.deterministic_dict()

print(
    "Planar contact task: reference_release=1, mechanics=3/3, "
    "random=3/3, local=3/3, replay_exact=1"
)
