import json
import math
from pathlib import Path
import tempfile

import simder
from simder_planning import (
    COMPARISON_SCHEMA_VERSION,
    PlanningTask,
    contact_free_task_a,
    run_comparison,
    spherical_shell_task_b,
)


task = spherical_shell_task_b((0, 1, 2))
assert task.mechanics.domain_type == "spherical_shell"
assert task.required_contact_event == "sustained"

# Calibrate the target independently using a known fixed-active continuation.
reference_session = simder.MechanicsSession(task.mechanics.to_binding())
reference_start = reference_session.solve_initial_state()
assert reference_start.active_contact_ids == [[14, 0]]
assert reference_start.multipliers[0] > 0.0
reference = reference_session.attempt_continuation(
    reference_start,
    simder.Actuation(5.0e-3, [0.1, 0.2, -0.25]),
    task.continuation.to_binding(),
)
assert reference.success
assert reference.contacts_added == 0
assert reference.contacts_released == 0
assert reference.state.active_contact_ids == [[14, 0]]
assert math.dist(reference.state.tip_position, task.target_tip) < 1.0e-14

# A v1 task written before domain fields existed must retain planar defaults.
legacy_values = contact_free_task_a((0,)).to_dict()
for name in (
    "domain_type",
    "shell_center",
    "shell_radius",
    "shell_minus_thickness",
    "shell_plus_thickness",
    "cavity_center",
    "cavity_radius",
    "obstacle_center",
    "obstacle_radius",
):
    legacy_values["mechanics"].pop(name)
legacy_values.pop("required_contact_event")
legacy_task = PlanningTask.from_dict(legacy_values)
assert legacy_task.mechanics.domain_type == "planar_slab"
assert legacy_task.required_contact_event == "none"

with tempfile.TemporaryDirectory() as temporary_directory:
    directory = Path(temporary_directory)
    task_path = directory / "spherical_shell_task.json"
    comparison_path = directory / "spherical_shell_comparison.json"
    task.save(task_path)
    loaded = PlanningTask.load(task_path)
    assert loaded.to_dict() == task.to_dict()

    first = run_comparison(loaded)
    first.save(comparison_path)
    stored = json.loads(comparison_path.read_text(encoding="utf-8"))
    assert stored["schema_version"] == COMPARISON_SCHEMA_VERSION

    for experiment in first.experiments:
        assert experiment.task == task.to_dict()
        assert all(run.success for run in experiment.runs)
        for run in experiment.runs:
            assert run.path_contact_transitions == 0
            assert run.minimum_contact_multiplier is not None
            assert run.minimum_contact_multiplier >= -1.0e-8
            assert run.minimum_body_gap >= -1.0e-8
            assert run.minimum_tip_clearance > 0.0
            assert run.minimum_stability_margin > 0.0
            assert all(state.tip_safe for state in run.path)
            assert all(state.active_contact_ids == ((14, 0),) for state in run.path)

    by_method = {
        aggregate.method: aggregate
        for aggregate in first.aggregates
    }
    assert by_method["mechanics_informed_rrt"].success_rate == 1.0
    assert by_method["random_actuation_rrt"].success_rate == 1.0
    assert by_method["local_goal_seeking"].success_rate == 1.0

    replay = run_comparison(PlanningTask.load(task_path))
    assert replay.deterministic_dict() == first.deterministic_dict()

print(
    "Spherical-shell Task B: sustained_contact=1, mechanics=3/3, "
    "random=3/3, local=3/3, replay_exact=1"
)
