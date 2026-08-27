import json
from pathlib import Path
import tempfile

from simder_planning import (
    EXPERIMENT_SCHEMA_VERSION,
    PlanningTask,
    contact_free_task_a,
    run_task,
)


task = contact_free_task_a((0, 1, 2))
with tempfile.TemporaryDirectory() as temporary_directory:
    directory = Path(temporary_directory)
    task_path = directory / "task_a.json"
    result_path = directory / "task_a_results.json"
    task.save(task_path)
    loaded = PlanningTask.load(task_path)
    assert loaded.to_dict() == task.to_dict()

    first = run_task(loaded)
    first.save(result_path)
    stored = json.loads(result_path.read_text(encoding="utf-8"))
    assert stored["schema_version"] == EXPERIMENT_SCHEMA_VERSION
    assert stored["task"] == json.loads(task_path.read_text(encoding="utf-8"))
    assert len(first.runs) == 3
    assert all(run.success for run in first.runs)
    assert all(run.final_tip_error <= task.planner.terminal_tolerance for run in first.runs)
    assert all(run.total_insertion >= 0.0 for run in first.runs)
    assert all(run.peak_field <= task.planner.maximum_field_magnitude for run in first.runs)
    assert all(run.minimum_stability_margin > 0.0 for run in first.runs)
    assert all(run.contact_transitions == 0 for run in first.runs)
    assert all(run.tree is None for run in first.runs)
    assert all(
        not state.active_contact_ids
        for run in first.runs
        for state in run.path
    )

    replay = run_task(PlanningTask.load(task_path))
    assert [run.deterministic_dict() for run in replay.runs] == [
        run.deterministic_dict() for run in first.runs
    ]

    full_tree = run_task(contact_free_task_a((0,)), include_full_tree=True)
    assert full_tree.runs[0].tree is not None
    assert len(full_tree.runs[0].tree) >= len(full_tree.runs[0].path)

print(
    "Task A experiment: seeds=3, successes=3, "
    "contact_transitions=0, replay_exact=1"
)
