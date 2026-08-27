import json
import math
from pathlib import Path
import tempfile

from simder_planning import (
    SWEEP_SCHEMA_VERSION,
    dead_end_sweep_cases,
    run_dead_end_sweep,
    wilson_score_interval,
)


full_cases = dead_end_sweep_cases((0,), profile="full")
assert len(full_cases) == 11
assert len({case.name for case in full_cases}) == len(full_cases)
assert all(case.task.seeds == (0,) for case in full_cases)
assert all(case.task.validated() == case.task for case in full_cases)

regression_cases = dead_end_sweep_cases((0,), profile="regression")
assert len(regression_cases) == 6
assert {case.parameter for case in regression_cases} == {
    "baseline",
    "obstacle_radius",
    "obstacle_center_z",
    "tip_safe_distance",
    "maximum_field_magnitude",
    "goal_bias",
}

zero_low, zero_high = wilson_score_interval(0, 10)
all_low, all_high = wilson_score_interval(10, 10)
assert zero_low == 0.0
assert math.isclose(zero_high, 0.2775327998628892)
assert math.isclose(all_low, 0.7224672001371107)
assert math.isclose(all_high, 1.0)

# Exercise persistence and deterministic replay on a single real mechanics case;
# the factory assertions above cover every one-factor parameterization.
first = run_dead_end_sweep(
    seeds=(0,),
    profile="regression",
    methods=("mechanics_informed_rrt",),
)
assert len(first.cases) == 6
for case in first.cases:
    assert len(case.summaries) == 1
    summary = case.summaries[0]
    assert summary.run_count == 1
    assert summary.success_count in (0, 1)
    assert 0.0 <= summary.success_rate_ci95_low <= summary.success_rate
    assert summary.success_rate <= summary.success_rate_ci95_high <= 1.0

with tempfile.TemporaryDirectory() as temporary_directory:
    output = Path(temporary_directory) / "dead_end_sweep.json"
    first.save(output)
    stored = json.loads(output.read_text(encoding="utf-8"))
    assert stored["schema_version"] == SWEEP_SCHEMA_VERSION
    assert len(stored["cases"]) == 6

replay = run_dead_end_sweep(
    seeds=(0,),
    profile="regression",
    methods=("mechanics_informed_rrt",),
)
assert replay.deterministic_dict() == first.deterministic_dict()

print("Dead-end sweep: cases=11, regression_cases=6, replay_exact=1")
