import json
import math
from pathlib import Path
import tempfile

import simder
from simder_planning import (
    PlanningTask,
    double_obstacle_contact_branching_task,
)


task = double_obstacle_contact_branching_task((0, 1, 2))
assert task.mechanics.domain_type == "double_spherical_obstacle"
assert task.required_contact_event == "any"
session = simder.MechanicsSession(task.mechanics.to_binding())
start = session.solve_initial_state()
assert not start.active_contact_ids
assert start.minimum_body_gap > 0.0

# The direct actuation is mechanically valid but becomes trapped in the narrow
# gate. It converges to a two-contact branch far from the calibrated goal.
target_actuation = simder.Actuation(5.0e-3, [0.0, 0.0, 0.3])
direct = session.attempt_continuation(
    start, target_actuation, task.continuation.to_binding()
)
assert direct.success
assert set(map(tuple, direct.state.active_contact_ids)) == {
    (14, 1),
    (14, 2),
}
assert math.dist(direct.state.tip_position, task.target_tip) > 1.0e-2

# The two routes select different contact histories, move around opposite
# sides of the gate, and converge to the same contact-free goal branch.
final_states = []
for sign, expected_contact in ((1.0, (14, 1)), (-1.0, (14, 2))):
    current = start
    history = []
    route = (
        simder.Actuation(5.0e-3, [0.0, sign * 0.05, 0.0]),
        simder.Actuation(5.0e-3, [0.0, sign * 0.8, 0.0]),
        simder.Actuation(5.0e-3, [0.0, sign * 0.8, 0.3]),
        target_actuation,
    )
    for waypoint in route:
        edge = session.attempt_continuation(
            current, waypoint, task.continuation.to_binding()
        )
        assert edge.success
        assert edge.state.tip_safe
        assert edge.state.minimum_body_gap >= -1.0e-8
        history.extend(map(tuple, edge.state.active_contact_ids))
        current = edge.state
    assert history[0] == expected_contact
    assert expected_contact in history
    assert not current.active_contact_ids
    assert math.dist(current.tip_position, task.target_tip) < 1.1e-8
    final_states.append(current)

assert math.dist(
    final_states[0].tip_position, final_states[1].tip_position
) < 1.1e-8

with tempfile.TemporaryDirectory() as temporary_directory:
    path = Path(temporary_directory) / "contact_branching_task.json"
    task.save(path)
    stored = json.loads(path.read_text(encoding="utf-8"))
    assert stored["mechanics"]["second_obstacle_radius"] == 1.0e-3
    loaded = PlanningTask.load(path)
    assert loaded.to_dict() == task.to_dict()

print(
    "Contact branching mechanics: direct_pinned=1, boundary1_route=1, "
    "boundary2_route=1, common_goal=1"
)
