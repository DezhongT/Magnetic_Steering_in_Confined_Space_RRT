import math

import simder
from simder_planning import RRTConfig, TipSpaceRRT


def make_session():
    config = simder.MechanicsConfig()
    config.gravity = [0.0, 0.0, 0.22]
    config.initial_actuation = simder.Actuation(1.0e-3, [0.0, 0.1, 0.0])
    return config, simder.MechanicsSession(config)


mechanics_config, goal_session = make_session()
goal_start = goal_session.solve_initial_state()
goal_edge = goal_session.attempt_continuation(
    goal_start, simder.Actuation(1.2e-3, [0.01, 0.11, 0.005])
)
assert goal_edge.success
goal = goal_edge.state.tip_position

success_config = RRTConfig(
    workspace_min=(-0.05, -0.15, -0.099),
    workspace_max=(1.05, 0.15, 0.099),
    goal_bias=1.0,
    local_tip_step=1.0e-3,
    goal_neighborhood_radius=1.0e-3,
    terminal_tolerance=1.0e-7,
    maximum_terminal_iterations=10,
    maximum_iterations=1,
    random_seed=11,
    maximum_insertion_step=2.0e-4,
    maximum_field_step=2.0e-2,
    maximum_field_magnitude=0.25,
)
success_session = simder.MechanicsSession(mechanics_config)
success_start = success_session.solve_initial_state()
success = TipSpaceRRT(success_session, success_config).plan(success_start, goal)
assert success.success
assert success.iterations == 0
assert success.terminal_diagnostics
successful_connector = success.terminal_diagnostics[0]
assert successful_connector.success
assert successful_connector.final_tip_error <= success_config.terminal_tolerance
assert successful_connector.steps
assert all(step.jacobian_rank == 3 for step in successful_connector.steps)
assert all(
    math.isfinite(step.jacobian_condition)
    for step in successful_connector.steps
)
assert len(success.path_indices()) == len(success.nodes)

restricted_config = RRTConfig(
    workspace_min=(-0.05, -0.15, -0.099),
    workspace_max=(1.05, 0.15, 0.099),
    goal_bias=1.0,
    local_tip_step=1.0e-3,
    goal_neighborhood_radius=1.0e-3,
    terminal_tolerance=1.0e-9,
    maximum_terminal_iterations=1,
    maximum_iterations=1,
    random_seed=11,
    maximum_insertion_step=1.0e-10,
    maximum_field_step=1.0e-10,
    maximum_field_magnitude=0.25,
    minimum_tip_progress=1.0e-12,
)
restricted_session = simder.MechanicsSession(mechanics_config)
restricted_start = restricted_session.solve_initial_state()
restricted = TipSpaceRRT(restricted_session, restricted_config).plan(
    restricted_start, goal
)
assert not restricted.success
assert restricted.goal_index is None
assert restricted.terminal_diagnostics
assert not restricted.terminal_diagnostics[0].success
assert restricted.terminal_diagnostics[0].source_index == 0
assert restricted.nodes[0].tip_position == tuple(restricted_start.tip_position)

print(
    "Terminal connector: "
    f"steps={len(successful_connector.steps)}, "
    f"final_error={successful_connector.final_tip_error}, "
    "restricted_failure_preserved=1"
)
