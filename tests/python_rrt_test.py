import math
from dataclasses import replace

import simder
from simder_planning import RRTConfig, RRTNode, TipSpaceRRT


config = simder.MechanicsConfig()
config.gravity = [0.0, 0.0, 0.22]
config.initial_actuation = simder.Actuation(1.0e-3, [0.0, 0.1, 0.0])
session = simder.MechanicsSession(config)
start = session.solve_initial_state()

# Construct a known reachable goal through the same nonlinear mechanics. The
# RRT receives only the resulting tip target, not the actuation that made it.
known_target = simder.Actuation(2.0e-3, [0.04, 0.16, 0.02])
known_edge = session.attempt_continuation(start, known_target)
assert known_edge.success
goal = known_edge.state.tip_position

rrt_config = RRTConfig(
    workspace_min=(-0.05, -0.15, -0.099),
    workspace_max=(1.05, 0.15, 0.099),
    goal_bias=1.0,
    local_tip_step=1.0e-3,
    goal_neighborhood_radius=1.0e-4,
    terminal_tolerance=1.0e-7,
    maximum_terminal_iterations=8,
    maximum_iterations=40,
    random_seed=7,
    maximum_insertion_step=2.0e-4,
    maximum_field_step=1.0e-2,
    maximum_field_magnitude=0.25,
)
planner = TipSpaceRRT(session, rrt_config)
result = planner.plan(start, goal)

assert result.success
assert result.goal_index is not None
assert result.iterations > 1
assert len(result.nodes) >= 5
assert result.path_indices()[0] == 0
assert result.path_indices()[-1] == result.goal_index
assert all(
    node.state.actuation.xi >= result.nodes[node.parent_index].state.actuation.xi
    for node in result.nodes[1:]
)
assert all(node.state.stability_margin > 0.0 for node in result.nodes)
assert any(diagnostic.accepted for diagnostic in result.edge_diagnostics)
assert result.terminal_diagnostics
successful_terminal = result.terminal_diagnostics[-1]
assert successful_terminal.success
assert successful_terminal.final_tip_error <= rrt_config.terminal_tolerance
assert successful_terminal.steps
assert all(step.jacobian_rank == 3 for step in successful_terminal.steps)
assert all(
    math.isfinite(step.jacobian_condition) and step.jacobian_condition > 0.0
    for step in successful_terminal.steps
)

# The fixed seed and goal-only sampling must reproduce the same tree outcome.
repeat_session = simder.MechanicsSession(config)
repeat_start = repeat_session.solve_initial_state()
repeat = TipSpaceRRT(repeat_session, rrt_config).plan(repeat_start, goal)
assert repeat.success
assert repeat.iterations == result.iterations
assert [node.tip_position for node in repeat.nodes] == [
    node.tip_position for node in result.nodes
]
assert result.nearest_neighbor_queries > 0
assert result.nearest_neighbor_distance_evaluations > 0
assert result.local_steering_evaluations > 0
assert result.shortlist_candidates_evaluated == 0
assert result.shortlist_non_nearest_selections == 0
assert result.shortlist_triggered_queries == 0

# The exact deterministic spatial index must preserve the brute-force tree,
# including index-based tie breaking.
linear_session = simder.MechanicsSession(config)
linear_start = linear_session.solve_initial_state()
linear = TipSpaceRRT(
    linear_session, replace(rrt_config, use_spatial_index=False)
).plan(linear_start, goal)
assert linear.success == result.success
assert linear.iterations == result.iterations
assert [node.tip_position for node in linear.nodes] == [
    node.tip_position for node in result.nodes
]
assert linear.nearest_neighbor_queries == result.nearest_neighbor_queries
assert (
    result.nearest_neighbor_distance_evaluations
    <= linear.nearest_neighbor_distance_evaluations
)

# With only the contact-free history populated, the optional shortlist is
# exactly equivalent to global-nearest selection and remains deterministic.
shortlist_session = simder.MechanicsSession(config)
shortlist_start = shortlist_session.solve_initial_state()
shortlisted = TipSpaceRRT(
    shortlist_session,
    replace(rrt_config, use_contact_history_shortlist=True),
).plan(shortlist_start, goal)
assert shortlisted.success
assert shortlisted.iterations == result.iterations
assert [node.tip_position for node in shortlisted.nodes] == [
    node.tip_position for node in result.nodes
]
assert shortlisted.shortlist_candidates_evaluated > 0
assert shortlisted.shortlist_non_nearest_selections == 0
assert shortlisted.shortlist_triggered_queries == shortlisted.nearest_neighbor_queries

# A periodic hybrid retains ordinary global-nearest steps between exact,
# deterministic shortlist queries.
periodic_session = simder.MechanicsSession(config)
periodic_start = periodic_session.solve_initial_state()
periodic = TipSpaceRRT(
    periodic_session,
    replace(rrt_config, contact_history_shortlist_period=3),
).plan(periodic_start, goal)
assert periodic.success
assert [node.tip_position for node in periodic.nodes] == [
    node.tip_position for node in result.nodes
]
assert periodic.shortlist_triggered_queries == periodic.iterations // 3
assert periodic.shortlist_candidates_evaluated == periodic.shortlist_triggered_queries

# Auxiliary history expansions use a separate random stream and never enter
# the backbone index. Consequently, the backbone remains an exact prefix of
# the unchanged deterministic baseline even while auxiliary nodes are added.
dual_session = simder.MechanicsSession(config)
dual_start = dual_session.solve_initial_state()
dual = TipSpaceRRT(
    dual_session,
    replace(
        rrt_config,
        use_dual_frontier=True,
        dual_frontier_auxiliary_period=3,
    ),
).plan(dual_start, goal)
assert dual.success
dual_backbone_tips = [
    node.tip_position for node in dual.nodes if node.frontier == "backbone"
]
assert dual_backbone_tips == [
    node.tip_position for node in result.nodes[:len(dual_backbone_tips)]
]
assert dual.backbone_expansion_queries > 0
assert dual.auxiliary_expansion_queries > 0
assert dual.shortlist_triggered_queries == dual.auxiliary_expansion_queries
assert dual.local_steering_evaluations > 0
assert any(node.frontier == "auxiliary" for node in dual.nodes)

# A nearby connector with deliberately tiny actuation steps and one connector
# iteration must fail cleanly while preserving the initial node for exploration.
restricted_config = RRTConfig(
    workspace_min=(-0.05, -0.15, -0.099),
    workspace_max=(1.05, 0.15, 0.099),
    goal_bias=1.0,
    local_tip_step=1.0e-3,
    goal_neighborhood_radius=1.0e-3,
    terminal_tolerance=1.0e-9,
    maximum_terminal_iterations=1,
    maximum_iterations=1,
    random_seed=7,
    maximum_insertion_step=1.0e-10,
    maximum_field_step=1.0e-10,
    maximum_field_magnitude=0.25,
    minimum_tip_progress=1.0e-12,
)
restricted_session = simder.MechanicsSession(config)
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

# Mechanics work is a deterministic soft cap: the current nonlinear edge is
# allowed to finish, but the tree cannot start another iteration afterward.
budget_session = simder.MechanicsSession(config)
budget_start = budget_session.solve_initial_state()
budget_config = replace(rrt_config, maximum_continuation_steps=1)
budgeted = TipSpaceRRT(budget_session, budget_config).plan(budget_start, goal)
assert not budgeted.success
assert budgeted.termination_reason == "mechanics_work_budget"
assert budgeted.continuation_attempted_steps >= 1
assert budgeted.continuation_attempted_steps == sum(
    diagnostic.attempted_steps for diagnostic in budgeted.edge_diagnostics
) + sum(
    step.attempted_steps
    for diagnostic in budgeted.terminal_diagnostics
    for step in diagnostic.steps
)

# A negative unconstrained insertion step activates dxi=0 and re-solves the
# field-only least-squares problem instead of clipping an inconsistent command.
class SyntheticState:
    actuation = simder.Actuation(0.5, [0.0, 0.0, 0.0])


diagnostic_planner = TipSpaceRRT(session, rrt_config)
active_set_diagnostic = diagnostic_planner._steer_with_diagnostics(
    SyntheticState(),
    (
        (1.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ),
    (-1.0e-3, 0.0, 0.0),
)
assert active_set_diagnostic is not None
active_set_command = active_set_diagnostic.command
assert active_set_command[0] == 0.0
assert active_set_command[1] < 0.0
assert active_set_diagnostic.insertion_constraint_active
assert active_set_diagnostic.projected_linear_residual_ratio < 1.0e-6

# If only retraction can produce the desired tip direction, projection onto
# dxi>=0 leaves the full residual and explicitly reports local infeasibility.
infeasible_diagnostic = diagnostic_planner._steer_with_diagnostics(
    SyntheticState(),
    (
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ),
    (-1.0e-3, 1.0e-6, 0.0),
)
assert infeasible_diagnostic is not None
assert infeasible_diagnostic.insertion_constraint_active
assert infeasible_diagnostic.projected_linear_residual_ratio > 0.999

# Across populated histories, projected feasible steering can intentionally
# select a slightly farther node when the spatially nearest node cannot move
# toward the sample without forbidden insertion retraction.
shortlist_nodes = [
    RRTNode(
        SyntheticState(),
        (0.0, 0.0, 0.0),
        None,
        None,
        tip_actuation_jacobian=(
            (-1.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        ),
        contact_history=(1,),
    ),
    RRTNode(
        SyntheticState(),
        (-5.0e-4, 0.0, 0.0),
        None,
        None,
        tip_actuation_jacobian=(
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        ),
        contact_history=(2,),
    ),
]
shortlist_selection = diagnostic_planner._select_shortlist_proposal(
    shortlist_nodes, (0, 1), (1.0, 1.0e-3, 0.0)
)
assert shortlist_selection is not None
assert shortlist_selection[0] == 1
assert shortlist_selection[2] == 2

print(
    "Baseline RRT: "
    f"iterations={result.iterations}, nodes={len(result.nodes)}, "
    f"edges={len(result.edge_diagnostics)}, "
    f"terminal_steps={len(successful_terminal.steps)}, deterministic=1"
)
