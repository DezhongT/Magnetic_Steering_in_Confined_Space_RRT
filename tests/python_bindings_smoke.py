import simder


config = simder.MechanicsConfig()
config.gravity = [0.0, 0.0, 0.22]
config.initial_actuation = simder.Actuation(1.0e-3, [0.0, 0.1, 0.0])

session = simder.MechanicsSession(config)
start = session.solve_initial_state()
assert start.tip_safe
assert start.tip_clearance > 0.0
assert start.minimum_body_gap > -1.0e-8
assert start.stability_margin > 0.0
start_configuration = list(start.configuration)

target = simder.Actuation(
    start.actuation.xi + 2.0e-4,
    [value + increment for value, increment in zip(
        start.actuation.field, [0.01, 0.01, 0.005]
    )],
)
options = simder.ContinuationOptions()
options.initial_step_fraction = 0.25
options.maximum_step_fraction = 0.4
options.easy_corrector_iterations = 10

accepted = session.attempt_continuation(start, target, options)
assert accepted.success
assert not accepted.rolled_back
assert accepted.stored_points >= 3
assert accepted.state.actuation.xi >= start.actuation.xi
assert accepted.state.stability_margin > 0.0
assert accepted.state.actuation.field == target.field

rejecting_options = simder.ContinuationOptions()
rejecting_options.stability_tolerance = start.stability_margin + 1.0
rejected = session.attempt_continuation(start, target, rejecting_options)
assert not rejected.success
assert rejected.rolled_back
assert rejected.state.configuration == start_configuration
assert rejected.state.actuation.xi == start.actuation.xi
assert rejected.state.actuation.field == start.actuation.field

recovered = session.attempt_continuation(start, target, options)
assert recovered.success

print(
    "Python mechanics session: "
    f"points={accepted.stored_points}, "
    f"attempts={accepted.attempted_steps}, rollback_verified=1"
)
