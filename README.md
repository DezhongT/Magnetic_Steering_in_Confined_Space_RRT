# simDER: rod in confined space

This repository contains the legacy discrete-elastic-rod simulator that will be
used as the mechanics foundation for the mechanics-informed RRT described in
[`doc/README_mechanics_informed_RRT.md`](doc/README_mechanics_informed_RRT.md).

For code review, including the governing equations, implementation map,
baseline definitions, complete test inventory, benchmark results, and failure
analysis, start with
[`doc/PLANNING_MATH_IMPLEMENTATION_REVIEW.md`](doc/PLANNING_MATH_IMPLEMENTATION_REVIEW.md).

The current refactor status, design decisions, completed work, and exact next
milestone are recorded in
[`doc/DEVELOPMENT_HANDOFF.md`](doc/DEVELOPMENT_HANDOFF.md).

## Build

The project requires a C++17 compiler, Eigen3, LAPACK, OpenGL, GLU, and GLUT.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

To build the optional Python planning interface, install pybind11 for the
selected Python interpreter and enable it explicitly:

```bash
python3 -m pip install pybind11
cmake -S . -B build-python -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DSIMDER_BUILD_PYTHON=ON
cmake --build build-python -j
ctest --test-dir build-python --output-on-failure
```

The module is produced in `build-python`; either install/package it later or
set `PYTHONPATH=build-python` while running development scripts.

## Run

Run the OpenGL demonstration from the repository root so the legacy relative
mesh paths resolve correctly:

```bash
./build/simDER option.txt
```

Command-line overrides follow `--`:

```bash
./build/simDER option.txt -- render 0 saveData 0 totalTime 0.05
```

## Test

```bash
ctest --test-dir build --output-on-failure
```

The current tests include a short headless regression check for the legacy
dynamic solver, a zero-load static-equilibrium residual check, and exact DER
state capture/restore checks. A contact-free static Newton test also verifies
convergence, preservation of clamped DOFs and physical time, and rollback after
a forced failure. Elastic finite differences, gravity-only relaxation
equivalence, an opt-in axial tip dipole, and equilibrium field sensitivity are
also validated. The independent geometry layer includes exact triangle-mesh
projection, a planar slab, and an analytic spherical shell with finite-
difference-tested clearance gradients and Hessians.
Body-contact candidate detection applies rod-radius clearance, stable IDs,
deterministic ordering, and a separate forbidden-tip safety check. A new
IPC-style planar barrier assembles its residual and Jacobian through that
detector and is finite-difference validated without changing the default
legacy contact model.

Static Newton now uses feasibility-preserving backtracking when
`planar_barrier` is selected. Nonpositive-gap trials are rejected, contact
steps must reduce the residual merit function, and failed solves restore their
exact input state. A loaded planar-contact equilibrium test verifies active
contact, positive final gap, tip safety, and backtracking behavior.

Contact selection is controlled by `contactModel`: `legacy` (default),
`planar_barrier`, `spherical_shell_barrier`, or `none`. The planar model treats
`thickness` as the half thickness on each side of `z = 0`. The spherical model
uses `shellCenter`, `shellRadius`, `shellMinusThickness`, and
`shellPlusThickness` for the admissible layer
`R-d_minus < ||x-center|| < R+d_plus`. `dBar`, `stiffness`, and
`tipSafeDistance` configure either analytic barrier domain.
The line search can be configured with `maxLineSearchIter`,
`lineSearchReduction`, and `lineSearchArmijo`.

A first fixed-active-set KKT layer can now be constructed from a solved planar
barrier equilibrium. It removes barrier terms from the physical mechanics,
maps active candidate normals into the free-DOF constraint Jacobian, initializes
positive compressive multipliers from the barrier force, and solves a symmetric
dense augmented system. This is currently a validated linearized correction,
and the planar implementation now iterates that correction to nonlinear KKT
equilibrium. Violated inactive contacts are inserted, negative-multiplier
contacts are released, unsafe tips are rejected, and failed solves roll back
exactly. The nonlinear result reports active IDs, multipliers, stationarity,
gap, complementarity, and feasibility diagnostics.

Converged planar KKT equilibria now support null-space reduced-Hessian
stability analysis and differentiated-KKT sensitivity to the three axial-tip
field components. Warm-started field re-solves preserve DER transported-frame
history and verify configuration and contact-multiplier derivatives while
requiring the active IDs to remain unchanged.

The same KKT mechanics now support the spherical shell. The contact tangent
includes `-lambda * Hessian(gap)` in the Lagrangian Hessian, and the barrier
predictor includes its corresponding curvature term. A curved fixed-active
regression validates equilibrium, positive multipliers, reduced-Hessian
stability, field sensitivity, combined-actuation continuation, and exact
rollback. Python planning configuration remains planar until Task B is
serialized and calibrated.

An adaptive magnetic-field predictor-corrector now uses those sensitivities to
advance planar KKT equilibria. It grows or shrinks field steps, rejects failed
or unstable corrections, records contact insertion/release events, recomputes
sensitivity after nonsmooth events, and rolls the complete state and field back
if the requested path cannot be completed.

Planning actuation is now represented explicitly as `(xi, B)`, and a
`PlannerState` owns the rod history, actuation, active contact IDs,
multipliers, minimum body gap, continuous tip clearance, tip-safety flag,
stability margin, and convergence diagnostics. The first
`InsertionModel` is deliberately provisional: `proximal_guide` applies a
compliant axial guide to the first free centerline vertex, with target offset
`xi` from its initial position. It provides analytic residual, tangent, and
`dF/dxi`; `none` remains the default so legacy behavior is unchanged.

The differentiated contact KKT system and adaptive continuation now support
the combined parameter ordering `(xi, Bx, By, Bz)`. Continuation can warm-start
directly from a `PlannerState`, enforces monotone insertion, preserves the
existing contact-event and stability checks, and rolls rod plus actuation back
on failure.

The optional `simder` Python module exposes a high-level `MechanicsSession`,
immutable planner-state snapshots, actuation and continuation options, and one
complete continuation-edge attempt. Newton/KKT iterations, contact updates,
and Eigen storage remain internal to C++. A minimal call sequence is:

```python
import simder

config = simder.MechanicsConfig()
config.gravity = [0.0, 0.0, 0.22]
config.initial_actuation = simder.Actuation(1e-3, [0.0, 0.1, 0.0])

session = simder.MechanicsSession(config)
start = session.solve_initial_state()
target = simder.Actuation(start.actuation.xi + 2e-4, [0.01, 0.11, 0.005])
result = session.attempt_continuation(start, target)
assert result.success
```

Set `config.domain_type = "spherical_shell"` and configure `shell_center`,
`shell_radius`, `shell_minus_thickness`, and `shell_plus_thickness` to use the
analytic curved domain. The default remains `"planar_slab"` for compatibility.

`MechanicsSession.evaluate_local_steering(state)` returns the distal-tip
position and the differentiated-KKT `3 x 4` Jacobian ordered as
`(xi, Bx, By, Bz)`. The pure-Python `simder_planning.TipSpaceRRT` uses this as
a proposal model, then validates every edge with full C++ continuation. The
baseline intentionally uses fixed-seed tip-space sampling, nearest-tip node
selection, bounded damped least-squares steering, parent indices, and
structured edge diagnostics. Nodes inside the goal neighborhood now invoke an
iterative terminal connector that relinearizes, continuation-validates every
step, reports Jacobian rank/condition and final tip error, and stops at an
explicit terminal tolerance. Failed connectors retain their nearby and
partially advanced valid nodes so exploration can continue. RRT*, rewiring,
bidirectional growth, and mechanics-aware nearest-neighbor ranking remain
deferred.

The contact-free Task-A benchmark and versioned JSON experiment logger can be
run directly after building the Python module:

```bash
PYTHONPATH=python:build-python python3 python/run_task_a.py \
  --output results/task_a.json --seeds 0 1 2
```

The default compact record stores the complete task configuration, seed,
summary metrics, and final parent-linked path. Add `--full-tree` only when all
explored nodes are needed. Task configuration files can be saved and reloaded
through `PlanningTask.save()` and `PlanningTask.load()` for deterministic
replay.

Run all three Task-A methods under identical mechanics, bounds, seeds, and
budgets with:

```bash
PYTHONPATH=python:build-python python3 python/run_task_a_comparison.py \
  --output results/task_a_comparison.json --seeds 0 1 2
```

The comparison includes mechanics-informed RRT, random-actuation RRT, and a
non-branching local goal-seeking controller. Task A is intentionally easy:
both Jacobian-guided methods currently succeed for all three seeds, while the
random-actuation method fails within the common iteration budget. Harder curved
or branching/dead-end tasks are still needed to distinguish global exploration
from local control.

The planar contact-release benchmark exercises the same methods from an
initial upper-wall contact to a contact-free target:

```bash
PYTHONPATH=python:build-python python3 python/run_planar_contact_comparison.py \
  --output results/planar_contact_comparison.json --seeds 0 1 2
```

Every successful selected path is checked for a recorded contact release,
nonnegative multipliers, nonpenetrating body gaps, nonnegative tip clearance,
tip safety, and positive reduced-Hessian stability. The current deterministic
regression succeeds for 3/3 mechanics-informed RRT seeds, 3/3 random-actuation
RRT seeds, and 3/3 local goal-seeking seeds. This validates contact transitions
and logging; it is not yet evidence of a global-planning advantage.

Run spherical-shell Task B with:

```bash
PYTHONPATH=python:build-python python3 python/run_spherical_shell_comparison.py \
  --output results/spherical_shell_comparison.json --seeds 0 1 2
```

Task B follows the inner sphere with one sustained body contact. All three
methods currently succeed for all three seeds. Median continuation-step counts
are 129 for mechanics-informed RRT, 204 for random-actuation RRT, and 21 for
local goal seeking. This validates curved equilibrium, sensitivity,
continuation, Python configuration, logging, and replay; a branching or
dead-end task is needed to test global exploration by itself.

The eccentric spherical-obstacle benchmark provides that dead end:

```bash
PYTHONPATH=python:build-python python3 python/run_dead_end_comparison.py \
  --output results/dead_end_comparison.json --seeds 0 1 2
```

The obstacle lies across the direct distal-tip path. Direct continuation and
local goal seeking stop with an explicit `tip_safety` failure, while a safe
reference route initially moves sideways. Mechanics-informed and random-
actuation RRT both succeed 3/3; local goal seeking fails 0/3. Median
continuation work is 288 steps for mechanics-informed RRT versus 675 for random
RRT, providing the first benchmark here where a global tree is necessary and
the mechanics-informed proposal reduces nonlinear mechanics work.

Run the ten-seed one-factor robustness sweep with:

```bash
PYTHONPATH=python:build-python python3 python/run_dead_end_sweep.py \
  --output results/dead_end_robustness_sweep.json \
  --seeds 0 1 2 3 4 5 6 7 8 9
```

The full profile varies obstacle radius and vertical offset, tip safety,
magnetic-field authority, and goal bias while retaining common per-case
budgets. It stores every run plus Wilson 95% success intervals and distribution
summaries. Across the current 11 cases, mechanics-informed RRT succeeds in
9--10/10 seeds. Random-actuation RRT succeeds in 7--10/10 and is most sensitive
to low goal bias; local goal seeking succeeds in 0/10 throughout.

The body-contact-history branching benchmark is run separately because its
1,000-iteration failure budget is substantially more expensive:

```bash
PYTHONPATH=python:build-python python3 \
  python/run_contact_branching_comparison.py \
  --output results/contact_branching_comparison.json --seeds 0 1 2
```

Two symmetric analytic obstacles form a narrow gate around rod vertex 14.
Direct goal actuation converges to a pinned two-contact equilibrium, while
calibrated routes around opposite sides insert and release boundary 1 or 2 and
reach the same contact-free goal. In the current three-seed comparison,
mechanics-informed RRT succeeds 3/3, random-actuation RRT succeeds 1/3, and
local goal seeking succeeds 0/3. Treat these counts as benchmark validation,
not a statistical ranking.

Run the deterministic mechanics-work study with:

```bash
PYTHONPATH=python:build-python python3 \
  python/run_contact_branching_budget_sweep.py \
  --output results/contact_branching_budget.json \
  --seeds 0 1 2 3 4 5 6 7 8 9
```

The soft work cap counts completed continuation attempts and records whether a
run ends at the goal, iteration limit, or mechanics-work limit. The runner
checkpoints after every case; use repeated `--case NAME` options to create
smaller deterministic shards. With corrected unilateral-insertion steering,
mechanics-informed RRT succeeds 4/10 at 3,000 work and 8/10 at both 7,000 and
15,000. The remaining failures separate into one terminal no-progress case and
one search miss, so increasing the raw budget alone is not justified.

Nearest-tip selection now uses an exact deterministic k-d tree with geometric
rebuilds and a linear insertion tail; disabling it reproduces the identical
tree through brute-force selection. Planner logs also include accumulated
contact-history classes and projected terminal residuals on the feasible
non-retracting actuation cone. In the two remaining high-budget failures,
seed 4 has projected residual ratio 1.0 and is locally unreachable without
retraction, while seed 5 allocates 541 of 651 nodes to mixed-contact history
and never enters the terminal neighborhood.

An optional contact-history-stratified mechanics shortlist can be enabled with
`--history-shortlist` in the budget-sweep runner. It queries the exact nearest
node in every populated accumulated-history class and selects the candidate
whose projected feasible endpoint is closest to the sample. On the same ten
seeds it scores 7/10, 8/10, 3/10, 7/10, and 8/10 across the five established
iteration/work cases, compared with 8/10, 8/10, 4/10, 8/10, and 8/10 for the
default global-nearest policy. At 15,000 work it rescues seeds 4 and 5 but
loses seeds 2 and 7, so it remains opt-in rather than becoming the default.

Hybrid modes are also available through `--history-shortlist-period`,
`--history-shortlist-stagnation`, `--history-shortlist-burst`, and
`--history-shortlist-after`. Single-iteration periodic and stagnation triggers
preserve baseline seeds 2 and 7 but do not recover seed 5. A five-iteration
burst every twenty iterations recovers all three critical seeds and reaches
seed 5 at error `2.62e-7`, but its full ten-seed 15,000-work result remains
8/10: seeds 1 and 4 replace the earlier failures. A continuous fallback
starting at iteration 300, 400, or 500 cannot recover seed 5 because its
decisive contact-history allocation occurs earlier. Global-nearest therefore
remains the default while these deterministic modes support controlled studies.

The dual-frontier planner, enabled with `--dual-frontier-period 2`, resolves
the single-tree tradeoff. Backbone expansions retain their own exact index and
the original sample stream; auxiliary nodes use a separate stream and cannot
change future backbone parents. Both frontiers share the continuation-work cap.
Across the established five cases it succeeds 10/10, 10/10, 5/10, 8/10, and
10/10, compared with global-nearest 8/10, 8/10, 4/10, 8/10, and 8/10. At
15,000 work every seed succeeds, median work is 3,174.5, and the original seed
4/5 failures reach errors `1.29e-7` and `8.46e-7`. Planner logs distinguish
backbone/auxiliary nodes and queries and count actual local-steering sensitivity
evaluations. Dual frontier is the recommended robust mode for this benchmark;
global nearest remains the simpler baseline and default.
