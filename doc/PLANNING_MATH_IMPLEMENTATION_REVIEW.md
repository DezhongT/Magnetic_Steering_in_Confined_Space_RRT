# Planning Mathematics, Implementation, Baselines, and Test Report

## 1. Purpose and source of truth

This document is the code-review guide for the mechanics and planning stack. It
answers four questions:

1. What mathematical problem is the code solving?
2. Where is each mathematical operation implemented?
3. What does each baseline mean, and what is a fair comparison?
4. What has been tested, what succeeded, and why did failures occur?

For current behavior, the source-of-truth order is:

1. the implementation and automated tests;
2. this review document;
3. [`DEVELOPMENT_HANDOFF.md`](DEVELOPMENT_HANDOFF.md), which also contains
   historical intermediate results;
4. [`README_mechanics_informed_RRT.md`](README_mechanics_informed_RRT.md), which
   began as the design specification and includes ideas not yet implemented.

This distinction matters. For example, an older handoff section records an
early planar-contact result of 2/3 for two RRT methods. The current automated
regression in [`python_planar_contact_task.py`](../tests/python_planar_contact_task.py)
requires 3/3 for all three methods. Current tests take precedence.

## 2. State, actuation, and admissibility

Let the free discrete-elastic-rod degrees of freedom be

\[
q\in\mathbb R^n.
\]

They contain vertex positions and twist variables after eliminating prescribed
degrees of freedom. A planner state is more than `q`: it also stores previous
configuration, velocities, transported directors, material directors,
tangents, reference twist, active contacts, multipliers, actuation, stability,
and safety diagnostics. This full state is required for exact restoration of
independent tree branches.

The actuation is

\[
u=\begin{bmatrix}\xi & B_x & B_y & B_z\end{bmatrix}^{T},
\]

where `xi` is the insertion coordinate and `B` is a uniform magnetic field.
Planning currently imposes

\[
\Delta\xi\ge 0,
\qquad
\|\Delta B\|\le \Delta B_{\max},
\qquad
\|B\|\le B_{\max}.
\]

The current insertion model is deliberately provisional. Vertex 2 is attached
to a compliant axial guide with energy

\[
E_{\rm ins}
=\frac{k_{\rm ins}}{2}
\left[a^T(x_2-x_{2,0})-\xi\right]^2.
\]

Therefore

\[
\frac{\partial r}{\partial\xi}=-k_{\rm ins}a.
\]

This is an architecture and derivative model, not yet a physical material-feed
or changing-deployed-length model. See
[`proximalGuideInsertionModel.cpp`](../insertion/proximalGuideInsertionModel.cpp).

## 3. Energy, residual, and tangent

The static potential used by the planning mechanics is conceptually

\[
\Pi(q;u)=E_{\rm stretch}+E_{\rm bend}+E_{\rm twist}
+E_{\rm gravity}+E_{\rm ins}-\mu\,t(q)^TB.
\]

The magnetic model is an axial distal-edge dipole, not yet a general
material-frame dipole. If the last edge is `e`, with length `l` and tangent
`t=e/l`, then

\[
E_{\rm mag}=-\mu t^TB.
\]

The residual and tangent are

\[
r(q;u)=\nabla_q\Pi(q;u),
\qquad
H(q;u)=\frac{\partial r}{\partial q}.
\]

The axial magnetic residual, exact Hessian, and field derivative are implemented
in [`tipMagneticForce.cpp`](../tipMagneticForce.cpp) and checked against finite
differences. Elastic residual/Jacobian assembly is similarly checked by
`elastic_jacobian_finite_difference`.

### Static Newton versus dynamic relaxation

Static Newton solves

\[
r(q;u)=0
\]

with line search. Dynamic relaxation integrates

\[
M\ddot q+C\dot q+r(q;u)=0
\]

until both velocity and the static residual are small. They are not different
equilibrium definitions: when relaxation reaches rest on the same stable branch,
both target `r=0`. Dynamic relaxation is useful as a robust recovery strategy;
static Newton/KKT is preferred for repeated planning solves and sensitivities.
The regression requires the two gravity-only equilibria to agree in
configuration within `1e-4`.

## 4. Contact mathematics

For every body-contact candidate `i`, the signed gap and normal reaction obey

\[
g_i(q)\ge0,
\qquad
\lambda_i\ge0,
\qquad
g_i(q)\lambda_i=0.
\]

Body contact is frictionless and allowed. Distal-tip contact is forbidden:

\[
g_{\rm tip}(q)\ge d_{\rm safe}.
\]

For a fixed active set `A`, define

\[
z=\begin{bmatrix}q\\\lambda_A\end{bmatrix},
\qquad
F_A(z;u)=
\begin{bmatrix}
r(q;u)-G_A(q)^T\lambda_A\\
g_A(q)
\end{bmatrix}.
\]

The nonlinear corrector solves `F_A=0`. Its tangent is

\[
K_A=
\begin{bmatrix}
H_L & -G_A^T\\
G_A & 0
\end{bmatrix},
\qquad
H_L=H-\sum_{i\in A}\lambda_i\nabla^2g_i.
\]

The curved-gap Hessian term is essential. Omitting it gives the wrong tangent,
stability, and sensitivity for spherical boundaries. Planar gaps have zero
Hessian, which is why the error is invisible in planar-only testing.

The barrier is used for contact prediction and KKT initialization. For
`0<g<d_hat`, the implemented energy is

\[
\phi(g)=-k(g-d_{\rm hat})^2\log(g/d_{\rm hat}).
\]

Exact unilateral validity still comes from the active-set KKT solution:
nonnegative active multipliers, nonnegative inactive gaps, and tip safety.

## 5. Stability, sensitivity, and local steering

Let `Z` span the null space of active contact constraints. The constrained
stability matrix is the reduced Hessian

\[
H_{\rm red}=Z^TH_LZ.
\]

The stored stability margin is its minimum eigenvalue. A valid planning state
requires a positive margin under the configured tolerance.

Differentiating fixed-active-set equilibrium gives

\[
K_A\frac{\partial z}{\partial u}
=-\frac{\partial F_A}{\partial u}.
\]

If `P_tip` extracts the distal position, the local actuation-to-tip Jacobian is

\[
J_{xu}=P_{\rm tip}\frac{\partial q}{\partial u}
\in\mathbb R^{3\times4}.
\]

This solve is exposed by
[`MechanicsSession::evaluateLocalSteering`](../planning/mechanicsSession.cpp).
The planner caches `J_xu` per node, because a node is immutable.

For desired tip displacement `d`, mechanics-informed steering solves the
regularized least-squares problem

\[
\min_{\Delta u}
\frac12\|J_{xu}\Delta u-d\|^2
+\frac{\rho_\xi}{2}(\Delta\xi-\Delta\xi_{\rm pref})^2
+\frac{\rho_B}{2}\|\Delta B\|^2,
\]

subject to `Delta xi >= 0` and the configured step/magnitude bounds.

The unilateral insertion constraint uses a one-constraint active set. If the
unconstrained solution requests retraction, the planner fixes
`Delta xi=0` and re-solves the three field variables. Clipping insertion after
an unconstrained solve is incorrect because the remaining field command is no
longer optimal.

The projected feasibility diagnostic is

\[
\rho_{\rm proj}
=\frac{\|J_{xu}\Delta u_{\rm feasible}-d\|}{\|d\|}.
\]

Values near zero mean the requested local displacement is representable by the
linearized feasible actuation. Values near one mean essentially none of it is.
This is a local statement at one equilibrium, not a proof that the global goal
is unreachable.

## 6. Nonlinear continuation and edge acceptance

The Jacobian command is only a proposal. Every tree edge is validated by
adaptive predictor-corrector continuation in `(xi,B)`:

1. restore the complete parent state;
2. predict the next equilibrium from differentiated KKT;
3. run the nonlinear KKT corrector;
4. add or release contacts when unilateral conditions demand it;
5. reject unsafe tips, unstable states, or failed correctors;
6. reduce the continuation step and retry;
7. grow the step after easy accepted corrections;
8. roll back exactly if the complete requested edge cannot be reached.

An edge is appended only when continuation succeeds and the tip moves by at
least `minimum_tip_progress`. The continuation-step budget is a deterministic
soft cap: an atomic edge may finish just beyond the cap, but no new tree
iteration begins afterward.

## 7. Planner implementations and baselines

### 7.1 Local goal seeking

`LocalGoalSeekingPlanner` is the no-global-search baseline. It repeatedly
steers one chain directly toward the goal with the same mechanics Jacobian,
actuation bounds, continuation, contact solver, terminal connector, and safety
rules. It is cheapest when the direct branch works, but cannot go around a
dead end because it neither samples nor branches.

### 7.2 Random-actuation RRT

This is the global-search/no-mechanics-steering baseline. It uses the same
workspace sampling distribution and goal-bias rule, nearest-node rule, tree,
continuation validator, safety checks, budgets, and terminal connector as
mechanics-informed RRT. Ordinary expansions choose random bounded nonnegative
insertion and random field increments instead of using `J_xu`.

The current random baseline uses one generator for workspace and actuation
draws. Consequently, the same seed does not imply the same sample-by-sample
workspace sequence as mechanics RRT after the first random actuation draw. The
comparison controls distributions and budgets, not paired samples. A future
paired-sample experiment should use separate sampling and proposal streams.

This baseline isolates the benefit of differentiated mechanics proposals from
the benefit of merely having a global tree.

### 7.3 Mechanics-informed global-nearest RRT

The principal single-tree baseline samples a tip target, selects the exact
Euclidean-nearest tree node in tip space, computes `d`, solves the local
least-squares steering problem, and validates it through continuation.

It is deliberately not RRT*: there is no rewiring, asymptotic optimality claim,
bidirectional growth, or learned model.

### 7.4 Linear scan versus deterministic k-d tree

These are computationally equivalent versions of the same baseline. The exact
three-dimensional k-d tree rebuilds geometrically and scans a linear insertion
tail. Equal distances select the lower node index. Regression tests require the
indexed and brute-force versions to produce identical nodes and iterations.

### 7.5 Contact-history shortlist

Every node accumulates boundary IDs visited along its ancestry. The always-on
shortlist finds the nearest node in each populated history class, computes a
feasible steering prediction for that small set, then orders candidates by:

1. predicted endpoint distance to the sample;
2. projected residual ratio;
3. node index.

Residual-first ordering was tested and rejected because it abandoned spatial
progress on 997/1000 seed-5 samples.

### 7.6 Periodic, stagnation, burst, and late-start hybrids

These single-tree variants invoke the history shortlist only periodically,
after insufficient goal progress, for a fixed burst, or after a chosen
iteration. They are diagnostic baselines for selector timing. Because their
auxiliary choices enter the same tree, they can rescue one seed while changing
and losing another deterministic branch.

### 7.7 Isolated dual-frontier RRT

The dual-frontier planner is the robust benchmark mode. It maintains:

- a backbone index containing only backbone nodes;
- the original seeded backbone sample stream;
- a separate auxiliary sample stream;
- history-stratified auxiliary expansions over all histories;
- one shared continuation-work cap.

Auxiliary nodes never enter the backbone index. Therefore auxiliary exploration
cannot alter later backbone parents. With auxiliary period 2, one auxiliary
attempt follows every two backbone iterations. The baseline remains the default
configuration; dual frontier is opt-in because broader perturbed testing is
still needed.

## 8. Terminal connector

When a node enters `goal_neighborhood_radius`, the connector repeatedly:

1. recomputes `J_xu` at the current equilibrium;
2. solves toward the remaining full tip error;
3. checks the rank and condition of `J_xu J_xu^T`;
4. validates the step through continuation;
5. accepts only strict nonlinear error reduction;
6. backtracks the desired correction if necessary.

Success requires final error no larger than `terminal_tolerance`. Valid partial
connector nodes remain in the tree after failure. This avoids falsely labeling
a near-goal tree node as an exact solution.

## 9. Fair-comparison rules

All planner comparisons use:

- the same mechanics parameters and initial equilibrium;
- the same target and workspace;
- the same actuation bounds;
- the same nonlinear continuation and contact rules;
- the same tip-safety and stability rules;
- the same iteration and continuation-work caps;
- explicit deterministic seeds;
- a fresh `MechanicsSession` per run.

A run is successful only if the planner reaches the numerical goal tolerance
and the selected path satisfies the task's required contact event (`none`,
`release`, `sustained`, or `any`). Wall time is recorded but removed from exact
replay comparisons. Continuation work and local-steering evaluation counts are
reported separately because a shortlist can spend sensitivity solves without
spending continuation steps.

## 10. Benchmark report

### 10.1 Three-seed task comparisons

| Task | Mechanics RRT | Random RRT | Local goal | What it tests |
|---|---:|---:|---:|---|
| Contact-free Task A | 3/3 | 0/3 | 3/3 | Steering, exact terminal connection, replay |
| Planar contact release | 3/3 | 3/3 | 3/3 | Initial contact, multiplier validity, release |
| Spherical-shell Task B | 3/3 | 3/3 | 3/3 | Curved KKT tangent and sustained contact |
| Tip-obstacle dead end | 3/3 | 3/3 | 0/3 | Need for branching around forbidden tip region |
| Double-obstacle contact branch | 3/3 | 1/3 | 0/3 | Alternative body-contact histories |

Interpretation:

- Task A is a sanity test, not evidence that global planning is needed. Local
  goal seeking is expected to work.
- The planar and shell tasks validate contact mechanics. They are not hard
  global-planning tasks because local goal seeking succeeds.
- The dead-end task establishes the need for a global tree. Mechanics steering
  uses fewer nonlinear continuation steps than random RRT.
- The double-obstacle task establishes that the relevant branches may be
  contact histories, not merely collision-free tip-space paths.

### 10.2 Ten-seed tip-obstacle robustness

For the baseline dead-end task, mechanics-informed RRT succeeds 9/10 with
median 291 continuation steps. Random RRT succeeds 10/10 with median 587.5
steps. Local goal seeking is 0/10. Across the eleven one-factor cases,
mechanics RRT remains 9--10/10, random RRT ranges from 7--10/10, and local goal
seeking remains 0/10.

The supported claim is not that ten seeds prove mechanics RRT has a larger
success probability. The supported claim is that global search is necessary,
and mechanics-informed proposals generally reduce nonlinear work and are less
sensitive to low goal bias.

### 10.3 Contact-branching budget matrix

The five cases are, in order: 500 iterations with 15,000 work, 1,000 iterations
with 15,000 work, and work caps 3,000, 7,000, and 15,000.

| Selector | 500 iter. | 1000 iter. | Work 3000 | Work 7000 | Work 15000 |
|---|---:|---:|---:|---:|---:|
| Global nearest | 8/10 | 8/10 | 4/10 | 8/10 | 8/10 |
| Always history shortlist | 7/10 | 8/10 | 3/10 | 7/10 | 8/10 |
| Period-20, burst-5 hybrid | not fully swept | not fully swept | not fully swept | not fully swept | 8/10 |
| Dual frontier, period 2 | **10/10** | **10/10** | **5/10** | **8/10** | **10/10** |

At 15,000 work, dual frontier has median continuation work 3,174.5 and a
Wilson 95% success interval of 0.722--1.000. Successful histories comprise five
boundary-1-only, two boundary-2-only, and three mixed paths. Former failures:

- seed 4 reaches error `1.29e-7` in 302 backbone iterations and 6,967 work;
- seed 5 reaches error `8.46e-7` in 182 iterations and 3,806 work.

At 3,000 work, dual-frontier failures are seeds 0, 1, 2, 4, and 5. At 7,000
work, failures are seeds 0 and 2. These are soft budget exhaustions, not KKT,
safety, or numerical failures.

## 11. Why runs fail

### Random RRT on Task A

The random proposal distribution does not reliably generate the small directed
actuation sequence needed within 100 iterations. Mechanics continuation is not
failing; the search simply does not reach the goal neighborhood.

### Local goal seeking at the tip obstacle

The direct branch enters the forbidden tip-clearance region. Continuation
returns `tip_safety` and rolls back. With no sampling or alternate parents,
local goal seeking cannot route around the obstacle.

### Direct actuation in the double-obstacle task

The solve succeeds mechanically but converges to a stable two-contact
equilibrium pinned between the obstacles, more than `1e-2` from the target.
This is an equilibrium-branch dead end, not a solver failure.

### Global-nearest seed 4

The tree reaches error `2.10674e-6`, just above the `2e-6` tolerance. Six
terminal trials activate `Delta xi=0`, and the final projected residual ratio
is approximately 1. The remaining correction is locally unavailable without
retraction at that equilibrium. This does not prove global infeasibility:
always-shortlist and dual-frontier runs reach the goal through another history.

### Global-nearest seed 5

The 15,000-work tree has 651 nodes: 541 mixed-history, 79 boundary-2-only, 30
boundary-1-only, and one contact-free. Its best error is about `0.00607`, and it
never invokes the terminal connector. The failure is history allocation and
node selection, not terminal accuracy or insufficient connector iterations.

### Always-shortlist and single-tree hybrids

Always-shortlist recovers seeds 4 and 5 but loses seeds 2 and 7. Period-20,
burst-5 preserves critical seeds 2, 5, and 7 but loses seeds 1 and 4 in the full
ten-seed run. Auxiliary choices alter the shared tree and therefore change all
later nearest parents. Aggregate success remains 8/10: failures are exchanged,
not eliminated.

### Sparse dual-frontier auxiliary cadence

Periods 4 and 8 recover seed 4 but not seed 5. Seed 5 needs sufficiently dense
early alternate-history growth. Period 2 supplies that density while keeping
all runs under the high work cap.

### General failure labels

| Label | Meaning |
|---|---|
| `goal_reached` | Planner reached terminal tolerance; experiment success also checks the task event |
| `iteration_budget` | Backbone iteration ceiling reached |
| `mechanics_work_budget` | Continuation-step soft cap exhausted |
| `local_stall` | Single-chain method cannot produce an improving edge |
| `tip_safety` | A continuation trial violates forbidden tip clearance |
| `no_progress` | Nonlinear terminal correction does not reduce error |
| continuation failure | Corrector/step adaptation cannot complete the edge |
| stability failure | Corrected equilibrium does not have acceptable reduced-Hessian margin |
| contact validity failure | Negative multiplier, penetration, or unresolved complementarity |

## 12. Automated test inventory

Both the normal build and strict `-Wall -Wextra -Wpedantic -Werror` build pass
36/36 tests.

| Tests | Coverage | Current result |
|---|---|---|
| 1--4 | Legacy smoke, zero-load equilibrium, state restoration, static Newton | Pass |
| 5--8 | Elastic derivatives, relaxation comparison, magnetic derivatives, field sensitivity | Pass |
| 9--13 | Triangle, slab, shell, one-obstacle, and two-obstacle geometry | Pass |
| 14--20 | Contact detection, barrier derivatives, planar equilibrium/KKT/stability/sensitivity | Pass |
| 21--25 | Curved contact, continuation, insertion derivatives, combined actuation, session facade | Pass |
| 26--29 | Python binding, baseline RRT, terminal connector, serialized experiment | Pass |
| 30--34 | Planner comparisons, contact tasks, dead-end task and robustness sweep | Pass |
| 35--36 | Contact-history mechanics and budget/dual-frontier integration | Pass |

### 12.1 Individual test contracts

| # | CTest name | Primary contract |
|---:|---|---|
| 1 | `legacy_headless_smoke` | Legacy executable completes a headless run |
| 2 | `zero_load_static_equilibrium` | Undeformed zero-load rod is an equilibrium |
| 3 | `rod_state_round_trip` | Full DER state restores exactly |
| 4 | `contact_free_static_newton` | Static Newton converges and failure rolls back |
| 5 | `elastic_jacobian_finite_difference` | Elastic tangent matches residual differences |
| 6 | `static_dynamic_relaxation_comparison` | Relaxation and Newton reach the same stable branch |
| 7 | `axial_tip_magnetic_validation` | Magnetic energy, residual, Hessian, and field derivative agree |
| 8 | `equilibrium_field_sensitivity` | Differentiated equilibrium matches central re-solves |
| 9 | `triangle_surface_projection` | Closest point and signed projection are consistent |
| 10 | `planar_slab_clearance` | Planar gaps, normals, and IDs are correct |
| 11 | `spherical_shell_clearance` | Inner/outer shell derivatives and IDs are correct |
| 12 | `spherical_obstacle_clearance` | Cavity/obstacle gap derivatives are correct |
| 13 | `double_spherical_obstacle_clearance` | Both obstacle boundaries retain stable IDs and derivatives |
| 14 | `contact_candidate_detection` | Candidate generation is deterministic and geometrically valid |
| 15 | `planar_barrier_contact_validation` | Barrier energy derivatives match finite differences |
| 16 | `loaded_planar_contact_equilibrium` | Loaded contact equilibrium is admissible |
| 17 | `planar_contact_kkt_linearization` | KKT matrix blocks and signs match numerical derivatives |
| 18 | `nonlinear_planar_contact_kkt` | Active-set nonlinear KKT converges with complementarity |
| 19 | `contact_constraint_analysis` | Reduced-Hessian stability and null space are valid |
| 20 | `contact_field_sensitivity` | Contact KKT sensitivity matches re-solves |
| 21 | `spherical_contact_kkt_mechanics` | Curved-gap Hessian enters equilibrium, stability, and tangent |
| 22 | `adaptive_field_continuation` | Predictor-corrector adapts steps and rolls back on failure |
| 23 | `proximal_guide_insertion_derivatives` | Insertion energy/residual derivatives are exact |
| 24 | `combined_actuation_continuation` | Coupled insertion/field continuation handles contact events |
| 25 | `mechanics_session_facade` | High-level session owns state, steering, continuation, rollback |
| 26 | `python_mechanics_session_smoke` | pybind11 types and session calls are usable from Python |
| 27 | `python_baseline_rrt` | Deterministic RRT, index parity, steering projection, dual isolation |
| 28 | `python_terminal_connector` | Exact connector succeeds/fails under controlled bounds |
| 29 | `python_task_a_experiment` | Task serialization, logging, and replay are exact |
| 30 | `python_task_a_comparison` | Current Task-A baseline success counts are enforced |
| 31 | `python_planar_contact_task` | Required release and all path safety/contact invariants hold |
| 32 | `python_spherical_shell_task` | Sustained curved contact remains valid for all methods |
| 33 | `python_dead_end_task` | Direct path is blocked and global methods route safely |
| 34 | `python_dead_end_sweep` | Sweep construction, Wilson intervals, persistence, replay work |
| 35 | `python_contact_branching_mechanics` | Direct pinned branch and both calibrated side routes are real |
| 36 | `python_contact_branching_budget` | Soft budgets, history metrics, dual metrics, and JSON agree |

Important invariants checked by this suite include:

- state capture/restore includes director and twist history;
- analytic residuals, Hessians, gap derivatives, and sensitivities agree with
  finite differences or independently re-solved equilibria;
- failed continuation rolls back exactly;
- active multipliers and inactive gaps satisfy unilateral contact tolerances;
- selected paths remain tip-safe and stable;
- insertion is monotone;
- indexed and linear nearest-neighbor searches build identical trees;
- fixed seeds replay exactly after excluding wall time;
- terminal success is based on nonlinear error, not a linear prediction;
- dual auxiliary nodes do not perturb the deterministic backbone sequence;
- planner, experiment, and sweep continuation-work totals agree.

Long ten-seed/five-case sweeps are intentionally not all part of routine CTest
because they take minutes to tens of minutes. Routine CTest exercises real
mechanics on reduced deterministic cases; the complete matrices above were run
with the checkpointed sweep CLI.

## 13. Reproduction commands

Build and run the complete automated suite:

```bash
cmake -S . -B build-python -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DSIMDER_BUILD_PYTHON=ON
cmake --build build-python -j
ctest --test-dir build-python --output-on-failure
```

Run the three-method task comparisons:

```bash
PYTHONPATH=python:build-python python3 python/run_task_a_comparison.py --output task_a.json
PYTHONPATH=python:build-python python3 python/run_planar_contact_comparison.py --output planar.json
PYTHONPATH=python:build-python python3 python/run_spherical_shell_comparison.py --output shell.json
PYTHONPATH=python:build-python python3 python/run_dead_end_comparison.py --output dead_end.json
PYTHONPATH=python:build-python python3 python/run_contact_branching_comparison.py --output branch.json
```

Run the recommended dual-frontier budget study:

```bash
PYTHONPATH=python:build-python \
python3 python/run_contact_branching_budget_sweep.py \
  --dual-frontier-period 2 \
  --seeds 0 1 2 3 4 5 6 7 8 9 \
  --output dual_frontier_budget.json
```

## 14. Code-review map

| Concern | Primary files |
|---|---|
| DER state and mechanics | [`elasticRod.cpp`](../elasticRod.cpp), [`world.cpp`](../world.cpp) |
| Axial magnetic model | [`tipMagneticForce.cpp`](../tipMagneticForce.cpp) |
| Insertion model | [`proximalGuideInsertionModel.cpp`](../insertion/proximalGuideInsertionModel.cpp) |
| Geometry | [`geometry/`](../geometry/) |
| Contact detection/barrier/KKT | [`contact/`](../contact/) |
| Predictor-corrector continuation | [`continuation/`](../continuation/) |
| C++ planning facade | [`mechanicsSession.cpp`](../planning/mechanicsSession.cpp) |
| RRT, selectors, terminal connector | [`rrt.py`](../python/simder_planning/rrt.py) |
| Serializable tasks | [`task.py`](../python/simder_planning/task.py) |
| Experiment records | [`experiment.py`](../python/simder_planning/experiment.py) |
| Robustness sweeps | [`sweep.py`](../python/simder_planning/sweep.py) |
| Tests | [`tests/`](../tests/) |

## 15. Known limitations and next review questions

1. The proximal-guide insertion coordinate is not yet a physical material-feed
   model.
2. The dipole is constrained to the distal edge tangent; a general material-
   frame magnet is not implemented.
3. Contact is frictionless; frictional stick/slip modes are absent.
4. Analytic slab/sphere domains use the modern KKT path; legacy triangle-mesh
   contact is not yet equivalent.
5. The continuation-work cap does not yet count Newton iterations, KKT
   factorizations, or all sensitivity-solve cost. Local-steering call count is
   logged separately as an interim metric.
6. Ten deterministic seeds have broad statistical confidence intervals.
7. The 10/10 dual-frontier result should be repeated for at least 30 seeds with
   small obstacle and target perturbations before making a general performance
   claim.
8. No RRT*, rewiring, bidirectional search, or path-cost optimality claim is
   present or implied.

The current supported conclusion is narrow but useful: for the implemented
frictionless analytic-contact model, isolated auxiliary contact-history
exploration removes the measured global-nearest branch failures without
changing the deterministic backbone, under the established 15,000-step
continuation budget.
