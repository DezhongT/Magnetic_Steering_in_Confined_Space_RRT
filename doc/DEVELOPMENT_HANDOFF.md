# Development Handoff

Last updated: 2026-08-26

## 1. Project purpose

This repository contains a C++ discrete-elastic-rod simulator (`simDER`) that
is being developed into the mechanics engine for the mechanics-informed RRT
described in
[`README_mechanics_informed_RRT.md`](README_mechanics_informed_RRT.md).

The original program models a clamped elastic rod with stretching, bending,
twisting, inertia, gravity, viscous damping, distributed magnetic loading, and
penalty/barrier contact against a triangle mesh. It can run headlessly or draw
the rod and mesh through OpenGL/GLUT.

The current repository should be regarded as a working legacy simulator plus
a tested mechanics-library refactor with planar equilibrium continuation. It
is not yet the full mechanics-informed RRT implementation.

## 2. Agreed technical direction

### Language split

The intended implementation is hybrid:

- **C++:** rod mechanics, geometry queries, contact, equilibrium/KKT solves,
  stability, sensitivity, and predictor-corrector continuation.
- **Python:** RRT, benchmark task definitions, configuration, experiment
  orchestration, logging, plotting, and paper-result analysis.
- **pybind11:** a narrow, high-level interface between the two layers. A full
  continuation attempt should execute inside C++; Python should not call into
  C++ once per force, contact, or Newton iteration.

If only one language can be used, C++ is preferred because equilibrium and
contact solves will dominate runtime and the existing analytic DER mechanics
are already implemented in C++.

### Solver strategy

Dynamic relaxation and static Newton/KKT will share exactly the same physical
static residual and Hessian.

The intended robust solve sequence is:

```text
warm-started static Newton/KKT
    -> if necessary, dynamic relaxation
    -> static Newton/KKT correction
    -> complementarity, clearance, and stability checks
```

Dynamic relaxation is valuable as a robust initial-equilibrium and recovery
method. The final static correction is required to obtain an accurate static
residual, KKT tangent, sensitivity, and stability margin.

The equations targeted by the two methods are:

```text
Static:              r(q; u) = 0
Dynamic relaxation:  M qddot + C qdot + r(q; u) = 0
```

If relaxation reaches rest, both target the same equilibrium. A converged
single implicit dynamic time step is not by itself proof of static equilibrium,
because inertia and damping can balance a nonzero static residual.

### Contact strategy

The current smooth penalty/barrier method can be retained as a fast contact
mode predictor and warm start. The planned exact contact workflow is:

```text
barrier solve
    -> infer candidate active contacts
    -> KKT correction
    -> add violated inactive contacts
    -> release contacts with negative multipliers
    -> accept a complementarity-consistent equilibrium
```

The barrier solution must not be treated directly as an exact KKT equilibrium:
its gaps are generally nonzero and depend on barrier stiffness.

### Magnetic model and symbolic derivatives

The planning specification uses a permanent magnetic dipole attached to the
distal tip in a spatially uniform field. The current distributed magnetic model
will eventually be replaced by a dedicated tip-dipole component.

SymEngine is intended for **offline** derivation and C++ code generation of the
local magnetic energy, residual, Hessian, and field derivative. SymEngine
should not be invoked inside runtime equilibrium solves. Generated C++ and its
generator script should both be committed and verified with finite-difference
tests.

For the tip material frame `(m1, m2, t)`, the target model is

```text
mu(q)       = mu1 m1(q) + mu2 m2(q) + mu3 t(q)
U_mag(q, B) = -mu(q)^T B
r_mag       = gradient_q U_mag
H_mag       = Hessian_q U_mag
F_B         = partial r_mag / partial B
```

## 3. Work completed

### 3.1 Build system

A CMake build was added in [`../CMakeLists.txt`](../CMakeLists.txt).

The source is now built as:

- `simder_core`: static C++ mechanics library containing the legacy mechanics.
- `simDER`: legacy command-line/OpenGL executable linked against `simder_core`.
- `static_equilibrium_test`: direct static-residual test executable.

CMake finds Eigen3, LAPACK, OpenGL, GLU, and GLUT instead of relying on the old
hard-coded Eigen path in `Makefile_Sample`.

### 3.2 Reproducible tests

CTest currently registers:

1. `legacy_headless_smoke`
   - Runs ten short legacy dynamic steps using
     [`../tests/smoke_options.txt`](../tests/smoke_options.txt).
   - Confirms the existing numerical execution path still builds and runs.

2. `zero_load_static_equilibrium`
   - Constructs the default straight rod with zero gravity and zero magnetic
     field.
   - Assembles only the static system.
   - Verifies that the residual is finite and below `1e-10`.
   - The measured residual is exactly `0` in the current environment.

3. `rod_state_round_trip`
   - Captures the complete mutable DER configuration and frame history.
   - Perturbs a free tip degree of freedom and verifies that the static
     residual becomes nonzero.
   - Restores the snapshot exactly and verifies that the original static
     residual is recovered.
   - Repeats exact capture/restore after several loaded dynamic steps.

4. `contact_free_static_newton`
   - Perturbs the free tip of an unloaded straight rod.
   - Solves static equilibrium without advancing time.
   - Verifies the clamped DOFs remain fixed and committed velocity is zero.
   - Forces a zero-iteration failure and verifies rollback to the exact input
     state.
   - Current convergence: residual `187.94` to `3.25898e-10` in 7 Newton
     iterations.

5. `elastic_jacobian_finite_difference`
   - Checks three deterministic random free-DOF directions at four central-
     difference step sizes.
   - Worst best relative directional error: `7.73313e-07`.

6. `static_dynamic_relaxation_comparison`
   - Applies the same modest gravity-only load to static Newton and the legacy
     damped dynamics.
   - Relaxation reaches the residual/velocity criteria at `t=5.782`.
   - Static versus relaxed configuration difference: `3.96063e-06`.

7. `axial_tip_magnetic_validation`
   - Validates tip-dipole energy gradient, Hessian, and field derivative.
   - Checks Hessian symmetry, zero net force, and the expected dipole torque.
   - Demonstrates contact-free magnetic bending with a final residual of
     `6.71812e-12`.

8. `equilibrium_field_sensitivity`
   - Computes `dq/dB = -J^{-1} F_B` from the equilibrium tangent.
   - Compares the full configuration sensitivity against re-solved central
     field perturbations.
   - Relative error: `4.87094e-12`.

9. `triangle_surface_projection`
   - Validates triangle interior, edge, and vertex closest-point regions.
   - Validates a degenerate-triangle segment fallback.
   - Checks oriented projection above, below, and outside a two-triangle patch.
   - Loads and queries the repository's 177-vertex, 312-triangle `mesh_1`.

10. `planar_slab_clearance`
    - Validates positive-inside and negative-outside clearance conventions.
    - Checks both slab boundaries, closest points, and boundary identifiers.
    - Verifies clearance normals by finite differences.
    - Checks subtraction of a rod-radius/safety margin.

11. `spherical_shell_clearance`
    - Implements the Task-B layer `R-d_minus < ||x-center|| < R+d_plus` with
      stable inner/outer boundary IDs, closest points, and positive-inside
      clearance.
    - Validates both analytic clearance gradients and Hessians against central
      finite differences, checks the curvature-aware barrier Jacobian, and
      rejects the derivative singularity at the center.

12. `spherical_obstacle_clearance`
    - Defines the smooth admissible intersection of an outer spherical cavity
      and the exterior of an eccentric excluded sphere.
    - Validates closest-boundary selection plus analytic gradient and Hessian
      for both boundaries and rejects obstacles not strictly inside the cavity.

13. `double_spherical_obstacle_clearance`
    - Selects the outer cavity or either excluded sphere using stable boundary
      IDs 0, 1, and 2.
    - Validates each analytic clearance gradient and Hessian with centered
      finite differences and rejects an obstacle outside the cavity.

14. `contact_candidate_detection`
    - Detects upper- and lower-boundary body candidates using rod-radius gap.
    - Verifies penetrations, activation ordering, stable vertex/boundary IDs,
      and inward gap normals.
    - Keeps the distal tip out of admissible body contacts and checks its
      independent safety margin.
    - Covers centered rods, a single tip, empty input, invalid parameters, and
      nonfinite vertices.

15. `planar_barrier_contact_validation`
    - Finite-difference checks the IPC-style scalar barrier energy against its
      analytic first and second gap derivatives.
    - Finite-difference checks the planar residual and Jacobian in a generic
      spatial direction.
    - Compares `contactModel planar_barrier` with `contactModel none` at the
      same rod state and verifies the exact residual/Jacobian contribution.
    - Checks cutoff behavior and explicit infeasibility of nonpositive gaps.

16. `loaded_planar_contact_equilibrium`
    - Solves a gravity-loaded planar-slab case to an active barrier-contact
      equilibrium.
    - Verifies residual `8.08352e-12`, minimum body gap `0.0146655`, positive
      tip clearance `0.0188737`, and four backtracks.
    - Uses an oversized-load, one-iteration solve to verify that infeasible
      wall-crossing trials are rejected and failure restores the input state.

17. `planar_contact_kkt_linearization`
    - Validates a manufactured symmetric KKT block system, solve residual,
      constraint equation, and positive-compression multiplier convention.
    - Builds a KKT seed from a solved planar barrier equilibrium after removing
      the barrier residual and Hessian from the physical mechanics.
    - Current seeded result: one contact, seed stationarity `2.52611e-10`, KKT
      linear residual `5.68664e-14`, predicted gap `1.73472e-18`, and corrected
      multiplier `2.83767e-4`.

18. `nonlinear_planar_contact_kkt`
    - Iterates physical-mechanics KKT corrections from the barrier warm start.
    - Retained-contact result: 10 nonlinear iterations, 17 backtracks,
      stationarity `1.92707e-9`, active gap `5.20417e-18`, and multiplier
      `2.66355e-4`.
    - Over-selected-contact result: releases the negative-multiplier contact,
      converges contact-free in 9 iterations, and finishes with minimum gap
      `0.00728365`.
    - Separately validates insertion of a violated inactive candidate and
      release of a tensile active candidate.
    - Rejects an unsafe distal-tip safety margin and restores the exact input
      state.

19. `contact_constraint_analysis`
    - Constructs null-space bases using a full-V SVD and checks `G Z = 0`.
    - Detects manufactured stable and unstable reduced-Hessian modes with
      minimum eigenvalues `2` and `-1`.
    - Validates constrained configuration and multiplier parameter derivatives,
      including the positive-compression sign convention.

20. `contact_field_sensitivity`
    - Analyzes the converged nonlinear planar-contact equilibrium on its fixed
      active manifold.
    - Current stability margin: minimum reduced-Hessian eigenvalue
      `0.00714882`.
    - Compares differentiated-KKT `dq/dB` and `dlambda/dB` against warm-started
      central re-solves with identical active IDs.
    - Relative configuration error: `6.32692e-13`; multiplier error:
      `7.745e-15`; differentiated-system residual: `1.02138e-14`.

21. `spherical_contact_kkt_mechanics`
    - Solves a fixed-active inner-sphere contact at vertex 14 with multiplier
      `0.0137523` and reduced-Hessian stability margin `0.00710001`.
    - Verifies the curved Lagrangian tangent through field-sensitivity central
      differences: relative configuration error `6.1552e-13` and multiplier
      error `6.23727e-15`.
    - Completes combined insertion/field continuation in 5 stored points and 4
      attempts, matches a direct endpoint solve, and checks exact rollback for
      a rejected stability threshold.

22. `adaptive_field_continuation`
    - Uses `dq/dB` and `dlambda/dB` to predict each field step, followed by the
      warm-started nonlinear KKT corrector and stability analysis.
    - Fixed-active path: 5 stored points, 4 accepted attempts, no rejections,
      and endpoint position error `2.14089e-12` versus a direct solve.
    - Contact-event path: 6 stored points, one explicitly recorded release,
      and successful continuation on the recomputed zero-contact manifold.
    - Verifies stability-threshold rejection restores the exact initial rod
      state and applied field.

23. `proximal_guide_insertion_derivatives`
    - Validates the provisional compliant proximal-guide energy gradient and
      configuration Hessian with central finite differences.
    - Validates the analytic residual derivative with respect to insertion
      coordinate `xi`.

24. `combined_actuation_continuation`
    - Validates fixed-active differentiated-KKT sensitivities in the parameter
      ordering `(xi, Bx, By, Bz)`.
    - Current insertion configuration error: `3.35871e-11`; multiplier error:
      `2.2482e-14`.
    - Verifies exact `PlannerState` rod/actuation round-trip and uses that state
      as the KKT continuation warm start.
    - Advances a simultaneous insertion/field path in 4 stored points while
      enforcing monotone `xi`; endpoint position error is `2.06616e-13` versus
      a direct target solve.

25. `mechanics_session_facade`
    - Exercises the high-level C++ facade used by future planner code.
    - Solves and captures an initial planner state, completes a combined
      continuation edge, and reports 4 stored points, 3 attempts, and minimum
      stability margin `0.00714881`.
    - Validates the bound tip position and `3 x 4` local actuation Jacobian
      against four central equilibrium finite differences; relative error is
      `3.11722e-9`.
    - Uses an impossible stability threshold to verify rejected-edge rollback
      and exact preservation of the input planner state.

26. `python_mechanics_session_smoke` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Runs the same solve/accepted-edge/rejected-edge workflow through the
      `simder` pybind11 module.
    - Verifies immutable state snapshots, monotone insertion, stable output,
      and exact rollback without exposing Newton or KKT internals.

27. `python_baseline_rrt` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Builds a tip-space RRT whose nodes own complete `PlannerState` snapshots.
    - Uses fixed-seed goal sampling, nearest-tip selection, the bound local
      steering Jacobian, bounded damped least-squares actuation proposals, and
      full nonlinear continuation validation.
    - Reaches a known reachable planar target in 11 RRT iterations, then uses 3
      final connector steps; the resulting tree has 23 stable nodes and is
      reproduced exactly in a second run.

28. `python_terminal_connector` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Starts inside the goal neighborhood and converges in 2 relinearized,
      continuation-validated steps to tip error `9.56482e-8`.
    - Requires full row rank and records a condition estimate consistent with
      the singular values of `J_tip_u`.
    - Verifies a one-iteration connector with `1e-10` actuator steps fails
      cleanly, preserves its nearby source node, and does not report a goal.

29. `python_task_a_experiment` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Serializes and reloads the complete versioned Task-A configuration, then
      executes contact-free seeds 0, 1, and 2.
    - All three seeds succeed with zero contact transitions; final tip errors
      are `4.13052e-9`, `6.95582e-9`, and `4.08580e-9`.
    - Writes compact JSON, verifies paths contain no active contacts, and tests
      optional full-tree inclusion.
    - Replays the saved configuration and exactly reproduces every path and
      non-timing summary field.

30. `python_task_a_comparison` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Runs mechanics-informed RRT, random-actuation RRT, and non-branching local
      goal seeking with identical Task-A mechanics, bounds, seeds, budgets,
      nonlinear continuation, terminal connector, and logging.
    - Mechanics-informed RRT succeeds 3/3 with median final error `4.13052e-9`
      and median 39 continuation attempts; local goal seeking succeeds 3/3
      with `2.53615e-9` and 27 attempts.
    - Random-actuation RRT fails 0/3 within 100 iterations, with median closest
      error `5.02652e-6` and 591 continuation attempts.
    - Verifies method labels, monotone insertion, common field bounds, contact-
      free paths, aggregate metrics, JSON output, and exact non-timing replay.

31. `python_planar_contact_task` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Calibrates a deterministic reference continuation from upper-wall
      contact to a contact-free target and requires a recorded release.
    - Runs the three shared planners for seeds 0--2: mechanics-informed RRT
      succeeds 3/3, random-actuation RRT succeeds 3/3, and local goal seeking
      succeeds 3/3.
    - Checks every selected path for the required contact transition,
      nonnegative multipliers, nonpenetrating minimum body gap, nonnegative
      continuous tip clearance, tip safety, and positive stability.
    - Serializes and reloads the task and exactly reproduces all non-timing
      comparison fields.

32. `python_spherical_shell_task` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Exposes analytic-domain selection and shell dimensions through
      `MechanicsConfig`/pybind11 and verifies old v1 task JSON retains planar
      defaults when those fields are absent.
    - Calibrates Task B from a known combined-actuation continuation and keeps
      the inner-sphere contact `(vertex 14, boundary 0)` active along every
      selected path state.
    - Mechanics-informed RRT, random-actuation RRT, and local goal seeking all
      succeed 3/3. Median continuation steps are 129, 204, and 21.
    - Checks multipliers, body gap, tip clearance/safety, stability, JSON
      output, sustained-contact success enforcement, and exact replay.

33. `python_dead_end_task` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Places a smooth eccentric spherical obstacle across the direct tip path;
      direct continuation is explicitly rejected with `tip_safety` after only
      23% progress and rolls back exactly.
    - Validates a three-edge reference route that first moves below the
      obstacle, advances laterally, and then approaches the target safely.
    - Mechanics-informed and random-actuation RRT both succeed 3/3, while
      deterministic local goal seeking fails 0/3 from tip-safety rejection.
    - Mechanics-informed RRT uses median 288 continuation steps versus 675 for
      random RRT; checks common budgets, safety, JSON, and exact replay.

34. `python_dead_end_sweep` (only with `SIMDER_BUILD_PYTHON=ON`)
    - Validates all 11 one-factor cases and both sweep profiles.
    - Checks Wilson interval edge cases, versioned JSON persistence, and exact
      non-timing replay of six real mechanics cases.

35. `python_contact_branching_mechanics` (only with
    `SIMDER_BUILD_PYTHON=ON`)
    - Verifies the direct actuation converges to the stable-ID pair
      `(vertex 14, boundary 1/2)` on a pinned equilibrium branch more than
      `1e-2` from the goal.
    - Calibrates two opposite routes that insert/release boundary 1 or 2,
      remain feasible and stable, and reach the same contact-free goal within
      `1.1e-8`.
    - Serializes and reloads all two-obstacle task parameters.

36. `python_contact_branching_budget` (only with
    `SIMDER_BUILD_PYTHON=ON`)
    - Validates all five iteration/work cases and the 3,000-step regression
      profile.
    - Exercises explicit mechanics-budget termination, aggregate reason
      counts, versioned JSON output, and the deterministic soft-cap contract.

Latest result:

```text
legacy_headless_smoke             Passed
zero_load_static_equilibrium      Passed
rod_state_round_trip              Passed
contact_free_static_newton        Passed
elastic_jacobian_finite_difference Passed
static_dynamic_relaxation_comparison Passed
axial_tip_magnetic_validation     Passed
equilibrium_field_sensitivity     Passed
triangle_surface_projection       Passed
planar_slab_clearance              Passed
spherical_shell_clearance          Passed
spherical_obstacle_clearance       Passed
double_spherical_obstacle_clearance Passed
contact_candidate_detection        Passed
planar_barrier_contact_validation  Passed
loaded_planar_contact_equilibrium  Passed
planar_contact_kkt_linearization   Passed
nonlinear_planar_contact_kkt       Passed
contact_constraint_analysis        Passed
contact_field_sensitivity          Passed
spherical_contact_kkt_mechanics    Passed
adaptive_field_continuation        Passed
proximal_guide_insertion_derivatives Passed
combined_actuation_continuation     Passed
mechanics_session_facade            Passed
python_mechanics_session_smoke      Passed
python_baseline_rrt                 Passed
python_terminal_connector           Passed
python_task_a_experiment             Passed
python_task_a_comparison             Passed
python_planar_contact_task           Passed
python_spherical_shell_task          Passed
python_dead_end_task                 Passed
python_dead_end_sweep                Passed
python_contact_branching_mechanics   Passed
python_contact_branching_budget      Passed

100% tests passed, 0 tests failed
```

### 3.3 Startup and ownership cleanup

The following bounded safety fixes were made without changing the intended rod
model:

- The executable checks that an option-file argument was supplied.
- Failure to open the option file returns a failure exit status.
- Configuration fields have deterministic defaults rather than indeterminate
  primitive values.
- Missing option lookup throws instead of dereferencing `map::end()`.
- The global `world` object is managed through `std::unique_ptr`, avoiding a
  shallow assignment from a temporary `world`.
- Copying `world` is disabled.
- `world`, `timeStepper`, and `elasticRod` release their owned allocations.
- Array allocations in `elasticRod` now use `delete[]` correctly.
- The LAPACK `dgbsv` status code is returned and checked.
- The obsolete `helixradius` and `helixpitch` entries were removed from
  `option.txt`; the corresponding helix construction code was already
  commented out.

### 3.4 Triangle-loop correction

The original contact and rendering loops stopped at `triangle_count - 1`, so
the final triangle was never processed. Both loops now include every triangle.

### 3.5 Static versus dynamic assembly extraction

The force/Jacobian calls previously embedded directly in
`world::updateTimeStep()` were separated into:

```cpp
void world::assembleStaticSystem();
void world::assembleDynamicSystem();
double world::computeResidualNorm() const;
double world::getStaticResidualNorm();
```

`assembleStaticSystem()` clears the current Newton system and assembles:

- stretching force and Jacobian;
- bending force and Jacobian;
- twisting force and Jacobian;
- gravity force and Jacobian;
- legacy magnetic force and Jacobian;
- penalty/barrier contact force and Jacobian.

`assembleDynamicSystem()` calls `assembleStaticSystem()` and then adds:

- inertia force and Jacobian;
- damping force and Jacobian.

The legacy implicit time-step loop now uses `assembleDynamicSystem()`. This is
the first shared mechanics seam needed by dynamic relaxation and static Newton.

`getStaticResidualNorm()` prepares the current rod geometry, assembles only the
static system, and returns the norm over unconstrained degrees of freedom. It
is a diagnostic interface, not yet a static equilibrium solver.

### 3.6 Explicit DER state capture and restoration

[`../rodState.h`](../rodState.h) defines a copyable `RodState` value containing:

- current configuration, previous configuration, and velocity;
- current and previous reference directors;
- current and previous material directors;
- current and previous tangents;
- current and previous reference twist;
- current edge lengths, curvature binormal, and curvature.

The director and twist-history fields are intentionally included. Restoring
only position/twist DOFs would mix a saved RRT branch with the parallel-
transport history of whichever branch was evaluated most recently.

`elasticRod` now provides:

```cpp
RodState captureState() const;
void restoreState(const RodState& state);
void applyFreeDofIncrement(const Eigen::VectorXd& increment);
```

`restoreState()` validates all dimensions before modifying the rod. `world`
exposes small forwarding methods and the number of free DOFs for tests and the
future solver interface.

### 3.7 Owned static evaluation and contact-free Newton solver

[`../staticEvaluation.h`](../staticEvaluation.h) stores an owned residual and
an owned copy of the LAPACK general-banded Jacobian, including its upper and
lower bandwidths. This is necessary because `dgbsv` overwrites both inputs.

`timeStepper` now provides:

```cpp
StaticEvaluation captureEvaluation() const;
int solveBandedSystem(
    const StaticEvaluation& evaluation,
    Eigen::VectorXd& solution) const;
```

`solveBandedSystem()` copies the evaluation into local LAPACK workspaces, so
the original evaluation remains available for diagnostics and future multiple
right-hand sides.

[`../equilibriumResult.h`](../equilibriumResult.h) records solver success,
initial/final residuals, Newton iteration count, LAPACK status, and the returned
`RodState`.

`world::solveStaticEquilibrium()` now:

1. snapshots the input rod state;
2. repeatedly evaluates the static residual/Jacobian;
3. solves `J dq = r` using copied banded storage;
4. applies `q <- q - dq` only to free DOFs;
5. commits a successful equilibrium with zero velocity and consistent frame
   history;
6. restores the exact input state after nonfinite residuals, linear-solver
   failures, or exhaustion of `maxIter`.

The solver does not advance `currentTime` or `timeStep`.

### 3.8 Elastic derivative and relaxation validation

`StaticEvaluation::multiplyJacobian()` applies the LAPACK banded Jacobian to a
vector without factorization or mutation. The elastic directional derivative
test restores identical DER history before every positive and negative
perturbation, preventing frame-history contamination.

A gravity-only comparison now demonstrates that dynamic relaxation and static
Newton approach the same stable equilibrium when relaxation is stopped using
both static-residual and velocity criteria.

### 3.9 Axial distal-tip dipole and field sensitivity

The legacy distributed magnetic model remains the default for backward
compatibility. Configuration now supports:

```text
magneticModel legacy
```

or the planning-oriented sanity model:

```text
magneticModel axial_tip
tipDipoleMoment <scalar moment magnitude>
baVector <Bx> <By> <Bz>
```

[`../tipMagneticForce.cpp`](../tipMagneticForce.cpp) implements an ideal dipole
aligned with the final edge:

```text
U = -mu t^T B
```

It provides an analytic local residual, Hessian, and residual derivative with
respect to the three field components. The resulting endpoint forces have zero
sum and produce torque `mu t x B`.

`world::computeConfigurationFieldSensitivity()` solves three copied banded
systems to obtain the full configuration derivative:

```text
dq/dB = -J^{-1} F_B
```

`world::setAppliedField()` supports re-solving nearby field values without
reconstructing the mechanics engine. This runtime update is currently enabled
only for `magneticModel=axial_tip`.

SymEngine and SymPy were not installed on the development workstation. The
axial model was derived analytically and independently finite-difference
validated. A general non-axial material-frame dipole remains deferred until
the terminal-frame convention and symbolic generation dependency are fixed.

### 3.10 Independent geometry query layer

The new [`../geometry/`](../geometry/) directory separates surface projection
from admissible-domain clearance:

```text
Surface::project(point)        -> SurfaceProjection
ConfinedDomain::query(point)   -> DomainQuery
```

`SurfaceProjection` reports unsigned distance, oriented signed offset, closest
point, surface normal, and primitive ID. `DomainQuery` reports clearance
(positive when admissible), its gradient direction, closest boundary point,
boundary ID, and admissibility.

`projectPointToTriangle()` implements exact Voronoi-region projection and a
longest-segment-style fallback for degenerate triangles.

`TriangleMeshSurface` validates indices and degeneracy, loads the existing text
mesh format, and performs exact brute-force closest-triangle queries. The
brute-force implementation is intentional until correctness is established;
BVH/AABB acceleration can be added without changing the interface.

`PlanarSlabDomain` represents an infinite oriented slab with independent minus
and plus thicknesses. Its clearance normal is explicitly the gradient of
clearance, so it is directly suitable for a future contact gap Jacobian.
`applyClearanceMargin()` converts centerline clearance to rod-body clearance or
applies a tip safety margin.

The geometry layer is connected to the new opt-in planar barrier contact. The
existing demonstration retains its prior behavior because `legacy` remains
the default contact model.

### 3.11 Contact-candidate and tip-safety detection

[`../contact/contactDetector.cpp`](../contact/contactDetector.cpp) consumes any
`ConfinedDomain` and produces a pure `ContactDetectionResult` without applying
forces or mutating rod state.

For every body vertex except the distal tip:

```text
gap = domain clearance - rod radius
```

Candidates with `gap <= activationDistance` retain the rod vertex, boundary
ID, centerline clearance, gap, center position, closest boundary point, and
clearance-gradient normal. Candidates are sorted deterministically by gap,
then vertex and boundary ID.

The final vertex is evaluated only as a forbidden-contact tip:

```text
tip safe <=> domain clearance >= tipSafeDistance
```

The result also records the minimum gap over all body vertices, even when no
candidate is inside the activation band. Invalid parameters, empty rods,
nonfinite vertices, and nonfinite domain responses are rejected explicitly.

### 3.12 Planar barrier contact assembly

[`../contact/barrierPotential.cpp`](../contact/barrierPotential.cpp) implements
the scalar IPC-style energy for `0 < g < dBar`:

```text
E(g) = -stiffness (g - dBar)^2 log(g / dBar)
```

The energy and its first two derivatives are zero for `g >= dBar`; a
nonpositive gap is explicitly marked infeasible. For the planar slab, whose
gap Hessian is zero, the contact component assembles

```text
residual = E'(g) n
Jacobian = E''(g) n n^T
```

where `n` is the inward clearance gradient from the candidate detector. The
distal tip remains excluded from body-contact forces and is only evaluated by
the separate safety query.

`world` accepts three contact modes: `legacy`, `planar_barrier`, and `none`.
The default is `legacy`, preserving existing runs. The planar model uses a
symmetric slab about `z = 0`, with `thickness` on each side, and is currently
restricted to `PlanarSlabDomain` so the zero-gap-Hessian assumption cannot be
silently reused for curved geometry.

### 3.13 Feasibility-preserving planar-contact Newton solve

For `contactModel planar_barrier`, `world::solveStaticEquilibrium()` now saves
the complete `RodState` before every Newton trial and uses configurable
backtracking. A nonpositive gap raises the dedicated
`InfeasibleContactError`, causing step reduction rather than solver exit. When
the current or trial configuration contains contact candidates, acceptance
uses Armijo sufficient decrease on `0.5 ||r||^2`.

Full Newton steps are retained while both ends of a trial are outside the
contact activation band. This preserves the previously validated convergence
of the legacy approximate tangent while applying strict merit globalization
where contact is involved. Every rejected trial restores configuration and
transported-frame history; exhaustion restores the solve's original input
state. `EquilibriumResult` reports total backtracks, infeasible-trial count,
minimum accepted step length, and line-search failure.

### 3.14 Fixed-active-set planar KKT linearization

[`../contact/contactKktSystem.cpp`](../contact/contactKktSystem.cpp) defines an
independent dense KKT builder and solver. The compressive multiplier convention
is

```text
r(q) - G(q)^T lambda = 0
g(q) >= 0
lambda >= 0
g_i lambda_i = 0
```

The symmetric augmented matrix uses the internal dual increment
`dmu = -dlambda`:

```text
[ H  G^T ] [dq ] = -[r - G^T lambda]
[ G   0  ] [dmu]    [g             ]
```

`world::buildPlanarContactKktSeed()` refreshes the stable candidate IDs at the
current barrier state, maps each planar normal into free DOFs, initializes
`lambda = -E'(g)`, and separately assembles the physical residual/Hessian with
the barrier disabled. Constrained-only candidates are omitted. The returned
seed is independent of the legacy LAPACK band storage except for converting an
owned `StaticEvaluation` into a dense mechanics Hessian.

This low-level API performs one fixed-contact linearized correction and does
not mutate `world`; the higher-level nonlinear corrector below applies and
relinearizes it and manages the active set.

### 3.15 Nonlinear planar KKT corrector

`world::solvePlanarContactKktEquilibrium()` now performs the full planar
barrier-to-KKT correction sequence within a rollback boundary:

```text
capture input state
  -> solve feasible barrier equilibrium
  -> seed stable active IDs and positive compression estimates
  -> assemble physical mechanics without barrier terms
  -> solve and globalize the dense KKT correction
  -> relinearize active boundary-specific gaps
  -> add violated inactive contacts
  -> release negative-multiplier contacts
  -> accept only after KKT, feasibility, complementarity, and tip checks
  -> otherwise restore the exact input state
```

`PlanarSlabDomain::queryBoundary()` evaluates a stable boundary ID even when
the nearest-boundary query could switch. Active-gap error is scaled into force
units with `characteristicForce / dBar` in the nonlinear merit function.

[`../contact/contactKktEquilibriumResult.h`](../contact/contactKktEquilibriumResult.h)
owns the accepted rod state, active candidate IDs, multipliers, iteration and
active-set counts, backtracks, stationarity, active-gap norm,
complementarity, minimum body gap, minimum multiplier, tip safety, and rollback
status. The tolerances are configurable through `kktGapTolerance`,
`kktMultiplierTolerance`, `kktComplementarityTolerance`, and
`kktMaxActiveSetUpdates`.

### 3.16 Contact stability and differentiated-KKT sensitivity

[`../contact/contactConstraintAnalysis.cpp`](../contact/contactConstraintAnalysis.cpp)
computes a full-V SVD null space `Z` of the active constraint Jacobian and the
symmetric reduced Hessian

```text
H_reduced = Z^T (H + H^T) Z / 2.
```

The result reports constraint rank, `Z`, `H_reduced`, its eigenpairs and
minimum eigenvalue, `||G Z||`, the original Hessian symmetry error, and a
stability classification. Empty active sets reduce to `Z = I`.

For the axial tip-dipole field derivative `F_B`, constrained sensitivity uses

```text
[ H  G^T ] [q_B ] = [-F_B]
[ G   0  ] [-l_B]   [  0 ]
```

and returns full-configuration `dq/dB` plus `dlambda/dB`.
`world::correctPlanarContactKktEquilibrium()` re-solves a nearby field from a
supplied KKT result without re-entering the singular barrier at zero gap. This
warm start also preserves reference-director history; cold-starting symmetric
plus/minus worlds can select different twist-coordinate gauges even when their
physical positions agree.

### 3.17 Adaptive magnetic-field predictor-corrector continuation

[`../continuation/fieldContinuation.h`](../continuation/fieldContinuation.h)
defines configurable step fractions, reduction/growth factors, stability
tolerance, easy-correction threshold, and attempt budget. Every accepted point
owns its field, nonlinear KKT result, stability analysis, path and step
fractions, configuration/multiplier predictor errors, and contact-event counts.

`world::continuePlanarContactField()` first establishes a stable planar KKT
equilibrium, predicts each trial with the differentiated KKT sensitivity,
corrects it from the predicted full `RodState`, and rejects failed, unsafe, or
insufficiently stable trials. Rejections restore the last accepted state and
reduce the step. Easy fixed-active corrections grow the step conservatively.
Contact insertion or release is recorded as a nonsmooth event; the following
step recomputes sensitivity on the new active manifold.

Multiplier predictor error is matched by stable contact ID, so different
vector dimensions across an event are handled explicitly. If the minimum step
or attempt budget is exhausted, both rod state and applied field roll back to
the continuation input.

### 3.18 Explicit planning actuation and state

[`../actuation.h`](../actuation.h) defines `Actuation` as insertion coordinate
`xi` plus applied field `B`. [`../plannerState.h`](../plannerState.h) defines a
self-contained planner-node value with the full `RodState`, actuation, active
contact IDs, multipliers, stability margin, and KKT convergence/feasibility
diagnostics. `world::capturePlannerState()` and `restorePlannerState()` provide
the tested boundary between mechanics and a future planner.

The continuation overload taking a `PlannerState` is important at contact: it
restores the known active set and multipliers and enters the nonlinear KKT
corrector directly. It therefore does not re-enter the open-gap barrier at an
exact zero-gap KKT state.

### 3.19 Provisional insertion model and combined continuation

[`../insertion/insertionModel.h`](../insertion/insertionModel.h) separates the
meaning of `xi` from the solver. The first implementation,
[`../insertion/proximalGuideInsertionModel.h`](../insertion/proximalGuideInsertionModel.h),
is intentionally provisional: it applies a compliant axial guide to vertex 2,
the first free centerline vertex. Its target is the vertex's initial position
plus `xi` along `insertionAxis`, with energy

```text
U_guide = 0.5 * insertionStiffness *
          (insertionAxis dot (x_2 - x_2_initial) - xi)^2.
```

This convention is useful for validating the architecture, derivatives, and
continuation, but it is not yet a material-feed or changing-deployed-length
model. `insertionModel none` remains the default. Selecting
`insertionModel proximal_guide` enables it; `insertionCoordinate`,
`insertionStiffness`, and `insertionAxis` configure the guide.

The model assembles its residual and tangent into both static and dynamic
systems and provides analytic `dF/dxi`. Differentiated KKT sensitivity now
returns four columns ordered `(xi, Bx, By, Bz)`. The combined adaptive
continuation uses the same nonlinear corrector, stability gate, active-set
events, step control, and rollback semantics as field continuation, while
rejecting any target with `xi_target < xi_start`.

### 3.20 High-level C++/Python planning interface

[`../planning/mechanicsSession.h`](../planning/mechanicsSession.h) defines the
value-oriented `MechanicsSession` facade. `MechanicsConfig` configures the
validated planar-barrier, axial-tip, proximal-guide scenario;
`solveInitialState()` returns a self-contained `PlannerState`; and
`attemptContinuation()` performs one complete edge attempt inside C++ and
returns its resulting state, step counts, contact-event totals, reached path
fraction, and minimum stability margin. A rejected attempt returns the exact
input state and leaves the session restored to it.

[`../python/bindings.cpp`](../python/bindings.cpp) exposes this facade as the
optional `simder` module. Python can access `MechanicsConfig`, `Actuation`,
`ContinuationOptions`, read-only `PlannerState` snapshots,
`ContinuationResult`, and `MechanicsSession`. Configuration vectors accept any
three-element Python sequence. Returned rod configurations and multipliers are
owned Python lists, avoiding borrowed Eigen storage and preventing Python from
mutating internal mechanics state.

`MechanicsSession::evaluateLocalSteering()` reconstructs the fixed-active KKT
equilibrium from a `PlannerState`, evaluates the existing four-column
actuation sensitivity, and extracts only the distal position rows. Python sees
the tip position, the `3 x 4` Jacobian ordered `(xi, Bx, By, Bz)`, and the
differentiated-system residual. It never receives the KKT matrix or DER degree-
of-freedom indexing. `PlannerState.tip_position` is also available as an owned
three-component value.

The module is disabled by default through `SIMDER_BUILD_PYTHON=OFF`; the C++
library and its tests do not require Python or pybind11. When enabled, CMake
first searches for pybind11 normally and then asks the selected interpreter
for `python -m pybind11 --cmakedir`. Newton, KKT, sensitivity, active-set, and
continuation loops remain entirely inside C++.

### 3.21 Deterministic baseline tip-space RRT

[`../python/simder_planning/rrt.py`](../python/simder_planning/rrt.py)
implements the deliberately minimal Phase-4 planner. Every `RRTNode` owns a
complete immutable `PlannerState`, cached tip position, optional cached local
Jacobian, parent index, and accepted edge command. Sampling uses a local
`random.Random` instance with an explicit seed; node selection is Euclidean
nearest-neighbor distance in tip space only.

Local inverse steering solves a four-variable damped normal equation with
partial-pivot Gaussian elimination, clamps insertion to a nonnegative bounded
increment, bounds field increments, and projects the applied field onto its
magnitude limit. The Jacobian is only a proposal: every edge is accepted only
after `MechanicsSession.attempt_continuation()` succeeds and produces minimum
tip progress. Per-edge diagnostics record the sample, parent, acceptance,
continuation attempts/rejections, contact events, and reached path fraction.

The implementation validates workspace bounds, goal and initial tip location,
actuation limits, regularization, iteration budget, and deterministic seed. It
returns the complete tree plus parent-based path reconstruction. It does not
include RRT*, rewiring, a mechanics-aware nearest-neighbor shortlist, partial-
edge insertion, or bidirectional growth.

### 3.22 Iterative exact terminal connector

When a node lies within `goal_neighborhood_radius`, the RRT repeatedly
recomputes its local tip Jacobian, solves the bounded steering problem toward
the remaining full tip error, validates the actuation through C++
continuation, and measures the true nonlinear error. It succeeds only below
`terminal_tolerance`, within `maximum_terminal_iterations`, and only accepts
connector steps that strictly reduce tip error while producing the configured
minimum tip motion.

Each valid connector equilibrium is appended as a normal parent-linked
`RRTNode`. If a later connector step fails or ceases to improve, the source and
all valid partial states remain in the tree and ordinary RRT exploration
continues. `TerminalDiagnostic` records source/final indices, initial/final tip
errors, success, and each step's rank, spectral condition estimate,
continuation-failure flag, and before/after errors.

The rank and condition estimates come from the eigenvalues of the `3 x 3`
Gram matrix `J_tip_u J_tip_u^T`, computed by a small symmetric Jacobi method
without adding NumPy as a runtime dependency. Full row rank is required for a
terminal attempt. The current reachable regression converges in two connector
steps to `9.56482e-8`; deliberately restrictive actuator bounds exercise clean
failure and source-node preservation.

### 3.23 Serializable Task A and experiment logging

[`../python/simder_planning/task.py`](../python/simder_planning/task.py)
defines versioned, JSON-serializable values for mechanics parameters,
continuation options, RRT configuration, target tip, and the explicit seed
set. `contact_free_task_a()` fixes a zero-gravity, zero-contact sanity case and
a reachable target generated from the validated mechanics. The task can be
saved and loaded without importing mutable C++ state; binding objects are
created only when a run begins.

[`../python/simder_planning/experiment.py`](../python/simder_planning/experiment.py)
runs every configured seed in a fresh `MechanicsSession`. Each versioned run
record includes success, final error, RRT/edge counts, failed continuations,
continuation step/rejection counts, terminal attempts/failures/steps, contact
transitions, peak field, total insertion, minimum path stability, wall time,
and the complete selected path. Compact mode stores no full tree; an explicit
flag adds every explored node. `RunRecord.deterministic_dict()` excludes only
wall time for exact replay comparison.

[`../python/run_task_a.py`](../python/run_task_a.py) is the command-line entry
point. The regression serializes Task A, reloads it, runs seeds 0--2, writes
and parses the result JSON, and replays all seeds exactly. The current runs use
10--14 RRT iterations, have no failed continuations or contact transitions,
and reach final tip errors below `7e-9`.

### 3.24 Task-A baseline comparisons

`TipSpaceRRT` now takes an explicit proposal policy. The
`mechanics_informed` policy retains differentiated-KKT inverse steering;
`random_actuation` samples nonnegative bounded insertion and a random field
increment, then applies the identical field-magnitude projection. It still
uses the same tip-space samples, nearest-node rule, nonlinear continuation,
terminal connector, and node/diagnostic structures, but does not evaluate
`J_tip_u` for ordinary expansion.

`LocalGoalSeekingPlanner` uses the same bounded Jacobian steering and
continuation along one parent-linked chain directed only at the goal. It never
samples or branches. Its seed is retained in the shared task record for fair
configuration, but does not affect its deterministic trajectory.

`run_comparison()` executes `mechanics_informed_rrt`,
`random_actuation_rrt`, and `local_goal_seeking` through the same Task-A
runner. Method identity is stored on each run. The versioned comparison record
contains all per-seed experiments plus success rate and median final error,
wall time, and continuation work for each method; its deterministic view
removes only timing fields. [`../python/run_task_a_comparison.py`](../python/run_task_a_comparison.py)
provides the CLI.

Task A intentionally does not establish a global-planning advantage: local
goal seeking succeeds 3/3 and is cheapest, mechanics-informed RRT succeeds
3/3, and random-actuation RRT fails 0/3 in the common 100-iteration budget.
This validates the comparison machinery and confirms that a contact/branching
task is required for the intended global-planning claim.

### 3.25 Planar contact-release benchmark and clearance logging

`PlannerState` and the nonlinear KKT result now carry the continuous distal
tip clearance in addition to `tip_safe`; the existing minimum body gap is also
bound and persisted. Experiment paths log both values at every node together
with the minimum active multiplier. Each run summarizes minimum body gap,
minimum tip clearance, minimum multiplier, and active-set transitions along
the selected path. A planner result counts as successful only if its selected
path also contains the event required by the serialized task.

`planar_contact_release_task()` defines a versioned planar-slab task with
gravity `(0, 0, 0.22)`, initial actuation `xi = 0.001` and
`B = (0, 0.1, 0)`, and a target calibrated by continuation to
`B = (0, 0.1, -0.25)`. The initial equilibrium has upper-wall contact with
multiplier `7.60705e-5`, body gap numerically zero, tip clearance
`0.00300218`, and stability margin `0.00714881`. The reference continuation
records one release and reaches a contact-free equilibrium.

[`../python/run_planar_contact_comparison.py`](../python/run_planar_contact_comparison.py)
runs the same mechanics-informed RRT, random-actuation RRT, and local
goal-seeking implementations under common bounds, budgets, and seeds. Current
success counts are 2/3, 2/3, and 3/3 respectively. Median final errors are
`6.35300e-6`, `6.72936e-6`, and `3.96599e-6`; median continuation step counts
are 280, 134, and 15. Every successful path contains the required release and
passes the multiplier, gap, tip-clearance, tip-safety, and stability checks.
Because local goal seeking remains best on this task, the result validates
contact-transition infrastructure but makes no global-planning claim.

### 3.26 Analytic spherical shell and curved KKT tangent

[`../geometry/sphericalShellDomain.cpp`](../geometry/sphericalShellDomain.cpp)
implements the Task-B admissible layer
`R-d_minus < ||x-center|| < R+d_plus`. Boundary ID 0 is the inner sphere and
ID 1 is the outer sphere. For radial unit vector `u`, the clearance gradients
are `u` and `-u`, while the Hessians are
`(I-u*u^T)/r` and `-(I-u*u^T)/r`. The center is rejected because these
derivatives are not uniquely defined there. `DomainQuery` and
`ContactCandidate` now carry this Hessian, and the detector preserves it with
the stable contact ID.

The analytic barrier Jacobian now contains both
`phi''(g) grad(g) grad(g)^T` and `phi'(g) Hessian(g)`. The exact-contact KKT
tangent uses the Lagrangian Hessian
`H_physical - sum(lambda_i Hessian(g_i))`; the same term is included in the
barrier-derived KKT seed. Planar behavior is unchanged because its gap Hessian
is identically zero.

The calibrated curved regression uses a shell centered at `(0, 0, 1)` with
reference radius `1`, inner thickness `0.1`, and outer thickness `1`, under
gravity `(0, 0, 0.5)`. It converges with one inner-boundary contact at vertex
14, multiplier `0.0137523`, positive tip clearance, and stability margin
`0.00710001`. Differentiated-KKT field sensitivity agrees with central
re-solves to relative configuration error `6.1552e-13` and multiplier error
`6.23727e-15`. Combined `(xi, B)` continuation reaches its direct-solve target
in five points and the forced instability case rolls rod and actuation back
exactly.

The older C++ class and method names containing `Planar` are retained for API
compatibility even though their analytic-domain KKT path now handles the
spherical shell. New generic method names remain a future API cleanup.

### 3.27 Serialized spherical-shell Task B

`MechanicsConfig` and its pybind11 interface select `planar_slab`,
`spherical_shell`, or `spherical_obstacle` and carry the corresponding domain
parameters. `MechanicsTaskConfig` persists the same fields. They have planar
defaults, and loading older schema-v1 JSON that omits them remains supported.

`spherical_shell_task_b()` uses the validated shell centered at `(0, 0, 1)`,
reference radius `1`, inner thickness `0.1`, outer thickness `1`, and gravity
`(0, 0, 0.5)`. Its target is independently calibrated from actuation
`xi = 0.005`, `B = (0.1, 0.2, -0.25)`. The selected-path contract
`required_contact_event = "sustained"` requires every path state to retain at
least one active contact; Task B retains `(vertex 14, inner boundary 0)`.

[`../python/run_spherical_shell_comparison.py`](../python/run_spherical_shell_comparison.py)
runs all three methods for seeds 0--2. All currently succeed 3/3. Median final
errors are `1.97980e-6` for mechanics-informed RRT, `1.75787e-6` for random
actuation, and `1.98591e-6` for local goal seeking. Median continuation-step
counts are 129, 204, and 21. Every selected path passes multiplier, gap,
tip-clearance, tip-safety, stability, and sustained-contact checks, and exact
non-timing replay is verified. Local goal seeking is still cheapest, so Task B
does not establish a global-planning advantage.

### 3.28 Eccentric spherical-obstacle dead end

[`../geometry/sphericalObstacleDomain.cpp`](../geometry/sphericalObstacleDomain.cpp)
defines the admissible intersection of an outer spherical cavity and the
exterior of an eccentric excluded sphere. Boundary 0 is the cavity and
boundary 1 is the obstacle. Both expose analytic signed clearance, gradient,
Hessian, closest point, and stable ID; central finite differences validate the
derivatives. The domain is exposed through C++, pybind11, and serialized task
configuration.

`spherical_obstacle_dead_end_task()` places a radius `2.5e-4` obstacle midway
across the direct distal-tip route and requires `1e-4` additional tip safety
clearance. Direct continuation toward the known target is rejected explicitly
with failure reason `tip_safety` after reaching path fraction `0.233248`, then
restores the exact start state. A calibrated safe route first changes the field
to move the tip below the obstacle, advances insertion and lateral field while
remaining below it, and finally approaches the target from the safe side.

[`../python/run_dead_end_comparison.py`](../python/run_dead_end_comparison.py)
uses the same bounds, budgets, seeds, terminal criteria, nonlinear continuation,
and safety rules for all methods. Mechanics-informed RRT succeeds 3/3 with
median final error `3.27307e-8` and 288 continuation steps. Random-actuation
RRT succeeds 3/3 with `6.50644e-8` and 675 steps. Deterministic local goal
seeking fails 0/3 with median closest error `0.00112438`; each run records an
explicit tip-safety failure. Thus the global tree is necessary on this task,
and differentiated-mechanics steering reduces median nonlinear work by about
57% relative to random RRT. This benchmark concerns tip-obstacle branching;
different body-contact histories remain a harder future extension.

### 3.29 Ten-seed dead-end robustness sweep

[`../python/simder_planning/sweep.py`](../python/simder_planning/sweep.py)
defines a versioned, deterministic sweep record containing raw comparisons,
per-method summaries, and Wilson 95% success intervals.
[`../python/run_dead_end_sweep.py`](../python/run_dead_end_sweep.py) runs 11
one-factor cases: the baseline and low/high changes to obstacle radius,
obstacle vertical offset, tip safety distance, field magnitude limit, and goal
bias. All methods retain identical budgets and safety criteria within a case.

For seeds 0--9, mechanics-informed RRT succeeds in 9--10/10 runs across every
case. Its baseline is 9/10 (95% Wilson interval 0.596--0.982) with median 291
continuation steps. Random-actuation RRT succeeds in 7--10/10; its baseline is
10/10 (0.722--1.000) with median 587.5 steps, but low goal bias reduces it to
7/10 (0.397--0.892) and raises median work to 1381 steps. Local goal seeking is
0/10 in all cases (upper confidence bound 0.278), always stopping at the direct
tip-safety dead end.

The mechanics-informed baseline miss is seed 4 at the fixed 300-iteration
budget; successful runs have median final errors near `3.27e-8`. No selected
path violates tip safety, body clearance, or stability. Obstacle-radius and
tip-safety perturbations with the same signed change produce identical results,
as expected: for this spherical tip obstacle they enter feasibility only via
the sum `obstacle_radius + tip_safe_distance`. These are sensitivity checks of
the same effective clearance, not independent statistical evidence.

The sweep supports the infrastructure and the usefulness of global planning,
but ten seeds still give broad confidence intervals and do not prove one RRT
has a higher success probability. The stronger supported claim is that
mechanics-informed proposals consistently reduce median continuation work and
are less sensitive to reduced goal bias in this benchmark.

### 3.30 Body-contact-history branching benchmark

[`../geometry/doubleSphericalObstacleDomain.cpp`](../geometry/doubleSphericalObstacleDomain.cpp)
defines an outer spherical cavity with two excluded spheres. Its closest
boundary query uses stable IDs 0/1/2, and every boundary supplies analytic
clearance, gradient, Hessian, and closest point. The domain and both obstacle
parameter sets are exposed through `MechanicsConfig`, pybind11, and versioned
task JSON.

`double_obstacle_contact_branching_task()` places two radius-`1e-3` obstacles
symmetrically around rod vertex 14. The gap between them is exactly one rod
diameter. Direct actuation to `B=(0,0,0.3)` converges successfully but remains
pinned on contacts `(14,1)` and `(14,2)`, leaving tip error above `1e-2`.
This is an equilibrium-branch dead end rather than a solver failure.

Two independently calibrated four-edge routes first select `B_y=+0.05` or
`-0.05`, inserting boundary 1 or 2. Increasing `|B_y|` releases that contact,
after which the route moves above the gate and returns to the common target.
Both routes finish contact-free and their tips agree within `1.1e-8`. The
regression verifies stable IDs, nonpenetration, tip safety, stability, and JSON
round-trip before relying on any planner result.

[`../python/run_contact_branching_comparison.py`](../python/run_contact_branching_comparison.py)
runs the current planner comparison. For seeds 0--2, mechanics-informed RRT
succeeds 3/3 with median error `1.92471e-8` and 5787 continuation steps;
random-actuation RRT succeeds 1/3 with median error `6.82704e-4` and 7497
steps; local goal seeking succeeds 0/3 with error `6.04901e-3` and 84 steps.
The successful mechanics paths contain a side-contact insertion and release.
These three-seed counts validate the benchmark but are not a robust
statistical comparison. Full failed trees are expensive because they consume
the 1,000-iteration budget, so the multi-seed comparison is intentionally not
part of routine CTest.

### 3.31 Budgeted contact-branching diagnosis

`RRTConfig.maximum_continuation_steps` is a deterministic soft mechanics-work
cap. A nonlinear continuation edge remains atomic and may finish just beyond
the cap, but the planner starts no further tree iteration. `RRTResult` and
`RunRecord` now distinguish `goal_reached`, `iteration_budget`,
`mechanics_work_budget`, and `local_stall`; logs also retain terminal-failure
classes and successful contact-history classes.

The diagnosis exposed a real steering defect. Local least squares treated
insertion as unconstrained, after which `_bounded_actuation()` clipped negative
insertion without recomputing the field command. A near-goal node at
`xi=0.0055` therefore repeatedly moved away from a target calibrated near
`xi=0.005`. Steering now uses a one-constraint active set: if the unconstrained
solution requests retraction, it fixes `dxi=0` and re-solves the three field
components. The terminal connector also has bounded step backtracking for
nonlinear no-progress trials. The previously failing seed 1 then reaches
`1.92489e-9` error in 5,787 continuation steps.

[`../python/run_contact_branching_budget_sweep.py`](../python/run_contact_branching_budget_sweep.py)
runs five cases over seeds 0--9 and checkpoints after each case. It can select
named cases with repeatable `--case` options for deterministic sharding. The
updated mechanics-informed results are:

- 500 iterations: 8/10 successes (95% Wilson interval 0.490--0.943).
- 1,000 iterations with a 15,000-work cap: 8/10 (0.490--0.943).
- 3,000 work: 4/10 (0.168--0.687).
- 7,000 work: 8/10 (0.490--0.943).
- 15,000 work: 8/10 (0.490--0.943).

Successful paths include five boundary-1-only and three mixed histories in the
high-budget cases. Seed 4 reaches `2.10674e-6`, narrowly outside the `2e-6`
goal tolerance, and reports terminal no progress. Seed 5 never reaches the
goal neighborhood and is a search/node-selection miss. Raising work from
7,000 to 15,000 rescues neither. The serial 50-run sweep took about 25 minutes;
linear nearest-node scans over growing trees are now a measured performance
bottleneck, and per-case checkpointing prevents losing completed cases.

### 3.32 Projected feasibility and deterministic nearest index

Tip-space nearest selection now uses an exact deterministic three-dimensional
k-d tree. Rebuilds occur geometrically, newly appended nodes remain in a
linear tail until the next rebuild, and equal distances select the lower node
index exactly as the brute-force implementation did. Tests run both versions
and require identical nodes and iterations. `RRTResult`/`RunRecord` report
nearest query and exact distance-evaluation counts.

Every `RRTNode` also accumulates boundary IDs visited along its ancestry.
Experiment records store node counts and closest goal distance for the
`none`, `boundary_1_only`, `boundary_2_only`, and `mixed` history classes; this
does not yet influence selection.

Terminal steps report the residual after projecting local least squares onto
the feasible unilateral-insertion cone, both as an absolute norm and relative
to the requested tip correction. Seed 4 at the 15,000-work cap ends at error
`2.1067405e-6`. All six terminal trials activate `dxi=0`; the final projected
residual is `1.31671e-7` with ratio `1.0`. Thus its remaining correction is
locally unavailable without insertion retraction. Relaxing the tolerance would
label it successful, but would hide this actuation-feasibility distinction.

Seed 5 uses 906 nearest queries, 148,012 exact distance evaluations, and builds
651 nodes. Its distribution is 541 mixed-history, 79 boundary-2-only, 30
boundary-1-only, and one contact-free node. Closest errors by class are
`0.00607287` mixed, `0.0159705` boundary 2, `0.00847577` boundary 1, and
`0.0177988` none. It never invokes the terminal connector. This confirms a
history-distribution/search failure rather than insufficient terminal
iterations. Seed 4 similarly uses 949 queries and 154,219 distance evaluations
for 697 nodes.

### 3.33 Contact-history-stratified mechanics shortlist

`RRTConfig.use_contact_history_shortlist` optionally partitions the exact
nearest-node query by accumulated contact-history tuple. Each populated class
contributes one exact nearest candidate. The planner evaluates feasible local
steering only for this small set and ranks candidates by projected endpoint
distance to the sampled tip, projected residual ratio, and node index. Separate
deterministic k-d trees serve the history classes, and brute-force mode retains
the same classwise semantics. `RRTResult` and `RunRecord` report candidates
evaluated and selections that intentionally differ from the global nearest.

The first residual-first ranking was rejected: it changed 997/1,000 seed-5
selections and increased final error to `0.0156253`. With projected endpoint
distance primary, seed 5 succeeds at error `1.63615e-7` after 621 iterations
and 8,952 continuation steps. Only 49 selections differ from global nearest;
its successful path has boundary-1-only accumulated history.

The full ten-seed sweep used the unchanged five cases and soft work caps. The
shortlist success counts are 7/10 at 500 iterations, 8/10 at 1,000 iterations,
3/10 at 3,000 work, 7/10 at 7,000 work, and 8/10 at 15,000 work. The matching
global-nearest counts are 8/10, 8/10, 4/10, 8/10, and 8/10. At 15,000 work the
shortlist rescues seeds 4 and 5 but loses seeds 2 and 7. Median candidate count
is 1,551.5 in the high-budget cases, with median 20.5 non-nearest selections
at 1,000 iterations/15,000 work. Total wall time for the five shortlist cases
was about 39 minutes versus about 28 minutes for the earlier global-nearest
sweep on this workstation. The policy is therefore useful as an opt-in branch
diagnostic, but is not a better default.

### 3.34 Deterministic hybrid shortlist triggers

`RRTConfig` now supports mutually exclusive periodic, stagnation, and late-start
hybrid modes in addition to the existing always-shortlist flag. Periodic and
stagnation triggers can request a deterministic consecutive-iteration burst.
Stagnation measures cumulative goal-distance reduction against a configurable
significant-progress threshold, and resets its clock after either significant
progress or a trigger. Logs report the exact number of shortlist queries.
The budget-sweep CLI exposes these controls as `--history-shortlist-period`,
`--history-shortlist-stagnation`, `--history-shortlist-progress`,
`--history-shortlist-burst`, and `--history-shortlist-after`.

Critical-seed tuning used the unchanged 1,000-iteration/15,000-work case.
Single-query periodic cadences 5, 10, and 20 all preserve seeds 2 and 7 but
fail seed 5. Stagnation intervals 25, 50, and 100 with progress threshold
`1e-4` do the same. A period-20, burst-5 policy succeeds on all three: seed 2
finishes in 365 iterations, seed 5 in 303 iterations at error `2.62e-7`, and
seed 7 in 117 iterations. Seed 5 uses 74 shortlist queries and ten non-nearest
selections.

On all ten seeds at 15,000 work, period-20/burst-5 still scores 8/10. It
succeeds on seeds 0, 2, 3, 5, 6, 7, 8, and 9, while seeds 1 and 4 exhaust the
work cap; seed 4 stops narrowly at `2.06327e-6`. Median continuation work drops
to 3,193.5, but total wall time is about 533 seconds because failed seeds still
evaluate multiple local sensitivities. A continuous shortlist fallback at
iteration 300, 400, or 500 does not rescue seed 5, showing that the decisive
history allocation occurs early. These hybrids remain experimental and the
global-nearest selector remains the default.

## 4. Build and run on another workstation

### Required dependencies

The build requires:

- C++17 compiler
- CMake 3.16 or newer
- Eigen3
- BLAS/LAPACK
- OpenGL
- GLU
- GLUT/freeglut
- pybind11 only when `SIMDER_BUILD_PYTHON=ON`

On Ubuntu, suitable packages are typically:

```bash
sudo apt install build-essential cmake libeigen3-dev libopenblas-dev \
    liblapack-dev libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev
```

Package names may differ on other systems.

### Configure and build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON
cmake --build build -j
```

Optional Python module and smoke test:

```bash
python3 -m pip install pybind11
cmake -S . -B build-python -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DBUILD_TESTING=ON -DSIMDER_BUILD_PYTHON=ON
cmake --build build-python -j
ctest --test-dir build-python --output-on-failure
```

### Run the tests

```bash
ctest --test-dir build --output-on-failure
```

To see the static residual test output:

```bash
ctest --test-dir build -V -R zero_load_static_equilibrium
```

### Run the simulator

The legacy mesh paths are still relative to the repository root, so run the
program from that directory:

```bash
./build/simDER option.txt
```

Short headless example:

```bash
./build/simDER option.txt -- render 0 saveData 0 totalTime 0.05
```

Command-line option overrides must appear after `--`.

## 5. Current execution architecture

The current legacy dynamic flow is:

```text
main
  -> load option file and command-line overrides
  -> construct world
  -> construct rod, stepper, and force objects
  -> updateTimeStep
       -> extrapolate dynamic guess
       -> prepare rod geometry
       -> assembleDynamicSystem
            -> assembleStaticSystem
            -> add inertia
            -> add damping
       -> solve banded Newton system with LAPACK dgbsv
       -> update unconstrained rod DOFs
       -> commit time-step state
```

The static diagnostic flow is:

```text
construct world and rod
  -> prepare rod geometry
  -> assembleStaticSystem
  -> compute residual norm over free DOFs
```

The static solve flow is:

```text
capture RodState
  -> evaluateStaticSystem into owned residual/banded Jacobian
  -> copy evaluation into LAPACK workspace
  -> solve J dq = r
  -> for planar barrier contact, backtrack q <- q - alpha dq
       -> reject nonpositive gaps
       -> require residual-merit decrease when contact is involved
       -> restore the exact pre-trial state before each retry
  -> otherwise apply the validated full Newton step
  -> repeat until absolute force tolerance
  -> commit with zero velocity, or restore input state on failure
```

Force classes still write directly into the raw residual and banded Jacobian
owned by `timeStepper`. That storage arrangement is temporary and will not be
sufficient for an augmented contact KKT system.

## 6. Known limitations that remain

### Mechanics and solver

- `PlannerState` captures rod configuration/history, actuation, contact state,
  and solver diagnostics, but it does not capture the legacy dynamic clock.
- Static Newton currently uses the legacy banded Jacobian and has only been
  validated for a small contact-free, zero-load perturbation.
- The planar barrier has feasibility-preserving line search, but the static
  solver does not yet have trust-region control or a dynamic-relaxation
  fallback. Legacy/contact-free modes retain their previously validated full
  Newton trajectory because their approximate tangent did not behave well with
  a strictly monotone residual-merit line search.
- The residual/Jacobian accumulator is still coupled to `timeStepper`.
- The legacy banded LAPACK layout cannot directly represent the augmented KKT
  system; the first KKT layer therefore uses an independent dense matrix. A
  sparse solver will be needed for larger systems.
- Static convergence, velocity convergence, and dynamic convergence are not
  yet reported independently.

### Magnetics

- The backward-compatible legacy model still applies loading over every rod
  edge, and `externalMagneticForce::computeJm()` remains empty.
- The new planning model currently supports only a dipole aligned with the
  terminal tangent; a general local material-frame dipole is not implemented.
- The legacy magnetic force evaluates parts of its state from the previous
  committed configuration, which is unsuitable for a consistent static Newton
  residual.

### Contact and geometry

- The planar slab, analytic spherical shell, and eccentric spherical-obstacle
  domain support the common KKT/barrier path; the legacy triangle mesh still
  uses only its original penalty/barrier force.
- The analytic barrier consumes the tested clearance query, Hessian, and
  candidate detector, but the legacy triangle-mesh contact force still does
  not.
- The nonlinear KKT result owns multipliers and stable candidate IDs, but they
  are exposed to Python planning only through the still-planar
  `MechanicsConfig`; C++ `PlannerState` also works for the analytic shell.
- Contact insertion/release and fixed-active curved continuation work for the
  analytic shell but have not been generalized to triangle-mesh contact.
- Triangle selection uses average distance to triangle vertices rather than an
  exact point-to-triangle closest-point query.
- Contact search is brute-force `O(number_of_rod_vertices * number_of_triangles)`.
- `mesh_1` paths are hard-coded in `elasticRod::readInputMesh()`.
- `mesh_2` is present but unused.
- Tip clearance is detected separately and enforced by the analytic-domain KKT solver,
  and both the Boolean `tip_safe` and continuous `tip_clearance` are exposed
  through `PlannerState` and recorded in experiment paths.

### Planning and experiment infrastructure

- The `proximal_guide` insertion model is an architecture/derivative
  validation model, not a physical material-feed or changing deployed-length
  implementation. The final experimental interpretation of `xi` remains to
  be selected and implemented behind `InsertionModel`.
- Combined `(xi, B)` sensitivity and continuation are implemented for the
  provisional guide, axial-tip magnetics, planar contact, and analytic
  spherical-shell contact.
- The planner is validated on a deterministic reachable-goal regression,
  multi-seed contact-free Task A, a planar contact-release task, spherical-
  shell Task B, a tip-obstacle dead end, a ten-seed one-factor robustness
  sweep, and a two-obstacle body-contact-history branch. It has an exact
  deterministic spatial index and an optional mechanics/history-aware
  shortlist, but no wall-clock timeout or hybrid fallback policy.
- `world::CoutData()` still does not write simulation state or metrics.
- Versioned planning logs exist, but Newton/KKT factorization counts are not
  yet exposed through `MechanicsSession`.

## 7. Exact next milestone

Always-on, low-cadence, stagnation, burst, and late-start shortlist policies are
implemented and measured. Period-20/burst-5 recovers critical seeds 2, 5, and
7, but the full ten-seed result remains 8/10 because seeds 1 and 4 are
displaced. A single shared tree therefore cannot add history exploration
without changing the deterministic global-nearest backbone.

Next, prototype a deterministic dual-frontier planner: keep a global-nearest
backbone whose index contains only backbone nodes, and schedule auxiliary
history-stratified expansions without allowing those nodes to perturb backbone
selection. Use an explicit shared continuation-work allocation and report
local-steering evaluations separately. First test whether it preserves all
eight baseline successes while adding seed 5; do not claim improvement unless
the ten-seed result exceeds 8/10 under the same total work cap. Continue to
defer RRT* and rewiring.

## 8. Planned order after the next milestone

1. ~~Explicit rod state and state capture/restore.~~ Completed.
2. ~~Owned static evaluation and contact-free static Newton solve.~~ Completed.
3. ~~Finite-difference validation of the elastic residual/Jacobian.~~ Completed.
4. ~~Dynamic-relaxation comparison using velocity and static residual.~~ Completed.
5. ~~Axial tip-dipole model and field sensitivity.~~ Completed; general
   material-frame dipole remains deferred.
6. ~~Signed-clearance `ConfinedDomain` geometry interface.~~ Completed.
7. ~~Barrier prediction followed by active-set/KKT contact correction.~~
   Completed for the planar slab and analytic spherical shell; triangle-mesh
   complementarity remains deferred.
8. ~~Reduced-Hessian stability and equilibrium sensitivity.~~ Completed for
   contact-free, planar-contact, and fixed-active spherical-contact changes.
9. ~~Predictor-corrector continuation.~~ Completed for magnetic-field paths,
   including a planar contact-release event and fixed-active curved contact.
10. ~~Explicit `Actuation`/`PlannerState` and insertion-coordinate model.~~
    Completed with a clearly provisional compliant proximal-guide convention
    and combined `(xi, B)` continuation.
11. ~~Narrow pybind11 interface.~~ Completed with an optional high-level
    `MechanicsSession` module and C++/Python rollback tests.
12. ~~Tip-steering facade and Python baseline RRT.~~ Completed with finite-
    difference validation and a deterministic reachable-goal tree.
13. ~~Iterative exact terminal connector.~~ Completed with full-rank/condition
    diagnostics and deterministic success/failure tests.
14. ~~Contact-free Task A and persistent experiment logging.~~ Completed for
    three replayed seeds with compact/full-tree JSON output.
15. ~~Random-actuation RRT and local goal-seeking Task-A comparisons.~~
    Completed with common limits, method labels, aggregates, and exact replay.
16. ~~Planar contact-transition benchmark with continuous clearance logging.~~
    Completed with an initial-contact/release task, three-method comparison,
    safety-margin checks, and exact replay.
17. ~~Analytic curved-domain gap Hessian and curved KKT mechanics validation.~~
    Completed for inner/outer sphere queries, barrier/KKT curvature,
    equilibrium, stability, sensitivity, continuation, and rollback.
18. ~~Serialized spherical-shell Task B and three-method comparison.~~
    Completed with backward-compatible domain configuration, sustained-contact
    enforcement, three replayed seeds per method, and safety checks.
19. ~~Branching/dead-end benchmark design and multi-seed comparison.~~
    Completed for an eccentric forbidden-tip obstacle: both RRT variants
    succeed 3/3, local goal seeking fails 0/3, failure reasons are logged, and
    mechanics-informed steering uses less continuation work.
20. ~~One-factor robustness sweep with more seeds and confidence intervals.~~
    Completed for 11 cases and ten seeds per method, with versioned raw logs,
    Wilson intervals, deterministic replay, and failure-seed analysis.
21. ~~Body-contact-history branching benchmark.~~ Completed with a
    two-obstacle analytic domain, two distinct insertion/release routes, a
    pinned direct equilibrium, and an initial three-seed planner comparison.
22. ~~Budgeted multi-seed robustness study for contact-history branching.~~
    Completed with deterministic soft work caps, ten seeds, Wilson intervals,
    failure classification, active-set steering correction, and checkpointed
    case shards.
23. ~~Projected terminal-feasibility diagnostic and deterministic spatial
    index.~~ Completed with exact brute-force equivalence, query metrics,
    ancestry history classes, and targeted seed-4/5 evidence.
24. ~~Optional contact-history-stratified mechanics shortlist and controlled
    baseline comparison.~~ Completed; it changes the high-budget successful
    seed set but not the 8/10 aggregate, so it remains opt-in.
25. ~~Deterministic low-cadence, stagnation, burst, and late-start hybrid
    shortlist policies.~~ Completed; period-20/burst-5 preserves critical
    seeds 2/5/7 but still swaps two failures in the full ten-seed case.
26. Separate deterministic global-backbone and auxiliary history frontiers
    under one explicit mechanics-work allocation.

Keep RRT*, bidirectional growth, and mechanics-aware node ranking deferred
until the baseline, terminal connector, and benchmark logging are validated.

## 9. Files added or materially changed

Added:

- `.gitignore`
- `CMakeLists.txt`
- `tests/smoke_options.txt`
- `tests/static_equilibrium_test.cpp`
- `tests/rod_state_test.cpp`
- `tests/static_solver_test.cpp`
- `tests/static_failure_options.txt`
- `rodState.h`
- `staticEvaluation.h`
- `equilibriumResult.h`
- `tipMagneticForce.h`
- `tipMagneticForce.cpp`
- `tests/elastic_jacobian_test.cpp`
- `tests/relaxation_comparison_test.cpp`
- `tests/relaxation_options.txt`
- `tests/tip_magnetic_test.cpp`
- `tests/tip_magnetic_options.txt`
- `tests/field_sensitivity_test.cpp`
- `geometry/surfaceProjection.h`
- `geometry/surface.h`
- `geometry/domainQuery.h`
- `geometry/confinedDomain.h`
- `geometry/triangleGeometry.h`
- `geometry/triangleGeometry.cpp`
- `geometry/triangleMeshSurface.h`
- `geometry/triangleMeshSurface.cpp`
- `geometry/planarSlabDomain.h`
- `geometry/planarSlabDomain.cpp`
- `tests/triangle_geometry_test.cpp`
- `tests/planar_slab_test.cpp`
- `contact/contactCandidate.h`
- `contact/contactDetector.h`
- `contact/contactDetector.cpp`
- `tests/contact_detector_test.cpp`
- `contact/barrierPotential.h`
- `contact/barrierPotential.cpp`
- `contact/planarBarrierContactForce.h`
- `contact/planarBarrierContactForce.cpp`
- `tests/barrier_contact_test.cpp`
- `tests/contact_equilibrium_test.cpp`
- `contact/contactKktSystem.h`
- `contact/contactKktSystem.cpp`
- `tests/contact_kkt_test.cpp`
- `contact/contactKktEquilibriumResult.h`
- `tests/contact_kkt_equilibrium_test.cpp`
- `contact/contactConstraintAnalysis.h`
- `contact/contactConstraintAnalysis.cpp`
- `tests/contact_constraint_analysis_test.cpp`
- `tests/contact_field_sensitivity_test.cpp`
- `continuation/fieldContinuation.h`
- `tests/field_continuation_test.cpp`
- `actuation.h`
- `plannerState.h`
- `insertion/insertionModel.h`
- `insertion/proximalGuideInsertionModel.h`
- `insertion/proximalGuideInsertionModel.cpp`
- `continuation/actuationContinuation.h`
- `tests/insertion_model_test.cpp`
- `tests/actuation_continuation_test.cpp`
- `planning/mechanicsSession.h`
- `planning/mechanicsSession.cpp`
- `python/bindings.cpp`
- `python/simder_planning/__init__.py`
- `python/simder_planning/rrt.py`
- `python/simder_planning/task.py`
- `python/simder_planning/experiment.py`
- `python/run_task_a.py`
- `python/run_task_a_comparison.py`
- `python/run_planar_contact_comparison.py`
- `tests/mechanics_session_test.cpp`
- `tests/python_bindings_smoke.py`
- `tests/python_rrt_test.py`
- `tests/python_terminal_connector_test.py`
- `tests/python_task_a_experiment_test.py`
- `tests/python_task_a_comparison_test.py`
- `tests/python_planar_contact_task_test.py`
- `geometry/sphericalShellDomain.h`
- `geometry/sphericalShellDomain.cpp`
- `tests/spherical_shell_test.cpp`
- `tests/spherical_contact_kkt_test.cpp`
- `python/run_spherical_shell_comparison.py`
- `tests/python_spherical_shell_task_test.py`
- `geometry/sphericalObstacleDomain.h`
- `geometry/sphericalObstacleDomain.cpp`
- `tests/spherical_obstacle_test.cpp`
- `python/run_dead_end_comparison.py`
- `tests/python_dead_end_task_test.py`
- `python/simder_planning/sweep.py`
- `python/run_dead_end_sweep.py`
- `tests/python_dead_end_sweep_test.py`
- `geometry/doubleSphericalObstacleDomain.h`
- `geometry/doubleSphericalObstacleDomain.cpp`
- `tests/double_spherical_obstacle_test.cpp`
- `python/run_contact_branching_comparison.py`
- `tests/python_contact_branching_mechanics_test.py`
- `python/run_contact_branching_budget_sweep.py`
- `tests/python_contact_branching_budget_test.py`
- `doc/DEVELOPMENT_HANDOFF.md`

Materially changed:

- `README.md`
- `main.cpp`
- `world.h`
- `world.cpp`
- `setInput.h`
- `setInput.cpp`
- `timeStepper.h`
- `timeStepper.cpp`
- `elasticRod.cpp`
- `externalContactForce.cpp`
- `option.txt`

## 10. GitHub transfer notes

Before committing, inspect the repository status and avoid committing local
build products unless intentionally desired. In particular, the following may
be generated local artifacts:

```text
build/
simDER
datafiles/
example_simulation.png
```

The screenshot is useful as documentation if intentionally included; the
compiled executable and generated data directory normally should not be
committed.

After cloning on the second workstation, start by running the configure, build,
and CTest commands in Section 4. Continue development only after both tests
pass there.
