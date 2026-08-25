# Development Handoff

Last updated: 2026-08-25

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
the beginning of a tested mechanics-library refactor. It is not yet the full
equilibrium-continuation or RRT implementation.

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

Latest result:

```text
legacy_headless_smoke             Passed
zero_load_static_equilibrium      Passed

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

Force classes still write directly into the raw residual and banded Jacobian
owned by `timeStepper`. That storage arrangement is temporary and will not be
sufficient for an augmented contact KKT system.

## 6. Known limitations that remain

### Mechanics and solver

- There is no `RodState` value type or safe snapshot/restore API.
- There is no static Newton solver yet.
- There is no `StaticEvaluation` object returning residual and Hessian.
- The residual/Jacobian accumulator is still coupled to `timeStepper`.
- The current banded LAPACK layout cannot directly represent the augmented KKT
  system with contact multipliers.
- Static convergence, velocity convergence, and dynamic convergence are not
  yet reported independently.

### Magnetics

- The current magnetic model applies loading over every rod edge rather than
  representing a distal-tip dipole.
- `externalMagneticForce::computeJm()` is empty.
- The legacy magnetic force evaluates parts of its state from the previous
  committed configuration, which is unsuitable for a consistent static Newton
  residual.

### Contact and geometry

- Contact is a penalty/barrier force, not complementarity/KKT contact.
- There is no explicit signed-gap query API.
- Contact multipliers and persistent active-contact identifiers do not exist.
- Activation and release are not handled explicitly.
- Triangle selection uses average distance to triangle vertices rather than an
  exact point-to-triangle closest-point query.
- Contact search is brute-force `O(number_of_rod_vertices * number_of_triangles)`.
- `mesh_1` paths are hard-coded in `elasticRod::readInputMesh()`.
- `mesh_2` is present but unused.
- Tip clearance is not treated separately from body contact.

### Planning and experiment infrastructure

- No insertion coordinate or `InsertionModel` exists.
- No equilibrium sensitivity, stability, or continuation implementation exists.
- No pybind11 module exists yet.
- No RRT or benchmark task implementation exists yet.
- `world::CoutData()` still does not write simulation state or metrics.
- Random-seed and structured experiment logging have not been added.

## 7. Exact next milestone

The next step is to introduce explicit, copyable mechanics state without yet
changing the force formulas.

Suggested initial types:

```cpp
struct Actuation {
    double xi = 0.0;
    Eigen::Vector3d B = Eigen::Vector3d::Zero();
};

struct RodState {
    Eigen::VectorXd q;
    Eigen::VectorXd velocity;
};

struct StaticEvaluation {
    Eigen::VectorXd residual;
    // Sparse Hessian will replace or supplement legacy banded storage.
};
```

Required state operations:

```cpp
RodState captureState() const;
void restoreState(const RodState& state);
```

The acceptance test for this milestone should:

1. capture the initial state;
2. perturb one or more free degrees of freedom;
3. confirm the static residual becomes nonzero;
4. restore the saved state;
5. confirm the static residual returns to zero.

After state capture/restore passes, implement the first contact-free static
Newton solver and test it on a perturbed, unloaded rod. Then compare its final
configuration and static residual with a dynamically relaxed result.

## 8. Planned order after the next milestone

1. Explicit rod state and state capture/restore.
2. Contact-free static Newton solve using the existing static assembly.
3. Finite-difference validation of the elastic residual/Jacobian.
4. Dynamic-relaxation stopping criteria based on velocity and static residual.
5. Tip-dipole magnetic model and SymEngine-generated derivatives.
6. Signed-gap `ConfinedDomain` geometry interface.
7. Barrier prediction followed by active-set/KKT contact correction.
8. Reduced-Hessian stability and equilibrium sensitivity.
9. Predictor-corrector continuation.
10. Narrow pybind11 interface.
11. Python baseline RRT and benchmark tasks.

Do not begin RRT implementation before the equilibrium, derivative, contact,
and continuation validation phases pass.

## 9. Files added or materially changed

Added:

- `CMakeLists.txt`
- `tests/smoke_options.txt`
- `tests/static_equilibrium_test.cpp`
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
