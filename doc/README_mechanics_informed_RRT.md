# Mechanics-Informed RRT for Magnetic Elastic-Rod Planning in Confined Space

## 1. Purpose

This document is an implementation specification for a simulation study of **mechanics-informed motion planning for a magnetically actuated elastic rod in confined space**.

The central idea is:

> **The planner does not search arbitrary rod shapes. It grows a tree of stable, mechanically valid equilibrium states. RRT provides global exploration; constrained rod mechanics provides local steering and validates every transition.**

The intended first implementation should stay as simple as possible. Do **not** add RRT*, bidirectional planning, learned surrogates, AtlasRRT, complicated adaptive sampling, or other refinements until the basic method is working and benchmarked.

---

# 2. Physical problem

We consider a slender elastic rod inserted into a three-dimensional confined domain

\[
\Omega \subset \mathbb{R}^3
\]

with curved boundary

\[
\partial \Omega.
\]

The framework should not assume that the domain is a conventional tube. Relevant examples include:

- a tubular neighborhood of a spatial curve;
- a volume swept by closed cross-sections;
- a finite-thickness neighborhood of a curved reference surface;
- a hollow spherical shell;
- a more general smooth confined cavity.

The implementation only needs a geometric interface capable of evaluating:

1. signed clearance/gap from the rod body to the boundary;
2. local surface normal at contact;
3. signed clearance of the magnetic tip;
4. optional nearest-surface point.

## 2.1 Body contact

Frictionless body contact with the confining surface is **allowed**.

For a body contact candidate \(i\),

\[
g_i(q) \ge 0,
\qquad
\lambda_i \ge 0,
\qquad
g_i(q)\lambda_i = 0.
\]

Here:

- \(g_i>0\): separated from boundary;
- \(g_i=0\): touching;
- \(\lambda_i>0\): active compressive contact reaction.

Tangential friction is neglected in the first version.

## 2.2 Tip constraint

The magnetized distal tip is not allowed to contact the boundary.

Require

\[
g_{\rm tip}(q) \ge d_{\rm safe} > 0.
\]

Body contact is part of the mechanics. Tip contact is a planning failure.

---

# 3. Actuation

The actuation vector is

\[
u = (\xi,\mathbf B),
\]

where:

- \(\xi\) is the insertion coordinate;
- \(\mathbf B\in\mathbb R^3\) is a spatially uniform magnetic field.

Allow both magnetic-field direction and magnitude to vary:

\[
\|\mathbf B\| \le B_{\max}.
\]

Insertion should be monotone nondecreasing in the first planner:

\[
\Delta \xi \ge 0.
\]

It is allowed to pause:

\[
\Delta \xi = 0.
\]

Retraction is not allowed in the first version.

## 3.1 Important unresolved modeling choice

Do **not** hard-code the exact insertion mechanics until the physical implementation is finalized.

Create an abstraction such as:

```text
InsertionModel
    apply(xi, rod_state)
    boundary_condition(xi)
    derivative_wrt_xi(...)
```

Possible interpretations include:

1. translating a proximal clamp;
2. feeding more material through a fixed guide;
3. changing deployed rod length;
4. a combination of the above.

The global planner does not depend on which interpretation is chosen, but the derivative \(F_\xi\) does.

---

# 4. Magnetic loading

A permanent magnetic dipole is rigidly attached to the distal tip.

Let the dipole in its local material frame be

\[
\boldsymbol\mu_0.
\]

If the tip orientation is \(R_{\rm tip}(q)\), then

\[
\boldsymbol\mu(q)
=
R_{\rm tip}(q)\boldsymbol\mu_0.
\]

For a spatially uniform magnetic field,

\[
U_{\rm mag}(q,\mathbf B)
=
-\boldsymbol\mu(q)^T\mathbf B.
\]

The magnetic torque is

\[
\boldsymbol\tau_{\rm mag}
=
\boldsymbol\mu \times \mathbf B.
\]

For an ideal dipole in a uniform field,

\[
\mathbf F_{\rm mag}=0.
\]

---

# 5. Equilibrium mechanics

Let \(q\in\mathbb R^n\) denote the free rod degrees of freedom.

The magneto-elastic potential is

\[
\Pi(q;u)
=
E_{\rm el}(q;\xi)
-
\boldsymbol\mu(q)^T\mathbf B.
\]

For an active body-contact set \(\mathcal A\), collect the active contact reactions into

\[
\lambda_{\mathcal A}.
\]

Define the augmented equilibrium variable

\[
z
=
\begin{bmatrix}
q\\
\lambda_{\mathcal A}
\end{bmatrix}.
\]

Let

\[
g_{\mathcal A}(q)
\]

be the active gap vector and

\[
G_{\mathcal A}
=
\frac{\partial g_{\mathcal A}}{\partial q}
\]

its Jacobian.

The equilibrium residual is

\[
F_{\mathcal A}(z;u)
=
\begin{bmatrix}
\nabla_q \Pi(q;u)
-
G_{\mathcal A}^{T}\lambda_{\mathcal A}\\[1mm]
g_{\mathcal A}(q)
\end{bmatrix}.
\]

A fixed-contact-mode equilibrium satisfies

\[
F_{\mathcal A}(z;u)=0.
\]

In addition, a physically valid unilateral-contact equilibrium must satisfy:

\[
\lambda_{\mathcal A}\ge 0,
\]

and every inactive contact candidate must satisfy

\[
g_j(q)\ge0.
\]

### Interpretation

\(F_{\mathcal A}=0\) defines the set of candidate equilibria for one contact mode. It does **not** globally define a unique configuration as a function of \(u\). Multiple equilibria may exist for the same \((\xi,\mathbf B)\).

The previous equilibrium determines which local equilibrium branch is followed.

---

# 6. Local equilibrium structure

At a known equilibrium,

\[
F_{\mathcal A}(z_k;u_k)=0,
\]

linearization gives

\[
F_z\,dz + F_u\,du = 0,
\]

where

\[
F_z
=
\frac{\partial F_{\mathcal A}}{\partial z},
\qquad
F_u
=
\frac{\partial F_{\mathcal A}}{\partial u}.
\]

If \(F_z\) is nonsingular,

\[
dz
=
-F_z^{-1}F_u\,du.
\]

This is the local tangent of the equilibrium branch.

For

\[
u=(\xi,B_x,B_y,B_z),
\]

the locally reachable equilibrium family is parameterized by four actuation coordinates even though \(q\) may be high-dimensional.

---

# 7. KKT tangent matrix

Define the constrained Lagrangian

\[
\mathcal L(q,\lambda;u)
=
\Pi(q;u)
-
\lambda_{\mathcal A}^{T}g_{\mathcal A}(q).
\]

Let

\[
H_{\mathcal L}
=
\nabla^2_{qq}\mathcal L.
\]

Then

\[
F_z
=
\begin{bmatrix}
H_{\mathcal L} & -G_{\mathcal A}^{T}\\
G_{\mathcal A} & 0
\end{bmatrix}.
\]

Use the same tangent/KKT operator for:

1. local sensitivity prediction;
2. Newton correction of equilibrium;
3. optionally local steering Jacobian computation.

Do not store the full KKT matrix at every RRT node. Rebuild/factorize it only when a node is expanded, unless profiling shows caching is beneficial.

---

# 8. Stability

Equilibrium does not imply stability.

For a fixed strongly active contact mode, admissible infinitesimal perturbations satisfy

\[
G_{\mathcal A}dq=0.
\]

Let \(Z\) span the nullspace of \(G_{\mathcal A}\):

\[
G_{\mathcal A}Z=0.
\]

Define the reduced Hessian

\[
H_{\rm red}
=
Z^T H_{\mathcal L} Z.
\]

A strict locally stable equilibrium satisfies

\[
H_{\rm red}\succ0.
\]

A useful scalar stability margin is

\[
\sigma_{\rm stab}
=
\lambda_{\min}(H_{\rm red}).
\]

In the first planner, reject or terminate continuation when

\[
\sigma_{\rm stab}\le\epsilon_{\rm stab}.
\]

Do **not** intentionally plan through snap-through or loss-of-stability transitions in Version 1.

At contact activation/release events where a multiplier is exactly zero, the simple equality-constrained reduced-Hessian test is only an approximation to the full unilateral second-order condition. Treat those events separately and evaluate stability once the new active mode is established.

---

# 9. Contact-mode transitions

The active contact set may change during continuation.

## Activation

An inactive contact \(j\) activates when

\[
g_j\rightarrow0.
\]

At the activation event,

\[
g_j=0,
\qquad
\lambda_j=0.
\]

Update

\[
\mathcal A^+
=
\mathcal A\cup\{j\}.
\]

The augmented state gains one multiplier.

## Release

An active contact \(i\) releases when

\[
\lambda_i\rightarrow0.
\]

At release,

\[
g_i=0,
\qquad
\lambda_i=0.
\]

Update

\[
\mathcal A^-
=
\mathcal A\setminus\{i\}.
\]

The augmented state loses one multiplier.

The physical rod configuration \(q\) should remain continuous across regular contact transitions, while its tangent may change.

---

# 10. Predictor-corrector equilibrium continuation

Given a known equilibrium

\[
(z_k,u_k),
\]

and a small proposed actuation increment

\[
\Delta u,
\]

predict

\[
\Delta z_{\rm pred}
=
-F_z^{-1}F_u\,\Delta u,
\]

and

\[
z_{k+1}^{\rm pred}
=
z_k+\Delta z_{\rm pred}.
\]

Then correct at

\[
u_{k+1}=u_k+\Delta u
\]

using Newton or the existing nonlinear equilibrium solver:

\[
F_z^{(j)}\delta z^{(j)}
=
-
F_{\mathcal A}(z^{(j)};u_{k+1}).
\]

Update until

\[
\|F_{\mathcal A}\|
<
\epsilon_{\rm eq}.
\]

During continuation:

1. detect contact activation;
2. detect contact release;
3. update the active set;
4. solve the new fixed-mode equilibrium;
5. check stability;
6. check tip clearance;
7. check actuator limits.

Use adaptive step size if needed.

### Important

Do not solve each new equilibrium from a cold start. The whole method relies on branch-following from the previous equilibrium.

---

# 11. Planning state

Each RRT tree node must correspond to one full mechanically valid equilibrium.

Conceptually:

\[
X_k
=
(q_k,\xi_k,\mathbf B_k).
\]

Attach mechanics metadata:

\[
\mathcal A_k,
\qquad
\lambda_k.
\]

Suggested node structure:

```python
Node:
    q                  # full rod equilibrium state
    xi                 # insertion coordinate
    B                  # 3-vector magnetic field
    active_set         # contact identifiers
    lambdas            # active contact reactions
    tip_position       # cached 3-vector
    J_tip_u            # optional cached 3x4 local steering Jacobian
    stability_margin   # optional cached scalar
    parent_index
    edge_command       # actuation increment or compact edge description
```

Do not store every temporary continuation state for every tree edge during planning.

After a final path is found, rerun high-resolution continuation along only the selected path to reconstruct the full continuous trajectory.

---

# 12. Local tip steering Jacobian

Let

\[
x_t=x_t(q)
\]

be the distal tip position and

\[
J_t
=
\frac{\partial x_t}{\partial q}.
\]

From

\[
dz=-F_z^{-1}F_u\,du,
\]

extract the \(dq\) block and define

\[
dx_t
=
J_{tu}\,du.
\]

Equivalently,

\[
J_{tu}
=
-J_t[I\;0]F_z^{-1}F_u.
\]

Split it as

\[
dx_t
=
J_{t\xi}\,d\xi
+
J_{tB}\,d\mathbf B.
\]

This Jacobian is the mechanics-derived local steering map.

It should be treated as a **proposal model**, not as the final nonlinear truth.

---

# 13. Core mechanics-informed RRT

The first implementation should use a simple RRT structure.

## 13.1 Sampling

Sample a desired tip-space exploration point

\[
x_{\rm rand}
\in
\Omega_{\rm tip},
\]

where

\[
\Omega_{\rm tip}
=
\{
x\in\Omega:
\operatorname{dist}(x,\partial\Omega)\ge d_{\rm safe}
\}.
\]

With small probability \(p_{\rm goal}\),

\[
x_{\rm rand}=x_{\rm goal}.
\]

Do not sample arbitrary rod configurations \(q\).

## 13.2 Node selection

### Baseline Version

Start with the simplest rule:

\[
X_{\rm near}
=
\arg\min_i
\|x_{t,i}-x_{\rm rand}\|.
\]

### Mechanics-informed Version

If the baseline repeatedly selects mechanically poor nodes:

1. find the \(K\) nearest nodes in tip space;
2. rank only those nodes by local mechanical steerability.

Do not make this refinement mandatory before the baseline is tested.

## 13.3 Desired local tip step

Given \(X_{\rm near}\),

\[
d
=
x_{\rm rand}-x_{t,\rm near}.
\]

Choose a local step length \(\delta_x\):

\[
\Delta x_{\rm des}
=
\delta_x
\frac{d}{\|d\|}.
\]

The sample can be global; the actual expansion should remain local.

## 13.4 Local inverse steering

Solve a small optimization in the four actuation increments

\[
\Delta u
=
(\Delta\xi,\Delta\mathbf B).
\]

A simple first form is

\[
\min_{\Delta u}
\left\|
J_{tu}\Delta u
-
\Delta x_{\rm des}
\right\|^2
+
\rho_B\|\Delta\mathbf B\|^2
+
\rho_\xi(\Delta\xi-\Delta\xi_{\rm pref})^2
\]

subject to

\[
0\le\Delta\xi\le\Delta\xi_{\max},
\]

\[
\|\Delta\mathbf B\|\le\Delta B_{\max},
\]

\[
\|\mathbf B+\Delta\mathbf B\|\le B_{\max}.
\]

Do not include the full nonlinear equilibrium equations as constraints in this small optimization.

The linearized equilibrium constraint is already embedded in \(J_{tu}\).

## 13.5 Nonlinear validation

Take the proposed

\[
\Delta u^\star
\]

and perform full predictor-corrector equilibrium continuation.

The actual new tip location will generally differ from the first-order prediction:

\[
x_{t,\rm new}
\neq
x_{t,\rm near}
+
J_{tu}\Delta u^\star.
\]

That is expected.

Only the nonlinear equilibrium result may become a tree node.

## 13.6 Edge validity

A candidate edge is valid only if every intermediate continuation state satisfies:

\[
F_{\mathcal A}=0,
\]

\[
H_{\rm red}\succ0,
\]

\[
g_{\rm body}\ge0,
\]

\[
\lambda\ge0,
\]

\[
g_{\rm tip}\ge d_{\rm safe},
\]

\[
\|\mathbf B\|\le B_{\max}.
\]

If continuation fails, reject the edge.

For Version 1, it is acceptable to reject the entire attempted edge. A later implementation may add the longest valid partial prefix as a new node.

## 13.7 Add node

If the continuation succeeds,

\[
X_{\rm new}
=
(q_{\rm new},\xi_{\rm new},\mathbf B_{\rm new})
\]

is added to the tree.

Store the parent and the actuation increment used to reach it.

Then sample again.

The next iteration may expand any node in the existing tree, not necessarily \(X_{\rm new}\).

This is essential for global exploration.

---

# 14. Optional mechanics-aware node ranking

Only implement this after the baseline RRT works.

For each of the \(K\) nearest tip-space candidates, define

\[
\hat d_i
=
\frac{x_{\rm rand}-x_{t,i}}
{\|x_{\rm rand}-x_{t,i}\|}.
\]

Set

\[
\Delta x_i^{\rm des}
=
\delta_x\hat d_i.
\]

Using the cached or recomputed \(J_{tu}^{(i)}\), solve the small local steering problem and define a score such as

\[
c_i
=
\min_{\Delta u}
\left\|
J_{tu}^{(i)}\Delta u
-
\Delta x_i^{\rm des}
\right\|^2
+
\rho\|W_u\Delta u\|^2.
\]

Choose the node with the smallest score among the geometric shortlist.

The design principle is:

> Geometry defines a local candidate neighborhood; mechanics chooses which equilibrium inside that neighborhood is easiest to steer toward the sample.

Do not rank every node only by mechanical cost, because that may weaken the usual RRT exploration bias.

---

# 15. Terminal goal connection

Ordinary RRT expansion does not need to land exactly at the target.

Once a node enters a goal neighborhood

\[
\|x_t-x_{\rm goal}\|
\le
r_{\rm goal},
\]

attempt a special local terminal connection.

Define

\[
e_k
=
x_{\rm goal}-x_{t,k}.
\]

Solve

\[
J_{tu,k}\Delta u_k
\approx
e_k
\]

subject to actuator limits.

Then perform full nonlinear equilibrium continuation, recompute the true error, relinearize, and repeat.

Terminate successfully when

\[
\|x_t-x_{\rm goal}\|
<
\epsilon_{\rm goal}.
\]

A useful local condition is

\[
\operatorname{rank}(J_{tu})=3.
\]

If the terminal connector fails from one nearby equilibrium, keep the node and continue the RRT. Another nearby equilibrium/contact history may provide a better-conditioned connection to the target.

---

# 16. Why the method is globally exploratory

The local inverse steering is a local optimization, but the entire planner is not a local optimizer.

RRT repeatedly:

1. samples globally;
2. searches the whole existing tree;
3. can expand any previously discovered equilibrium;
4. preserves alternative contact/equilibrium histories.

A dead branch does not stop the tree.

Example:

```text
X0
├── X1
│   ├── X2
│   │   └── X3   # dead end
│   └── Y2
│       └── Y3
│           └── Goal
└── Z1
```

The planner does not need to physically "escape" from \(X_3\) during offline planning. It simply explores a different stored branch that did not enter that dead end.

---

# 17. Minimal algorithm pseudocode

```python
tree = [initial_equilibrium_node]

while not timeout:

    # 1. Global exploration sample
    x_rand = sample_tip_space(goal_bias=True)

    # 2. Choose node to expand
    X_near = select_near_node(tree, x_rand)

    # 3. Define local desired tip motion
    dx_des = local_tip_step(X_near.tip_position, x_rand)

    # 4. Obtain local mechanics
    J_tip_u = get_or_compute_tip_actuation_jacobian(X_near)

    # 5. Cheap local steering proposal
    du = solve_local_steering(J_tip_u, dx_des, constraints)

    if du is None:
        continue

    # 6. Expensive nonlinear mechanics validation
    result = continue_equilibrium(
        X_near,
        du,
        adaptive_step=True,
        handle_contact_changes=True,
        require_stability=True,
        require_tip_clearance=True,
    )

    if not result.success:
        continue

    X_new = result.endpoint
    X_new.parent = X_near
    tree.append(X_new)

    # 7. Exact terminal connection only near the goal
    if distance(X_new.tip_position, x_goal) <= r_goal:

        terminal = connect_to_goal_iteratively(
            X_new,
            x_goal,
            equilibrium_continuation=True,
        )

        if terminal.success:
            return reconstruct_tree_path(terminal.goal_node)

return FAILURE
```

---

# 18. Suggested software architecture

```text
project/
├── geometry/
│   ├── confined_domain.py
│   ├── sphere_shell.py
│   ├── surface_layer.py
│   └── gap_queries.py
│
├── rod/
│   ├── rod_model.py
│   ├── elastic_energy.py
│   ├── magnetic_tip.py
│   └── insertion_model.py
│
├── contact/
│   ├── contact_detection.py
│   ├── active_set.py
│   └── contact_jacobian.py
│
├── mechanics/
│   ├── equilibrium_residual.py
│   ├── equilibrium_solver.py
│   ├── kkt.py
│   ├── stability.py
│   ├── sensitivity.py
│   └── continuation.py
│
├── planning/
│   ├── node.py
│   ├── sampler.py
│   ├── nearest.py
│   ├── local_steer.py
│   ├── terminal_connector.py
│   └── mechanics_rrt.py
│
├── tasks/
│   ├── task_free_space.py
│   ├── task_spherical_shell.py
│   ├── task_surface_layer.py
│   └── task_dead_end.py
│
├── tests/
│   ├── test_equilibrium.py
│   ├── test_contact.py
│   ├── test_sensitivity.py
│   ├── test_continuation.py
│   └── test_planner.py
│
└── scripts/
    ├── run_mechanics_validation.py
    ├── run_planner_benchmark.py
    └── plot_results.py
```

Adapt this structure to the existing codebase rather than forcing a complete rewrite.

---

# 19. Implementation order

Do not implement the entire planner at once.

## Phase 1 — Mechanics sanity checks

Goal: verify that the equilibrium solver is trustworthy.

Required tests:

1. zero-field relaxed rod;
2. contact-free magnetic bending;
3. finite-difference check of magnetic generalized force;
4. single wall contact;
5. multiple contacts;
6. contact activation and release.

### Contact-free pure-moment sanity check

For a straight isotropic rod with an axially aligned tip magnet, uniform magnetic field produces a tip torque with no translational magnetic force.

The free rod should reproduce the expected constant-moment behavior qualitatively.

This is a useful mechanics sanity test before introducing planning.

---

## Phase 2 — Sensitivity validation

At a known equilibrium, compute

\[
J_{tu}.
\]

For random small \(du\):

1. predict

\[
\Delta x_t^{\rm pred}
=
J_{tu}du;
\]

2. solve the full nonlinear equilibrium after applying \(du\);
3. measure

\[
e_{\rm pred}
=
\left\|
x_t^{\rm actual}
-
x_t^{\rm pred}
\right\|.
\]

Repeat over decreasing step sizes.

For smooth fixed-contact-mode cases, verify approximately second-order prediction error:

\[
e_{\rm pred}
=
O(\|du\|^2).
\]

This is an important validation figure for the paper.

---

## Phase 3 — Continuation validation

Test a prescribed smooth actuation path without RRT.

Verify:

- branch following;
- warm-start convergence;
- contact activation;
- contact release;
- continuous \(q\);
- changes in \(\lambda\);
- positive stability margin;
- adaptive step-size behavior.

Plot:

\[
g_i,
\qquad
\lambda_i,
\qquad
\sigma_{\rm stab},
\qquad
\|F\|,
\qquad
x_{\rm tip}.
\]

---

## Phase 4 — Baseline RRT

Implement the simplest planner:

- sample one \(x_{\rm rand}\) per iteration;
- nearest node by tip distance only;
- mechanics-informed local steering via \(J_{tu}\);
- nonlinear continuation validation;
- no mechanics-aware shortlist ranking yet;
- no RRT*;
- no rewiring;
- no batching.

This is the clean baseline method.

---

## Phase 5 — Terminal connector

Add the iterative exact goal connection near the target.

Measure convergence rate and rank/conditioning of \(J_{tu}\).

---

## Phase 6 — Mechanics-aware node selection

Only if needed, add the \(K\)-nearest tip-space shortlist followed by mechanical steerability ranking.

Compare it against pure tip-distance nearest neighbor.

---

# 20. Benchmark tasks

Use tasks that isolate different algorithmic capabilities.

## Task A — Contact-free sanity case

Purpose:

- verify steering;
- verify goal connection;
- benchmark computational cost without active-set changes.

The target should be reachable without wall contact.

---

## Task B — Simple spherical shell

Reference surface:

\[
\|\mathbf x\|=R.
\]

Admissible shell:

\[
R-d_- < \|\mathbf x\| < R+d_+.
\]

Purpose:

- simple curved confinement;
- body contact on inner/outer surfaces;
- easy-to-visualize contact mechanics;
- test surface-following equilibria.

---

## Task C — Curved surface layer

Create a non-spherical curved reference surface and a finite-thickness neighborhood.

Purpose:

- demonstrate that the method is not tube-specific;
- test changing surface normals;
- test several contact-mode transitions.

---

## Task D — Dead-end / branching environment

Construct a confined geometry where:

- moving greedily toward the target leads to a dead end;
- a successful route initially moves sideways or away from the target;
- different body-contact histories create different future reachability.

Purpose:

- demonstrate why a global tree is necessary;
- compare RRT against local goal-seeking control.

This is one of the most important planning demonstrations.

---

## Task E — Narrow mechanically reachable route

Construct a case where only a relatively small region of equilibrium space leads to success.

Purpose:

- assess RRT exploration difficulty;
- test goal bias;
- test optional mechanics-aware node ranking.

---

## Task F — Magnetic-authority sweep

Run the same planning problem over

\[
B_{\max}
\]

or a nondimensional magnetic authority parameter.

Measure:

- success rate;
- planning time;
- number of equilibrium solves;
- final peak field;
- reachable target region.

This can distinguish true mechanical unreachability from planner failure.

---

# 21. Baselines

At minimum compare against:

## Baseline 1 — Local goal-seeking controller

At every step, locally minimize distance to the final target using the same mechanics Jacobian.

No tree. No global random exploration.

Purpose:

> show that local mechanics alone can enter dead ends.

---

## Baseline 2 — Random-actuation RRT

Use the same equilibrium solver and tree structure, but do not use \(J_{tu}\) for steering.

Sample admissible \((\Delta\xi,\Delta B)\) directly.

Purpose:

> isolate the benefit of mechanics-informed local steering.

---

## Baseline 3 — Tip-space RRT + mechanics validation

Optional intermediate baseline.

Use geometric tip direction to choose a simple actuation heuristic, then validate using full mechanics.

Purpose:

> separate the benefit of the equilibrium tangent from the benefit of global tree search.

---

# 22. Metrics

For each task and method, record:

### Planning quality

- success/failure;
- final tip error;
- total insertion;
- peak magnetic field;
- integrated magnetic-field magnitude;
- number of contact-mode transitions.

### Computational cost

- total planning wall time;
- number of RRT iterations;
- number of attempted edges;
- number of successful edges;
- number of failed equilibrium continuations;
- total number of nonlinear equilibrium corrections;
- total Newton iterations;
- total KKT factorizations;
- mean continuation steps per accepted edge.

### Mechanics quality

- maximum equilibrium residual;
- minimum tip clearance;
- minimum stability margin;
- number of active-set corrections;
- local sensitivity prediction error.

These metrics will reveal whether the main bottleneck is RRT exploration or mechanics continuation.

---

# 23. Reproducibility

All planning experiments should support fixed random seeds.

For every benchmark, run multiple seeds, e.g.

```text
seed = 0, 1, 2, ...
```

Do not report only one successful tree.

Store:

```text
config.yaml
random_seed
solver_tolerances
geometry parameters
rod parameters
magnet parameters
planner parameters
summary metrics
final path
```

---

# 24. Logging recommendations

For each attempted edge, log at least:

```text
parent_node_id
sampled_tip_target
desired_tip_step
proposed_du
predicted_tip
actual_tip
prediction_error
equilibrium_success
num_newton_iterations
num_continuation_steps
contact_mode_changes
minimum_tip_clearance
minimum_stability_margin
failure_reason
```

Recommended failure reasons:

```text
NO_LOCAL_STEER
FIELD_LIMIT
EQUILIBRIUM_DIVERGENCE
CONTACT_INCONSISTENCY
TIP_COLLISION
STABILITY_LOSS
STEP_TOO_SMALL
MAX_CONTINUATION_STEPS
```

This will make debugging and ablation studies much easier.

---

# 25. Plots for the first paper-quality results

## Mechanics figures

1. equilibrium branch with contact activation/release;
2. contact gaps \(g_i\);
3. contact forces \(\lambda_i\);
4. stability margin \(\lambda_{\min}(H_{\rm red})\);
5. local sensitivity prediction error vs. step size.

## Planning figures

1. RRT tree in tip space;
2. several representative full rod equilibrium configurations along different branches;
3. failed/dead-end branch vs. successful branch;
4. final equilibrium path;
5. magnetic field \(B(\xi)\);
6. tip position trajectory;
7. body contact locations along the final path.

## Benchmark plots

1. success rate;
2. planning time;
3. number of nonlinear equilibrium solves;
4. mechanics-informed vs. random steering;
5. local goal-seeking vs. global RRT.

---

# 26. Important implementation principles

### Principle 1

Do not sample arbitrary rod configurations.

### Principle 2

Every stored tree node must correspond to a validated stable equilibrium.

### Principle 3

Use the previous equilibrium to follow the same local branch.

### Principle 4

The local Jacobian is a proposal mechanism, not a replacement for the nonlinear mechanics solver.

### Principle 5

Body contact is mechanically admissible and may be useful.

### Principle 6

Tip contact is always forbidden.

### Principle 7

Do not allow instability/snap-through in Version 1.

### Principle 8

Preserve global exploration. Do not let mechanics-aware steering become a purely greedy goal optimizer.

### Principle 9

Measure the number of expensive equilibrium solves. This is likely to be the main computational bottleneck.

### Principle 10

Start simple. Add algorithmic refinements only when a measured failure mode justifies them.

---

# 27. Definition of Done for Version 1

Version 1 is successful when all of the following are true:

- [ ] equilibrium solver passes contact-free and contact tests;
- [ ] sensitivity prediction agrees with finite nonlinear perturbations;
- [ ] predictor-corrector follows a stable branch;
- [ ] contact activation/release is handled automatically;
- [ ] RRT stores complete equilibrium nodes;
- [ ] RRT can grow a tree using mechanics-informed steering;
- [ ] failed branches do not stop exploration;
- [ ] terminal connector reaches a target to tolerance;
- [ ] planner succeeds on at least one geometry where local goal seeking fails;
- [ ] mechanics-informed steering is compared against random-actuation steering;
- [ ] computational statistics are logged;
- [ ] final path can be reconstructed at high resolution.

---

# 28. Questions that must remain configurable

Do not bake these into the architecture:

1. Exact insertion implementation.
2. DER discretization density.
3. Contact candidate representation.
4. RRT tip step \(\delta_x\).
5. Goal bias \(p_{\rm goal}\).
6. Goal neighborhood \(r_{\rm goal}\).
7. Terminal tolerance \(\epsilon_{\rm goal}\).
8. Maximum field \(B_{\max}\).
9. Maximum local field step \(\Delta B_{\max}\).
10. Maximum insertion step \(\Delta\xi_{\max}\).
11. Stability threshold \(\epsilon_{\rm stab}\).
12. Whether to use mechanics-aware shortlist ranking.
13. Number \(K\) of shortlist nodes.
14. Whether to keep a valid partial edge after an attempted extension fails.

These should be exposed as configuration parameters and investigated experimentally.

---

# 29. Core conceptual summary

The mechanics layer defines the equilibrium structure:

\[
F_{\mathcal A}(z;u)=0.
\]

The tangent of that structure is

\[
dz
=
-F_z^{-1}F_u\,du.
\]

The corresponding local tip response is

\[
dx_t
=
J_{tu}\,du.
\]

The planner uses this local response to propose physically meaningful actuation changes, while full nonlinear equilibrium continuation verifies the actual motion.

Thus:

\[
\boxed{
\text{RRT decides where to explore;}
\quad
\text{mechanics decides how an equilibrium can move;}
\quad
\text{nonlinear continuation decides what the rod actually does.}
}
\]

The intended contribution is **not** “RRT plus a rod simulator.” The intended contribution is a global planner that searches through **stable, contact-supported equilibrium histories**, with local steering derived from the mechanics of the equilibrium set.
