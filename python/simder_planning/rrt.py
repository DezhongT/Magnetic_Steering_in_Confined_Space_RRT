from dataclasses import dataclass
import math
import random
from typing import List, Optional, Sequence, Tuple

import simder


Vector3 = Tuple[float, float, float]


def _vector3(values: Sequence[float], name: str) -> Vector3:
    if len(values) != 3:
        raise ValueError(f"{name} must have three components")
    result = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must be finite")
    return result


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _distance(left: Sequence[float], right: Sequence[float]) -> float:
    return _norm([a - b for a, b in zip(left, right)])


def _squared_distance(left: Sequence[float], right: Sequence[float]) -> float:
    return sum((a - b) * (a - b) for a, b in zip(left, right))


@dataclass(frozen=True)
class _TipIndexNode:
    index: int
    axis: int
    left: Optional["_TipIndexNode"]
    right: Optional["_TipIndexNode"]


class _DeterministicTipIndex:
    """Exact 3-D k-d tree with deterministic ties and a linear insertion tail."""

    def __init__(self, rebuild_minimum: int = 32) -> None:
        self.rebuild_minimum = rebuild_minimum
        self.indices: List[int] = []
        self.indexed_count = 0
        self.root: Optional[_TipIndexNode] = None

    def add(self, index: int) -> None:
        self.indices.append(index)

    @staticmethod
    def _build(
        nodes: Sequence["RRTNode"], indices: List[int], depth: int = 0
    ) -> Optional[_TipIndexNode]:
        if not indices:
            return None
        axis = depth % 3
        indices.sort(key=lambda index: (nodes[index].tip_position[axis], index))
        middle = len(indices) // 2
        return _TipIndexNode(
            index=indices[middle],
            axis=axis,
            left=_DeterministicTipIndex._build(
                nodes, indices[:middle], depth + 1
            ),
            right=_DeterministicTipIndex._build(
                nodes, indices[middle + 1:], depth + 1
            ),
        )

    def _synchronize(self, nodes: Sequence["RRTNode"]) -> None:
        count = len(self.indices)
        rebuild_threshold = max(
            self.rebuild_minimum,
            2 * self.indexed_count,
        )
        if self.indexed_count == 0 or count >= rebuild_threshold:
            self.root = self._build(nodes, self.indices[:])
            self.indexed_count = count

    def nearest(
        self, nodes: Sequence["RRTNode"], sample: Vector3
    ) -> Tuple[int, int]:
        if not self.indices:
            raise ValueError("Cannot query an empty tip index")
        self._synchronize(nodes)
        first_index = self.indices[0]
        best = (
            _squared_distance(nodes[first_index].tip_position, sample),
            first_index,
        )
        evaluations = 0

        def visit(tree_node: Optional[_TipIndexNode]) -> None:
            nonlocal best, evaluations
            if tree_node is None:
                return
            candidate_index = tree_node.index
            candidate = (
                _squared_distance(nodes[candidate_index].tip_position, sample),
                candidate_index,
            )
            evaluations += 1
            if candidate < best:
                best = candidate
            axis = tree_node.axis
            difference = sample[axis] - nodes[candidate_index].tip_position[axis]
            near, far = (
                (tree_node.left, tree_node.right)
                if difference < 0.0
                else (tree_node.right, tree_node.left)
            )
            visit(near)
            if difference * difference <= best[0]:
                visit(far)

        visit(self.root)
        for index in self.indices[self.indexed_count:]:
            candidate = (_squared_distance(nodes[index].tip_position, sample), index)
            evaluations += 1
            if candidate < best:
                best = candidate
        return best[1], evaluations


def _solve_dense(matrix: List[List[float]], rhs: List[float]) -> Optional[List[float]]:
    size = len(rhs)
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1.0e-14:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        for entry in range(column, size + 1):
            augmented[column][entry] /= scale
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(column, size + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][size] for row in range(size)]


def _jacobian_rank_condition(
    jacobian: Sequence[Sequence[float]],
) -> Tuple[int, float]:
    gram = [
        [sum(jacobian[row][column] * jacobian[other][column] for column in range(4))
         for other in range(3)]
        for row in range(3)
    ]
    for _ in range(32):
        row, column = max(
            ((0, 1), (0, 2), (1, 2)),
            key=lambda pair: abs(gram[pair[0]][pair[1]]),
        )
        if abs(gram[row][column]) < 1.0e-24:
            break
        off_diagonal = gram[row][column]
        tau = (
            gram[column][column] - gram[row][row]
        ) / (2.0 * off_diagonal)
        tangent = (
            1.0 if tau >= 0.0 else -1.0
        ) / (abs(tau) + math.sqrt(1.0 + tau * tau))
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        diagonal_row = gram[row][row]
        diagonal_column = gram[column][column]
        for other in range(3):
            if other == row or other == column:
                continue
            row_value = gram[other][row]
            column_value = gram[other][column]
            gram[other][row] = gram[row][other] = (
                cosine * row_value - sine * column_value
            )
            gram[other][column] = gram[column][other] = (
                sine * row_value + cosine * column_value
            )
        gram[row][row] = diagonal_row - tangent * off_diagonal
        gram[column][column] = diagonal_column + tangent * off_diagonal
        gram[row][column] = gram[column][row] = 0.0
    eigenvalues = sorted(max(0.0, gram[index][index]) for index in range(3))
    maximum = eigenvalues[-1]
    if maximum == 0.0:
        return 0, math.inf
    threshold = maximum * 1.0e-12
    rank = sum(value > threshold for value in eigenvalues)
    condition = (
        math.sqrt(maximum / eigenvalues[0])
        if rank == 3 and eigenvalues[0] > 0.0
        else math.inf
    )
    return rank, condition


@dataclass(frozen=True)
class RRTConfig:
    workspace_min: Vector3 = (-0.05, -0.15, -0.099)
    workspace_max: Vector3 = (1.05, 0.15, 0.099)
    goal_bias: float = 0.15
    local_tip_step: float = 0.01
    goal_neighborhood_radius: float = 0.01
    terminal_tolerance: float = 1.0e-5
    maximum_terminal_iterations: int = 10
    maximum_terminal_backtracks: int = 4
    terminal_backtrack_reduction: float = 0.5
    maximum_iterations: int = 200
    maximum_continuation_steps: Optional[int] = None
    random_seed: int = 0
    maximum_insertion_step: float = 2.0e-4
    maximum_field_step: float = 2.0e-2
    maximum_field_magnitude: float = 0.25
    insertion_regularization: float = 1.0e-8
    field_regularization: float = 1.0e-8
    preferred_insertion_step: float = 0.0
    minimum_tip_progress: float = 1.0e-8
    use_spatial_index: bool = True
    spatial_index_rebuild_minimum: int = 32
    use_contact_history_shortlist: bool = False
    contact_history_shortlist_start_iteration: Optional[int] = None
    contact_history_shortlist_period: Optional[int] = None
    contact_history_shortlist_stagnation_iterations: Optional[int] = None
    contact_history_shortlist_stagnation_progress: float = 1.0e-4
    contact_history_shortlist_burst_iterations: int = 1
    use_dual_frontier: bool = False
    dual_frontier_auxiliary_period: int = 4

    def validated(self) -> "RRTConfig":
        lower = _vector3(self.workspace_min, "workspace_min")
        upper = _vector3(self.workspace_max, "workspace_max")
        if any(a >= b for a, b in zip(lower, upper)):
            raise ValueError("workspace_min must be below workspace_max")
        if not 0.0 <= self.goal_bias <= 1.0:
            raise ValueError("goal_bias must lie in [0, 1]")
        positive = (
            self.local_tip_step,
            self.goal_neighborhood_radius,
            self.terminal_tolerance,
            self.maximum_insertion_step,
            self.maximum_field_step,
            self.maximum_field_magnitude,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in positive):
            raise ValueError("RRT step, goal, and actuation bounds must be positive")
        if self.maximum_iterations <= 0 or self.random_seed < 0:
            raise ValueError("RRT iteration count must be positive and seed nonnegative")
        if (
            self.maximum_continuation_steps is not None
            and self.maximum_continuation_steps <= 0
        ):
            raise ValueError("Continuation-step budget must be positive")
        if self.maximum_terminal_iterations <= 0:
            raise ValueError("Terminal iteration count must be positive")
        if self.maximum_terminal_backtracks < 0:
            raise ValueError("Terminal backtrack count must be nonnegative")
        if not 0.0 < self.terminal_backtrack_reduction < 1.0:
            raise ValueError("Terminal backtrack reduction must lie in (0, 1)")
        if self.terminal_tolerance > self.goal_neighborhood_radius:
            raise ValueError("Terminal tolerance must not exceed goal neighborhood")
        if self.insertion_regularization <= 0.0 or self.field_regularization <= 0.0:
            raise ValueError("Local-steering regularization must be positive")
        if not 0.0 <= self.preferred_insertion_step <= self.maximum_insertion_step:
            raise ValueError("preferred insertion must respect insertion bounds")
        if self.minimum_tip_progress < 0.0:
            raise ValueError("minimum_tip_progress must be nonnegative")
        if self.spatial_index_rebuild_minimum <= 0:
            raise ValueError("Spatial-index rebuild minimum must be positive")
        if (
            self.contact_history_shortlist_start_iteration is not None
            and self.contact_history_shortlist_start_iteration <= 0
        ):
            raise ValueError("Contact-history shortlist start must be positive")
        if (
            self.contact_history_shortlist_period is not None
            and self.contact_history_shortlist_period <= 0
        ):
            raise ValueError("Contact-history shortlist period must be positive")
        if (
            self.contact_history_shortlist_stagnation_iterations is not None
            and self.contact_history_shortlist_stagnation_iterations <= 0
        ):
            raise ValueError("Contact-history stagnation interval must be positive")
        if self.contact_history_shortlist_stagnation_progress <= 0.0:
            raise ValueError("Contact-history stagnation progress must be positive")
        if self.contact_history_shortlist_burst_iterations <= 0:
            raise ValueError("Contact-history shortlist burst must be positive")
        if self.dual_frontier_auxiliary_period <= 0:
            raise ValueError("Dual-frontier auxiliary period must be positive")
        shortlist_modes = sum((
            self.use_contact_history_shortlist,
            self.contact_history_shortlist_start_iteration is not None,
            self.contact_history_shortlist_period is not None,
            self.contact_history_shortlist_stagnation_iterations is not None,
        ))
        if shortlist_modes > 1:
            raise ValueError("Only one contact-history shortlist mode may be active")
        if self.use_dual_frontier and shortlist_modes > 0:
            raise ValueError("Dual-frontier and single-tree shortlist modes conflict")
        return self


@dataclass
class RRTNode:
    state: object
    tip_position: Vector3
    parent_index: Optional[int]
    edge_command: Optional[Tuple[float, float, float, float]]
    tip_actuation_jacobian: Optional[Tuple[Tuple[float, ...], ...]] = None
    contact_history: Tuple[int, ...] = ()
    frontier: str = "backbone"


def _state_contact_boundaries(state: object) -> Tuple[int, ...]:
    return tuple(sorted({int(contact[1]) for contact in state.active_contact_ids}))


def _extended_contact_history(parent: RRTNode, state: object) -> Tuple[int, ...]:
    return tuple(sorted(set(parent.contact_history) | set(_state_contact_boundaries(state))))


@dataclass(frozen=True)
class EdgeDiagnostic:
    iteration: int
    parent_index: int
    accepted: bool
    continuation_failed: bool
    sampled_tip: Vector3
    attempted_steps: int
    rejected_steps: int
    contacts_added: int
    contacts_released: int
    reached_path_fraction: float
    failure_reason: str = ""
    frontier: str = "backbone"


@dataclass(frozen=True)
class TerminalStepDiagnostic:
    connector_iteration: int
    parent_index: int
    accepted: bool
    tip_error_before: float
    tip_error_after: float
    jacobian_rank: int
    jacobian_condition: float
    continuation_failed: bool
    attempted_steps: int
    rejected_steps: int
    contacts_added: int
    contacts_released: int
    failure_reason: str = ""
    projected_linear_residual_norm: float = math.inf
    projected_linear_residual_ratio: float = math.inf
    insertion_constraint_active: bool = False


@dataclass(frozen=True)
class TerminalDiagnostic:
    source_index: int
    success: bool
    final_index: int
    initial_tip_error: float
    final_tip_error: float
    steps: Tuple[TerminalStepDiagnostic, ...]


@dataclass(frozen=True)
class _TerminalConnection:
    success: bool
    final_index: int
    diagnostic: TerminalDiagnostic


@dataclass(frozen=True)
class _SteeringProposal:
    target: simder.Actuation
    command: Tuple[float, float, float, float]
    projected_linear_residual_norm: float
    projected_linear_residual_ratio: float
    insertion_constraint_active: bool


@dataclass
class RRTResult:
    success: bool
    nodes: List[RRTNode]
    edge_diagnostics: List[EdgeDiagnostic]
    terminal_diagnostics: List[TerminalDiagnostic]
    goal_index: Optional[int]
    iterations: int
    termination_reason: str = ""
    continuation_attempted_steps: int = 0
    nearest_neighbor_queries: int = 0
    nearest_neighbor_distance_evaluations: int = 0
    shortlist_candidates_evaluated: int = 0
    shortlist_non_nearest_selections: int = 0
    shortlist_triggered_queries: int = 0
    backbone_expansion_queries: int = 0
    auxiliary_expansion_queries: int = 0
    local_steering_evaluations: int = 0

    def path_indices(self) -> List[int]:
        if self.goal_index is None:
            return []
        path = []
        index = self.goal_index
        while index is not None:
            path.append(index)
            index = self.nodes[index].parent_index
        path.reverse()
        return path


class TipSpaceRRT:
    def __init__(
        self,
        mechanics_session: object,
        config: RRTConfig = RRTConfig(),
        continuation_options: Optional[object] = None,
        proposal_policy: str = "mechanics_informed",
    ) -> None:
        self.session = mechanics_session
        self.config = config.validated()
        self.continuation_options = continuation_options or simder.ContinuationOptions()
        if proposal_policy not in ("mechanics_informed", "random_actuation"):
            raise ValueError("Unsupported RRT proposal policy")
        self.proposal_policy = proposal_policy

    def _sample(self, generator: random.Random, goal: Vector3) -> Vector3:
        if generator.random() < self.config.goal_bias:
            return goal
        return tuple(
            generator.uniform(lower, upper)
            for lower, upper in zip(self.config.workspace_min, self.config.workspace_max)
        )

    @staticmethod
    def _nearest(nodes: Sequence[RRTNode], sample: Vector3) -> int:
        return min(
            range(len(nodes)),
            key=lambda index: _distance(nodes[index].tip_position, sample),
        )

    def _desired_step(self, source: Vector3, sample: Vector3) -> Optional[Vector3]:
        direction = tuple(target - current for current, target in zip(source, sample))
        distance = _norm(direction)
        if distance == 0.0:
            return None
        length = min(distance, self.config.local_tip_step)
        return tuple(length * value / distance for value in direction)

    def _steer(
        self,
        state: object,
        jacobian: Sequence[Sequence[float]],
        desired_step: Vector3,
    ) -> Optional[Tuple[simder.Actuation, Tuple[float, float, float, float]]]:
        proposal = self._steer_with_diagnostics(state, jacobian, desired_step)
        if proposal is None:
            return None
        return proposal.target, proposal.command

    def _steer_with_diagnostics(
        self,
        state: object,
        jacobian: Sequence[Sequence[float]],
        desired_step: Vector3,
    ) -> Optional[_SteeringProposal]:
        normal = [[0.0 for _ in range(4)] for _ in range(4)]
        rhs = [0.0 for _ in range(4)]
        for row in range(3):
            for left in range(4):
                rhs[left] += jacobian[row][left] * desired_step[row]
                for right in range(4):
                    normal[left][right] += jacobian[row][left] * jacobian[row][right]
        normal[0][0] += self.config.insertion_regularization
        rhs[0] += (
            self.config.insertion_regularization
            * self.config.preferred_insertion_step
        )
        for index in range(1, 4):
            normal[index][index] += self.config.field_regularization
        increment = _solve_dense(normal, rhs)
        if increment is None:
            return None
        insertion_constraint_active = increment[0] < 0.0

        # Insertion is unilateral. If the unconstrained least-squares step
        # requests retraction, activate dxi=0 and re-solve the remaining field
        # variables; merely clipping dxi would leave field components that are
        # no longer optimal for the feasible problem.
        if insertion_constraint_active:
            field_normal = [[0.0 for _ in range(3)] for _ in range(3)]
            field_rhs = [0.0 for _ in range(3)]
            for row in range(3):
                for left in range(3):
                    field_rhs[left] += jacobian[row][left + 1] * desired_step[row]
                    for right in range(3):
                        field_normal[left][right] += (
                            jacobian[row][left + 1] * jacobian[row][right + 1]
                        )
            for index in range(3):
                field_normal[index][index] += self.config.field_regularization
            field_increment = _solve_dense(field_normal, field_rhs)
            if field_increment is None:
                return None
            increment = [0.0] + field_increment

        bounded = self._bounded_actuation(state, increment)
        if bounded is None:
            return None
        target, command = bounded
        predicted = tuple(
            sum(jacobian[row][column] * command[column] for column in range(4))
            for row in range(3)
        )
        projected_residual = _distance(predicted, desired_step)
        return _SteeringProposal(
            target=target,
            command=command,
            projected_linear_residual_norm=projected_residual,
            projected_linear_residual_ratio=(
                projected_residual / _norm(desired_step)
                if _norm(desired_step) > 0.0 else 0.0
            ),
            insertion_constraint_active=insertion_constraint_active,
        )

    def _select_shortlist_proposal(
        self,
        nodes: Sequence[RRTNode],
        candidates: Sequence[int],
        sample: Vector3,
    ) -> Optional[Tuple[int, _SteeringProposal, int, int]]:
        """Choose the most feasible local prediction, with deterministic ties."""
        best: Optional[Tuple[Tuple[float, float, int], int, _SteeringProposal]] = None
        evaluated = 0
        steering_evaluations = 0
        for parent_index in sorted(candidates):
            parent = nodes[parent_index]
            desired = self._desired_step(parent.tip_position, sample)
            if desired is None:
                continue
            if parent.tip_actuation_jacobian is None:
                steering = self.session.evaluate_local_steering(parent.state)
                steering_evaluations += 1
                parent.tip_actuation_jacobian = tuple(
                    tuple(row) for row in steering.tip_actuation_jacobian
                )
            evaluated += 1
            proposal = self._steer_with_diagnostics(
                parent.state, parent.tip_actuation_jacobian, desired
            )
            if proposal is None:
                continue
            predicted_tip = tuple(
                parent.tip_position[row]
                + sum(
                    parent.tip_actuation_jacobian[row][column]
                    * proposal.command[column]
                    for column in range(4)
                )
                for row in range(3)
            )
            score = (
                _squared_distance(predicted_tip, sample),
                proposal.projected_linear_residual_ratio,
                parent_index,
            )
            candidate = (score, parent_index, proposal)
            if best is None or candidate[0] < best[0]:
                best = candidate
        if best is None:
            return None
        return best[1], best[2], evaluated, steering_evaluations

    def _bounded_actuation(
        self,
        state: object,
        increment: List[float],
    ) -> Optional[Tuple[simder.Actuation, Tuple[float, float, float, float]]]:

        increment[0] = min(
            self.config.maximum_insertion_step,
            max(0.0, increment[0]),
        )
        field_increment_norm = _norm(increment[1:])
        if field_increment_norm > self.config.maximum_field_step:
            scale = self.config.maximum_field_step / field_increment_norm
            increment[1:] = [value * scale for value in increment[1:]]

        current_field = list(state.actuation.field)
        target_field = [
            current + change for current, change in zip(current_field, increment[1:])
        ]
        target_field_norm = _norm(target_field)
        if target_field_norm > self.config.maximum_field_magnitude:
            scale = self.config.maximum_field_magnitude / target_field_norm
            target_field = [value * scale for value in target_field]
            increment[1:] = [
                target - current for target, current in zip(target_field, current_field)
            ]
        if _norm(increment) < 1.0e-14:
            return None

        target = simder.Actuation(
            state.actuation.xi + increment[0],
            target_field,
        )
        return target, tuple(increment)

    def _random_actuation(
        self,
        state: object,
        generator: random.Random,
    ) -> Optional[Tuple[simder.Actuation, Tuple[float, float, float, float]]]:
        increment = [
            generator.uniform(0.0, self.config.maximum_insertion_step),
            generator.uniform(-1.0, 1.0),
            generator.uniform(-1.0, 1.0),
            generator.uniform(-1.0, 1.0),
        ]
        return self._bounded_actuation(state, increment)

    def _validate_problem(
        self,
        initial_state: object,
        goal: Sequence[float],
    ) -> Tuple[Vector3, Vector3]:
        goal_vector = _vector3(goal, "goal")
        initial_tip = _vector3(initial_state.tip_position, "initial tip")
        for name, point in (("initial tip", initial_tip), ("goal", goal_vector)):
            if any(
                value < lower or value > upper
                for value, lower, upper in zip(
                    point, self.config.workspace_min, self.config.workspace_max
                )
            ):
                raise ValueError(f"{name} lies outside the sampling workspace")
        if _norm(initial_state.actuation.field) > self.config.maximum_field_magnitude:
            raise ValueError("Initial field exceeds maximum_field_magnitude")
        return goal_vector, initial_tip

    def _connect_terminal(
        self,
        nodes: List[RRTNode],
        source_index: int,
        goal: Vector3,
    ) -> _TerminalConnection:
        current_index = source_index
        initial_error = _distance(nodes[current_index].tip_position, goal)
        steps: List[TerminalStepDiagnostic] = []
        if initial_error <= self.config.terminal_tolerance:
            diagnostic = TerminalDiagnostic(
                source_index, True, current_index, initial_error, initial_error, tuple()
            )
            return _TerminalConnection(True, current_index, diagnostic)

        for connector_iteration in range(1, self.config.maximum_terminal_iterations + 1):
            current = nodes[current_index]
            error_before = _distance(current.tip_position, goal)
            steering = self.session.evaluate_local_steering(current.state)
            jacobian = tuple(tuple(row) for row in steering.tip_actuation_jacobian)
            current.tip_actuation_jacobian = jacobian
            rank, condition = _jacobian_rank_condition(jacobian)
            if rank < 3:
                steps.append(TerminalStepDiagnostic(
                    connector_iteration=connector_iteration,
                    parent_index=current_index,
                    accepted=False,
                    tip_error_before=error_before,
                    tip_error_after=error_before,
                    jacobian_rank=rank,
                    jacobian_condition=condition,
                    continuation_failed=False,
                    attempted_steps=0,
                    rejected_steps=0,
                    contacts_added=0,
                    contacts_released=0,
                ))
                break
            desired = self._desired_step(current.tip_position, goal)
            proposal = self._steer_with_diagnostics(
                current.state, jacobian, desired
            )
            if proposal is None:
                steps.append(TerminalStepDiagnostic(
                    connector_iteration=connector_iteration,
                    parent_index=current_index,
                    accepted=False,
                    tip_error_before=error_before,
                    tip_error_after=error_before,
                    jacobian_rank=rank,
                    jacobian_condition=condition,
                    continuation_failed=False,
                    attempted_steps=0,
                    rejected_steps=0,
                    contacts_added=0,
                    contacts_released=0,
                ))
                break
            target, command = proposal.target, proposal.command
            edge = self.session.attempt_continuation(
                current.state, target, self.continuation_options
            )
            if not edge.success:
                steps.append(TerminalStepDiagnostic(
                    connector_iteration=connector_iteration,
                    parent_index=current_index,
                    accepted=False,
                    tip_error_before=error_before,
                    tip_error_after=error_before,
                    jacobian_rank=rank,
                    jacobian_condition=condition,
                    continuation_failed=True,
                    attempted_steps=edge.attempted_steps,
                    rejected_steps=edge.rejected_steps,
                    contacts_added=edge.contacts_added,
                    contacts_released=edge.contacts_released,
                    failure_reason=edge.failure_reason,
                    projected_linear_residual_norm=(
                        proposal.projected_linear_residual_norm
                    ),
                    projected_linear_residual_ratio=(
                        proposal.projected_linear_residual_ratio
                    ),
                    insertion_constraint_active=(
                        proposal.insertion_constraint_active
                    ),
                ))
                break
            tip = _vector3(edge.state.tip_position, "terminal tip")
            error_after = _distance(tip, goal)
            accepted = (
                error_after < error_before
                and _distance(tip, current.tip_position)
                >= self.config.minimum_tip_progress
            )
            steps.append(TerminalStepDiagnostic(
                connector_iteration=connector_iteration,
                parent_index=current_index,
                accepted=accepted,
                tip_error_before=error_before,
                tip_error_after=error_after,
                jacobian_rank=rank,
                jacobian_condition=condition,
                continuation_failed=False,
                attempted_steps=edge.attempted_steps,
                rejected_steps=edge.rejected_steps,
                contacts_added=edge.contacts_added,
                contacts_released=edge.contacts_released,
                projected_linear_residual_norm=(
                    proposal.projected_linear_residual_norm
                ),
                projected_linear_residual_ratio=(
                    proposal.projected_linear_residual_ratio
                ),
                insertion_constraint_active=proposal.insertion_constraint_active,
            ))
            backtrack = 0
            while not accepted and backtrack < self.config.maximum_terminal_backtracks:
                backtrack += 1
                desired = tuple(
                    self.config.terminal_backtrack_reduction * value
                    for value in desired
                )
                proposal = self._steer_with_diagnostics(
                    current.state, jacobian, desired
                )
                if proposal is None:
                    break
                target, command = proposal.target, proposal.command
                edge = self.session.attempt_continuation(
                    current.state, target, self.continuation_options
                )
                if edge.success:
                    tip = _vector3(edge.state.tip_position, "terminal tip")
                    error_after = _distance(tip, goal)
                    accepted = (
                        error_after < error_before
                        and _distance(tip, current.tip_position)
                        >= self.config.minimum_tip_progress
                    )
                else:
                    tip = current.tip_position
                    error_after = error_before
                    accepted = False
                steps.append(TerminalStepDiagnostic(
                    connector_iteration=connector_iteration,
                    parent_index=current_index,
                    accepted=accepted,
                    tip_error_before=error_before,
                    tip_error_after=error_after,
                    jacobian_rank=rank,
                    jacobian_condition=condition,
                    continuation_failed=not edge.success,
                    attempted_steps=edge.attempted_steps,
                    rejected_steps=edge.rejected_steps,
                    contacts_added=edge.contacts_added,
                    contacts_released=edge.contacts_released,
                    failure_reason=edge.failure_reason,
                    projected_linear_residual_norm=(
                        proposal.projected_linear_residual_norm
                    ),
                    projected_linear_residual_ratio=(
                        proposal.projected_linear_residual_ratio
                    ),
                    insertion_constraint_active=(
                        proposal.insertion_constraint_active
                    ),
                ))
            if not accepted:
                break
            nodes.append(RRTNode(
                edge.state,
                tip,
                current_index,
                command,
                contact_history=_extended_contact_history(current, edge.state),
                frontier=current.frontier,
            ))
            current_index = len(nodes) - 1
            if error_after <= self.config.terminal_tolerance:
                diagnostic = TerminalDiagnostic(
                    source_index, True, current_index,
                    initial_error, error_after, tuple(steps)
                )
                return _TerminalConnection(True, current_index, diagnostic)

        final_error = _distance(nodes[current_index].tip_position, goal)
        diagnostic = TerminalDiagnostic(
            source_index, False, current_index,
            initial_error, final_error, tuple(steps)
        )
        return _TerminalConnection(False, current_index, diagnostic)

    def _plan_dual_frontier(
        self, initial_state: object, goal: Sequence[float]
    ) -> RRTResult:
        if self.proposal_policy != "mechanics_informed":
            raise ValueError("Dual-frontier planning requires mechanics steering")
        goal_vector, initial_tip = self._validate_problem(initial_state, goal)
        backbone_generator = random.Random(self.config.random_seed)
        auxiliary_generator = random.Random(
            self.config.random_seed ^ 0x9E3779B9
        )
        nodes = [RRTNode(
            initial_state,
            initial_tip,
            None,
            None,
            contact_history=_state_contact_boundaries(initial_state),
            frontier="backbone",
        )]
        diagnostics: List[EdgeDiagnostic] = []
        terminal_diagnostics: List[TerminalDiagnostic] = []
        continuation_work = 0
        local_steering_evaluations = 0
        nearest_queries = 0
        nearest_evaluations = 0
        shortlist_candidates = 0
        shortlist_non_nearest = 0
        shortlist_queries = 0
        backbone_queries = 0
        auxiliary_queries = 0
        backbone_index = _DeterministicTipIndex(
            self.config.spatial_index_rebuild_minimum
        )
        backbone_indices: List[int] = []
        history_indices = {}
        registered_node_count = 0

        def register_new_nodes() -> None:
            nonlocal registered_node_count
            while registered_node_count < len(nodes):
                index = registered_node_count
                node = nodes[index]
                if node.frontier == "backbone":
                    backbone_index.add(index)
                    backbone_indices.append(index)
                if node.contact_history not in history_indices:
                    history_indices[node.contact_history] = _DeterministicTipIndex(
                        self.config.spatial_index_rebuild_minimum
                    )
                history_indices[node.contact_history].add(index)
                registered_node_count += 1

        register_new_nodes()

        def budget_exhausted() -> bool:
            return (
                self.config.maximum_continuation_steps is not None
                and continuation_work >= self.config.maximum_continuation_steps
            )

        def make_result(
            success: bool,
            goal_index: Optional[int],
            iterations: int,
            reason: str,
        ) -> RRTResult:
            return RRTResult(
                success,
                nodes,
                diagnostics,
                terminal_diagnostics,
                goal_index,
                iterations,
                reason,
                continuation_work,
                nearest_queries,
                nearest_evaluations,
                shortlist_candidates,
                shortlist_non_nearest,
                shortlist_queries,
                backbone_queries,
                auxiliary_queries,
                local_steering_evaluations,
            )

        def attempt_expansion(
            sample: Vector3,
            iteration: int,
            frontier: str,
            use_shortlist: bool,
        ) -> Optional[_TerminalConnection]:
            nonlocal continuation_work
            nonlocal nearest_queries, nearest_evaluations
            nonlocal shortlist_candidates, shortlist_non_nearest
            nonlocal shortlist_queries, backbone_queries, auxiliary_queries
            nonlocal local_steering_evaluations

            if frontier == "backbone":
                backbone_queries += 1
            else:
                auxiliary_queries += 1

            shortlist_proposal: Optional[_SteeringProposal] = None
            if use_shortlist:
                shortlist_queries += 1
                if self.config.use_spatial_index:
                    candidates = []
                    evaluations = 0
                    for history in sorted(history_indices):
                        candidate, history_evaluations = history_indices[
                            history
                        ].nearest(nodes, sample)
                        candidates.append(candidate)
                        evaluations += history_evaluations
                else:
                    nearest_by_history = {}
                    for index, node in enumerate(nodes):
                        candidate = (
                            _squared_distance(node.tip_position, sample), index
                        )
                        previous = nearest_by_history.get(node.contact_history)
                        if previous is None or candidate < previous:
                            nearest_by_history[node.contact_history] = candidate
                    candidates = [
                        candidate[1]
                        for _, candidate in sorted(nearest_by_history.items())
                    ]
                    evaluations = len(nodes)
                global_nearest = min(
                    candidates,
                    key=lambda index: (
                        _squared_distance(nodes[index].tip_position, sample), index
                    ),
                )
                selection = self._select_shortlist_proposal(
                    nodes, candidates, sample
                )
                nearest_queries += 1
                nearest_evaluations += evaluations
                if selection is None:
                    return None
                (
                    parent_index,
                    shortlist_proposal,
                    evaluated,
                    new_steering_evaluations,
                ) = selection
                shortlist_candidates += evaluated
                local_steering_evaluations += new_steering_evaluations
                shortlist_non_nearest += int(parent_index != global_nearest)
            else:
                if self.config.use_spatial_index:
                    parent_index, evaluations = backbone_index.nearest(
                        nodes, sample
                    )
                else:
                    evaluations = len(backbone_indices)
                    parent_index = min(
                        backbone_indices,
                        key=lambda index: (
                            _squared_distance(nodes[index].tip_position, sample),
                            index,
                        ),
                    )
                nearest_queries += 1
                nearest_evaluations += evaluations

            parent = nodes[parent_index]
            desired = self._desired_step(parent.tip_position, sample)
            if desired is None:
                return None
            if shortlist_proposal is not None:
                proposal = shortlist_proposal.target, shortlist_proposal.command
            else:
                if parent.tip_actuation_jacobian is None:
                    steering = self.session.evaluate_local_steering(parent.state)
                    local_steering_evaluations += 1
                    parent.tip_actuation_jacobian = tuple(
                        tuple(row) for row in steering.tip_actuation_jacobian
                    )
                proposal = self._steer(
                    parent.state, parent.tip_actuation_jacobian, desired
                )
                if proposal is None:
                    return None
            target, command = proposal
            edge = self.session.attempt_continuation(
                parent.state, target, self.continuation_options
            )
            continuation_work += edge.attempted_steps
            accepted = edge.success
            if accepted:
                tip = _vector3(edge.state.tip_position, "continued tip")
                accepted = (
                    _distance(tip, parent.tip_position)
                    >= self.config.minimum_tip_progress
                )
            diagnostics.append(EdgeDiagnostic(
                iteration,
                parent_index,
                accepted,
                not edge.success,
                sample,
                edge.attempted_steps,
                edge.rejected_steps,
                edge.contacts_added,
                edge.contacts_released,
                edge.reached_path_fraction,
                edge.failure_reason,
                frontier,
            ))
            if not accepted:
                return None
            nodes.append(RRTNode(
                edge.state,
                tip,
                parent_index,
                command,
                contact_history=_extended_contact_history(parent, edge.state),
                frontier=frontier,
            ))
            register_new_nodes()
            if _distance(tip, goal_vector) > self.config.goal_neighborhood_radius:
                return None
            terminal = self._connect_terminal(
                nodes, len(nodes) - 1, goal_vector
            )
            register_new_nodes()
            terminal_diagnostics.append(terminal.diagnostic)
            continuation_work += sum(
                step.attempted_steps for step in terminal.diagnostic.steps
            )
            local_steering_evaluations += len({
                step.connector_iteration for step in terminal.diagnostic.steps
            })
            return terminal if terminal.success else None

        if _distance(initial_tip, goal_vector) <= self.config.goal_neighborhood_radius:
            terminal = self._connect_terminal(nodes, 0, goal_vector)
            register_new_nodes()
            terminal_diagnostics.append(terminal.diagnostic)
            continuation_work += sum(
                step.attempted_steps for step in terminal.diagnostic.steps
            )
            local_steering_evaluations += len({
                step.connector_iteration for step in terminal.diagnostic.steps
            })
            if terminal.success:
                return make_result(True, terminal.final_index, 0, "goal_reached")
            if budget_exhausted():
                return make_result(False, None, 0, "mechanics_work_budget")

        for iteration in range(1, self.config.maximum_iterations + 1):
            if budget_exhausted():
                return make_result(
                    False, None, iteration - 1, "mechanics_work_budget"
                )
            backbone_terminal = attempt_expansion(
                self._sample(backbone_generator, goal_vector),
                iteration,
                "backbone",
                False,
            )
            if backbone_terminal is not None:
                return make_result(
                    True, backbone_terminal.final_index, iteration, "goal_reached"
                )
            if budget_exhausted():
                return make_result(
                    False, None, iteration, "mechanics_work_budget"
                )
            if iteration % self.config.dual_frontier_auxiliary_period == 0:
                auxiliary_terminal = attempt_expansion(
                    self._sample(auxiliary_generator, goal_vector),
                    iteration,
                    "auxiliary",
                    True,
                )
                if auxiliary_terminal is not None:
                    return make_result(
                        True,
                        auxiliary_terminal.final_index,
                        iteration,
                        "goal_reached",
                    )

        return make_result(
            False, None, self.config.maximum_iterations, "iteration_budget"
        )

    def plan(self, initial_state: object, goal: Sequence[float]) -> RRTResult:
        if self.config.use_dual_frontier:
            return self._plan_dual_frontier(initial_state, goal)
        goal_vector, initial_tip = self._validate_problem(initial_state, goal)
        generator = random.Random(self.config.random_seed)
        nodes = [RRTNode(
            initial_state,
            initial_tip,
            None,
            None,
            contact_history=_state_contact_boundaries(initial_state),
        )]
        diagnostics: List[EdgeDiagnostic] = []
        terminal_diagnostics: List[TerminalDiagnostic] = []
        continuation_work = 0
        nearest_queries = 0
        nearest_evaluations = 0
        shortlist_candidates = 0
        shortlist_non_nearest = 0
        shortlist_queries = 0
        local_steering_evaluations = 0
        last_significant_progress_iteration = 0
        last_significant_progress_distance = _distance(initial_tip, goal_vector)
        last_stagnation_trigger_iteration = 0
        shortlist_burst_remaining = 0
        tip_index = _DeterministicTipIndex(
            self.config.spatial_index_rebuild_minimum
        )
        history_indices = {}
        registered_node_count = 0

        def register_new_nodes() -> None:
            nonlocal registered_node_count
            while registered_node_count < len(nodes):
                index = registered_node_count
                tip_index.add(index)
                history = nodes[index].contact_history
                if history not in history_indices:
                    history_indices[history] = _DeterministicTipIndex(
                        self.config.spatial_index_rebuild_minimum
                    )
                history_indices[history].add(index)
                registered_node_count += 1

        register_new_nodes()

        def budget_exhausted() -> bool:
            return (
                self.config.maximum_continuation_steps is not None
                and continuation_work >= self.config.maximum_continuation_steps
            )

        if _distance(initial_tip, goal_vector) <= self.config.goal_neighborhood_radius:
            terminal = self._connect_terminal(nodes, 0, goal_vector)
            register_new_nodes()
            terminal_diagnostics.append(terminal.diagnostic)
            continuation_work += sum(
                step.attempted_steps for step in terminal.diagnostic.steps
            )
            local_steering_evaluations += len({
                step.connector_iteration for step in terminal.diagnostic.steps
            })
            if terminal.success:
                return RRTResult(
                    True, nodes, diagnostics, terminal_diagnostics,
                    terminal.final_index, 0, "goal_reached", continuation_work,
                    nearest_queries, nearest_evaluations,
                    shortlist_candidates, shortlist_non_nearest,
                    shortlist_queries, nearest_queries, 0,
                    local_steering_evaluations
                )
            if budget_exhausted():
                return RRTResult(
                    False, nodes, diagnostics, terminal_diagnostics,
                    None, 0, "mechanics_work_budget", continuation_work,
                    nearest_queries, nearest_evaluations,
                    shortlist_candidates, shortlist_non_nearest,
                    shortlist_queries, nearest_queries, 0,
                    local_steering_evaluations
                )

        for iteration in range(1, self.config.maximum_iterations + 1):
            if budget_exhausted():
                return RRTResult(
                    False, nodes, diagnostics, terminal_diagnostics,
                    None, iteration - 1, "mechanics_work_budget",
                    continuation_work, nearest_queries, nearest_evaluations,
                    shortlist_candidates, shortlist_non_nearest,
                    shortlist_queries, nearest_queries, 0,
                    local_steering_evaluations
                )
            sample = self._sample(generator, goal_vector)
            periodic_trigger = (
                shortlist_burst_remaining == 0
                and
                self.config.contact_history_shortlist_period is not None
                and iteration % self.config.contact_history_shortlist_period == 0
            )
            stagnation_anchor = max(
                last_significant_progress_iteration,
                last_stagnation_trigger_iteration,
            )
            stagnation_trigger = (
                shortlist_burst_remaining == 0
                and
                self.config.contact_history_shortlist_stagnation_iterations
                is not None
                and iteration - stagnation_anchor
                >= self.config.contact_history_shortlist_stagnation_iterations
            )
            hybrid_trigger = self.proposal_policy == "mechanics_informed" and (
                periodic_trigger or stagnation_trigger
            )
            if hybrid_trigger:
                shortlist_burst_remaining = (
                    self.config.contact_history_shortlist_burst_iterations
                )
            shortlisted = self.proposal_policy == "mechanics_informed" and (
                self.config.use_contact_history_shortlist
                or (
                    self.config.contact_history_shortlist_start_iteration
                    is not None
                    and iteration
                    >= self.config.contact_history_shortlist_start_iteration
                )
                or shortlist_burst_remaining > 0
            )
            shortlist_proposal: Optional[_SteeringProposal] = None
            if shortlisted:
                shortlist_queries += 1
                if stagnation_trigger:
                    last_stagnation_trigger_iteration = iteration
                if shortlist_burst_remaining > 0:
                    shortlist_burst_remaining -= 1
                if self.config.use_spatial_index:
                    candidates = []
                    evaluations = 0
                    for history in sorted(history_indices):
                        candidate, history_evaluations = history_indices[
                            history
                        ].nearest(nodes, sample)
                        candidates.append(candidate)
                        evaluations += history_evaluations
                else:
                    nearest_by_history = {}
                    for index, node in enumerate(nodes):
                        candidate = (
                            _squared_distance(node.tip_position, sample), index
                        )
                        previous = nearest_by_history.get(node.contact_history)
                        if previous is None or candidate < previous:
                            nearest_by_history[node.contact_history] = candidate
                    candidates = [
                        candidate[1]
                        for _, candidate in sorted(nearest_by_history.items())
                    ]
                    evaluations = len(nodes)
                global_nearest = min(
                    candidates,
                    key=lambda index: (
                        _squared_distance(nodes[index].tip_position, sample),
                        index,
                    ),
                )
                selection = self._select_shortlist_proposal(
                    nodes, candidates, sample
                )
                nearest_queries += 1
                nearest_evaluations += evaluations
                if selection is None:
                    continue
                (
                    parent_index,
                    shortlist_proposal,
                    evaluated,
                    new_steering_evaluations,
                ) = selection
                shortlist_candidates += evaluated
                local_steering_evaluations += new_steering_evaluations
                shortlist_non_nearest += int(parent_index != global_nearest)
            else:
                if self.config.use_spatial_index:
                    parent_index, evaluations = tip_index.nearest(nodes, sample)
                else:
                    parent_index = self._nearest(nodes, sample)
                    evaluations = len(nodes)
                nearest_queries += 1
                nearest_evaluations += evaluations
            parent = nodes[parent_index]
            desired = self._desired_step(parent.tip_position, sample)
            if desired is None:
                continue
            if self.proposal_policy == "mechanics_informed":
                if shortlist_proposal is not None:
                    proposal = (
                        shortlist_proposal.target,
                        shortlist_proposal.command,
                    )
                else:
                    if parent.tip_actuation_jacobian is None:
                        steering = self.session.evaluate_local_steering(parent.state)
                        local_steering_evaluations += 1
                        parent.tip_actuation_jacobian = tuple(
                            tuple(row) for row in steering.tip_actuation_jacobian
                        )
                    proposal = self._steer(
                        parent.state, parent.tip_actuation_jacobian, desired
                    )
            else:
                proposal = self._random_actuation(parent.state, generator)
            if proposal is None:
                continue
            target, command = proposal
            edge = self.session.attempt_continuation(
                parent.state, target, self.continuation_options
            )
            continuation_work += edge.attempted_steps
            accepted = edge.success
            if accepted:
                tip = _vector3(edge.state.tip_position, "continued tip")
                accepted = (
                    _distance(tip, parent.tip_position)
                    >= self.config.minimum_tip_progress
                )
            diagnostics.append(
                EdgeDiagnostic(
                    iteration,
                    parent_index,
                    accepted,
                    not edge.success,
                    sample,
                    edge.attempted_steps,
                    edge.rejected_steps,
                    edge.contacts_added,
                    edge.contacts_released,
                    edge.reached_path_fraction,
                    edge.failure_reason,
                )
            )
            if not accepted:
                continue
            node = RRTNode(
                edge.state,
                tip,
                parent_index,
                command,
                contact_history=_extended_contact_history(parent, edge.state),
            )
            nodes.append(node)
            register_new_nodes()
            goal_distance = _distance(tip, goal_vector)
            if (
                goal_distance
                <= last_significant_progress_distance
                - self.config.contact_history_shortlist_stagnation_progress
            ):
                last_significant_progress_distance = goal_distance
                last_significant_progress_iteration = iteration
            if _distance(tip, goal_vector) <= self.config.goal_neighborhood_radius:
                terminal = self._connect_terminal(
                    nodes, len(nodes) - 1, goal_vector
                )
                register_new_nodes()
                terminal_diagnostics.append(terminal.diagnostic)
                continuation_work += sum(
                    step.attempted_steps for step in terminal.diagnostic.steps
                )
                local_steering_evaluations += len({
                    step.connector_iteration for step in terminal.diagnostic.steps
                })
                if terminal.success:
                    return RRTResult(
                        True, nodes, diagnostics, terminal_diagnostics,
                        terminal.final_index, iteration, "goal_reached",
                        continuation_work, nearest_queries, nearest_evaluations,
                        shortlist_candidates, shortlist_non_nearest,
                        shortlist_queries, nearest_queries, 0,
                        local_steering_evaluations
                    )
                if budget_exhausted():
                    return RRTResult(
                        False, nodes, diagnostics, terminal_diagnostics,
                        None, iteration, "mechanics_work_budget",
                        continuation_work, nearest_queries, nearest_evaluations,
                        shortlist_candidates, shortlist_non_nearest,
                        shortlist_queries, nearest_queries, 0,
                        local_steering_evaluations
                    )

        return RRTResult(
            False,
            nodes,
            diagnostics,
            terminal_diagnostics,
            None,
            self.config.maximum_iterations,
            "iteration_budget",
            continuation_work,
            nearest_queries,
            nearest_evaluations,
            shortlist_candidates,
            shortlist_non_nearest,
            shortlist_queries,
            nearest_queries,
            0,
            local_steering_evaluations,
        )


class LocalGoalSeekingPlanner(TipSpaceRRT):
    def __init__(
        self,
        mechanics_session: object,
        config: RRTConfig = RRTConfig(),
        continuation_options: Optional[object] = None,
    ) -> None:
        super().__init__(
            mechanics_session,
            config,
            continuation_options,
            proposal_policy="mechanics_informed",
        )

    def plan(self, initial_state: object, goal: Sequence[float]) -> RRTResult:
        goal_vector, initial_tip = self._validate_problem(initial_state, goal)
        nodes = [RRTNode(
            initial_state,
            initial_tip,
            None,
            None,
            contact_history=_state_contact_boundaries(initial_state),
        )]
        diagnostics: List[EdgeDiagnostic] = []
        terminal_diagnostics: List[TerminalDiagnostic] = []
        continuation_work = 0
        local_steering_evaluations = 0

        def budget_exhausted() -> bool:
            return (
                self.config.maximum_continuation_steps is not None
                and continuation_work >= self.config.maximum_continuation_steps
            )

        if _distance(initial_tip, goal_vector) <= self.config.goal_neighborhood_radius:
            terminal = self._connect_terminal(nodes, 0, goal_vector)
            terminal_diagnostics.append(terminal.diagnostic)
            continuation_work += sum(
                step.attempted_steps for step in terminal.diagnostic.steps
            )
            local_steering_evaluations += len({
                step.connector_iteration for step in terminal.diagnostic.steps
            })
            if terminal.success:
                return RRTResult(
                    True, nodes, diagnostics, terminal_diagnostics,
                    terminal.final_index, 0, "goal_reached", continuation_work,
                    local_steering_evaluations=local_steering_evaluations,
                )
            if budget_exhausted():
                return RRTResult(
                    False, nodes, diagnostics, terminal_diagnostics,
                    None, 0, "mechanics_work_budget", continuation_work,
                    local_steering_evaluations=local_steering_evaluations,
                )

        for iteration in range(1, self.config.maximum_iterations + 1):
            if budget_exhausted():
                return RRTResult(
                    False, nodes, diagnostics, terminal_diagnostics,
                    None, iteration - 1, "mechanics_work_budget",
                    continuation_work,
                    local_steering_evaluations=local_steering_evaluations,
                )
            parent_index = len(nodes) - 1
            parent = nodes[parent_index]
            desired = self._desired_step(parent.tip_position, goal_vector)
            if desired is None:
                break
            steering = self.session.evaluate_local_steering(parent.state)
            local_steering_evaluations += 1
            parent.tip_actuation_jacobian = tuple(
                tuple(row) for row in steering.tip_actuation_jacobian
            )
            proposal = self._steer(
                parent.state, parent.tip_actuation_jacobian, desired
            )
            if proposal is None:
                break
            target, command = proposal
            edge = self.session.attempt_continuation(
                parent.state, target, self.continuation_options
            )
            continuation_work += edge.attempted_steps
            accepted = edge.success
            if accepted:
                tip = _vector3(edge.state.tip_position, "goal-seeking tip")
                accepted = (
                    _distance(tip, goal_vector)
                    < _distance(parent.tip_position, goal_vector)
                    and _distance(tip, parent.tip_position)
                    >= self.config.minimum_tip_progress
                )
            diagnostics.append(EdgeDiagnostic(
                iteration=iteration,
                parent_index=parent_index,
                accepted=accepted,
                continuation_failed=not edge.success,
                sampled_tip=goal_vector,
                attempted_steps=edge.attempted_steps,
                rejected_steps=edge.rejected_steps,
                contacts_added=edge.contacts_added,
                contacts_released=edge.contacts_released,
                reached_path_fraction=edge.reached_path_fraction,
                failure_reason=edge.failure_reason,
            ))
            if not accepted:
                break
            nodes.append(RRTNode(
                edge.state,
                tip,
                parent_index,
                command,
                contact_history=_extended_contact_history(parent, edge.state),
            ))
            if _distance(tip, goal_vector) <= self.config.goal_neighborhood_radius:
                terminal = self._connect_terminal(
                    nodes, len(nodes) - 1, goal_vector
                )
                terminal_diagnostics.append(terminal.diagnostic)
                continuation_work += sum(
                    step.attempted_steps for step in terminal.diagnostic.steps
                )
                local_steering_evaluations += len({
                    step.connector_iteration for step in terminal.diagnostic.steps
                })
                if terminal.success:
                    return RRTResult(
                        True, nodes, diagnostics, terminal_diagnostics,
                        terminal.final_index, iteration, "goal_reached",
                        continuation_work,
                        local_steering_evaluations=local_steering_evaluations,
                    )
                break

        reason = "mechanics_work_budget" if budget_exhausted() else "local_stall"
        return RRTResult(
            False,
            nodes,
            diagnostics,
            terminal_diagnostics,
            None,
            len(diagnostics),
            reason,
            continuation_work,
            local_steering_evaluations=local_steering_evaluations,
        )
