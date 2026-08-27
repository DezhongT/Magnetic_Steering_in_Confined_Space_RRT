import argparse
from pathlib import Path

from simder_planning import (
    PLANNING_METHODS,
    RobustnessSweepRecord,
    contact_branching_budget_cases,
    run_sweep,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sweep deterministic planning budgets on the contact branch"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(10)))
    parser.add_argument(
        "--profile", choices=("full", "regression"), default="full"
    )
    parser.add_argument(
        "--methods", nargs="+", choices=PLANNING_METHODS,
        default=["mechanics_informed_rrt"],
    )
    parser.add_argument(
        "--case", action="append", dest="cases",
        help="Run only a named case; repeat to create a deterministic shard",
    )
    shortlist_group = parser.add_mutually_exclusive_group()
    shortlist_group.add_argument(
        "--history-shortlist",
        action="store_true",
        help=(
            "Select among the nearest nodes from each accumulated contact "
            "history using feasible local mechanics"
        ),
    )
    shortlist_group.add_argument(
        "--history-shortlist-period",
        type=int,
        help="Invoke the history shortlist every N planner iterations",
    )
    shortlist_group.add_argument(
        "--history-shortlist-after",
        type=int,
        help="Use the history shortlist on every iteration starting at N",
    )
    shortlist_group.add_argument(
        "--history-shortlist-stagnation",
        type=int,
        help="Invoke the history shortlist after N iterations without progress",
    )
    parser.add_argument(
        "--history-shortlist-progress",
        type=float,
        default=1.0e-4,
        help="Goal-distance reduction that resets the stagnation trigger",
    )
    parser.add_argument(
        "--history-shortlist-burst",
        type=int,
        default=1,
        help="Number of consecutive shortlist iterations after a hybrid trigger",
    )
    arguments = parser.parse_args()

    cases = contact_branching_budget_cases(
        seeds=arguments.seeds,
        profile=arguments.profile,
        use_contact_history_shortlist=arguments.history_shortlist,
        contact_history_shortlist_start_iteration=(
            arguments.history_shortlist_after
        ),
        contact_history_shortlist_period=arguments.history_shortlist_period,
        contact_history_shortlist_stagnation_iterations=(
            arguments.history_shortlist_stagnation
        ),
        contact_history_shortlist_stagnation_progress=(
            arguments.history_shortlist_progress
        ),
        contact_history_shortlist_burst_iterations=(
            arguments.history_shortlist_burst
        ),
    )
    if arguments.cases:
        requested = set(arguments.cases)
        available = {case.name for case in cases}
        unknown = requested - available
        if unknown:
            parser.error(f"unknown case(s): {', '.join(sorted(unknown))}")
        cases = tuple(case for case in cases if case.name in requested)

    completed = []
    for case in cases:
        completed.extend(run_sweep((case,), arguments.methods).cases)
        # Checkpoint every completed case so an interrupted long sweep retains
        # a valid prefix. Separate --case invocations can be run as shards.
        RobustnessSweepRecord(cases=tuple(completed)).save(arguments.output)
    sweep = RobustnessSweepRecord(cases=tuple(completed))
    for case in sweep.cases:
        print(f"{case.name} ({case.parameter}={case.value:g}):")
        for summary in case.summaries:
            reasons = dict(summary.termination_reason_counts)
            histories = dict(summary.successful_contact_history_counts)
            print(
                f"  {summary.method}: {summary.success_count}/"
                f"{summary.run_count} success "
                f"(95% CI {summary.success_rate_ci95_low:.3f}-"
                f"{summary.success_rate_ci95_high:.3f}), "
                f"median_steps={summary.median_continuation_attempted_steps:.1f}, "
                f"termination={reasons}, successful_histories={histories}"
            )
    print(f"Contact-branching budget sweep results={arguments.output}")


if __name__ == "__main__":
    main()
