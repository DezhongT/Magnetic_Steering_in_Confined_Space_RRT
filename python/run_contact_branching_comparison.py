import argparse
from pathlib import Path

from simder_planning import (
    double_obstacle_contact_branching_task,
    run_comparison,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare planners on the two-obstacle contact branch"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--full-tree", action="store_true")
    arguments = parser.parse_args()

    comparison = run_comparison(
        double_obstacle_contact_branching_task(arguments.seeds),
        include_full_tree=arguments.full_tree,
    )
    comparison.save(arguments.output)
    for aggregate in comparison.aggregates:
        print(
            f"{aggregate.method}: success_rate={aggregate.success_rate:.3f}, "
            f"median_error={aggregate.median_final_tip_error:.6g}, "
            f"median_continuation_steps="
            f"{aggregate.median_continuation_attempted_steps:.1f}"
        )
    print(f"Contact-branching comparison results={arguments.output}")


if __name__ == "__main__":
    main()
