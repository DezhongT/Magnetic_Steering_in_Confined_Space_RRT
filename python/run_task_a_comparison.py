import argparse
from pathlib import Path

from simder_planning import contact_free_task_a, run_comparison


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare simDER Task-A planners")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--full-tree", action="store_true")
    arguments = parser.parse_args()

    comparison = run_comparison(
        contact_free_task_a(arguments.seeds),
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
    print(f"Comparison results={arguments.output}")


if __name__ == "__main__":
    main()
