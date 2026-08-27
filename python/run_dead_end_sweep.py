import argparse
from pathlib import Path

from simder_planning import PLANNING_METHODS, run_dead_end_sweep


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run one-factor robustness sweeps on the obstacle dead end"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(10)))
    parser.add_argument(
        "--profile", choices=("full", "regression"), default="full"
    )
    parser.add_argument(
        "--methods", nargs="+", choices=PLANNING_METHODS,
        default=list(PLANNING_METHODS),
    )
    arguments = parser.parse_args()

    sweep = run_dead_end_sweep(
        seeds=arguments.seeds,
        profile=arguments.profile,
        methods=arguments.methods,
    )
    sweep.save(arguments.output)
    for case in sweep.cases:
        print(f"{case.name} ({case.parameter}={case.value}):")
        for summary in case.summaries:
            print(
                f"  {summary.method}: {summary.success_count}/"
                f"{summary.run_count} success "
                f"(95% CI {summary.success_rate_ci95_low:.3f}-"
                f"{summary.success_rate_ci95_high:.3f}), "
                f"median_steps={summary.median_continuation_attempted_steps:.1f}"
            )
    print(f"Robustness sweep results={arguments.output}")


if __name__ == "__main__":
    main()
