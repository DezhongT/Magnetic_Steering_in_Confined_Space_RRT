import argparse
from pathlib import Path

from simder_planning import contact_free_task_a, run_task


def main() -> None:
    parser = argparse.ArgumentParser(description="Run simDER contact-free Task A")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--full-tree", action="store_true")
    arguments = parser.parse_args()

    task = contact_free_task_a(arguments.seeds)
    record = run_task(task, include_full_tree=arguments.full_tree)
    record.save(arguments.output)
    successes = sum(run.success for run in record.runs)
    print(
        f"Task A: {successes}/{len(record.runs)} successful; "
        f"results={arguments.output}"
    )


if __name__ == "__main__":
    main()
