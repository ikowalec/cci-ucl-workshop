import argparse
import os
import sys

# Make local project imports work when executed from the repo root.
sys.path.append(rf"{os.getcwd()}")

from codebase.rematch.memory_rematch import get_unique_db


def main():
    parser = argparse.ArgumentParser(
        description="Prune similar structures from an ASE database using SOAP REMatch."
    )
    parser.add_argument(
        "--db-path",
        default="codebase/data/prescreened_structures.db",
        help="Input ASE database containing candidate structures.",
    )
    parser.add_argument(
        "--db-out-path",
        default="unique.db",
        help="Output ASE database that will store unique structures.",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.9999,
        help="SOAP REMatch similarity threshold above which structures are treated as duplicates.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=4,
        help="Number of CPUs to be used for similarity check. 1 is serial",
    )

    args = parser.parse_args()

    get_unique_db(
        db_path=args.db_path,
        db_out_path=args.db_out_path,
        similarity_threshold=args.similarity_threshold,
        n_jobs=args.n_jobs
    )


if __name__ == "__main__":
    import time

    start = time.time()
    main()
    print(f"Finished after {time.time()-start} s.")

