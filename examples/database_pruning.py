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
        default="2-layer-set-1.db", #"codebase/data/prescreened_structures.db",
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
        "--rematch-alpha",
        type=float,
        default=1.0,
        help="non-zero positive float - weighting between between the best match of local environments and the averaging strategy",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=4,
        help="Number of CPUs to be used for similarity check. 1 is serial",
    )
    parser.add_argument(
        "--slab-shave-distance",
        type=float,
        default=3.0,
        help="Provide the z-coordinate for slab atoms to be retained. e.g. 3.0 is about 2 layers, -1 is cluster-only.",
    )

    args = parser.parse_args()

    get_unique_db(
        db_path=args.db_path,
        db_out_path=args.db_out_path,
        similarity_threshold=args.similarity_threshold,
        rematch_alpha=args.rematch_alpha,
        slab_shave_distance=args.slab_shave_distance,
        n_jobs=args.n_jobs
    )


if __name__ == "__main__":
    import time

    start = time.time()
    main()
    print(f"Finished after {time.time()-start} s.")

