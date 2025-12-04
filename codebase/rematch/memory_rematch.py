import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from ase.db import connect
from sklearn.preprocessing import normalize
from codebase.rematch.get_rematch import get_cached_rematch_kernel
from codebase.run.evaluate import shave_slab


def get_soap(atoms):
    from dscribe.descriptors import SOAP

    species = ["Ce", "O", "Au", "Cu"]
    r_cut = 6.0
    n_max = 8
    l_max = 6

    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=True,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
    )

    return soap.create(atoms)


def filter_similar_structures(db_path, 
                              similarity_threshold=0.9999, 
                              kernel_gamma=1.0, 
                              kernel_alpha=1.0, 
                              kernel_threshold=1e-6):
    # First pass: load and normalize SOAP descriptors (O(n * descriptor_size) memory).
    soap_list = []
    keep_mask = None
    for row in connect(db_path).select():
        atoms = row.toatoms()
        if keep_mask is None:
            # shaved atoms for speed - increase for better accuracy
            keep_mask = shave_slab(atoms, threshold=3.0, fix=["Ce", "O"])[0]
        soap_list.append(normalize(get_soap(atoms[keep_mask])))

    rematch_kernel = get_cached_rematch_kernel(
        gamma=kernel_gamma, alpha=kernel_alpha, threshold=kernel_threshold
    )

    active_indices = list(range(len(soap_list)))
    keep_indices = []
    i = 0

    # Iteratively prune: compare item i to all remaining items j>i and remove near-duplicates.
    while i < len(active_indices):
        current_index = active_indices[i]
        current_soap = soap_list[current_index]
        max_score = None

        j = i + 1
        while j < len(active_indices):
            other_index = active_indices[j]
            score = rematch_kernel.create([current_soap, soap_list[other_index]])[0, 1]
            if max_score is None or score > max_score:
                max_score = score

            if score >= similarity_threshold:
                print(
                    f"Pruning structure {other_index} (SOAP {score:.6f}) "
                    f"as duplicate of {current_index} (>= {similarity_threshold})."
                )
                active_indices.pop(j)
            else:
                j += 1

        keep_indices.append(current_index)
        similarity_msg = (
            f"{max_score:.6f}" if max_score is not None else "N/A (only structure remaining)"
        )
        print(
            f"Structure {current_index} max SOAP REMatch similarity: {similarity_msg}"
        )
        i += 1

    return keep_indices

def get_unique_db(db_path, db_out_path, similarity_threshold=0.9999):
    '''Take a database of Atoms objects and save unique structures
    in a new database'''

    keep_indices = filter_similar_structures(db_path, similarity_threshold)
    print(f"Keeping {len(keep_indices)} unique structures.")

    print(keep_indices)
    '''MAKE SURE TO SAVE THE KEEP INDICES SOMEWHERE SAFE'''
    with connect(db_path) as db:
        with connect(db_out_path) as db_out:
            for i in keep_indices:
                row = db.get(i + 1)  # ASE DB indices are 1-based
                db_out.write(row.toatoms())

if __name__ == "__main__":
    # Convenience entrypoint for ad-hoc runs; avoid running on import.
    get_unique_db(
        db_path=r"codebase\data\prescreened_structures.db",
        db_out_path="test_pruning.db",
    )
