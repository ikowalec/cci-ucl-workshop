import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from ase.db import connect
from codebase.rematch.get_rematch import soap_rematch_kernel_matrix
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

    # Extract SOAP to memory
    soap_list = []

    # initialise the list 
    keep = []
    # Retrieve trimmed structures from database
    for row in connect(db_path).select():
        if not keep:
            # shaved atoms for speed - increase for better accuracy
            keep = shave_slab(row.toatoms(), threshold=3.0, fix=["Ce", "O"])[0]

        # INITIALISE AN ARRAY AND ACCESS ARRAY ELEMENT INSTEAD OF APPENDING
        soap_list.append(get_soap(row.toatoms()[keep]))

    keep_indices = []

    similarity_matrix = soap_rematch_kernel_matrix(
        soap_list, gamma=kernel_gamma, alpha=kernel_alpha, threshold=kernel_threshold
    )

    active_indices = list(range(len(soap_list)))
    i = 0

    while i < len(active_indices):
        current_index = active_indices[i]
        max_score = None
        duplicates_to_remove = []

        for j, other_index in enumerate(active_indices):
            if current_index == other_index:
                continue

            score = similarity_matrix[current_index, other_index]
            if max_score is None or score > max_score:
                max_score = score

            # Only prune entries that appear later in the active list so previously
            # processed structures remain untouched.
            if j > i and score >= similarity_threshold:
                duplicates_to_remove.append((other_index, score))

        # Remove duplicates now that the similarity scores are known.
        if duplicates_to_remove:
            for duplicate_index, duplicate_score in duplicates_to_remove:
                active_indices.remove(duplicate_index)
                print(
                    f"Pruning structure {duplicate_index} (SOAP {duplicate_score:.6f}) "
                    f"as duplicate of {current_index} (>= {similarity_threshold})."
                )

        keep_indices.append(current_index)
        similarity_msg = (
            f"{max_score:.6f}" if max_score is not None else "N/A (only structure remaining)"
        )
        print(
            f"Structure {current_index} max SOAP REMatch similarity: {similarity_msg}"
        )
        i += 1

    return keep_indices

"""
keep_indices = filter_similar_structures("codebase/data/prescreened_structures.db")
print(f"Keeping {len(keep_indices)} unique structures.")
"""

'''MAKE SURE TO SAVE THE KEEP INDICES SOMEWHERE SAFE
with connect("codebase/data/prescreened_structures.db") as db:
    with connect("codebase/data/unique_structures.db") as db_out:
        for i in keep_indices:
            row = db.get(i + 1)  # ASE DB indices are 1-based
            db_out.write(row.toatoms())
'''