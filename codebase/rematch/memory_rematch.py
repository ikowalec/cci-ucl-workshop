from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from ase.db import connect
from sklearn.preprocessing import normalize
from codebase.rematch.get_rematch import get_cached_rematch_kernel
from codebase.rematch.get_rematch import soap_rematch_similarity
from codebase.run.evaluate import shave_slab
from ase import Atoms
from joblib import Parallel, delayed
import numpy as np
import time

def get_soap(atoms):
    from dscribe.descriptors import SOAP

    assert isinstance(atoms, Atoms), f"Must pass an Atoms object to get_soap, not {type(atoms)}."
    
    species = list(set(atoms.symbols))
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
                              slab_shave_distance=3.0,
                              similarity_threshold=0.9999, 
                              kernel_gamma=1.0, 
                              kernel_alpha=1.0, 
                              kernel_threshold=1e-6,
                              n_jobs=1,
                              backend="threading"):
    # First pass: load and normalize SOAP descriptors (O(n * descriptor_size) memory).
    soap_list = []
    keep_mask = None
    atoms_list = []
    for row in connect(db_path).select():
        atoms = row.toatoms()
        if keep_mask is None:
            # shaved atoms for speed - increase for better accuracy
            # TODO: These species should not be hardcoded
            keep_mask = shave_slab(atoms, threshold=slab_shave_distance, fix=["Ce", "Ti", "O"])[0]
        atoms_list.append(atoms[keep_mask])
    
    def process_atoms_to_soap(atoms):
        '''Helper function to save memory in Parallel() output list'''
        return normalize(get_soap(atoms))

    soap_list = Parallel(n_jobs=n_jobs, backend=backend)(delayed(process_atoms_to_soap)(atoms) for atoms in atoms_list)
    soap_list = np.array(soap_list)

    # clear memory
    del atoms_list

    print("SOAP memory usage is", soap_list.nbytes/1e6, "Mb.")

    rematch_kernel = get_cached_rematch_kernel(
        gamma=kernel_gamma, alpha=kernel_alpha, threshold=kernel_threshold
    )

    active_indices = np.array(list(range(len(soap_list))))
    keep_indices = []
    
    while active_indices.size != 0:

        # compare current_soap to the active subset
        targets = [soap_list[int(idx)] for idx in active_indices]
        s = time.time() # keep track of the time
        current_index = active_indices[0]
        print(f"Now processing: {current_index}")
        current_soap = targets[0]
        
        # Parallelise similarity calculation over cpus
        scores = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(rematch_kernel.create)([current_soap, soap])
            for soap in targets
        )
        scores = [round(score[0,1], 10) for score in scores] 
        scores = np.array(scores)

        # Find the indices to keep
        mask = scores < similarity_threshold
        active_indices = active_indices[mask]
        keep_indices.append(int(current_index + 1))
        
        # Print the timings
        print(time.time()-s, "s per iteration.")
   
    return keep_indices


def get_unique_db(db_path, db_out_path, slab_shave_distance=3.0, 
                  rematch_alpha=1.0, similarity_threshold=0.9999, n_jobs=1):
    '''Take a database of Atoms objects and save unique structures
    in a new database'''

    keep_indices = filter_similar_structures(
        db_path, kernel_alpha=rematch_alpha, similarity_threshold=similarity_threshold, n_jobs=n_jobs
    )
    print(f"Keeping {len(keep_indices)} unique structures.")

    print(keep_indices)
    '''MAKE SURE TO SAVE THE KEEP INDICES SOMEWHERE SAFE'''
    with connect(db_path) as db:
        with connect(db_out_path) as db_out:
            for i in keep_indices:
                row = db.get(i)  # ASE DB indices are 1-based
                db_out.write(
                    row.toatoms(),
                    data=dict(row.data),
                    key_value_pairs=dict(row.key_value_pairs),
                )

