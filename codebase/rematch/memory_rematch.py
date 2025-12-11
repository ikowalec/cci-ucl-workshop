from ase.db import connect
from sklearn.preprocessing import normalize
from codebase.rematch.get_rematch import get_cached_rematch_kernel
from codebase.run.evaluate import shave_slab
from ase import Atoms
from joblib import Parallel, delayed
import numpy as np
import time
from codebase.rematch.prescreen import split_by_1st_layer
from tqdm import tqdm

def get_soap(atoms):
    from dscribe.descriptors import SOAP

    assert isinstance(atoms, Atoms), f"Must pass an Atoms object to get_soap, not {type(atoms)}."
    
    species = list(set(atoms.symbols))
    r_cut = 6.0
    n_max = 8
    l_max = 6

    if atoms.cell:
        periodic = True
    else:
        periodic = False 
    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=periodic,
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
                              kernel_threshold=1e-6):
    # First pass: load and normalize SOAP descriptors (O(n * descriptor_size) memory).
    # Shaved typically 1 MB per structure
    soap_list = []
    keep_mask = None
    atoms_list = []

    # With statement avoids repeated I/O calls
    with connect(db_path) as db:
        for row in db.select():
            atoms = row.toatoms()
            if keep_mask is None:
                # shaved atoms for speed - increase for better accuracy
                # TODO: These species should not be hardcoded
                keep_mask = shave_slab(atoms, threshold=slab_shave_distance, fix=["Ce", "Ti", "O"])[0]
            atoms_list.append(atoms[keep_mask])
    
    def process_atoms_to_soap(atoms):
        '''Helper function to save memory in Parallel() output list'''
        return normalize(get_soap(atoms))


    soap_list = np.array(list(map(process_atoms_to_soap, atoms_list)))

    #soap_list = Parallel(n_jobs=n_jobs, backend=backend)(delayed(process_atoms_to_soap)(atoms) for atoms in atoms_list)
    #soap_list = np.array(soap_list)

    # clear memory
    del atoms_list

    print("SOAP memory usage is", soap_list.nbytes/1e6, "Mb.")

    rematch_kernel = get_cached_rematch_kernel(
        gamma=kernel_gamma, alpha=kernel_alpha, threshold=kernel_threshold
    )

    # Split the datasets to more manageable subsets 
    data_subsets = split_by_1st_layer(connect(db_path), layer_tolerance=1.0)
    # initialise the container for row ids to retain
    keep_indices = []
    
    for _,subset in enumerate(tqdm(data_subsets)):

        active_indices = np.array(data_subsets[subset])
        print(f"Starting subset {_} of {len(data_subsets)}, with size {len(active_indices)}.")
        s = time.time() # keep track of the time

        while active_indices.size != 0:

            # compare current_soap to the active subset
            targets = [soap_list[int(idx - 1)] for idx in active_indices]

            current_index = active_indices[0]
            current_soap = targets[0]
            
            def get_soap_score(soap):
                return round(rematch_kernel.create([current_soap, soap])[0,1], 8)

            scores = np.array(list(map(get_soap_score, targets)))

            # Find the indices to keep
            active_indices = active_indices[scores < similarity_threshold]
            keep_indices.append(int(current_index))
            
        # Print the timings
        print(f"{time.time()-s} s per subset.")
   
    return keep_indices


def get_unique_db(db_path, db_out_path, slab_shave_distance=3.0, 
                  rematch_alpha=1.0, similarity_threshold=0.9999):
    '''Take a database of Atoms objects and save unique structures
    in a new database'''

    keep_indices = filter_similar_structures(
        db_path, slab_shave_distance=slab_shave_distance, 
        kernel_alpha=rematch_alpha, similarity_threshold=similarity_threshold
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

