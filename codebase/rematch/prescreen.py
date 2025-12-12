"""
# DFT protocol -> get metal oxide surfaces
# Prepare the MASTER.GIN file
    - put a slab of metal on the metal oxide surface
    - save as a .gin
# KLMC cluster generation
# Pre-screening for similar and unstable structures -> metal cluster with selected sizes
# metal atoms only -> connectivity matrix -> make sure the no. connected atoms = n


# remove the "multiple cluster" system
# SOS/Rematch stuff (Yunyu)
# Run Mace calculations on the sctructures
# Identify the most stable structures and similar structures
# Take a sample -> run DFT on selected top X structures
# Compare ordering from MACE and DFT
# Once happy - choose clusters
# Retrain MACE on the DFT data
######
# Cluster decoration
# Cluster Expansion (for getting ground states of alloyed structures)>
"""

def calculate_molecules(atoms,  mult=1, print_output=False):
    '''
    Returns information on the molecules in the system. Needs refinement!
    
    
    Parameters:
    atoms: Atoms object
        Input structure from which to calculate molecular information
    Mult: a multiplier for the cutoffs
        Set to 1 as default, but can be adjusted depending on application
    print_output: Boolean 
        Whether to print any information whilst working out molecules

    Returns:
    n_molecules: int
        Number of molecules in the system
   
    '''

    from ase.neighborlist import natural_cutoffs, NeighborList
    from scipy import sparse
    cutOff = natural_cutoffs(atoms, mult=mult)
    neighborList = NeighborList(cutOff, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    n_molecules, component_list = sparse.csgraph.connected_components(matrix)

    return n_molecules

def match_cluster_size(size: int, slab, species: list):

    # TODO: use the size also
    stripped_cluster = slab.copy()
    stripped_cluster = stripped_cluster[[atom.index for atom in slab if atom.symbol in species]]
    is_single_cluster = calculate_molecules(stripped_cluster, mult=0.9) == 1

    return is_single_cluster

def is_linked_periodically(size: int, slab, species: list):
    if not (slab.pbc).all():
        return False
    
    stripped_cluster = slab.copy()
    stripped_cluster = stripped_cluster[[atom.index for atom in slab if atom.symbol in species]]
    stripped_cluster = stripped_cluster.repeat((2,2,2))
    is_linked_periodically =  calculate_molecules(stripped_cluster, mult=1) !=  8

    return is_linked_periodically

def split_by_1st_layer(db, species=["Au", "Cu"], layer_tolerance=0.5):
    import numpy as np
    atoms_list = [row.toatoms()[[atom.index for atom in row.toatoms() if atom.symbol in species]] for row in db.select()]
    no_atoms = [len([atom.index for atom in atoms if atom.z < np.amin(atoms.positions[:,2]) + layer_tolerance]) for atoms in atoms_list]
    
    split = {n:[] for n in set(no_atoms)}
    print(split)

    for index, count in enumerate(no_atoms):
        split[count].append(index+1)

    # therefore we get {layer_atom_count:database_indices}

    return split