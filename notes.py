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
######
# Cluster decoration
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

    stripped_cluster = slab.copy()
    stripped_cluster = stripped_cluster[[atom.index for atom in slab if atom.symbol in species]]
    is_single_cluster = calculate_molecules(stripped_cluster, mult=1) == 1

    return is_single_cluster

from ase.build import fcc100
slab = fcc100('Pt', size=(4,4,4), vacuum=10.0)
print(match_cluster_size(16, slab, ['Pt']))

