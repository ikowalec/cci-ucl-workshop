import sys
import os
import numpy as np
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from codebase.run.evaluate import evaluate_structure, get_mace_energy
from codebase.alloys.decoration import decorate_clusters
from ase.io import read

# generate structures or read - if generated - prescreen using connectivity matrix
atoms =  read("codebase/data/prescreened_structures.db@0")

species = ["Au", "Cu"]

# cluster structure can be optimised
# atoms = evaluate_structure(atoms, index=0)

# cluster alloy can be decorated for a number of different substitutions
atoms_list = decorate_clusters(atoms, to_substitute=["Au"], substitute_with=species, 
                      mutations={0:1, 1:5, 2:5, 3:15, 4:5, 5:5, 12:1})

from ase.visualize import view
# reconstruct
slab = atoms[[atom.index for atom in atoms if atom.symbol in ['Ce', 'O']]]
slabs = [slab + cluster for cluster in atoms_list]

# identical structures will have identical energy - filter these out
# Ideally we can use rematch/soap similarity here - how to define thresholds?
"""
unique_slabs = []
energies = []
for _, slab_structure in enumerate(slabs):
    print(_)
    energy = get_mace_energy(slab_structure[[atom.index for atom in slab_structure if atom.symbol in species]])
    if any(np.isclose(energy, existing) for existing in energies):
        continue
    unique_slabs.append(slab_structure)
    energies.append(energy)

slabs = unique_slabs
"""
