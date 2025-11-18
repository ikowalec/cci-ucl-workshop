from ase.db import connect
from mace.calculators import mace_mp
from ase.optimize import FIRE2
from ase.constraints import FixAtoms
import os
from ase.visualize import view
import numpy as np

PARAMS = {'model': "small",
          'dispersion': False,
          'default_dtype': 'float32',
          'device': 'cpu'}

calc =  mace_mp(**PARAMS)
calc.directory = 'scratch'
os.makedirs(calc.directory, exist_ok=True)


def evaluate_structure(atoms, index, fix=["Ce", "O"]):
    atoms.calc = calc
    positions = atoms.get_positions()
    oxide_z_coordinate = float(np.max([[position[2] for _, position in enumerate(positions) if atoms[_].symbol in fix]]))

    threshold = oxide_z_coordinate - 2.0
    symbols = np.array(atoms.get_chemical_symbols())
    fix_mask = np.isin(symbols, fix)
    bottom_mask = positions[:, 2] <= threshold
    keep_indices = np.nonzero(~(fix_mask & bottom_mask))[0].tolist()
    if not keep_indices:
        return
    temp_atoms = atoms.copy()[keep_indices]
    
    # Run preoptimisation on reduced slab
    c = FixAtoms(indices=[atom.index for atom in temp_atoms if atom.symbol in fix])
    temp_atoms.set_constraint(c)
    temp_atoms.calc = calc
    opt = FIRE2(temp_atoms, trajectory=f"{calc.directory}/pre_{index}.traj", logfile=f"{calc.directory}/{index}.log")
    opt.run(fmax=0.05)

    # Restore positions to original atoms object
    counter = 0
    for atom in atoms:
        if atom.index in keep_indices and not atom.symbol in fix:
            cluster = temp_atoms[[temp_atom.index for temp_atom in temp_atoms if not temp_atom.symbol in fix]]
            atom.position = cluster[counter].position
            counter += 1
    
    # Reevaluate
    c = FixAtoms(indices=[atom.index for atom in atoms if atom.symbol in fix])
    atoms.set_constraint(c)
    atoms.calc = calc
    opt = FIRE2(atoms, trajectory=f"{calc.directory}/{index}.traj", logfile=f"{calc.directory}/{index}.log")
    opt.run(fmax=0.01)

    return atoms

def get_mace_energy(atoms):
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    return energy

"""
with connect("codebase/data/prescreened_structures.db") as db:
    for row in db.select():
        
        view(evaluate_structure(row.toatoms(), row.id))
        if row.id == 1:
            break
"""
