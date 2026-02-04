from ase.optimize import FIRE2
from ase.constraints import FixAtoms
import os
import numpy as np


path_to_default_model = r"C:\Users\c1528354\GitHub\ucl-cci\cci-ucl-workshop\codebase\run\mace-omat-0-medium.model"

def __init_calc(model=path_to_default_model):
    # Lazy import to avoid slow startup times
    from mace.calculators import mace_mp
    PARAMS = {'model': model,
            'dispersion': False,
            'default_dtype': 'float64',
            'device': 'cpu'}

    calc =  mace_mp(**PARAMS)
    calc.directory = 'scratch'
    os.makedirs(calc.directory, exist_ok=True)

    return calc

# TODO: add user input for metal oxide / support slab species


def indentify_metal_oxide_atoms(atoms, threshold=1.0, fix=["Ce", "Ti", "O"]):
    '''Identify atoms below a certain z-coordinate threshold, only for specified species'''
    positions = atoms.get_positions()
    oxide_z_coordinate = float(np.max([[position[2] for _, position in enumerate(positions) if atoms[_].symbol in fix]]))

    threshold = oxide_z_coordinate - threshold
    symbols = np.array(atoms.get_chemical_symbols())
    fix_mask = np.isin(symbols, fix)
    bottom_mask = positions[:, 2] <= threshold
    indices = np.nonzero(~(fix_mask & bottom_mask))[0].tolist()

    return indices


def shave_slab(atoms, threshold=1.0, fix=["Ce", "Ti", "O"]):
    ''' Remove atoms below a certain z-coordinate threshold, only for specified species to fix'''
    keep_indices = indentify_metal_oxide_atoms(atoms, threshold, fix)
    if not keep_indices:
        return [], atoms
    temp_atoms = atoms.copy()[keep_indices]
    return keep_indices, temp_atoms

def evaluate_structure(atoms, calc, threshold=None, fix=["Ce", "Ti", "O"]):
    ''' Optimise slab, fixing a specified portion of the slab if requested '''

    # Reevaluate
    # Fix specified amount or just bottom layer
    if threshold:
        constraint = FixAtoms(indices=indentify_metal_oxide_atoms(atoms, threshold, fix))
    else:
        # Fix bottom layer
        z_cut = min(a.position[2] for a in atoms) + 2.0
        constraint = FixAtoms(indices=[a.index for a in atoms if a.position[2] < z_cut])

    atoms.set_constraint(constraint)
    atoms.calc = calc
    opt = FIRE2(atoms)
    opt.run(fmax=0.01)

    return atoms

def get_mace_energy(atoms):
    atoms.calc = __init_calc()
    energy = atoms.get_potential_energy()
    return energy

