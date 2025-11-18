import sys
import os
sys.path.append(rf"{os.getcwd()}")
from codebase.io.klmc_parser import gin_files_to_db
from codebase.rematch.prescreen import match_cluster_size

gin_files_to_db("codebase/data")

from ase.db import connect

prescreened_file = "codebase/data/prescreened_structures.db"
if os.path.exists(prescreened_file):
    with connect("codebase/data/structures.db") as db:
        with connect(prescreened_file) as db_out:
            for row in db.select():
                is_single_cluster = match_cluster_size(slab=row.toatoms(), size=12, species=["Au"])
                print(f"ID: {row.id}, Formula: {row.formula}, Single cluster {is_single_cluster}")
                
                if is_single_cluster:
                    db_out.write(row.toatoms())
           
from ase.visualize import view
from ase.io import read
# view(read(f"{prescreened_file}@:"))


# Calculate the energy of each structure and save in a database
# for image in images:
from mace.calculators import mace_mp
from ase.optimize import BFGS

PARAMS = {'model': "medium",
          'dispersion': True,
          'default_dtype': 'float64',
          'device': 'cpu'}

def evaluate_structure(atoms, index):
    atoms.calc = mace_mp(**PARAMS)
    opt = BFGS(atoms, trajectory=f"{index}.traj", logfile=f"{index}.log")
    opt.run(fmax=0.05)

with connect("codebase/data/prescreened_structures.db") as db:
    for row in db.select():
        evaluate_structure(row.toatoms(), row.id)