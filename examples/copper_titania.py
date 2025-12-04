import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from codebase.rematch.prescreen import match_cluster_size
from ase.db import connect

initial_db = "adjacent.db"
final_db = "adjacent_screened.db"
"""
import time
start = time.time()

with connect(initial_db) as idb:
    with connect(final_db) as fdb:
        for row in idb.select():
            atoms = row.toatoms()
            if match_cluster_size(size=13, slab=atoms, species='Cu'):
                fdb.write(atoms)

elapsed = time.time() - start
print(f"Screening of 25001 structures took {elapsed} seconds.")
"""

'''Note that KLMC structures are already unique'''

from codebase.run.evaluate import __init_calc, evaluate_structure

calc = __init_calc()
db_optimised = "adjacent_optimised.db"

with connect(db_optimised) as odb:
    with connect(final_db) as fdb:
        for row in fdb.select():
            atoms = row.toatoms()
            atoms = evaluate_structure(atoms, calc=calc, index=row.id, fix=["Ti", "O"])
            odb.write(atoms)
    