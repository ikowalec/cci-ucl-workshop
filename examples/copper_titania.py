import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules
from codebase.rematch.prescreen import match_cluster_size, is_linked_periodically
from codebase.run.evaluate import evaluate_structure, __init_calc
from ase.db import connect

#initial_db = "adjacent.db"
#final_db = "adjacent_screened.db"

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


calc = __init_calc()
with connect("unique_combined.db") as db:
    with connect("unique_combined_result.db") as db2:
        for row in db.select():
            atoms = evaluate_structure(row.toatoms(), calc=calc, index=row.id)
            db2.write(atoms)
"""
'''Following optimisation:'''
with connect("scratch/CuTiO2_unique.db") as db:
    rows = db.select()
    from ase.visualize import view
    rows = sorted(rows, key=lambda k: k.energy)
    atoms_list = [row.toatoms() for row in rows]

atoms_intact = []
for atoms in atoms_list:
    if match_cluster_size(size=atoms.symbols.count("Cu"), 
                          slab=atoms, 
                          species=["Cu"]):
        
        atoms_intact.append(atoms)

separate_clusters = []
for atoms in atoms_intact:
    metal_indices = [atom.index for atom in atoms if atom.symbol in ["Cu"]]
    if not is_linked_periodically(13, atoms, species=["Cu"]):
        separate_clusters.append(atoms)

with connect("scratch/Cu13TiO2_optimised_screened.db") as db:
    for atoms in separate_clusters:
        db.write(atoms)

