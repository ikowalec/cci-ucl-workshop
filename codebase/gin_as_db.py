import sys
import os


'''THIS SNIPPET OF CODE POPULATES A DATABASE WITH STRUCTURES
FILTER TO ENSURE ONLY SINGLE CLUSTERS OF A GIVEN SIZE ARE INCLUDED'''

sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules
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
           
