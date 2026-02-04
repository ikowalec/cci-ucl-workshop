from ase.io import read
from ase.db import connect
import os

def vasp_outputs_to_db(folder, output_db):
    """
    Convert a list of VASP output files (e.g., OUTCAR or CONTCAR) to an ASE database.

    Parameters:
    vasp_output_files (list of str): List of paths to VASP output files.
    output_db (str): Path to the output ASE database file.
    """
    filenames = os.listdir(folder)

    with connect(output_db) as db:
        for vasp_file in filenames:
            try:
                atoms = read(os.path.join(folder, vasp_file))
                db.write(atoms, filename=vasp_file)
                print(f"Added {vasp_file} to {output_db}")
            except Exception as e:
                print(f"Error processing {vasp_file}: {e}")

if __name__ == "__main__":
    folder = "path/to/vasp/outputs"  # Replace with your folder path
    output_db = "vasp_structures.db"  # Replace with your desired output database name
    vasp_outputs_to_db(folder, output_db)