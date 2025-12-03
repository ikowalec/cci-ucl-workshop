import os
from ase import Atoms
from ase.db import connect

def parse_gin(gin_path):
    symbols = []
    positions = []
    cell = None
    reading_atoms = False

    with open(gin_path, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if not line or line.lower().startswith('name'):
                i += 1
                continue

            if line.lower().startswith('scell'):
                i += 1
                while i < len(lines) and not lines[i].strip():
                    i += 1
                if i >= len(lines):
                    raise ValueError("Expected cell parameters after 'scell'")
                parts = lines[i].strip().split()
                if len(parts) < 3:
                    raise ValueError("Not enough values in scell line")
                a, b, c = map(float, parts[:3])
                cell = [[a, 0, 0], [0, b, 0], [0, 0, c]]
                i += 1
                continue

            if line.lower().startswith('cartesian'):
                reading_atoms = True
                i += 1
                continue

            if reading_atoms:
                parts = line.split()
                if len(parts) < 5:
                    i += 1
                    continue
                if parts[1].lower() != 'core':
                    i += 1
                    continue

                symbol = parts[0]
                try:
                    x, y, z = map(float, parts[2:5])
                except ValueError:
                    raise ValueError(f"Cannot parse coordinates in line: {line}")

                symbols.append(symbol)
                positions.append([x, y, z])

            i += 1

    if cell is None:
        raise ValueError("No 'scell' cell parameters found in the .gin file.")

    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return atoms

def gin_files_to_db(folder_path, output_name="structures.db"):
    
    input_folder = os.path.abspath(folder_path)

    gin_files = [f for f in os.listdir(input_folder) if f.endswith(".gin")]

    if not gin_files:
        print(f"No .gin files found in {input_folder}")
        return

    with connect(os.path.join(input_folder, output_name)) as db:
        for gin_file in gin_files:
            input_path = os.path.join(input_folder, gin_file)

            try:
                atoms = parse_gin(input_path)
                db.write(atoms)

            except Exception as e:
                print(f"Failed to process {gin_file}: {e}")
    

import os
import tarfile
from io import BytesIO
from urllib.request import urlopen

from ase import Atoms
from ase.db import connect


def _parse_gin_from_lines(lines):
    symbols = []
    positions = []
    cell = None
    reading_atoms = False

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if not line or line.lower().startswith('name'):
            i += 1
            continue

        if line.lower().startswith('scell'):
            i += 1
            # skip blank lines
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines):
                raise ValueError("Expected cell parameters after 'scell'")
            parts = lines[i].strip().split()
            if len(parts) < 3:
                raise ValueError("Not enough values in scell line")
            a, b, c = map(float, parts[:3])
            cell = [[a, 0, 0], [0, b, 0], [0, 0, c]]
            i += 1
            continue

        if line.lower().startswith('cartesian'):
            reading_atoms = True
            i += 1
            continue

        if reading_atoms:
            parts = line.split()
            if len(parts) < 5:
                i += 1
                continue
            if parts[1].lower() != 'core':
                i += 1
                continue

            symbol = parts[0]
            try:
                x, y, z = map(float, parts[2:5])
            except ValueError:
                raise ValueError(f"Cannot parse coordinates in line: {line}")

            symbols.append(symbol)
            positions.append([x, y, z])

        i += 1

    if cell is None:
        raise ValueError("No 'scell' cell parameters found in the .gin file.")

    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return atoms


def parse_gin_path(gin_path):
    """Backward-compatible: parse from a file on disk."""
    with open(gin_path, 'r') as f:
        lines = f.readlines()
    return _parse_gin_from_lines(lines)


def parse_gin_bytes(gin_bytes, encoding="utf-8"):
    """New: parse from bytes (e.g. from tarfile.extractfile)."""
    text = gin_bytes.decode(encoding)
    lines = text.splitlines(keepends=False)
    return _parse_gin_from_lines(lines)


def gin_tar_url_to_db(url, output_db="structures.db", subdir="run/"):
    """
    Download a .tar.gz from `url`, read all *.gin files (optionally only
    those under `subdir` inside the tar), parse to Atoms, and write to ASE DB.
    No .gin files are written to disk.
    """
    # Download tar.gz into memory
    with urlopen(url) as r:
        tar_bytes = r.read()

    with tarfile.open(fileobj=BytesIO(tar_bytes), mode="r:gz") as tf, \
         connect(output_db) as db:

        for member in tf.getmembers():
            if not member.isfile():
                continue
            if not member.name.endswith(".gin"):
                continue
            # if you only want run/x0.gin, run/x1.gin, ...:
            if subdir and not member.name.startswith(subdir):
                continue

            fobj = tf.extractfile(member)
            if fobj is None:
                continue

            try:
                atoms = parse_gin_bytes(fobj.read())
            except Exception as e:
                print(f"Failed to process {member.name}: {e}")
                continue

            # store path inside tar as metadata if useful
            db.write(atoms, data={"gin_path_in_tar": member.name})


# If you have a local tar.gz instead of a URL:
def gin_tar_file_to_db(tar_path, output_db="structures.db", subdir="run/"):
    with tarfile.open(tar_path, mode="r:gz") as tf, \
         connect(output_db) as db:

        for member in tf.getmembers():
            if not member.isfile():
                continue
            if not member.name.endswith(".gin"):
                continue
            if subdir and not member.name.startswith(subdir):
                continue

            fobj = tf.extractfile(member)
            if fobj is None:
                continue

            try:
                atoms = parse_gin_bytes(fobj.read())
            except Exception as e:
                print(f"Failed to process {member.name}: {e}")
                continue

            db.write(atoms, data={"gin_path_in_tar": member.name})
