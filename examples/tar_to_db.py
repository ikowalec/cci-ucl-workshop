import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from codebase.io.klmc_parser import gin_tar_file_to_db

tarball = r"C:\Users\c1528354\GitHub\ucl-cci\cci-ucl-workshop\scratch\parallel.tar.gz"
gin_tar_file_to_db(tarball, output_db="parallel.db", subdir="run/")
