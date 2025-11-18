import sys
import os
import numpy as np
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

from codebase.rematch.get_rematch import soap_rematch_similarity 

def get_soap(atoms):
    from dscribe.descriptors import SOAP

    species = ["Ce", "O", "Au", "Cu"]
    r_cut = 6.0
    n_max = 8
    l_max = 6

    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=True,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
    )

    return soap.create(atoms)

from ase.io import read 
atoms = read(r'codebase\data\prescreened_structures.db', index=0)
atoms2 = read(r'codebase\data\prescreened_structures.db', index=50)
soap = get_soap(atoms)
soap2 = get_soap(atoms2)
score = soap_rematch_similarity(soap, soap2, gamma=1.0, alpha=1.0)
print("SOAP REMatch similarity:", score)