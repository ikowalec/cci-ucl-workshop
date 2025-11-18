
def strip_metal_oxide(atoms, keep_elements=['Au', 'Cu']):
    """
    Strip away all atoms that are not in keep_elements list.
    """
    from ase import Atoms

    indices_to_keep = [i for i, atom in enumerate(atoms) if atom.symbol in keep_elements]
    stripped_atoms = atoms[indices_to_keep]
    atoms.cell = None # Remove periodic boundary conditions
    return stripped_atoms


def decorate_clusters(atoms, 
                      to_substitute=["Au"], 
                      substitute_with=['Au', 'Cu'], 
                      mutations={1:5}):
    
    from acat.ga.group_operators import GroupSubstitute
    atoms = strip_metal_oxide(atoms, keep_elements=to_substitute)
    
    # Validations
    assert len(to_substitute) == 1, "Currently only single element substitution is supported."
    assert to_substitute[0] in atoms.get_chemical_symbols(), f"No atoms of type {to_substitute[0]} found in the structure."
    assert len(substitute_with) > 0, "substitute_with list cannot be empty."


    # Define groups: here each atom is its own group
    groups = [[atom.index] for atom in atoms if atom.symbol in to_substitute]

    atoms_list = []
    for num_mutations in mutations:
        for _ in range(mutations[num_mutations]):# Define which elements you allow in the alloy
            op = GroupSubstitute(
                groups=groups,
                elements=substitute_with,  # whatever alloy elements you want
                num_muts=num_mutations     # number of group substitutions per call
            )

            # Simple one-off mutation on this Atoms object
            atoms_list.append(op.substitute(atoms))   # does the actual substitution

    return atoms_list
