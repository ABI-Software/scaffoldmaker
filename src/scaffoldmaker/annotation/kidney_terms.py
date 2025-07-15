"""
Common resource for kidney annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
kidney_terms = [
    ("core", ""),
    ("kidney capsule", "UBERON:0002015", "ILX:0733912"),
    ("major calyx", "UBERON:0001226", "ILX:0730785"),
    ("minor calyx", "UBERON:0001227", "ILX:0730473"),
    ("renal pelvis", "UBERON:0001224", "ILX:0723968"),
    ("renal pyramid", "UBERON:0004200", "ILX:0727514"),
    ("shell", "")
    ]

def get_kidney_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in kidney_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Kidney annotation term '" + name + "' not found.")
