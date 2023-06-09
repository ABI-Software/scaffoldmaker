"""
Common resource for cecum annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
cecum_terms = [
    ("caecum", "UBERON:0001153", "FMA:14541", "ILX:0732270")
    ]


def get_cecum_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in cecum_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Cecum annotation term '" + name + "' not found.")
