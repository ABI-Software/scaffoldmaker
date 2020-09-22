"""
Common resource for stellate annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stellate_terms = [
    ( "cervicothoracic ganglion", "UBERON:2441", "FMA:6469", "ILX:733799")
    ]

def get_stellate_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in stellate_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Stellate annotation term '" + name + "' not found.")
