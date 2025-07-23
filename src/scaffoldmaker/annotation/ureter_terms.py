"""
Common resource for ureter annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
ureter_terms = [
    ("ureter", "UBERON:0000056", "ILX:0728080")
    ]

def get_ureter_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in ureter_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Body annotation term '" + name + "' not found.")