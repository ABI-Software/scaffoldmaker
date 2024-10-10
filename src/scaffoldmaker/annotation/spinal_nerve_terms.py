"""
Common resource for spinal nerve annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
spinal_nerve_terms = [
    ("dorsal root ganglion", "UBERON:0000044", "FMA:5888", "ILX:0103471"),
    ("dorsal root of spinal cord", "UBERON:0002261", "FMA:5980", "ILX:0726708"),
    ("spinal nerve", "UBERON:0001780", "FMA:5858", "ILX:0729912"),
    ("ventral root of spinal cord", "UBERON:0002260", "FMA:5979", "ILX:0724498")
    ]

def get_spinal_nerve_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in spinal_nerve_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Spinal nerve annotation term '" + name + "' not found.")
