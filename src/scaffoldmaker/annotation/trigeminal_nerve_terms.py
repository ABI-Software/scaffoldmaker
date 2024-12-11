"""
Common resource for trigeminal nerve annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
trigeminal_nerve_terms = [
    ("mandibular nerve", "UBERON:0000375", "FMA:52996", "ILX:0732443"),
    ("maxillary nerve", "UBERON:0000377", "FMA:52724", "ILX:0728987"),
    ("ophthalmic nerve", "UBERON:0000348", "FMA:52621", "ILX:0734030"),
    ("trigeminal ganglion", "UBERON:0001675", "FMA:52618", "ILX:0724630"),
    ("trigeminal nerve", "UBERON:0001645", "FMA:50866", "ILX:0736097"),
    ("trigeminal nerve root", "UBERON:0004673", "FMA:52610", "ILX:0111966")
    ]

def get_trigeminal_nerve_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in trigeminal_nerve_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Trigeminal nerve annotation term '" + name + "' not found.")
