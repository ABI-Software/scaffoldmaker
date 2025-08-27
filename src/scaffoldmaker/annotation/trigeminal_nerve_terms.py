"""
Common resource for trigeminal nerve annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
trigeminal_nerve_terms = [
    ("mandibular nerve", "UBERON:0000375", "ILX:0732443", "FMA:52996"),
    ("maxillary nerve", "UBERON:0000377", "ILX:0728987", "FMA:52724"),
    ("ophthalmic nerve", "UBERON:0000348", "ILX:0734030", "FMA:52621"),
    ("trigeminal ganglion", "UBERON:0001675", "ILX:0724630", "FMA:52618"),
    ("trigeminal nerve", "UBERON:0001645", "ILX:0736097", "FMA:50866"),
    ("trigeminal nerve root", "UBERON:0004673", "ILX:0111966", "FMA:52610")
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
