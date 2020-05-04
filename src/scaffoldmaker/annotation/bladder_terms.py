"""
Common resource for bladder annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
bladder_terms = [
    ( "urinary bladder", "FMA:15900", "UBERON:0001255" ),
    ( "neck of urinary bladder", "FMA:15912", "UBERON:0001258"),
    ( "Dome of the Bladder", None),  # needs to be replaced with an actual term e.g. urinary bladder (which includes neck)
    ( "Serosa of body of urinary bladder", None),
    ( "Lumen of body of urinary bladder", None),
    ("Serosa of neck of urinary bladder", None),
    ("Lumen of neck of urinary bladder", None)
   ]

def get_bladder_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in bladder_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Bladder annotation term '" + name + "' not found.")
