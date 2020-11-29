"""
Common resource for bladder annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
bladder_terms = [
    ("urinary bladder", "UBERON:0001255", "FMA:15900"),
    ("neck of urinary bladder", "UBERON:0001258", "FMA:15912"),
    ("Dome of the Bladder", "ILX:0738433"),
    ("serosa of urinary bladder", "UBERON:0001260"),
    ("bladder lumen", "UBERON:0009958"),
    ("serosa of body of urinary bladder", "ILX:0739276"),
    ("lumen of body of urinary bladder", "ILX:0739252"),
    ("serosa of neck of urinary bladder", "ILX:0739277"),
    ("lumen of neck of urinary bladder", "ILX:0739256"),
    ("urethra", "UBERON:0000057", "ILX:0733022"),
    ("lumen of urethra", "ILX:0736762"),
    ("serosa of urethra", "ILX:0739282")
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
