"""
Common resource for body annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
body_terms = [
    ( "abdomen", "UBERON:0000916", "ILX:0725977" ),
    ( "thorax", "ILX:0742178" ),
    ( "neck", "UBERON:0000974", "ILX:0733967" ),
    ( "head", "UBERON:0000033", "ILX:0104909" ),
    ( "neck core", "" ),
    ( "head core", "" ),
    ( "skin epidermis", "UBERON:0001003", "ILX:0728574" ),
    ( "diaphragm", "UBERON:0001103", "ILX:0103194" ),
    ( "spinal cord", "UBERON:0002240", "ILX:0110909" ),
    ( "body", "UBERON:0000468", "ILX:0101370" ),
    ( "core", "" ),
    ( "non core", "" ),
    ( "core boundary", "" ),
    ]

def get_body_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in body_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Body annotation term '" + name + "' not found.")
