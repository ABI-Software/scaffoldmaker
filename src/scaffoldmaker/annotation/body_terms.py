"""
Common resource for body annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
body_terms = [
    ("abdomen", "UBERON:0000916", "ILX:0725977"),
    ("abdominal cavity", "UBERON:0003684"),
    ("arm", "UBERON:0001460"),
    ("left arm", "FMA:24896"),
    ("right arm", "FMA:24895"),
    ("body", "UBERON:0000468", "ILX:0101370"),
    ("core", ""),
    ("core boundary", ""),
    ("head", "UBERON:0000033", "ILX:0104909"),
    ("head core", ""),
    ("diaphragm", "UBERON:0001103", "ILX:0103194"),
    ("hand", "FMA:9712"),
    ("leg", "UBERON:0000978"),
    ("left leg", "FMA:24981"),
    ("right leg", "FMA:24980"),
    ("foot", "FMA:9664"),
    ("neck", "UBERON:0000974", "ILX:0733967"),
    ("neck core", ""),
    ("non core", ""),
    ("shell", ""),
    ("skin epidermis", "UBERON:0001003", "ILX:0728574"),
    ("spinal cord", "UBERON:0002240", "ILX:0110909"),
    ("thoracic cavity", "UBERON:0002224"),
    ("thorax", "ILX:0742178"),
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
