"""
Common resource for kidney annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
kidney_terms = [
    ("anterior surface of kidney", "UBERON:0035368", "ILX:0724840"),
    ("cortex of kidney", "UBERON:0001225", "ILX:0726853"),
    ("dorsal surface of kidney", ""),
    ("hilum of kidney", "UBERON:0008716", "ILX:0731719"),
    ("kidney", "UBERON:0002113", "ILX:0735723"),
    ("kidney capsule", "UBERON:0002015", "ILX:0733912"),
    ("lateral edge of kidney", ""),
    ("lateral surface of kidney", ""),
    ("major calyx", "UBERON:0001226", "ILX:0730785"),
    ("medial edge of kidney", ""),
    ("medial surface of kidney", ""),
    ("minor calyx", "UBERON:0001227", "ILX:0730473"),
    ("renal medulla", "UBERON:0000362", "ILX:0729114"),
    ("renal pelvis", "UBERON:0001224", "ILX:0723968"),
    ("renal pyramid", "UBERON:0004200", "ILX:0727514"),
    ("posterior surface of kidney", "UBERON:0035471", "ILX:0724479"),
    ("ventral surface of kidney", "")
    ]

def get_kidney_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in kidney_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Kidney annotation term '" + name + "' not found.")
