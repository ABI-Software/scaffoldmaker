"""
Common resource for kidney annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
kidney_terms = [
    ("anterior surface of kidney", "UBERON:0035368", "ILX:0724840"),
    ("anterior surface of left kidney", ""),
    ("anterior surface of right kidney", ""),
    ("cortex of kidney", "UBERON:0001225", "ILX:0726853"),
    ("cortex of left kidney", "ILX:0791219"),
    ("cortex of right kidney", "ILX:0791182"),
    ("dorsal surface of kidney", ""),
    ("dorsal surface of left kidney", ""),
    ("dorsal surface of right kidney", ""),
    ("hilum of kidney", "UBERON:0008716", "ILX:0731719"),
    ("hilum of left kidney", ""),
    ("hilum of right kidney", ""),
    ("inferior pole of left kidney", "FMA:15609"),
    ("inferior pole of right kidney", "FMA:15608"),
    ("inner medulla of left kidney", "ILX:0784932"),
    ("inner medulla of right kidney", "ILX:0791193"),
    ("juxtamedullary cortex", "UBERON:0005271", "ILX:0730126"),
    ("juxtamedullary cortex surface of kidney", ""),
    ("juxtamedullary cortex surface of left kidney", ""),
    ("juxtamedullary cortex surface of right kidney", ""),
    ("kidney", "UBERON:0002113", "ILX:0735723"),
    ("kidney capsule", "UBERON:0002015", "ILX:0733912"),
    ("lateral edge of kidney", ""),
    ("lateral edge of left kidney", ""),
    ("lateral edge of right kidney", ""),
    ("lateral surface of kidney", ""),
    ("lateral surface of left kidney", ""),
    ("lateral surface of right kidney", ""),
    ("left kidney", "UBERON:0004538", "ILX:0725163"),
    ("left kidney capsule", ""),
    ("major calyx", "UBERON:0001226", "ILX:0730785"),
    ("medial edge of kidney", ""),
    ("medial edge of left kidney", ""),
    ("medial edge of right kidney", ""),
    ("medial surface of kidney", ""),
    ("medial surface of left kidney", ""),
    ("medial surface of right kidney", ""),
    ("medulla of left kidney", ""),
    ("medulla of right kidney", ""),
    ("minor calyx", "UBERON:0001227", "ILX:0730473"),
    ("outer medulla of left kidney", ""),
    ("outer medulla of right kidney", ""),
    ("renal medulla", "UBERON:0000362", "ILX:0729114"),
    ("renal pelvis", "UBERON:0001224", "ILX:0723968"),
    ("renal pyramid", "UBERON:0004200", "ILX:0727514"),
    ("right kidney", "UBERON:0004539", "ILX:0735697"),
    ("right kidney capsule", ""),
    ("posterior surface of kidney", "UBERON:0035471", "ILX:0724479"),
    ("posterior surface of left kidney", ""),
    ("posterior surface of right kidney", ""),
    ("superior pole of left kidney", "FMA:15607"),
    ("superior pole of right kidney", "FMA:15606"),
    ("ventral surface of kidney", ""),
    ("ventral surface of left kidney", ""),
    ("ventral surface of right kidney", "")
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
