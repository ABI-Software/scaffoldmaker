"""
Common resource for lungs annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
lung_terms = [
    ("apex of right lung accessory lobe", "ILX:0778119"),
    ("apex of left lung", "ILX:0778112"),
    ("apex of right lung", "ILX:0778113"),
    ("dorsal base of right lung accessory lobe", "ILX:0778125"),
    ("dorsal base of left lung", "ILX:0778126"),
    ("dorsal base of right lung", "ILX:0778127"),
    ("horizontal fissure of right lung", "ILX:0746327"),
    ("laterodorsal tip of middle lobe of right lung", "ILX:0778124"),
    ("left lung", "UBERON:0002168", "ILX:0733450"),
    ("lower lobe of left lung", "UBERON:0008953", "ILX:0735534"),
    ("lower lobe of right lung", "UBERON:0002171", "ILX:0725712"),
    ("lung", "UBERON:0002048", "ILX:0726937"),
    ("middle lobe of right lung", "UBERON:0002174", "ILX:0733737"),
    ("medial base of left lung", "ILX:0778120"),
    ("medial base of right lung", "ILX:0778121"),
    ("oblique fissure of left lung", "ILX:0778115"),
    ("oblique fissure of right lung", "ILX:0778114"),
    ("right lung", "UBERON:0002167", "ILX:0729582"),
    ("right lung accessory lobe", "UBERON:0004890", "ILX:0728696"),
    ("upper lobe of left lung", "UBERON:0008952", "ILX:0735339"),
    ("upper lobe of right lung", "UBERON:0002170", "ILX:0728821"),
    ("ventral base of right lung accessory lobe", "ILX:0778123"),
    ("ventral base of left lung", "ILX:0778118"),
    ("ventral base of right lung", "ILX:0778122")
]

def get_lung_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in lung_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Lung annotation term '" + name + "' not found.")
