"""
Common resource for lungs annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
lung_terms = [
    ("apex of left lung", "None"),
    ("apex of right lung", "None"),
    ("horizontal fissure of right lung", "None"),
    ("lateral side of right lung", "None"),
    ("left lung", "UBERON:0002168", "ILX:0733450"),
    ("lower lobe of left lung", "UBERON:0008953", "ILX:0735534"),
    ("lower lobe of right lung", "UBERON:0002171", "ILX:0725712"),
    ("lung", "UBERON:0002048", "ILX:0726937"),
    ("middle lobe of right lung", "UBERON:0002174", "ILX:0733737"),
    ("medial base of left lung", "None"),
    ("oblique fissure of left lung", "None"),
    ("oblique fissure of right lung", "None"),
    ("right lung", "UBERON:0002167", "ILX:0729582"),
    ("right lung accessory lobe", "UBERON:0004890", "ILX:0728696"),
    ("upper lobe of left lung", "UBERON:0008952", "ILX:0735339"),
    ("upper lobe of right lung", "UBERON:0002170", "ILX:0728821"),
    ("ventral base of left lung", "None"),
    ("ventral base of right lung", "None")
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
