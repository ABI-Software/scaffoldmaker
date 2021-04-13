"""
Common resource for lungs annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
lung_terms = [
    ("left lung", "UBERON:0002168", "ILX:0733450"),
    ("left lung bronchiole", "UBERON:0003539", "ILX:0734592"),
    ("left lung lower lobe bronchiole", "UBERON:0012055", "ILX:0736185"),
    ("left lung upper lobe bronchiole", "UBERON:0012056", "ILX:0727728"),
    ("left main bronchus", "UBERON:0002178", "ILX:0726074"),
    ("lower lobe of left lung", "UBERON:0008953", "ILX:0735534"),
    ("lower lobe of right lung", "UBERON:0002171", "ILX:0725712"),
    ("lung", "UBERON:0002048", "ILX:0726937"),
    ("middle lobe of right lung", "UBERON:0002174", "ILX:0733737"),
    ("respiratory airway", "UBERON:0001005", "ILX:0728555"),
    ("right lung", "UBERON:0002167", "ILX:0729582"),
    ("right lung accessory lobe", "UBERON:0004890", "ILX:0728696"),
    ("right lung accessory lobe bronchiole", "UBERON:0005682", "ILX:0736892"),
    ("right lung bronchioles", "UBERON:0003538", "ILX:0736755"),
    ("right lung lower lobe bronchiole", "UBERON:0012059", "ILX:0733237"),
    ("right lung middle lobe bronchiole", "UBERON:0012068", "ILX:0728473"),
    ("right lung upper lobe bronchiole", "UBERON:0005681", "ILX:0732909"),
    ("right main bronchus", "UBERON:0002177", "ILX:0728716"),
    ("trachea", "UBERON:0003126", "ILX:0729183"),
    ("upper lobe of left lung", "UBERON:0008952", "ILX:0735339"),
    ("upper lobe of right lung", "UBERON:0002170", "ILX:0728821"),
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
