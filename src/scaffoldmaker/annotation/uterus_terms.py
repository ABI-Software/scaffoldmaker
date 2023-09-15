"""
Common resource for uterus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
uterus_terms = [
    ("uterus", "UBERON:0000995"),
    ("serosa of uterus", "UBERON:0001297"),
    ("body of uterus", "UBERON:0009853"),
    ("cervix of uterus", "ILX:0745917 "),
    ("uterine cervix", "UBERON:0000002"),
    ("uterine wall", "UBERON:0000459"),
    ("uterine horn", "UBERON:000224"),
    ("uterine lumen", "UBERON:0013769"),
    ("left uterine horn", "UBERON:0009020"),
    ("right uterine horn", "UBERON:0009022"),
    ("serosa of uerine cervix", "None"),
    ("lumen of uerine cervix", "None"),
    ("serosa of right horn", "None"),
    ("lumen of right horn", "None"),
    ("serosa of left horn", "None"),
    ("lumen of left horn", "None"),
    ("serosa of uterus", "None"),
    ("lumen of uterus", "None"),
    ("fundus", "None"),
    ("ventral cervix junction with vagina", "None"),
    ("dorsal cervix junction with vagina", "None"),
    ("dorsal top right horn", "None"),
    ("ventral top right horn", "None"),
    ("dorsal top left horn", "None"),
    ("ventral top left horn", "None"),
    ("vagina", "None")
]

def get_uterus_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in uterus_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Uterus annotation term '" + name + "' not found.")