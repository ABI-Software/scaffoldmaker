"""
Common resource for uterus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
uterus_terms = [
    ("body of uterus", "UBERON:0009853", "FMA:17739", "ILX:0730129"),
    ("broad ligament of uterus", "UBERON:0012332", "FMA:16516", "ILX:0733266"),
    ("dorsal cervix junction with vagina", "None"),
    ("dorsal top left horn", "None"),
    ("dorsal top right horn", "None"),
    ("external cervical os", "UBERON:0013760", "FMA:76836", "ILX:0736534"),
    ("fallopian tube", "UBERON:0003889", "FMA:18245", "ILX:0724851"),
    ("fundus", "None"),
    ("internal cervical os", "UBERON:0013759", "FMA:17747", "ILX:0729495"),
    ("left broad ligament of uterus", "None"),
    ("left round ligament of uterus", "None"),
    ("left transverse cervical ligament", "None"),
    ("left uterine horn", "UBERON:0009020"),
    ("lumen of body of uterus", "None"),
    ("lumen of fallopian tube", "None"),
    ("lumen of left horn", "None"),
    ("lumen of right horn", "None"),
    ("lumen of uterine cervix", "None"),
    ("lumen of uterus", "None"),
    ("lumen of vagina", "None"),
    ("myometrium", "UBERON:0001296", "FMA:17743", "	ILX:0735601"),
    ("pubocervical ligament (TA98)", "ILX:0743760"),
    ("right broad ligament of uterus", "None"),
    ("right round ligament of uterus", "None"),
    ("right transverse cervical ligament", "None"),
    ("right uterine horn", "UBERON:0009022"),
    ("round ligament of uterus", "UBERON:0006589", "ILX:0725512"),
    ("serosa of body of uterus", "None"),
    ("serosa of fallopian tube", "None"),
    ("serosa of left horn", "None"),
    ("serosa of right horn", "None"),
    ("serosa of uterine cervix", "None"),
    ("serosa of uterus", "UBERON:0001297"),
    ("serosa of vagina", "None"),
    ("uterine cervix", "UBERON:0000002","FMA:17740", "ILX:0724162"),
    ("uterine horn", "UBERON:000224"),
    ("uterine lumen", "UBERON:0013769"),
    ("uterine wall", "UBERON:0000459", "FMA:17560", "ILX:0735839"),
    ("uterus", "UBERON:0000995", "FMA:17558", "ILX:0726002"),
    ("vagina", "UBERON:0000996", "FMA:19949", "ILX:0736016"),
    ("vagina orifice", "UBERON:0012317", "FMA:19984", "ILX:0729556"),
    ("ventral cervix junction with vagina", "None"),
    ("ventral top left horn", "None"),
    ("ventral top right horn", "None")]

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