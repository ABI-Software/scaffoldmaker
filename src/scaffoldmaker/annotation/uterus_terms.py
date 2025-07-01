"""
Common resource for uterus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
uterus_terms = [
    ("body junction", "None"),
    ("body of uterus", "UBERON:0009853", "ILX:0730129", "FMA:17739"),
    ("broad ligament of uterus", "UBERON:0012332", "ILX:0733266", "FMA:16516"),
    ("dorsal cervix junction with vagina", "None"),
    ("dorsal top left horn", "None"),
    ("dorsal top right horn", "None"),
    ("dorsal uterus", "None"),
    ("external cervical os", "UBERON:0013760", "ILX:0736534", "FMA:76836"),
    ("fundus of uterus", "None"),
    ("internal cervical os", "UBERON:0013759", "ILX:0729495", "FMA:17747"),
    ("junction of left round ligament with uterus", "None"),
    ("junction of right round ligament with uterus", "None"),
    ("junction of left uterosacral ligament with uterus", "None"),
    ("junction of right uterosacral ligament with uterus", "None"),
    ("left broad ligament of uterus", "None"),
    ("left transverse cervical ligament", "None"),
    ("left uterine horn", "UBERON:0009020"),
    ("left uterine tube", "UBERON:0001303", "ILX:0734218", "FMA:18484"),
    ("left uterus", "None"),
    ("lumen of body of uterus", "None"),
    ("lumen of fallopian tube", "None"),
    ("lumen of fundus of uterus", "None"),
    ("lumen of left horn", "None"),
    ("lumen of left uterine tube", "None"),
    ("lumen of right horn", "None"),
    ("lumen of right uterine tube", "None"),
    ("lumen of uterine cervix", "None"),
    ("lumen of uterus", "None"),
    ("lumen of vagina", "None"),
    ("myometrium", "UBERON:0001296", "ILX:0735601", "FMA:17743"),
    ("pubocervical ligament (TA98)", "ILX:0743760"),
    ("right broad ligament of uterus", "None"),
    ("right transverse cervical ligament", "None"),
    ("right uterine horn", "UBERON:0009022"),
    ("right uterine tube", "UBERON:0001302", "ILX:0724908", "FMA:18483"),
    ("right uterus", "None"),
    ("serosa of body of uterus", "None"),
    ("serosa of fundus of uterus", "None"),
    ("serosa of left uterine tube", "None"),
    ("serosa of left horn", "None"),
    ("serosa of right horn", "None"),
    ("serosa of right uterine tube", "None"),
    ("serosa of uterine cervix", "None"),
    ("serosa of uterus", "UBERON:0001297"),
    ("serosa of vagina", "None"),
    ("uterine cervix", "UBERON:0000002", "ILX:0724162", "FMA:17740"),
    ("uterine horn", "UBERON:000224"),
    ("uterine lumen", "UBERON:0013769"),
    ("uterine wall", "UBERON:0000459", "ILX:0735839", "FMA:17560"),
    ("uterus", "UBERON:0000995", "ILX:0726002", "FMA:17558"),
    ("vagina", "UBERON:0000996", "ILX:0736016", "FMA:19949"),
    ("vagina orifice", "UBERON:0012317", "ILX:0729556", "FMA:19984"),
    ("ventral cervix junction with vagina", "None"),
    ("ventral top left horn", "None"),
    ("ventral top right horn", "None"),
    ("ventral uterus", "None")]


def get_uterus_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return: ( preferred name, preferred id )
    """
    for term in uterus_terms:
        if name in term:
            return term[0], term[1]
    raise NameError("Uterus annotation term '" + name + "' not found.")
