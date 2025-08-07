"""
Common resource for uterus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
uterus_terms = [
    ("body of uterus", "UBERON:0009853", "ILX:0730129", "FMA:17739"),
    ("broad ligament of uterus", "UBERON:0012332", "ILX:0733266", "FMA:16516"),
    ("external cervical os", "UBERON:0013760", "ILX:0736534", "FMA:76836"),
    ("fundus of uterus", "ILX:0743898"),
    ("internal cervical os", "UBERON:0013759", "ILX:0729495", "FMA:17747"),
    ("junction of left round ligament with uterus", ""),
    ("junction of left uterosacral ligament with uterus", ""),
    ("junction of right round ligament with uterus", ""),
    ("junction of right uterosacral ligament with uterus", ""),
    ("left broad ligament of uterus", ""),
    ("left lumen of uterus, cervix and vagina", ""),
    ("left lumen of uterus and cervix", ""),
    ("left oviduct", "UBERON:0001303", "ILX:0734218", "FMA:18484"),  # left uterine tube
    ("left serosa of uterus and vagina", ""),
    ("left transverse cervical ligament", ""),
    ("left uterine horn", "UBERON:0009020"),
    ("lumen of body of uterus", ""),
    ("lumen of fundus of uterus", ""),
    ("lumen of left uterine horn", ""),
    ("lumen of left oviduct", ""),
    ("lumen of right uterine horn", ""),
    ("lumen of right oviduct", ""),
    ("lumen of uterus, cervix and vagina", ""),
    ("lumen of uterus and cervix", ""),
    ("myometrium", "UBERON:0001296", "ILX:0735601", "FMA:17743"),
    ("pubocervical ligament", "ILX:0743760"),
    ("right broad ligament of uterus", ""),
    ("right lumen of uterus, cervix and vagina", ""),
    ("right lumen of uterus and cervix", ""),
    ("right oviduct", "UBERON:0001302", "ILX:0724908", "FMA:18483"),  # right uterine tube
    ("right serosa of uterus and vagina", ""),
    ("right transverse cervical ligament", ""),
    ("right uterine horn", "UBERON:0009022"),
    ("serosa of body of uterus", ""),
    ("serosa of fundus of uterus", ""),
    ("serosa of left oviduct", ""),
    ("serosa of left uterine horn", ""),
    ("serosa of right uterine horn", ""),
    ("serosa of right oviduct", ""),
    ("serosa of uterus", "UBERON:0001297"),
    ("serosa of uterus and vagina", ""),
    ("serosa of vagina", ""),
    ("uterine cervix", "UBERON:0000002", "ILX:0724162", "FMA:17740"),
    ("uterine cavity", "UBERON:0013769", "ILX:0727680", "FMA:17745"),
    ("uterine horn", "UBERON:000224"),
    ("uterine wall", "UBERON:0000459", "ILX:0735839", "FMA:17560"),
    ("uterus", "UBERON:0000995", "ILX:0726002", "FMA:17558"),
    ("vagina", "UBERON:0000996", "ILX:0736016", "FMA:19949"),
    ("vaginal canal", "UBERON:0011894", "ILX:0735924", "FMA:19982"),
    ("vagina orifice", "UBERON:0012317", "ILX:0729556", "FMA:19984")]


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
