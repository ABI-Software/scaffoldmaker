"""
Common resource for colon annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
colon_terms = [
    ("ascending colon", "UBERON:0001156", "ILX:0734393", "FMA:14545"),
    ("descending colon", "UBERON:0001158", "ILX:0724444", "FMA:14547"),
    ("circular muscle layer of colon", "ILX:0772428"),
    ("colon", "UBERON:0001155", "ILX:0736005", "FMA:14543"),
    ("colonic mucosa", "UBERON:0000317", "ILX:0731046", "FMA:14984"),
    ("left colon", "UBERON:0008971", "ILX:0727523"),
    ("longitudinal muscle layer of colon", "ILX:0775554"),
    ("luminal surface of the colonic mucosa", "ILX:0793083"),
    ("mesenteric zone", ""),
    ("non-mesenteric zone", ""),
    ("right colon", "UBERON:0008972", "ILX:0733240"),
    ("serosa of colon", "UBERON:0003335", "ILX:0736932", "FMA:14990"),
    ("spiral colon", "UBERON:0010239", "ILX:0735018"),
    ("submucosa of colon", "UBERON:0003331", "ILX:0728041", "FMA:14985"),
    ("taenia coli", "UBERON:0012419", "ILX:0731555", "FMA:15041"),
    ("taenia libera", "ILX:0739285"),
    ("taenia mesocolica", "ILX:0739284"),
    ("taenia omentalis", "ILX:0739286"),
    ("transverse colon", "UBERON:0001157", "ILX:0728767", "FMA:14546")
    ]


def get_colon_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in colon_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Colon annotation term '" + name + "' not found.")
