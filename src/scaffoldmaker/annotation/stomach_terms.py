"""
Common resource for stomach annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stomach_terms = [
    ( "body of stomach", "UBERON:0001161", " FMA:14560", "ILX:0724929"),
    ( "cardia of stomach", "UBERON:0001162", " FMA:14561", "ILX:0729096"),
    ( "duodenum", "UBERON:0002114", " FMA:7206", "ILX:0726125"),
    ( "esophagus", "UBERON:0001043", "FMA: 7131", "ILX:0735017"),
    ( "esophagogastric junction", "UBERON:0007650", "FMA: 9434", "ILX:0733910"),
    ( "forestomach-glandular stomach junction", "UBERON:0012270", "ILX:0729974"),
    ( "forestomach-glandular stomach junction on inner wall", None),
    ( "forestomach-glandular stomach junction on outer wall", None),
    ( "fundus of stomach", "UBERON:0001160", " FMA:14559", "ILX:0724443"),
    ( "pyloric antrum", "UBERON:0001165", " FMA:14579", "ILX:0728672"),
    ( "pylorus", "UBERON:0001166", " FMA:14581", "ILX:0734150"),
    ( "stomach", "UBERON:0000945", "FMA:7148", "ILX:0736697")
    ]

def get_stomach_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in stomach_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Stomach annotation term '" + name + "' not found.")
