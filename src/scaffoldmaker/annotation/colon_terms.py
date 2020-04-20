"""
Common resource for colon annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
colon_terms = [
    ( "colon", "UBERON:0001155", "FMA:14543" ),
    ( "mesenteric zone", None ),
    ( "non-mesenteric zone", None ),
    ( "taenia coli", "UBERON:0012419", "FMA:15041" ),
    ( "tenia libera", None ),
    ( "tenia mesocolica", None ),
    ( "tenia omentalis", None )
    ]

def get_colon_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in colon_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Colon annotation term '" + name + "' not found.")
