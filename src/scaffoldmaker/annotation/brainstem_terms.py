"""
Common resource for testing annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
brainstem_terms = [
                    ("medulla oblongata", "FMA:62004", "UBERON:0001896"),
                    ("pons", "FMA:67943", "UBERON:0000988"),
                    ("midbrain", "FMA:61993", "UBERON:0001891"),
                    ("diencephalon", "FMA:62001", "UBERON:0001894"),
                    ("brainstem", "FMA:79876", "UBERON:0002298")
                ]

def get_brainstem_annotation_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in brainstem_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Brainstem annotation term '" + name + "' not found.")
