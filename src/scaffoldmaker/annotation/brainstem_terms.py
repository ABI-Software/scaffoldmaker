"""
Common resource for testing annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
brainstem_terms = [ # Landmarks and groups
                    ("brainstem", "UBERON:0002298", "ILX:0101444"),
                    ("central canal of spinal cord", "UBERON:0002291", "ILX:0724457"),
                    ("cerebral aqueduct", "UBERON:0002289", "ILX:0101977"),
                    ("diencephalon", "UBERON:0001894", "ILX:0103217"),
                    ("foramen caecum of medulla oblongata", "ILX:0746371"),
                    ("medulla oblongata", "UBERON:0001896", "ILX:0106736"),
                    ("midbrain", "UBERON:0001891", "ILX:0106935"),
                    ("middle cerebellar peduncle", "UBERON:0002152", "ILX:0106956"),
                    ("obex", "ILX:0107862"),
                    ("pons", "UBERON:0000988", "ILX:0109019"),

                    # Geometric markers
                    ("brainstem dorsal midline caudal point", "ILX:0778144"),
                    ("brainstem ventral midline caudal point", "ILX:0778145"),
                    ("brainstem dorsal midline cranial point", "ILX:0778146 "),
                    ("brainstem ventral midline cranial point", "ILX:0778147"),
                    ("brainstem dorsal midline pons-medulla junction", "ILX:0778148"),
                    ("brainstem ventral midline pons-medulla junction", "ILX:0778149"),
                    ("brainstem dorsal midline midbrain-pons junction", "ILX:0778150"),
                    ("brainstem ventral midline midbrain-pons junction", "ILX:0778151"),

                    # Surface
                    ("brainstem exterior", "ILX:0778157"),
                    ("midbrain exterior", "ILX:0778158"),
                    ("medulla oblongata exterior", "ILX:0778159"),
                    ("pons exterior", "ILX:0778160"),
                    ("brainstem-spinal cord interface", "ILX:0778162"),
                    ("thalamus-brainstem interface", "ILX:0778163")

]

def get_brainstem_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in brainstem_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Brainstem annotation term '" + name + "' not found.")
