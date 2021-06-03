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
                    ("brainstem dorsal midline caudal point", "None"),
                    ("brainstem ventral midline caudal point", "None"),
                    ("brainstem dorsal midline cranial point", "None"),
                    ("brainstem ventral midline cranial point", "None"),
                    ("brainstem dorsal midline pons-medulla junction", "None"),
                    ("brainstem ventral midline pons-medulla junction", "None"),
                    ("brainstem dorsal midline midbrain-pons junction", "None"),
                    ("brainstem ventral midline midbrain-pons junction", "None")
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
