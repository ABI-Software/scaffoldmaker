"""
Common resource for testing annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
bone_terms = [
                    #General bones
                    ("Cancellous bone", "ILX:0491986"),
                    ("Cortical bone", "ILX:0745808"),
                    # Geometric markers for bones
                    ("Olecranon", "UBERON:0006810","ILX:0735843"),
                    ("Coronoid process of ulna", "UBERON:0010994","ILX:0735423"),
                    ("Greater tubercle of radius","ILX:0745221"),
                    ("Medial epicondyle of radius","ILX:0748959"),
                    ("Lateral epicondyle of radius", "ILX:0748301"),
                    ("Intercondylar eminence of tibia", "ILX:0743525"),
                    ("Inferior articular surface of tibia", "ILX:0748099"),
                    ("Medial surface of tibia", "ILX:0747893"),
                    ("Lateral surface of tibia", "ILX:0742857"),

]


def get_bone_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in bone_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Bone annotation term '" + name + "' not found.")
