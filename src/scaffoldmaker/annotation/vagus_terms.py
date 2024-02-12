"""
Common resource for vagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
vagus_terms = [
    ("centroid of level of superior border of the jugular foramen on the vagal trunk", "None"),
    ("centroid of level of inferior border of the jugular foramen on the vagal trunk", "None"),
    ("centroid of level of C1 transverse process on the vagal trunk", "None"),
    ("centroid of level of angle of mandible on the vagal trunk", "None"),
    ("centroid of level of tubercles of the greater horn of hyoid bone on the vagal trunk", "None"),
    ("centroid of level of carotid bifurcation on the vagal trunk", "None"),
    ("centroid of level of laryngeal prominence on the vagal trunk", "None"),
    ("centroid of level of superior border of clavicle on the vagal trunk", "None"),
    ("centroid of level of jugular notch on the vagal trunk", "None"),
    ("centroid of level of sternal angle on the vagal trunk", "None"),
    ("centroid of 1 cm superior to esophageal plexus on the vagal trunk", "None"),
    ("centroid of level of esophageal hiatus on the vagal trunk", "None"),
    ("centroid of level of aortic hiatus on the vagal trunk", "None"),
    ]

def get_vagus_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in vagus_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Vagus annotation term '" + name + "' not found.")