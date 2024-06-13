"""
Common resource for vagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
vagus_terms = [
    # marker names
    ("level of superior border of jugular foramen on the vagus nerve", "ILX:0794617"),
    ("right level of superior border of jugular foramen on the vagus nerve", "ILX:0794618"),
    ("left level of superior border of jugular foramen on the vagus nerve", "ILX:0794619"),

    ("level of inferior border of jugular foramen on the vagus nerve", "ILX:0794620"),
    ("right level of inferior border of jugular foramen on the vagus nerve", "ILX:0794621"),
    ("left level of inferior border of jugular foramen on the vagus nerve", "ILX:0794622"),

    ("level of inferior border of cranium on the vagus nerve", "ILX:0794623"),
    ("right level of inferior border of cranium on the vagus nerve", "ILX:0794624"),
    ("left level of inferior border of cranium on the vagus nerve", "ILX:0794625"),

    ("level of C1 transverse process on the vagus nerve", "ILX:0794626"),
    ("right level of C1 transverse process on the vagus nerve", "ILX:0794627"),
    ("left level of C1 transverse process on the vagus nerve", "ILX:0794628"),

    ("level of greater horn of hyoid on the vagus nerve", "ILX:0794629"),
    ("right level of greater horn of hyoid on the vagus nerve", "ILX:0794630"),
    ("left level of greater horn of hyoid on the vagus nerve", "ILX:0794631"),

    ("level of laryngeal prominence on the vagus nerve", "ILX:0794632"),
    ("right level of laryngeal prominence on the vagus nerve", "ILX:0794633"),
    ("left level of laryngeal prominence on the vagus nerve","ILX:0794634"),

    ("level of angle of the mandible on the vagus nerve", "ILX:0794635"),
    ("right level of angle of the mandible on the vagus nerve", "ILX:0794636"),
    ("left level of angle of the mandible on the vagus nerve", "ILX:0794637"),

    ("level of carotid bifurcation on the vagus nerve", "ILX:0794638"),
    ("right level of carotid bifurcation on the vagus nerve", "ILX:0794639"),
    ("left level of carotid bifurcation on the vagus nerve", "ILX:0794640"),

    ("level of superior border of the clavicle on the vagus nerve", "ILX:0794641"),
    ("right level of superior border of the clavicle on the vagus nerve", "ILX:0794642"),
    ("left level of superior border of the clavicle on the vagus nerve", "ILX:0794643"),

    ("level of jugular notch on the vagus nerve", "ILX:0794644"),
    ("right level of jugular notch on the vagus nerve", "ILX:0794645"),
    ("left level of jugular notch on the vagus nerve", "ILX:0794646"),

    ("level of sternal angle on the vagus nerve", "ILX:0794647"),
    ("right level of sternal angle on the vagus nerve", "ILX:0794648"),
    ("left level of sternal angle on the vagus nerve", "ILX:0794649"),

    ("1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794650"), # !!!
    ("right 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794651"), # !!!
    ("left 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794652"), # !!!

    ("level of esophageal hiatus on the vagus nerve", "ILX:0794653"),
    ("right level of esophageal hiatus on the vagus nerve", "ILX:0794654"),
    ("left level of esophageal hiatus on the vagus nerve", "ILX:0794655"),

    ("level of aortic hiatus on the vagus nerve", "ILX:0794656"),
    ("right level of aortic hiatus on the vagus nerve", "ILX:0794657"),
    ("left level of aortic hiatus on the vagus nerve", "ILX:0794658"),

    # branch names - in use at the moment
    ()
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