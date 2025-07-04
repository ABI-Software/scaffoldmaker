"""
Common resource for spinal nerve annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
spinal_nerve_terms = [
    ("C1 spinal nerve",             "ILX:0794592"),
    ("C2 spinal nerve",             "ILX:0794595"),
    ("C3 spinal nerve",             "ILX:0794598"),
    ("C4 spinal nerve",             "ILX:0794601"),
    ("C5 spinal nerve",             "ILX:0794604"),
    ("C6 spinal nerve",             "ILX:0794607"),
    ("C7 spinal nerve",             "ILX:0794610"),
    ("C8 spinal nerve",             "ILX:0794613"),
    ("cervical spinal nerve",       "UBERON:0000962",   "FMA:5859", "ILX:0733509"),
    ("coccygeal nerve",             "UBERON:0009629",   "FMA:5863", "ILX:0723795"),
    ("dorsal root ganglion",        "UBERON:0000044",   "FMA:5888", "ILX:0103471"),
    ("dorsal root of spinal cord",  "UBERON:0002261",   "FMA:5980", "ILX:0726708"),
    ("eighth thoracic nerve",       "ILX:0792754",      "FMA:6302"),
    ("eleventh thoracic nerve",     "ILX:0791454",      "FMA:6312"),
    ("fifth lumbar nerve",          "ILX:0787261",      "FMA:6416"),
    ("fifth sacral nerve",          "ILX:0787432",      "FMA:6427"),
    ("fifth thoracic nerve",        "ILX:0792589",      "FMA:6293"),
    ("first lumbar nerve",          "ILX:0788972",      "FMA:6172"),
    ("first sacral nerve",          "ILX:0790504",      "FMA:6423"),
    ("first thoracic nerve",        "ILX:0738397",      "FMA:6037"),
    ("fourth lumbar nerve",         "ILX:0788896",      "FMA:6415"),
    ("fourth sacral nerve",         "ILX:0792327",      "FMA:6426"),
    ("fourth thoracic nerve",       "ILX:0784755",      "FMA:6290"),
    ("lumbar nerve",                "UBERON:0009624",   "FMA:5861", "ILX:0735627"),
    ("ninth thoracic nerve",        "ILX:0789545",      "FMA:6306"),
    ("sacral nerve",                "UBERON:0009625",   "FMA:5862", "ILX:0726145"),
    ("second lumbar nerve",         "ILX:0786362",      "FMA:6176"),
    ("second sacral nerve",         "ILX:0787528",      "FMA:6424"),
    ("second thoracic nerve",       "ILX:0739243",      "FMA:6169"),
    ("seventh thoracic nerve",      "ILX:0789008",      "FMA:6299"),
    ("sixth thoracic nerve",        "ILX:0791122",      "FMA:6296"),
    ("spinal nerve",                "UBERON:0001780",   "FMA:5858", "ILX:0729912"),
    ("tenth thoracic nerve",        "ILX:0792933",      "FMA:6309"),
    ("third lumbar nerve",          "ILX:0789374",      "FMA:6414"),
    ("third sacral nerve",          "ILX:0791971",      "FMA:6425"),
    ("third thoracic nerve",        "ILX:0739246",      "FMA:6171"),
    ("thoracic nerve",              "UBERON:0003726",   "FMA:5860", "ILX:0735852"),
    ("twelfth thoracic nerve",      "ILX:0791046",      "FMA:6167"),
    ("ventral root of spinal cord", "UBERON:0002260",   "FMA:5979", "ILX:0724498")
    ]

def get_spinal_nerve_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in spinal_nerve_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Spinal nerve annotation term '" + name + "' not found.")
