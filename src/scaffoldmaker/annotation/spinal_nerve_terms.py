"""
Common resource for spinal nerve annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
spinal_nerve_terms = [
    ("C1 spinal nerve", "ILX:0794592"),
    ("C2 spinal nerve", "ILX:0794595"),
    ("C3 spinal nerve", "ILX:0794598"),
    ("C4 spinal nerve", "ILX:0794601"),
    ("C5 spinal nerve", "ILX:0794604"),
    ("C6 spinal nerve", "ILX:0794607"),
    ("C7 spinal nerve", "ILX:0794610"),
    ("C8 spinal nerve", "ILX:0794613"),
    ("cervical spinal nerve", "UBERON:0000962", "FMA:5859", "ILX:0733509"),
    ("coccygeal nerve", "UBERON:0009629", "FMA:5863", "ILX:0723795"),
    ("dorsal root ganglion", "UBERON:0000044", "FMA:5888", "ILX:0103471"),
    ("dorsal root of spinal cord", "UBERON:0002261", "FMA:5980", "ILX:0726708"),
    ("eighth thoracic nerve", "FMA:6302", "ILX:0792754"),
    ("eleventh thoracic nerve", "FMA:6312", "ILX:0791454"),
    ("fifth lumbar nerve", "FMA:6416", "ILX:0787261"),
    ("fifth sacral nerve", "FMA:6427", "ILX:0787432"),
    ("fifth thoracic nerve", "FMA:6293", "ILX:0792589"),
    ("first lumbar nerve", "FMA:6172", "ILX:0788972"),
    ("first sacral nerve", "FMA:6423", "ILX:0790504"),
    ("first thoracic nerve", "FMA:6037", "ILX:0738397"),
    ("fourth lumbar nerve", "FMA:6415", "ILX:0788896"),
    ("fourth sacral nerve", "FMA:6426", "ILX:0792327"),
    ("fourth thoracic nerve", "FMA:6290", "ILX:0784755"),
    ("lumbar nerve", "UBERON:0009624", "FMA:5861", "ILX:0735627"),
    ("ninth thoracic nerve", "FMA:6306", "ILX:0789545"),
    ("sacral nerve", "UBERON:0009625", "FMA:5862", "ILX:0726145"),
    ("second lumbar nerve", "FMA:6176", "ILX:0786362"),
    ("second sacral nerve", "FMA:6424", "ILX:0787528"),
    ("second thoracic nerve", "FMA:6169", "ILX:0739243"),
    ("seventh thoracic nerve", "FMA:6299", "ILX:0789008"),
    ("sixth thoracic nerve", "FMA:6296", "ILX:0791122"),
    ("spinal nerve", "UBERON:0001780", "FMA:5858", "ILX:0729912"),
    ("tenth thoracic nerve", "FMA:6309", "ILX:0792933"),
    ("third lumbar nerve", "FMA:6414", "ILX:0789374"),
    ("third sacral nerve", "FMA:6425", "ILX:0791971"),
    ("third thoracic nerve", "FMA:6171", "ILX:0739246"),
    ("thoracic nerve", "UBERON:0003726", "FMA:5860", "ILX:0735852"),
    ("twelfth thoracic nerve", "FMA:6167", "ILX:0791046"),
    ("ventral root of spinal cord", "UBERON:0002260", "FMA:5979", "ILX:0724498")
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
