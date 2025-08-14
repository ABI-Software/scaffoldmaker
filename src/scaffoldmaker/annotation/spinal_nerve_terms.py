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
    ("T1 spinal nerve", "ILX:0738397", "FMA:6037"),
    ("T2 spinal nerve", "ILX:0739243", "FMA:6169"),
    ("T3 spinal nerve", "ILX:0739246", "FMA:6171"),
    ("T4 spinal nerve", "ILX:0784755", "FMA:6290"),
    ("T5 spinal nerve", "ILX:0792589", "FMA:6293"),
    ("T6 spinal nerve", "ILX:0791122", "FMA:6296"),
    ("T7 spinal nerve", "ILX:0789008", "FMA:6299"),
    ("T8 spinal nerve", "ILX:0792754", "FMA:6302"),
    ("T9 spinal nerve", "ILX:0789545", "FMA:6306"),
    ("T10 spinal nerve", "ILX:0792933", "FMA:6309"),
    ("T11 spinal nerve", "ILX:0791454", "FMA:6312"),
    ("T12 spinal nerve", "ILX:0791046", "FMA:6167"),
    ("L1 spinal nerve", "ILX:0788972", "FMA:6172"),
    ("L2 spinal nerve", "ILX:0786362", "FMA:6176"),
    ("L3 spinal nerve", "ILX:0789374", "FMA:6414"),
    ("L4 spinal nerve", "ILX:0788896", "FMA:6415"),
    ("L5 spinal nerve", "ILX:0787261", "FMA:6416"),
    ("S1 spinal nerve", "ILX:0790504", "FMA:6423"),
    ("S2 spinal nerve", "ILX:0787528", "FMA:6424"),
    ("S3 spinal nerve", "ILX:0791971", "FMA:6425"),
    ("S4 spinal nerve", "ILX:0792327", "FMA:6426"),
    ("S5 spinal nerve", "ILX:0787432", "FMA:6427"),
    ("cervical spinal nerve", "UBERON:0000962", "ILX:0733509", "FMA:5859"),
    ("coccygeal nerve", "UBERON:0009629", "ILX:0723795", "FMA:5863"),
    ("dorsal root ganglion", "UBERON:0000044", "ILX:0103471", "FMA:5888"),
    ("dorsal root of spinal cord", "UBERON:0002261", "ILX:0726708", "FMA:5980"), 
    ("lumbar nerve", "UBERON:0009624", "ILX:0735627", "FMA:5861"),
    ("sacral nerve", "UBERON:0009625", "ILX:0726145", "FMA:5862"),
    ("spinal nerve", "UBERON:0001780", "ILX:0729912", "FMA:5858"),
    ("thoracic nerve", "UBERON:0003726", "ILX:0735852", "FMA:5860"),
    ("ventral root of spinal cord", "UBERON:0002260", "ILX:0724498", "FMA:5979")
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
