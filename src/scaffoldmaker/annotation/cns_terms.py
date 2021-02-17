"""
Common resource for CNS (brain, brainstem, and spinal cord) annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
spinalcord_terms = [
                     ("Cervical spinal cord", "UBERON:0002726", "FMA:71166", "ILX:0102009"),
                     ("Thoracic spinal cord", "UBERON:0003038", "FMA:71167", "ILX:0111710"),
                     ("Lumbar spinal cord", "UBERON:0002792", "FMA:71168", "ILX:0106391"),
                     ("Lumbar spinal cord", "UBERON:0005843", "FMA:256623", "ILX:0110295"),
                     ("caudal segment of spinal cord", "UBERON:0005845", "ILX:0725272")
                 ],


def get_spinalcord_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in spinalcord_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Spinal cord annotation term '" + name + "' not found.")
