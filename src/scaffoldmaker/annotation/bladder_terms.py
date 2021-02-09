"""
Common resource for bladder annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
bladder_terms = [
    ("urinary bladder", "UBERON:0001255", "FMA:15900"),
    ("neck of urinary bladder", "UBERON:0001258", "FMA:15912"),
    ("Dome of the Bladder", "ILX:0738433"),
    ("serosa of urinary bladder", "UBERON:0001260"),
    ("bladder lumen", "UBERON:0009958"),
    ("serosa of body of urinary bladder", "ILX:0739276"),
    ("lumen of body of urinary bladder", "ILX:0739252"),
    ("serosa of neck of urinary bladder", "ILX:0739277"),
    ("lumen of neck of urinary bladder", "ILX:0739256"),
    ("urethra", "UBERON:0000057", "ILX:0733022"),
    ("lumen of urethra", "UBERON:0010390", "ILX:0736762"),
    ("serosa of urethra", "ILX:0739282"),
    ("ureter", "UBERON:0000056", "ILX:0728080"),
    ("Dorsal part of serosa of urinary bladder", "ILX:0739248"),
    ("Ventral part of serosa of urinary bladder", "ILX:0739249"),
    ("dorsal part of bladder lumen", "None"),
    ("ventral part of bladder lumen", "None"),
    ("dorsal part of serosa of body of urinary bladder", "ILX:0739278"),
    ("ventral part of serosa of body of urinary bladder", "ILX:0739279"),
    ("Dorsal part of lumen of body of urinary bladder", "ILX:0739250"),
    ("Ventral part of lumen of body of urinary bladder", "ILX:0739251"),
    ("dorsal part of serosa of neck of urinary bladder", "ILX:0739280"),
    ("ventral part of serosa of neck of urinary bladder", "ILX:0739281"),
    ("Dorsal part of lumen of neck of urinary bladder", "ILX:0739255"),
    ("Ventral part of lumen of neck of urinary bladder", "ILX:0739257"),
    ("dorsal part of serosa of urethra", "ILX:0739283"),
    ("ventral part of serosa of urethra", "ILX:0739306"),
    ("Dorsal part of lumen of urethra", "ILX:0739260"),
    ("Ventral part of lumen of urethra", "ILX:0739261"),
    ("dorsal part of the bladder", "None"),
    ("ventral part of the bladder", "None"),
    ("Dorsal part of urethra", "ILX:0739258"),
    ("Ventral part of urethra", "ILX:0739259"),
    ("dorsal part of the scaffold", "None"),
    ("ventral part of the scaffold", "None")
]

def get_bladder_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in bladder_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Bladder annotation term '" + name + "' not found.")
