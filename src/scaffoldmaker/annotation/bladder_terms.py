"""
Common resource for bladder annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
bladder_terms = [
    ("bladder lumen", "UBERON:0009958"),
    ("dome of the bladder", "ILX:0738433"),
    ("dorsal part of urinary bladder lumen", "ILX:0777128"),
    ("dorsal part of lumen of body of urinary bladder", "ILX:0739250"),
    ("dorsal part of lumen of neck of urinary bladder", "ILX:0739255"),
    ("dorsal part of lumen of urethra", "ILX:0739260"),
    ("dorsal part of bladder", "ILX:0778155"),
    ("dorsal part of urethra", "ILX:0739258"),
    ("dorsal part of serosa of body of urinary bladder", "ILX:0739278"),
    ("dorsal part of serosa of neck of urinary bladder", "ILX:0739280"),
    ("dorsal part of serosa of urethra", "ILX:0739283"),
    ("dorsal part of serosa of urinary bladder", "ILX:0739248"),
    ("left hemi-bladder", "ILX:0793661"),
    ("lumen of body of urinary bladder", "ILX:0739252"),
    ("lumen of neck of urinary bladder", "ILX:0739256"),
    ("lumen of urethra", "UBERON:0010390", "ILX:0736762"),
    ("neck of urinary bladder", "UBERON:0001258", "FMA:15912"),
    ("right hemi-bladder", "ILX:0793662"),
    ("serosa of body of urinary bladder", "ILX:0739276"),
    ("serosa of neck of urinary bladder", "ILX:0739277"),
    ("serosa of urethra", "ILX:0739282"),
    ("serosa of urinary bladder", "UBERON:0001260"),
    ("ureter", "UBERON:0000056", "ILX:0728080"),
    ("urethra", "UBERON:0000057", "ILX:0733022"),
    ("urinary bladder", "UBERON:0001255", "FMA:15900"),
    ("ventral part of urinary bladder lumen", "ILX:0777129"),
    ("ventral part of lumen of body of urinary bladder", "ILX:0739251"),
    ("ventral part of lumen of neck of urinary bladder", "ILX:0739257"),
    ("ventral part of lumen of urethra", "ILX:0739261"),
    ("ventral part of serosa of body of urinary bladder", "ILX:0739279"),
    ("ventral part of serosa of neck of urinary bladder", "ILX:0739281"),
    ("ventral part of serosa of urethra", "ILX:0739306"),
    ("ventral part of serosa of urinary bladder", "ILX:0739249"),
    ("ventral part of bladder", "ILX:0778156 "),
    ("ventral part of urethra", "ILX:0739259"),
    ("apex of urinary bladder", "ILX:0774405"),
    ("left ureter junction with bladder", "ILX:0778142"),
    ("right ureter junction with bladder", "ILX:0778143"),
    ("urethra junction of dorsal bladder neck", "ILX:0778154"),
    ("urethra junction of ventral bladder neck", "ILX:0738410")
]


def get_bladder_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in bladder_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Bladder annotation term '" + name + "' not found.")
