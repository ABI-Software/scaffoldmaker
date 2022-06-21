"""
Common resource for lungs annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
lung_terms = [
    ("anterior border of left lung", "ILX:0793130"),
    ("anterior border of right lung", "ILX:0793129"),
    ("anterior mediastinum of left lung", "ILX:0793183"),
    ("anterior mediastinum of right lung", "ILX:0793184"),
    ("apex of left lung", "ILX:0778112"),
    ("apex of right lung", "ILX:0778113"),
    ("apex of right lung accessory lobe", "ILX:0778119"),
    ("base of left lung surface", "ILX:0793187"),
    ("base of lower lobe of left lung surface", "ILX:0793577"),
    ("base of upper lobe of left lung surface", "ILX:0793578"),
    ("base of right lung accessory lobe surface", "ILX:0793189"),
    ("base of right lung surface", "ILX:0793188"),
    ("base of lower lobe of right lung surface", "ILX:0793579"),
    ("base of middle lobe of right lung surface", "ILX:0793580"),
    ("dorsal apex of right lung accessory lobe", "ILX:0793603"),
    ("dorsal base of right lung accessory lobe", "ILX:0778125"),
    ("dorsal base of left lung", "ILX:0778126"),
    ("dorsal base of right lung", "ILX:0778127"),
    ("dorsal surface of right lung accessory lobe", "ILX:0793609"),
    ("horizontal fissure of right lung", "ILX:0746327"),
    ("horizontal fissure of lower lobe of right lung", "ILX:0793581"),
    ("horizontal fissure of middle lobe of right lung", "ILX:0793582"),
    ("horizontal fissure of upper lobe of right lung", "ILX:0793583"),
    ("laterodorsal tip of middle lobe of right lung", "ILX:0778124"),
    ("lateral surface of left lung", "ILX:0793601"),
    ("lateral surface of right lung", "ILX:0793602"),
    ("lateral surface of lower lobe of left lung", "ILX:0793584"),
    ("lateral surface of lower lobe of right lung", "ILX:0793585"),
    ("lateral surface of middle lobe of right lung", "ILX:0793586"),
    ("lateral surface of upper lobe of left lung", "ILX:0793587"),
    ("lateral surface of upper lobe of right lung", "ILX:0793588"),
    ("left dorsal base of right lung accessory lobe", "ILX:0793607"),
    ("left lung", "UBERON:0002168", "ILX:0733450"),
    ("left lung surface", "ILX:0793186"),
    ("left surface of right lung accessory lobe", "ILX:0793611"),
    ("left ventral base of right lung accessory lobe", "ILX:0793605"),
    ("lower lobe of left lung", "UBERON:0008953", "ILX:0735534"),
    ("lower lobe of left lung surface", "ILX:0793192"),
    ("lower lobe of right lung", "UBERON:0002171", "ILX:0725712"),
    ("lower lobe of right lung surface", "ILX:0793191"),
    ("lung", "UBERON:0002048", "ILX:0726937"),
    ("medial surface of left lung", "ILX:0793599"),
    ("medial surface of right lung", "ILX:0793600"),
    ("medial surface of lower lobe of left lung", "ILX:0793589"),
    ("medial surface of lower lobe of right lung", "ILX:0793590"),
    ("medial surface of middle lobe of right lung", "ILX:0793591"),
    ("medial surface of upper lobe of left lung", "ILX:0793592"),
    ("medial surface of upper lobe of right lung", "ILX:0793593"),
    ("middle lobe of right lung", "UBERON:0002174", "ILX:0733737"),
    ("middle lobe of right lung surface", "ILX:0793193"),
    ("medial base of left lung", "ILX:0778120"),
    ("medial base of right lung", "ILX:0778121"),
    ("anterior mediastinum of left lung", "UBERON:0008820", "ILX:0725455"),
    ("anterior mediastinum of right lung", "UBERON:0008820", "ILX:0725455"),
    ("oblique fissure of left lung", "ILX:0778115"),
    ("oblique fissure of lower lobe of left lung", "ILX:0793594"),
    ("oblique fissure of upper lobe of left lung", "ILX:0793595"),
    ("oblique fissure of right lung", "ILX:0778114"),
    ("oblique fissure of lower lobe of right lung", "ILX:0793596"),
    ("oblique fissure of middle lobe of right lung", "ILX:0793597"),
    ("oblique fissure of upper lobe of right lung", "ILX:0793598"),
    ("right dorsal base of right lung accessory lobe", "ILX:0793608"),
    ("right lung", "UBERON:0002167", "ILX:0729582"),
    ("right lung surface", "ILX:0793185"),
    ("right lung accessory lobe", "UBERON:0004890", "ILX:0728696"),
    ("right lung accessory lobe surface", "ILX:0793190"),
    ("right surface of right lung accessory lobe", "ILX:0793612"),
    ("right ventral base of right lung accessory lobe", "ILX:0793606"),
    ("upper lobe of left lung", "UBERON:0008952", "ILX:0735339"),
    ("upper lobe of left lung surface", "ILX:0793194"),
    ("upper lobe of right lung", "UBERON:0002170", "ILX:0728821"),
    ("upper lobe of right lung surface", "ILX:0793195"),
    ("ventral apex of right lung accessory lobe", "ILX:0793604"),
    ("ventral base of right lung accessory lobe", "ILX:0778123"),
    ("ventral base of left lung", "ILX:0778118"),
    ("ventral base of right lung", "ILX:0778122"),
    ("ventral surface of right lung accessory lobe", "ILX:0793610")

]

def get_lung_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in lung_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Lung annotation term '" + name + "' not found.")
