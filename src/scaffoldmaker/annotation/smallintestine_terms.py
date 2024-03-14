"""
Common resource for small intestine annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
smallintestine_terms = [
    ("circular-longitudinal muscle interface of first segment of the duodenum along the gastric-omentum attachment", "ILX:0793090"),
    ("circular muscle layer of ileum", "ILX:0776561"),
    ("circular muscle layer of small intestine", "ILX:0772669"),
    ("duodenum", "UBERON:0002114", "FMA:7206", "ILX:0726125"),
    ("ileocecal junction", "UBERON:0001073", "FMA:11338", "ILX:0730012"),
    ("ileum", "UBERON:0002116", "FMA:7208", "ILX:0728151"),
    ("jejunum", "UBERON:0002115", "FMA:7207", "ILX:0724224"),
    ("longitudinal muscle layer of ileum", "ILX:0770304"),
    ("longitudinal muscle layer of small intestine", "ILX:0772125"),
    ("luminal surface of duodenum", "ILX:0793121"),
    ("mucosa of ileum", "ILX:0770578"),
    ("mucosa of small intestine", "UBERON:0001204", "FMA:14933", "ILX:0770578"),
    ("serosa of duodenum", "UBERON:0003336", "FMA:14948", "ILX:0732373"),
    ("serosa of ileum", "ILX:0774472"),
    ("serosa of small intestine", "UBERON:0001206", "FMA:14938", "ILX:0727465"),
    ("small intestine", "UBERON:0002108", "FMA:7200", "ILX:0726770"),
    ("submucosa of ileum", "UBERON:0004946", "FMA:14957", "ILX:0734297"),
    ("submucosa of small intestine", "UBERON:0001205", "FMA:14934", "ILX:0735609")
    ]


def get_smallintestine_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in smallintestine_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Small intestine annotation term '" + name + "' not found.")
