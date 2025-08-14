"""
Common resource for esophagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
esophagus_terms = [
    ("abdominal part of esophagus", "UBERON:0035177", "ILX:0735274", "FMA:9397"),
    ("circular-longitudinal muscle interface of esophagus along the cut margin", "ILX:0793091"),
    ("cervical part of esophagus", "UBERON:0035450", "ILX:0734725", "FMA:9395"),
    ("distal point of lower esophageal sphincter serosa on the greater curvature of stomach", "ILX:0793179"),
    ("distal point of lower esophageal sphincter serosa on the lesser curvature of stomach", "ILX:0793180"),
    ("esophagus", "UBERON:0001043", "ILX:0735017", "FMA:7131"),
    ("esophagus mucosa", "UBERON:0002469", "ILX:0725079", "FMA:62996"),
    ("esophagus smooth muscle circular layer", "UBERON:0009960", "ILX:0727608", "FMA:67605"),
    ("esophagus smooth muscle longitudinal layer", "UBERON:0009961", "ILX:0735992", "FMA:63573"),
    ("luminal surface of esophagus", "ILX:0793122"),
    ("proximodorsal midpoint on serosa of upper esophageal sphincter", "ILX:0793181"),
    ("proximoventral midpoint on serosa of upper esophageal sphincter", "ILX:0793182"),
    ("serosa of esophagus", "UBERON:0001975", "ILX:0725745", "FMA:63057"),
    ("submucosa of esophagus", "UBERON:0001972", "ILX:0728662", "FMA:62997"),
    ("thoracic part of esophagus", "UBERON:0035216", "ILX:0732442", "FMA:9396"),
    ]

def get_esophagus_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in esophagus_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Esophagus annotation term '" + name + "' not found.")
