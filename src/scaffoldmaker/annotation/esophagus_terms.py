"""
Common resource for esophagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
esophagus_terms = [
    ("abdominal part of esophagus", "UBERON:0035177", "FMA:9397", "ILX:0735274"),
    ("cervical part of esophagus", "UBERON:0035450", "FMA:9395", "ILX:0734725"),
    ("distal point of lower esophageal sphincter serosa on the greater curvature of stomach", "ILX:0793179"),
    ("distal point of lower esophageal sphincter serosa on the lesser curvature of stomach", "ILX:0793180"),
    ("esophagus", "UBERON:0001043", "FMA:7131", "ILX:0735017"),
    ("esophagus mucosa", "UBERON:0002469", "FMA:62996", "ILX:0725079"),
    ("esophagus smooth muscle circular layer", "UBERON:0009960", "FMA:67605", "ILX:0727608"),
    ("esophagus smooth muscle longitudinal layer", "UBERON:0009961", "FMA:63573", "ILX:0735992"),
    ("proximodorsal midpoint on serosa of upper esophageal sphincter", "ILX:0793181"),
    ("proximoventral midpoint on serosa of upper esophageal sphincter", "ILX:0793182"),
    ("serosa of esophagus", "UBERON:0001975", "FMA:63057", "ILX:0725745"),
    ("submucosa of esophagus", "UBERON:0001972", "FMA:62997", "ILX:0728662"),
    ("thoracic part of esophagus", "UBERON:0035216", "FMA:9396", "ILX:0732442"),
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
