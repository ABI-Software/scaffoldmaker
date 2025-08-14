"""
Common resource for stomach annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stomach_terms = [
    ("body of stomach", "UBERON:0001161", "ILX:0724929", " FMA:14560"),
    ("body-antrum junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793127"),
    ("body-antrum junction along the greater curvature on luminal surface", "ILX:0793128"),
    ("body-antrum junction along the greater curvature on serosa", "ILX:0793087"),
    ("cardia of stomach", "UBERON:0001162", "ILX:0729096", " FMA:14561"),
    ("circular muscle layer of stomach", "ILX:0774731"),
    ("circular-longitudinal muscle interface of body of stomach along the gastric-omentum attachment", "ILX:0793088"),
    ("circular-longitudinal muscle interface of dorsal stomach", "ILX:0793089"),
    ("circular-longitudinal muscle interface of fundus of stomach along the greater curvature", "ILX:0793092"),
    ("circular-longitudinal muscle interface of gastroduodenal junction", "ILX:0793093"),
    ("circular-longitudinal muscle interface of pyloric antrum along the greater curvature", "ILX:0793135"),
    ("circular-longitudinal muscle interface of pyloric antrum along the lesser curvature", "ILX:0793136"),
    ("circular-longitudinal muscle interface of pyloric canal along the greater curvature", "ILX:0793137"),
    ("circular-longitudinal muscle interface of pyloric canal along the lesser curvature", "ILX:0793138"),
    ("circular-longitudinal muscle interface of stomach", "ILX:0793096"),
    ("circular-longitudinal muscle interface of ventral stomach", "ILX:0793097"),
    ("distal point of lower esophageal sphincter serosa on the greater curvature of stomach", "ILX:0793179"),
    ("distal point of lower esophageal sphincter serosa on the lesser curvature of stomach", "ILX:0793180"),
    ("dorsal stomach", "ILX:0793086"),
    ("duodenum part of stomach", "UBERON:0002114", "ILX:0726125", "FMA:7206"),
    ("esophagogastric junction", "UBERON:0007650", "ILX:0733910", "FMA: 9434"),
    ("esophagogastric junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793098"),
    ("esophagogastric junction along the greater curvature on luminal surface", "ILX:0793099"),
    ("esophagogastric junction along the greater curvature on serosa", "ILX:0793100"),
    ("esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface", "ILX:0793101"),
    ("esophagogastric junction along the lesser curvature on luminal surface", "ILX:0793102"),
    ("esophagogastric junction along the lesser curvature on serosa", "ILX:0793103"),
    ("esophagus part of stomach", "UBERON:0001043", "ILX:0735017", "FMA:7131"),
    ("forestomach-glandular stomach junction", "UBERON:0012270", "ILX:0729974"),
    ("fundus of stomach", "UBERON:0001160", "ILX:0724443", " FMA:14559"),
    ("fundus-body junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793104"),
    ("fundus-body junction along the greater curvature on luminal surface", "ILX:0793105"),
    ("fundus-body junction along the greater curvature on serosa", "ILX:0793106"),
    ("gastroduodenal junction", "UBERON:0012650", "ILX:0725406", "FMA:17046"),
    ("gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793107"),
    ("gastroduodenal junction along the greater curvature on luminal surface", "ILX:0793108"),
    ("gastroduodenal junction along the greater curvature on serosa", "ILX:0793109"),
    ("gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface", "ILX:0793110"),
    ("gastroduodenal junction along the lesser curvature on luminal surface", "ILX:0793111"),
    ("gastroduodenal junction along the lesser curvature on serosa", "ILX:0793112"),
    ("greater curvature of stomach", "UBERON:0001164", "ILX:0724395", "FMA:14574"),
    ("lesser curvature of stomach", "UBERON:0001163", "ILX:0733753", "FMA:14572"),
    ("limiting ridge at the greater curvature on the circular-longitudinal muscle interface", "ILX:0793113"),
    ("limiting ridge at the greater curvature on the luminal surface", "ILX:0793114"),
    ("limiting ridge at the greater curvature on serosa", "ILX:0793115"),
    ("limiting ridge on circular-longitudinal muscle interface", "ILX:0793116"),
    ("limiting ridge on luminal surface", "ILX:0793117"),
    ("limiting ridge on serosa", "ILX:0793118"),
    ("longitudinal muscle layer of stomach", "ILX:0772619"),
    ("luminal surface of body of stomach", "ILX:0793119"),
    ("luminal surface of cardia of stomach", "ILX:0793120"),
    ("luminal surface of fundus of stomach", "ILX:0793123"),
    ("luminal surface of pyloric antrum", "ILX:0793124"),
    ("luminal surface of pyloric canal", "ILX:0793125"),
    ("luminal surface of stomach", "ILX:0793126"),
    ("mucosa of stomach", "UBERON:0001199", "ILX:0736669", "FMA:14907"),
    ("pyloric antrum", "UBERON:0001165", "ILX:0728672", " FMA:14579"),
    ("pyloric canal", "UBERON:0008858", "ILX:0735898", "FMA:14580"),
    ("serosa of body of stomach", "ILX:0771402"),
    ("serosa of cardia of stomach", "ILX:0776646"),
    ("serosa of fundus of stomach", "UBERON:0012503", "ILX:0726906", "FMA:17073"),
    ("serosa of pyloric antrum", "ILX:0777005"),
    ("serosa of pyloric canal", "ILX:0775898"),
    ("serosa of stomach", "UBERON:0001201", "ILX:0735818", "FMA:14914"),
    ("stomach", "UBERON:0000945", "ILX:0736697", "FMA:7148"),
    ("submucosa of stomach", "UBERON:0001200", "ILX:0732950", "FMA:14908"),
    ("ventral stomach", "ILX:0793085")
    ]


def get_stomach_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in stomach_terms:
        if name in term:
            return (term[0], term[1])
    raise NameError("Stomach annotation term '" + name + "' not found.")
