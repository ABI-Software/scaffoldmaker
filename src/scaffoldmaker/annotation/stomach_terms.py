"""
Common resource for stomach annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stomach_terms = [
    ("body of stomach", "UBERON:0001161", " FMA:14560", "ILX:0724929"),
    ("body-antrum junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793127"),
    ("body-antrum junction along the greater curvature on luminal surface", "ILX:0793128"),
    ("body-antrum junction along the greater curvature on serosa", "ILX:0793087"),
    ("cardia of stomach", "UBERON:0001162", " FMA:14561", "ILX:0729096"),
    ("circular muscle layer of stomach", "ILX:0774731"),
    ("circular-longitudinal muscle interface of body of stomach along the gastric-omentum attachment", "ILX:0793088"),
    ("circular-longitudinal muscle interface of dorsal stomach", "ILX:0793089"),
    ("circular-longitudinal muscle interface of first segment of the duodenum along the gastric-omentum attachment", "ILX:0793090"),
    ("circular-longitudinal muscle interface of esophagus along the cut margin", "ILX:0793091"),
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
    ("duodenum", "UBERON:0002114", " FMA:7206", "ILX:0726125"),
    ("esophagogastric junction", "UBERON:0007650", "FMA: 9434", "ILX:0733910"),
    ("esophagogastric junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793098"),
    ("esophagogastric junction along the greater curvature on luminal surface", "ILX:0793099"),
    ("esophagogastric junction along the greater curvature on serosa", "ILX:0793100"),
    ("esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface", "ILX:0793101"),
    ("esophagogastric junction along the lesser curvature on luminal surface", "ILX:0793102"),
    ("esophagogastric junction along the lesser curvature on serosa", "ILX:0793103"),
    ("esophagus", "UBERON:0001043", "FMA:7131", "ILX:0735017"),
    ("esophagus mucosa", "UBERON:0002469", "FMA:62996", "ILX:0725079"),
    ("esophagus smooth muscle circular layer", "UBERON:0009960", "FMA:67605", "ILX:0735992"),
    ("esophagus smooth muscle longitudinal layer", "UBERON:0009961", "FMA:63573", "ILX:0727608"),
    ("forestomach-glandular stomach junction", "UBERON:0012270", "ILX:0729974"),
    ("fundus of stomach", "UBERON:0001160", " FMA:14559", "ILX:0724443"),
    ("fundus-body junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793104"),
    ("fundus-body junction along the greater curvature on luminal surface", "ILX:0793105"),
    ("fundus-body junction along the greater curvature on serosa", "ILX:0793106"),
    ("gastroduodenal junction", "UBERON:0012650", "FMA:17046", "ILX:0725406"),
    ("gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface", "ILX:0793107"),
    ("gastroduodenal junction along the greater curvature on luminal surface", "ILX:0793108"),
    ("gastroduodenal junction along the greater curvature on serosa", "ILX:0793109"),
    ("gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface", "ILX:0793110"),
    ("gastroduodenal junction along the lesser curvature on luminal surface", "ILX:0793111"),
    ("gastroduodenal junction along the lesser curvature on serosa", "ILX:0793112"),
    ("greater curvature of stomach", "UBERON:0001164", "FMA:14574", "ILX:0724395"),
    ("lesser curvature of stomach", "UBERON:0001163", "FMA:14572", "ILX:0733753"),
    ("limiting ridge at the greater curvature on the circular-longitudinal muscle interface", "ILX:0793113"),
    ("limiting ridge at the greater curvature on the luminal surface", "ILX:0793114"),
    ("limiting ridge at the greater curvature on serosa", "ILX:0793115"),
    ("limiting ridge on circular-longitudinal muscle interface", "ILX:0793116"),
    ("limiting ridge on luminal surface", "ILX:0793117"),
    ("limiting ridge on serosa", "ILX:0793118"),
    ("longitudinal muscle layer of stomach", "ILX:0772619"),
    ("luminal surface of body of stomach", "ILX:0793119"),
    ("luminal surface of cardia of stomach", "ILX:0793120"),
    ("luminal surface of duodenum", "ILX:0793121"),
    ("luminal surface of esophagus", "ILX:0793122"),
    ("luminal surface of fundus of stomach", "ILX:0793123"),
    ("luminal surface of pyloric antrum", "ILX:0793124"),
    ("luminal surface of pyloric canal", "ILX:0793125"),
    ("luminal surface of stomach", "ILX:0793126"),
    ("mucosa of stomach", "UBERON:0001199", "FMA:14907", "ILX:0736669"),
    ("pyloric antrum", "UBERON:0001165", " FMA:14579", "ILX:0728672"),
    ("pyloric canal", "UBERON:0008858", "FMA:14580", "ILX:0735898"),
    ("serosa of body of stomach", "ILX:0771402"),
    ("serosa of cardia of stomach", "ILX:0776646"),
    ("serosa of duodenum", "UBERON:0003336", "FMA:14948", "ILX:0732373"),
    ("serosa of esophagus", "UBERON:0001975", "FMA:63057", "ILX:0725745"),
    ("serosa of fundus of stomach", "UBERON:0012503", "FMA:17073", "ILX:0726906"),
    ("serosa of pyloric antrum", "ILX:0777005"),
    ("serosa of pyloric canal", "ILX:0775898"),
    ("serosa of stomach", "UBERON:0001201", "FMA:14914", "ILX:0735818"),
    ("stomach", "UBERON:0000945", "FMA:7148", "ILX:0736697"),
    ("submucosa of esophagus", "UBERON:0001972", "FMA:62997", "ILX:0728662"),
    ("submucosa of stomach", "UBERON:0001200", "FMA:14908", "ILX:0732950"),
    ("ventral stomach", "ILX:0793085"),
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
