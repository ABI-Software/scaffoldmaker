"""
Common resource for stomach annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stomach_terms = [
    ( "body of stomach", "UBERON:0001161", " FMA:14560", "ILX:0724929"),
    ( "body-antrum junction along the greater curvature on circular-longitudinal muscle interface", "None"),
    ( "body-antrum junction along the greater curvature on luminal surface", "None"),
    ( "body-antrum junction along the greater curvature on serosa", "None"),
    ( "cardia of stomach", "UBERON:0001162", " FMA:14561", "ILX:0729096"),
    ( "circular muscle layer of stomach", "ILX:0774731"),
    ( "circular-longitudinal muscle interface of body of stomach along the gastric-omentum attachment", "None"),
    ( "circular-longitudinal muscle interface of dorsal stomach", "None"),
    ( "circular-longitudinal muscle interface of duodenum along the gastric-omentum attachment", "None"),
    ( "circular-longitudinal muscle interface of esophagus along the cut margin", "None"),
    ( "circular-longitudinal muscle interface of fundus of stomach along the gastric-omentum attachment", "None"),
    ( "circular-longitudinal muscle interface of gastroduodenal junction", "None"),
    ( "circular-longitudinal muscle interface of pyloric antrum along the gastric-omentum attachment", "None"),
    ( "circular-longitudinal muscle interface of pyloric canal along the gastric-omentum attachment", "None"),
    ( "circular-longitudinal muscle interface of stomach", "None"),
    ( "circular-longitudinal muscle interface of ventral stomach", "None"),
    ( "dorsal stomach", "ILX:0793086"),
    ( "duodenum", "UBERON:0002114", " FMA:7206", "ILX:0726125"),
    ( "esophagogastric junction", "UBERON:0007650", "FMA: 9434", "ILX:0733910"),
    ( "esophagogastric junction along the greater curvature on circular-longitudinal muscle interface", "None"),
    ( "esophagogastric junction along the greater curvature on luminal surface", "None"),
    ( "esophagogastric junction along the greater curvature on serosa", "None"),
    ( "esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface", "None"),
    ( "esophagogastric junction along the lesser curvature on luminal surface", "None"),
    ( "esophagogastric junction along the lesser curvature on serosa", "None"),
    ( "esophagus", "UBERON:0001043", "FMA:7131", "ILX:0735017"),
    ( "esophagus mucosa", "UBERON:0002469", "FMA:62996", "ILX:0725079"),
    ( "esophagus smooth muscle circular layer", "UBERON:0009960", "FMA:67605", "ILX:0735992"),
    ( "esophagus smooth muscle longitudinal layer", "UBERON:0009961", "FMA:63573", "ILX:0727608"),
    ( "forestomach-glandular stomach junction", "UBERON:0012270", "ILX:0729974"),
    ( "fundus of stomach", "UBERON:0001160", " FMA:14559", "ILX:0724443"),
    ( "fundus-body junction along the greater curvature on circular-longitudinal muscle interface", "None"),
    ( "fundus-body junction along the greater curvature on luminal surface", "None"),
    ( "fundus-body junction along the greater curvature on serosa", "None"),
    ( "gastroduodenal junction", "UBERON:0012650", "FMA:17046", "ILX:0725406"),
    ( "gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface", "None"),
    ( "gastroduodenal junction along the greater curvature on luminal surface", "None"),
    ( "gastroduodenal junction along the greater curvature on serosa", "None"),
    ( "gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface", "None"),
    ( "gastroduodenal junction along the lesser curvature on luminal surface", "None"),
    ( "gastroduodenal junction along the lesser curvature on serosa", "None"),
    ( "limiting ridge along the greater curvature on circular-longitudinal muscle interface", "None"),
    ( "limiting ridge along the greater curvature on luminal surface", "None"),
    ( "limiting ridge along the greater curvature on serosa", "None"),
    ( "limiting ridge on circular-longitudinal muscle interface", "None"),
    ( "limiting ridge on luminal surface", "None"),
    ( "limiting ridge on serosa", "None"),
    ( "longitudinal muscle layer of stomach", "ILX:0772619"),
    ( "luminal surface of body of stomach", "None"),
    ( "luminal surface of cardia of stomach", "None"),
    ( "luminal surface of duodenum", "None"),
    ( "luminal surface of esophagus", "None"),
    ( "luminal surface of fundus of stomach", "None"),
    ( "luminal surface of pyloric antrum", "None"),
    ( "luminal surface of pyloric canal", "None"),
    ( "luminal surface of stomach", "None"),
    ( "mucosa of stomach", "UBERON:0001199", "FMA:14907", "ILX:0736669"),
    ( "pyloric antrum", "UBERON:0001165", " FMA:14579", "ILX:0728672"),
    ( "pyloric canal", "UBERON:0008858", "FMA:14580", "ILX:0735898"),
    ( "serosa of body of stomach", "ILX:0771402"),
    ( "serosa of cardia of stomach", "ILX:0776646"),
    ( "serosa of duodenum", "UBERON:0003336", "FMA:14948", "ILX:0732373"),
    ( "serosa of esophagus", "UBERON:0001975", "FMA:63057", "ILX:0725745"),
    ( "serosa of fundus of stomach", "UBERON:0012503", "FMA:17073", "ILX:0726906"),
    ( "serosa of pyloric antrum", "ILX:0777005"),
    ( "serosa of pyloric canal", "ILX:0775898"),
    ( "serosa of stomach", "UBERON:0001201", "FMA:14914", "ILX:0735818"),
    ( "stomach", "UBERON:0000945", "FMA:7148", "ILX:0736697"),
    ( "submucosa of esophagus", "UBERON:0001972", "FMA:62997", "ILX:0728662"),
    ( "submucosa of stomach", "UBERON:0001200", "FMA:14908", "ILX:0732950"),
    ( "ventral stomach", "ILX:0793085"),
    ]

def get_stomach_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in stomach_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Stomach annotation term '" + name + "' not found.")
