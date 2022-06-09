"""
Common resource for heart annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
heart_terms = [
    # heart - volume terms
    ("heart", "UBERON:0000948", "FMA:7088"),  # group of the entire heart
    ("epicardial fat", "UBERON:0015129"),  # not used
    ("epicardium", "UBERON:0002348", "FMA:9461"),  # volumetric layer outside myocardium to pericardial cavity
    # heart - surface terms
    ("outer surface of epicardium", "ILX:0793548"),
    ("outer surface of myocardium", "ILX:0793530"),

    # ventricles - volume terms
    ("conus arteriosus", "UBERON:0003983"),
    ("endocardium of left ventricle", "UBERON:0009713", "FMA:9559"),
    ("endocardium of right ventricle", "UBERON:0009712", "FMA:9536"),
    ("heart left ventricle", "UBERON:0002084"),  # the whole left ventricle
    ("heart right ventricle", "UBERON:0002080"),  # the whole right ventricle
    ("interventricular septum", "UBERON:0002094", "FMA:7133"),
    ("left fibrous ring", "FMA:77124"),
    ("left ventricle myocardium", "UBERON:0006566", "FMA:9558"),
    ("right fibrous ring", "FMA:77125"),
    ("right ventricle myocardium", "UBERON:0006567", "FMA:9535"),
    # ventricles - surface terms
    ("luminal surface of left ventricle", "ILX:0793537"),
    ("luminal surface of right ventricle", "ILX:0793538"),
    ("outer surface of myocardium of left ventricle", "ILX:0793533"),
    ("outer surface of myocardium of right ventricle", "ILX:0793534"),

    # atria - volume terms
    ("left atrium endocardium", "UBERON:0034903", "FMA:7286"),
    ("endocardium of left auricle", "UBERON:0011006", "FMA:13236"),
    ("right atrium endocardium", "UBERON:0009129", "FMA:7281"),
    ("endocardium of right auricle", "UBERON:0011007", "FMA:13235"),
    ("epicardium of left auricle", "FMA:13233"),
    ("epicardium of right auricle", "FMA:13232"),
    ("fossa ovalis", "UBERON:0003369", "FMA:9246"),
    ("inferior vena cava", "UBERON:0001072", "FMA:10951", "posterior vena cava"),
    ("inferior vena cava inlet", "ILX:0738358"),
    ("interatrial septum", "UBERON:0002085", "FMA:7108"),
    ("left atrium myocardium", "FMA:7285"),
    ("left auricle", "UBERON:0006630", "FMA:7219"),
    ("left cardiac atrium", "UBERON:0002079"),  # the whole left atrium
    ("left inferior pulmonary vein", "FMA:49913"),
    ("left pulmonary vein", "UBERON:0009030"),
    ("left superior pulmonary vein", "FMA:49916"),
    ("middle pulmonary vein", "ILX:0739222"),  # in mouse, rat, rabbit
    ("pulmonary vein", "UBERON:0002016", "FMA:66643"),
    ("right atrium myocardium", "FMA:7282"),
    ("right auricle", "UBERON:0006631", "FMA:7218"),
    ("right cardiac atrium", "UBERON:0002078"),  # the whole right atrium
    ("right inferior pulmonary vein", "FMA:49911"),
    ("right pulmonary vein", "UBERON:0009032"),
    ("right superior pulmonary vein", "FMA:49914"),
    ("superior vena cava", "UBERON:0001585", "FMA:4720", "anterior vena cava"),
    ("superior vena cava inlet", "ILX:0738367"),
    # atria - surface terms
    ("luminal surface of left atrium", "ILX:0793535"),
    ("luminal surface of right atrium", "ILX:0793536"),
    ("luminal surface of inferior vena cava", "ILX:0793542"),
    ("luminal surface of left pulmonary vein", "ILX:0793539"),
    ("luminal surface of middle pulmonary vein", "ILX:0793540"),
    ("luminal surface of right pulmonary vein", "ILX:0793541"),
    ("luminal surface of superior vena cava", "ILX:0793543"),
    ("outer surface of myocardium of left atrium", "ILX:0793531"),
    ("outer surface of myocardium of right atrium", "ILX:0793532"),

    # arterial valves and great vessels - volume terms
    ("root of aorta", "FMA:3740"),
    ("posterior cusp of aortic valve", "FMA:7253"),
    ("right cusp of aortic valve", "FMA:7252"),
    ("left cusp of aortic valve", "FMA:7251"),
    ("root of pulmonary trunk", "FMA:8612"),
    ("right cusp of pulmonary valve", "FMA:7250"),
    ("anterior cusp of pulmonary valve", "FMA:7249"),
    ("left cusp of pulmonary valve", "FMA:7247"),
    # arterial valves and great vessels - surface terms
    ("luminal surface of aorta", "ILX:0793544"),
    ("luminal surface of pulmonary trunk", "ILX:0793546"),
    ("luminal surface of root of aorta", "ILX:0793545"),
    ("luminal surface of root of pulmonary trunk", "ILX:0793547"),

    # fiducial markers
    ("apex of heart", "UBERON:0002098", "FMA:7164"),
    # point on posterior surface where the four chambers meet on interventricular, A-V and interatrial sulci
    ("crux cordis", "ILX:0777104", "FMA:7220"),
    # point at centre of pulmonary vein ostia on left atrium epicardium, on anterior/ventral side for rodents
    ("left atrium epicardium venous midpoint", "ILX:0778116"),
    # point at centre of inferior & superior vena cavae on right atrium epicardium
    ("right atrium epicardium venous midpoint", "ILX:0778117")
]

def get_heart_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in heart_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Heart annotation term '" + name + "' not found.")
