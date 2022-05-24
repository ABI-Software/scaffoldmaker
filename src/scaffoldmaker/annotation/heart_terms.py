"""
Common resource for heart annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
heart_terms = [
    ( "heart", "UBERON:0000948", "FMA:7088" ),
    # ventricles
    ( "left ventricle myocardium", "UBERON:0006566", "FMA:9558" ),
    ( "right ventricle myocardium", "UBERON:0006567", "FMA:9535" ),
    ( "interventricular septum", "UBERON:0002094", "FMA:7133" ),
    ( "endocardium of left ventricle", "UBERON:0009713", "FMA:9559" ),
    ( "endocardium of right ventricle", "UBERON:0009712", "FMA:9536" ),
    ( "epicardial fat", "UBERON:0015129"),
    ( "epicardium", "UBERON:0002348", "FMA:9461"),
    ( "outer surface of epicardium", "None"),
    #( "epicardium of ventricle", "FMA:12150", "UBERON:0001082" ),
    # ventricles with base
    ( "conus arteriosus", "UBERON:0003983" ),
    ( "left fibrous ring", "FMA:77124" ),
    ( "right fibrous ring", "FMA:77125" ),
    # atria
    ( "left atrium myocardium", "FMA:7285" ),
    ( "right atrium myocardium", "FMA:7282" ),
    ( "endocardium of left atrium", "UBERON:0034903", "FMA:7286" ),
    ( "endocardium of right atrium", "UBERON:0009129", "FMA:7281" ),
    ( "interatrial septum", "UBERON:0002085", "FMA:7108" ),
    ( "fossa ovalis", "UBERON:0003369", "FMA:9246" ),
    ( "left auricle", "UBERON:0006630", "FMA:7219" ),  # uncertain if just the tissue like myocardium
    ( "right auricle", "UBERON:0006631", "FMA:7218" ),  # uncertain if just the tissue like myocardium
    ( "endocardium of left auricle", "UBERON:0011006", "FMA:13236" ),
    ( "endocardium of right auricle", "UBERON:0011007", "FMA:13235" ),
    ( "epicardium of left auricle", "FMA:13233" ),
    ( "epicardium of right auricle", "FMA:13232" ),
    ( "pulmonary vein", "UBERON:0002016", "FMA:66643" ),
    ( "left pulmonary vein", "UBERON:0009030" ),
    ( "left inferior pulmonary vein", "FMA:49913" ),
    ( "left superior pulmonary vein", "FMA:49916" ),
    ( "middle pulmonary vein", "ILX:0739222" ),  # in mouse, rat, rabbit
    ( "right pulmonary vein", "UBERON:0009032" ),
    ( "right inferior pulmonary vein", "FMA:49911" ),
    ( "right superior pulmonary vein", "FMA:49914" ),
    ( "inferior vena cava", "UBERON:0001072", "FMA:10951", "posterior vena cava" ),
    ( "inferior vena cava inlet", "ILX:0738358" ),
    ( "superior vena cava", "UBERON:0001585", "FMA:4720", "anterior vena cava" ),
    ( "superior vena cava inlet", "ILX:0738367" ),
    # arterial root
    ( "root of aorta", "FMA:3740" ),
    ( "posterior cusp of aortic valve", "FMA:7253" ),
    ( "right cusp of aortic valve", "FMA:7252" ),
    ( "left cusp of aortic valve", "FMA:7251" ),
    ( "root of pulmonary trunk", "FMA:8612" ),
    ( "right cusp of pulmonary valve", "FMA:7250" ),
    ( "anterior cusp of pulmonary valve", "FMA:7249" ),
    ( "left cusp of pulmonary valve", "FMA:7247" ),
    # fiducial markers
    ( "apex of heart", "UBERON:0002098", "FMA:7164"),
    ( "crux of heart", "ILX:0777104", "FMA:7220" ),  # point on posterior surface where the four chambers meet on interventricular, atrio-ventricular and interatrial sulci
    ( "left atrium epicardium venous midpoint", "ILX:0778116"),  # point at centre of pulmonary vein ostia on left atrium epicardium, on anterior/ventral side for rodents
    ( "right atrium epicardium venous midpoint", "ILX:0778117")  # point at centre of inferior & superior vena cavae on right atrium epicardium
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
