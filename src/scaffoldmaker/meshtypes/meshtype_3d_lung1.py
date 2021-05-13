'''
Generates 3D lung surface mesh.
'''

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_lung1(Scaffold_base):
    '''
    3D lung scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Lung 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        if parameterSetName == 'Default':
            parameterSetName = 'Mouse 1'
        options['Base parameter set'] = parameterSetName
        options['Refine'] = False
        options['Refine number of elements'] = 4
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Refine',
            'Refine number of elements'
            ]
        return optionNames

    @classmethod
    def checkOptions(cls, options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        for key in [
            'Refine number of elements']:
            if options[key] < 1:
                options[key] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''
        parameterSetName = options['Base parameter set']
        isMouse = 'Mouse 1' in parameterSetName
        isHuman = 'Human 1' in parameterSetName

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fm.findMeshByDimension(3)

        eftfactory = eftfactory_tricubichermite(mesh, None)
        eftRegular = eftfactory.createEftBasic()

        elementtemplateRegular = mesh.createElementtemplate()
        elementtemplateRegular.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateRegular.defineField(coordinates, -1, eftRegular)

        elementtemplateCustom = mesh.createElementtemplate()
        elementtemplateCustom.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Annotation groups
        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        annotationGroups = [leftLungGroup, lungGroup]

        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLungGroup)
        lowerRightLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of right lung"))
        lowerRightLungMeshGroup = lowerRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(lowerRightLungGroup)
        upperRightLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of right lung"))
        upperRightLungMeshGroup = upperRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(upperRightLungGroup)
        middleRightLungGroup = AnnotationGroup(region, get_lung_term("middle lobe of right lung"))
        middleRightLungMeshGroup = middleRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(middleRightLungGroup)

        # Marker points/groups
        leftApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_lung_term("apex of left lung"))
        rightApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_lung_term("apex of right lung"))
        leftDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_lung_term("dorsal base of left lung"))
        rightDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_lung_term("dorsal base of right lung"))
        leftVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_lung_term("ventral base of left lung"))
        rightVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_lung_term("ventral base of right lung"))

        if isHuman:
            # Annotation groups
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(upperLeftLungGroup)

            # Marker points/groups
            leftMedialGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_lung_term("medial base of left lung"))
            rightMedialGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_lung_term("medial base of right lung"))
            rightLateralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_lung_term(
                                                                       "laterodorsal tip of middle lobe of right lung"))

        elif isMouse:
            # Annotation groups
            diaphragmaticLungGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
            diaphragmaticLungMeshGroup = diaphragmaticLungGroup.getMeshGroup(mesh)
            annotationGroups.append(diaphragmaticLungGroup)

            accessoryApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                    get_lung_term("apex of accessory lung"))
            accessoryVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term("ventral base of accessory lung"))
            accessoryDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                      get_lung_term("dorsal base of accessory lung"))

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        cache = fm.createFieldcache()

        # common element field templates
        eftWedgeCollapseXi1_15 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
        eftWedgeCollapseXi1_37 = eftfactory.createEftWedgeCollapseXi1Quadrant([3, 7])
        eftWedgeCollapseXi1_57 = eftfactory.createEftWedgeCollapseXi1Quadrant([5, 7])
        eftWedgeCollapseXi2_56 = eftfactory.createEftWedgeCollapseXi2Quadrant([5, 6])
        eftWedgeCollapseXi2_78 = eftfactory.createEftWedgeCollapseXi2Quadrant([7, 8])
        eftTetCollapseXi1Xi2_82 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(8, 2)
        eftTetCollapseXi1Xi2_63 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(6, 3)

        # common parameters in species
        generateParameters = False
        leftLung = 0
        rightLung = 1

        # The number of the elements in the generic lungs
        # These counts are only values that work for nodeFieldParameters (KEEP THEM FIXED)
        lElementsCount1 = 2
        lElementsCount2 = 4
        lElementsCount3 = 3

        uElementsCount1 = 2
        uElementsCount2 = 4
        uElementsCount3 = 4

        if isHuman:
            #valueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ]
            nodeFieldParameters = [
                (   1, [ [238.96169, 226.06812, -343.43996], [-61.11589,   7.39465,  11.89523], [ -30.76581, -29.91568, 35.20122], [-8.03835,15.14892,94.83977] ] ),
                (   2, [ [285.20366, 190.31182, -350.07243], [-54.32449, -11.37468,  71.07193], [  26.97968, -42.41533,  1.60136], [-27.94843,44.55751,77.87808] ] ),
                (   3, [ [223.34100, 188.44400, -308.23200], [-63.07335,   8.96391,   4.33079], [   2.85232, -41.99340, 26.47692], [ 6.50901,20.99596,51.37272] ] ),
                (   4, [ [176.12078, 187.95290, -334.30467], [-12.21420,  11.49438, -33.57606], [   2.99473, -26.22842, 19.79406], [-6.63367,-13.93528,91.81281] ] ),
                (   5, [ [292.99000, 153.58700, -341.67400], [-28.36359,   4.94768,  56.56986], [  -0.87774, -47.75565, 11.50698], [ 7.56668,36.67408,57.88686] ] ),
                (   6, [ [241.78700, 149.54300, -291.57800], [-71.60526, -13.46081,  38.76136], [  14.61976, -40.54216,  6.56251], [ 0.33257,36.38556,37.23453] ] ),
                (   7, [ [155.35122, 127.69892, -273.82518], [-96.50597, -31.65008,  -1.92059], [  -1.01701, -23.72904,  6.49136], [18.59949,33.07221,32.67652] ] ),
                (   8, [ [279.34600,  98.45500, -327.71700], [-12.81487,  17.66496,  37.55641], [ -21.29912, -41.68125,  8.70509], [31.03294,48.25749,54.30009] ] ),
                (   9, [ [251.88700, 110.97900, -294.25900], [-40.86100,   5.67079,  25.71928], [   0.10022, -35.84751, -6.32927], [ 6.49329,41.85543,52.16744] ] ),
                (  10, [ [203.26300, 108.03400, -281.64700], [-52.69189, -10.80320,  -0.46282], [  36.78814, -19.74531, -6.60712], [16.37854,35.98275,29.77094] ] ),
                (  11, [ [256.41200,  71.15200, -323.52500], [- 9.85470,   9.24288,  22.09696], [ -17.61383, -17.51492, 14.32986], [42.65656,30.14688,73.65583] ] ),
                (  12, [ [243.93500,  80.99900, -302.87200], [-15.00801,  10.36550,  19.00436], [  -7.36567, -20.80071, -0.79453], [30.12027,31.34641,61.43762] ] ),
                (  13, [ [226.62800,  91.70200, -285.89200], [-19.43479,  10.94410,  14.82505], [  13.45744, -17.45727, -6.70345], [19.48569,35.01710,31.24560] ] ),
                (  14, [ [232.07411, 235.69876, -247.76428], [-54.92565, - 7.81817,  16.63557], [ -13.17593, -32.29020, -3.12488], [-4.93492, 4.95235,90.50376] ] ),
                (  15, [ [267.56700, 218.53800, -268.60800], [-33.25224, -18.69676,  24.30420], [  34.37980, -29.86829,-19.35106], [-6.43345,10.47365,82.56681] ] ),
                (  16, [ [227.62600, 202.77300, -250.31600], [-45.70013, -12.31054,  11.60033], [   4.91341, -32.00831, -1.82826], [ 1.96535, 7.35354,63.70444] ] ),
                (  17, [ [178.21000, 194.84000, -246.53300], [-51.96221, - 3.47719,  -3.94552], [ -30.90750, -49.20223,  4.73225], [-2.96190, 1.29501,75.13714] ] ),
                (  18, [ [296.25000, 178.15400, -283.77300], [-44.05805,   6.07426,  39.82616], [  14.21016, -42.15834, -1.64208], [-1.21954,11.62212,56.59275] ] ),
                (  19, [ [240.70600, 174.73100, -251.29800], [-65.09937, -13.18642,  23.37871], [  11.73316, -27.40803,  1.68665], [-2.50679,12.65356,41.95740] ] ),
                (  20, [ [170.03600, 151.29900, -240.51000], [-74.19583, -32.77433,  -1.75436], [  22.63784, -36.07139, -0.93095], [ 9.78386,16.37566,32.19940] ] ),
                (  21, [ [297.50200, 143.35500, -275.67900], [-41.75575,  11.05269,  40.73634], [  -3.98149, -37.43511, 15.98202], [ 4.76819,40.74808,48.88201] ] ),
                (  22, [ [250.97800, 148.43100, -247.19500], [-49.40725, - 1.39965,  14.39268], [   8.88604, -29.27432,  6.34914], [-8.37384,32.64539,41.45804] ] ),
                (  23, [ [204.68000, 141.63600, -246.51300], [-40.28449, -11.37060, -12.15256], [  30.73245, -12.09131,  2.29128], [-14.89368,28.25726,38.04476] ] ),
                (  24, [ [287.46400, 106.72600, -251.65500], [-27.81682,   7.82344,  19.33390], [ -20.50545, -33.64894, 21.60006], [18.86081,40.58653,69.07123] ] ),
                (  25, [ [257.92200, 116.60700, -238.39300], [-30.53741,  11.73331,   6.68288], [  -1.06420, -33.38460,  7.35631], [-3.23290,38.73873,65.30394] ] ),
                (  26, [ [228.11000, 129.47200, -238.39100], [-28.24731,  13.59282,  -6.48616], [  20.41224, -25.44354,  8.38236], [-18.63835,36.71911,60.36230] ] ),
                (  27, [ [228.63372, 235.81606, -159.90821], [-63.49559, - 4.67022,   9.11393], [  -7.67720, -35.93077,-17.25352], [-2.77465,-4.84501,80.63342] ] ),
                (  28, [ [271.47900, 212.59800, -191.07500], [-40.60440, -15.55496,  12.31944], [  32.83245, -30.47300,-38.31859], [ 4.66882,-9.67209,77.93521] ] ),
                (  29, [ [226.88600, 201.94300, -182.15400], [-47.95595, - 5.51536,   5.33273], [   4.39046, -30.83852,-26.76900], [ 0.49085,-6.34474,65.50401] ] ),
                (  30, [ [176.81200, 202.10800, -180.83300], [-51.54777,   5.77320,  -2.65752], [ -27.38298, -44.89569,-30.87234], [-0.54598,-4.15468,66.37461] ] ),
                (  31, [ [291.42800, 178.26800, -232.36000], [-47.54245,   3.79596,  28.32713], [  13.09667, -35.07916,-42.86933], [-0.96373,-14.42396,51.38954] ] ),
                (  32, [ [237.26800, 175.71200, -212.16500], [-59.34096, - .02266,   11.20691], [  12.02454, -26.89842,-32.58142], [-1.75287,-9.65403,39.79137] ] ),
                (  33, [ [175.73200, 160.13200, -211.46600], [-62.30361, -21.64151,  -9.58922], [  16.03710, -32.97626,-36.43165], [ 8.00308, 4.42352,35.78965] ] ),
                (  34, [ [226.47258, 226.91239, - 86.94871], [-70.21971, - 6.56204,   0.09607], [  -1.07781, -34.80919,-20.76585], [-4.86046,-17.18503,59.39277] ] ),
                (  35, [ [276.87000, 199.24100, -113.45500], [-45.99412, -10.81828, -11.43359], [  41.56431, -41.60880,-53.00759], [-8.05186,-15.58240,70.33218] ] ),
                (  36, [ [228.51200, 190.56400, -119.93200], [-50.29281, - 6.43480,  -1.41375], [   5.18484, -36.97720,-44.65763], [-2.95235,-11.14150,60.43018] ] ),
                (  37, [ [177.14300, 186.58700, -116.04100], [-52.01698, - 1.50679,   9.12067], [ -22.15506, -40.80430,-51.77319], [-1.01643,-12.85378,56.83227] ] ),
                (  38, [ [294.61800, 150.47300, -187.48500], [-59.45644,   7.40936,  14.28534], [ -17.23438, -53.48293, 39.36971], [-4.44658,47.14659,70.33578] ] ),
                (  39, [ [237.51300, 155.97400, -176.52600], [-54.63582,   3.57797,   7.60437], [ - 2.25011, -44.18651, 26.88604], [-14.46436,37.15769,59.59554] ] ),
                (  40, [ [185.65700, 157.79900, -171.99400], [-48.93933,   0.07183,   1.45556], [  24.27252, -31.58739, 30.71017], [-24.55990,30.11494,63.74141] ] ),
                (  41, [ [246.49100,  63.88000, -307.11300], [- 6.82060,   6.94095,   5.85056], [ -13.61923,  -4.42194, 13.68915], [20.44407,20.84230,77.79503] ] ),
                (  42, [ [239.32500,  70.69500, -300.57100], [- 7.49968,   6.67712,   7.22339], [  -6.61967,  -9.56751,  2.88190], [22.65399,13.69900,68.94088] ] ),
                (  43, [ [231.51300,  77.17800, -292.65600], [- 8.11343,   6.28045,   8.59508], [   2.10276, -15.03856, -5.78031], [22.50001,14.45866,56.93531] ] ),
                (  44, [ [230.97800,  62.44400, -297.24500], [- 9.80500,   9.97500,   8.60100], [  -9.87027,  -6.79403,  3.69374], [ 6.57943, 3.56483,54.87536] ] ),
                (  45, [ [258.32900,  78.76600, -234.21000], [- 8.86050,   3.99010,   1.14251], [ -35.41421, -21.81730,  5.57717], [ 3.05080, 8.74505,67.32174] ] ),
                (  46, [ [248.79900,  84.34400, -233.10400], [-10.13806,   7.13823,   1.06157], [ -19.14901, -27.52461, -1.11793], [-4.21696,13.29002,64.43814] ] ),
                (  47, [ [238.27500,  93.15700, -232.13600], [-10.87080,  10.45014,   0.87130], [  -3.95002, -35.41035, -1.74361], [-10.04643,16.81149,61.39606] ] ),
                (  48, [ [223.57300,  66.18700, -240.08000], [-17.02700,  20.66500,  -0.83100], [ -28.35511,  -7.96167,-11.62546], [-21.66337, 3.77275,57.16987] ] ),
                (  49, [ [254.22300,  82.22600, -174.23700], [-21.92261,  15.19138,  -0.67785], [ -63.01724, -44.54553, -0.87789], [-11.27104, 4.28773,58.22535] ] ),
                (  50, [ [232.66900,  96.60200, -174.81200], [-21.18048,  13.55722,  -0.47200], [ -30.28502, -49.68869, -6.62427], [-20.10583, 8.25059,53.92831] ] ),
                (  51, [ [211.88800, 109.35800, -175.18600], [-20.37604,  11.95157,  -0.27593], [  -1.09208, -52.71729, -9.55424], [-29.02890, 9.80012,50.83595] ] ),
                (  52, [ [187.82100,  69.71300, -187.14000], [-17.21900,  27.27500,   2.31300], [ -50.52789,  -3.47788,-15.33564], [-30.98173, 6.57516,49.90838] ] ),
                (  53, [ [220.21825, 205.66197, - 42.33880], [-54.99046,   0.02684,  -6.40726], [   0.05371, -19.15208,-11.14849], [-13.50643,-29.31064,32.37773] ] ),
                (  54, [ [258.13000, 182.77700, - 53.57100], [-32.14671, - 2.77922,  -9.15111], [  31.14741, -41.46945,-24.39708], [-41.29860,-13.74292,49.20956] ] ),
                (  55, [ [221.27200, 179.75700, - 61.79100], [-41.50471, - 3.25519,  -7.27050], [   2.05338, -32.51214,-27.67109], [-14.47662,-7.86193,47.81114] ] ),
                (  56, [ [175.16700, 176.30000, - 67.69800], [-50.65118, - 3.65490,  -4.53865], [ -16.89522, -36.08558,-37.61351], [13.47706,-6.68229,49.47788] ] ),
                (  57, [ [270.01700, 129.27200, - 88.09600], [-47.89412,  16.13724, -13.75119], [ -12.34485, -52.58830,-36.48171], [-64.68053,-9.72582,74.23944] ] ),
                (  58, [ [224.62600, 141.72000, - 98.40600], [-42.69048,   8.69225,  -6.81213], [  -6.07947, -40.46491,-33.05834], [-28.11511,-12.34519,60.89072] ] ),
                (  59, [ [185.27400, 147.07700, -102.14500], [-35.75214,   2.00707,  -0.66103], [   3.23605, -32.39972,-32.70151], [ 2.59920,-14.06669,59.62659] ] ),
                (  60, [ [236.41700,  87.16000, -119.82500], [-26.26780,  14.71701,  -6.84809], [ -58.16650, -30.54529,-30.45656], [-47.65768,14.07330,55.97662] ] ),
                (  61, [ [209.60500, 101.12400, -126.12100], [-27.33472,  13.19895,  -5.73831], [ -33.86432, -34.96476,-23.72738], [-35.46433, 6.56261,46.91283] ] ),
                (  62, [ [181.79200, 113.53600, -131.29200], [-28.26820,  11.61556,  -4.59993], [ -12.80399, -36.54972,-20.58378], [-27.06173, 0.31634,48.70083] ] ),
                (  63, [ [161.72100,  78.67100, -141.17900], [-21.72600,  42.02900, -15.24000], [ -57.07212,  -9.16533, -5.88999], [-15.06849,20.16875,54.10475] ] ),
                (  64, [ [203.02800, 174.61900, - 25.28600], [-57.48190,  -3.93560, -15.18085], [ -19.17222, -44.19620, -4.68998], [-21.18458,-2.32339,24.25232] ] ),
                (  65, [ [189.72900, 132.31300, - 57.38600], [-60.34501,  19.97808, -33.24085], [ -22.34843, -33.20915,-31.02907], [-37.99491,-5.89703,19.27990] ] ),
                (  66, [ [162.05800, 109.62300, - 84.65900], [-50.08517,  33.10155, -44.20013], [ -18.24852, -28.95492,-43.61351], [-58.32053,10.20628,35.22056] ] ),
                (  67, [ [ 73.77101, 232.02781, -353.40564], [-61.96117,  -3.47476,  -3.27214], [  13.20375, -37.14290, 36.57091], [ 5.36168,14.44806,101.87448] ] ),
                (  68, [ [128.29817, 193.86368, -339.86060], [-34.19083,  -7.61523,  35.65625], [  33.29531, -50.12321, 17.87942], [ 6.64332, 1.01572,116.83948] ] ),
                (  69, [ [ 80.16891, 198.25628, -320.32901], [-56.82543,  -1.57657,  -2.12030], [  -0.54056, -30.02715, 29.21509], [-0.45140,-6.25962,70.62403] ] ),
                (  70, [ [ 27.23616, 201.07814, -346.07232], [-44.71641,   6.58371, -45.01390], [ -29.41621, -46.29326, 10.61852], [-0.85876,13.89331,79.48348] ] ),
                (  71, [ [154.51637, 136.17472, -271.15466], [-81.81813,  63.44949, -10.60706], [ -33.74301, -18.21468, 11.94521], [-1.99676,16.79498,40.66719] ] ),
                (  72, [ [ 74.44600, 172.92500, -295.84700], [-74.82363,   3.65865, -36.43927], [  -5.91351, -36.69837, 20.49153], [ 3.81517,-6.94447,56.03773] ] ),
                (  73, [ [ 19.59100, 150.81500, -334.75500], [-30.18953, -41.43263, -35.80608], [   3.90168, -55.34070, 11.86376], [-10.54738,31.63248,72.90707] ] ),
                (  74, [ [110.04100, 132.22200, -267.11200], [-37.03161,   5.81778,  -4.18020], [ -32.80840, -30.61837, -1.05827], [17.73640, 4.95587,25.87416] ] ),
                (  75, [ [ 69.15100, 126.20800, -283.82700], [-40.80611, -18.46512, -28.80478], [   5.17211, -46.40825,  6.63929], [-2.57390, 8.97880,43.15763] ] ),
                (  76, [ [ 36.29500,  94.58600, -323.31800], [-23.88850, -42.94968, -48.12749], [  27.25153, -41.37917, 20.59743], [-37.53451,49.79384,66.69272] ] ),
                (  77, [ [102.46800,  96.61600, -272.12400], [-19.87993, -12.23683,  -7.75418], [   8.92165, -24.92582, -2.23952], [ 5.77902,23.22023,35.98234] ] ),
                (  78, [ [ 83.59900,  83.99900, -282.74000], [-17.66199, -12.87648, -13.40134], [  24.80146, -22.88787,  1.66119], [-23.09932,18.37536,19.14698] ] ),
                (  79, [ [ 67.54200,  71.10600, -298.77800], [-14.32693, -12.79780, -18.51304], [  32.26623, -11.30821, 15.96647], [-50.74255,24.27259,13.77249] ] ),
                (  80, [ [ 78.39418, 242.54684, -252.84156], [-64.56190,   2.22180, -16.49452], [   6.93614, -53.89662,  1.26171], [ 3.87961, 6.57643,99.15800] ] ),
                (  81, [ [133.75682, 202.97920, -239.93939], [-50.31990, -15.26410,  -7.57279], [  38.21559, -50.16420, 10.20990], [ 4.73278,-1.82704,82.44790] ] ),
                (  82, [ [ 81.14882, 197.21536, -250.60239], [-53.64818,   4.11497, -13.56541], [  -1.45720, -36.53061,  3.21111], [ 2.41302, 4.20259,68.54914] ] ),
                (  83, [ [ 28.99369, 211.44117, -266.54450], [-49.58059,  23.81715, -17.92775], [ -31.71814, -39.27233, -8.97282], [ 4.37554, 6.80519,79.41450] ] ),
                (  84, [ [152.42300, 152.23200, -233.45500], [-80.34343,  24.38261,  -7.96455], [ -16.42111, -38.93772, -0.32957], [-4.55484,13.45907,36.26878] ] ),
                (  85, [ [ 77.26500, 169.88100, -247.13100], [-69.35399,  10.72766, -19.32613], [  -5.49530, -27.18993,  3.51498], [ 1.81127, 0.87750,41.22460] ] ),
                (  86, [ [ 14.89600, 174.39000, -270.71000], [-54.70185,  -1.68861, -27.48907], [  -7.00255, -35.65020,  1.52682], [ 1.21927,15.33190,54.75512] ] ),
                (  87, [ [128.46100, 142.53200, -236.96600], [-58.18538,   0.70524,   0.71469], [ -21.55380, -16.26582,  1.15824], [18.90423,15.60842,34.12699] ] ),
                (  88, [ [ 70.17900, 142.99000, -243.59200], [-57.76264,   0.20330, -13.97426], [  -4.78554, -27.59853,  2.69204], [ 4.70758,24.31423,36.00992] ] ),
                (  89, [ [ 14.15300, 142.92900, -264.77700], [-53.75475,  -0.32209, -28.11612], [   4.35738, -30.96406,  8.11414], [-6.08703,46.01332,49.21219] ] ),
                (  90, [ [112.85800, 122.96800, -231.66900], [-45.05147, -10.89848,  -8.59796], [ -10.64309, -27.65369,  1.69910], [14.98001,29.39949,44.79706] ] ),
                (  91, [ [ 67.80400, 115.23000, -241.80600], [-44.88763,  -4.53666, -11.64381], [   7.96662, -30.41510,  4.31772], [-6.02814,42.12770,60.67983] ] ),
                (  92, [ [ 23.42400, 113.89200, -254.85000], [-43.70265,   1.85346, -14.38831], [  21.73880, -35.91609, 12.13603], [-29.81552,57.62669,71.99957] ] ),
                (  93, [ [ 81.54483, 245.30030, -155.33935], [-70.83114,   6.11577, -11.63960], [   8.52360, -41.09975,-27.17989], [ 2.43537,-5.74056,80.18700] ] ),
                (  94, [ [137.51453, 200.35571, -174.84086], [-54.19592,   4.96162, -11.67151], [  31.24860, -46.20180,-26.76120], [ 4.18297,-13.17289,62.43592] ] ),
                (  95, [ [ 84.92800, 206.45900, -183.82300], [-50.88567,   7.23658,  -6.27308], [  -1.83217, -36.22166,-29.54855], [ 6.29304,-0.56475,67.58513] ] ),
                (  96, [ [ 35.94235, 214.70120, -187.54303], [-46.98626,   9.22830,  -1.16451], [ -31.26538, -33.47569,-38.69700], [ 6.58466, 2.16199,76.79664] ] ),
                (  97, [ [146.83400, 163.96300, -199.82700], [-73.67951,   9.29057, -14.41012], [  -4.82075, -31.73641,-33.95790], [-12.39683,-0.54149,28.06956] ] ),
                (  98, [ [ 78.26400, 173.39400, -213.78500], [-63.45119,   9.57025, -13.50406], [  -7.40980, -31.80247,-29.97100], [ 1.45578,-1.29769,43.62631] ] ),
                (  99, [ [ 19.91500, 182.97600, -226.74800], [-53.23611,   9.59182, -12.41944], [ -11.16675, -36.28244,-39.16367], [ 4.61389, 4.05126,58.15681] ] ),
                ( 100, [ [ 83.50224, 234.64926, - 93.34096], [-58.79918, -13.27188,  -2.98609], [  11.68399, -31.50287, -9.03831], [ 9.58529,-18.59117,60.22329] ] ),
                ( 101, [ [145.80563, 174.07081, -118.02793], [-46.74854,  17.81685,  -2.30639], [ -13.34685, -78.32083,-40.55773], [ 8.54408,-16.04052,56.07402] ] ),
                ( 102, [ [ 93.58800, 196.02800, -117.20400], [-49.95875,  19.16507,   1.86099], [  -0.43197, -37.89186,-36.52477], [ 8.86393,-14.36488,62.11610] ] ),
                ( 103, [ [ 41.92640, 212.44694, -112.86108], [-52.85215,  20.39957,   6.28767], [ -18.68472, -20.86182,-37.31912], [15.09654,-6.79844,67.94684] ] ),
                ( 104, [ [133.05500, 155.14700, -183.01400], [-49.38179,   7.00844,  26.28397], [ -18.24900, -30.11614, 18.53991], [15.50630,28.60785,58.61571] ] ),
                ( 105, [ [ 80.12000, 164.70100, -163.61300], [-55.81089,  12.00344,  12.15752], [   5.48197, -51.45307, 35.70974], [13.86825,40.73158,62.33716] ] ),
                ( 106, [ [ 22.83900, 178.91500, -159.27800], [-58.03270,  16.22372,  -3.44488], [  16.80473, -64.31852, 51.09538], [13.86011,52.54326,71.32694] ] ),
                ( 107, [ [113.14500,  88.04600, -272.21500], [- 8.43016,  -8.66642,  -9.02426], [  12.25223,  -9.54323,  0.97035], [-6.92901,-1.83778,25.00206] ] ),
                ( 108, [ [104.28100,  79.01900, -281.41300], [- 9.29690,  -9.38662,  -9.37074], [  22.32816,  -3.11395,  6.37341], [-23.63112, 7.26560,38.06333] ] ),
                ( 109, [ [ 94.54500,  69.27200, -290.94400], [-10.17418, -10.10647,  -9.69039], [  30.55827,   2.92727, 14.09682], [-40.00238, 2.49956,37.65681] ] ),
                ( 110, [ [126.00500,  78.32100, -270.12000], [- 1.52300,  -8.67600,  -5.12400], [  20.52272,   1.66938, 15.75421], [13.02504,-7.49930,24.26177] ] ),
                ( 111, [ [110.77900,  90.71500, -234.92900], [-20.82535,   0.51561,   2.85764], [  12.90491, -29.42784, -6.21718], [ 2.43409, 7.23866,48.71450] ] ),
                ( 112, [ [ 87.11100,  87.49900, -235.17900], [-25.90332,  -6.96265,  -3.44097], [  34.70113, -25.23976,  1.19083], [-10.12621, 9.51526,53.46615] ] ),
                ( 113, [ [ 60.27700,  76.35300, -242.36200], [-27.36242, -15.10726, -10.76674], [  59.10038, -26.02342,  8.21833], [-27.36205,11.58923,58.40431] ] ),
                ( 114, [ [134.24700,  72.87000, -242.12600], [- 0.24000, -10.72900,   1.39900], [  56.38595,  -3.80341,-14.27833], [ 3.04781,-3.16598,30.96038] ] ),
                ( 115, [ [120.54600, 104.05100, -176.46100], [-30.22927,  -6.78939,   0.39580], [  -0.52830, -49.98930,-17.24829], [17.65083, 7.33038,47.96876] ] ),
                ( 116, [ [ 86.01400,  97.68300, -176.26300], [-38.80687,  -5.94035,  -0.00016], [  31.73477, -53.29160,-28.21754], [11.83818, 5.54348,52.97363] ] ),
                ( 117, [ [ 42.94400,  92.42700, -176.50800], [-47.31026,  -4.56944,  -0.48960], [  65.57180, -66.08095,-30.85236], [ 6.26843, 8.95760,67.95456] ] ),
                ( 118, [ [131.42000,  72.41200, -210.02400], [- 0.22600, -19.16800,   1.86300], [  50.24854,   2.33869,-33.43067], [10.88159,-0.41785,46.58427] ] ),
                ( 119, [ [100.41654, 208.82619, - 37.66502], [-57.67642,  14.49144,  -1.30042], [  10.40468, -25.81947,-16.62988], [21.99927,-34.77575,35.39534] ] ),
                ( 120, [ [153.82600, 168.08500, - 64.52900], [-56.06588,   2.38347,   3.63052], [   8.70418, -48.67971,-45.75992], [-11.85076,-1.56973,54.58364] ] ),
                ( 121, [ [102.49800, 178.38300, - 60.04200], [-45.57797,  18.16950,   5.27794], [  -6.74504, -33.81801,-27.31969], [15.64260,-9.75675,47.66143] ] ),
                ( 122, [ [ 64.18900, 202.49300, - 54.33900], [-30.20692,  29.24394,   5.96358], [ -37.74103, -27.04993,-24.31724], [42.26290,-21.94720,50.07548] ] ),
                ( 123, [ [123.54000, 138.16800, -101.48500], [-34.04147,   3.25521,  11.07149], [  -5.75333, -35.47229,-43.69763], [-8.48470,-11.05867,65.87765] ] ),
                ( 124, [ [ 85.02100, 144.01300, - 90.69500], [-42.90106,   8.42566,  10.47746], [  -0.49183, -43.01317,-39.19981], [23.93840,-17.02025,56.10322] ] ),
                ( 125, [ [ 37.90900, 155.36900, - 80.88000], [-51.24508,  14.26467,   9.13865], [  -0.99332, -62.08434,-35.75497], [29.37213,-35.88090,67.46481] ] ),
                ( 126, [ [141.97200, 106.50800, -140.39200], [-33.02569,  -8.26429,   4.12525], [  21.54894, -35.18838,-28.54772], [15.82392,-0.81452,47.47247] ] ),
                ( 127, [ [107.36900,  98.97900, -132.63300], [-35.96992,  -6.74106,  11.36647], [  42.06196, -38.61797,-35.69944], [30.05861, 1.85532,44.88732] ] ),
                ( 128, [ [ 70.49800,  93.22900, -117.45400], [-37.59767,  -4.73696,  18.90384], [  66.64553, -50.02738,-42.81734], [55.68621, 4.60237,52.62383] ] ),
                ( 129, [ [163.74900,  72.26500, -157.01400], [- 9.56000, -38.17200,  14.94100], [  66.25687, -13.87968,-12.24198], [11.42222,14.68756,67.86727] ] ),
                ( 130, [ [121.25700, 174.97900, - 23.48600], [-59.48007,  21.83152,  13.42335], [  12.10669, -47.30221, -8.25033], [20.84002, 2.80918,24.24598] ] ),
                ( 131, [ [117.05700, 132.47900, - 55.76400], [-59.63637,   8.29556,  37.61918], [  13.62644, -39.62645,-35.02889], [34.45597,-5.19218,11.81235] ] ),
                ( 132, [ [146.43500, 101.41700, - 88.48000], [-57.99731, -13.25263,  40.93618], [  26.27834, -32.59566,-51.11229], [47.63012, 2.99282,43.01833] ] )
                ]

            # Create nodes
            nodeIndex = 0
            nodeIdentifier = 1
            lowerLeftNodeIds = []
            upperLeftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []

            # Left lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(leftLung, cache, coordinates, generateParameters,
                nodes, nodetemplate, nodeFieldParameters,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerLeftNodeIds, upperLeftNodeIds, nodeIndex, nodeIdentifier)

            # Right lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(rightLung, cache, coordinates, generateParameters,
                nodes, nodetemplate, nodeFieldParameters,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, nodeIndex, nodeIdentifier)

            # Create elements
            elementIdentifier = 1

            # Left lung elements
            elementIdentifier = getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                leftLungMeshGroup, lowerLeftLungMeshGroup, None, upperLeftLungMeshGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerLeftNodeIds, upperLeftNodeIds, elementIdentifier)

            # Right lung elements
            getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                rightLungMeshGroup, lowerRightLungMeshGroup, middleRightLungMeshGroup, upperRightLungMeshGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, elementIdentifier)

            # Marker points
            lungNodesetGroup = lungGroup.getNodesetGroup(nodes)
            markerList = []

            lowerLeftElementCount = (lElementsCount1 * (lElementsCount2-1) * lElementsCount3 + lElementsCount1)

            idx = lowerLeftElementCount + (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)
            markerList.append({ "group" : leftApexGroup, "elementId" : idx, "xi" : [0.0, 0.0, 1.0] })

            idx = 1
            markerList.append({"group": leftDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 // 2)
            markerList.append({"group": leftMedialGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 - 1) + 1
            markerList.append({"group": leftVentralGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            upperLeftElementCount = (uElementsCount1 * uElementsCount2 * (uElementsCount3-1))
            leftLungElementCount = lowerLeftElementCount + upperLeftElementCount
            lowerRightElementCount = lowerLeftElementCount

            idx = leftLungElementCount + 1
            markerList.append({"group": rightDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = leftLungElementCount + lowerRightElementCount + \
                  (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)
            markerList.append({"group": rightApexGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx = leftLungElementCount + lElementsCount1 + 1
            markerList.append({"group": rightMedialGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            idx = leftLungElementCount + lElementsCount1 * lElementsCount2
            markerList.append({"group": rightVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * (lElementsCount3 - 1)) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 1.0, 1.0]})

        elif isMouse:
            # valueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ]
            nodeFieldParameters = [
                (   1, [ [   0.07403,  -9.56577,   3.60357], [   0.95906,   0.71429,   2.93871], [  -1.67031,  -0.43146,   0.60812], [  -1.76769,  -4.15710,   3.30552] ] ),
                (   2, [ [   0.00308,  -9.73991,   6.23603], [  -0.91221,  -0.88041,   1.92742], [  -1.86919,   0.54982,   0.87706], [  -0.58548,  -0.52947,   2.21052] ] ),
                (   3, [ [  -1.92122, -11.21142,   1.24012], [   0.20451,   1.97277,   3.02946], [  -2.80724,  -0.35177,  -1.27853], [  -2.64529,  -2.41687,   4.96326] ] ),
                (   4, [ [  -1.78634,  -9.65609,   4.22937], [   0.06380,   1.12389,   2.92754], [  -2.00998,   0.26127,   0.62875], [  -1.98142,  -2.36558,   1.94931] ] ),
                (   5, [ [  -1.78517,  -8.93248,   7.02409], [  -0.06077,   0.31970,   2.63203], [  -1.68448,   1.05833,   0.68834], [  -0.39630,   0.49892,   2.33932] ] ),
                (   6, [ [  -4.29865, -10.50146,   1.20789], [   0.37882,   1.55423,   4.04373], [  -2.63582,   1.48617,   0.80382], [  -2.39199,  -1.12119,   5.24368] ] ),
                (   7, [ [  -3.85980,  -9.00517,   4.82763], [   0.49854,   1.43694,   3.19209], [  -2.06718,   1.18962,   0.78973], [  -1.21391,  -1.24641,   1.50229] ] ),
                (   8, [ [  -3.32678,  -7.65254,   7.59759], [   0.56646,   1.26599,   2.34353], [  -1.57620,   1.59817,   0.72268], [   0.28745,   0.97206,   2.28447] ] ),
                (   9, [ [  -6.57368,  -8.33292,   2.96552], [   0.78921,   0.84268,   2.86361], [  -1.63553,   2.60797,   2.67981], [  -0.81873,  -0.02233,   5.10575] ] ),
                (  10, [ [  -5.75096,  -7.25801,   5.78564], [   0.85324,   1.30395,   2.76578], [  -1.66434,   1.74323,   0.77869], [  -0.75552,  -0.82191,   1.46800] ] ),
                (  11, [ [  -4.87296,  -5.72876,   8.47415], [   0.89962,   1.74845,   2.60216], [  -1.18379,   2.18438,   0.84461], [   0.85546,   0.88778,   2.03581] ] ),
                (  12, [ [  -7.18669,  -5.55843,   6.39603], [   0.55788,   2.60065,   3.29819], [  -1.20131,   1.64796,   0.43996], [   0.52133,   0.38157,   2.81760] ] ),
                (  13, [ [  -5.65157,  -3.36472,   9.25580], [   2.37223,   1.68711,   2.28630], [  -0.36630,   2.49508,   0.70495], [   1.25584,   1.27257,   2.59212] ] ),
                (  14, [ [  -1.46058, -12.47910,   7.13881], [   2.13302,   1.46092,   1.58836], [  -1.88471,   1.14866,  -0.03002], [  -1.24286,  -1.53158,   3.65525] ] ),
                (  15, [ [  -0.43063, -10.07412,   8.71137], [  -0.06459,   2.95816,   1.37507], [  -1.52906,   1.62194,   0.84465], [  -0.27598,  -0.13356,   2.71766] ] ),
                (  16, [ [  -4.02723, -12.66079,   5.97470], [   0.49565,   0.81981,   0.68090], [  -2.55416,   0.85551,  -0.61402], [  -1.51986,  -0.43904,   4.41795] ] ),
                (  17, [ [  -3.22600, -11.20300,   7.18000], [   1.10652,   2.09529,   1.72930], [  -1.63687,   1.39790,   0.11255], [  -0.66214,  -0.44676,   3.72000] ] ),
                (  18, [ [  -1.83711,  -8.46276,   9.43760], [   1.67115,   3.38496,   2.78572], [  -1.28118,   1.59789,   0.60630], [   0.29795,   0.43356,   2.45509] ] ),
                (  19, [ [  -6.00294, -11.06987,   5.98787], [   0.97269,   0.66609,   0.76324], [  -1.71493,   2.34334,   0.60678], [  -0.99047,  -0.00338,   4.25902] ] ),
                (  20, [ [  -4.71987,  -9.70115,   7.36101], [   1.57616,   2.05951,   1.96948], [  -1.50108,   1.76558,   0.35819], [  -0.33627,   0.02896,   3.35413] ] ),
                (  21, [ [  -2.99367,  -6.89314,   9.93095], [   1.87184,   3.54817,   3.16296], [  -1.13503,   1.70437,   0.56225], [   0.37722,   0.54150,   2.36993] ] ),
                (  22, [ [  -7.05528,  -8.21578,   7.23623], [  -0.04743,   0.21677,   0.14801], [  -0.35683,   3.03360,   1.61413], [  -0.14227,   0.25667,   3.42193] ] ),
                (  23, [ [  -6.17292,  -7.68518,   7.91870], [   1.67759,   1.45942,   1.63685], [  -1.01458,   2.30610,   0.89443], [  -0.01928,   0.04277,   2.66382] ] ),
                (  24, [ [  -4.09024,  -5.06093,  10.56252], [   2.47439,   3.76869,   3.63115], [  -0.77293,   2.14754,   0.90306], [   0.70426,   0.44195,   2.12733] ] ),
                (  25, [ [  -6.63445,  -5.23571,   9.11277], [   1.29479,   2.91361,   2.31569], [   0.08865,   2.51141,   1.44679], [   0.58295,   0.26372,   2.61475] ] ),
                (  26, [ [  -4.41971,  -2.69228,  11.72475], [   3.04409,   2.11044,   2.82422], [   0.11158,   2.53492,   1.39130], [   1.17715,   0.04117,   2.28236] ] ),
                (  27, [ [  -2.33590, -12.83832,  10.40150], [   2.88381,   1.20100,   1.09925], [  -0.18173,   2.66677,   0.70385], [  -0.87332,   0.16147,   3.18425] ] ),
                (  28, [ [  -0.51618,  -9.97224,  11.62323], [   0.65861,   3.94935,   1.17161], [  -0.27951,   1.97737,   0.12068], [  -0.79293,  -0.21003,   2.63273] ] ),
                (  29, [ [  -4.99078, -12.32037,   9.79621], [   2.17352,   1.67823,   1.22241], [  -2.16881,   1.27754,  -0.44518], [  -0.48233,   0.66527,   3.61343] ] ),
                (  30, [ [  -2.94100, -10.41126,  10.92996], [   1.91384,   2.13057,   1.03823], [  -1.02412,   2.12341,   0.33619], [  -0.11625,   0.59407,   3.35135] ] ),
                (  31, [ [  -1.18625,  -8.08002,  11.85981], [   1.58710,   2.51833,   0.81707], [  -1.05238,   1.74869,   0.34892], [   0.05191,  -0.11824,   2.50982] ] ),
                (  32, [ [  -6.41236, -10.67331,   9.55112], [   2.38573,   1.91452,   1.67719], [  -1.12592,   2.29408,  -0.04315], [  -0.05981,   0.67603,   3.43748] ] ),
                (  33, [ [  -4.25751,  -8.67962,  11.08139], [   1.91673,   2.06705,   1.37826], [  -1.36250,   1.71786,  -0.03767], [   0.25731,   0.69329,   3.22903] ] ),
                (  34, [ [  -2.57791,  -6.57308,  12.30631], [   1.43636,   2.13694,   1.06705], [  -1.15376,   1.64565,   0.42215], [   0.22890,  -0.09508,   2.33575] ] ),
                (  35, [ [  -6.97952,  -7.90413,   9.77554], [   0.86359,   0.30450,   0.66060], [   0.11900,   3.06937,   1.03036], [   0.39443,   0.52278,   2.69307] ] ),
                (  36, [ [  -5.63978,  -7.00879,  10.85483], [   1.79715,   1.47957,   1.48364], [  -0.94905,   1.95983,   0.28992], [   0.69943,   0.67916,   2.82106] ] ),
                (  37, [ [  -3.47246,  -4.84263,  12.69322], [   2.53100,   2.84545,   2.18753], [  -0.42103,   1.80597,   0.72699], [   0.55941,  -0.14862,   2.08769] ] ),
                (  38, [ [  -6.02540,  -5.02789,  11.62196], [   1.85963,   2.51542,   2.04771], [   0.16536,   1.86189,   1.15728], [   0.77715,   0.08202,   2.22181] ] ),
                (  39, [ [  -3.37913,  -3.15743,  13.67816], [   3.29322,   1.17563,   1.98068], [   0.57374,   1.47703,   1.17346], [   0.86448,  -0.80650,   1.86683] ] ),
                (  40, [ [  -3.17111, -12.20100,  13.38085], [   2.12652,   0.55899,   0.44433], [   0.16750,   2.33744,   0.63169], [  -0.85030,   1.34507,   2.74098] ] ),
                (  41, [ [  -1.85859, -10.43779,  13.75973], [   1.44952,   1.72894,   0.60352], [   0.24785,   1.72556,   0.51281], [  -1.01790,  -0.02591,   2.50059] ] ),
                (  42, [ [  -5.03857, -11.40132,  13.08586], [   1.70052,   1.42937,   0.86783], [  -1.58299,   1.32180,  -0.32245], [   0.56034,   1.07102,   3.01045] ] ),
                (  43, [ [  -3.35343,  -9.99978,  13.79183], [   1.66537,   1.37001,   0.54186], [  -0.53576,   2.01439,   0.17660], [  -0.56268,   0.01490,   2.40321] ] ),
                (  44, [ [  -1.71989,  -8.67015,  14.17671], [   1.59657,   1.28512,   0.22717], [  -0.35050,   1.77131,   0.38853], [  -0.62823,  -0.95180,   2.32388] ] ),
                (  45, [ [  -6.14401,  -9.75429,  12.77928], [   2.16121,   1.63557,   1.09963], [  -0.71076,   2.21275,  -0.29315], [   0.92132,   1.20642,   3.09774] ] ),
                (  46, [ [  -4.15956,  -8.25973,  13.75872], [   1.80752,   1.35343,   0.85917], [  -0.72593,   1.83964,  -0.12732], [  -0.23415,   0.14089,   2.27387] ] ),
                (  47, [ [  -2.53185,  -7.04795,  14.50314], [   1.44769,   1.06998,   0.62958], [  -0.63860,   1.69106,   0.24785], [   0.29451,  -0.80052,   2.67644] ] ),
                (  48, [ [  -6.22965,  -7.15528,  12.54983], [   1.25476,   0.72701,   0.90323], [   0.54615,   2.54332,   0.40531], [   1.05654,   0.75533,   2.60184] ] ),
                (  49, [ [  -4.79229,  -6.33533,  13.53300], [   1.61960,   0.91268,   1.06285], [  -0.48367,   1.61129,  -0.10092], [   0.24131,   0.20255,   2.14045] ] ),
                (  50, [ [  -2.98701,  -5.33022,  14.66882], [   1.99066,   1.09737,   1.20860], [  -0.04573,   1.48484,   0.52507], [   0.62367,  -0.56781,   2.05117] ] ),
                (  51, [ [  -5.14160,  -5.05188,  13.52129], [   1.85459,   1.66867,   1.50671], [  -0.21383,   0.95066,   0.07710], [   0.60599,  -0.78027,   1.80907] ] ),
                (  52, [ [  -2.71728,  -4.25018,  15.38349], [   2.81886,  -0.06146,   2.08793], [   0.52113,   0.60132,   0.80528], [   0.59496,  -1.42790,   1.58518] ] ),
                (  53, [ [  -3.97000, -10.25800,  15.67100], [   1.60402,   0.70083,   1.68144], [  -0.82192,   2.28449,   0.98993], [  -0.64717,  -0.51289,   1.30805] ] ),
                (  54, [ [  -2.42877,  -9.98303,  16.42373], [   1.15748,  -0.11813,  -0.13777], [   0.04611,   1.36745,   2.12860], [  -0.78095,  -1.65578,   2.14658] ] ),
                (  55, [ [  -4.60043,  -8.32257,  15.55367], [   2.24339,   0.78867,   2.52215], [  -0.47795,   1.97306,  -0.38370], [  -0.62324,  -0.25655,   1.26655] ] ),
                (  56, [ [  -1.90227,  -8.23856,  17.57678], [   2.98703,  -0.58799,   1.44387], [   0.14674,   2.33566,   0.26025], [   0.95946,  -1.57221,   3.45219] ] ),
                (  57, [ [  -4.90484,  -6.41106,  14.92762], [   2.03469,   0.62463,   2.22346], [  -0.28161,   1.70116,  -1.06322], [  -0.37351,  -0.28350,   0.51957] ] ),
                (  58, [ [  -2.21713,  -5.98027,  16.77811], [   3.24505,   0.23016,   1.43520], [  -0.42259,   2.04103,  -1.13906], [   0.91436,  -0.73090,   2.16332] ] ),
                (  59, [ [   3.21000, -14.40000,   2.61000], [   2.45000,  -0.85500,   0.01500], [   0.34400,   2.19000,   1.23000], [   0.00400,  -0.80200,   2.52000] ] ),
                (  60, [ [   0.88200, -11.40000,   2.59000], [   2.06000,  -0.29900,   1.50000], [  -1.23000,   2.97000,   0.75300], [  -0.22000,  -1.37000,   1.92000] ] ),
                (  61, [ [   3.47000, -12.00000,   3.89000], [   3.07000,  -1.43000,   0.64300], [   0.57500,   1.72000,   2.00000], [   0.02000,  -1.12000,   1.70000] ] ),
                (  62, [ [   6.56000, -14.00000,   3.58000], [   2.97000,  -2.53000,  -1.24000], [   3.13000,   1.39000,   1.62000], [  -0.06700,  -0.64000,   2.48000] ] ),
                (  63, [ [   0.77700, -10.50000,   6.79000], [   3.07000,  -0.42700,  -0.51900], [   1.10000,   0.84100,   2.57000], [  -1.04000,  -1.14000,   1.95000] ] ),
                (  64, [ [   4.14000, -11.00000,   6.27000], [   4.18000,  -0.53400,  -0.58200], [   0.53400,   1.27000,   2.35000], [  -0.11100,  -0.23400,   1.97000] ] ),
                (  65, [ [   8.84000, -11.70000,   5.54000], [   5.22000,  -0.85200,  -0.88100], [   2.29000,   3.24000,   2.57000], [  -0.40200,  -0.68200,   1.53000] ] ),
                (  66, [ [   1.69000,  -9.41000,   9.46000], [   1.74000,   0.27200,  -0.64900], [  -0.32900,   0.88800,   2.23000], [  -0.55100,  -1.14000,   0.63600] ] ),
                (  67, [ [   4.49000,  -9.40000,   8.57000], [   3.82000,   0.49200,  -0.79600], [   0.08200,   2.48000,   2.40000], [  -0.02500,  -0.13800,   0.93500] ] ),
                (  68, [ [   9.35000,  -8.58000,   8.03000], [   3.28000,   1.25700,   0.50500], [  -0.55000,   2.58000,   1.57000], [  -0.33400,  -1.54000,   0.45000] ] ),
                (  69, [ [   1.02000,  -8.07000,  11.90000], [   2.34000,   2.47000,  -1.51000], [  -0.99100,   3.00000,   1.34000], [   0.18900,  -0.82600,   1.25000] ] ),
                (  70, [ [   4.18000,  -6.03000,  10.70000], [   3.39000,   1.34000,  -1.05000], [  -1.25000,   3.11000,   1.17000], [   0.22400,  -1.86000,   0.43500] ] ),
                (  71, [ [   7.92000,  -5.44000,   9.89000], [   3.64500,   0.55700,  -0.59800], [  -0.83300,   1.42000,   0.60700], [   0.64600,  -1.31000,   0.80000] ] ),
                (  72, [ [   3.31000, -15.10000,   5.27000], [   2.59000,  -0.53600,  -0.02700], [   0.15400,   2.25000,   0.50800], [   0.26500,  -0.37600,   3.06000] ] ),
                (  73, [ [   0.91900, -13.40000,   6.67000], [   2.53100,   0.43500,   0.04200], [  -1.55000,   1.75000,   2.17000], [   0.17000,  -0.55500,   3.21000] ] ),
                (  74, [ [   3.66000, -13.07900,   6.59200], [   2.92900,  -0.58200,  -0.17400], [   0.44000,   1.98000,   1.54000], [   0.28000,  -0.79900,   2.90000] ] ),
                (  75, [ [   6.47000, -14.30000,   6.16000], [   2.66000,  -2.23000,  -0.02900], [   1.75000,   1.43000,   0.96800], [   0.24600,   0.20600,   2.69000] ] ),
                (  76, [ [  -0.19700, -11.60000,   8.90000], [   4.38000,   0.93300,  -0.40500], [   1.03000,   0.73100,   1.20000], [  -0.25400,  -0.21500,   2.42000] ] ),
                (  77, [ [   4.19400, -11.28500,   8.31500], [   4.41300,  -0.56100,  -0.89600], [   0.43300,   1.59000,   1.78000], [   0.15900,  -0.27800,   2.03000] ] ),
                (  78, [ [   8.29000, -12.30000,   7.11000], [   2.56500,  -1.43000,  -0.76300], [   1.03000,   1.97000,   0.91300], [  -0.40400,  -0.27900,   1.46000] ] ),
                (  79, [ [   0.90300, -10.80000,  10.40000], [   2.72000,   1.16000,  -0.46300], [   0.46700,   1.34000,   2.46000], [  -0.94100,  -1.12300,   0.97300] ] ),
                (  80, [ [   4.53300,  -9.71400,   9.97600], [   4.48800,   0.06600,  -1.13000], [   0.19000,   1.69000,   1.54000], [  -0.01500,  -1.18000,   1.14000] ] ),
                (  81, [ [   9.03000, -10.10000,   8.41000], [   2.04200,  -1.07500,  -1.05600], [   0.08600,   1.44000,   1.87000], [  -0.52600,  -1.73000,   0.70400] ] ),
                (  82, [ [   1.40200,  -9.18300,  13.23100], [   2.30100,   1.11800,  -1.69000], [   0.89900,   1.61000,   1.22000], [  -0.21600,  -1.29000,   1.83000] ] ),
                (  83, [ [   4.54000,  -7.87000,  11.30000], [   3.58100,   1.06000,  -0.76100], [  -0.31900,   2.31000,   1.22900], [   0.48000,  -1.81000,   0.71500] ] ),
                (  84, [ [   8.33000,  -6.82000,  10.60000], [   4.01100,   1.10400,  -0.40400], [  -0.97400,   2.10000,   1.31000], [   0.30800,  -2.36000,   1.31000] ] ),
                (  85, [ [   3.67000, -15.20000,   8.43000], [   3.66000,  -0.19900,   0.11900], [   0.44400,   1.71000,   1.59000], [   0.16900,   0.35900,   2.69000] ] ),
                (  86, [ [   1.24000, -13.70000,  10.10000], [   2.30000,   0.53400,   0.14700], [  -1.13000,   0.95700,   1.01000], [   0.38000,  -0.15800,   3.00000] ] ),
                (  87, [ [   4.04000, -13.40000,   9.71000], [   2.78000,   0.10800,  -0.90100], [   0.44400,   1.85000,   0.91400], [   0.63100,   0.42500,   3.05000] ] ),
                (  88, [ [   6.82000, -13.90000,   8.87000], [   2.89000,  -0.77800,  -0.24100], [   2.12000,   1.57000,   0.27200], [   0.38800,   0.64900,   2.08000] ] ),
                (  89, [ [  -0.21500, -11.89600,  11.60400], [   3.40000,   0.10000,  -1.14600], [   1.05800,   1.23800,  -1.21900], [   0.21900,  -0.02100,   1.92900] ] ),
                (  90, [ [   4.46000, -11.50000,  10.30000], [   4.10000,   0.28400,  -1.23000], [   0.23100,   1.93000,   0.05800], [   0.34200,   0.08800,   1.66900] ] ),
                (  91, [ [   8.43000, -11.60000,   9.04000], [   2.10300,  -0.39700,  -0.46200], [   0.83200,   1.39000,  -0.62600], [   0.42400,   1.62000,   2.11000] ] ),
                (  92, [ [   4.07100, -14.59800,  11.20400], [   2.31200,   0.61700,  -0.27000], [   0.65900,   1.60300,   0.60700], [   0.40800,   0.67700,   3.13000] ] ),
                (  93, [ [   1.60000, -13.50000,  13.30000], [   3.15200,   0.71300,  -1.53500], [  -0.98900,   1.18900,   1.04700], [   0.67100,   0.45200,   2.11000] ] ),
                (  94, [ [   4.61000, -12.80000,  11.90000], [   2.43900,   0.41200,  -0.57300], [   0.51900,   2.15000,   0.56800], [  -0.40500,  -0.00900,   1.83700] ] ),
                (  95, [ [   7.18000, -13.00000,  11.10000], [   1.24900,  -0.68100,  -0.50900], [   1.62000,   2.25000,   0.54700], [  -0.20100,   0.76100,   2.90000] ] ),
                (  96, [ [   0.19500, -11.60000,  14.40000], [   4.67000,   1.97000,  -2.39000], [   0.52800,   0.71700,   1.76000], [  -0.16000,  -0.52100,   0.39100] ] ),
                (  97, [ [   4.95300,  -9.85200,  12.34500], [   4.34200,   0.68900,  -0.96600], [   0.31300,   1.53000,   2.23000], [   0.19600,  -2.32000,   0.49900] ] ),
                (  98, [ [   8.57700,  -9.27300,  11.71000], [   1.67600,   0.01000,  -0.41500], [  -0.58900,   1.41000,   2.69000], [  -0.52300,  -3.64000,   0.34200] ] ),
                (  99, [ [   0.17500,  -4.79000,  11.40000], [   2.19000,   1.07000,  -0.47700], [  -1.17000,   3.24000,  -0.42000], [   2.95000,  -1.42000,   2.27000] ] ),
                ( 100, [ [   2.42000,  -3.62000,  11.00000], [   2.29000,   1.28000,  -0.45800], [  -2.62000,   2.34000,   0.24600], [   1.85000,  -0.53100,   1.67000] ] ),
                ( 101, [ [   4.36000,  -1.84000,  10.70000], [   2.08000,   1.65000,  -0.35600], [  -4.71000,   2.20000,   0.65100], [   1.87000,  -0.62300,   1.67000] ] ),
                ( 102, [ [  -1.02000,  -1.69000,  11.00000], [   0.77400,   1.81000,  -0.26700], [  -1.41000,   0.43000,   0.19700], [   0.91000,   1.90000,   2.73000] ] ),
                ( 103, [ [   2.50200,  -6.10300,  13.80300], [   1.32800,   1.43600,  -0.84000], [   0.26600,   3.66000,   0.88300], [   1.03000,  -1.59000,   2.47000] ] ),
                ( 104, [ [   3.88300,  -4.54900,  12.89500], [   1.71000,   1.65000,  -0.78500], [  -1.79000,   3.19000,   1.46000], [   0.70800,  -1.36000,   1.72300] ] ),
                ( 105, [ [   6.06000,  -2.90000,  12.40000], [   1.02600,   0.45300,  -0.08800], [  -3.10000,   2.21000,   1.63000], [   0.95000,  -1.10000,   1.62000] ] ),
                ( 106, [ [   1.07000,  -1.86000,  14.10000], [   0.56000,   2.90000,   0.59200], [  -3.00000,   0.21000,   1.58000], [   0.70100,  -0.98000,   1.70900] ] ),
                ( 107, [ [   2.37000,  -7.49000,  15.90000], [   1.75700,   1.17000,  -0.77700], [   0.63200,   3.39800,  -0.01300], [  -0.58000,  -1.14000,   2.10000] ] ),
                ( 108, [ [   4.35000,  -6.18000,  15.00000], [   2.03000,   1.57000,  -1.03000], [  -1.98000,   3.34000,   2.01000], [  -0.07100,  -1.66000,   2.17000] ] ),
                ( 109, [ [   6.58000,  -4.52000,  13.90000], [   0.63600,   1.38500,  -0.17500], [  -1.88000,   1.91000,   1.25000], [  -0.63500,  -1.92000,   2.68000] ] ),
                ( 110, [ [   1.86000,  -3.31000,  16.00000], [   0.37500,   2.65100,   1.04300], [  -2.62000,   0.83400,   0.85500], [  -0.15700,  -1.09000,   1.66000] ] ),
                ( 111, [ [   4.42000, -13.70000,  14.20000], [   2.33000,   0.71600,   0.00700], [  -0.22400,   0.74100,   0.29300], [   0.14000,   1.05000,   1.94000] ] ),
                ( 112, [ [   2.64000, -12.70000,  15.30000], [   0.94900,  -0.20000,  -0.35200], [  -1.44000,   1.35000,   1.27000], [   1.14000,   0.58600,   1.28000] ] ),
                ( 113, [ [   4.27000, -12.70000,  14.60000], [   2.09400,   0.41000,  -0.44700], [  -0.31000,   1.48000,   0.66300], [   0.04800,   0.24200,   1.17500] ] ),
                ( 114, [ [   6.77000, -12.20000,  14.00000], [   1.14000,   0.37800,   0.29900], [   1.21000,   1.53000,   0.27000], [  -1.19000,   0.02700,   1.81000] ] ),
                ( 115, [ [   1.07000, -10.90000,  16.80000], [   2.55200,   0.17900,  -1.11500], [  -0.66800,   2.40000,   1.38000], [   1.79000,   0.15800,   1.91000] ] ),
                ( 116, [ [   3.87000, -10.70000,  15.60000], [   3.16200,   0.52100,  -1.11600], [  -0.10800,   2.70300,   1.43300], [  -0.50900,  -0.45800,   2.86000] ] ),
                ( 117, [ [   7.33000, -10.10000,  14.40000], [   3.17000,   0.05600,   0.01400], [  -0.09600,   3.07000,   1.14000], [  -2.59000,  -0.76300,   3.56000] ] ),
                ( 118, [ [   1.34000,  -8.28000,  17.80000], [   2.52200,   0.57800,  -0.72600], [   0.32100,   2.86000,   0.62500], [   0.36100,  -0.63500,   1.89000] ] ),
                ( 119, [ [   3.76000,  -7.78000,  17.10000], [   2.18000,   0.81800,  -0.17600], [  -1.05000,   2.81000,   1.59000], [  -0.92300,  -0.82100,   1.77300] ] ),
                ( 120, [ [   5.60000,  -6.96000,  17.00000], [   1.52000,   0.93200,   0.13400], [  -2.95000,   2.55000,   2.23000], [  -1.90000,  -2.20000,   2.72000] ] ),
                ( 121, [ [   1.73000,  -5.15000,  18.50000], [   2.92000,   1.23000,   0.96600], [  -2.87000,   2.35000,   1.05000], [   0.35800,  -2.13000,   1.39000] ] ),
                ( 122, [ [   4.23000, -12.40000,  16.00000], [   2.26000,  -0.31500,  -0.71700], [  -0.40600,   1.32000,   1.51000], [   0.02500,   0.14000,   1.42500] ] ),
                ( 123, [ [   3.87000, -10.70000,  17.60000], [   3.79000,   0.42100,  -0.94000], [  -0.81100,   1.93000,   1.48000], [   0.44100,   0.35000,   1.01000] ] ),
                ( 124, [ [   2.99000,  -8.50000,  18.90000], [   2.71000,   0.49000,  -0.17200], [  -1.57000,   3.81000,   1.44000], [  -0.14300,  -0.36100,   1.17300] ] ),
                ( 125, [ [   2.49300,  -9.65700,   8.51800], [  -1.35500,   0.80500,   0.25900], [  -2.09500,  -2.28700,  -5.51800], [  -0.14100,  -1.15900,   1.75800] ] ),
                ( 126, [ [   1.82900,  -8.51000,   9.08400], [   0.62000,   1.16800,   0.77800], [  -0.05400,  -0.02400,  -0.02000], [  -0.67000,  -0.79900,   1.19500] ] ),
                ( 127, [ [   3.26354,  -8.27035,   9.45759], [   2.53800,  -0.40100,   0.43300], [  -1.64700,   2.16100,   0.67400], [  -1.60100,   0.06100,   1.77100] ] ),
                ( 128, [ [   4.10608,  -5.54413,  10.83025], [  -1.93100,   3.85000,   2.14600], [  -1.36100,   1.34600,   0.44700], [  -3.02000,  -1.80300,   0.23600] ] ),
                ( 129, [ [   0.06500,  -11.43700,  7.21500], [   0.18400,   2.21200,  -0.19300], [  -0.54200,  -0.39700,  -0.29100], [   0.10900,   0.04700,   2.24600] ] ),
                ( 130, [ [  -0.39049,  -9.30579,   6.36823], [   1.13854,   2.34541,   3.26477], [  -1.09706,   0.57877,  -0.80982], [  -1.15689,  -0.02495,   3.62477] ] ),
                ( 131, [ [   1.19500,  -6.75800,   9.68000], [   1.40000,   2.49800,   2.16700], [  -2.14900,   0.85000,  -0.15400], [  -1.51700,  -1.03900,   0.74200] ] ),
                ( 132, [ [   2.49300,  -4.64600,  11.05500], [   1.15800,   1.67200,   0.56600], [  -1.65900,   0.66700,   0.08800], [  -2.05100,  -1.69200,   0.35800] ] ),
                ( 133, [ [  -2.16643,  -8.69021,   6.53180], [   1.24044,   1.63402,   2.28323], [  -2.64042,   2.09831,   1.43965], [  -0.10634,   0.81744,   3.47189] ] ),
                ( 134, [ [  -0.79200,  -6.61800,   9.20000], [   1.48600,   2.37600,   2.30700], [  -1.85100,   0.70100,  -0.13300], [  -0.68500,  -0.48100,   1.57900] ] ),
                ( 135, [ [   0.87700,  -4.12500,  10.97600], [   1.81900,   2.56400,   1.22300], [  -1.73500,   0.43000,  -0.13400], [  -1.77200,  -1.55600,   0.42900] ] ),
                ( 136, [ [  -3.56733,  -7.10525,   8.05865], [   1.29300,   1.54400,   1.32600], [  -1.32900,   1.77000,   1.27100], [   0.40200,   0.69000,   2.47200] ] ),
                ( 137, [ [  -2.27700,  -5.49600,   9.40000], [   1.31400,   1.64800,   1.36800], [  -1.33600,   1.34600,   0.31300], [  -0.46000,  -0.52000,   1.27200] ] ),
                ( 138, [ [  -0.95000,  -3.76000,  10.74000], [   1.33200,   1.75200,   1.40900], [  -1.63800,   0.73000,  -0.22100], [  -1.09700,  -1.31800,   0.38500] ] ),
                ( 139, [ [  -4.69069,  -5.37016,   8.95014], [   1.32600,   1.46500,   0.91600], [  -0.79700,   1.90500,   0.71200], [   0.31300,   0.30300,   1.46900] ] ),
                ( 140, [ [  -3.40400,  -3.97400,   9.81500], [   1.19400,   1.28900,   0.79400], [  -0.86800,   1.93100,   0.39500], [  -0.44800,  -0.51400,   0.82400] ] ),
                ( 141, [ [  -2.26000,  -2.71430,  10.60630], [   1.05900,   1.11500,   0.67500], [  -0.92400,   1.77400,   0.34600], [  -1.01000,  -0.90200,   0.52200] ] ),
                ( 142, [ [  -5.26917,  -3.01247,   9.67344], [   1.25017,   1.16785,   0.36178], [  -0.40957,   2.74678,   0.70837], [   0.09849,   0.11854,   1.06835] ] ),
                ( 143, [ [  -3.88800,  -1.70800,  10.15300], [   1.51100,   1.44000,   0.59700], [  -0.09800,   2.54100,   0.27400], [  -0.42400,  -0.07600,   0.80100] ] ),
                ( 144, [ [  -2.21419,  -0.52570,  10.63808], [   2.02900,   0.25900,   0.44900], [  -1.61300,   2.07400,  -0.10700], [  -1.05300,  -0.66900,   1.14800] ] ),
                ( 145, [ [   1.75100,  -10.58100, 10.76100], [  -0.97700,   1.55100,  -0.18500], [  -2.85000,  -1.04200,  -1.37200], [  -1.33400,  -0.62600,   2.63200] ] ),
                ( 146, [ [   1.22000,  -9.13500,  10.80000], [  -0.03300,   1.25800,   0.27400], [  -2.33700,  -0.38600,  -0.53500], [  -0.52300,  -0.41900,   2.19300] ] ),
                ( 147, [ [   1.51600,  -8.25800,  11.17000], [   0.26000,   0.83600,   0.41900], [  -1.73100,   0.31600,  -0.13700], [  -1.20200,  -0.48000,   2.11000] ] ),
                ( 148, [ [   1.73700,  -7.47700,  11.62600], [   0.18100,   0.72300,   0.49100], [  -0.88900,   1.03000,   0.36900], [  -1.66600,  -1.94500,   1.42400] ] ),
                ( 149, [ [  -0.11300,  -11.26300,  9.86400], [  -0.57200,   2.27600,   0.46100], [  -0.87900,  -0.32100,  -0.42300], [  -0.46500,   0.30100,   3.03200] ] ),
                ( 150, [ [  -0.91900,  -9.25600,  10.19700], [   0.88500,   1.55700,   0.86400], [  -1.03300,   1.70500,   0.19900], [  -0.13400,   0.50100,   3.47100] ] ),
                ( 151, [ [  -0.07900,  -7.79200,  11.03900], [   0.79700,   1.37000,   0.82000], [  -1.44500,   0.61400,  -0.12300], [  -0.92600,  -0.95700,   1.92600] ] ),
                ( 152, [ [   0.67300,  -6.51400,  11.83400], [   0.70700,   1.18500,   0.76900], [  -1.22100,   0.87800,   0.03900], [  -1.53900,  -2.00300,   1.19100] ] ),
                ( 153, [ [  -2.08900,  -7.88200,  10.27500], [   0.72600,   0.57800,   0.57600], [  -1.14200,   1.39600,   0.12000], [   0.24500,   0.67100,   3.27300] ] ),
                ( 154, [ [  -1.35600,  -7.06100,  10.92600], [   0.73000,   1.05600,   0.71900], [  -1.31100,   0.94800,  -0.09300], [  -0.43900,  -0.40200,   1.86600] ] ),
                ( 155, [ [  -0.67100,  -5.75900,  11.69300], [   0.63700,   1.53900,   0.80900], [  -1.33200,   0.75700,  -0.19600], [  -1.29400,  -1.68600,   0.99800] ] ),
                ( 156, [ [  -3.19900,  -6.46900,  10.43600], [   0.46100,   0.43500,   0.34300], [  -1.08300,   1.48100,   0.16100], [   0.35900,   0.55500,   2.29400] ] ),
                ( 157, [ [  -2.66200,  -5.88700,  10.86100], [   0.61400,   0.72800,   0.50500], [  -1.20000,   1.37800,   0.01100], [  -0.30300,  -0.25500,   1.63200] ] ),
                ( 158, [ [  -1.98200,  -5.00400,  11.44300], [   0.74300,   1.03700,   0.66000], [  -1.29100,   1.07400,  -0.18500], [  -0.99200,  -1.13500,   0.95800] ] ),
                ( 159, [ [  -4.25000,  -4.92100,  10.59600], [   0.52600,   0.53800,   0.33800], [  -0.86500,   1.66300,   0.12300], [   0.51300,   0.55600,   1.80000] ] ),
                ( 160, [ [  -3.71400,  -4.32700,  10.95400], [   0.54700,   0.64900,   0.37600], [  -0.82300,   1.98900,   0.12000], [  -0.14600,  -0.16400,   1.40600] ] ),
                ( 161, [ [  -3.16000,  -3.62200,  11.34900], [   0.56100,   0.76100,   0.41300], [  -0.74200,   1.91400,   0.07300], [  -0.68900,  -0.74300,   1.01200] ] ),
                ( 162, [ [  -4.80632,  -2.73141,  10.82287], [   0.67475,   1.45292,   0.16666], [  -0.55718,   1.84148,   0.21418], [   0.59165,   0.14617,   1.35023] ] ),
                ( 163, [ [  -4.16900,  -1.97000,  11.09800], [   0.89300,   0.89100,   0.49400], [  -0.08600,   2.67300,   0.16500], [  -0.11900,  -0.44500,   1.05400] ] ),
                ( 164, [ [  -3.22900,  -1.40800,  11.61400], [   0.94000,   0.22100,   0.51200], [   0.56900,   2.37200,   0.43200], [  -0.95100,  -1.06000,   0.79300] ] ),
                ( 165, [ [  -0.21700,  -10.72500, 13.59000], [   1.15600,   1.58500,  -0.23100], [  -0.98100,  -0.11900,  -0.49200], [  -2.54900,   0.33100,   2.96300] ] ),
                ( 166, [ [   0.88000,  -9.22100,  13.37100], [   1.03800,   1.42300,  -0.20700], [  -1.37100,   0.43200,   0.23800], [  -0.06900,  -1.36500,   2.16900] ] ),
                ( 167, [ [  -0.93000,  -10.81200, 13.23200], [   0.43900,   2.22800,   0.14800], [  -0.44500,  -0.05400,  -0.22300], [  -1.16400,   0.59900,   3.69000] ] ),
                ( 168, [ [  -0.48300,  -8.53800,  13.38300], [   0.45600,   2.31900,   0.15400], [  -1.30900,   0.92000,  -0.22200], [   0.11300,  -0.51400,   2.64700] ] ),
                ( 169, [ [  -1.65200,  -7.41300,  12.91600], [   0.63100,   0.43300,   0.19800], [  -1.20000,   1.28300,  -0.38200], [  -0.15200,  -0.30000,   2.10300] ] ),
                ( 170, [ [  -2.86200,  -5.97800,  12.63700], [   0.71200,   0.49800,   0.26600], [  -1.00400,   1.61400,  -0.19300], [  -0.09600,   0.07200,   1.90500] ] ),
                ( 171, [ [  -3.63227,  -4.45910,  12.49220], [   1.04500,   0.89800,   0.69600], [  -0.68456,   1.67452,  -0.18977], [   0.30547,  -0.09891,   1.64883] ] ),
                ( 172, [ [  -4.20869,  -2.63800,  12.25506], [   0.79000,   0.79000,   0.97500], [  -0.54977,   1.47235,  -0.32850], [   0.03926,  -0.88294,   1.24872] ] )
            ]
            # The number of the elements in the left mouse lung
            elementsCount1 = 2
            elementsCount2 = 4
            elementsCount3 = 4

            # The number of the elements in the diaphragmatic animal lung
            diaphragmaticElementsCount1 = 3
            diaphragmaticElementsCount2 = 5
            diaphragmaticElementsCount3 = 2

            # Create nodes
            nodeIndex = 0
            nodeIdentifier = 1
            leftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []
            diaphragmaticNodeIds = []

            # Left lung nodes
            d1 = [0.5, 0.0, 0.0]
            d2 = [0.0, 0.5, 0.0]
            d3 = [0.0, 0.0, 1.0]
            for n3 in range(elementsCount3 + 1):
                leftNodeIds.append([])
                for n2 in range(elementsCount2 + 1):
                    leftNodeIds[n3].append([])
                    for n1 in range(elementsCount1 + 1):
                        leftNodeIds[n3][n2].append(None)
                        if n3 < elementsCount3:
                            if (n1 == 0) and ((n2 == 0) or (n2 == elementsCount2)):
                                continue
                        else:
                            if (n2 == 0) or (n2 == elementsCount2) or (n1 == 0):
                                continue
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        if generateParameters:
                            x = [0.5 * (n1 - 1), 0.5 * (n2 - 1), 1.0 * n3]
                        else:
                            nodeParameters = nodeFieldParameters[nodeIndex]
                            nodeIndex += 1
                            assert nodeIdentifier == nodeParameters[0]
                            x, d1, d2, d3 = nodeParameters[1]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        leftNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Right lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(rightLung, cache, coordinates, generateParameters,
                 nodes, nodetemplate, nodeFieldParameters,
                 lElementsCount1, lElementsCount2, lElementsCount3,
                 uElementsCount1, uElementsCount2, uElementsCount3,
                 lowerRightNodeIds, upperRightNodeIds, nodeIndex, nodeIdentifier)

            # Diaphragm lung nodes
            nodeIndex, nodeIdentifier = getDiaphragmaticLungNodes(cache, coordinates, generateParameters,
                 nodes, nodetemplate, nodeFieldParameters,
                 diaphragmaticElementsCount1, diaphragmaticElementsCount2, diaphragmaticElementsCount3,
                 diaphragmaticNodeIds, nodeIndex, nodeIdentifier)

            # Create elements
            elementIdentifier = 1

            # Left lung elements
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            leftNodeIds[e3    ][e2][e1], leftNodeIds[e3    ][e2][e1 + 1], leftNodeIds[e3    ][e2 + 1][e1], leftNodeIds[e3    ][e2 + 1][e1 + 1],
                            leftNodeIds[e3 + 1][e2][e1], leftNodeIds[e3 + 1][e2][e1 + 1], leftNodeIds[e3 + 1][e2 + 1][e1], leftNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                        if (e3 < elementsCount3 - 1):
                            if (e2 == 0) and (e1 == 0):
                                # Back wedge elements
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(0)
                                eft = eftWedgeCollapseXi1_15
                            elif (e2 == elementsCount2 - 1) and (e1 == 0):
                                # Front wedge elements
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(2)
                                eft = eftWedgeCollapseXi1_37
                        else:
                            if (e2 == 0) and (e1 == 1):
                                # Top back wedge elements
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(4)
                                eft = eftWedgeCollapseXi2_56
                            elif (e2 == elementsCount2 - 1) and (e1 == 1):
                                # Top front wedge elements
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                eft = eftWedgeCollapseXi2_78
                            elif ((0 < e2 < (elementsCount2 - 1))) and (e1 == 0):
                                # Top middle back wedge element
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                eft = eftWedgeCollapseXi1_57
                            elif (e2 == 0) and (e1 == 0):
                                # Top back tetrahedron element
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(0)
                                eft = eftTetCollapseXi1Xi2_82
                            elif (e2 == elementsCount2 - 1) and (e1 == 0):
                                # Top front tetrahedron element
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(2)
                                eft = eftTetCollapseXi1Xi2_63

                        if eft is eftRegular:
                            element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                        else:
                            elementtemplateCustom.defineField(coordinates, -1, eft)
                            element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                        element.setNodesByIdentifier(eft, nodeIdentifiers)
                        if eft.getNumberOfLocalScaleFactors() == 1:
                            element.setScaleFactors(eft, [-1.0])
                        elementIdentifier += 1
                        leftLungMeshGroup.addElement(element)
                        lungMeshGroup.addElement(element)

            # Right lung elements
            elementIdentifier = getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                rightLungMeshGroup, lowerRightLungMeshGroup, middleRightLungMeshGroup, upperRightLungMeshGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, elementIdentifier)

            # Diaphragm lung elements
            getDiaphragmaticLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                rightLungMeshGroup, diaphragmaticLungMeshGroup,
                diaphragmaticElementsCount1, diaphragmaticElementsCount2, diaphragmaticElementsCount3,
                diaphragmaticNodeIds, elementIdentifier)

            # Marker points
            lungNodesetGroup = lungGroup.getNodesetGroup(nodes)
            markerList = []

            idx = elementsCount1 * elementsCount2 * (elementsCount3 - 1) + elementsCount1 * (elementsCount2 // 2)
            markerList.append({"group": leftApexGroup, "elementId": idx, "xi": [1.0, 1.0, 1.0]})

            idx = elementsCount1 * elementsCount2
            markerList.append({"group": leftVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = elementsCount1 // 2 + elementsCount2
            markerList.append({"group": leftDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            lowerRightLungElementsCount = lElementsCount1 * lElementsCount2 * lElementsCount3
            leftLungElementsCount = elementsCount1 * elementsCount2 * elementsCount3

            idx_1 = uElementsCount1 * uElementsCount2 * (uElementsCount3 - 2) + (uElementsCount1 // 2)
            idx = leftLungElementsCount + lowerRightLungElementsCount + idx_1
            markerList.append({"group": rightApexGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})

            idx = leftLungElementsCount + lowerRightLungElementsCount - (lElementsCount1 // 2)
            markerList.append({"group": rightVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = leftLungElementsCount + (lElementsCount1 // 2)
            markerList.append({"group": rightDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            upperRightLungElementsCount = (uElementsCount1 - 1) * uElementsCount2 * (uElementsCount3 + 1)
            rightLungElementsCount = lowerRightLungElementsCount + upperRightLungElementsCount

            idx_1 = diaphragmaticElementsCount1 * (diaphragmaticElementsCount2 - 1) * (diaphragmaticElementsCount3 - 1)
            idx = rightLungElementsCount + leftLungElementsCount + idx_1
            markerList.append({"group": accessoryApexGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx_1 = diaphragmaticElementsCount1 * (diaphragmaticElementsCount2 - 1) * (
                        diaphragmaticElementsCount3 - 1) - 1
            idx = rightLungElementsCount + leftLungElementsCount + idx_1
            markerList.append({"group": accessoryVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = rightLungElementsCount + leftLungElementsCount + diaphragmaticElementsCount1 + 1

            markerList.append({"group": accessoryDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

        for marker in markerList:
            annotationGroup = marker["group"]
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            cache.setNode(markerPoint)
            markerLocation.assignMeshLocation(cache, mesh.findElementByIdentifier(marker["elementId"]),
                                              marker["xi"])
            markerName.assignString(cache, annotationGroup.getName())
            annotationGroup.getNodesetGroup(nodes).addNode(markerPoint)
            lungNodesetGroup.addNode(markerPoint)
            nodeIdentifier += 1

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCount = options['Refine number of elements']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        parameterSetName = options['Base parameter set']
        isMouse = 'Mouse 1' in parameterSetName
        isHuman = 'Human 1' in parameterSetName

        # create fissure groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        if isHuman:
            upperLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of left lung"))
            lowerLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of left lung"))

            is_upperLeftGroup = upperLeftGroup.getFieldElementGroup(mesh2d)
            is_lowerLeftGroup = lowerLeftGroup.getFieldElementGroup(mesh2d)

            is_obliqueLeftGroup = fm.createFieldAnd(is_upperLeftGroup, is_lowerLeftGroup)
            obliqueLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of left lung"))
            obliqueLeftGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueLeftGroup)

        if isHuman or isMouse:
            upperRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of right lung"))
            middleRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("middle lobe of right lung"))
            lowerRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of right lung"))

            is_upperRightGroup = upperRightGroup.getFieldElementGroup(mesh2d)
            is_middleRightGroup = middleRightGroup.getFieldElementGroup(mesh2d)
            is_lowerRightGroup = lowerRightGroup.getFieldElementGroup(mesh2d)

            is_obliqueRightGroup = fm.createFieldAnd(fm.createFieldOr(is_middleRightGroup, is_upperRightGroup),
                                                     is_lowerRightGroup)
            is_horizontalRightGroup = fm.createFieldAnd(is_upperRightGroup, is_middleRightGroup)

            obliqueRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of right lung"))
            obliqueRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueRightGroup)
            horizontalRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("horizontal fissure of right lung"))
            horizontalRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_horizontalRightGroup)


def getLungNodes(lungSide, cache, coordinates, generateParameters, nodes, nodetemplate, nodeFieldParameters,
                 lElementsCount1, lElementsCount2, lElementsCount3,
                 uElementsCount1, uElementsCount2, uElementsCount3,
                 lowerNodeIds, upperNodeIds, nodeIndex, nodeIdentifier):
    """
    :param lowerNodeIds: nodeIdentifier array in the lower lobe filled by this function
        including indexing by [lElementsCount3 + 1][lElementsCount2 + 1][lElementsCount1 + 1]
    :param upperNodeIds: nodeIdentifier array in the upper lobe filled by this function
        including indexing by [uElementsCount3 + 1][uElementsCount2 + 1][uElementsCount1 + 1]
    :return: nodeIndex, nodeIdentifier
    """
    leftLung = 0

    # Initialise parameters
    d1 = [1.0, 0.0, 0.0]
    d2 = [0.0, 1.0, 0.0]
    d3 = [0.0, 0.0, 1.0]

    # Offset
    xMirror = 0 if lungSide == leftLung else 150

    # Lower lobe nodes
    for n3 in range(lElementsCount3 + 1):
        lowerNodeIds.append([])
        for n2 in range(lElementsCount2 + 1):
            lowerNodeIds[n3].append([])
            for n1 in range(lElementsCount1 + 1):
                lowerNodeIds[n3][n2].append(None)
                if ((n1 == 0) or (n1 == lElementsCount1)) and (n2 == 0):
                    continue
                if (n3 > (lElementsCount3 - 2)) and (n2 > (lElementsCount2 - 2)):
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                if generateParameters:
                    x = [1.0 * (n1 - 1) + xMirror, 1.0 * (n2 - 1), 1.0 * n3]
                else:
                    nodeParameters = nodeFieldParameters[nodeIndex]
                    nodeIndex += 1
                    assert nodeIdentifier == nodeParameters[0]
                    x, d1, d2, d3 = nodeParameters[1]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                lowerNodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

    # Upper lobe nodes
    for n3 in range(uElementsCount3 + 1):
        upperNodeIds.append([])
        for n2 in range(uElementsCount2 + 1):
            upperNodeIds[n3].append([])
            for n1 in range(uElementsCount1 + 1):
                upperNodeIds[n3][n2].append(None)
                if ((n1 == 0) or (n1 == uElementsCount1)) and ((n2 == 0) or (n2 == uElementsCount2)):
                    continue
                if (n2 < (uElementsCount2 - 2)) and (n3 < (uElementsCount3 - 2)):
                    continue
                if ((n2 == 0) or (n2 == uElementsCount2)) and (n3 == uElementsCount3):
                    continue
                if ((n1 == 0) or (n1 == uElementsCount1)) and (n3 == uElementsCount3):
                    continue

                # Oblique fissure nodes
                if (n2 == (uElementsCount2 - 2)) and (n3 < (uElementsCount3 - 2)):
                    upperNodeIds[n3][n2][n1] = lowerNodeIds[n3][lElementsCount2][n1]
                    continue
                elif (n2 < (uElementsCount2 - 1)) and (n3 == (uElementsCount3 - 2)):
                    upperNodeIds[n3][n2][n1] = lowerNodeIds[lElementsCount3][n2][n1]
                    continue

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                if generateParameters:
                    x = [1.0 * (n1 - 1) + xMirror, 1.0 * (n2 - 1) + 2.5, 1.0 * n3 + 2.0]
                else:
                    nodeParameters = nodeFieldParameters[nodeIndex]
                    nodeIndex += 1
                    assert nodeIdentifier == nodeParameters[0]
                    x, d1, d2, d3 = nodeParameters[1]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                upperNodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

    return nodeIndex, nodeIdentifier

def getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular, elementtemplateCustom, mesh,
                    lungMeshGroup, lungSideMeshGroup, lowerLobeMeshGroup, middleLobeMeshGroup, upperLobeMeshGroup,
                    lElementsCount1, lElementsCount2, lElementsCount3,
                    uElementsCount1, uElementsCount2, uElementsCount3,
                    lowerNodeIds, upperNodeIds, elementIdentifier):
    """
    :param lowerNodeIds: Indexing by [lElementsCount3 + 1][lElementsCount2 + 1][lElementsCount1 + 1]
    :param upperNodeIds: Indexing by [uElementsCount3 + 1][uElementsCount2 + 1][uElementsCount1 + 1]
    :return: elementIdentifier
    """

    eftWedgeCollapseXi1_15 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
    eftWedgeCollapseXi1_26 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
    eftWedgeCollapseXi1_57 = eftfactory.createEftWedgeCollapseXi1Quadrant([5, 7])
    eftWedgeCollapseXi1_68 = eftfactory.createEftWedgeCollapseXi1Quadrant([6, 8])
    eftWedgeCollapseXi2_78 = eftfactory.createEftWedgeCollapseXi2Quadrant([7, 8])
    eftTetCollapseXi1Xi2_71 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(7, 1)
    eftTetCollapseXi1Xi2_82 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(8, 2)

    # Lower lobe elements
    for e3 in range(lElementsCount3):
        for e2 in range(lElementsCount2):
            for e1 in range(lElementsCount1):
                eft = eftRegular
                nodeIdentifiers = [
                    lowerNodeIds[e3][e2][e1], lowerNodeIds[e3][e2][e1 + 1], lowerNodeIds[e3][e2 + 1][e1],
                    lowerNodeIds[e3][e2 + 1][e1 + 1],
                    lowerNodeIds[e3 + 1][e2][e1], lowerNodeIds[e3 + 1][e2][e1 + 1],
                    lowerNodeIds[e3 + 1][e2 + 1][e1], lowerNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e2 == 0) and (e1 == 0):
                    # Back wedge elements
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftWedgeCollapseXi1_15
                elif (e2 == 0) and (e1 == (lElementsCount1 - 1)):
                    # Back wedge elements
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(1)
                    eft = eftWedgeCollapseXi1_26
                elif (e3 == 1) and (e2 == (lElementsCount2 - 2)):
                    # Middle wedge
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    eft = eftWedgeCollapseXi2_78
                elif (e3 == (lElementsCount3 - 1)) and (e2 == (lElementsCount2 - 3)):
                    # Remapped cube element 1
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                elif (e3 == (lElementsCount3 - 1)) and (e2 == (lElementsCount2 - 2)):
                    # Remapped cube element 2
                    nodeIdentifiers[2] = lowerNodeIds[e3 - 1][e2 + 1][e1]
                    nodeIdentifiers[3] = lowerNodeIds[e3 - 1][e2 + 1][e1 + 1]
                    nodeIdentifiers[6] = lowerNodeIds[e3 - 1][e2 + 2][e1]
                    nodeIdentifiers[7] = lowerNodeIds[e3 - 1][e2 + 2][e1 + 1]
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, [])])
                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if lowerLobeMeshGroup:
                    lowerLobeMeshGroup.addElement(element)

    # Upper lobe elements
    for e3 in range(uElementsCount3):
        for e2 in range(uElementsCount2):
            for e1 in range(uElementsCount1):
                eft = eftRegular
                nodeIdentifiers = [
                    upperNodeIds[e3][e2][e1], upperNodeIds[e3][e2][e1 + 1], upperNodeIds[e3][e2 + 1][e1],
                    upperNodeIds[e3][e2 + 1][e1 + 1],
                    upperNodeIds[e3 + 1][e2][e1], upperNodeIds[e3 + 1][e2][e1 + 1],
                    upperNodeIds[e3 + 1][e2 + 1][e1], upperNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                    # Distal-front wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(2)
                    eft = eftfactory.createEftBasic()
                    nodes = [3, 4, 7, 8]
                    collapseNodes = [3, 7]
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    # Distal-back wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(3)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [3, 4, 7, 8]
                    collapseNodes = [4, 8]
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e3 == (uElementsCount3 - 2)) and (e2 == 0) and (e1 == 0):
                    # Medial-front wedge elements
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftWedgeCollapseXi1_15
                elif (e3 == (uElementsCount3 - 2)) and (e2 == 0) and (e1 == (uElementsCount1 - 1)):
                    # Medial-back wedge elements
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(1)
                    eft = eftWedgeCollapseXi1_26
                elif (e3 == (uElementsCount3 - 1)) and (0 < e2 < (uElementsCount2 - 1)) and (e1 == 0):
                    # Top-front wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    eft = eftWedgeCollapseXi1_57
                elif (e3 == (uElementsCount3 - 1)) and (0 < e2 < (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    # Top-back wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    eft = eftWedgeCollapseXi1_68
                elif (e3 == (uElementsCount3 - 1)) and (e2 == 0) and (e1 == 0):
                    # Top-front-medial tetrahedron wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftTetCollapseXi1Xi2_82
                elif (e3 == (uElementsCount3 - 1)) and (e2 == 0) and (e1 == (uElementsCount1 - 1)):
                    # Top-back-medial tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(1)
                    eft = eftTetCollapseXi1Xi2_71
                elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                    # Top-front-distal tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(2)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    # remap parameters on xi3 = 1 before collapsing nodes
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, [3], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 4, 4, 4]
                    remapEftLocalNodes(eft, 4, ln_map)

                elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    # Top-front-distal tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(3)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    # remap parameters on xi3 = 1 before collapsing nodes
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, [4], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 4, 4, 4]
                    remapEftLocalNodes(eft, 4, ln_map)

                elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 3)):
                    # Remapped cube element 1
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 2)):
                    # Remapped cube element 2
                    eft = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft, [1, 2], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                elementIdentifier += 1
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if middleLobeMeshGroup and (e3 < (uElementsCount3 - 2)):
                    middleLobeMeshGroup.addElement(element)
                elif upperLobeMeshGroup:
                    upperLobeMeshGroup.addElement(element)

    return elementIdentifier

def getDiaphragmaticLungNodes(cache, coordinates, generateParameters, nodes, nodetemplate, nodeFieldParameters,
                 elementsCount1, elementsCount2, elementsCount3,
                 nodeIds, nodeIndex, nodeIdentifier):
    """
    :parameter:
    :return: nodeIndex, nodeIdentifier
    """

    # Initialise parameters
    d1 = [1.0, 0.0, 0.0]
    d2 = [0.0, 1.0, 0.0]
    d3 = [0.0, 0.0, 1.0]

    # Offset
    xMirror = 75

    # Diaphragmatic lobe nodes
    for n3 in range(elementsCount3 + 1):
        nodeIds.append([])
        for n2 in range(elementsCount2 + 1):
            nodeIds[n3].append([])
            for n1 in range(elementsCount1 + 1):
                nodeIds[n3][n2].append(None)
                if ((n1 == elementsCount1) or (n1 == 1)) and (n3 == elementsCount3):
                    continue
                if (n2 > 1) and (n1 == 0):
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                if generateParameters:
                    if n1 == 0:
                        # The side boxes node
                        x = [1.0 * (n1 - 1) + xMirror, 1.0 * (n2 - 1), 1.0 * n3 + 0.5]
                    elif n3 == (elementsCount3 - 1):
                        # Middle row
                        x = [0.5 * (n1 - 1) + 0.5 + xMirror, 1.0 * (n2 - 1), 1.0 * n3]
                    else:
                        x = [1.0 * (n1 - 1) + xMirror, 1.0 * (n2 - 1), 1.0 * n3]
                else:
                    nodeParameters = nodeFieldParameters[nodeIndex]
                    nodeIndex += 1
                    assert nodeIdentifier == nodeParameters[0]
                    x, d1, d2, d3 = nodeParameters[1]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                nodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

    return nodeIndex, nodeIdentifier

def getDiaphragmaticLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular, elementtemplateCustom,
                    mesh, lungMeshGroup, lungSideMeshGroup, diaphragmaticLobeMeshGroup,
                    elementsCount1, elementsCount2, elementsCount3,
                    NodeIds, elementIdentifier):
    """
    :parameter:
    :return: elementIdentifier
    """

    # Diaphragmatic lobe elements
    for e3 in range(elementsCount3):
        for e2 in range(elementsCount2):
            for e1 in range(elementsCount1):
                eft = eftRegular
                nodeIdentifiers = [
                    NodeIds[e3][e2][e1], NodeIds[e3][e2][e1 + 1], NodeIds[e3][e2 + 1][e1],
                    NodeIds[e3][e2 + 1][e1 + 1],
                    NodeIds[e3 + 1][e2][e1], NodeIds[e3 + 1][e2][e1 + 1],
                    NodeIds[e3 + 1][e2 + 1][e1], NodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e1 == 1) and (e3 == (elementsCount3 - 1)):
                    # wedge elements along crest
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    eft = eftfactory.createEftBasic()
                    nodes = [5, 6, 7, 8]
                    collapseNodes = [5, 7]
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                    if e2 == 0:
                        setEftScaleFactorIds(eft, [1], [])
                        remapEftNodeValueLabel(eft, [3], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [1])])
                    ln_map = [1, 2, 3, 4, 5, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e1 == 1) and (e2 == 0) and (e3 < (elementsCount3 - 1)):
                    # Remapping the elements
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [1])])

                elif (e1 == elementsCount1 - 1) and (e3 == (elementsCount3 - 1)):
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    collapseNodes = [6, 8]
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    ln_map = [1, 2, 3, 4, 5, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e1 == 0) and (e2 == 0):
                    # Remapping the elements
                    if e3 == 0:
                        eft = eftfactory.createEftBasic()
                        setEftScaleFactorIds(eft, [1], [])
                        remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [])])
                    elif (e3 == (elementsCount3 - 1)):
                        nodeIdentifiers[7] = NodeIds[e3 + 1][e2 + 1][e1 + 2]
                        nodeIdentifiers[5] = NodeIds[e3 + 1][e2][e1 + 2]
                        eft = eftfactory.createEftBasic()
                        setEftScaleFactorIds(eft, [1], [])
                        collapseNodes = [6, 8]
                        remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS3,
                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                        remapEftNodeValueLabel(eft, [4], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft, [4], Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [])])

                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                diaphragmaticLobeMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)

    return elementIdentifier
