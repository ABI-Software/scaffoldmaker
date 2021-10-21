'''
Generates 3D lung surface mesh.
'''
import copy
import math
from sympy import symbols, nsolve, solve, Eq
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, remapEftNodeValueLabelsVersion, setEftScaleFactorIds
from scaffoldmaker.utils.cylindermesh import createEllipsePerimeter
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.utils.zinc.general import ChangeManager
from scaffoldmaker.utils.geometry import createEllipsoidPoints, getEllipseRadiansToX, getEllipseArcLength, getApproximateEllipsePerimeter, \
    updateEllipseAngleByArcLength
from scaffoldmaker.utils.interpolation import DerivativeScalingMode
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.vector import magnitude, setMagnitude, crossproduct3, normalise
from opencmiss.utils.zinc.field import Field, findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from opencmiss.zinc.element import Element
from opencmiss.zinc.node import Node




class MeshType_3d_lung1(Scaffold_base):
    '''
    3D lung scaffold.
    '''

    """
        Generates an ellipsoid with a tear-shaped base for the lung mathematically,
         with x, y, z length and the position, angle of the oblique fissure. 
         Regions and markers of the lung are annotated.
    """

    @staticmethod
    def getName():
        return '3D Lung 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        if parameterSetName == 'Default':
            parameterSetName = 'Human 1'
        options['Base parameter set'] = parameterSetName
        options['Refine'] = False
        options['Refine number of elements'] = 4
        options['Discontinuity on the core boundary'] = True
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Refine',
            'Refine number of elements',
            'Discontinuity on the core boundary'
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
        isRat = 'Rat 1' in parameterSetName

        discontinuity = options['Discontinuity on the core boundary']

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
        rightLateralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_lung_term("laterodorsal tip of middle lobe of right lung"))

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

        elif isMouse or isRat:
            # Annotation groups
            diaphragmaticLungGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
            diaphragmaticLungMeshGroup = diaphragmaticLungGroup.getMeshGroup(mesh)
            annotationGroups.append(diaphragmaticLungGroup)

            accessoryApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                    get_lung_term("apex of right lung accessory lobe"))
            accessoryVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term("ventral base of right lung accessory lobe"))
            accessoryDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                      get_lung_term("dorsal base of right lung accessory lobe"))

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
        generateParameters = True
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
            elementIdentifier, leftUpperLobeElementID, leftLowerLobeElementID = \
                getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
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

            if discontinuity:
                # create discontinuity in d3 on the core boundary
                nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
                elementtemplate = mesh.createElementtemplate()
                undefineNodetemplate = nodes.createNodetemplate()
                undefineNodetemplate.undefineField(coordinates)
                nodetemplate = nodes.createNodetemplate()
                fieldcache = fm.createFieldcache()

                with ChangeManager(fm):
                    # localNodeIndexes = [1, 2, 5, 6]
                    valueLabel = Node.VALUE_LABEL_D_DS2

                    for e3 in range(len(leftUpperLobeElementID)):
                        for e2 in range(len(leftUpperLobeElementID[e3])):
                            for e1 in range(len(leftUpperLobeElementID[e3][e2])):

                                # Oblique fissure upper lobe Element
                                if ((e3 < 2) and (e2 == 2)) or ((e3 == 2) and (e2 < 3)):
                                    elementIdentifier = leftUpperLobeElementID[e3][e2][e1]

                                    if e3 < 2:
                                        localNodeIndexes = [1, 2, 5, 6]
                                        valueLabel = Node.VALUE_LABEL_D_DS2
                                    elif (e3 > 1) and (e2 == 0):
                                        localNodeIndexes = [1, 2, 3]
                                        valueLabel = Node.VALUE_LABEL_D_DS3
                                    elif (e3 > 1) and (e2 < 2):
                                        localNodeIndexes = [1, 2]
                                        valueLabel = Node.VALUE_LABEL_D_DS3
                                    elif (e3 > 1) and (e2 < 3):
                                        localNodeIndexes = [1, 2]
                                        valueLabel = Node.VALUE_LABEL_D_DS3


                                    element = mesh.findElementByIdentifier(elementIdentifier)
                                    eft = element.getElementfieldtemplate(coordinates, -1)
                                    nodeIds = get_element_node_identifiers(element, eft)

                                    print('e3',e3,'e2',e2,'e1',e1,'||', elementIdentifier, nodeIds)

                                    for localNodeIndex in localNodeIndexes:
                                        node = element.getNode(eft, localNodeIndex)
                                        nodetemplate.defineFieldFromNode(coordinates, node)
                                        versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1,
                                                                                              valueLabel)
                                        if versionsCount == 1:
                                            if e3 < 2:
                                                fieldcache.setNode(node)
                                                result0, x = coordinates.getNodeParameters(fieldcache, -1,
                                                                                           Node.VALUE_LABEL_VALUE, 1, 3)
                                                result0, d1 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            Node.VALUE_LABEL_D_DS1, 1, 3)
                                                result0, d2 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            valueLabel, 1, 3)
                                                result0, d3 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                                                            3)
                                                result1 = node.merge(undefineNodetemplate)
                                                result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel,
                                                                                                2)
                                                result3 = node.merge(nodetemplate)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_VALUE, 1, x)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS1, 1, d1)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS3, 1, d3)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d2)
                                                result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d2)
                                            elif (e3 > 1) and (e2 < 2):
                                                fieldcache.setNode(node)
                                                result0, x = coordinates.getNodeParameters(fieldcache, -1,
                                                                                           Node.VALUE_LABEL_VALUE, 1, 3)
                                                result0, d1 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            Node.VALUE_LABEL_D_DS1, 1, 3)
                                                result0, d2 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            Node.VALUE_LABEL_D_DS2, 1, 3)
                                                result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1,
                                                                                            3)
                                                result1 = node.merge(undefineNodetemplate)
                                                result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel,
                                                                                                2)
                                                result3 = node.merge(nodetemplate)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_VALUE, 1, x)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS1, 1, d1)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS2, 1, d2)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                                result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                            elif (e3 > 1) and (e2 < 3):
                                                fieldcache.setNode(node)
                                                result0, x = coordinates.getNodeParameters(fieldcache, -1,
                                                                                           Node.VALUE_LABEL_VALUE, 1, 3)
                                                result0, d1 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            Node.VALUE_LABEL_D_DS1, 1, 3)
                                                result0, d2 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                            Node.VALUE_LABEL_D_DS2, 1, 3)
                                                result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1,
                                                                                            3)
                                                result1 = node.merge(undefineNodetemplate)
                                                result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel,
                                                                                                2)
                                                result3 = node.merge(nodetemplate)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_VALUE, 1, x)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS1, 1, d1)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS2, 1, d2)
                                                result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                        Node.VALUE_LABEL_D_DS2, 2, d2)
                                                result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                                result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                            else:
                                                continue

                                    remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)

                                    result1 = elementtemplate.defineField(coordinates, -1, eft)
                                    result2 = element.merge(elementtemplate)
                                    result3 = element.setNodesByIdentifier(eft, nodeIds)

                                    if (e3 == 2) and (e2 < 3):
                                        element.setScaleFactors(eft, [-1.0])

                                else:
                                    fieldcache = fm.createFieldcache()

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

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * lElementsCount3) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 0.0, 1.0]})

        elif isMouse or isRat:
            # valueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ]
            if isMouse:
                nodeFieldParameters = [
                (   1, [ [   0.07403,  -9.56577,   3.60357], [   0.86121,   0.91663,   3.13808], [  -2.20273,  -0.23540,   0.46445], [  -1.94940,  -4.27711,   3.05223] ] ),
                (   2, [ [  -0.26840,  -9.51435,   6.62504], [  -1.32766,  -0.69883,   2.49447], [  -1.88606,   0.44731,  -0.04715], [  -0.69427,  -0.26800,   2.50448] ] ),
                (   3, [ [  -1.92122, -11.21142,   1.24012], [  -0.28342,   2.14970,   2.92579], [  -2.81056,  -0.34804,  -1.27573], [  -2.55210,  -2.28852,   5.04724] ] ),
                (   4, [ [  -2.11315,  -9.49992,   4.16902], [  -0.09808,   1.25545,   2.90775], [  -2.14191,   0.37028,   0.66019], [  -1.47041,  -2.64230,   2.03708] ] ),
                (   5, [ [  -2.12755,  -8.67391,   6.97346], [   0.06838,   0.39146,   2.66632], [  -1.73858,   1.21136,   0.74634], [   0.25644,   0.37400,   2.97672] ] ),
                (   6, [ [  -4.29865, -10.50146,   1.20789], [  -0.07732,   1.64159,   3.92748], [  -2.46054,   1.42565,   0.66974], [  -1.93502,  -1.59827,   5.10560] ] ),
                (   7, [ [  -4.14920,  -8.83455,   4.90568], [   0.37643,   1.68784,   3.45762], [  -1.80302,   1.03452,   0.80364], [  -0.64093,  -1.37872,   1.23650] ] ),
                (   8, [ [  -3.57171,  -7.14590,   8.10918], [   0.77574,   1.68338,   2.93877], [  -1.29608,   1.66826,   0.98992], [   0.15957,   0.71843,   2.11887] ] ),
                (   9, [ [  -6.33487,  -8.48704,   2.61154], [   0.55076,   0.44025,   2.97780], [  -1.73125,   2.62991,   2.59931], [  -1.05677,  -0.95380,   4.51872] ] ),
                (  10, [ [  -5.64451,  -7.49459,   5.73593], [   0.82111,   1.53758,   3.22310], [  -1.53921,   1.63772,   0.76315], [  -1.03865,  -0.57333,   1.61890] ] ),
                (  11, [ [  -4.70306,  -5.37742,   8.94242], [   1.05156,   2.67077,   3.15912], [  -0.89806,   2.08222,   0.80487], [   0.37686,   0.48294,   1.66467] ] ),
                (  12, [ [  -7.18669,  -5.55843,   6.39603], [   0.59266,   2.89465,   3.74833], [  -1.53648,   2.22206,   0.55392], [   0.54053,   0.26298,   3.06847] ] ),
                (  13, [ [  -5.30307,  -3.02550,   9.68332], [   2.97978,   2.03799,   2.65285], [  -0.29788,   2.58614,   0.66778], [   0.94911,   0.77526,   2.35180] ] ),
                (  14, [ [  -1.46058, -12.47910,   7.13881], [   2.09808,   1.60645,   1.90038], [  -1.91804,   1.09320,   0.05112], [  -1.01520,  -1.31998,   3.85443] ] ),
                (  15, [ [  -0.66496,  -9.64655,   9.32022], [  -0.45511,   3.64447,   2.21117], [  -1.11840,   1.37152,   0.44962], [  -0.09056,   0.00680,   2.85604] ] ),
                (  16, [ [  -3.84500, -12.69000,   6.04000], [   0.26136,   0.69836,   0.43532], [  -2.37651,   0.65682,  -0.60958], [  -1.25919,  -0.63612,   4.48080] ] ),
                (  17, [ [  -3.22600, -11.20300,   7.18000], [   0.97561,   2.27288,   1.84296], [  -1.59753,   1.45029,   0.03086], [  -0.58171,  -0.45193,   3.74440] ] ),
                (  18, [ [  -1.88919,  -8.18409,   9.76756], [   1.69766,   3.76414,   3.33147], [  -1.32971,   1.55297,   0.44493], [   0.21986,   0.60501,   2.60653] ] ),
                (  19, [ [  -5.74928, -11.42016,   5.96798], [   0.89023,   1.09247,   0.39442], [  -1.74757,   1.95743,   0.30639], [  -0.94426,  -0.22098,   4.35656] ] ),
                (  20, [ [  -4.63302,  -9.60659,   7.20025], [   1.30646,   2.49071,   2.05426], [  -1.52303,   1.81789,   0.31218], [  -0.20767,   0.09073,   3.12297] ] ),
                (  21, [ [  -3.32590,  -6.54328,  10.20546], [   1.30123,   3.61772,   3.93635], [  -1.19682,   1.64930,   0.44198], [   0.33162,   0.48490,   2.06808] ] ),
                (  22, [ [  -7.01142,  -8.92140,   6.69950], [   0.23567,   0.57897,   0.36479], [  -0.67829,   3.29557,   1.55266], [  -0.28477,   0.09552,   3.60778] ] ),
                (  23, [ [  -6.22043,  -7.60863,   7.85901], [   1.34296,   2.03838,   1.94904], [  -1.05423,   2.27614,   0.99210], [  -0.00580,   0.40454,   2.45985] ] ),
                (  24, [ [  -4.28573,  -4.92425,  10.64091], [   2.52517,   3.32869,   3.61294], [  -0.65090,   2.02589,   0.78353], [   0.45760,   0.42313,   1.73144] ] ),
                (  25, [ [  -6.63445,  -5.23571,   9.11277], [   1.35015,   3.27219,   2.64912], [   0.21711,   2.37051,   1.45456], [   0.56325,   0.38213,   2.36108] ] ),
                (  26, [ [  -4.43830,  -2.62681,  11.78349], [   2.92245,   1.86906,   2.58639], [   0.33621,   2.49812,   1.46020], [   0.76889,   0.01269,   1.81992] ] ),
                (  27, [ [  -1.93100, -12.43600,  10.55900], [   2.55948,   1.63705,   1.27258], [  -1.15116,   2.23905,   0.69133], [  -0.63533,   0.34578,   3.55330] ] ),
                (  28, [ [  -0.41906,  -9.48854,  12.25446], [   0.42218,   3.87081,   1.92577], [  -1.25698,   2.34270,  -0.16936], [  -0.92031,  -0.18621,   2.96052] ] ),
                (  29, [ [  -4.52000, -12.65600,   9.95600], [   1.35217,   2.21656,   0.83086], [  -2.40411,   0.82849,  -0.43856], [  -0.46652,   0.64974,   3.64498] ] ),
                (  30, [ [  -3.13300, -10.25400,  10.92300], [   1.42083,   2.58580,   1.10253], [  -1.23997,   2.09990,   0.02894], [  -0.18554,   0.74076,   3.42815] ] ),
                (  31, [ [  -1.68905,  -7.48625,  12.16850], [   1.46636,   2.94825,   1.38780], [  -1.27716,   1.65104,  -0.00177], [  -0.07484,   0.05607,   2.61833] ] ),
                (  32, [ [  -6.25192, -11.09845,   9.75297], [   2.03540,   3.10903,   0.47581], [  -1.36274,   2.20607,  -0.13783], [  -0.09509,   0.79342,   3.51115] ] ),
                (  33, [ [  -4.37800,  -8.28200,  10.63900], [   1.68764,   2.48600,   1.29046], [  -1.20779,   1.81108,  -0.31319], [   0.13704,   0.85273,   3.27166] ] ),
                (  34, [ [  -2.91467,  -6.17486,  12.23057], [   1.21597,   1.69613,   1.85746], [  -1.06873,   1.49397,   0.11655], [   0.27563,   0.08593,   2.12064] ] ),
                (  35, [ [  -6.98986,  -8.44416,   9.71097], [   1.17990,   1.47658,  -0.09107], [  -0.07160,   3.36699,   0.60791], [   0.36234,   0.70327,   2.96964] ] ),
                (  36, [ [  -5.53200,  -6.64300,  10.30600], [   1.67267,   2.04670,   1.28602], [  -0.91425,   1.84457,   0.24689], [   0.68684,   0.84626,   2.70679] ] ),
                (  37, [ [  -3.78791,  -4.53368,  12.40070], [   1.78494,   2.13535,   2.85449], [  -0.46849,   1.73952,   0.54313], [   0.51030,   0.08286,   1.66357] ] ),
                (  38, [ [  -6.08345,  -4.82242,  11.11703], [   1.72607,   2.92044,   1.90180], [  -0.17690,   1.68471,   1.28954], [   0.70579,   0.15651,   2.08713] ] ),
                (  39, [ [  -3.80003,  -2.88835,  13.25652], [   2.66118,   0.88778,   2.22689], [   0.41935,   1.46422,   1.10304], [   0.90439,  -0.64196,   1.70368] ] ),
                (  40, [ [  -2.72404, -11.77560,  14.14633], [   1.55340,   1.17580,   1.01170], [  -0.80753,   2.16476,  -0.12068], [  -1.17968,   1.29992,   2.61646] ] ),
                (  41, [ [  -2.42583, -10.01198,  14.66971], [  -0.68264,   1.67734,   0.02500], [   0.88132,   1.72091,   0.32221], [  -0.93777,   0.07446,   2.84817] ] ),
                (  42, [ [  -4.78400, -11.50200,  13.19700], [   1.12674,   1.88575,   0.73723], [  -1.79954,   1.09346,  -0.69608], [   0.33288,   1.24262,   2.91702] ] ),
                (  43, [ [  -3.53400,  -9.71000,  13.96000], [   1.36979,   1.69244,   0.78650], [  -0.81167,   1.96451,  -0.25187], [  -0.46943,  -0.15270,   2.45281] ] ),
                (  44, [ [  -2.05298,  -8.12746,  14.76513], [   1.58734,   1.46811,   0.82122], [  -0.19471,   1.93276,  -0.15298], [  -0.28410,  -0.98265,   2.58715] ] ),
                (  45, [ [  -5.98782,  -9.93135,  12.86221], [   1.66750,   2.39142,   0.80255], [  -0.84633,   2.06795,  -0.33823], [   0.83745,   1.49286,   2.98998] ] ),
                (  46, [ [  -4.34300,  -7.84900,  13.64900], [   1.61760,   1.76679,   0.76885], [  -0.66575,   1.88525,  -0.36179], [  -0.14614,  -0.03047,   2.50738] ] ),
                (  47, [ [  -2.78800,  -6.38400,  14.38400], [   1.48492,   1.15738,   0.69764], [  -0.62649,   1.73742,  -0.43431], [   0.41608,  -0.63445,   2.75821] ] ),
                (  48, [ [  -6.29328,  -7.53087,  12.56359], [   1.34012,   1.73997,   0.67057], [   0.34309,   2.67619,   0.14977], [   1.06474,   1.03413,   2.64452] ] ),
                (  49, [ [  -4.86400,  -5.95300,  13.24000], [   1.51296,   1.40866,   0.67950], [  -0.46940,   1.45903,  -0.17125], [   0.24780,  -0.04317,   2.43130] ] ),
                (  50, [ [  -3.28712,  -4.71867,  13.91557], [   1.63291,   1.05492,   0.66842], [   0.16640,   1.63261,   0.55187], [   0.66240,  -0.36288,   1.98688] ] ),
                (  51, [ [  -5.23459,  -4.94899,  13.22187], [   1.99474,   2.05185,   1.35243], [  -0.26306,   0.53137,   0.13065], [   0.63126,  -0.84942,   2.04003] ] ),
                (  52, [ [  -2.63217,  -4.00133,  15.05523], [   2.94186,  -0.14344,   2.12092], [   0.81897,  -0.14175,   1.23719], [   0.63217,  -1.29051,   1.59520] ] ),
                (  53, [ [  -3.97000, -10.25800,  15.67100], [   1.35917,   1.05958,   2.09582], [  -1.01355,   2.02429,   0.72368], [  -0.35939,  -0.84212,   0.86523] ] ),
                (  54, [ [  -2.24601,  -9.44164,  17.27413], [   1.99351,   0.54698,   1.05978], [   0.29630,   1.69013,   1.65367], [  -0.10099,  -1.63017,   2.40788] ] ),
                (  55, [ [  -4.59066,  -8.17597,  15.56314], [   2.11029,   1.21320,   2.47664], [  -0.46995,   1.96813,  -0.40168], [  -0.33349,  -0.59545,   1.26153] ] ),
                (  56, [ [  -1.97927,  -7.59115,  17.60644], [   2.96636,  -0.04152,   1.53435], [  -0.06400,   2.29997,  -0.39224], [   1.19050,  -1.76373,   3.65326] ] ),
                (  57, [ [  -4.90484,  -6.41106,  14.92762], [   1.97629,   1.14920,   1.97034], [  -0.33841,   1.71414,  -1.18187], [  -0.28114,  -0.74487,   0.80545] ] ),
                (  58, [ [  -2.46755,  -5.31416,  16.35299], [   2.79825,   1.00856,   0.85000], [  -0.31319,   1.79797,  -1.34572], [   0.97558,  -0.82711,   2.88452] ] ),
                (  59, [ [   3.21000, -14.40000,   2.61000], [   3.47061,  -1.46960,   0.62900], [   0.01454,   2.77649,   0.53967], [  -0.02362,  -0.96063,   2.35693] ] ),
                (  60, [ [   0.88200, -11.40000,   2.59000], [   1.94636,   0.11461,   1.79858], [  -1.73324,   2.72998,   2.71221], [  -0.13041,  -2.87473,   4.22435] ] ),
                (  61, [ [   3.47000, -12.00000,   3.89000], [   3.04343,  -1.32557,   0.62935], [   0.50438,   1.81716,   1.98022], [   0.09308,  -1.41074,   2.42143] ] ),
                (  62, [ [   6.56000, -14.00000,   3.58000], [   2.98922,  -2.54879,  -1.19066], [   3.06738,   1.43351,   1.57229], [  -0.30031,  -0.62940,   2.43173] ] ),
                (  63, [ [   0.77700, -10.50000,   6.79000], [   2.69445,  -0.40008,  -0.41523], [   0.53427,   1.08623,   3.53337], [  -1.34972,  -1.39050,   1.63022] ] ),
                (  64, [ [   4.14000, -11.00000,   6.27000], [   4.03155,  -0.59992,  -0.62477], [   0.52237,   1.31148,   2.37710], [  -0.05390,  -0.31750,   2.05428] ] ),
                (  65, [ [   8.84000, -11.70000,   5.54000], [   5.36845,  -0.80008,  -0.83523], [   1.48966,   2.82967,   2.32541], [  -0.74699,  -1.06591,   1.06379] ] ),
                (  66, [ [   1.69000,  -9.41000,   9.46000], [   1.68345,  -0.30698,  -0.92692], [   0.11705,   1.29073,   2.70861], [  -0.58039,  -1.45068,   0.75685] ] ),
                (  67, [ [   4.49000,  -9.40000,   8.57000], [   3.85907,   0.33746,  -0.82143], [   0.08264,   2.48591,   2.37816], [   0.07492,   0.46490,   1.26661] ] ),
                (  68, [ [   9.35000,  -8.58000,   8.03000], [   5.82607,   1.29479,  -0.25703], [  -0.49429,   3.27510,   2.27152], [  -0.17912,  -1.49909,   0.25232] ] ),
                (  69, [ [   1.02000,  -8.07000,  11.90000], [   2.69283,   2.66072,  -1.33110], [  -0.93296,   2.77332,   1.35610], [   0.80579,  -0.25500,   0.91207] ] ),
                (  70, [ [   4.18000,  -6.03000,  10.70000], [   3.55302,   1.34601,  -1.03225], [  -1.24668,   3.09689,   1.19257], [   0.32782,  -1.73556,   0.38030] ] ),
                (  71, [ [   7.92000,  -5.44000,   9.89000], [   3.81564,  -0.16130,  -0.57109], [  -2.49577,   3.54450,   1.49256], [   0.43340,  -0.81824,   0.48104] ] ),
                (  72, [ [   3.31000, -15.10000,   5.27000], [   3.34298,  -0.57249,  -0.33456], [   0.25534,   2.10578,   1.10841], [   0.22390,  -0.42791,   2.93496] ] ),
                (  73, [ [   0.91900, -13.40000,   6.67000], [   2.42639,   1.01308,   0.09650], [  -1.82871,   1.84847,   1.92951], [   0.20658,  -1.07744,   3.86536] ] ),
                (  74, [ [   3.66000, -13.07900,   6.59200], [   2.91851,  -0.42832,  -0.25795], [   0.44294,   1.92202,   1.52812], [   0.28553,  -0.72620,   2.94642] ] ),
                (  75, [ [   6.47000, -14.30000,   6.16000], [   2.59687,  -1.93569,  -0.58258], [   2.59812,   1.55550,   0.98686], [   0.12548,   0.04025,   2.68635] ] ),
                (  76, [ [  -0.19700, -11.60000,   8.90000], [   4.37405,   0.97246,  -0.25094], [   0.26249,   1.43921,   2.15124], [  -0.53157,  -0.74077,   2.50920] ] ),
                (  77, [ [   4.19400, -11.28500,   8.31500], [   4.32540,  -0.36082,  -0.91432], [   0.43356,   1.68386,   1.69726], [   0.16201,  -0.25187,   2.03166] ] ),
                (  78, [ [   8.29000, -12.30000,   7.11000], [   3.79245,  -1.63717,  -1.46700], [   1.30230,   2.17472,   1.16973], [  -0.25546,   0.00512,   1.93728] ] ),
                (  79, [ [   0.90300, -10.80000,  10.40000], [   2.93644,   1.64797,   0.11953], [   0.95889,   1.22587,   2.21386], [  -0.98597,  -1.31022,   1.11318] ] ),
                (  80, [ [   4.53300,  -9.71400,   9.97600], [   4.19252,   0.45049,  -0.97287], [   0.17330,   1.72458,   1.50513], [  -0.01153,  -1.23322,   1.16309] ] ),
                (  81, [ [   9.03000, -10.10000,   8.41000], [   4.69558,  -1.19553,  -2.11151], [   0.17676,   2.83868,   1.78605], [  -0.45937,  -1.52825,   0.50555] ] ),
                (  82, [ [   1.40200,  -9.18300,  13.23100], [   2.68026,   1.38653,  -2.45675], [   0.89270,   2.62359,   1.91527], [  -0.25806,  -1.90256,   1.50513] ] ),
                (  83, [ [   4.54000,  -7.87000,  11.30000], [   3.53639,   1.20877,  -1.35085], [  -0.26225,   2.58994,   1.53879], [   0.39022,  -1.93410,   0.81744] ] ),
                (  84, [ [   8.33000,  -6.82000,  10.60000], [   3.96122,   0.87307,  -0.04815], [  -1.45545,   3.69541,   2.08716], [   0.37949,  -1.92834,   0.93107] ] ),
                (  85, [ [   3.67000, -15.20000,   8.43000], [   3.36962,  -0.15600,  -0.77104], [   0.33116,   1.68642,   1.59980], [   0.38621,   0.27271,   2.99273] ] ),
                (  86, [ [   1.24000, -13.70000,  10.10000], [   2.71440,   0.68249,  -0.15581], [  -2.00373,   1.75166,   1.66351], [   0.34382,  -0.04162,   3.33665] ] ),
                (  87, [ [   4.04000, -13.40000,   9.71000], [   2.83798,  -0.09447,  -0.62146], [   0.40411,   1.88945,   0.93731], [   0.51175,   0.21715,   2.70557] ] ),
                (  88, [ [   6.82000, -13.90000,   8.87000], [   2.67979,  -0.89148,  -1.04212], [   2.46763,   1.99209,   0.31201], [   0.36024,   0.67514,   2.48872] ] ),
                (  89, [ [  -0.21500, -11.89600,  11.60400], [   5.01150,   0.66113,  -1.31507], [   0.05959,   2.49326,  -0.11161], [   0.19416,  -0.01108,   2.83486] ] ),
                (  90, [ [   4.46000, -11.50000,  10.30000], [   4.33017,   0.12977,  -1.29074], [   0.25110,   1.93747,   0.11788], [   0.40643,   0.65398,   2.24849] ] ),
                (  91, [ [   8.43000, -11.60000,   9.04000], [   3.60023,  -0.32889,  -1.22599], [   1.08492,   1.98884,  -0.35887], [   0.15797,   1.43960,   2.44040] ] ),
                (  92, [ [   4.07100, -14.59800,  11.20400], [   3.29405,   0.27004,  -1.32010], [   0.57485,   1.17505,   0.74144], [   0.37806,   0.74667,   2.89387] ] ),
                (  93, [ [   1.60000, -13.50000,  13.30000], [   3.13550,   1.18083,  -1.68222], [  -2.02670,   1.67769,   1.66322], [   0.79135,   0.57396,   2.65542] ] ),
                (  94, [ [   4.61000, -12.80000,  11.90000], [   2.84705,   0.20507,  -1.09769], [   0.48779,   2.38956,   0.63076], [   0.15375,   0.38108,   2.50891] ] ),
                (  95, [ [   7.18000, -13.00000,  11.10000], [   2.23932,  -0.59092,  -0.49056], [   2.55415,   2.87602,   0.25575], [   0.02058,   0.88556,   2.61369] ] ),
                (  96, [ [   0.19500, -11.60000,  14.40000], [   4.67000,   1.97000,  -2.39000], [   0.52800,   0.71700,   1.76000], [  -0.16000,  -0.52100,   0.39100] ] ),
                (  97, [ [   4.95300,  -9.85200,  12.34500], [   4.29342,   1.10731,  -1.27324], [   0.09660,   2.65815,   2.51776], [   0.09472,  -2.62528,   0.43948] ] ),
                (  98, [ [   8.57700,  -9.27300,  11.71000], [   2.87364,   0.04930,   0.00315], [  -0.76016,   3.64135,   2.77137], [  -0.46615,  -3.35150,   0.47526] ] ),
                (  99, [ [   0.17500,  -4.79000,  11.40000], [   2.34176,   0.84241,  -0.43968], [  -1.04463,   3.26038,  -0.45954], [   3.44665,  -1.02970,   2.16091] ] ),
                ( 100, [ [   2.42000,  -3.62000,  11.00000], [   2.11991,   1.48740,  -0.35500], [  -2.62542,   2.32119,   0.17888], [   1.83921,  -0.51835,   1.65213] ] ),
                ( 101, [ [   4.36000,  -1.84000,  10.70000], [   1.74102,   2.05015,  -0.24235], [  -4.86044,   2.09150,   0.61269], [   2.19254,  -0.66239,   1.65955] ] ),
                ( 102, [ [-  1.02000,  -1.69000,  11.00000], [   0.77400,   1.81000,  -0.26700], [  -4.16096,   1.50495,  -0.17494], [   2.63934,   0.66375,   3.45319] ] ),
                ( 103, [ [   2.50200,  -6.10300,  13.80300], [   0.95882,   1.44092,  -1.05373], [   0.02231,   3.94853,   0.50210], [   0.99787,  -1.53371,   2.51376] ] ),
                ( 104, [ [   3.88300,  -4.54900,  12.89500], [   1.78322,   1.63709,  -0.74034], [  -1.78301,   3.19172,   1.48926], [   1.01418,  -1.31919,   2.07264] ] ),
                ( 105, [ [   6.06000,  -2.90000,  12.40000], [   2.54207,   1.64236,  -0.24687], [  -3.93726,   2.82979,   1.93979], [   1.12555,  -1.43286,   1.67845] ] ),
                ( 106, [ [   1.07000,  -1.86000,  14.10000], [   0.56000,   2.90000,   0.59200], [  -3.66070,   2.08257,   0.87707], [   1.44454,  -1.02792,   2.62106] ] ),
                ( 107, [ [   2.37000,  -7.49000,  15.90000], [   1.84900,   1.13800,  -0.80053], [   0.80494,   4.52235,   0.81871], [  -0.63021,  -1.13204,   2.09764] ] ),
                ( 108, [ [   4.35000,  -6.18000,  15.00000], [   2.10984,   1.48129,  -0.99896], [  -1.72372,   3.46252,   1.88835], [  -0.06378,  -1.65703,   2.15657] ] ),
                ( 109, [ [   6.58000,  -4.52000,  13.90000], [   2.34926,   1.83801,  -1.20057], [  -3.80548,   3.32178,   2.41313], [  -0.02160,  -2.17176,   2.35359] ] ),
                ( 110, [ [   1.86000,  -3.31000,  16.00000], [   0.37500,   2.65100,   1.04300], [  -3.02069,   2.11270,   0.10357], [   0.39487,  -1.69293,   2.25976] ] ),
                ( 111, [ [   4.42000, -13.70000,  14.20000], [   2.48958,   0.16078,  -0.84608], [  -0.03545,   0.48061,   0.11896], [   0.03557,   1.19311,   2.41818] ] ),
                ( 112, [ [   2.64000, -12.70000,  15.30000], [   1.15268,  -0.20740,  -0.69147], [  -1.75251,   1.41923,   1.33386], [   1.50450,   0.57562,   1.40754] ] ),
                ( 113, [ [   4.27000, -12.70000,  14.60000], [   2.08040,   0.21224,  -0.69238], [  -0.26416,   1.51416,   0.67974], [  -0.16006,   0.25770,   2.06198] ] ),
                ( 114, [ [   6.77000, -12.20000,  14.00000], [   2.89987,   0.78244,  -0.50419], [   1.53892,   2.08249,   0.15167], [  -1.60075,   0.35686,   2.75576] ] ),
                ( 115, [ [   1.07000, -10.90000,  16.80000], [   2.45679,   0.01483,  -1.18169], [  -0.70547,   2.40458,   1.35924], [   2.09623,   0.54531,   1.93302] ] ),
                ( 116, [ [   3.87000, -10.70000,  15.60000], [   3.13502,   0.38512,  -1.21437], [  -0.29331,   2.47238,   1.25316], [  -0.43796,  -0.34293,   2.74080] ] ),
                ( 117, [ [   7.33000, -10.10000,  14.40000], [   3.77802,   0.81338,  -1.18345], [  -0.26994,   3.07283,   1.43474], [  -2.29983,  -0.80189,   3.14117] ] ),
                ( 118, [ [   1.34000,  -8.28000,  17.80000], [   2.64869,   0.29638,  -1.01268], [   0.33185,   2.90989,   0.87613], [   0.50231,  -0.70803,   2.12300] ] ),
                ( 119, [ [   3.76000,  -7.78000,  17.10000], [   2.16042,   0.70016,  -0.37551], [  -1.07609,   2.92854,   1.52888], [  -0.71475,  -1.14019,   1.99526] ] ),
                ( 120, [ [   5.60000,  -6.96000,  17.00000], [   1.48343,   0.91748,   0.17134], [  -2.98647,   2.66755,   2.20947], [  -1.96033,  -2.09494,   2.62983] ] ),
                ( 121, [ [   1.73000,  -5.15000,  18.50000], [   2.92000,   1.23000,   0.96600], [  -2.85335,   2.22944,   1.21550], [   0.58134,  -2.88134,   1.73760] ] ),
                ( 122, [ [   4.23000, -12.40000,  16.00000], [   2.54690,   0.34868,  -0.35124], [  -0.27522,   1.50866,   1.71904], [   0.07752,   0.33143,   0.71459] ] ),
                ( 123, [ [   3.87000, -10.70000,  17.60000], [   3.80757,   0.44270,  -0.93879], [  -0.61445,   1.97382,   1.49213], [   0.39916,   0.31255,   1.14764] ] ),
                ( 124, [ [   2.99000,  -8.50000,  18.90000], [   2.88202,   0.63405,  -0.04194], [  -1.13454,   2.92919,   0.62640], [  -0.80580,  -0.29274,   1.56691] ] ),
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
            elif isRat:
                nodeFieldParameters = [
                (   1, [ [-19.28248,43.39184,-72.39267], [-6.81883,-0.76909,-0.82548], [-0.15439,-3.87742, 2.41863], [-0.30227, 0.15351, 9.41919] ] ),
                (   2, [ [-25.65539,41.36925,-72.92428], [-5.10506,-1.63402, 0.17299], [-1.91151,-2.22103,-0.44048], [ 1.72416, 2.20000, 3.94218] ] ),
                (   3, [ [-12.96631,44.61706,-71.42559], [-4.26267,-6.74508, 2.16272], [ 5.87508,-0.19389, 0.53558], [-4.43226,-7.47402, 9.31360] ] ),
                (   4, [ [-19.75896,39.40753,-70.92949], [-8.79009,-2.83129,-1.44071], [-1.01458,-3.90603, 0.33240], [-4.24313,-1.24331, 4.19554] ] ),
                (   5, [ [-27.54855,38.90758,-73.64003], [-8.62111, 3.45887,-4.73442], [-2.10229,-3.18282,-1.16764], [ 0.95101, 0.67480, 5.68807] ] ),
                (   6, [ [-8.24626,43.25155,-71.30254], [-12.71817,-11.22346, 1.91151], [ 1.71778,-4.25743,-1.34540], [-11.44486,-12.73333, 6.65367] ] ),
                (   7, [ [-21.05045,35.93955,-71.59126], [-12.24577,-2.83183,-2.58580], [-1.60752,-3.52968,-1.01410], [-4.34432,-0.99332, 2.06949] ] ),
                (   8, [ [-29.91127,36.11379,-74.54629], [-6.24348, 4.00417,-4.34221], [-2.36866,-3.57722,-1.90147], [-0.95150, 0.32507, 3.26945] ] ),
                (   9, [ [-11.21335,33.78873,-72.64563], [-11.56976,-2.16560, 2.13990], [-9.22051,-9.35360,-2.54093], [-13.58337,-5.09777, 3.24788] ] ),
                (  10, [ [-22.96489,32.44201,-72.96980], [-11.54077,-0.45437,-2.86084], [-2.41561,-3.68940,-2.10526], [-4.63166,-0.94433, 0.73123] ] ),
                (  11, [ [-33.36703,32.74014,-78.02511], [-9.24213, 1.33599,-7.26902], [-2.47000,-2.97566,-4.54060], [-1.70444, 2.82825, 1.45321] ] ),
                (  12, [ [-25.87378,28.75083,-75.88490], [-12.98867,-0.82493,-6.55124], [-3.37369,-3.66204,-3.69376], [-8.57326,-0.57690,-0.93028] ] ),
                (  13, [ [-35.25635,31.45622,-83.95748], [-5.13848, 5.54699,-8.53429], [-0.94920, 0.07170,-6.62137], [-4.21806, 2.84354,-0.47125] ] ),
                (  14, [ [-21.02558,42.06394,-64.20965], [-3.98901, 2.26598,-1.77711], [-4.24254,-4.20261,-2.89632], [-3.39012,-2.79862, 6.36332] ] ),
                (  15, [ [-24.86894,41.90825,-67.16180], [-3.13373,-2.18427,-3.49772], [-4.55818,-1.12695,-0.44482], [-0.61689,-1.49228, 6.85251] ] ),
                (  16, [ [-18.34908,38.32479,-64.19348], [-6.08605,-1.34198,-3.22096], [ 1.29757,-4.51231,-1.20353], [-6.23950,-4.95237, 4.95352] ] ),
                (  17, [ [-24.23247,38.12881,-66.92476], [-5.54823, 0.97923,-2.17145], [-2.13541,-3.63217,-2.50944], [-4.69719,-1.31217, 3.80731] ] ),
                (  18, [ [-29.15176,40.04655,-68.48561], [-4.14261, 2.75789,-0.91753], [-3.81914,-2.54987,-2.18441], [-1.60487,-0.35500, 5.56713] ] ),
                (  19, [ [-18.76814,34.08153,-66.36278], [-6.62271, 0.34816,-2.92925], [-2.12190,-4.49056,-3.15189], [-9.45396,-5.44544, 3.14158] ] ),
                (  20, [ [-25.46210,34.97410,-69.10687], [-6.73725, 1.43551,-2.54656], [-1.79686,-3.08863,-2.34006], [-4.46586,-0.93458, 2.89304] ] ),
                (  21, [ [-32.18658,36.94750,-71.43195], [-6.68455, 2.50112,-2.09509], [-3.24840,-2.39513,-3.37255], [-1.92413,-0.06207, 4.58901] ] ),
                (  22, [ [-22.79552,30.25432,-70.18097], [-3.57717, 1.06874,-0.38473], [-7.20966,-3.18431,-5.08520], [-9.53308,-1.95308, 1.66999] ] ),
                (  23, [ [-27.83499,32.05697,-71.54908], [-6.45580, 2.52283,-2.34655], [-3.74704,-3.09834,-3.35785], [-5.01558, 0.19321, 2.09554] ] ),
                (  24, [ [-35.53206,35.33892,-75.10899], [-8.91628, 4.03109,-4.76148], [-3.77937,-1.29053,-5.12030], [-2.17581, 1.78496, 4.51906] ] ),
                (  25, [ [-33.09514,29.25775,-75.75029], [-8.90281, 2.61684,-6.26016], [-6.70324,-2.47426,-4.99242], [-5.60883, 1.60827, 1.22778] ] ),
                (  26, [ [-39.47414,34.68000,-81.71061], [-3.52348, 7.51977,-5.17345], [-4.06682,-0.02705,-8.00822], [-3.42306, 3.06845, 5.05374] ] ),
                (  27, [ [-25.17136,38.47093,-60.17458], [-2.10562, 3.15120, 0.48303], [-5.33777,-0.79346,-3.10293], [-5.19258,-3.54021, 3.68313] ] ),
                (  28, [ [-28.12093,40.42136,-60.43981], [-3.31623, 0.65535,-0.88597], [-4.23939,-0.48540,-3.08947], [-4.06590,-1.21596, 5.70877] ] ),
                (  29, [ [-24.52913,34.82634,-61.34773], [-5.32700, 1.50880,-2.93256], [-0.30970,-3.79568,-2.35926], [-6.52953,-2.67733, 2.31139] ] ),
                (  30, [ [-29.14234,36.78669,-63.32945], [-3.78815, 2.38039,-0.96964], [-2.34028,-2.53577,-3.05341], [-4.71837,-1.31337, 3.81152] ] ),
                (  31, [ [-31.93696,39.18917,-63.49316], [-1.68592, 2.26952, 0.60117], [-3.31383,-1.96996,-2.95978], [-2.78541,-1.07349, 5.15423] ] ),
                (  32, [ [-25.96802,31.57611,-64.73225], [-3.58155, 2.46147,-1.34201], [-2.73887,-2.92413,-4.08061], [-6.83837,-1.43898, 1.67464] ] ),
                (  33, [ [-29.94295,34.08261,-65.79517], [-4.35148, 2.53996,-0.77751], [-1.89028,-2.13602,-2.89743], [-4.79528,-0.45542, 3.47994] ] ),
                (  34, [ [-34.65083,36.62302,-66.24347], [-5.04713, 2.53226,-0.11869], [-2.89918,-1.64552,-3.00835], [-2.57245,-0.18697, 5.44332] ] ),
                (  35, [ [-30.03021,29.46471,-69.17341], [-1.71310, 3.36324, 0.78334], [-5.64054,-0.45915,-5.00552], [-6.25633, 0.39180, 1.72188] ] ),
                (  36, [ [-32.81018,32.89477,-68.77757], [-3.78693, 3.37924,-0.01906], [-3.47746,-1.42243,-4.08670], [-4.26984, 0.50306, 2.74373] ] ),
                (  37, [ [-37.56010,35.92441,-69.32334], [-5.61948, 2.63623,-1.05495], [-3.24513,-0.11524,-4.15152], [-2.35156, 0.82690, 5.26770] ] ),
                (  38, [ [-36.76597,31.30014,-74.03135], [-5.95296, 4.26313,-2.73459], [-4.42749,-1.76418,-6.41126], [-3.00306, 1.87518, 3.08653] ] ),
                (  39, [ [-40.88461,36.66219,-74.49937], [-1.98272, 5.60795, 1.56109], [-3.36313, 1.57175,-6.12628], [-1.49782, 1.09399, 6.74011] ] ),
                (  40, [ [-31.29964,35.14734,-57.04794], [-0.71115, 3.37305, 2.18144], [-2.54959, 1.16452,-1.68981], [-6.27721,-1.97946, 1.83340] ] ),
                (  41, [ [-32.67316,39.49667,-55.87883], [-1.94711, 5.09336, 0.14994], [-1.71380,-1.71877,-2.28749], [-4.84906,-1.27534, 3.94664] ] ),
                (  42, [ [-31.17133,32.99811,-59.60737], [-3.11394, 2.37816,-0.18963], [-0.41435,-2.11960,-3.08764], [-6.58628,-0.20556, 2.60348] ] ),
                (  43, [ [-33.65882,35.50518,-59.31332], [-1.77342, 2.56906, 0.78307], [-1.99100,-0.53004,-2.72312], [-3.97962,-1.20281, 3.45864] ] ),
                (  44, [ [-34.70181,37.89624,-58.20712], [-0.29619, 2.09720, 1.35450], [-2.33449,-1.47307,-2.35707], [-2.89756,-0.65441, 5.48021] ] ),
                (  45, [ [-32.21542,31.08525,-63.09201], [-3.01438, 3.26697, 0.51969], [-1.95221,-1.26580,-3.89220], [-7.24248, 1.09167, 2.68557] ] ),
                (  46, [ [-35.01993,34.11512,-62.17701], [-2.57229, 2.76854, 1.30646], [-1.38845,-1.23985,-3.43618], [-4.91431, 0.06716, 3.20905] ] ),
                (  47, [ [-37.32647,36.58752,-60.55346], [-2.01932, 2.15338, 1.92022], [-2.85985,-0.62049,-3.25626], [-3.74977, 0.48245, 7.07252] ] ),
                (  48, [ [-35.03483,30.72582,-67.04592], [-0.22263, 1.45697, 0.49358], [-3.40618, 0.84487,-3.65119], [-5.04087, 2.30419, 2.67275] ] ),
                (  49, [ [-36.37566,33.12141,-66.15938], [-2.44578, 3.24758, 1.25013], [-1.83684,-0.69686,-3.93334], [-3.51992, 0.62857, 2.45551] ] ),
                (  50, [ [-40.15215,36.94891,-64.61725], [-5.07912, 4.38319, 1.82404], [-2.53334, 0.19117,-3.94704], [-3.11243, 1.10132, 5.35555] ] ),
                (  51, [ [-38.66820,32.74574,-69.95618], [-4.25291, 3.52802,-0.97680], [-2.71936,-0.05392,-3.62181], [-1.63009, 1.46782, 5.10388] ] ),
                (  52, [ [-42.39457,36.98693,-68.41099], [-2.87053, 4.44447, 3.64861], [-1.94774,-0.11492,-3.63340], [-1.55993, 0.66816, 7.94747] ] ),
                (  53, [ [-37.09195,34.40484,-56.39785], [-3.44662, 3.06359, 4.17545], [-4.70811,-0.49934,-1.77992], [-2.88548,-0.99747, 2.37133] ] ),
                (  54, [ [-37.69270,37.91658,-52.60956], [ 1.82617, 3.22097, 2.76649], [-5.18328,-0.81757, 1.77987], [-3.06163, 0.69000, 5.67306] ] ),
                (  55, [ [-39.74261,34.21219,-59.36023], [-5.36829, 3.65830, 5.94397], [-1.56161,-0.15557,-4.06315], [-4.52478, 0.12680, 2.42116] ] ),
                (  56, [ [-42.25963,37.81043,-52.24554], [-1.14371, 4.26825, 1.15547], [-3.35380, 0.80208,-2.10884], [-6.10045, 1.95819, 9.51821] ] ),
                (  57, [ [-39.80735,34.12657,-63.89066], [-4.50306, 3.73948, 4.20816], [ 0.45949,-0.65640,-5.40058], [-3.32105, 1.37248, 2.06796] ] ),
                (  58, [ [-43.81423,38.10534,-58.63538], [-3.47183, 4.17137, 6.23263], [-0.37501,-0.30033,-8.37088], [-4.21039, 1.21116, 6.60610] ] ),
                (  59, [ [-9.42530,47.43204,-70.79379], [ 0.57719, 4.76517,-0.20173], [-9.51838, 2.83456, 2.44664], [-4.20508, 1.52435, 4.07203] ] ),
                (  60, [ [-17.91019,45.82126,-70.99168], [-0.88384, 2.84358, 1.59093], [-7.80242,-1.69516,-0.35164], [-0.88061, 0.62507, 2.13826] ] ),
                (  61, [ [-18.41467,49.14604,-70.11403], [-0.07888, 3.65721, 0.08114], [-8.17492, 0.50843,-1.16049], [-0.28615, 0.94395, 1.33018] ] ),
                (  62, [ [-18.02598,52.80329,-70.90142], [ 0.82222, 3.51189,-1.59009], [-8.39128, 3.42554,-1.28845], [-0.94817,-1.30658, 3.10231] ] ),
                (  63, [ [-24.98494,44.64776,-71.58747], [-0.14441, 3.85323,-1.01033], [-6.15487,-2.75398,-1.75322], [-1.04961,-1.25122, 2.56244] ] ),
                (  64, [ [-25.19678,48.65712,-72.61913], [-0.42914, 5.12982,-0.76627], [-5.97065,-1.24710,-2.39105], [-0.62364, 1.12870, 3.50757] ] ),
                (  65, [ [-24.93136,54.36629,-72.78067], [ 0.94523, 6.19185, 0.43638], [-7.38512, 0.54168,-2.74144], [-1.09870,-0.14879, 3.54450] ] ),
                (  66, [ [-29.52826,40.45800,-74.47788], [-0.02291, 6.02040, 0.23701], [-3.66954, 0.05781,-1.19788], [-2.32408, 1.73507, 4.84870] ] ),
                (  67, [ [-30.25395,46.86216,-74.81230], [-1.42818, 6.71047,-0.90890], [-5.39172,-1.13197,-1.94084], [-3.77532, 2.34172, 3.01557] ] ),
                (  68, [ [-32.43734,53.71238,-76.34687], [-2.91328, 6.92972,-2.14162], [-6.19492,-2.97902,-2.70034], [-2.54344, 0.23562, 4.61368] ] ),
                (  69, [ [-34.25227,44.07147,-75.36504], [-1.76757, 2.00629,-0.79570], [-5.17808,-0.15220,-2.44558], [-2.25310, 0.36918, 3.45912] ] ),
                (  70, [ [-35.78348,46.44279,-76.42578], [-1.26699, 2.70474,-1.31322], [-4.89602,-3.11388,-2.89897], [-1.96501, 0.32783, 3.14246] ] ),
                (  71, [ [-36.69687,49.43041,-77.98112], [-0.55416, 3.23764,-1.77940], [-3.99916,-6.07215,-2.73032], [-3.01557,-0.39567, 3.34298] ] ),
                (  72, [ [-14.52547,48.39932,-67.38809], [ 0.31544, 7.09591, 0.70058], [-5.61349, 2.30454,-0.54694], [-5.87905, 0.36809, 2.62683] ] ),
                (  73, [ [-20.40379,45.73571,-68.49881], [-0.36133, 5.03486, 0.43705], [-6.16727,-2.76813,-0.81707], [-3.97937,-0.88647, 2.53858] ] ),
                (  74, [ [-20.53314,49.77990,-68.13399], [ 0.10314, 3.04639, 0.29196], [-6.29842, 0.41416,-0.93479], [-3.86384, 0.03692, 2.22569] ] ),
                (  75, [ [-20.32928,51.81531,-67.92738], [ 0.29945, 1.00717, 0.11922], [-6.58663, 2.93166,-1.08067], [-3.57130,-0.54931, 2.56069] ] ),
                (  76, [ [-26.41363,43.06192,-68.98285], [-0.67781, 6.74215,-0.26423], [-5.31491,-1.94080,-0.93496], [-3.25082,-1.11327, 2.10407] ] ),
                (  77, [ [-26.92290,49.16899,-69.23480], [-0.34058, 5.47051,-0.23961], [-5.92941,-0.61246,-1.35286], [-2.77725,-0.19790, 2.97231] ] ),
                (  78, [ [-27.12173,53.99812,-69.45858], [-0.05706, 4.18585,-0.20785], [-6.70776, 0.91113,-1.79972], [-3.22140,-0.57934, 2.90406] ] ),
                (  79, [ [-30.88480,41.78279,-70.19224], [-1.69434, 7.59902,-0.56631], [-5.12341, 0.59375,-1.23486], [-0.35133, 0.88639, 3.64400] ] ),
                (  80, [ [-32.36687,48.56378,-70.78886], [-1.26956, 5.96188,-0.62683], [-5.02031,-1.19807,-1.52008], [ 0.07832, 0.73349, 4.60890] ] ),
                (  81, [ [-33.43277,53.70501,-71.41772], [-0.86195, 4.31912,-0.63068], [-6.09293,-2.33953,-2.14637], [ 0.74082,-0.26779, 4.90317] ] ),
                (  82, [ [-35.71183,44.45981,-71.20044], [-0.88500, 2.29362,-0.91199], [-5.10898,-0.33068,-2.05892], [-0.58077, 0.39353, 4.73923] ] ),
                (  83, [ [-36.86494,46.87422,-72.23453], [-1.41826, 2.52753,-1.15313], [-4.06757,-3.17004,-1.71812], [-0.08215, 0.51570, 5.05493] ] ),
                (  84, [ [-38.35478,49.39724,-73.51100], [-1.98138, 2.68954,-1.38495], [-4.01075,-5.59524,-1.79491], [-0.46140, 0.54039, 5.31785] ] ),
                (  85, [ [-20.94080,48.11145,-65.72082], [ 0.12914, 6.33903, 0.17324], [-5.29893, 1.05302,-0.70181], [-6.21864,-0.38611, 1.94144] ] ),
                (  86, [ [-25.73088,43.76476,-66.60596], [-1.09664, 5.54917, 0.16140], [-5.17326,-2.95178,-1.00482], [-4.85761,-1.91747, 1.92836] ] ),
                (  87, [ [-26.00667,48.65601,-66.43733], [ 0.56223, 4.14649, 0.17334], [-4.80286, 0.03015,-0.72725], [-4.97385,-0.94351, 1.77423] ] ),
                (  88, [ [-24.90920,51.89351,-66.27686], [ 1.55791, 2.22183, 0.14084], [-5.47474, 2.85464,-0.87298], [-4.88424,-0.01728, 1.84976] ] ),
                (  89, [ [-30.57055,42.24859,-67.58091], [ 0.62940, 6.39412, 0.77120], [-2.63689,-1.15868,-2.83813], [-4.76971, 0.76499, 1.46439] ] ),
                (  90, [ [-30.42882,48.22974,-67.15170], [-0.34953, 5.53184, 0.08284], [-3.62917,-0.04105,-2.55687], [-4.47397,-0.92304, 1.40170] ] ),
                (  91, [ [-31.18463,53.23749,-67.36390], [-1.15121, 4.44165,-0.50249], [-4.58541, 0.97266,-3.25779], [-5.39850,-1.72656, 1.02822] ] ),
                (  92, [ [-26.92969,47.63327,-63.52805], [-0.44991, 6.33276, 0.15801], [-2.59658, 0.52399,-0.61481], [-5.09656,-0.59708, 2.49457] ] ),
                (  93, [ [-30.09697,41.92876,-64.68329], [-0.81200, 6.93981,-0.18633], [-6.18975,-1.78683,-2.06690], [-4.60082,-1.74380, 3.14111] ] ),
                (  94, [ [-30.45491,47.88330,-64.63811], [ 0.09899, 4.94465, 0.27735], [-4.42411,-0.02994,-1.59826], [-4.71360,-0.74800, 2.63499] ] ),
                (  95, [ [-30.06270,51.76820,-64.22987], [ 0.67629, 2.78740, 0.53193], [-6.09486, 2.53823,-2.51577], [-5.28403,-0.64647, 2.39407] ] ),
                (  96, [ [-35.22505,44.81750,-66.30665], [ 0.17126, 2.16318,-0.02729], [-5.86104, 0.30294, 0.20940], [ 3.32567,-1.41150, 4.51742] ] ),
                (  97, [ [-35.62722,47.42796,-66.75754], [-0.98358, 2.95690,-0.87322], [-5.56317,-2.15276,-0.16258], [ 3.68907, 0.57837, 4.34537] ] ),
                (  98, [ [-37.48956,50.60156,-68.03821], [-2.32327, 3.13913,-1.80509], [-5.60467,-4.61630,-0.52488], [ 4.28957, 1.30111, 5.56668] ] ),
                (  99, [ [-37.13421,40.88386,-77.83753], [-1.23430,-0.13641,-2.23973], [-2.70493,-4.18931,-3.98722], [-2.26612, 0.36583, 3.41846] ] ),
                ( 100, [ [-38.40227,41.03416,-79.98313], [-1.29001, 0.43833,-2.03004], [-1.61074,-5.33401,-3.51330], [-1.78889, 1.28351, 5.62558] ] ),
                ( 101, [ [-39.68601,41.73791,-81.86016], [-1.26351, 0.95859,-1.70520], [-0.89889,-7.10270,-2.62863], [-1.80254, 1.20759, 7.23027] ] ),
                ( 102, [ [-39.07545,36.01159,-83.29566], [-2.64042,-0.38794,-2.71184], [ 0.26028,-4.63784,-3.06335], [-2.37516, 0.79706, 8.96019] ] ),
                ( 103, [ [-39.20754,41.47258,-73.47469], [-0.52971, 1.00791,-0.46373], [-2.91236,-4.08168,-1.81735], [-1.85255, 0.80711, 5.26500] ] ),
                ( 104, [ [-39.95620,42.38006,-74.06185], [-0.95351, 0.78021,-0.69825], [-2.25057,-5.24082,-1.24150], [-1.31606, 1.40621, 6.20786] ] ),
                ( 105, [ [-41.09110,42.97218,-74.84851], [-1.29321, 0.39695,-0.85973], [-1.34704,-6.73871,-0.56217], [-1.00486, 1.25909, 6.78188] ] ),
                ( 106, [ [-41.21210,36.72630,-74.63147], [-4.69466, 0.27130,-3.24406], [-0.25337,-5.88425, 0.09919], [-1.89759, 0.63217, 8.36612] ] ),
                ( 107, [ [-40.65986,42.51801,-67.34039], [-0.05377, 1.56232, 0.08330], [-4.22966,-4.25995,-0.11935], [-1.11299, 0.39202, 6.24481] ] ),
                ( 108, [ [-41.00679,43.84009,-67.59459], [-0.63511, 0.93735,-0.59940], [-4.03495,-5.51020, 0.07056], [-1.04220, 0.16418, 6.33567] ] ),
                ( 109, [ [-41.72041,44.24509,-68.32667], [-0.66934,-0.10762,-0.73072], [-3.00257,-7.10457, 0.82621], [-0.95777, 0.49216, 6.17055] ] ),
                ( 110, [ [-42.88231,37.28006,-66.57116], [-4.94741,-0.03310,-2.28429], [ 0.26173,-7.01547, 1.82193], [-1.45119, 0.74470, 7.45725] ] ),
                ( 111, [ [-31.12982,46.96638,-60.89436], [-1.36265, 6.75404,-0.83401], [-4.28404, 0.59745,-0.07587], [-4.58429,-0.30934, 2.97256] ] ),
                ( 112, [ [-34.46256,40.50541,-60.35173], [-1.21190, 8.35222,-0.91170], [-2.02811,-0.10563,-0.39397], [-4.20563, 3.04947, 4.90517] ] ),
                ( 113, [ [-35.24458,47.20606,-61.10836], [-0.35128, 5.04315,-0.60091], [-3.92238,-0.12130,-0.35173], [-2.68630,-0.38190, 4.05762] ] ),
                ( 114, [ [-35.32122,50.56765,-61.52977], [ 0.19617, 1.66456,-0.23969], [-5.20365, 0.72073,-0.86514], [-3.38251,-2.78750, 3.90268] ] ),
                ( 115, [ [-37.86389,44.68804,-60.97787], [-1.24781, 2.09463,-0.54185], [-5.13777, 0.66537,-0.49730], [-2.29286, 0.98629, 5.60903] ] ),
                ( 116, [ [-38.87840,46.76104,-61.57470], [-0.77454, 2.04019,-0.64892], [-3.83694,-2.16607,-0.34534], [-1.93508,-0.36437, 5.64065] ] ),
                ( 117, [ [-39.42268,48.73852,-62.26042], [-0.31178, 1.90117,-0.71737], [-4.12880,-3.32466,-0.61036], [-1.19209,-1.99948, 6.43758] ] ),
                ( 118, [ [-41.41939,42.25915,-61.09277], [-0.36856, 0.24421,-0.36235], [-3.38307,-3.43611, 0.58587], [-1.35009, 0.52323, 7.12082] ] ),
                ( 119, [ [-41.99411,42.81891,-61.66582], [-0.77607, 0.87212,-0.77902], [-2.73276,-4.43819, 0.88663], [-1.14152,-0.17716, 6.74252] ] ),
                ( 120, [ [-42.94079,44.04136,-62.62320], [-1.11559, 1.57038,-1.13400], [-2.62531,-5.68322, 1.22406], [-0.95016,-0.28418, 7.37566] ] ),
                ( 121, [ [-44.11520,38.17068,-59.79145], [-2.82601, 3.17743,-2.57947], [-1.46775,-4.72411, 2.78308], [-0.29194, 3.33108, 6.78326] ] ),
                ( 122, [ [-36.05698,47.08828,-57.59312], [-0.33087, 6.49415,-1.52348], [-4.22413,-0.19610, 2.45611], [ 0.86437, 0.11917, 2.42078] ] ),
                ( 123, [ [-39.45486,46.70341,-55.90009], [-0.52891, 2.95603,-0.83710], [-3.83060,-1.56861, 1.91682], [ 0.75162, 0.23938, 5.48570] ] ),
                ( 124, [ [-43.26144,43.67862,-54.30818], [-2.11455, 1.11413,-0.62146], [-3.58216,-5.45812,-1.65953], [-1.38549, 1.91028, 7.94585] ] ),
                ( 125, [ [-27.80895,42.29657,-71.76398], [-2.39908,-2.21192,-1.83394], [ 2.66008,-2.38403,-0.11581], [-0.62274,-0.65278, 1.68398] ] ),
                ( 126, [ [-29.57739,40.41327,-74.41851], [-2.32766, 0.20361,-0.92860], [ 1.02254,-1.40708, 0.12967], [-0.89774, 0.95716, 2.98675] ] ),
                ( 127, [ [-31.53743,42.08706,-74.75581], [-1.79202, 1.61708,-0.20859], [-0.60802,-2.71699,-0.55044], [ 0.13308,-1.13913, 2.56261] ] ),
                ( 128, [ [-33.07261,43.56021,-74.84381], [-1.27451, 1.32523, 0.03249], [-1.23146, 1.13905,-0.20501], [ 0.18233,-1.49460, 3.76854] ] ),
                ( 129, [ [-26.34498,40.88118,-71.47592], [-2.83460,-3.50334,-1.51836], [ 0.15850,-0.14205,-0.00690], [-2.30608,-0.42110, 3.09443] ] ),
                ( 130, [ [-28.15167,38.16528,-73.75299], [-4.45232,-0.20121,-2.45714], [-2.21565,-2.76126,-1.15169], [-1.71227, 0.39194, 1.70072] ] ),
                ( 131, [ [-32.51213,39.59458,-75.66792], [-3.86135, 3.07822,-1.14798], [-1.32468,-2.19335,-1.25866], [-0.27050,-0.64795, 1.67969] ] ),
                ( 132, [ [-35.18203,43.76194,-75.84548], [-1.35322, 4.81128, 0.72570], [-2.48689,-1.19853,-1.71501], [-0.02335,-3.92115, 1.96130] ] ),
                ( 133, [ [-31.04068,34.87847,-75.47543], [-3.10366, 2.81366,-2.26134], [-2.46397,-2.63993,-2.07332], [-1.51093, 1.62692, 1.31693] ] ),
                ( 134, [ [-34.12141,37.81376,-77.21110], [-3.03149, 3.03308,-1.19085], [-1.48388,-1.70356,-1.73119], [-0.92808,-1.09941, 1.30692] ] ),
                ( 135, [ [-37.04515,40.87618,-77.86950], [-2.78686, 3.05978,-0.12465], [-1.32162,-2.28164,-1.94407], [ 0.67002,-2.74647, 2.34099] ] ),
                ( 136, [ [-33.29146,32.74370,-78.03313], [-1.69664, 3.63708,-1.45299], [-1.76486,-1.61931,-2.80031], [-1.44159, 2.12996, 2.07405] ] ),
                ( 137, [ [-35.46357,36.20698,-79.10730], [-2.45178, 3.19038,-0.91276], [-1.19015,-1.49411,-2.16061], [-0.97405,-0.32188, 1.55115] ] ),
                ( 138, [ [-37.88900,39.21238,-79.55922], [-2.37657, 2.79396, 0.00884], [-0.80895,-1.89919,-2.06378], [ 0.47536,-2.26322, 2.59743] ] ),
                ( 139, [ [-34.47167,31.73662,-80.95570], [-1.89776, 3.53869,-0.56964], [-1.12990,-0.63629,-3.22982], [-2.28249, 3.29646, 2.58848] ] ),
                ( 140, [ [-36.47679,34.85381,-81.50793], [-2.10384, 2.67958,-0.53223], [-1.02133,-1.27183,-2.32907], [-1.77285, 0.67554, 2.54875] ] ),
                ( 141, [ [-38.61104,37.10000,-82.01016], [-2.14794, 1.79880,-0.46858], [-0.56405,-1.64124,-1.89578], [-0.08156,-0.10228, 3.46853] ] ),
                ( 142, [ [-35.29802,31.56232,-84.03479], [-2.52798, 1.97099, 0.19535], [-0.51301, 0.28231,-2.87352], [-5.42833, 3.79002, 1.00935] ] ),
                ( 143, [ [-37.50222,33.66238,-83.76215], [-1.86147, 2.21435, 0.34847], [-1.02924,-1.11071,-2.17874], [-2.49882, 2.04923, 1.65780] ] ),
                ( 144, [ [-39.01443,35.93252,-83.35268], [-1.15121, 2.30243, 0.46572], [-0.24273,-0.69371,-0.78924], [-0.59296, 0.43710, 2.01947] ] ),
                ( 145, [ [-28.89715,42.05184,-69.19026], [-2.65658,-1.46782,-1.04197], [ 0.60494,-1.76289, 0.79032], [-1.29990, 0.10470, 3.01137] ] ),
                ( 146, [ [-30.43873,41.36221,-71.17016], [-1.64770,-0.19798,-0.66127], [ 0.81326,-3.40475,-1.08727], [-0.84583, 0.92548, 3.49209] ] ),
                ( 147, [ [-31.39295,41.49878,-71.55841], [-1.21593, 0.56347,-0.02935], [-0.47735,-2.67446,-1.28525], [ 0.15131, 0.00171, 3.74414] ] ),
                ( 148, [ [-32.51085,42.47326,-71.03994], [-0.94278, 1.28077, 0.98569], [-1.80398,-2.44667,-1.84750], [ 0.93880,-0.65971, 3.78979] ] ),
                ( 149, [ [-28.07969,40.55486,-69.18155], [-1.74206,-2.73513,-1.70182], [ 0.20922,-0.60969, 0.27333], [-1.35261,-0.13534, 2.53856] ] ),
                ( 150, [ [-29.90161,38.37547,-71.60989], [-2.57835,-0.04563,-2.18332], [-2.35314,-2.43925,-1.99758], [-1.79146, 0.12567, 2.52558] ] ),
                ( 151, [ [-32.44716,38.92879,-73.09666], [-2.35217, 1.05436,-0.62351], [-1.61858,-2.39556,-1.75762], [ 0.40667,-0.66871, 3.42414] ] ),
                ( 152, [ [-34.26246,40.22361,-72.97721], [-1.13101, 1.35825, 0.76297], [-1.69543,-2.04746,-2.02314], [ 1.86391,-2.91723, 3.65605] ] ),
                ( 153, [ [-32.52232,36.10170,-73.74075], [-2.17993, 0.06395,-1.74026], [-2.31706,-1.74079,-2.07316], [-1.29996, 0.40785, 2.36583] ] ),
                ( 154, [ [-34.61412,36.85084,-75.00972], [-1.79846, 1.26383,-0.65242], [-1.85325,-1.49457,-1.99000], [ 0.03627,-0.71553, 2.96401] ] ),
                ( 155, [ [-35.89447,38.37865,-75.06796], [-0.67773, 1.59315, 0.47653], [-1.44549,-1.51928,-1.91846], [ 1.61813,-2.19443, 3.21593] ] ),
                ( 156, [ [-34.79444,34.99180,-75.91677], [-1.66778, 0.63109,-1.85512], [-2.05605,-0.95213,-2.33691], [-1.19817, 0.75442, 2.01100] ] ),
                ( 157, [ [-36.11893,35.88082,-76.93059], [-1.42448, 1.29543,-0.49544], [-1.62249,-0.73084,-1.92223], [-0.29531,-0.31677, 2.73641] ] ),
                ( 158, [ [-37.14464,37.16301,-76.78428], [-0.54302, 1.09909, 0.68257], [-1.35877,-0.86864,-1.87048], [ 1.00957,-1.81749, 2.93175] ] ),
                ( 159, [ [-36.58526,34.39201,-78.42312], [-1.41708, 0.83705,-0.67530], [-2.48985,-0.14623,-3.06540], [-1.95359, 1.92906, 2.52147] ] ),
                ( 160, [ [-37.83354,35.39756,-78.82538], [-1.01126, 1.20983,-0.17296], [-1.75536,-0.29613,-2.24756], [-0.91555, 0.40240, 2.78027] ] ),
                ( 161, [ [-38.57230,36.67275,-78.75432], [-0.44691, 1.28492, 0.30200], [-1.27527,-0.49600,-2.12731], [ 0.15971,-0.75139, 3.01472] ] ),
                ( 162, [ [-39.47397,34.72704,-81.75620], [-0.10048, 0.49489, 0.30598], [-3.23871, 0.87092,-3.63881], [-2.25174, 1.66786, 3.18703] ] ),
                ( 163, [ [-39.57581,35.33660,-81.41714], [-0.10304, 0.72344, 0.37165], [-1.71837, 0.17313,-2.91762], [-1.53676, 1.20769, 2.95818] ] ),
                ( 164, [ [-39.67552,36.17645,-81.02003], [-0.09632, 0.95562, 0.42229], [-0.92690,-0.49432,-2.39307], [-0.72528, 0.04786, 2.63242] ] ),
                ( 165, [ [-30.36982,42.42619,-66.24109], [-1.10379,-0.21005,-1.50004], [ 1.39510,-4.44932,-0.86350], [-2.60222, 0.63441, 4.78763] ] ),
                ( 166, [ [-31.24582,42.25949,-67.43156], [-0.64821,-0.12335,-0.88090], [ 0.77823,-4.30120,-0.53797], [ 0.14073, 1.49608, 4.43944] ] ),
                ( 167, [ [-29.64218,40.14305,-66.66589], [-3.42687,-3.30561,-3.99791], [ 0.00778,-0.02482,-0.00482], [-1.17406,-0.66008, 2.60229] ] ),
                ( 168, [ [-31.54090,38.36864,-68.89869], [-0.34167,-0.32959,-0.39861], [-1.46530,-3.06898,-2.29517], [-0.84087,-0.24505, 1.92152] ] ),
                ( 169, [ [-33.64994,36.50360,-71.44255], [-0.39332, 0.82355,-0.16046], [-2.19761,-1.36920,-2.46011], [-0.41980,-0.16235, 1.39866] ] ),
                ( 170, [ [-35.94629,35.83690,-73.80384], [-0.33050, 1.22614,-0.12953], [-2.27997,-0.51240,-2.29570], [-0.93819, 0.15378, 1.43130] ] ),
                ( 171, [ [-38.30971,35.65559,-76.04625], [-0.59359, 1.17190, 0.42562], [-2.29060, 0.19141,-2.20138], [-1.11098, 0.14031, 1.37798] ] ),
                ( 172, [ [-40.44931,35.97948,-78.12884], [-0.02214, 0.94644, 0.08406], [-1.98489, 0.45552,-1.96014], [-0.20060, 0.07449, 3.45244] ] )
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

            idx = leftLungElementsCount + + (lElementsCount1 * lElementsCount2 * lElementsCount3) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 0.0, 1.0]})

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
        isRat = 'Rat 1' in parameterSetName

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

        if isHuman or isMouse or isRat:
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

    if generateParameters is False:
        # lower lobe
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
    else:
        # Initialise parameters
        leftLung = 0
        d1 = [1.0, 0.0, 0.0]
        d2 = [0.0, 1.0, 0.0]
        d3 = [0.0, 0.0, 1.0]

        # Input
        lengthUnit = 6
        widthUnit = 3
        heightUnit = 10
        fissureAngle = 45 # degree
        obliqueProportion = 0.8
        centre = [0.0, 0.0, 0.0]
        centreRight = [0.0, 25.0, 0.0]
        xAxis = [1.0, 0.0, 0.0]
        yAxis = [0.0, 1.0, 0.0]
        zAxis = [0.0, 0.0, 1.0]

        majorXAxis = [xAxis[i] * lengthUnit for i in range(3)]
        minorYAxis = [yAxis[i] * widthUnit for i in range(3)]
        minorZAxis = [zAxis[i] * heightUnit for i in range(3)]

        fissureRadian = fissureAngle/180 * math.pi
        tanFissureRadian = math.tan(fissureRadian)

        # Create 5 points at the centre of 2D Ellipse base
        p1 = [xAxis[i] * lengthUnit for i in range(3)]
        p2 = [xAxis[i] * -lengthUnit for i in range(3)]
        p3 = [[p1[i] - (p1[i] - p2[i]) * obliqueProportion for i in range(3)]] # point of oblique fissure on diaphragm

        C = -tanFissureRadian * p3[0][0]
        X = symbols('x')
        Z = symbols('z')
        eq1 = Eq((X / lengthUnit) ** 2 + (Z / heightUnit) ** 2, 1)
        eq2 = Eq(Z - (tanFissureRadian * X), C)
        result = solve([eq1, eq2], (X, Z))
        p4 = [float(result[1][0]), 0.0, float(result[1][1])] # point of oblique fissure on the curve

        # points on 2D Ellipse base (two on the curve and one centre)
        x = p3[0][0]
        y = math.sqrt(1 - (x/lengthUnit) ** 2) * widthUnit
        p3.append([p3[0][0], y, p3[0][2]])
        p3.insert(0, [p3[0][0], -y, p3[0][2]])

        # points around on 2D Ellipse base
        # --------------------------------------------------- Lower lobe -------------------------------------------------
        lower_row1 = []
        lower_row1_d2 = []
        lower_row2 = []
        lower_row2_d2 = []
        obl = []
        obl_d2 = []
        lower_edge = []
        lower_edge_d2 = []
        lowerObl = []
        lowerObl_d2 = []

        for n2 in range(lElementsCount2 + 1):
            if n2 == 0:
                # row 1
                dx = p3[2][0]
                vec1 = majorXAxis
                vec2 = minorYAxis
                planeCentre = centre
                elementsCountAround = 4
                nodesCountAround = elementsCountAround + 1
                x = lower_row1
                nd2 = lower_row1_d2
            elif n2 == 1:
                # oblique
                dx = 0.0
                vec1 = ([(p4[i] - p3[1][i]) for i in range(3)])
                vec2 = minorYAxis
                elementsCountAround = 4
                planeCentre = p3[1]
                nodesCountAround = elementsCountAround + 1
                x = obl
                nd2 = obl_d2
            elif n2 == 2:
                # lower lobe
                dx = obl[0][0]
                elementsCountAround = 3
                vec1 = majorXAxis
                vec2 = [minorZAxis[i] * -1 for i in range(3)]
                planeCentre = centre
                nodesCountAround = elementsCountAround + 1
                x = lower_edge
                nd2 = lower_edge_d2
            elif n2 == 3:
                # 2nd oblique
                row1_temp = [lower_row1[3][0], 0.0, lower_row1[3][2]]
                vec1 = ([lower_edge[2][i] - row1_temp[i] for i in range(3)])
                vec2 = minorYAxis
                dx = 0.0
                planeCentre = row1_temp
                elementsCountAround = 4
                nodesCountAround = elementsCountAround
                x = lowerObl
                nd2 = lowerObl_d2
            elif n2 == 4:
                # row 2
                y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (lower_edge[1][2] / heightUnit) ** 2) * widthUnit
                row1_centre = [lowerObl[3][0], 0.0, lowerObl[3][2]]
                row1_temp = [0.0, 0.0, lower_edge[1][2]]
                row1_ellipse = [0.0, y, lower_edge[1][2]]
                vec1 = ([lower_edge[1][i] - row1_centre[i] for i in range(3)])
                vec2 = ([row1_ellipse[i] - row1_temp[i] for i in range(3)])
                dx = 0.0
                planeCentre = row1_centre
                elementsCountAround = 3
                nodesCountAround = elementsCountAround + 1
                x = lower_row2
                nd2 = lower_row2_d2

            magMajorAxis = magnitude(vec1)
            magMinorAxis = magnitude(vec2)
            unitMajorAxis = normalise(vec1)
            unitMinorAxis = normalise(vec2)

            radians = 0.0
            totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, dx, 0.5 * math.pi)
            arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
            elementArcLength = arclengthAround / elementsCountAround
            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

            x_temp = []
            nd2_temp = []
            for n1 in range(nodesCountAround):
                cosRadiansAround = math.cos(radians)
                sinRadiansAround = math.sin(radians)

                temp = ([(planeCentre[c] + cosRadiansAround * vec1[c] + sinRadiansAround * vec2[c]) for c in range(3)])
                x_temp.append(temp)

                ndab = setMagnitude([-sinRadiansAround * magMajorAxis, cosRadiansAround * magMinorAxis], -elementArcLength)
                nd2_temp.append(
                    [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])
                radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, -elementArcLength)

            # Gradual changes in node spacing
            if n2 == 0:
                # row 1
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 1, arcLengthDerivatives = True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 1:
                # oblique
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 1, arcLengthDerivatives = True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 2:
                # lower lobe
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2, lengthFractionStart = 0.5, lengthFractionEnd = 0.5, arcLengthDerivatives = True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 3:
                # 2nd oblique
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2,  arcLengthDerivatives = True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 4:
                # rows 2
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2, arcLengthDerivatives = True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            else:
                # even spacing
                [x.append(x_temp[i]) for i in range(len(x_temp))]
                [nd2.append(nd2_temp[i]) for i in range(len(nd2_temp))]

        # complete lower_row2 and lowerObl list
        lower_row2.append(obl[3])
        lower_row2_d2.append(obl_d2[3])
        lowerObl.append(lower_row1[3])
        lowerObl_d2.append(lower_row1[3])

        # smooth derivatives
        tx_d2 = []
        td2 = []
        tx_d3 = []
        td3 = []
        md2 = []
        md3 = []

        for n3 in range(lElementsCount3 + 1):
            tx_d2.append([])
            td2.append([])
            if n3 == 0:
                x = lower_row1
                d2 = lower_row1
            elif n3 == 1:
                x = lower_row2
                d2 = lower_row2
            elif n3 == 2:
                x = lowerObl
                d2 = lowerObl
            elif n3 == 3:
                x = obl
                d2 = obl
            for n2 in range(lElementsCount2 + 1):
                tx_d2[n3].append(x[n2])
                td2[n3].append(d2[n2])
                if n3 == 0:
                    if n2 == 4:
                        reverse_obl = [obl[-i] for i in range(1, len(obl)+1)]
                        tx_d3.append(reverse_obl)
                        td3.append(reverse_obl)
                    else:
                        tx_d3.append([lower_row1[n2], lower_row2[n2], lowerObl[n2], obl[n2]])
                        td3.append([lower_row1[n2], lower_row2[n2], lowerObl[n2], obl[n2]])
                    md3.append(smoothCubicHermiteDerivativesLine(tx_d3[n2], td3[n2]))
            md2.append(smoothCubicHermiteDerivativesLine(tx_d2[n3], td2[n3]))

        # create nodes
        for n3 in range(lElementsCount3 + 1):
            lowerNodeIds.append([])
            for n2 in range(lElementsCount2 + 1):
                lowerNodeIds[n3].append([])
                for n1 in range(lElementsCount1 + 1):
                    lowerNodeIds[n3][n2].append(None)
                    d2 = md2[n3][n2]
                    d3 = md3[n2][n3]

                    # each i row
                    if n3 == 0:
                        x = copy.deepcopy(lower_row1[n2])
                        if n2 < 4:
                            next_xd2 = copy.deepcopy(lower_row1[n2 + 1])
                    elif n3 == 1:
                        x = copy.deepcopy(lower_row2[n2])
                        if n2 < 4:
                            next_xd2 = copy.deepcopy(lower_row2[n2 + 1])
                    elif (n3 == 2) and (n2 < 3):
                        x = copy.deepcopy(lowerObl[n2])
                        if n2 < 4:
                            next_xd2 = copy.deepcopy(lowerObl[n2 + 1])
                    elif (n3 == 3) and (n2 < 3):
                        x = copy.deepcopy(obl[n2])
                        if n2 < 2:
                            next_xd2 = copy.deepcopy(obl[n2 + 1])
                        else:
                            # 3-way point
                            d3 = [-md2[n3][n2][0], -md2[n3][n2][1], -md2[n3][n2][2]]
                            d2 = [md3[n2][n3][0], md3[n2][n3][1], md3[n2][n3][2]]
                    else:
                        continue

                    # skipping the first and last column and row - edge.
                    if (n2 == 0) and (n1 != 1):
                        continue

                    # symmetry
                    if n1 == 0:
                        next_xd1 = [x[0], 0.0, x[2]]
                        d1 = [next_xd1[i] - x[i] for i in range(3)]
                    elif n1 == 1:
                        x = [x[0], 0.0, x[2]]
                        if n2 == 0:
                            d1 = [-d2[i] for i in range(3)]
                        # 3-way point
                        if (n3 == 3) and (n2 == 2):
                            d3 = [-md2[n3][n2][0], 0.0, -md2[n3][n2][2]]
                        elif (n2 < 4):
                            next_xd2[1] = 0.0
                            d2 = [next_xd2[i] - x[i] for i in range(3)]
                            previous_d2 = d2
                        else:
                            d2 = previous_d2
                    else:
                        x = [x[0], -x[1], x[2]]
                        d2 = [d2[0], -d2[1], d2[2]]
                        d3 = [d3[0], -d3[1], d3[2]]

                    # translate right lung to the defined centre of the right lung
                    if lungSide != leftLung:
                        x = [x[i] + centreRight[i] for i in range(3)]

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    lowerNodeIds[n3][n2][n1] = nodeIdentifier
                    nodeIdentifier += 1

        # ------------------------------------------------------ Upper lobe ----------------------------------------------
        upper_row1 = []
        upper_row2 = []
        upper_row3 = []
        upper_row4 = []
        upper_edge = []

        for n2 in range(uElementsCount1):
            if n2 == 0:
                # Upper lobe along the edge
                dx = 0.0
                vec1 = [majorXAxis[i] * -1 for i in range(3)]
                vec2 = minorZAxis
                planeCentre = centre
                elementsCountAround = 5
                nodesCountAround = elementsCountAround + 1
                x = upper_edge
            else:
                # Upper lobe along the edge on the lower lobe
                dx = obl[0][2]
                vec1 = minorZAxis
                vec2 = majorXAxis
                planeCentre = centre
                elementsCountAround = 3
                nodesCountAround = elementsCountAround + 1
                x = upper_edge

            radians = 0.0
            totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, dx, 0.5 * math.pi)
            arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
            elementArcLength = arclengthAround / elementsCountAround
            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

            for n1 in range(nodesCountAround):
                cosRadians = math.cos(radians)
                sinRadians = math.sin(radians)
                x_temp = ([(planeCentre[c] + cosRadians * vec1[c] + sinRadians * vec2[c]) for c in range(3)])
                radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, elementArcLength)
                if (n2 == 1) and (n1 == 0):
                    continue
                x.insert(0, x_temp)

        # Interpolate the nodes between the upper lobe edge and fissures
        for n2 in range(uElementsCount3):
            if n2 == 0:
                # row 1
                dx = -p3[2][0]
                vec2 = minorYAxis
                vec1 = [majorXAxis[i] * -1 for i in range(3)]
                planeCentre = centre
                elementsCountAround = 2
                nodesCountAround = elementsCountAround + 1
                x = upper_row1
            elif n2 == 1:
                # row 2
                y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (upper_edge[-2][2] / heightUnit) ** 2) * widthUnit
                obl_temp = [0.0, 0.0, upper_edge[-2][2]]
                obl_ellipse = [0.0, y, upper_edge[-2][2]]
                dx = 0.0
                vec1 = [upper_edge[-2][i] - obl_temp[i] for i in range(3)]
                vec2 = [obl_ellipse[i] - obl_temp[i] for i in range(3)]
                planeCentre = obl_temp
                elementsCountAround = 2
                nodesCountAround = elementsCountAround + 1
                x = upper_row2
            elif n2 == 2:
                # row 3
                y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (upper_edge[-3][2] / heightUnit) ** 2) * widthUnit
                obl_temp = [0.0, 0.0, upper_edge[-3][2]]
                obl_ellipse = [0.0, y, upper_edge[-3][2]]
                dx = 0.0
                vec1 = [upper_edge[-3][i] - obl_temp[i] for i in range(3)]
                vec2 = [obl_ellipse[i] - obl_temp[i] for i in range(3)]
                planeCentre = obl_temp
                elementsCountAround = 2
                nodesCountAround = elementsCountAround + 1
                x = upper_row3
            else:
                # row 4
                p6_1 = [(upper_edge[-4][i] + upper_edge[1][i]) / 2 for i in range(3)]
                z = math.sqrt(1 - (p6_1[0] / lengthUnit) ** 2 - (p6_1[2] / heightUnit) ** 2) * widthUnit
                p6_0 = [p6_1[0], z, p6_1[2]]
                radius = [upper_edge[-4][i] - p6_1[i] for i in range(3)]
                dx = -magnitude(radius)
                vec1 = [(upper_edge[-4][i] - upper_edge[1][i]) / 2 for i in range(3)]
                vec2 = [p6_0[i] - p6_1[i] for i in range(3)]
                planeCentre = p6_1
                elementsCountAround = 4
                nodesCountAround = elementsCountAround + 1
                x = upper_row4

            radians = 0.0
            totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, dx, 0.5 * math.pi)
            arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
            elementArcLength = arclengthAround / elementsCountAround
            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

            x_temp = []
            nd2_temp = []
            for n1 in range(nodesCountAround):
                cosRadians = math.cos(radians)
                sinRadians = math.sin(radians)

                temp = ([(planeCentre[c] + cosRadians * vec1[c] + sinRadians * vec2[c]) for c in range(3)])
                x_temp.insert(0, temp)

                ndab = setMagnitude([-sinRadiansAround * magMajorAxis, cosRadiansAround * magMinorAxis], -elementArcLength)
                nd2_temp.append(
                    [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])

                radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, -elementArcLength)

            # Gradual changes in node spacing
            if n2 == 0:
                # row 1
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                                 elementLengthStartEndRatio=1, arcLengthDerivatives=True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 1:
                # row 2
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                                 elementLengthStartEndRatio=1.5, arcLengthDerivatives=True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 2:
                # row 3
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                                 elementLengthStartEndRatio=2, lengthFractionStart=0.5,
                                                                 lengthFractionEnd=1, arcLengthDerivatives=True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            elif n2 == 3:
                # row 4
                px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                                 elementLengthStartEndRatio= 0.75 , arcLengthDerivatives=True)
                [x.append(px[i]) for i in range(len(px))]
                [nd2.append(pd1[i]) for i in range(len(pd1))]
            else:
                # even spacing
                [x.append(x_temp[i]) for i in range(len(x_temp))]
                [nd2.append(nd2_temp[i]) for i in range(len(nd2_temp))]

        # smooth derivatives
        tx_d2 = []
        tx_d3 = []
        md2 = []
        md3 = []
        for n3 in range(uElementsCount3 + 1):
            if n3 == 0:
                tx_d2 = upper_row1
                tx_d3 = [upper_edge[-i] for i in range(1,6)]
            elif n3 == 1:
                tx_d2 = upper_row2
                tx_d3 = [upper_row1[-n3-1], upper_row2[-n3-1], upper_row3[-n3-1], upper_row4[-n3-1], upper_edge[-5]]
            elif n3 == 2:
                tx_d2 = upper_row3
                tx_d3 = [upper_row4[n3], upper_edge[n3 + 1]]
            elif n3 == 3:
                tx_d2 = upper_row4
                tx_d3 = [upper_row4[n3-2], upper_edge[n3 - 1]]
            elif n3 == 4:
                # Apex row
                tx_d2 = [upper_edge[i] for i in range(1, 6)]
                tx_d3 = [upper_edge[i] for i in range(1,3)]

            md2.append(smoothCubicHermiteDerivativesLine(tx_d2, tx_d2))
            md3.append(smoothCubicHermiteDerivativesLine(tx_d3, tx_d3))

        # create nodes
        for i in range(uElementsCount3 + 1):
            upperNodeIds.append([])
            for j in range(uElementsCount2 + 1):
                upperNodeIds[i].append([])
                for k in range(uElementsCount1 + 1):
                    upperNodeIds[i][j].append(None)

                    # Oblique fissure nodes
                    if (i < 2) and (j == 2):
                        upperNodeIds[i][j][k] = lowerNodeIds[i][-1][k]
                    elif (i == 2) and (j < 4):
                        upperNodeIds[i][j][k] = lowerNodeIds[-1][j][k]

                    # each i row
                    if (i == 0) and (j > 2):
                        x = copy.deepcopy(upper_row1[j - 2])
                        d2 = md2[i][j-2]
                        idx = j-(j//2)*2
                        d3 = md3[idx][i]
                        if j < 4:
                            next_xd2 = copy.deepcopy(upper_row1[j - 1])
                    elif (i == 1) and (j > 2):
                        x = copy.deepcopy(upper_row2[j - 2])
                        d2 = md2[i][j-2]
                        idx = j-(j//2)*2
                        d3 = md3[idx][i]
                        if j < 4:
                            next_xd2 = copy.deepcopy(upper_row2[j - 1])
                    elif (i == 2) and (j > 2):
                        x = copy.deepcopy(upper_row3[j - 2])
                        d2 = md2[i][j-2]
                        idx = j-(j//2)*2
                        d3 = md3[idx][i]
                        if j < 4:
                            next_xd2 = copy.deepcopy(upper_row3[j - 1])
                    elif (i == 3):
                        x = copy.deepcopy(upper_row4[j])
                        d2 = md2[i][j]
                        d3 = md3[-j-1][i] if j > 2 else md3[-j-1][i-3]
                        # d3 = md3[i][j]
                        if j < 4:
                            next_xd2 = copy.deepcopy(upper_row4[j+1])
                    elif (i == 4) and (j > 0) and (j < 4):
                        x = copy.deepcopy(upper_edge[j+1])
                        d2 = md2[i][j]
                        d3 = md3[-j-1][i] if j == 3 else md3[-j-1][i-3]
                        if j < 4:
                            next_xd2 = copy.deepcopy(upper_edge[j+2])
                    else:
                        continue

                    # skipping the first and last, k.
                    if (((j == 0) or (j == 4)) and (k != 1)) or ((i == 4) and (k != 1)):
                        continue

                    # symmetry
                    if k == 0:
                        # d2 = md2[i][j]
                        pass
                    elif k == 1:
                        if j == 4:
                            # Ridges
                            d2 = previous_d2
                        elif i == 4:
                            # Apex row
                            d3 = [d3[0], 0.0, d3[2]]
                        else:
                            x = [x[0], 0.0, x[2]]
                            next_xd2 = [next_xd2[0], 0.0, next_xd2[2]]
                            d2 = [next_xd2[i] - x[i] for i in range(3)]
                            previous_d2 = d2
                            d3 = [d3[0], 0.0, d3[2]]
                    else:
                        x = [x[0], -x[1], x[2]]
                        d2 = [d2[0], -d2[1], d2[2]]
                        d3 = [d3[0], -d3[1], d3[2]]

                    # translate right lung to the defined centre of the right lung
                    if lungSide != leftLung:
                        x = [x[i] + centreRight[i] for i in range(3)]

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    upperNodeIds[i][j][k] = nodeIdentifier
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

    lowerLobeElementID = []
    upperLobeElementID = []
    # Lower lobe elements
    for e3 in range(lElementsCount3):
        lowerLobeElementID.append([])
        for e2 in range(lElementsCount2):
            lowerLobeElementID[e3].append([])
            for e1 in range(lElementsCount1):
                lowerLobeElementID[e3][e2].append(None)

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
                lowerLobeElementID[e3][e2][e1] = elementIdentifier
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if lowerLobeMeshGroup:
                    lowerLobeMeshGroup.addElement(element)

    # Upper lobe elements
    for e3 in range(uElementsCount3):
        upperLobeElementID.append([])
        for e2 in range(uElementsCount2):
            upperLobeElementID[e3].append([])
            for e1 in range(uElementsCount1):
                upperLobeElementID[e3][e2].append(None)

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
                upperLobeElementID[e3][e2][e1] = elementIdentifier
                elementIdentifier += 1
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if middleLobeMeshGroup and (e3 < (uElementsCount3 - 2)):
                    middleLobeMeshGroup.addElement(element)
                elif upperLobeMeshGroup:
                    upperLobeMeshGroup.addElement(element)

    return elementIdentifier, upperLobeElementID, lowerLobeElementID

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
