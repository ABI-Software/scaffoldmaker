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

        if isHuman:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(upperLeftLungGroup)
        elif isMouse:
            diaphragmaticLungGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
            diaphragmaticLungMeshGroup = diaphragmaticLungGroup.getMeshGroup(mesh)
            annotationGroups.append(diaphragmaticLungGroup)


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
                (   1, [ [ 204.214, 217.387,-341.252], [ -62.947, -35.911,   0.688], [  25.814, -29.305,  29.820], [  -4.512,  17.903,  91.032] ] ),
                (   2, [ [ 269.306, 209.072,-354.889], [ -20.103,  -4.908,  79.411], [  47.047, -49.139,  -2.543], [   0.175,  15.720,  89.961] ] ),
                (   3, [ [ 223.341, 188.444,-308.232], [ -60.426, -24.274,  -1.081], [  20.503, -36.477,  21.902], [   3.351,  18.019,  52.584] ] ),
                (   4, [ [ 175.176, 162.285,-317.267], [ -27.384, -29.319, -26.704], [   9.020, -29.937,  61.149], [   2.480,  21.158,  82.660] ] ),
                (   5, [ [ 292.990, 153.587,-341.674], [ -10.831,  26.074,  66.575], [ -10.247, -58.090,  17.240], [  10.595,  29.035,  61.452] ] ),
                (   6, [ [ 241.787, 149.543,-291.578], [ -78.934, -21.344,  11.789], [  22.813, -35.233,  12.622], [  -6.632,  28.612,  43.054] ] ),
                (   7, [ [ 155.702, 126.322,-273.198], [ -89.521, -38.973,  20.606], [  16.492, -32.096,   4.354], [  14.895,  29.074,  37.867] ] ),
                (   8, [ [ 279.346,  98.455,-327.717], [ -21.666,  16.370,  35.485], [ -18.452, -43.354,   8.934], [  18.541,  53.843,  54.860] ] ),
                (   9, [ [ 251.887, 110.979,-294.259], [ -46.884,  -0.667,  13.029], [  -6.640, -34.923,  -9.542], [  -1.793,  34.831,  57.261] ] ),
                (  10, [ [ 203.263, 108.034,-281.647], [ -46.333, -22.115,  13.236], [  35.945, -18.664,  -6.836], [   5.249,  39.630,  30.209] ] ),
                (  11, [ [ 256.412,  71.152,-323.525], [ -14.537,  11.023,  15.628], [  -7.850, -10.251,  25.280], [  35.613,  35.978,  71.913] ] ),
                (  12, [ [ 243.935,  80.999,-302.872], [ -16.394,   8.598,  18.165], [  -8.074, -19.651,  -6.343], [  10.839,  34.046,  66.123] ] ),
                (  13, [ [ 226.628,  91.702,-285.892], [ -23.230,  -0.285,  15.358], [   7.839, -18.769, -11.448], [   6.974,  29.647,  42.160] ] ),
                (  14, [ [ 217.057, 233.615,-251.001], [ -56.387, -12.798,  -2.074], [   6.551, -32.468,   4.928], [   7.180,   7.996,  92.202] ] ),
                (  15, [ [ 267.567, 218.538,-268.608], [ -31.271, -28.559,  17.780], [  28.703, -41.920, -22.876], [  -2.499,   2.835,  82.368] ] ),
                (  16, [ [ 227.626, 202.773,-250.316], [ -46.421,  -9.158,  11.317], [   7.297, -31.033,   0.920], [  -0.016,   8.237,  63.668] ] ),
                (  17, [ [ 178.210, 194.840,-246.533], [ -50.393,  -5.414,  -8.376], [ -22.308, -44.954,  12.222], [   3.296,  11.647,  70.649] ] ),
                (  18, [ [ 296.250, 178.154,-283.773], [ -52.959,  -4.397,  27.035], [  10.998, -43.061,  -0.027], [  -2.037,   9.722,  56.957] ] ),
                (  19, [ [ 240.706, 174.731,-251.298], [ -65.503, -16.663,  18.653], [  12.413, -26.875,   3.862], [  -0.209,   7.605,  43.189] ] ),
                (  20, [ [ 170.036, 151.299,-240.510], [ -77.888, -18.667,   9.104], [  21.815, -36.197,   2.313], [  11.396,  18.147,  30.385] ] ),
                (  21, [ [ 297.502, 143.355,-275.679], [ -48.044,   9.944,  32.993], [  -5.929, -36.823,  16.652], [  -0.988,  42.077,  47.842] ] ),
                (  22, [ [ 250.978, 148.431,-247.195], [ -50.687,   1.846,   9.041], [   9.032, -28.662,   9.417], [  -4.776,  32.784,  41.932] ] ),
                (  23, [ [ 204.680, 141.636,-246.513], [ -35.493, -20.244,  16.094], [  31.013, -11.624,   3.008], [ -11.195,  25.677,  39.181] ] ),
                (  24, [ [ 287.464, 106.726,-251.655], [ -25.621,   6.219,  23.134], [ -19.294, -29.795,  27.926], [  17.692,  37.852,  69.883] ] ),
                (  25, [ [ 257.922, 116.607,-238.393], [ -31.830,   9.651,   4.808], [  -9.432, -33.031,   3.881], [ -10.300,  42.470,  62.047] ] ),
                (  26, [ [ 228.110, 129.472,-238.391], [ -30.191,   5.166,  -9.731], [  19.291, -27.684,   5.002], [ -37.186,  41.031,  47.384] ] ),
                (  27, [ [ 219.598, 234.911,-158.376], [ -59.865, -18.569, -15.474], [   6.365, -34.542, -21.886], [  -2.376,  -5.948,  82.683] ] ),
                (  28, [ [ 271.479, 212.598,-191.075], [ -45.292, -11.794,   7.859], [  30.495, -31.862, -46.294], [   1.874, -12.175,  77.537] ] ),
                (  29, [ [ 226.886, 201.943,-182.154], [ -46.036,  -0.006,   5.281], [   5.759, -30.320, -27.476], [   4.237, -12.188,  64.641] ] ),
                (  30, [ [ 176.812, 202.108,-180.833], [ -45.198,  -1.262,   0.548], [  -3.376, -48.006, -27.616], [  -1.725,  -4.852,  66.277] ] ),
                (  31, [ [ 291.428, 178.268,-232.360], [ -43.644, -20.853,  24.989], [   9.889, -35.193, -43.486], [   4.353, -16.910,  50.697] ] ),
                (  32, [ [ 237.268, 175.712,-212.165], [ -56.954, -16.756,  10.424], [  13.032, -26.067, -32.916], [   0.506, -11.344,  40.699] ] ),
                (  33, [ [ 175.732, 160.132,-211.466], [ -60.962, -16.488,  -9.264], [  19.400, -26.292, -38.189], [   2.855,   9.850,  36.322] ] ),
                (  34, [ [ 219.695, 225.898, -85.992], [ -73.383,  -7.270,   0.075], [   9.466, -34.972, -20.253], [  -0.443, -16.457,  58.378] ] ),
                (  35, [ [ 276.870, 199.241,-113.455], [ -50.589,  -4.004,   2.780], [  27.986, -43.217, -64.463], [  -2.551, -14.345,  71.169] ] ),
                (  36, [ [ 228.512, 190.564,-119.932], [ -48.632,  -2.241,   2.901], [  10.553, -37.072, -44.027], [  -7.211, -10.078,  60.219] ] ),
                (  37, [ [ 177.143, 186.587,-116.041], [ -45.024, -10.633,   6.490], [  36.691, -40.359, -41.741], [  -7.404, -14.171,  56.958] ] ),
                (  38, [ [ 294.618, 150.473,-187.485], [ -61.912,   3.854,  13.211], [  -7.435, -61.052,  31.811], [   4.313,  46.520,  70.475] ] ),
                (  39, [ [ 237.513, 155.974,-176.526], [ -52.879,   5.296,   5.774], [  -7.421, -52.768,   1.120], [  -9.440,  39.508,  59.010] ] ),
                (  40, [ [ 185.657, 157.799,-171.994], [ -42.801,   6.979,   3.264], [  17.078, -26.547,   0.836], [  -6.732,  20.365,  70.432] ] ),
                (  41, [ [ 246.491,  63.880,-307.113], [  -7.559,   4.781,   7.258], [ -10.994,  -7.365,  11.231], [  30.221,  28.933,  71.221] ] ),
                (  42, [ [ 239.325,  70.695,-300.571], [  -8.415,   6.741,   6.479], [  -7.094,  -9.314,   3.351], [  10.779,  13.380,  72.272] ] ),
                (  43, [ [ 231.513,  77.178,-292.656], [ -12.525,   4.255,   3.319], [  -0.928, -16.243,  -4.197], [  11.706,  13.532,  59.617] ] ),
                (  44, [ [ 230.978,  62.444,-297.245], [  -9.805,   9.975,   8.601], [  -8.962,  -8.354,   1.170], [  -5.556,   2.237,  55.703] ] ),
                (  45, [ [ 258.329,  78.766,-234.210], [  -8.603,   4.498,   0.874], [ -37.867, -20.422,   2.056], [   3.424,   4.618,  67.449] ] ),
                (  46, [ [ 248.799,  84.344,-233.104], [ -10.725,   6.321,   0.445], [ -24.858, -21.841,  -4.453], [ -12.855,  -3.519,  64.769] ] ),
                (  47, [ [ 238.275,  93.157,-232.136], [ -12.861,   7.971,  -1.034], [  -9.789, -33.602, -12.451], [  -9.411,  14.946,  61.906] ] ),
                (  48, [ [ 223.573,  66.187,-240.080], [ -17.027,  20.665,  -0.831], [ -28.759,  -7.073,  -9.354], [ -34.523,   0.750,  50.708] ] ),
                (  49, [ [ 254.223,  82.226,-174.237], [ -21.821,  14.346,  -0.313], [ -70.819, -36.828,  -4.473], [  -6.474,   1.765,  58.654] ] ),
                (  50, [ [ 232.669,  96.602,-174.812], [ -20.714,  13.936,   0.272], [ -45.067, -35.411,  -8.030], [ -18.872,  10.884,  54.020] ] ),
                (  51, [ [ 211.888, 109.358,-175.186], [ -20.256,  12.504,   0.821], [ -16.313, -50.413, -13.321], [ -35.462,  14.924,  45.481] ] ),
                (  52, [ [ 187.821,  69.713,-187.140], [ -17.219,  27.275,   2.313], [ -45.037, -16.749, -14.126], [ -30.897,   6.798,  50.061] ] ),
                (  53, [ [ 213.425, 207.382, -42.148], [ -56.500,   0.342,  -5.827], [   6.048, -18.275, -16.938], [  -5.756, -22.958,  37.324] ] ),
                (  54, [ [ 258.130, 182.777, -53.571], [ -32.759,  -3.828,  -5.952], [  46.842, -36.257, -14.249], [ -42.970, -18.327,  45.780] ] ),
                (  55, [ [ 221.272, 179.757, -61.791], [ -41.743,  -3.435,  -5.875], [  10.486, -36.897, -21.690], [  -5.754, -11.017,  49.078] ] ),
                (  56, [ [ 175.167, 176.300, -67.698], [ -50.920,  -3.892,   0.663], [  -2.971, -33.698, -41.085], [  -2.018, -13.036,  51.511] ] ),
                (  57, [ [ 270.017, 129.272, -88.096], [ -48.699,  18.376,  -7.516], [ -17.418, -51.841, -36.718], [ -50.518, -29.109,  80.611] ] ),
                (  58, [ [ 224.626, 141.720, -98.406], [ -43.872,   3.149,  -4.298], [  -5.587, -42.256, -31.773], [   2.711, -18.020,  68.031] ] ),
                (  59, [ [ 185.274, 147.077,-102.145], [ -35.411,  -3.106,  -4.228], [  15.191, -29.940, -31.756], [ -14.714,  -1.454,  64.340] ] ),
                (  60, [ [ 236.417,  87.160,-119.825], [ -26.717,  14.046,  -6.516], [ -56.297, -42.646, -20.424], [ -33.135,   2.045,  67.489] ] ),
                (  61, [ [ 209.605, 101.124,-126.121], [ -27.728,  12.727,  -4.885], [ -42.756, -25.066, -21.644], [ -36.638,   1.272,  45.800] ] ),
                (  62, [ [ 181.792, 113.536,-131.292], [ -27.851,  13.168,  -2.607], [  -7.595, -34.516,  -8.836], [ -30.082,  -2.456,  33.978] ] ),
                (  63, [ [ 161.721,  78.671,-141.179], [ -21.726,  42.029, -15.240], [ -51.129, -20.611, -10.957], [  -8.563,  21.673,  52.565] ] ),
                (  64, [ [ 203.028, 174.619, -25.286], [ -60.155,  -2.415,  -3.955], [ -17.934, -44.396,   2.254], [  -5.864,   0.361,  32.179] ] ),
                (  65, [ [ 189.729, 132.313, -57.386], [ -66.731,  15.839, -24.611], [  -7.400, -14.578, -13.799], [ -31.717,   2.116,  31.478] ] ),
                (  66, [ [ 162.058, 109.623, -84.659], [ -29.742,  25.246, -50.572], [ -39.636, -28.589, -27.407], [ -49.298,   7.984,  46.787] ] ),
                (  67, [ [ 112.805, 220.636,-344.408], [ -57.668,  31.639, -21.158], [ -27.490,  -7.728,  31.544], [  -7.261,  25.118,  94.160] ] ),
                (  68, [ [ 138.804, 176.487,-317.842], [ -42.283,  17.815,  25.150], [  35.114, -29.696,  45.703], [ -14.588,  30.925,  88.491] ] ),
                (  69, [ [  91.579, 200.374,-316.049], [ -54.312,  26.071, -18.610], [ -21.624, -23.063,  24.318], [   4.200,  -7.127,  69.610] ] ),
                (  70, [ [  45.375, 218.344,-353.416], [ -33.266,   4.881, -54.504], [ -61.732, -41.224,   6.887], [   0.041,  11.357,  87.726] ] ),
                (  71, [ [ 157.132, 141.529,-272.796], [ -81.449,  61.031, -11.075], [ -12.618, -42.064,  37.035], [  -3.164,  13.541,  40.904] ] ),
                (  72, [ [  74.446, 172.925,-295.847], [ -76.526,   2.015, -34.054], [ -12.428, -38.556,  16.527], [   4.041,  -7.852,  55.914] ] ),
                (  73, [ [  19.591, 150.815,-334.755], [  -3.346, -58.150, -34.776], [  -2.241, -66.932,  16.771], [ -37.214,  32.111,  65.155] ] ),
                (  74, [ [ 110.041, 132.222,-267.112], [ -38.045,  -0.257,   4.674], [ -17.685, -35.955,   7.777], [  13.262,   5.765,  25.693] ] ),
                (  75, [ [  69.151, 126.208,-283.827], [ -43.252, -15.831, -25.640], [   5.343, -46.939,   6.574], [  -2.573,   8.979,  43.158] ] ),
                (  76, [ [  36.295,  94.586,-323.318], [ -23.915, -46.011, -46.866], [  28.030, -41.990,  18.732], [ -35.404,  53.450,  65.721] ] ),
                (  77, [ [ 102.468,  96.616,-272.124], [ -19.736, -13.518,  -6.858], [   7.724, -16.668,   0.398], [   8.984,  23.773,  34.657] ] ),
                (  78, [ [  83.599,  83.999,-282.740], [ -14.132, -14.812, -15.431], [  32.211, -10.811,   2.136], [ -23.147,  18.638,  18.839] ] ),
                (  79, [ [  67.542,  71.106,-298.778], [ -14.808, -11.987, -18.636], [  32.574, -12.596,  14.591], [ -50.592,  24.213,  13.600] ] ),
                (  80, [ [ 109.300, 234.171,-248.531], [ -52.302,  25.979, -18.010], [ -19.303, -36.966,   0.740], [  -4.582,   1.836,  96.421] ] ),
                (  81, [ [ 135.331, 190.757,-238.135], [ -39.821,   0.430,  -4.178], [  22.077, -42.173,   7.556], [  -3.889,   9.995,  73.763] ] ),
                (  82, [ [  91.699, 199.576,-247.800], [ -46.192,  17.407, -15.151], [ -16.031, -32.155,   0.697], [  -3.291,   2.421,  66.468] ] ),
                (  83, [ [  46.055, 227.254,-268.328], [ -44.924,  35.895, -26.293], [ -53.545, -36.498, -12.333], [   1.708,   5.118,  82.513] ] ),
                (  84, [ [ 152.423, 152.232,-233.455], [ -80.326,  24.447,  -7.950], [  -7.571, -34.911,  -1.433], [  -4.926,  12.512,  36.525] ] ),
                (  85, [ [  77.265, 169.881,-247.131], [ -69.352,  10.546, -19.469], [ -10.609, -28.577,   2.252], [   1.874,   0.636,  41.232] ] ),
                (  86, [ [  14.896, 174.390,-270.710], [ -54.454,  -2.347, -27.972], [ -12.870, -45.688,   3.505], [  -0.893,  15.664,  55.197] ] ),
                (  87, [ [ 128.461, 142.532,-236.966], [ -58.190,   0.601,   0.855], [ -23.184, -13.581,  -4.697], [  19.504,  21.536,  33.741] ] ),
                (  88, [ [  70.179, 142.990,-243.592], [ -57.694,   0.305, -14.313], [  -4.787, -27.614,   2.693], [   4.708,  24.314,  36.010] ] ),
                (  89, [ [  14.153, 142.929,-264.777], [ -53.191,  -0.026, -29.220], [   3.625, -31.174,   8.060], [   4.661,  44.056,  51.260] ] ),
                (  90, [ [ 112.858, 122.968,-231.669], [ -45.051, -10.899,  -8.599], [ -10.735, -28.480,   2.607], [  14.630,  28.758,  45.554] ] ),
                (  91, [ [  67.804, 115.230,-241.806], [ -44.888,  -4.537, -11.644], [   7.971, -30.404,   4.318], [  -6.116,  42.327,  60.537] ] ),
                (  92, [ [  23.424, 113.892,-254.850], [ -43.703,   1.854, -14.388], [  21.699, -35.829,  12.761], [ -29.770,  57.684,  72.277] ] ),
                (  93, [ [ 104.207, 234.698,-152.733], [ -57.745,  26.895,  -7.074], [ -25.137, -25.018, -30.805], [   0.558,  -2.469,  77.755] ] ),
                (  94, [ [ 135.749, 192.199,-173.993], [ -57.597,   9.559, -12.340], [  21.057, -36.953, -26.041], [   2.245,  -7.457,  61.353] ] ),
                (  95, [ [  84.928, 206.459,-183.823], [ -43.517,  18.872,  -7.207], [ -13.067, -31.107, -30.940], [   0.766,  -1.586,  66.478] ] ),
                (  96, [ [  49.256, 228.018,-188.624], [ -27.295,  23.783,  -2.350], [ -46.031, -28.715, -40.634], [   2.797,  -2.517,  77.864] ] ),
                (  97, [ [ 146.834, 163.963,-199.827], [ -73.680,   9.289, -14.410], [  -1.766, -30.206, -31.659], [ -15.120,   1.557,  27.005] ] ),
                (  98, [ [  78.264, 173.394,-213.785], [ -63.451,   9.570, -13.504], [  -7.411, -31.806, -29.975], [   1.456,  -1.298,  43.628] ] ),
                (  99, [ [  19.915, 182.976,-226.748], [ -53.236,   9.593, -12.420], [ -16.164, -43.819, -39.572], [   3.936,   4.415,  60.578] ] ),
                ( 100, [ [ 108.125, 230.310, -93.528], [ -60.591,  28.132,   1.067], [ -13.931, -33.044, -12.007], [   4.596, -12.946,  57.694] ] ),
                ( 101, [ [ 141.968, 177.527,-116.980], [ -50.186,  17.132,  -3.637], [  13.607, -45.829, -52.443], [   9.726, -10.010,  55.555] ] ),
                ( 102, [ [  93.588, 196.028,-117.204], [ -45.032,  22.969,   1.979], [ -14.635, -34.317, -34.909], [   8.867, -14.370,  62.135] ] ),
                ( 103, [ [  51.988, 222.000,-112.818], [ -38.646,  27.605,   7.146], [ -48.820, -28.128, -36.207], [   1.902,  -5.868,  69.749] ] ),
                ( 104, [ [ 133.055, 155.147,-183.014], [ -50.398,   9.084,  23.516], [  -9.198, -33.156,  20.315], [  12.146,  20.340,  62.470] ] ),
                ( 105, [ [  80.120, 164.701,-163.613], [ -55.883,  12.515,  11.024], [   5.477, -51.428,  35.601], [  13.858,  40.686,  62.263] ] ),
                ( 106, [ [  22.839, 178.915,-159.278], [ -58.296,  15.237,  -1.951], [ -15.341, -83.376,  10.208], [  25.833,  57.492,  69.483] ] ),
                ( 107, [ [ 113.145,  88.046,-272.215], [  -9.350,  -8.812,  -7.836], [  13.248,  -7.657,  -1.039], [  -0.655,   3.227,  25.179] ] ),
                ( 108, [ [ 104.281,  79.019,-281.413], [  -8.550,  -9.748,  -9.808], [  22.001,   2.416,   6.648], [ -24.411,   8.417,  37.367] ] ),
                ( 109, [ [  94.545,  69.272,-290.944], [  -5.290, -12.555, -10.931], [  29.050,   6.005,  16.206], [ -30.419,   6.714,  45.052] ] ),
                ( 110, [ [ 126.005,  78.321,-270.120], [  -1.523,  -8.676,  -5.124], [  23.157,  -2.063,  11.752], [  14.097, -10.214,  23.163] ] ),
                ( 111, [ [ 110.779,  90.715,-234.929], [ -20.626,  -3.554,   0.196], [  25.282, -19.816,  -8.626], [  -7.947,   4.073,  48.804] ] ),
                ( 112, [ [  87.111,  87.499,-235.179], [ -25.846,  -6.325,  -4.452], [  34.640, -25.610,   1.322], [ -10.371,   9.996,  53.349] ] ),
                ( 113, [ [  60.277,  76.353,-242.362], [ -29.944,  -9.378, -10.303], [  59.866, -22.816,  11.218], [ -24.802,  14.723,  58.770] ] ),
                ( 114, [ [ 134.247,  72.870,-242.126], [  -0.240, -10.729,   1.399], [  58.738,   2.608,   2.449], [   2.955,  -0.060,  31.056] ] ),
                ( 115, [ [ 120.546, 104.051,-176.461], [ -29.942,  -7.987,   0.166], [   1.555, -50.037, -14.328], [  26.018,   2.579,  44.916] ] ),
                ( 116, [ [  86.014,  97.683,-176.263], [ -39.074,  -4.278,   0.229], [  32.793, -54.376, -24.428], [  11.840,   5.546,  52.991] ] ),
                ( 117, [ [  42.944,  92.427,-176.508], [ -47.600,   1.553,   0.293], [  65.220, -66.919, -26.935], [   7.926,   7.876,  68.237] ] ),
                ( 118, [ [ 131.420,  72.412,-210.024], [  -0.226, -19.168,   1.863], [  53.830,  -1.229, -27.395], [  -0.892,  -2.235,  47.338] ] ),
                ( 119, [ [ 113.261, 209.394, -39.187], [ -58.015,  16.946,   3.696], [  -7.484, -28.183, -16.526], [   8.287, -32.069,  37.482] ] ),
                ( 120, [ [ 153.826, 168.085, -64.529], [ -56.165,  -3.576,   3.482], [   4.306, -46.196, -41.906], [ -10.396,  -1.517,  54.580] ] ),
                ( 121, [ [ 102.498, 178.383, -60.042], [ -46.388,  16.516,   4.961], [ -14.042, -33.069, -26.005], [  16.284, -10.358,  47.201] ] ),
                ( 122, [ [  64.189, 202.493, -54.339], [ -28.979,  30.763,   5.289], [ -43.499, -29.808, -24.779], [  33.688, -24.925,  55.097] ] ),
                ( 123, [ [ 123.540, 138.168,-101.485], [ -34.841,  -0.986,   9.056], [  -6.374, -36.724, -42.407], [  -7.235, -12.421,  65.846] ] ),
                ( 124, [ [  85.021, 144.013, -90.695], [ -43.409,   7.248,   9.410], [  -0.486, -43.000, -39.187], [  23.904, -17.493,  55.956] ] ),
                ( 125, [ [  37.909, 155.369, -80.880], [ -50.463,  16.948,   9.056], [   0.425, -61.681, -36.664], [  34.615, -21.375,  78.761] ] ),
                ( 126, [ [ 141.972, 106.508,-140.392], [ -33.103,  -6.698,   5.569], [  24.138, -37.596, -25.097], [  17.276,  -3.374,  47.967] ] ),
                ( 127, [ [ 107.369,  98.979,-132.633], [ -36.268,  -6.344,  10.645], [  42.344, -39.637, -37.346], [  30.705,   0.158,  44.660] ] ),
                ( 128, [ [  70.498,  93.229,-117.454], [ -37.708,  -2.027,  19.331], [  68.482, -57.114, -32.013], [  58.629, -15.041,  46.951] ] ),
                ( 129, [ [ 163.749,  72.265,-157.014], [  -9.560, -38.172,  14.941], [  64.175, -23.344,   7.170], [ -12.951,  22.136,  66.011] ] ),
                ( 130, [ [ 121.257, 174.979, -23.486], [ -61.541,  10.197,  19.011], [  -3.863, -47.244,  -6.568], [  14.695,  -1.216,  28.119] ] ),
                ( 131, [ [ 117.057, 132.479, -55.764], [ -61.489,   4.844,  37.672], [  11.585, -39.956, -35.077], [  34.195,  -6.085,  11.996] ] ),
                ( 132, [ [ 146.435, 101.417, -88.480], [ -12.117, -28.636,  11.354], [  25.548, -28.865, -55.831], [  50.980,  -2.219,  39.663] ] )
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
            lowerLeftElementCount = (lElementsCount1 * (lElementsCount2-1) * lElementsCount3 + lElementsCount1)
            idx = lowerLeftElementCount + (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'apex of left lung')
            markerLocation.assignMeshLocation(cache, element1, [0.0, 0.0, 1.0])

            idx = lElementsCount1 * (lElementsCount2 - 1) + 1
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'ventral base of left lung')
            markerLocation.assignMeshLocation(cache, element1, [0.0, 1.0, 0.0])

            idx = lElementsCount1 * (lElementsCount2 // 2)
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'medial base of left lung')
            markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 0.0])

            upperLeftElementCount = (uElementsCount1 * uElementsCount2 * (uElementsCount3-1))
            leftLungElementCount = lowerLeftElementCount + upperLeftElementCount

            lowerRightElementCount = lowerLeftElementCount

            idx = leftLungElementCount + lowerRightElementCount + \
                  (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)

            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'apex of right lung')
            markerLocation.assignMeshLocation(cache, element1, [0.0, 0.0, 1.0])

            idx = leftLungElementCount + lElementsCount1 * lElementsCount2
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'ventral base of right lung')
            markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 0.0])

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * (lElementsCount3 - 1)) + lElementsCount1
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'lateral side of right lung')
            markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 1.0])

        elif isMouse:
            # valueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ]
            nodeFieldParameters = [
                (   1, [ [  -0.06500,  -9.79600,   3.32500], [   0.93700,   0.96700,   2.47500], [  -1.72900,  -0.55700,   0.27900], [  -1.97500,  -3.88300,   3.53100] ] ),
                (   2, [ [  -0.00500,  -9.22800,   5.76600], [  -0.74000,   0.15300,   2.17800], [  -1.51700,   0.36800,  -0.21400], [  -0.71900,  -2.10500,   2.41600] ] ),
                (   3, [ [  -1.87500, -11.05600,   1.28600], [  -0.21900,   0.91700,   2.49900], [  -2.39900,  -0.36900,  -1.06700], [  -2.61300,  -2.51900,   4.94800] ] ),
                (   4, [ [  -1.93400, -10.00200,   3.65200], [   0.10200,   1.18600,   2.22000], [  -1.96200,   0.15900,   0.36800], [  -1.82400,  -2.13500,   2.94900] ] ),
                (   5, [ [  -1.68200,  -8.70900,   5.70400], [   0.40000,   1.39200,   1.87100], [  -1.82200,   0.66700,   0.09200], [  -0.58800,  -0.21200,   2.94500] ] ),
                (   6, [ [  -3.99400, -10.67100,   1.23900], [  -0.01000,   0.98000,   3.19800], [  -2.14700,   0.89400,   0.42200], [  -1.85400,  -1.87000,   4.95800] ] ),
                (   7, [ [  -3.89800,  -9.45400,   4.04500], [   0.20100,   1.44300,   2.38200], [  -1.94000,   0.83500,   0.33100], [  -0.94500,  -0.58800,   3.06100] ] ),
                (   8, [ [  -3.61900,  -7.88600,   5.97900], [   0.35100,   1.66300,   1.45900], [  -1.91700,   1.03900,   0.16000], [  -0.23900,   0.30600,   2.67600] ] ),
                (   9, [ [  -5.86600,  -9.33300,   2.12700], [   0.01700,   0.59900,   2.31000], [  -1.70800,   1.80100,   1.54800], [  -1.04600,  -0.54700,   4.42900] ] ),
                (  10, [ [  -5.76000,  -8.34000,   4.30200], [   0.19500,   1.37500,   1.98800], [  -1.60500,   1.16300,   0.13500], [  -0.13300,   0.32500,   2.99100] ] ),
                (  11, [ [  -5.48200,  -6.63600,   6.01500], [   0.35200,   1.99000,   1.40800], [  -1.65900,   1.25900,  -0.05100], [   0.25300,   0.26400,   2.50200] ] ),
                (  12, [ [  -7.10400,  -7.17500,   4.33400], [  -0.44700,   2.06800,   1.94900], [  -1.07500,   1.15600,  -0.07000], [  -0.35600,   0.07900,   2.65000] ] ),
                (  13, [ [  -6.93700,  -5.39200,   5.89100], [   0.73400,   1.40700,   1.09300], [  -1.24500,   1.22400,  -0.19700], [   0.68600,   1.13000,   2.17900] ] ),
                (  14, [ [  -1.55200, -12.41900,   7.12200], [   1.89000,   1.23100,   1.52300], [  -1.81000,   0.92100,   0.01900], [  -0.91700,  -1.20200,   3.91600] ] ),
                (  15, [ [  -0.54600, -10.44300,   8.76100], [   0.11000,   2.46200,   1.58800], [  -1.25400,   2.16800,   0.00900], [  -0.30400,  -0.15400,   3.37700] ] ),
                (  16, [ [  -3.84500, -12.69000,   6.04000], [   0.31200,   0.90800,   0.87300], [  -2.16400,   0.54100,  -0.59300], [  -1.28400,  -0.70700,   4.47700] ] ),
                (  17, [ [  -3.22600, -11.20300,   7.18000], [   0.92300,   2.05600,   1.39900], [  -1.50800,   1.49400,   0.09700], [  -0.63100,  -0.11700,   3.89800] ] ),
                (  18, [ [  -1.95900,  -8.56200,   8.76100], [   1.60800,   3.22200,   1.76000], [  -1.55900,   1.57300,  -0.00900], [   0.04300,   0.51000,   3.11900] ] ),
                (  19, [ [  -5.48600, -11.63400,   5.93600], [   1.00200,   2.17400,   1.35300], [  -1.53700,   1.70600,   0.08800], [  -1.09700,  -0.02200,   4.34400] ] ),
                (  20, [ [  -4.51300,  -9.46500,   7.31500], [   0.94300,   2.16400,   1.40500], [  -1.27500,   1.72900,   0.06700], [  -0.25700,   0.58500,   3.38600] ] ),
                (  21, [ [  -3.60000,  -7.30600,   8.74500], [   0.88400,   2.15200,   1.45500], [  -1.54100,   1.29600,  -0.07900], [   0.27900,   0.85000,   2.82400] ] ),
                (  22, [ [  -6.64300,  -9.35400,   6.25900], [   0.93200,   1.51100,   0.93700], [  -0.87500,   2.49400,   0.58500], [  -0.49800,   0.51200,   3.78400] ] ),
                (  23, [ [  -5.77100,  -7.75200,   7.31600], [   0.81100,   1.68900,   1.17300], [  -1.37400,   1.36400,  -0.12000], [   0.11200,   0.84800,   3.01900] ] ),
                (  24, [ [  -5.03400,  -5.98300,   8.60700], [   0.66100,   1.84700,   1.40700], [  -1.28000,   1.58200,  -0.26400], [   0.64000,   1.03900,   2.64400] ] ),
                (  25, [ [  -7.17300,  -6.74200,   7.10200], [   0.25600,   2.72800,   1.02500], [  -1.39400,   0.63900,  -0.30000], [   0.22600,   0.78500,   2.83700] ] ),
                (  26, [ [  -6.09900,  -4.17500,   8.21500], [   1.80500,   2.29700,   1.14800], [  -0.84000,   2.00500,  -0.51300], [   0.98900,   1.30200,   2.46800] ] ),
                (  27, [ [  -1.93100, -12.43600,  10.55900], [   2.49000,   1.68500,   1.18300], [  -1.15200,   2.23900,   0.69100], [  -0.54300,   0.32500,   3.55900] ] ),
                (  28, [ [  -0.56200,  -9.47900,  12.04200], [   0.22400,   3.80500,   1.60300], [  -0.78100,   1.89400,  -0.16000], [  -0.44600,   0.30600,   3.16500] ] ),
                (  29, [ [  -4.52000, -12.65600,   9.95600], [   1.30600,   2.32000,   0.98900], [  -2.41100,   0.92900,  -0.48400], [  -0.46700,   0.65000,   3.64600] ] ),
                (  30, [ [  -3.13300, -10.25400,  10.92300], [   1.46800,   2.48500,   0.94400], [  -1.24000,   2.10000,   0.02900], [  -0.18500,   0.74000,   3.42200] ] ),
                (  31, [ [  -1.58200,  -7.68800,  11.83900], [   1.63300,   2.64600,   0.88900], [  -1.25200,   1.66600,  -0.24300], [  -0.06200,   0.20100,   3.11400] ] ),
                (  32, [ [  -6.19700, -10.92400,   9.69400], [   2.06600,   2.90800,   0.94500], [  -1.25600,   2.20900,  -0.19200], [  -0.23700,   1.12300,   3.54400] ] ),
                (  33, [ [  -4.37800,  -8.28200,  10.63900], [   1.57100,   2.37300,   0.94500], [  -1.20500,   1.80600,  -0.31200], [   0.08200,   0.79000,   3.19500] ] ),
                (  34, [ [  -3.04500,  -6.18300,  11.56000], [   1.09400,   1.82300,   0.89500], [  -1.31100,   1.58200,  -0.31600], [   0.41300,   0.44900,   2.91100] ] ),
                (  35, [ [  -6.88900,  -8.42400,   9.59400], [   1.35800,   1.61900,   0.61600], [  -0.24700,   2.73500,   0.10200], [   0.26400,   1.11000,   3.22300] ] ),
                (  36, [ [  -5.53200,  -6.64300,  10.30600], [   1.35200,   1.94000,   0.80700], [  -1.15000,   1.33500,  -0.37700], [   0.46000,   0.90100,   2.97900] ] ),
                (  37, [ [  -4.19800,  -4.54300,  11.21100], [   1.31400,   2.25600,   1.00300], [  -0.95600,   1.70800,  -0.37900], [   0.81800,   0.65100,   2.72300] ] ),
                (  38, [ [  -6.63100,  -5.59800,   9.91100], [   0.94700,   2.93800,   0.60800], [  -1.03500,   0.74500,  -0.40700], [   0.98500,   1.06600,   2.93400] ] ),
                (  39, [ [  -4.95200,  -2.79500,  10.80900], [   2.34400,   2.59400,   1.15400], [  -0.54700,   1.77500,  -0.42200], [   1.20400,   0.55900,   2.59700] ] ),
                (  40, [ [  -2.63300, -11.75400,  14.13900], [   2.07100,   1.01000,   1.01300], [  -0.94500,   2.13200,  -0.11000], [  -1.19500,   1.29000,   2.62900] ] ),
                (  41, [ [  -1.34700,  -9.76300,  14.88200], [   0.43300,   2.56600,   0.40900], [  -0.67800,   1.57400,  -0.00900], [  -0.85400,   0.07400,   2.62100] ] ),
                (  42, [ [  -4.78400, -11.50200,  13.19700], [   1.12400,   1.88500,   0.74200], [  -1.93100,   1.31700,  -0.71500], [   0.33300,   1.24300,   2.91700] ] ),
                (  43, [ [  -3.53400,  -9.71000,  13.96000], [   1.37200,   1.69300,   0.78200], [  -0.85500,   1.95400,  -0.24900], [  -0.46900,  -0.15300,   2.45300] ] ),
                (  44, [ [  -2.04800,  -8.12700,  14.75500], [   1.59400,   1.46900,   0.80600], [  -0.72300,   1.69400,  -0.24600], [  -0.31400,  -0.89000,   2.71200] ] ),
                (  45, [ [  -5.99600,  -9.47500,  12.89200], [   1.70200,   1.70700,   0.76600], [  -0.70000,   2.26100,  -0.31000], [   0.83600,   1.57400,   3.01000] ] ),
                (  46, [ [  -4.34300,  -7.84900,  13.64900], [   1.60500,   1.54500,   0.74600], [  -0.66600,   1.88500,  -0.36200], [  -0.13900,   0.17400,   2.48300] ] ),
                (  47, [ [  -2.78800,  -6.38400,  14.38400], [   1.50500,   1.38500,   0.72400], [  -0.70800,   1.74800,  -0.46700], [   0.53500,  -0.43200,   3.00600] ] ),
                (  48, [ [  -6.14400,  -7.17700,  12.60200], [   1.21900,   1.17800,   0.66400], [   0.37800,   2.50500,  -0.03200], [   1.15600,   1.37100,   2.63400] ] ),
                (  49, [ [  -4.86400,  -5.95300,  13.24000], [   1.34200,   1.27000,   0.61200], [  -0.41900,   1.59500,  -0.39800], [   0.39900,   0.41100,   2.26200] ] ),
                (  50, [ [  -3.46100,  -4.63700,  13.82300], [   1.46400,   1.36000,   0.55300], [  -0.50700,   1.71100,  -0.60900], [   0.85100,  -0.09300,   2.33400] ] ),
                (  51, [ [  -5.18900,  -4.65900,  12.87300], [   1.22600,   2.11800,   0.29700], [  -0.23000,   0.99100,  -0.33700], [   1.02700,  -0.31500,   2.68800] ] ),
                (  52, [ [  -3.80900,  -2.97900,  13.17500], [   1.49600,   1.21300,   0.30100], [  -0.18800,   1.59800,  -0.68200], [   1.26300,  -0.91100,   2.64300] ] ),
                (  53, [ [  -3.97000, -10.25800,  15.67100], [   1.39000,   1.09500,   2.04300], [  -1.07900,   2.17800,   0.78400], [  -0.36000,  -0.84200,   0.86500] ] ),
                (  54, [ [  -2.21000,  -9.38800,  17.16100], [   2.01800,   0.61000,   0.88900], [  -0.38600,   1.79200,   1.71400], [  -0.00900,  -1.59600,   2.05300] ] ),
                (  55, [ [  -4.58000,  -7.85600,  15.56900], [   2.11100,   1.24600,   2.39500], [  -0.34700,   2.32300,  -0.47500], [  -0.32900,  -0.18600,   1.33500] ] ),
                (  56, [ [  -1.94700,  -7.07000,  17.50600], [   3.03100,   0.31400,   1.42200], [  -0.14400,   2.64400,  -0.63500], [   1.14000,  -0.93300,   3.21700] ] ),
                (  57, [ [  -4.66700,  -5.75200,  14.78800], [   1.90000,   1.26100,   1.65700], [  -0.32900,   1.75000,  -1.46200], [  -0.00600,  -0.00900,   0.81900] ] ),
                (  58, [ [  -2.54100,  -4.72700,  15.86200], [   2.25100,   0.75500,   0.47100], [  -0.93900,   2.14800,  -2.20600], [   0.98300,  -0.08500,   1.73600] ] ),
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
                ( 127, [ [   3.15700,  -8.31000,   9.39800], [   2.53800,  -0.40100,   0.43300], [  -1.64700,   2.16100,   0.67400], [  -1.60100,   0.06100,   1.77100] ] ),
                ( 128, [ [   4.11100,  -5.52400,  10.77500], [  -1.93100,   3.85000,   2.14600], [  -1.36100,   1.34600,   0.44700], [  -3.02000,  -1.80300,   0.23600] ] ),
                ( 129, [ [   0.06500, -11.43700,   7.21500], [   0.18400,   2.21200,  -0.19300], [  -0.54200,  -0.39700,  -0.29100], [   0.10900,   0.04700,   2.24600] ] ),
                ( 130, [ [  -0.16200,  -9.48800,   6.61900], [   1.30200,   2.93700,   3.92300], [  -1.02400,   0.49500,  -1.87400], [  -1.35400,  -0.03500,   3.61400] ] ),
                ( 131, [ [   1.19500,  -6.75800,   9.68000], [   1.40000,   2.49800,   2.16700], [  -2.14900,   0.85000,  -0.15400], [  -1.51700,  -1.03900,   0.74200] ] ),
                ( 132, [ [   2.49300,  -4.64600,  11.05500], [   1.15800,   1.67200,   0.56600], [  -1.65900,   0.66700,   0.08800], [  -2.05100,  -1.69200,   0.35800] ] ),
                ( 133, [ [  -2.02800,  -8.78200,   6.40500], [   0.97200,   1.92400,   3.23400], [  -2.03500,   1.31500,   0.66900], [  -0.36500,   1.12500,   4.45200] ] ),
                ( 134, [ [  -0.79200,  -6.61800,   9.20000], [   1.48600,   2.37600,   2.30700], [  -1.85100,   0.70100,  -0.13300], [  -0.68500,  -0.48100,   1.57900] ] ),
                ( 135, [ [   0.87700,  -4.12500,  10.97600], [   1.81900,   2.56400,   1.22300], [  -1.73500,   0.43000,  -0.13400], [  -1.77200,  -1.55600,   0.42900] ] ),
                ( 136, [ [  -3.58000,  -7.09200,   8.05300], [   1.29300,   1.54400,   1.32600], [  -1.32900,   1.77000,   1.27100], [   0.40200,   0.69000,   2.47200] ] ),
                ( 137, [ [  -2.27700,  -5.49600,   9.40000], [   1.31400,   1.64800,   1.36800], [  -1.33600,   1.34600,   0.31300], [  -0.46000,  -0.52000,   1.27200] ] ),
                ( 138, [ [  -0.95000,  -3.76000,  10.74000], [   1.33200,   1.75200,   1.40900], [  -1.63800,   0.73000,  -0.22100], [  -1.09700,  -1.31800,   0.38500] ] ),
                ( 139, [ [  -4.66400,  -5.35100,   8.96000], [   1.32600,   1.46500,   0.91600], [  -0.79700,   1.90500,   0.71200], [   0.31300,   0.30300,   1.46900] ] ),
                ( 140, [ [  -3.40400,  -3.97400,   9.81500], [   1.19400,   1.28900,   0.79400], [  -0.86800,   1.93100,   0.39500], [  -0.44800,  -0.51400,   0.82400] ] ),
                ( 141, [ [  -2.29900,  -2.77300,  10.58400], [   1.05900,   1.11500,   0.67500], [  -0.92400,   1.77400,   0.34600], [  -1.01000,  -0.90200,   0.52200] ] ),
                ( 142, [ [  -5.16100,  -3.36500,   9.46100], [   1.02300,   1.84500,   0.77500], [  -0.19200,   2.02200,   0.28300], [  -0.02400,   0.01700,   0.98200] ] ),
                ( 143, [ [  -3.88800,  -1.70800,  10.15300], [   1.51100,   1.44000,   0.59700], [  -0.09800,   2.54100,   0.27400], [  -0.42400,  -0.07600,   0.80100] ] ),
                ( 144, [ [  -2.21800,  -0.53800,  10.63300], [   2.02900,   0.25900,   0.44900], [  -1.61300,   2.07400,  -0.10700], [  -1.05300,  -0.66900,   1.14800] ] ),
                ( 145, [ [   1.75100, -10.58100,  10.76100], [  -0.97700,   1.55100,  -0.18500], [  -2.85000,  -1.04200,  -1.37200], [  -1.33400,  -0.62600,   2.63200] ] ),
                ( 146, [ [   1.22000,  -9.13500,  10.80000], [  -0.03300,   1.25800,   0.27400], [  -2.33700,  -0.38600,  -0.53500], [  -0.52300,  -0.41900,   2.19300] ] ),
                ( 147, [ [   1.51600,  -8.25800,  11.17000], [   0.26000,   0.83600,   0.41900], [  -1.73100,   0.31600,  -0.13700], [  -1.20200,  -0.48000,   2.11000] ] ),
                ( 148, [ [   1.73700,  -7.47700,  11.62600], [   0.18100,   0.72300,   0.49100], [  -0.88900,   1.03000,   0.36900], [  -1.66600,  -1.94500,   1.42400] ] ),
                ( 149, [ [  -0.11300, -11.26300,   9.86400], [  -0.57200,   2.27600,   0.46100], [  -0.87900,  -0.32100,  -0.42300], [  -0.46500,   0.30100,   3.03200] ] ),
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
                ( 162, [ [  -4.91500,  -3.17200,  10.67900], [   0.58300,   1.47100,   0.33300], [  -0.45900,   1.81800,   0.04300], [   0.51700,   0.36900,   1.41700] ] ),
                ( 163, [ [  -4.16900,  -1.97000,  11.09800], [   0.89300,   0.89100,   0.49400], [  -0.08600,   2.67300,   0.16500], [  -0.11900,  -0.44500,   1.05400] ] ),
                ( 164, [ [  -3.22900,  -1.40800,  11.61400], [   0.94000,   0.22100,   0.51200], [   0.56900,   2.37200,   0.43200], [  -0.95100,  -1.06000,   0.79300] ] ),
                ( 165, [ [  -0.21700, -10.72500,  13.59000], [   1.15600,   1.58500,  -0.23100], [  -0.98100,  -0.11900,  -0.49200], [  -2.54900,   0.33100,   2.96300] ] ),
                ( 166, [ [   0.88000,  -9.22100,  13.37100], [   1.03800,   1.42300,  -0.20700], [  -1.37100,   0.43200,   0.23800], [  -0.06900,  -1.36500,   2.16900] ] ),
                ( 167, [ [  -0.93000, -10.81200,  13.23200], [   0.43900,   2.22800,   0.14800], [  -0.44500,  -0.05400,  -0.22300], [  -1.16400,   0.59900,   3.69000] ] ),
                ( 168, [ [  -0.48300,  -8.53800,  13.38300], [   0.45600,   2.31900,   0.15400], [  -1.30900,   0.92000,  -0.22200], [   0.11300,  -0.51400,   2.64700] ] ),
                ( 169, [ [  -1.65200,  -7.41300,  12.91600], [   0.63100,   0.43300,   0.19800], [  -1.20000,   1.28300,  -0.38200], [  -0.15200,  -0.30000,   2.10300] ] ),
                ( 170, [ [  -2.86200,  -5.97800,  12.63700], [   0.71200,   0.49800,   0.26600], [  -1.00400,   1.61400,  -0.19300], [  -0.09600,   0.07200,   1.90500] ] ),
                ( 171, [ [  -3.63100,  -4.22700,  12.53600], [   1.04500,   0.89800,   0.69600], [  -0.60500,   1.70400,  -0.23300], [   0.30200,   0.35300,   1.71400] ] ),
                ( 172, [ [  -4.07700,  -2.59900,  12.18800], [   0.79000,   0.79000,   0.97500], [  -0.28400,   1.54000,  -0.45900], [   0.29600,  -0.79400,   1.10100] ] )
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
            diaNodeIndex, nodeIdentifier = getDiaphragmaticLungNodes(cache, coordinates, generateParameters,
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

            # Apex annotation point
            idx = elementsCount1 * elementsCount2 * (elementsCount3 - 1) + elementsCount1 * (elementsCount2 // 2)
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'apex of left lung')
            markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 1.0])

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
