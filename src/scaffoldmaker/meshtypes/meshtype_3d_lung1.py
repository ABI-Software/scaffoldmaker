'''
Generates 3D lung surface mesh.
'''

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
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
        isMouse = 'Mouse' in parameterSetName
        isHuman = 'Human' in parameterSetName

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

        if isHuman:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(upperLeftLungGroup)
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
        eftWedgeCollapseXi1_26 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
        eftWedgeCollapseXi1_37 = eftfactory.createEftWedgeCollapseXi1Quadrant([3, 7])
        eftWedgeCollapseXi1_48 = eftfactory.createEftWedgeCollapseXi1Quadrant([4, 8])
        eftWedgeCollapseXi1_57 = eftfactory.createEftWedgeCollapseXi1Quadrant([5, 7])
        eftWedgeCollapseXi1_68 = eftfactory.createEftWedgeCollapseXi1Quadrant([6, 8])
        eftWedgeCollapseXi2_56 = eftfactory.createEftWedgeCollapseXi2Quadrant([5, 6])
        eftWedgeCollapseXi2_78 = eftfactory.createEftWedgeCollapseXi2Quadrant([7, 8])
        eftTetCollapseXi1Xi2_71 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(7, 1)
        eftTetCollapseXi1Xi2_82 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(8, 2)
        eftTetCollapseXi1Xi2_63 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(6, 3)
        eftTetCollapseXi1Xi2_53 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(5, 3)

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

            generateParameters = False

            lElementsCount1 = 2
            lElementsCount2 = 4
            lElementsCount3 = 3

            uElementsCount1 = 2
            uElementsCount2 = 4
            uElementsCount3 = 4

            # Create nodes
            nodeIdentifier = 1
            lowerLeftNodeIds = []
            upperLeftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []
            d1 = [1.0, 0.0, 0.0]
            d2 = [0.0, 1.0, 0.0]
            d3 = [0.0, 0.0, 1.0]
            nodeIndex = 0
            leftLung = 0
            rightLung = 1

            for lung in (leftLung, rightLung):
                lowerNodeIds = lowerLeftNodeIds if (lung == leftLung) else lowerRightNodeIds
                xMirror = 0 if (lung == leftLung) else 150 # Offset

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
                upperNodeIds = upperLeftNodeIds if (lung == leftLung) else upperRightNodeIds

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

            # Create elements
            elementIdentifier = 1
            for lung in (leftLung, rightLung):
                lowerNodeIds = lowerLeftNodeIds if (lung == leftLung) else lowerRightNodeIds

                # Lower lobe elements
                for e3 in range(lElementsCount3):
                    for e2 in range(lElementsCount2):
                        for e1 in range(lElementsCount1):
                            eft = eftRegular
                            nodeIdentifiers = [
                                lowerNodeIds[e3    ][e2][e1], lowerNodeIds[e3    ][e2][e1 + 1], lowerNodeIds[e3    ][e2 + 1][e1], lowerNodeIds[e3    ][e2 + 1][e1 + 1],
                                lowerNodeIds[e3 + 1][e2][e1], lowerNodeIds[e3 + 1][e2][e1 + 1], lowerNodeIds[e3 + 1][e2 + 1][e1], lowerNodeIds[e3 + 1][e2 + 1][e1 + 1]]

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
                                nodeIdentifiers[2] = lowerNodeIds[e3 - 1][e2 + 1][e1    ]
                                nodeIdentifiers[3] = lowerNodeIds[e3 - 1][e2 + 1][e1 + 1]
                                nodeIdentifiers[6] = lowerNodeIds[e3 - 1][e2 + 2][e1    ]
                                nodeIdentifiers[7] = lowerNodeIds[e3 - 1][e2 + 2][e1 + 1]
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
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
                            if lung == leftLung:
                                lowerLeftLungMeshGroup.addElement(element)
                                leftLungMeshGroup.addElement(element)
                            else:
                                rightLungMeshGroup.addElement(element)
                                lowerRightLungMeshGroup.addElement(element)

                # Upper lobe elements
                upperNodeIds = upperLeftNodeIds if (lung == leftLung) else upperRightNodeIds
                for e3 in range(uElementsCount3):
                    for e2 in range(uElementsCount2):
                        for e1 in range(uElementsCount1):
                            eft = eftRegular
                            nodeIdentifiers = [
                                upperNodeIds[e3    ][e2][e1], upperNodeIds[e3    ][e2][e1 + 1], upperNodeIds[e3    ][e2 + 1][e1], upperNodeIds[e3    ][e2 + 1][e1 + 1],
                                upperNodeIds[e3 + 1][e2][e1], upperNodeIds[e3 + 1][e2][e1 + 1], upperNodeIds[e3 + 1][e2 + 1][e1], upperNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                            if (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                                # Distal-front wedge elements / lung edge (inside)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(2)
                                # eft = eftWedgeCollapseXi1_37
                                # ------------------------------
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                nodes = [3, 4, 7, 8]
                                collapseNodes = [3, 7]
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                                # zero cross derivative parameters
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D3_DS1DS2DS3, [])
                                remapEftLocalNodes(eft, 6, ln_map)

                            elif (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == (uElementsCount1 - 1)):
                                # Distal-back wedge elements / lung edge (outside)
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(3)
                                # eft = eftWedgeCollapseXi1_48
                                # ------------------------------
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                nodes = [3, 4, 7, 8]
                                collapseNodes = [4, 8]
                                remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                                ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                                # zero cross derivative parameters
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D3_DS1DS2DS3, [])
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
                            elif (e3 == (uElementsCount3 - 1)) and (0 < e2 < (uElementsCount2 - 1)) and (e1 == (uElementsCount1 - 1)):
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
                                # Top-front-distal tetrahedron wedge elements / lung edge (inside)
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(2)
                                #eft = eftTetCollapseXi1Xi2_63
                                # ------------------------------
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                nodes = [5, 6, 7, 8]

                                # remap parameters on xi3 = 1 before collapsing nodes
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                                remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft, [3], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])

                                ln_map = [1, 2, 3, 3, 4, 4, 4, 4]

                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS2DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D3_DS1DS2DS3, [])
                                remapEftLocalNodes(eft, 4, ln_map)

                            elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == (uElementsCount1 - 1)):
                                # Top-front-distal tetrahedron wedge elements / lung edge (outside)
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(3)
                                # eft = eftTetCollapseXi1Xi2_53
                                # ------------------------------
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                nodes = [5, 6, 7, 8]

                                # remap parameters on xi3 = 1 before collapsing nodes
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                                remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])

                                ln_map = [1, 2, 3, 3, 4, 4, 4, 4]

                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS2, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS1DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D2_DS2DS3, [])
                                remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D3_DS1DS2DS3, [])
                                remapEftLocalNodes(eft, 4, ln_map)


                            elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 3)):
                                # Remapped cube element 1
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 2)):
                                # Remapped cube element 2
                                eft = eftfactory.createEftBasic()
                                setEftScaleFactorIds(eft, [1], [])
                                remapEftNodeValueLabel(eft, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
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
                            if lung == leftLung:
                                leftLungMeshGroup.addElement(element)
                                upperLeftLungMeshGroup.addElement(element)
                            else:
                                rightLungMeshGroup.addElement(element)
                                if e3 < (uElementsCount3 - 2):
                                    middleRightLungMeshGroup.addElement(element)
                                else:
                                    upperRightLungMeshGroup.addElement(element)

        elif isMouse:
            # valueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ]
            nodeFieldParameters = [
                ( 1, [ [ -0.06510, -9.79594,  3.32525], [  0.93722,  0.96743,  2.47497], [ -1.72856, -0.55693,  0.27902], [ -1.97460, -3.88256,  3.53085] ] ),
                ( 2, [ [ -0.00541, -9.22791,  5.76568], [ -0.74021,  0.15263,  2.17755], [ -1.51717,  0.36826, -0.21360], [ -0.71889, -2.10495,  2.41646] ] ),
                ( 3, [ [ -1.87483,-11.05622,  1.28598], [ -0.21863,  0.91720,  2.49890], [ -2.39917, -0.36948, -1.06683], [ -2.61253, -2.51865,  4.94809] ] ),
                ( 4, [ [ -1.93372,-10.00227,  3.65228], [  0.10207,  1.18560,  2.21981], [ -1.96249,  0.15916,  0.36760], [ -1.82422, -2.13455,  2.94947] ] ),
                ( 5, [ [ -1.68151, -8.70879,  5.70416], [  0.39968,  1.39208,  1.87146], [ -1.82198,  0.66682,  0.09240], [ -0.58759, -0.21221,  2.94457] ] ),
                ( 6, [ [ -3.99360,-10.67083,  1.23895], [ -0.00974,  0.98009,  3.19786], [ -2.14704,  0.89450,  0.42247], [ -1.85405, -1.87022,  4.95808] ] ),
                ( 7, [ [ -3.89795, -9.45423,  4.04517], [  0.20114,  1.44312,  2.38197], [ -1.93969,  0.83545,  0.33089], [ -0.94518, -0.58799,  3.06084] ] ),
                ( 8, [ [ -3.61875, -7.88583,  5.97946], [  0.35072,  1.66260,  1.45935], [ -1.91709,  1.03902,  0.16009], [ -0.23855,  0.30601,  2.67575] ] ),
                ( 9, [ [ -5.86561, -9.33331,  2.12747], [  0.01662,  0.59905,  2.31040], [ -1.70810,  1.80053,  1.54826], [ -1.04603, -0.54686,  4.42926] ] ),
                ( 10, [ [ -5.75955, -8.33957,  4.30205], [  0.19514,  1.37528,  1.98808], [ -1.60492,  1.16263,  0.13548], [ -0.13285,  0.32516,  2.99109] ] ),
                ( 11, [ [ -5.48232, -6.63587,  6.01506], [  0.35197,  1.99049,  1.40847], [ -1.65938,  1.25924, -0.05112], [  0.25309,  0.26378,  2.50204] ] ),
                ( 12, [ [ -7.10444, -7.17483,  4.33436], [ -0.44741,  2.06805,  1.94872], [ -1.07465,  1.15589, -0.07020], [ -0.35618,  0.07890,  2.64986] ] ),
                ( 13, [ [ -6.93708, -5.39157,  5.89076], [  0.73437,  1.40695,  1.09300], [ -1.24473,  1.22404, -0.19662], [  0.68595,  1.12974,  2.17928] ] ),
                ( 14, [ [ -1.55172,-12.41852,  7.12162], [  1.88969,  1.23070,  1.52327], [ -1.81005,  0.92137,  0.01914], [ -0.91712, -1.20231,  3.91613] ] ),
                ( 15, [ [ -0.54610,-10.44301,  8.76054], [  0.11002,  2.46228,  1.58814], [ -1.25406,  2.16768,  0.00905], [ -0.30398, -0.15393,  3.37660] ] ),
                ( 16, [ [ -3.84472,-12.68973,  6.03981], [  0.31220,  0.90844,  0.87326], [ -2.16380,  0.54075, -0.59336], [ -1.28386, -0.70655,  4.47742] ] ),
                ( 17, [ [ -3.22579,-11.20330,  7.18002], [  0.92288,  2.05636,  1.39941], [ -1.50774,  1.49362,  0.09734], [ -0.63137, -0.11709,  3.89816] ] ),
                ( 18, [ [ -1.95902, -8.56165,  8.76085], [  1.60817,  3.22197,  1.75954], [ -1.55877,  1.57254, -0.00853], [  0.04254,  0.51009,  3.11885] ] ),
                ( 19, [ [ -5.48617,-11.63414,  5.93599], [  1.00208,  2.17448,  1.35335], [ -1.53695,  1.70630,  0.08816], [ -1.09664, -0.02164,  4.34388] ] ),
                ( 20, [ [ -4.51347, -9.46458,  7.31507], [  0.94319,  2.16434,  1.40462], [ -1.27505,  1.72935,  0.06740], [ -0.25714,  0.58514,  3.38599] ] ),
                ( 21, [ [ -3.60002, -7.30604,  8.74482], [  0.88359,  2.15244,  1.45468], [ -1.54078,  1.29572, -0.07879], [  0.27874,  0.85007,  2.82435] ] ),
                ( 22, [ [ -6.64349, -9.35382,  6.25941], [  0.93225,  1.51110,  0.93725], [ -0.87547,  2.49448,  0.58513], [ -0.49781,  0.51207,  3.78420] ] ),
                ( 23, [ [ -5.77052, -7.75206,  7.31586], [  0.81146,  1.68880,  1.17340], [ -1.37368,  1.36400, -0.12041], [  0.11169,  0.84799,  3.01929] ] ),
                ( 24, [ [ -5.03360, -5.98263,  8.60719], [  0.66118,  1.84670,  1.40670], [ -1.27950,  1.58159, -0.26354], [  0.64045,  1.03863,  2.64367] ] ),
                ( 25, [ [ -7.17256, -6.74200,  7.10165], [  0.25582,  2.72824,  1.02490], [ -1.39392,  0.63939, -0.30016], [  0.22637,  0.78534,  2.83697] ] ),
                ( 26, [ [ -6.09924, -4.17507,  8.21541], [  1.80545,  2.29701,  1.14832], [ -0.83983,  2.00500, -0.51273], [  0.98915,  1.30231,  2.46819] ] ),
                ( 27, [ [ -1.93094,-12.43564, 10.55902], [  2.48980,  1.68514,  1.18313], [ -1.15156,  2.23902,  0.69098], [ -0.54345,  0.32541,  3.55936] ] ),
                ( 28, [ [ -0.56152, -9.47868, 12.04160], [  0.22408,  3.80497,  1.60343], [ -0.78094,  1.89408, -0.16015], [ -0.44578,  0.30576,  3.16546] ] ),
                ( 29, [ [ -4.52036,-12.65650,  9.95617], [  1.30558,  2.31961,  0.98943], [ -2.41122,  0.92866, -0.48408], [ -0.46669,  0.65040,  3.64612] ] ),
                ( 30, [ [ -3.13315,-10.25385, 10.92288], [  1.46839,  2.48489,  0.94365], [ -1.24000,  2.09955,  0.02903], [ -0.18490,  0.73991,  3.42247] ] ),
                ( 31, [ [ -1.58226, -7.68819, 11.83927], [  1.63292,  2.64568,  0.88888], [ -1.25193,  1.66600, -0.24275], [ -0.06247,  0.20128,  3.11434] ] ),
                ( 32, [ [ -6.19695,-10.92411,  9.69362], [  2.06583,  2.90846,  0.94515], [ -1.25610,  2.20894, -0.19218], [ -0.23740,  1.12349,  3.54366] ] ),
                ( 33, [ [ -4.37803, -8.28227, 10.63922], [  1.57052,  2.37312,  0.94537], [ -1.20473,  1.80608, -0.31249], [  0.08243,  0.78962,  3.19547] ] ),
                ( 34, [ [ -3.04476, -6.18307, 11.56019], [  1.09442,  1.82261,  0.89526], [ -1.31144,  1.58163, -0.31638], [  0.41299,  0.44880,  2.91121] ] ),
                ( 35, [ [ -6.88857, -8.42366,  9.59418], [  1.35762,  1.61893,  0.61571], [ -0.24666,  2.73500,  0.10203], [  0.26350,  1.11037,  3.22349] ] ),
                ( 36, [ [ -5.53243, -6.64271, 10.30584], [  1.35243,  1.94031,  0.80660], [ -1.15016,  1.33550, -0.37674], [  0.45979,  0.90102,  2.97864] ] ),
                ( 37, [ [ -4.19832, -4.54309, 11.21099], [  1.31429,  2.25634,  1.00256], [ -0.95618,  1.70784, -0.37890], [  0.81785,  0.65091,  2.72280] ] ),
                ( 38, [ [ -6.63094, -5.59808,  9.91147], [  0.94728,  2.93781,  0.60812], [ -1.03451,  0.74487, -0.40713], [  0.98453,  1.06649,  2.93435] ] ),
                ( 39, [ [ -4.95196, -2.79494, 10.80914], [  2.34379,  2.59443,  1.15428], [ -0.54699,  1.77512, -0.42163], [  1.20381,  0.55853,  2.59722] ] ),
                ( 40, [ [ -2.63336,-11.75381, 14.13924], [  2.07109,  1.01048,  1.01275], [ -0.94505,  2.13211, -0.10995], [ -1.19500,  1.29010,  2.62884] ] ),
                ( 41, [ [ -1.34679, -9.76263, 14.88246], [  0.43346,  2.56584,  0.40897], [ -0.67824,  1.57386, -0.00913], [ -0.85439,  0.07357,  2.62111] ] ),
                ( 42, [ [ -4.78393,-11.50152, 13.19654], [  1.12434,  1.88464,  0.74238], [ -1.93124,  1.31679, -0.71499], [  0.33287,  1.24325,  2.91749] ] ),
                ( 43, [ [ -3.53394, -9.70992, 13.95975], [  1.37215,  1.69271,  0.78173], [ -0.85538,  1.95405, -0.24894], [ -0.46947, -0.15258,  2.45261] ] ),
                ( 44, [ [ -2.04827, -8.12657, 14.75512], [  1.59422,  1.46941,  0.80650], [ -0.72281,  1.69383, -0.24553], [ -0.31388, -0.89001,  2.71195] ] ),
                ( 45, [ [ -5.99636, -9.47494, 12.89247], [  1.70153,  1.70689,  0.76557], [ -0.70039,  2.26076, -0.31029], [  0.83634,  1.57429,  3.01028] ] ),
                ( 46, [ [ -4.34314, -7.84906, 13.64854], [  1.60474,  1.54469,  0.74649], [ -0.66560,  1.88514, -0.36153], [ -0.13861,  0.17375,  2.48261] ] ),
                ( 47, [ [ -2.78818, -6.38426, 14.38389], [  1.50499,  1.38473,  0.72413], [ -0.70779,  1.74800, -0.46650], [  0.53530, -0.43185,  3.00637] ] ),
                ( 48, [ [ -6.14433, -7.17716, 12.60211], [  1.21857,  1.17774,  0.66374], [  0.37835,  2.50490, -0.03228], [  1.15580,  1.37115,  2.63362] ] ),
                ( 49, [ [ -4.86378, -5.95285, 13.24022], [  1.34191,  1.27028,  0.61214], [ -0.41950,  1.59513, -0.39811], [  0.39946,  0.41073,  2.26153] ] ),
                ( 50, [ [ -3.46060, -4.63749, 13.82273], [  1.46385,  1.35988,  0.55265], [ -0.50746,  1.71135, -0.60949], [  0.85079, -0.09290,  2.33369] ] ),
                ( 51, [ [ -5.18891, -4.65904, 12.87258], [  1.22599,  2.11794,  0.29701], [ -0.23039,  0.99094, -0.33664], [  1.02691, -0.31517,  2.68807] ] ),
                ( 52, [ [ -3.80902, -2.97851, 13.17526], [  1.49629,  1.21273,  0.30081], [ -0.18838,  1.59811, -0.68182], [  1.26323, -0.91091,  2.64350] ] ),
                ( 53, [ [ -3.97013,-10.25768, 15.67056], [  1.39027,  1.09513,  2.04273], [ -1.07864,  2.17783,  0.78384], [ -0.35969, -0.84181,  0.86507] ] ),
                ( 54, [ [ -2.20978, -9.38789, 17.16098], [  2.01794,  0.61041,  0.88858], [ -0.38615,  1.79197,  1.71376], [ -0.00894, -1.59603,  2.05268] ] ),
                ( 55, [ [ -4.57962, -7.85650, 15.56861], [  2.11148,  1.24645,  2.39525], [ -0.34685,  2.32344, -0.47473], [ -0.32884, -0.18552,  1.33517] ] ),
                ( 56, [ [ -1.94683, -7.06974, 17.50631], [  3.03089,  0.31430,  1.42233], [ -0.14412,  2.64353, -0.63523], [  1.13996, -0.93302,  3.21747] ] ),
                ( 57, [ [ -4.66713, -5.75185, 14.78763], [  1.90006,  1.26128,  1.65660], [ -0.32912,  1.75003, -1.46178], [ -0.00605, -0.00857,  0.81862] ] ),
                ( 58, [ [ -2.54099, -4.72653, 15.86220], [  2.25081,  0.75532,  0.47130], [ -0.93933,  2.14774, -2.20592], [  0.98340, -0.08475,  1.73639] ] )
            ]

            generateParameters = False

            elementsCount1 = 2
            elementsCount2 = 4
            elementsCount3 = 4

            # Create nodes
            nodeIdentifier = 1
            lNodeIds = []
            nodeIndex = 0
            d1 = [0.5, 0.0, 0.0]
            d2 = [0.0, 0.5, 0.0]
            d3 = [0.0, 0.0, 1.0]
            for n3 in range(elementsCount3 + 1):
                lNodeIds.append([])
                for n2 in range(elementsCount2 + 1):
                    lNodeIds[n3].append([])
                    for n1 in range(elementsCount1 + 1):
                        lNodeIds[n3][n2].append(None)
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
                        lNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Create elements
            elementIdentifier = 1
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            lNodeIds[e3    ][e2][e1], lNodeIds[e3    ][e2][e1 + 1], lNodeIds[e3    ][e2 + 1][e1], lNodeIds[e3    ][e2 + 1][e1 + 1],
                            lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1]]

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

