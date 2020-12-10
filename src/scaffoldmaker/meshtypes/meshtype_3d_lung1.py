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
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion


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
        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        annotationGroups = [leftLungGroup, lungGroup, rightLungGroup]

        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)

        if isHuman:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            lowerRightLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of right lung"))
            lowerRightLungMeshGroup = lowerRightLungGroup.getMeshGroup(mesh)
            upperRightLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of right lung"))
            upperRightLungMeshGroup = upperRightLungGroup.getMeshGroup(mesh)
            middleRightLungGroup = AnnotationGroup(region, get_lung_term("middle lobe of right lung"))
            middleRightLungMeshGroup = middleRightLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            annotationGroups.append(upperLeftLungGroup)
            annotationGroups.append(lowerRightLungGroup)
            annotationGroups.append(upperRightLungGroup)
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
                (   7, [ [ 155.646, 121.977,-271.628], [ -83.355, -36.289,  19.187], [  35.892, -31.652,   7.476], [  14.994,  29.267,  38.119] ] ),
                (   8, [ [ 279.346,  98.455,-327.717], [ -21.666,  16.370,  35.485], [ -18.452, -43.354,   8.934], [  18.541,  53.843,  54.860] ] ),
                (   9, [ [ 251.887, 110.979,-294.259], [ -46.884,  -0.667,  13.029], [  -6.640, -34.923,  -9.542], [  -1.793,  34.831,  57.261] ] ),
                (  10, [ [ 203.263, 108.034,-281.647], [ -46.008, -21.960,  13.143], [  34.048, -17.679,  -6.475], [   5.003,  37.773,  28.793] ] ),
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
                (  44, [ [ 230.978,  62.444,-297.245], [ -14.698,  11.188,  -2.359], [  -8.962,  -8.354,   1.170], [   0.931,   7.312,  54.971] ] ),
                (  45, [ [ 258.329,  78.766,-234.210], [  -8.603,   4.498,   0.874], [ -37.867, -20.422,   2.056], [   3.424,   4.618,  67.449] ] ),
                (  46, [ [ 248.799,  84.344,-233.104], [ -10.725,   6.321,   0.445], [ -24.858, -21.841,  -4.453], [ -12.855,  -3.519,  64.769] ] ),
                (  47, [ [ 238.275,  93.157,-232.136], [ -12.861,   7.971,  -1.034], [  -9.789, -33.602, -12.451], [  -9.411,  14.946,  61.906] ] ),
                (  48, [ [ 223.573,  66.187,-240.080], [ -19.654,  31.273,  -7.973], [ -20.558, -22.685,  -0.982], [ -25.963,  -0.396,  55.488] ] ),
                (  49, [ [ 254.223,  82.226,-174.237], [ -21.821,  14.346,  -0.313], [ -70.819, -36.828,  -4.473], [  -6.474,   1.765,  58.654] ] ),
                (  50, [ [ 232.669,  96.602,-174.812], [ -20.714,  13.936,   0.272], [ -45.067, -35.411,  -8.030], [ -18.872,  10.884,  54.020] ] ),
                (  51, [ [ 211.888, 109.358,-175.186], [ -20.256,  12.504,   0.821], [ -16.313, -50.413, -13.321], [ -35.462,  14.924,  45.481] ] ),
                (  52, [ [ 187.821,  69.713,-187.140], [ -13.187,  61.646,   5.702], [ -34.715, -33.733, -12.517], [ -31.504,  10.855,  49.018] ] ),
                (  53, [ [ 213.425, 207.382, -42.148], [ -56.500,   0.342,  -5.827], [   6.048, -18.275, -16.938], [  -5.756, -22.958,  37.324] ] ),
                (  54, [ [ 258.130, 182.777, -53.571], [ -32.759,  -3.828,  -5.952], [  46.842, -36.257, -14.249], [ -42.970, -18.327,  45.780] ] ),
                (  55, [ [ 221.272, 179.757, -61.791], [ -41.743,  -3.435,  -5.875], [  10.486, -36.897, -21.690], [  -5.754, -11.017,  49.078] ] ),
                (  56, [ [ 175.167, 176.300, -67.698], [ -50.920,  -3.892,   0.663], [  -2.971, -33.698, -41.085], [  -2.018, -13.036,  51.511] ] ),
                (  57, [ [ 270.017, 129.272, -88.096], [ -48.699,  18.376,  -7.516], [ -17.418, -51.841, -36.718], [ -50.518, -29.109,  80.611] ] ),
                (  58, [ [ 224.626, 141.720, -98.406], [ -43.872,   3.149,  -4.298], [  -5.587, -42.256, -31.773], [   2.711, -18.020,  68.031] ] ),
                (  59, [ [ 185.274, 147.077,-102.145], [ -35.411,  -3.106,  -4.228], [  15.191, -29.940, -31.756], [ -14.714,  -1.454,  64.340] ] ),
                (  60, [ [ 236.417,  87.160,-119.825], [ -26.717,  14.046,  -6.516], [ -56.297, -42.646, -20.424], [ -33.135,   2.045,  67.489] ] ),
                (  61, [ [ 209.605, 101.124,-126.121], [ -27.728,  12.727,  -4.885], [ -42.756, -25.066, -21.644], [ -36.638,   1.272,  45.800] ] ),
                (  62, [ [ 181.792, 113.536,-131.292], [ -27.851,  13.168,  -2.607], [ -10.444, -42.179, -14.896], [ -23.513,  15.200,  47.902] ] ),
                (  63, [ [ 162.311,  78.408,-140.383], [  -5.741,  48.954, -11.944], [ -48.917, -26.754,   6.663], [ -11.741,  16.148,  53.528] ] ),
                (  64, [ [ 203.028, 174.619, -25.286], [ -60.155,  -2.415,  -3.955], [ -17.934, -44.396,   2.254], [  -5.864,   0.361,  32.179] ] ),
                (  65, [ [ 189.331, 131.538, -54.161], [ -73.019,  16.272, -14.723], [ -17.980, -26.665, -45.176], [  -2.116, -19.030,  46.007] ] ),
                (  66, [ [ 157.938, 103.099, -87.963], [ -38.581,  50.234, -40.028], [ -29.507, -44.106, -22.230], [ -56.230,  -1.913,  36.555] ] ),
                (  67, [ [ 112.805, 220.636,-344.408], [ -57.668,  31.639, -21.158], [ -27.490,  -7.728,  31.544], [  -7.261,  25.118,  94.161] ] ),
                (  68, [ [ 144.624, 170.181,-319.593], [ -45.954,  32.756,  25.839], [  29.801, -31.427,  46.872], [ -14.588,  30.925,  88.492] ] ),
                (  69, [ [  91.579, 200.374,-316.049], [ -54.312,  26.071, -18.610], [ -21.624, -23.063,  24.318], [   4.200,  -7.127,  69.610] ] ),
                (  70, [ [  45.375, 218.344,-353.416], [ -33.267,   4.881, -54.505], [ -61.732, -41.224,   6.887], [   0.041,  11.357,  87.727] ] ),
                (  71, [ [ 153.816, 137.480,-272.108], [ -79.344,  59.453, -10.789], [ -31.650, -32.032,  32.484], [  -3.217,  13.768,  41.588] ] ),
                (  72, [ [  74.446, 172.925,-295.847], [ -76.526,   2.015, -34.054], [ -12.428, -38.556,  16.527], [   4.041,  -7.852,  55.915] ] ),
                (  73, [ [  19.591, 150.815,-334.755], [  -3.346, -58.150, -34.776], [  -2.241, -66.932,  16.771], [ -37.214,  32.111,  65.155] ] ),
                (  74, [ [ 108.352, 123.286,-268.077], [ -31.450,  15.693,   5.225], [ -33.972, -20.910,  -4.547], [  22.886,  10.749,  29.999] ] ),
                (  75, [ [  69.151, 126.208,-283.827], [ -43.253, -15.831, -25.640], [   5.343, -46.935,   6.573], [  -2.573,   8.979,  43.158] ] ),
                (  76, [ [  36.295,  94.586,-323.318], [ -23.915, -46.011, -46.866], [  28.031, -41.992,  18.733], [ -35.404,  53.450,  65.721] ] ),
                (  77, [ [ 102.378,  96.563,-272.210], [ -19.527, -13.375,  -6.786], [   9.534, -20.075,   2.002], [   9.033,  23.903,  34.847] ] ),
                (  78, [ [  83.599,  83.999,-282.740], [ -14.132, -14.812, -15.431], [  32.104, -10.775,   2.129], [ -23.146,  18.638,  18.839] ] ),
                (  79, [ [  67.542,  71.106,-298.778], [ -14.808, -11.987, -18.636], [  32.907, -12.725,  14.740], [ -50.592,  24.213,  13.600] ] ),
                (  80, [ [ 109.300, 234.171,-248.531], [ -52.302,  25.979, -18.010], [ -19.304, -36.967,   0.740], [  -4.582,   1.836,  96.421] ] ),
                (  81, [ [ 135.331, 190.757,-238.135], [ -39.822,   0.430,  -4.178], [  22.077, -42.173,   7.556], [  -3.889,   9.995,  73.763] ] ),
                (  82, [ [  91.699, 199.576,-247.800], [ -46.192,  17.407, -15.151], [ -16.031, -32.155,   0.697], [  -3.291,   2.421,  66.468] ] ),
                (  83, [ [  46.055, 227.254,-268.328], [ -44.923,  35.895, -26.293], [ -53.545, -36.498, -12.333], [   1.708,   5.118,  82.513] ] ),
                (  84, [ [ 152.423, 152.232,-233.455], [ -80.327,  24.447,  -7.950], [  -7.571, -34.911,  -1.433], [  -4.926,  12.512,  36.525] ] ),
                (  85, [ [  77.265, 169.881,-247.131], [ -69.352,  10.546, -19.469], [ -10.609, -28.576,   2.252], [   1.874,   0.636,  41.232] ] ),
                (  86, [ [  14.896, 174.390,-270.710], [ -54.453,  -2.347, -27.972], [ -12.870, -45.688,   3.505], [  -0.893,  15.664,  55.196] ] ),
                (  87, [ [ 128.461, 142.532,-236.966], [ -58.190,   0.601,   0.855], [ -23.174, -13.575,  -4.695], [  19.504,  21.536,  33.741] ] ),
                (  88, [ [  70.179, 142.990,-243.592], [ -57.695,   0.305, -14.313], [  -4.787, -27.614,   2.693], [   4.708,  24.314,  36.010] ] ),
                (  89, [ [  14.153, 142.929,-264.777], [ -53.192,  -0.026, -29.221], [   3.625, -31.170,   8.059], [   4.661,  44.056,  51.260] ] ),
                (  90, [ [ 112.858, 122.968,-231.669], [ -45.051, -10.899,  -8.599], [ -10.184, -27.018,   2.473], [  14.630,  28.758,  45.554] ] ),
                (  91, [ [  67.804, 115.230,-241.806], [ -44.888,  -4.537, -11.644], [   7.968, -30.392,   4.316], [  -6.116,  42.327,  60.537] ] ),
                (  92, [ [  23.424, 113.892,-254.850], [ -43.703,   1.854, -14.388], [  21.202, -35.008,  12.469], [ -29.770,  57.684,  72.277] ] ),
                (  93, [ [ 104.207, 234.698,-152.733], [ -57.745,  26.895,  -7.074], [ -25.137, -25.018, -30.805], [   0.558,  -2.469,  77.754] ] ),
                (  94, [ [ 135.749, 192.199,-173.993], [ -57.597,   9.559, -12.340], [  21.057, -36.953, -26.041], [   2.245,  -7.457,  61.353] ] ),
                (  95, [ [  84.928, 206.459,-183.823], [ -43.517,  18.872,  -7.207], [ -13.067, -31.107, -30.940], [   0.766,  -1.586,  66.477] ] ),
                (  96, [ [  49.256, 228.018,-188.624], [ -27.295,  23.783,  -2.350], [ -46.031, -28.715, -40.634], [   2.799,  -2.519,  77.925] ] ),
                (  97, [ [ 146.834, 163.963,-199.827], [ -73.679,   9.289, -14.410], [  -1.766, -30.207, -31.660], [ -15.120,   1.557,  27.005] ] ),
                (  98, [ [  78.264, 173.394,-213.785], [ -63.451,   9.570, -13.504], [  -7.411, -31.806, -29.975], [   1.456,  -1.298,  43.628] ] ),
                (  99, [ [  19.915, 182.976,-226.748], [ -53.237,   9.593, -12.420], [ -16.164, -43.819, -39.572], [   3.936,   4.415,  60.584] ] ),
                ( 100, [ [ 108.125, 230.310, -93.528], [ -60.591,  28.132,   1.067], [ -13.931, -33.044, -12.007], [   4.596, -12.946,  57.696] ] ),
                ( 101, [ [ 141.968, 177.527,-116.980], [ -50.185,  17.132,  -3.637], [  13.607, -45.829, -52.443], [   9.725, -10.009,  55.550] ] ),
                ( 102, [ [  93.588, 196.028,-117.204], [ -45.032,  22.969,   1.979], [ -14.635, -34.317, -34.909], [   8.867, -14.370,  62.135] ] ),
                ( 103, [ [  51.988, 222.000,-112.818], [ -38.646,  27.605,   7.146], [ -48.820, -28.128, -36.207], [  13.793,  -1.876,  68.626] ] ),
                ( 104, [ [ 133.055, 155.147,-183.014], [ -50.399,   9.084,  23.517], [  -9.191, -33.132,  20.300], [  12.146,  20.340,  62.470] ] ),
                ( 105, [ [  80.120, 164.701,-163.613], [ -55.883,  12.515,  11.024], [   5.475, -51.414,  35.591], [  13.858,  40.686,  62.263] ] ),
                ( 106, [ [  22.839, 178.915,-159.278], [ -58.295,  15.237,  -1.951], [ -15.355, -83.455,  10.218], [  25.833,  57.492,  69.483] ] ),
                ( 107, [ [ 113.145,  88.046,-272.215], [  -9.070,  -8.548,  -7.601], [   9.970,  -9.918,  -1.262], [  -0.776,   3.825,  29.841] ] ),
                ( 108, [ [ 104.281,  79.019,-281.413], [  -8.787, -10.019, -10.080], [  18.983,  -1.180,   4.253], [ -24.412,   8.417,  37.368] ] ),
                ( 109, [ [  95.507,  67.421,-291.581], [  -8.328, -11.367, -12.477], [  28.245,   3.608,  12.840], [ -42.549,   9.803,  37.381] ] ),
                ( 110, [ [ 120.464,  78.319,-274.878], [  -4.724, -18.291, -14.002], [  13.797,   2.389,   6.929], [  17.827,  -6.047,  29.184] ] ),
                ( 111, [ [ 111.260,  92.585,-233.087], [ -21.088,  -3.634,   0.201], [  16.044, -25.578,  -6.389], [   7.310,  10.462,  47.435] ] ),
                ( 112, [ [  87.111,  87.499,-235.179], [ -27.022,  -6.613,  -4.654], [  33.506, -24.771,   1.279], [ -10.371,   9.996,  53.349] ] ),
                ( 113, [ [  58.498,  77.030,-243.404], [ -31.745,  -9.942, -10.923], [  61.050, -27.308,   9.195], [ -25.237,  14.981,  59.801] ] ),
                ( 114, [ [ 132.087,  72.782,-241.575], [  36.023, -40.841, -12.928], [  52.732, -11.634,  -8.495], [   5.541,  -5.631,  36.289] ] ),
                ( 115, [ [ 120.546, 104.051,-176.461], [ -29.943,  -7.987,   0.166], [   1.491, -47.972, -13.737], [  17.209,   6.934,  47.579] ] ),
                ( 116, [ [  86.014,  97.683,-176.263], [ -39.074,  -4.278,   0.229], [  32.208, -53.406, -23.992], [  11.839,   5.546,  52.988] ] ),
                ( 117, [ [  42.944,  92.427,-176.508], [ -47.599,   1.553,   0.293], [  68.016, -69.787, -28.090], [   7.893,   7.843,  67.952] ] ),
                ( 118, [ [ 131.468,  70.803,-204.192], [  26.047, -61.656, -37.983], [  50.646,  -1.156, -25.775], [  17.435, -17.925,  40.372] ] ),
                ( 119, [ [ 113.261, 209.394, -39.187], [ -58.016,  16.946,   3.696], [  -7.484, -28.184, -16.526], [   8.333, -32.245,  37.688] ] ),
                ( 120, [ [ 153.826, 168.085, -64.529], [ -56.166,  -3.576,   3.482], [   4.306, -46.192, -41.902], [ -10.285,  -1.501,  53.996] ] ),
                ( 121, [ [ 102.498, 178.383, -60.042], [ -46.388,  16.516,   4.961], [ -14.042, -33.069, -26.005], [  16.324, -10.383,  47.317] ] ),
                ( 122, [ [  64.189, 202.493, -54.339], [ -28.978,  30.762,   5.289], [ -43.522, -29.824, -24.792], [  33.818, -25.021,  55.310] ] ),
                ( 123, [ [ 123.540, 138.168,-101.485], [ -34.568,  -0.978,   8.985], [  -6.301, -36.302, -41.920], [  -7.205, -12.370,  65.576] ] ),
                ( 124, [ [  85.021, 144.013, -90.695], [ -43.680,   7.293,   9.469], [  -0.486, -42.997, -39.184], [  23.771, -17.396,  55.644] ] ),
                ( 125, [ [  37.425, 155.294, -80.316], [ -51.251,  17.213,   9.197], [   0.428, -62.132, -36.932], [  34.650, -21.397,  78.842] ] ),
                ( 126, [ [ 142.299, 107.091,-139.264], [ -33.152,  -9.595,   2.863], [  20.194, -35.241, -29.898], [  14.718,  -2.715,  45.939] ] ),
                ( 127, [ [ 107.369,  98.979,-132.633], [ -36.394,  -6.366,  10.682], [  41.976, -39.293, -37.022], [  29.911,   0.154,  43.505] ] ),
                ( 128, [ [  70.498,  93.229,-117.454], [ -37.588,  -2.021,  19.270], [  74.771, -47.522, -46.139], [  50.105, -10.662,  56.630] ] ),
                ( 129, [ [ 164.309,  71.103,-160.266], [  33.650, -70.498, -26.496], [  68.512, -15.582, -15.674], [  11.380,  14.458,  65.190] ] ),
                ( 130, [ [ 122.301, 174.072, -24.143], [ -59.195,  19.083,  18.719], [   3.799, -47.122,  -6.384], [  19.328,   1.610,  25.335] ] ),
                ( 131, [ [ 116.697, 132.839, -56.273], [ -61.128,   4.816,  37.451], [  11.875, -40.955, -35.954], [  33.215,  -5.911,  11.652] ] ),
                ( 132, [ [ 145.531,  97.776, -91.141], [ -48.613, -14.290,  47.805], [  40.807, -39.044, -38.399], [  47.670,  -2.075,  37.088] ] )
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
            lowerNodeIdsLeft = []
            upperNodeIdsLeft = []
            lowerNodeIdsRight = []
            upperNodeIdsRight = []
            d1 = [1.0, 0.0, 0.0]
            d2 = [0.0, 1.0, 0.0]
            d3 = [0.0, 0.0, 1.0]
            nodeIndex = 0
            leftLung = 0
            rightLung = 1

            for lung in (leftLung, rightLung):
                nodeIds = lowerNodeIdsLeft if (lung == leftLung) else lowerNodeIdsRight
                xMirror = 0 if (lung == leftLung) else 150 # Offset

                # Lower lobe nodes
                for n3 in range(lElementsCount3 + 1):
                    nodeIds.append([])
                    for n2 in range(lElementsCount2 + 1):
                        nodeIds[n3].append([])
                        for n1 in range(lElementsCount1 + 1):
                            nodeIds[n3][n2].append(None)
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
                            nodeIds[n3][n2][n1] = nodeIdentifier
                            nodeIdentifier += 1

                # Upper lobe nodes
                nodeIds = upperNodeIdsLeft if (lung == leftLung) else upperNodeIdsRight
                lowerNodeIds = lowerNodeIdsLeft if (lung == leftLung) else lowerNodeIdsRight

                for n3 in range(uElementsCount3 + 1):
                    nodeIds.append([])
                    for n2 in range(uElementsCount2 + 1):
                        nodeIds[n3].append([])
                        for n1 in range(uElementsCount1 + 1):
                            nodeIds[n3][n2].append(None)
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
                                nodeIds[n3][n2][n1] = lowerNodeIds[n3][lElementsCount2][n1]
                                continue
                            elif (n2 < (uElementsCount2 - 1)) and (n3 == (uElementsCount3 - 2)):
                                nodeIds[n3][n2][n1] = lowerNodeIds[lElementsCount3][n2][n1]
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
                            nodeIds[n3][n2][n1] = nodeIdentifier
                            nodeIdentifier += 1

            # Create elements
            elementIdentifier = 1
            for lung in (leftLung, rightLung):
                nodeIds = lowerNodeIdsLeft if (lung == leftLung) else lowerNodeIdsRight
                lowerNodeIds = lowerNodeIdsLeft if (lung == leftLung) else lowerNodeIdsRight

                # Lower lobe elements
                for e3 in range(lElementsCount3):
                    for e2 in range(lElementsCount2):
                        for e1 in range(lElementsCount1):
                            eft = eftRegular
                            nodeIdentifiers = [
                                nodeIds[e3    ][e2][e1], nodeIds[e3    ][e2][e1 + 1], nodeIds[e3    ][e2 + 1][e1], nodeIds[e3    ][e2 + 1][e1 + 1],
                                nodeIds[e3 + 1][e2][e1], nodeIds[e3 + 1][e2][e1 + 1], nodeIds[e3 + 1][e2 + 1][e1], nodeIds[e3 + 1][e2 + 1][e1 + 1]]

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
                nodeIds = upperNodeIdsLeft if (lung == leftLung) else upperNodeIdsRight
                for e3 in range(uElementsCount3):
                    for e2 in range(uElementsCount2):
                        for e1 in range(uElementsCount1):
                            eft = eftRegular
                            nodeIdentifiers = [
                                nodeIds[e3][e2][e1], nodeIds[e3][e2][e1 + 1], nodeIds[e3][e2 + 1][e1],
                                nodeIds[e3][e2 + 1][e1 + 1],
                                nodeIds[e3 + 1][e2][e1], nodeIds[e3 + 1][e2][e1 + 1], nodeIds[e3 + 1][e2 + 1][e1],
                                nodeIds[e3 + 1][e2 + 1][e1 + 1]]

                            if (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                                # Distal-front wedge elements
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(2)
                                eft = eftWedgeCollapseXi1_37
                            elif (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == (uElementsCount1 - 1)):
                                # Distal-back wedge elements
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(3)
                                eft = eftWedgeCollapseXi1_48
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
                                # Top-front-distal tetrahedron wedge elements
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(2)
                                eft = eftTetCollapseXi1Xi2_63
                            elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == (uElementsCount1 - 1)):
                                # Top-front-distal tetrahedron wedge elements
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(3)
                                eft = eftTetCollapseXi1Xi2_53
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
            elementsCount1 = 2
            elementsCount2 = 4
            elementsCount3 = 4

            # Create nodes
            nodeIdentifier = 1
            lowerNodeIdsLeft = []
            d1 = [0.5, 0.0, 0.0]
            d2 = [0.0, 0.5, 0.0]
            d3 = [0.0, 0.0, 1.0]
            for n3 in range(elementsCount3 + 1):
                lowerNodeIdsLeft.append([])
                for n2 in range(elementsCount2 + 1):
                    lowerNodeIdsLeft[n3].append([])
                    for n1 in range(elementsCount1 + 1):
                        lowerNodeIdsLeft[n3][n2].append(None)
                        if n3 < elementsCount3:
                            if (n1 == 0) and ((n2 == 0) or (n2 == elementsCount2)):
                                continue
                        else:
                            if (n2 == 0) or (n2 == elementsCount2) or (n1 == 0):
                                continue
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        x = [0.5 * (n1 - 1), 0.5 * (n2 - 1), 1.0 * n3]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        lowerNodeIdsLeft[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Create elements
            elementIdentifier = 1
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            lowerNodeIdsLeft[e3    ][e2][e1], lowerNodeIdsLeft[e3    ][e2][e1 + 1], lowerNodeIdsLeft[e3    ][e2 + 1][e1], lowerNodeIdsLeft[e3    ][e2 + 1][e1 + 1],
                            lowerNodeIdsLeft[e3 + 1][e2][e1], lowerNodeIdsLeft[e3 + 1][e2][e1 + 1], lowerNodeIdsLeft[e3 + 1][e2 + 1][e1], lowerNodeIdsLeft[e3 + 1][e2 + 1][e1 + 1]]

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

