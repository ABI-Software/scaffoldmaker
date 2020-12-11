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
        #rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        annotationGroups = [leftLungGroup, lungGroup]

        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        #rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)

        if isHuman:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            annotationGroups.append(upperLeftLungGroup)

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
                (  1, [ [  204.214,  217.387, -341.252 ], [  -62.947,  -35.911,    0.688 ], [   25.814,  -29.305,   29.820 ], [   -4.512,   17.903,   91.032 ] ] ),
                (  2, [ [  269.306,  209.072, -354.889 ], [  -20.103,   -4.908,   79.411 ], [   47.047,  -49.139,   -2.543 ], [    0.175,   15.720,   89.961 ] ] ),
                (  3, [ [  223.341,  188.444, -308.232 ], [  -60.426,  -24.274,   -1.081 ], [   20.503,  -36.477,   21.902 ], [    3.351,   18.019,   52.584 ] ] ),
                (  4, [ [  175.499,  161.680, -317.300 ], [  -27.384,  -29.319,  -26.704 ], [    9.020,  -29.937,   61.149 ], [    2.480,   21.158,   82.660 ] ] ),
                (  5, [ [  292.990,  153.587, -341.674 ], [  -10.831,   26.074,   66.575 ], [  -10.247,  -58.090,   17.240 ], [   10.595,   29.035,   61.452 ] ] ),
                (  6, [ [  241.787,  149.543, -291.578 ], [  -78.934,  -21.344,   11.789 ], [   22.813,  -35.233,   12.622 ], [   -6.632,   28.612,   43.054 ] ] ),
                (  7, [ [  159.624,  124.828, -273.695 ], [  -83.355,  -36.289,   19.187 ], [   42.385,  -36.600,    4.581 ], [   14.994,   29.267,   38.119 ] ] ),
                (  8, [ [  279.346,   98.455, -327.717 ], [  -21.666,   16.370,   35.485 ], [  -18.452,  -43.354,    8.934 ], [   18.541,   53.843,   54.860 ] ] ),
                (  9, [ [  251.887,  110.979, -294.259 ], [  -46.884,   -0.667,   13.029 ], [   -6.640,  -34.923,   -9.542 ], [   -1.793,   34.831,   57.261 ] ] ),
                ( 10, [ [  203.673,  108.786, -280.769 ], [  -46.008,  -21.960,   13.143 ], [   34.048,  -17.679,   -6.475 ], [    5.003,   37.773,   28.793 ] ] ),
                ( 11, [ [  256.681,   71.412, -322.006 ], [  -14.537,   11.023,   15.628 ], [   -7.850,  -10.251,   25.280 ], [   35.613,   35.978,   71.913 ] ] ),
                ( 12, [ [  243.935,   80.999, -302.872 ], [  -16.394,    8.598,   18.165 ], [   -8.074,  -19.651,   -6.343 ], [   10.839,   34.046,   66.123 ] ] ),
                ( 13, [ [  226.628,   91.702, -285.892 ], [  -23.230,   -0.285,   15.358 ], [    7.839,  -18.769,  -11.448 ], [    6.974,   29.647,   42.160 ] ] ),
                ( 14, [ [  217.057,  233.615, -251.001 ], [  -56.387,  -12.798,   -2.074 ], [    6.551,  -32.468,    4.928 ], [    7.180,    7.996,   92.202 ] ] ),
                ( 15, [ [  267.567,  218.538, -268.608 ], [  -31.271,  -28.559,   17.780 ], [   28.703,  -41.920,  -22.876 ], [   -2.499,    2.835,   82.368 ] ] ),
                ( 16, [ [  227.626,  202.773, -250.316 ], [  -46.421,   -9.158,   11.317 ], [    7.297,  -31.033,    0.920 ], [   -0.016,    8.237,   63.668 ] ] ),
                ( 17, [ [  178.790,  194.631, -246.555 ], [  -50.393,   -5.414,   -8.376 ], [  -22.308,  -44.954,   12.222 ], [    3.296,   11.647,   70.649 ] ] ),
                ( 18, [ [  296.250,  178.154, -283.773 ], [  -52.959,   -4.397,   27.035 ], [   10.998,  -43.061,   -0.027 ], [   -2.037,    9.722,   56.957 ] ] ),
                ( 19, [ [  240.706,  174.731, -251.298 ], [  -65.503,  -16.663,   18.653 ], [   12.413,  -26.875,    3.862 ], [   -0.209,    7.605,   43.189 ] ] ),
                ( 20, [ [  170.036,  151.299, -240.510 ], [  -77.888,  -18.667,    9.104 ], [   21.815,  -36.197,    2.313 ], [   11.396,   18.147,   30.385 ] ] ),
                ( 21, [ [  297.502,  143.355, -275.679 ], [  -48.044,    9.944,   32.993 ], [   -5.929,  -36.823,   16.652 ], [   -0.988,   42.077,   47.842 ] ] ),
                ( 22, [ [  250.978,  148.431, -247.195 ], [  -50.687,    1.846,    9.041 ], [    9.032,  -28.662,    9.417 ], [   -4.776,   32.784,   41.932 ] ] ),
                ( 23, [ [  204.680,  141.636, -246.513 ], [  -35.493,  -20.244,   16.094 ], [   31.013,  -11.624,    3.008 ], [  -11.195,   25.677,   39.181 ] ] ),
                ( 24, [ [  287.464,  106.726, -251.655 ], [  -25.621,    6.219,   23.134 ], [  -19.294,  -29.795,   27.926 ], [   17.692,   37.852,   69.883 ] ] ),
                ( 25, [ [  257.922,  116.607, -238.393 ], [  -31.830,    9.651,    4.808 ], [   -9.432,  -33.031,    3.881 ], [  -10.300,   42.470,   62.047 ] ] ),
                ( 26, [ [  228.110,  129.472, -238.391 ], [  -30.191,    5.166,   -9.731 ], [   19.291,  -27.684,    5.002 ], [  -37.186,   41.031,   47.384 ] ] ),
                ( 27, [ [  219.598,  234.911, -158.376 ], [  -59.865,  -18.569,  -15.474 ], [    6.365,  -34.542,  -21.886 ], [   -2.376,   -5.948,   82.683 ] ] ),
                ( 28, [ [  271.479,  212.598, -191.075 ], [  -45.292,  -11.794,    7.859 ], [   30.495,  -31.862,  -46.294 ], [    1.874,  -12.175,   77.537 ] ] ),
                ( 29, [ [  226.886,  201.943, -182.154 ], [  -46.036,   -0.006,    5.281 ], [    5.759,  -30.320,  -27.476 ], [    4.237,  -12.188,   64.641 ] ] ),
                ( 30, [ [  181.157,  202.415, -182.448 ], [  -45.198,   -1.262,    0.548 ], [   -3.376,  -48.006,  -27.616 ], [   -1.725,   -4.852,   66.277 ] ] ),
                ( 31, [ [  291.428,  178.268, -232.360 ], [  -43.644,  -20.853,   24.989 ], [    9.889,  -35.193,  -43.486 ], [    4.353,  -16.910,   50.697 ] ] ),
                ( 32, [ [  237.268,  175.712, -212.165 ], [  -56.954,  -16.756,   10.424 ], [   13.032,  -26.067,  -32.916 ], [    0.506,  -11.344,   40.699 ] ] ),
                ( 33, [ [  177.122,  161.633, -212.853 ], [  -60.962,  -16.488,   -9.264 ], [   19.400,  -26.292,  -38.189 ], [    2.855,    9.850,   36.322 ] ] ),
                ( 34, [ [  219.695,  225.898,  -85.992 ], [  -73.383,   -7.270,    0.075 ], [    9.466,  -34.972,  -20.253 ], [   -0.443,  -16.457,   58.378 ] ] ),
                ( 35, [ [  276.870,  199.241, -113.455 ], [  -50.589,   -4.004,    2.780 ], [   27.986,  -43.217,  -64.463 ], [   -2.551,  -14.345,   71.169 ] ] ),
                ( 36, [ [  228.512,  190.564, -119.932 ], [  -48.632,   -2.241,    2.901 ], [   10.553,  -37.072,  -44.027 ], [   -7.211,  -10.078,   60.219 ] ] ),
                ( 37, [ [  181.214,  185.829, -116.571 ], [  -45.024,  -10.633,    6.490 ], [   36.691,  -40.359,  -41.741 ], [   -7.404,  -14.171,   56.958 ] ] ),
                ( 38, [ [  294.618,  150.473, -187.485 ], [  -61.912,    3.854,   13.211 ], [   -7.435,  -61.052,   31.811 ], [    4.313,   46.520,   70.475 ] ] ),
                ( 39, [ [  237.513,  155.974, -176.526 ], [  -52.879,    5.296,    5.774 ], [   -7.421,  -52.768,    1.120 ], [   -9.440,   39.508,   59.010 ] ] ),
                ( 40, [ [  189.392,  160.334, -172.880 ], [  -42.801,    6.979,    3.264 ], [   21.841,  -37.859,   28.056 ], [   -6.732,   20.365,   70.432 ] ] ),
                ( 41, [ [  247.653,   64.318, -306.449 ], [   -7.559,    4.781,    7.258 ], [  -10.993,   -6.381,   15.496 ], [   30.221,   28.933,   71.221 ] ] ),
                ( 42, [ [  239.325,   70.695, -300.571 ], [   -8.415,    6.741,    6.479 ], [   -7.094,   -9.314,    3.351 ], [   10.779,   13.380,   72.272 ] ] ),
                ( 43, [ [  231.513,   77.178, -292.656 ], [  -12.525,    4.255,    3.319 ], [   -0.928,  -16.243,   -4.197 ], [   11.706,   13.532,   59.617 ] ] ),
                ( 44, [ [  230.978,   62.444, -297.245 ], [  -14.698,   11.188,   -2.359 ], [   -8.962,   -8.354,    1.170 ], [    0.931,    7.312,   54.971 ] ] ),
                ( 45, [ [  258.329,   78.766, -234.210 ], [   -8.603,    4.498,    0.874 ], [  -37.867,  -20.422,    2.056 ], [    3.424,    4.618,   67.449 ] ] ),
                ( 46, [ [  248.799,   84.344, -233.104 ], [  -10.725,    6.321,    0.445 ], [  -24.858,  -21.841,   -4.453 ], [  -12.855,   -3.519,   64.769 ] ] ),
                ( 47, [ [  238.275,   93.157, -232.136 ], [  -12.861,    7.971,   -1.034 ], [   -9.789,  -33.602,  -12.451 ], [   -9.411,   14.946,   61.906 ] ] ),
                ( 48, [ [  223.573,   66.187, -240.080 ], [  -19.654,   31.273,   -7.973 ], [  -20.558,  -22.685,   -0.982 ], [  -25.963,   -0.396,   55.488 ] ] ),
                ( 49, [ [  253.968,   82.516, -174.339 ], [  -21.821,   14.346,   -0.313 ], [  -70.819,  -36.828,   -4.473 ], [   -6.474,    1.765,   58.654 ] ] ),
                ( 50, [ [  232.669,   96.602, -174.812 ], [  -20.714,   13.936,    0.272 ], [  -45.067,  -35.411,   -8.030 ], [  -18.872,   10.884,   54.020 ] ] ),
                ( 51, [ [  211.888,  109.358, -175.186 ], [  -20.256,   12.504,    0.821 ], [  -16.313,  -50.413,  -13.321 ], [  -35.462,   14.924,   45.481 ] ] ),
                ( 52, [ [  187.821,   69.713, -187.140 ], [  -13.187,   61.646,    5.702 ], [  -34.715,  -33.733,  -12.517 ], [  -31.504,   10.855,   49.018 ] ] ),
                ( 53, [ [  213.425,  207.382,  -42.148 ], [  -56.500,    0.342,   -5.827 ], [    6.048,  -18.275,  -16.938 ], [   -5.756,  -22.958,   37.324 ] ] ),
                ( 54, [ [  258.130,  182.777,  -53.571 ], [  -32.759,   -3.828,   -5.952 ], [   46.842,  -36.257,  -14.249 ], [  -42.970,  -18.327,   45.780 ] ] ),
                ( 55, [ [  221.272,  179.757,  -61.791 ], [  -41.743,   -3.435,   -5.875 ], [   10.486,  -36.897,  -21.690 ], [   -5.754,  -11.017,   49.078 ] ] ),
                ( 56, [ [  175.167,  176.300,  -67.698 ], [  -50.920,   -3.892,    0.663 ], [   -2.971,  -33.698,  -41.085 ], [   -2.018,  -13.036,   51.511 ] ] ),
                ( 57, [ [  270.017,  129.272,  -88.096 ], [  -48.699,   18.376,   -7.516 ], [  -17.418,  -51.841,  -36.718 ], [  -50.518,  -29.109,   80.611 ] ] ),
                ( 58, [ [  224.626,  141.720,  -98.406 ], [  -43.872,    3.149,   -4.298 ], [   -5.587,  -42.256,  -31.773 ], [    2.711,  -18.020,   68.031 ] ] ),
                ( 59, [ [  185.274,  147.077, -102.145 ], [  -35.411,   -3.106,   -4.228 ], [   15.191,  -29.940,  -31.756 ], [  -14.714,   -1.454,   64.340 ] ] ),
                ( 60, [ [  236.417,   87.160, -119.825 ], [  -26.717,   14.046,   -6.516 ], [  -56.297,  -42.646,  -20.424 ], [  -33.135,    2.045,   67.489 ] ] ),
                ( 61, [ [  209.605,  101.124, -126.121 ], [  -27.728,   12.727,   -4.885 ], [  -42.756,  -25.066,  -21.644 ], [  -36.638,    1.272,   45.800 ] ] ),
                ( 62, [ [  181.792,  113.536, -131.292 ], [  -27.851,   13.168,   -2.607 ], [  -10.444,  -42.179,  -14.896 ], [  -23.513,   15.200,   47.902 ] ] ),
                ( 63, [ [  162.311,   78.408, -140.383 ], [  -34.319,   54.554,   -5.830 ], [  -48.917,  -26.754,    6.663 ], [  -11.741,   16.148,   53.528 ] ] ),
                ( 64, [ [  203.028,  174.619,  -25.286 ], [  -60.155,   -2.415,   -3.955 ], [  -17.934,  -44.396,    2.254 ], [   -5.864,    0.361,   32.179 ] ] ),
                ( 65, [ [  189.331,  131.538,  -54.161 ], [  -73.019,   16.272,  -14.723 ], [  -17.980,  -26.665,  -45.176 ], [   -2.116,  -19.030,   46.007 ] ] ),
                ( 66, [ [  159.408,  102.276,  -88.784 ], [  -38.581,   50.234,  -40.028 ], [  -29.507,  -44.106,  -22.230 ], [  -56.230,   -1.913,   36.555 ] ] )
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
            lNodeIds = []
            uNodeIds = []
            d1 = [1.0, 0.0, 0.0]
            d2 = [0.0, 1.0, 0.0]
            d3 = [0.0, 0.0, 1.0]
            nodeIndex = 0
            # Lower lobe nodes
            for n3 in range(lElementsCount3 + 1):
                lNodeIds.append([])
                for n2 in range(lElementsCount2 + 1):
                    lNodeIds[n3].append([])
                    for n1 in range(lElementsCount1 + 1):
                        lNodeIds[n3][n2].append(None)
                        if ((n1 == 0) or (n1 == lElementsCount1)) and (n2 == 0):
                            continue
                        if (n3 > (lElementsCount3 - 2)) and (n2 > (lElementsCount2 - 2)):
                            continue

                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        if generateParameters:
                            x = [1.0 * (n1 - 1), 1.0 * (n2 - 1), 1.0 * n3]
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

            # Upper lobe nodes
            for n3 in range(uElementsCount3 + 1):
                uNodeIds.append([])
                for n2 in range(uElementsCount2 + 1):
                    uNodeIds[n3].append([])
                    for n1 in range(uElementsCount1 + 1):
                        uNodeIds[n3][n2].append(None)
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
                            uNodeIds[n3][n2][n1] = lNodeIds[n3][lElementsCount2][n1]
                            #print('uNodeIds', uNodeIds[n3][n2][n1])
                            continue
                        elif (n2 < (uElementsCount2 - 1)) and (n3 == (uElementsCount3 - 2)):
                            uNodeIds[n3][n2][n1] = lNodeIds[lElementsCount3][n2][n1]
                            #print('uNodeIds', uNodeIds[n3][n2][n1])
                            continue

                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        if generateParameters:
                            x = [1.0 * (n1 - 1), 1.0 * (n2 - 1) + 2.5, 1.0 * n3 + 2.0]
                        else:
                            nodeParameters = nodeFieldParameters[nodeIndex]
                            nodeIndex += 1
                            assert nodeIdentifier == nodeParameters[0]
                            x, d1, d2, d3 = nodeParameters[1]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        uNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Create elements
            # Lower lobe elements
            elementIdentifier = 1
            for e3 in range(lElementsCount3):
                for e2 in range(lElementsCount2):
                    for e1 in range(lElementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            lNodeIds[e3    ][e2][e1], lNodeIds[e3    ][e2][e1 + 1], lNodeIds[e3    ][e2 + 1][e1], lNodeIds[e3    ][e2 + 1][e1 + 1],
                            lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1]]

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
                            nodeIdentifiers[2] = lNodeIds[e3 - 1][e2 + 1][e1    ]
                            nodeIdentifiers[3] = lNodeIds[e3 - 1][e2 + 1][e1 + 1]
                            nodeIdentifiers[6] = lNodeIds[e3 - 1][e2 + 2][e1    ]
                            nodeIdentifiers[7] = lNodeIds[e3 - 1][e2 + 2][e1 + 1]
                            eft = eftfactory.createEftBasic()
                            setEftScaleFactorIds(eft, [1], [])
                            remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif None in nodeIdentifiers:
                            continue

                        # print('element ', elementIdentifier, '|| ', nodeIdentifiers)
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
                        lowerLeftLungMeshGroup.addElement(element)

            # Upper lobe elements
            for e3 in range(uElementsCount3):
                for e2 in range(uElementsCount2):
                    for e1 in range(uElementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            uNodeIds[e3][e2][e1], uNodeIds[e3][e2][e1 + 1], uNodeIds[e3][e2 + 1][e1],
                            uNodeIds[e3][e2 + 1][e1 + 1],
                            uNodeIds[e3 + 1][e2][e1], uNodeIds[e3 + 1][e2][e1 + 1], uNodeIds[e3 + 1][e2 + 1][e1],
                            uNodeIds[e3 + 1][e2 + 1][e1 + 1]]

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

                        # print('element ', elementIdentifier, '|| ', nodeIdentifiers)
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
                        upperLeftLungMeshGroup.addElement(element)

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
