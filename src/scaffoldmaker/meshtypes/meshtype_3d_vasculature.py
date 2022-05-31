"""
Generates a 3-D unit box mesh with variable numbers of elements in 3 directions.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldFiniteElement
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup


class MeshType_3d_vasculature1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Vasculature 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements 1' : 1,
            'Number of elements 2' : 1,
            'Number of elements 3' : 1,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements 1' : 1,
            'Refine number of elements 2' : 1,
            'Refine number of elements 3' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        # elementsCount1 = options['Number of elements 1']
        # elementsCount2 = options['Number of elements 2']
        # elementsCount3 = options['Number of elements 3']
        # useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        radius = findOrCreateFieldFiniteElement(fm, "radius", components_count=1, managed=True)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.defineField(radius)
        nodetemplate.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)

        nodetemplate1 = nodes.createNodetemplate()
        nodetemplate1.defineField(coordinates)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 3)
        nodetemplate1.defineField(radius)
        nodetemplate1.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        linearBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft1.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 1)
        eft2 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft2.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 2)
        eft3 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft3.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 3)
        eft32 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft32.setTermNodeParameter(2, 1, 2, Node.VALUE_LABEL_D_DS1, 3)
        eftRadius = mesh.createElementfieldtemplate(linearBasis)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        elementtemplate1.defineField(radius, -1, eftRadius)

        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate2.defineField(coordinates, -1, eft2)
        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate3.defineField(coordinates, -1, eft3)
        elementtemplate2.defineField(radius, -1, eftRadius)
        elementtemplate3.defineField(radius, -1, eftRadius)



        cache = fm.createFieldcache()

        vasculature = False
        vagus = True
        if vasculature:
            x_list = [
            [20.9901, -55.2915, 1237.12],
            [21.7547, -56.1276, 1268.6],
            [24.838, -74.4412, 1299.32],
            [21.3563, -91.2, 1307.97],
            [12.9058, -111.015, 1302.16],
            [25.4696, -95.415, 1344.72],
            [39.5394, -109.4, 1435.59],
            [-10.5303, -102.14, 1331.9],
            [-27.2225, -91.4847, 1430.7],
            [-27.5775, -97.4023, 1345.7],
            [5.157545857793905, -73.7144264439337, 1067.611422937166],
            [0.7123922184414442, -76.12139177722604, 949.3444291866757],
            [1.480653991761836, -73.62566826677576, 815.1125925698087],
            [-19.50824359035452, -72.86872429654028, 929.2051634762531],
            [-20.28614848561935, -73.51853769979229, 899.0024026991548],
            [-21.61102738650197, -74.8636341463615, 841.6568853077786],
            [-18.32152657377993, -80.08693448724975, 975.6309361438332],
            [82.33929143728056, -75.38569430920019, 965.0794458269447],
            [-99.29494952824054, -84.197838408155, 963.3384151864503],
            [-14.63217490185608, -80.76093647541154, 1026.190525855144],
            [-100.6218523253879, -81.37468570139932, 1070.890294941285],
            [17.07417461634516, -120.3755886862945, 1249.284411678494]
                ]

            d1_list = [
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
            ]

            nodes_w_ver = [5,8,4,3,11,20,17]
        elif vagus:
            x_list = [
                [-3.067919428038835, -58.88831896656854, 1485.94100956215],
                [-0.6281586078592581, -62.69296504830071, 1481.347190265109],
                [19.14045007878925, -71.17249683579628, 1471.653806089609],
                [35.50375246183219, -78.35854544883568, 1463.384930602326],
                [36.38036257645483, -80.58804084795787, 1457.154198693773],
                [35.70333609041023, -82.09528266416706, 1448.653619505636],
                [34.97165431464403, -83.39378105979046, 1441.867216455742],
                [34.26489569026836, -84.45087019096185, 1435.143450297613],
                [32.29029472400944, -92.73679568270747, 1396.092960178755],
                [28.92067653338183, -109.6287737218008, 1316.75069003233],
                [28.32112681266167, -112.2390790222721, 1305.054388871628],
                [6.058014106179172, -96.19324278756504, 1168.221198705996],
                [4.49003139871628, -120.3327063642019, 982.41199971608],
                [-4.484831774677777, -76.72363699970431, 919.0809110751771],
                [-9.096316816507812, -18.3837646869268, 858.4382217409927],
                [17.463, -71.807, 1476.875],
                [19.059, -76.211, 1480.118],
                [15.862, -66.876, 1461.220],
                [14.086, -75.333, 1461.325],
                [4.653, -62.902, 1463.956],
                [12.273, -59.511, 1466.982],
                [-6.224, -58.028, 1472.469],
                [-57.275, -70.875, 1519.911],
                [-65.468, -70.375, 1520.633],
                [60.913, -52.711, 1519.712],
                [0.526, -56.094, 1488.768],
                [11.232, -116.557, 1452.669],
                [9.582, -154.877, 1443.576],
                [-14.000, -152.355, 1446.593],
                [1.878, -56.261, 1489.300],
                [-3.468, -73.077, 1344.169],
                [2.607, -56.401, 1489.594],
                [20.440, -92.164, 1438.379],
                [14.327, -97.349, 1437.414],
                [19.268, -95.333, 1427.899],
                [13.106, -99.879, 1430.559],
                [6.677, -105.285, 1436.500],
                [-0.539, -131.251, 1434.710],
                [9.361, -130.529, 1433.801],
                [5.960, -114.822, 1414.415],
                [1.823, -112.755, 1413.751],
                [0.392, -105.992, 1406.523],
                [4.033, -107.383, 1404.800],
                [0.755, -107.608, 1406.119],
                [1.376, -109.143, 1406.395],
                [12.367, -99.729, 1421.523],
                [4.832, -136.024, 1448.806],
                [19.826, -96.819, 1421.789],
                [8.188, -105.664, 1400.628],
                [7.465, -106.606, 1402.082],
                [6.623, -108.164, 1402.941],
                [6.818, -105.579, 1401.325],
                [2.098, -109.068, 1405.412],
                [29.099, -94.640, 1421.734],
                [27.685, -97.760, 1415.520],
                [25.584, -101.952, 1409.817],
                [30.962, -99.435, 1395.981],
                [30.626, -102.783, 1394.211],
                [30.621, -102.116, 1398.184],
                [12.808, -107.604, 1396.148],
                [12.660, -104.603, 1398.421],
                [13.545, -109.239, 1392.497],
                [13.373, -107.842, 1393.326],
                [12.293, -108.004, 1394.644],
                [11.744, -108.913, 1398.904],
                [19.598, -110.678, 1389.008],
                [20.305, -109.662, 1385.525],
                [13.402, -106.487, 1391.751],
                [17.812, -111.228, 1390.543],
                [16.942, -106.518, 1385.462],
                [18.305, -110.622, 1389.761],
                [18.880, -110.329, 1387.243],
                [19.183, -111.350, 1388.633],
                [2.263, -68.435, 1345.906],
                [0.501, -67.230, 1351.551],
                [3.296, -68.465, 1338.180],
                [20.571, -47.974, 1214.357],
                [13.336, -56.091, 1200.123],
                [4.126, -91.315, 1143.098],
                [-0.849, -78.626, 1059.708],
                [-0.497, -88.117, 1055.124],
                [4.051, -80.134, 952.868],
                [2.251, -72.056, 934.224],
                [84.559, -64.284, 856.328],
                [15.735, -106.477, 1387.002],
            ]

            d1_list = [
                [-5.129, -2.099, -5.383],
                [25.147, 6.953, 2.966],
                [3.230, -2.449, -7.906],
                [-1.057, -4.299, -8.305],
                [-0.118, -2.112, -6.641],
                [-6.773, -2.898, -16.891],
                [0.902, 1.073, -8.921],
                [-3.332, -7.049, -20.250],
                [0.321, -5.982, -10.694],
                [-3.914, -3.308, -13.084],
                [0.365, -9.088, -15.532],
                [-7.486, 1.904, -38.190],
                [-5.290, 48.180, -94.925],
                [2.729, 53.763, -37.471],
                [73.681, 23.211, -65.919],
                [-6.734, 0.139, 6.414],
                [-6.259, 0.228, 3.865],
                [-16.113, -0.110, 1.188],
                [3.242, -4.337, -0.252],
                [-7.282, 8.493, 9.077],
                [7.606, -0.578, -1.022],
                [-26.479, -6.517, 43.508],
                [-65.695, -23.698, 36.302],
                [-80.971, -22.695, 37.571],
                [113.397, 7.643, 25.844],
                [-3.962, 3.013, 3.485],
                [-17.740, -38.752, -6.999],
                [12.923, -33.908, -10.012],
                [-32.062, -32.179, -5.049],
                [-3.416, 3.366, 3.014],
                [-68.816, 21.122, -164.179],
                [-3.918, 3.780, 3.375],
                [-5.057, -3.352, -3.807],
                [-5.856, -6.144, -4.017],
                [-2.347, -4.675, -8.715],
                [-7.534, -5.737, -20.761],
                [-6.059, -19.197, -11.926],
                [-7.680, -30.025, 7.655],
                [9.952, -27.252, 5.685],
                [3.841, 0.102, -26.783],
                [-3.000, 3.499, -27.600],
                [-17.761, -6.441, -27.108],
                [-10.611, -9.270, -30.754],
                [-17.080, -9.671, -27.975],
                [-15.831, -12.715, -27.404],
                [1.399, 2.283, -5.055],
                [-15.246, -52.924, 50.184],
                [0.141, 1.614, -5.649],
                [-12.738, -7.253, -31.429],
                [-14.226, -9.153, -28.678],
                [-15.945, -12.281, -27.056],
                [-15.423, -7.075, -30.010],
                [-14.437, -12.602, -29.442],
                [-3.226, -6.683, -10.800],
                [-1.761, -3.668, -5.999],
                [-2.883, -5.866, -8.726],
                [-4.261, -2.824, 1.454],
                [-0.863, -1.635, -4.163],
                [-5.481, -5.016, 2.894],
                [-21.763, -7.091, 1.085],
                [-21.174, -1.106, 5.361],
                [-20.354, -10.352, -6.156],
                [-20.648, -7.563, -4.500],
                [-22.842, -7.903, -1.896],
                [-23.743, -9.631, 6.505],
                [-10.296, -9.173, -9.205],
                [-8.697, -6.995, -15.815],
                [-20.209, -4.791, -7.477],
                [-13.768, -10.203, -6.101],
                [-15.129, -0.845, -15.793],
                [-12.874, -9.060, -7.701],
                [-11.712, -8.465, -12.709],
                [-11.115, -10.505, -9.945],
                [-55.339, 23.836, -191.394],
                [-58.220, 25.910, -179.434],
                [-0.341, -1.935, -12.395],
                [-11.614, 110.759, -121.987],
                [-4.754, 7.795, -8.567],
                [-6.426, 7.025, -34.377],
                [-5.782, 17.092, -99.868],
                [-5.067, -1.880, -108.785],
                [-0.119, 8.718, 0.525],
                [1.907, -1.596, 6.514],
                [82.082, -83.081, 44.569],
                [-17.503, -0.767, -12.83],

                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                # [10.0, 0.0, 0.0],
                ]


            nodes_w_ver = []
        # create nodes
        nodeIdentifier = 1
        for n3 in range(len(x_list)):
            if n3+1 in nodes_w_ver:
                nt = nodetemplate1
            else:
                nt = nodetemplate
            node = nodes.createNode(nodeIdentifier, nt)
            cache.setNode(node)

            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x_list[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1_list[n3])
            if n3+1 in nodes_w_ver:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1_list[n3])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, d1_list[n3])

            radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 6)
            nodeIdentifier = nodeIdentifier + 1

        if vasculature:
            elem_list = [
                [2, 1],
                [3, 2],
                [4, 3],
                [5, 4],
                [3, 6],
                [4, 7],
                [5, 8],
                [8, 9],
                [8, 10],
                [1, 11],
                [11, 20],
                [20, 21],
                [20, 17],
                [17, 18],
                [17, 19],
                [11, 12],
                [12, 13],
                [17, 14],
                [14, 15],
                [15, 16],
                [22, 5],
            ]
        elif vagus:
            elem_list = [
                [1, 2],
                [2, 3],
                [3, 4],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 8],
                [8, 9],
                [9, 10],
                [10, 11],
                [11, 12],
                [12, 13],
                [13, 14],
                [14, 15],
                [3, 16],
                [3, 17],
                [4, 18],
                [18, 20],
                [20, 22],
                [18, 19],
                [20, 21],
                [22, 23],
                [22, 25],
                [22, 24],
                [5, 27],
                [2, 26],
                [27, 29],
                [27, 28],
                [6, 31],
                [2, 30],
                [7, 33],
                [2, 32],
                [33, 34],
                [33, 35],
                [34, 36],
                [34, 37],
                [37, 38],
                [37, 39],
                [37, 40],
                [37, 41],
                [36, 42],
                [36, 43],
                [36, 44],
                [36, 45],
                [36, 46],
                [35, 47],
                [35, 48],
                [35, 49],
                [35, 50],
                [35, 51],
                [35, 52],
                [8, 54],
                [36, 53],
                [54, 55],
                [55, 56],
                [54, 56],
                [9, 57],
                [57, 58],
                [57, 59],
                [57, 60],
                [57, 62],
                [57, 61],
                [57, 65],
                [57, 64],
                [57, 63],
                [57, 68],
                [59, 69],
                [59, 67],
                [59, 70],
                [59, 66],
                [59, 73],
                [59, 72],
                [59, 71],
                [3, 74],
                [3, 75],
                [2, 76],
                [10, 77],
                [11, 78],
                [12, 79],
                [79, 80],
                [79, 81],
                [13, 82],
                [14, 83],
                [15, 84],
                [59, 85],
            ]

        # create elements
        elementIdentifier = 1
        for e3 in range(len(elem_list)):
            if vasculature:
                if elementIdentifier in [7, 9, 6, 5, 16, 13, 14]:
                    elementtemplate = elementtemplate2
                    eft = eft2
                elif elementIdentifier in [8, 2, 11, 12, 15]:
                    elementtemplate = elementtemplate3
                    eft = eft3
                else:
                    elementtemplate = elementtemplate1
                    eft = eft1
            elif vagus:
                elementtemplate = elementtemplate1
                eft = eft1
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result = element.setNodesByIdentifier(eft, elem_list[e3])
            result = element.setNodesByIdentifier(eftRadius, elem_list[e3])
            elementIdentifier = elementIdentifier + 1

        vagusGroup = AnnotationGroup(region, ("n_5", ""))
        annotationGroups = [vagusGroup]

        # Groups of different parts of the vagus
        is_vagus = fm.createFieldConstant(1)
        vagusMeshGroup = vagusGroup.getMeshGroup(mesh)
        vagusMeshGroup.addElementsConditional(is_vagus)


        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCount1 = options['Refine number of elements 1']
        refineElementsCount2 = options['Refine number of elements 2']
        refineElementsCount3 = options['Refine number of elements 3']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount1, refineElementsCount2, refineElementsCount3)
