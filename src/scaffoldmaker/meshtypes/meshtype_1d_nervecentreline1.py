"""
Generates a 1-D mesh for nerve centrelines.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldCoordinates,\
    findOrCreateFieldFiniteElement, findOrCreateFieldStoredString
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation import nerve_terms


class MeshType_1d_nervecentreline1(Scaffold_base):
    '''
    Generates a 1-D mesh for nerve centrelines
    '''
    @staticmethod
    def getName():
        return '1D Nerve Centreline 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements 1' : 1,
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements 1',
        ]

    @staticmethod
    def checkOptions(options):
        for key in ['Number of elements 1']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base cubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        # elementsCount1 = options['Number of elements 1']
        # elementsCount2 = options['Number of elements 2']
        # elementsCount3 = options['Number of elements 3']
        # useCrossDerivatives = options['Use cross derivatives']

        # markerName.assignString(cache, x0name)

        # group
        # is_all = fm.createFieldConstant(1)
        # markerGroup = findOrCreateFieldGroup(fm, 'centreline')
        # markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        # markerPoints.addNodesConditional(is_all)
        #
        # elementGroup = findOrCreateFieldGroup(fm, 'centreline')
        # meshGroup = elementGroup.createFieldElementGroup(mesh).getMeshGroup()
        # meshGroup.addElementsConditional(is_all)

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        radius = findOrCreateFieldFiniteElement(fm, "radius", components_count=1, managed=True)
        markerName = findOrCreateFieldStoredString(fm, name="marker_name", managed=True)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.defineField(radius)
        nodetemplate.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.defineField(markerName)

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

        # markers group, coordinates and derivatives
        def marker_record(name, x, d1, radius_size):
            return {"group": nerve_terms.get_nerve_term(name), "x": x, "d1": d1, "r": radius_size}

        nerveMarkerPoints = [
            marker_record("brain_44-1", [-3.068, -58.888, 1485.941], [-5.129, -2.099, -5.383], 6.0),
            marker_record("brain_12-1", [0.526, -56.094, 1488.768], [-3.962, 3.013, 3.485], 6.0),
            marker_record("brain_34-1", [1.878, -56.261, 1489.3], [-3.416, 3.366, 3.014], 6.0),
            marker_record("brain_33-1", [2.607, -56.401, 1489.594], [-3.918, 3.78, 3.375], 6.0),
            marker_record("point_1", [-0.628, -62.693, 1481.347], [25.147, 6.953, 2.966], 6.0),
            marker_record("ardell_4_branching_point", [3.296, -68.465, 1338.18], [-0.341, -1.935, -12.395], 6.0),
            marker_record("ganglion_1-1", [19.14, -71.172, 1471.654], [3.23, -2.449, -7.906], 6.0),
            marker_record("bolser_X_1", [2.263, -68.435, 1345.906], [-55.339, 23.836, -191.394], 6.0),
            marker_record("ardell_3_branching_point", [0.501, -67.23, 1351.551], [-58.22, 25.91, -179.434], 6.0),
            marker_record("point_3", [17.463, -71.807, 1476.875], [-6.734, 0.139, 6.414], 6.0),
            marker_record("label_1-1", [19.059, -76.211, 1480.118], [-6.259, 0.228, 3.865], 6.0),
            marker_record("ganglion_2-1", [35.504, -78.359, 1463.385], [-1.057, -4.299, -8.305], 6.0),
            marker_record("point_6", [15.862, -66.876, 1461.22], [-16.113, -0.11, 1.188], 6.0),
            marker_record("point_7", [14.086, -75.333, 1461.325], [3.242, -4.337, -0.252], 6.0),
            marker_record("point_8", [4.653, -62.902, 1463.956], [-7.282, 8.493, 9.077], 6.0),
            marker_record("point_9", [12.273, -59.511, 1466.982], [7.606, -0.578, -1.022], 6.0),
            marker_record("point_10", [-6.224, -58.028, 1472.469], [-26.479, -6.517, 43.508], 6.0),
            marker_record("label_4-1", [60.913, -52.711, 1519.712], [113.397, 7.643, 25.844], 6.0),
            marker_record("label_3-1", [-57.275, -70.875, 1519.911], [-65.695, -23.698, 36.302], 6.0),
            marker_record("label_2-1", [-65.468, -70.375, 1520.633], [-80.971, -22.695, 37.571], 6.0),
            marker_record("point_11", [36.38, -80.588, 1457.154], [-0.118, -2.112, -6.641], 6.0),
            marker_record("point_13", [11.232, -116.557, 1452.669], [-17.74, -38.752, -6.999], 6.0),
            marker_record("label_6-1", [-14.0, -152.355, 1446.593], [-32.062, -32.179, -5.049], 6.0),
            marker_record("label_5-1", [9.582, -154.877, 1443.576], [12.923, -33.908, -10.012], 6.0),
            marker_record("point_14", [35.703, -82.095, 1448.654], [-6.773, -2.898, -16.891], 6.0),
            marker_record("cardio_4-1", [-3.468, -73.077, 1344.169], [-68.816, 21.122, -164.179], 6.0),
            marker_record("point_16", [34.972, -83.394, 1441.867], [0.902, 1.073, -8.921], 6.0),
            marker_record("point_18", [20.44, -92.164, 1438.379], [-5.057, -3.352, -3.807], 6.0),
            marker_record("point_19", [14.327, -97.349, 1437.414], [-5.856, -6.144, -4.017], 6.0),
            marker_record("point_20", [6.677, -105.285, 1436.5], [-6.059, -19.197, -11.926], 6.0),
            marker_record("label_7-1", [-0.539, -131.251, 1434.71], [-7.68, -30.025, 7.655], 6.0),
            marker_record("label_8-1", [9.361, -130.529, 1433.801], [9.952, -27.252, 5.685], 6.0),
            marker_record("respiratory_11-1", [1.823, -112.755, 1413.751], [-3.0, 3.499, -27.6], 6.0),
            marker_record("respiratory_5-1", [5.96, -114.822, 1414.415], [3.841, 0.102, -26.783], 6.0),
            marker_record("point_21", [13.106, -99.879, 1430.559], [-7.534, -5.737, -20.761], 6.0),
            marker_record("label_10-1", [2.098, -109.068, 1405.412], [-14.437, -12.602, -29.442], 6.0),
            marker_record("label_11-1", [0.755, -107.608, 1406.119], [-17.08, -9.671, -27.975], 6.0),
            marker_record("label_12-1", [0.392, -105.992, 1406.523], [-17.761, -6.441, -27.108], 6.0),
            marker_record("point_37", [12.367, -99.729, 1421.523], [1.399, 2.283, -5.055], 6.0),
            marker_record("bolser_X-5", [4.033, -107.383, 1404.8], [-10.611, -9.27, -30.754], 6.0),
            marker_record("bolser_X-4", [1.376, -109.143, 1406.395], [-15.831, -12.715, -27.404], 6.0),
            marker_record("external_laryngeal_n_branching_point", [19.268, -95.333, 1427.899], [-2.347, -4.675, -8.715],
                          6.0),
            marker_record("plexus_2-1", [4.832, -136.024, 1448.806], [-15.246, -52.924, 50.184], 6.0),
            marker_record("bolser_X-14", [7.465, -106.606, 1402.082], [-14.226, -9.153, -28.678], 6.0),
            marker_record("label_15-1", [6.818, -105.579, 1401.325], [-15.423, -7.075, -30.01], 6.0),
            marker_record("label_14-1", [6.623, -108.164, 1402.941], [-15.945, -12.281, -27.056], 6.0),
            marker_record("bolser_X-15", [8.188, -105.664, 1400.628], [-12.738, -7.253, -31.429], 6.0),
            marker_record("point_41", [19.826, -96.819, 1421.789], [0.141, 1.614, -5.649], 6.0),
            marker_record("point_22", [34.265, -84.451, 1435.143], [-3.332, -7.049, -20.25], 6.0),
            marker_record("point_24", [29.099, -94.64, 1421.734], [-3.226, -6.683, -10.8], 6.0),
            marker_record("point_25", [27.685, -97.76, 1415.52], [-1.761, -3.668, -5.999], 6.0),
            marker_record("plexus_3-1", [25.584, -101.952, 1409.817], [-2.883, -5.866, -8.726], 6.0),
            marker_record("point_26", [32.29, -92.737, 1396.093], [0.321, -5.982, -10.694], 6.0),
            marker_record("point_28", [30.962, -99.435, 1395.981], [-4.261, -2.824, 1.454], 6.0),
            marker_record("label_9-1", [12.293, -108.004, 1394.644], [-22.842, -7.903, -1.896], 6.0),
            marker_record("point_29", [30.621, -102.116, 1398.184], [-5.481, -5.016, 2.894], 6.0),
            marker_record("label_15-1", [15.735, -106.477, 1387.002], [-17.503, -0.767, -12.83], 6.0),
            marker_record("label_16-1", [16.942, -106.518, 1385.462], [-15.129, -0.845, -15.793], 6.0),
            marker_record("label_17-1", [18.305, -110.622, 1389.761], [-12.874, -9.06, -7.701], 6.0),
            marker_record("label_18-1", [17.812, -111.228, 1390.543], [-13.768, -10.203, -6.101], 6.0),
            marker_record("label_19-1", [20.305, -109.662, 1385.525], [-8.697, -6.995, -15.815], 6.0),
            marker_record("label_20-1", [19.183, -111.35, 1388.633], [-11.115, -10.505, -9.945], 6.0),
            marker_record("label_21-1", [18.88, -110.329, 1387.243], [-11.712, -8.465, -12.709], 6.0),
            marker_record("label_22-1", [19.598, -110.678, 1389.008], [-10.296, -9.173, -9.205], 6.0),
            marker_record("label_13-1", [12.808, -107.604, 1396.148], [-21.763, -7.091, 1.085], 6.0),
            marker_record("respiratory_8-1", [11.744, -108.913, 1398.904], [-23.743, -9.631, 6.505], 6.0),
            marker_record("respiratory_11-1", [12.66, -104.603, 1398.421], [-21.174, -1.106, 5.361], 6.0),
            marker_record("point_38", [30.626, -102.783, 1394.211], [-0.863, -1.635, -4.163], 6.0),
            marker_record("bolser_X-2", [13.545, -109.239, 1392.497], [-20.354, -10.352, -6.156], 6.0),
            marker_record("bolser_X-3", [13.402, -106.487, 1391.751], [-20.209, -4.791, -7.477], 6.0),
            marker_record("point_30", [28.921, -109.629, 1316.751], [-3.914, -3.308, -13.084], 6.0),
            marker_record("plexus_4-1", [20.571, -47.974, 1214.357], [-11.614, 110.759, -121.987], 6.0),
            marker_record("point_31", [28.321, -112.239, 1305.054], [0.365, -9.088, -15.532], 6.0),
            marker_record("respiratory_17-1", [13.336, -56.091, 1200.123], [-4.754, 7.795, -8.567], 6.0),
            marker_record("point_32", [6.058, -96.193, 1168.221], [-7.486, 1.904, -38.19], 6.0),
            marker_record("point_33", [4.126, -91.315, 1143.098], [-6.426, 7.025, -34.377], 6.0),
            marker_record("plexus_5-1", [-0.497, -88.117, 1055.124], [-5.067, -1.88, -108.785], 6.0),
            marker_record("digestive_12-1", [-0.849, -78.626, 1059.708], [-5.782, 17.092, -99.868], 6.0),
            marker_record("point_34", [4.49, -120.333, 982.412], [-5.29, 48.18, -94.925], 6.0),
            marker_record("digestive_9-1", [4.051, -80.134, 952.868], [-0.119, 8.718, 0.525], 6.0),
            marker_record("point_35", [-4.485, -76.724, 919.081], [2.729, 53.763, -37.471], 6.0),
            marker_record("plexus_6-1", [2.251, -72.056, 934.224], [1.907, -1.596, 6.514], 6.0),
            marker_record("point_36", [-9.096, -18.384, 858.438], [73.681, 23.211, -65.919], 6.0),
            marker_record("urinary_13-1", [84.559, -64.284, 856.328], [82.082, -83.081, 44.569], 6.0),
            marker_record("brain_D_1", [13.373, -107.842, 1393.326], [-20.648, -7.563, -4.5], 6.0),
            # {"group": get_term("brain_44-1"), "x": [-3.068, -58.888, 1485.941], "d1": [-5.129, -2.099, -5.383],
            #  "r": 6.0},
            # {"group": get_term("brain_12-1"), "x": [0.526, -56.094, 1488.768], "d1": [-3.962, 3.013, 3.485], "r": 6.0},
            # {"group": get_term("brain_34-1"), "x": [1.878, -56.261, 1489.3], "d1": [-3.416, 3.366, 3.014], "r": 6.0},
            # {"group": get_term("brain_33-1"), "x": [2.607, -56.401, 1489.594], "d1": [-3.918, 3.78, 3.375], "r": 6.0},
            # {"group": get_term("point_1"), "x": [-0.628, -62.693, 1481.347], "d1": [25.147, 6.953, 2.966], "r": 6.0},
            # {"group": get_term("ardell_4_branching_point"), "x": [3.296, -68.465, 1338.18], "d1": [-0.341, -1.935, -12.395],
            #  "r": 6.0},
            # {"group": get_term("ganglion_1-1"), "x": [19.14, -71.172, 1471.654], "d1": [3.23, -2.449, -7.906],
            #  "r": 6.0},
            # {"group": get_term("bolser_X_1"), "x": [2.263, -68.435, 1345.906], "d1": [-55.339, 23.836, -191.394],
            #  "r": 6.0},
            # {"group": get_term("ardell_3_branching_point"), "x": [0.501, -67.23, 1351.551], "d1": [-58.22, 25.91, -179.434],
            #  "r": 6.0},
            # {"group": get_term("point_3"), "x": [17.463, -71.807, 1476.875], "d1": [-6.734, 0.139, 6.414], "r": 6.0},
            # {"group": get_term("label_1-1"), "x": [19.059, -76.211, 1480.118], "d1": [-6.259, 0.228, 3.865], "r": 6.0},
            # {"group": get_term("ganglion_2-1"), "x": [35.504, -78.359, 1463.385], "d1": [-1.057, -4.299, -8.305],
            #  "r": 6.0},
            # {"group": get_term("point_6"), "x": [15.862, -66.876, 1461.22], "d1": [-16.113, -0.11, 1.188], "r": 6.0},
            # {"group": get_term("point_7"), "x": [14.086, -75.333, 1461.325], "d1": [3.242, -4.337, -0.252],
            #  "r": 6.0},
            # {"group": get_term("point_8"), "x": [4.653, -62.902, 1463.956], "d1": [-7.282, 8.493, 9.077], "r": 6.0},
            # {"group": get_term("point_9"), "x": [12.273, -59.511, 1466.982], "d1": [7.606, -0.578, -1.022],
            #  "r": 6.0},
            # {"group": get_term("point_10"), "x": [-6.224, -58.028, 1472.469], "d1": [-26.479, -6.517, 43.508],
            #  "r": 6.0},
            # {"group": get_term("label_4-1"), "x": [60.913, -52.711, 1519.712], "d1": [113.397, 7.643, 25.844],
            #  "r": 6.0},
            # {"group": get_term("label_3-1"), "x": [-57.275, -70.875, 1519.911], "d1": [-65.695, -23.698, 36.302],
            #  "r": 6.0},
            # {"group": get_term("label_2-1"), "x": [-65.468, -70.375, 1520.633], "d1": [-80.971, -22.695, 37.571],
            #  "r": 6.0},
            # {"group": get_term("point_11"), "x": [36.38, -80.588, 1457.154], "d1": [-0.118, -2.112, -6.641],
            #  "r": 6.0},
            # {"group": get_term("point_13"), "x": [11.232, -116.557, 1452.669], "d1": [-17.74, -38.752, -6.999],
            #  "r": 6.0},
            # {"group": get_term("label_6-1"), "x": [-14.0, -152.355, 1446.593], "d1": [-32.062, -32.179, -5.049],
            #  "r": 6.0},
            # {"group": get_term("label_5-1"), "x": [9.582, -154.877, 1443.576], "d1": [12.923, -33.908, -10.012],
            #  "r": 6.0},
            # {"group": get_term("point_14"), "x": [35.703, -82.095, 1448.654], "d1": [-6.773, -2.898, -16.891],
            #  "r": 6.0},
            # {"group": get_term("cardio_4-1"), "x": [-3.468, -73.077, 1344.169], "d1": [-68.816, 21.122, -164.179],
            #  "r": 6.0},
            # {"group": get_term("point_16"), "x": [34.972, -83.394, 1441.867], "d1": [0.902, 1.073, -8.921], "r": 6.0},
            # {"group": get_term("point_18"), "x": [20.44, -92.164, 1438.379], "d1": [-5.057, -3.352, -3.807],
            #  "r": 6.0},
            # {"group": get_term("point_19"), "x": [14.327, -97.349, 1437.414], "d1": [-5.856, -6.144, -4.017],
            #  "r": 6.0},
            # {"group": get_term("point_20"), "x": [6.677, -105.285, 1436.5], "d1": [-6.059, -19.197, -11.926],
            #  "r": 6.0},
            # {"group": get_term("label_7-1"), "x": [-0.539, -131.251, 1434.71], "d1": [-7.68, -30.025, 7.655],
            #  "r": 6.0},
            # {"group": get_term("label_8-1"), "x": [9.361, -130.529, 1433.801], "d1": [9.952, -27.252, 5.685],
            #  "r": 6.0},
            # {"group": get_term("respiratory_11-1"), "x": [1.823, -112.755, 1413.751], "d1": [-3.0, 3.499, -27.6], "r": 6.0},
            # {"group": get_term("respiratory_5-1"), "x": [5.96, -114.822, 1414.415], "d1": [3.841, 0.102, -26.783], "r": 6.0},
            # {"group": get_term("point_21"), "x": [13.106, -99.879, 1430.559], "d1": [-7.534, -5.737, -20.761],
            #  "r": 6.0},
            # {"group": get_term("label_10-1"), "x": [2.098, -109.068, 1405.412], "d1": [-14.437, -12.602, -29.442],
            #  "r": 6.0},
            # {"group": get_term("label_11-1"), "x": [0.755, -107.608, 1406.119], "d1": [-17.08, -9.671, -27.975],
            #  "r": 6.0},
            # {"group": get_term("label_12-1"), "x": [0.392, -105.992, 1406.523], "d1": [-17.761, -6.441, -27.108],
            #  "r": 6.0},
            # {"group": get_term("point_37"), "x": [12.367, -99.729, 1421.523], "d1": [1.399, 2.283, -5.055], "r": 6.0},
            # {"group": get_term("bolser_X-5"), "x": [4.033, -107.383, 1404.8], "d1": [-10.611, -9.27, -30.754],
            #  "r": 6.0},
            # {"group": get_term("bolser_X-4"), "x": [1.376, -109.143, 1406.395], "d1": [-15.831, -12.715, -27.404],
            #  "r": 6.0},
            # {"group": get_term("external_laryngeal_n_branching_point"), "x": [19.268, -95.333, 1427.899], "d1": [-2.347, -4.675, -8.715],
            #  "r": 6.0},
            # {"group": get_term("plexus_2-1"), "x": [4.832, -136.024, 1448.806], "d1": [-15.246, -52.924, 50.184],
            #  "r": 6.0},
            # {"group": get_term("bolser_X-14"), "x": [7.465, -106.606, 1402.082], "d1": [-14.226, -9.153, -28.678],
            #  "r": 6.0},
            # {"group": get_term("label_15-1"), "x": [6.818, -105.579, 1401.325], "d1": [-15.423, -7.075, -30.01],
            #  "r": 6.0},
            # {"group": get_term("label_14-1"), "x": [6.623, -108.164, 1402.941], "d1": [-15.945, -12.281, -27.056],
            #  "r": 6.0},
            # {"group": get_term("bolser_X-15"), "x": [8.188, -105.664, 1400.628], "d1": [-12.738, -7.253, -31.429],
            #  "r": 6.0},
            # {"group": get_term("point_41"), "x": [19.826, -96.819, 1421.789], "d1": [0.141, 1.614, -5.649], "r": 6.0},
            # {"group": get_term("point_22"), "x": [34.265, -84.451, 1435.143], "d1": [-3.332, -7.049, -20.25],
            #  "r": 6.0},
            # {"group": get_term("point_24"), "x": [29.099, -94.64, 1421.734], "d1": [-3.226, -6.683, -10.8], "r": 6.0},
            # {"group": get_term("point_25"), "x": [27.685, -97.76, 1415.52], "d1": [-1.761, -3.668, -5.999], "r": 6.0},
            # {"group": get_term("plexus_3-1"), "x": [25.584, -101.952, 1409.817], "d1": [-2.883, -5.866, -8.726],
            #  "r": 6.0},
            # {"group": get_term("point_26"), "x": [32.29, -92.737, 1396.093], "d1": [0.321, -5.982, -10.694],
            #  "r": 6.0},
            # {"group": get_term("point_28"), "x": [30.962, -99.435, 1395.981], "d1": [-4.261, -2.824, 1.454],
            #  "r": 6.0},
            # {"group": get_term("label_9-1"), "x": [12.293, -108.004, 1394.644], "d1": [-22.842, -7.903, -1.896],
            #  "r": 6.0},
            # {"group": get_term("point_29"), "x": [30.621, -102.116, 1398.184], "d1": [-5.481, -5.016, 2.894],
            #  "r": 6.0},
            # {"group": get_term("label_15-1"), "x": [15.735, -106.477, 1387.002], "d1": [-17.503, -0.767, -12.83],
            #  "r": 6.0},
            # {"group": get_term("label_16-1"), "x": [16.942, -106.518, 1385.462], "d1": [-15.129, -0.845, -15.793],
            #  "r": 6.0},
            # {"group": get_term("label_17-1"), "x": [18.305, -110.622, 1389.761], "d1": [-12.874, -9.06, -7.701],
            #  "r": 6.0},
            # {"group": get_term("label_18-1"), "x": [17.812, -111.228, 1390.543], "d1": [-13.768, -10.203, -6.101],
            #  "r": 6.0},
            # {"group": get_term("label_19-1"), "x": [20.305, -109.662, 1385.525], "d1": [-8.697, -6.995, -15.815],
            #  "r": 6.0},
            # {"group": get_term("label_20-1"), "x": [19.183, -111.35, 1388.633], "d1": [-11.115, -10.505, -9.945],
            #  "r": 6.0},
            # {"group": get_term("label_21-1"), "x": [18.88, -110.329, 1387.243], "d1": [-11.712, -8.465, -12.709],
            #  "r": 6.0},
            # {"group": get_term("label_22-1"), "x": [19.598, -110.678, 1389.008], "d1": [-10.296, -9.173, -9.205],
            #  "r": 6.0},
            # {"group": get_term("label_13-1"), "x": [12.808, -107.604, 1396.148], "d1": [-21.763, -7.091, 1.085],
            #  "r": 6.0},
            # {"group": get_term("respiratory_8-1"), "x": [11.744, -108.913, 1398.904], "d1": [-23.743, -9.631, 6.505],
            #  "r": 6.0},
            # {"group": get_term("respiratory_11-1"), "x": [12.66, -104.603, 1398.421], "d1": [-21.174, -1.106, 5.361],
            #  "r": 6.0},
            # {"group": get_term("point_38"), "x": [30.626, -102.783, 1394.211], "d1": [-0.863, -1.635, -4.163],
            #  "r": 6.0},
            # {"group": get_term("bolser_X-2"), "x": [13.545, -109.239, 1392.497], "d1": [-20.354, -10.352, -6.156],
            #  "r": 6.0},
            # {"group": get_term("bolser_X-3"), "x": [13.402, -106.487, 1391.751], "d1": [-20.209, -4.791, -7.477],
            #  "r": 6.0},
            # {"group": get_term("point_30"), "x": [28.921, -109.629, 1316.751], "d1": [-3.914, -3.308, -13.084],
            #  "r": 6.0},
            # {"group": get_term("plexus_4-1"), "x": [20.571, -47.974, 1214.357], "d1": [-11.614, 110.759, -121.987],
            #  "r": 6.0},
            # {"group": get_term("point_31"), "x": [28.321, -112.239, 1305.054], "d1": [0.365, -9.088, -15.532],
            #  "r": 6.0},
            # {"group": get_term("respiratory_17-1"), "x": [13.336, -56.091, 1200.123], "d1": [-4.754, 7.795, -8.567],
            #  "r": 6.0},
            # {"group": get_term("point_32"), "x": [6.058, -96.193, 1168.221], "d1": [-7.486, 1.904, -38.19], "r": 6.0},
            # {"group": get_term("point_33"), "x": [4.126, -91.315, 1143.098], "d1": [-6.426, 7.025, -34.377],
            #  "r": 6.0},
            # {"group": get_term("plexus_5-1"), "x": [-0.497, -88.117, 1055.124], "d1": [-5.067, -1.88, -108.785],
            #  "r": 6.0},
            # {"group": get_term("digestive_12-1"), "x": [-0.849, -78.626, 1059.708], "d1": [-5.782, 17.092, -99.868],
            #  "r": 6.0},
            # {"group": get_term("point_34"), "x": [4.49, -120.333, 982.412], "d1": [-5.29, 48.18, -94.925], "r": 6.0},
            # {"group": get_term("brain_44-1"), "x": [4.051, -80.134, 952.868], "d1": [-0.119, 8.718, 0.525], "r": 6.0},
            # {"group": get_term("point_35"), "x": [-4.485, -76.724, 919.081], "d1": [2.729, 53.763, -37.471],
            #  "r": 6.0},
            # {"group": get_term("plexus_6-1"), "x": [2.251, -72.056, 934.224], "d1": [1.907, -1.596, 6.514], "r": 6.0},
            # {"group": get_term("point_36"), "x": [-9.096, -18.384, 858.438], "d1": [73.681, 23.211, -65.919],
            #  "r": 6.0},
            # {"group": get_term("urinary_13-1"), "x": [84.559, -64.284, 856.328], "d1": [82.082, -83.081, 44.569],
            #  "r": 6.0},
            # {"group": get_term("brain_D_1"), "x": [13.373, -107.842, 1393.326], "d1": [-20.648, -7.563, -4.5],
            #  "r": 6.0},
            ]

        # x_list = [
        #     [-3.068, -58.888, 1485.941],
        #     [-0.628, -62.693, 1481.347],
        #     [19.140, -71.172, 1471.654],
        #     [35.504, -78.359, 1463.385],
        #     [36.380, -80.588, 1457.154],
        #     [35.703, -82.095, 1448.654],
        #     [34.972, -83.394, 1441.867],
        #     [34.265, -84.451, 1435.143],
        #     [32.290, -92.737, 1396.093],
        #     [28.921, -109.629, 1316.751],
        #     [28.321, -112.239, 1305.054],
        #     [6.058, -96.193, 1168.221],
        #     [4.490, -120.333, 982.412],
        #     [-4.485, -76.724, 919.081],
        #     [-9.096, -18.384, 858.438],
        #     [17.463, -71.807, 1476.875],
        #     [19.059, -76.211, 1480.118],
        #     [15.862, -66.876, 1461.220],
        #     [14.086, -75.333, 1461.325],
        #     [4.653, -62.902, 1463.956],
        #     [12.273, -59.511, 1466.982],
        #     [-6.224, -58.028, 1472.469],
        #     [-57.275, -70.875, 1519.911],
        #     [-65.468, -70.375, 1520.633],
        #     [60.913, -52.711, 1519.712],
        #     [0.526, -56.094, 1488.768],
        #     [11.232, -116.557, 1452.669],
        #     [9.582, -154.877, 1443.576],
        #     [-14.000, -152.355, 1446.593],
        #     [1.878, -56.261, 1489.300],
        #     [-3.468, -73.077, 1344.169],
        #     [2.607, -56.401, 1489.594],
        #     [20.440, -92.164, 1438.379],
        #     [14.327, -97.349, 1437.414],
        #     [19.268, -95.333, 1427.899],
        #     [13.106, -99.879, 1430.559],
        #     [6.677, -105.285, 1436.500],
        #     [-0.539, -131.251, 1434.710],
        #     [9.361, -130.529, 1433.801],
        #     [5.960, -114.822, 1414.415],
        #     [1.823, -112.755, 1413.751],
        #     [0.392, -105.992, 1406.523],
        #     [4.033, -107.383, 1404.800],
        #     [0.755, -107.608, 1406.119],
        #     [1.376, -109.143, 1406.395],
        #     [12.367, -99.729, 1421.523],
        #     [4.832, -136.024, 1448.806],
        #     [19.826, -96.819, 1421.789],
        #     [8.188, -105.664, 1400.628],
        #     [7.465, -106.606, 1402.082],
        #     [6.623, -108.164, 1402.941],
        #     [6.818, -105.579, 1401.325],
        #     [2.098, -109.068, 1405.412],
        #     [29.099, -94.640, 1421.734],
        #     [27.685, -97.760, 1415.520],
        #     [25.584, -101.952, 1409.817],
        #     [30.962, -99.435, 1395.981],
        #     [30.626, -102.783, 1394.211],
        #     [30.621, -102.116, 1398.184],
        #     [12.808, -107.604, 1396.148],
        #     [12.660, -104.603, 1398.421],
        #     [13.545, -109.239, 1392.497],
        #     [13.373, -107.842, 1393.326],
        #     [12.293, -108.004, 1394.644],
        #     [11.744, -108.913, 1398.904],
        #     [19.598, -110.678, 1389.008],
        #     [20.305, -109.662, 1385.525],
        #     [13.402, -106.487, 1391.751],
        #     [17.812, -111.228, 1390.543],
        #     [16.942, -106.518, 1385.462],
        #     [18.305, -110.622, 1389.761],
        #     [18.880, -110.329, 1387.243],
        #     [19.183, -111.350, 1388.633],
        #     [2.263, -68.435, 1345.906],
        #     [0.501, -67.230, 1351.551],
        #     [3.296, -68.465, 1338.180],
        #     [20.571, -47.974, 1214.357],
        #     [13.336, -56.091, 1200.123],
        #     [4.126, -91.315, 1143.098],
        #     [-0.849, -78.626, 1059.708],
        #     [-0.497, -88.117, 1055.124],
        #     [4.051, -80.134, 952.868],
        #     [2.251, -72.056, 934.224],
        #     [84.559, -64.284, 856.328],
        #     [15.735, -106.477, 1387.002],
        # ]
        #
        # d1_list = [
        #     [-5.129, -2.099, -5.383],
        #     [25.147, 6.953, 2.966],
        #     [3.230, -2.449, -7.906],
        #     [-1.057, -4.299, -8.305],
        #     [-0.118, -2.112, -6.641],
        #     [-6.773, -2.898, -16.891],
        #     [0.902, 1.073, -8.921],
        #     [-3.332, -7.049, -20.250],
        #     [0.321, -5.982, -10.694],
        #     [-3.914, -3.308, -13.084],
        #     [0.365, -9.088, -15.532],
        #     [-7.486, 1.904, -38.190],
        #     [-5.290, 48.180, -94.925],
        #     [2.729, 53.763, -37.471],
        #     [73.681, 23.211, -65.919],
        #     [-6.734, 0.139, 6.414],
        #     [-6.259, 0.228, 3.865],
        #     [-16.113, -0.110, 1.188],
        #     [3.242, -4.337, -0.252],
        #     [-7.282, 8.493, 9.077],
        #     [7.606, -0.578, -1.022],
        #     [-26.479, -6.517, 43.508],
        #     [-65.695, -23.698, 36.302],
        #     [-80.971, -22.695, 37.571],
        #     [113.397, 7.643, 25.844],
        #     [-3.962, 3.013, 3.485],
        #     [-17.740, -38.752, -6.999],
        #     [12.923, -33.908, -10.012],
        #     [-32.062, -32.179, -5.049],
        #     [-3.416, 3.366, 3.014],
        #     [-68.816, 21.122, -164.179],
        #     [-3.918, 3.780, 3.375],
        #     [-5.057, -3.352, -3.807],
        #     [-5.856, -6.144, -4.017],
        #     [-2.347, -4.675, -8.715],
        #     [-7.534, -5.737, -20.761],
        #     [-6.059, -19.197, -11.926],
        #     [-7.680, -30.025, 7.655],
        #     [9.952, -27.252, 5.685],
        #     [3.841, 0.102, -26.783],
        #     [-3.000, 3.499, -27.600],
        #     [-17.761, -6.441, -27.108],
        #     [-10.611, -9.270, -30.754],
        #     [-17.080, -9.671, -27.975],
        #     [-15.831, -12.715, -27.404],
        #     [1.399, 2.283, -5.055],
        #     [-15.246, -52.924, 50.184],
        #     [0.141, 1.614, -5.649],
        #     [-12.738, -7.253, -31.429],
        #     [-14.226, -9.153, -28.678],
        #     [-15.945, -12.281, -27.056],
        #     [-15.423, -7.075, -30.010],
        #     [-14.437, -12.602, -29.442],
        #     [-3.226, -6.683, -10.800],
        #     [-1.761, -3.668, -5.999],
        #     [-2.883, -5.866, -8.726],
        #     [-4.261, -2.824, 1.454],
        #     [-0.863, -1.635, -4.163],
        #     [-5.481, -5.016, 2.894],
        #     [-21.763, -7.091, 1.085],
        #     [-21.174, -1.106, 5.361],
        #     [-20.354, -10.352, -6.156],
        #     [-20.648, -7.563, -4.500],
        #     [-22.842, -7.903, -1.896],
        #     [-23.743, -9.631, 6.505],
        #     [-10.296, -9.173, -9.205],
        #     [-8.697, -6.995, -15.815],
        #     [-20.209, -4.791, -7.477],
        #     [-13.768, -10.203, -6.101],
        #     [-15.129, -0.845, -15.793],
        #     [-12.874, -9.060, -7.701],
        #     [-11.712, -8.465, -12.709],
        #     [-11.115, -10.505, -9.945],
        #     [-55.339, 23.836, -191.394],
        #     [-58.220, 25.910, -179.434],
        #     [-0.341, -1.935, -12.395],
        #     [-11.614, 110.759, -121.987],
        #     [-4.754, 7.795, -8.567],
        #     [-6.426, 7.025, -34.377],
        #     [-5.782, 17.092, -99.868],
        #     [-5.067, -1.880, -108.785],
        #     [-0.119, 8.718, 0.525],
        #     [1.907, -1.596, 6.514],
        #     [82.082, -83.081, 44.569],
        #     [-17.503, -0.767, -12.83],
        #     ]
        nodes_w_ver = []

        annotationGroups = []
        marker_node_dict = {}
        # create nodes
        nodeIdentifier = 1
        for n3 in range(len(nerveMarkerPoints)):
            if n3+1 in nodes_w_ver:
                nt = nodetemplate1
            else:
                nt = nodetemplate
            node = nodes.createNode(nodeIdentifier, nt)
            cache.setNode(node)

            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x_list[n3])
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1_list[n3])
            marker_node_dict[nerveMarkerPoints[n3]["group"][0]] = nodeIdentifier
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nerveMarkerPoints[n3]["x"])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nerveMarkerPoints[n3]["d1"])
            markerName.assignString(cache, nerveMarkerPoints[n3]["group"][0])
            annotationGroup = AnnotationGroup(region, nerveMarkerPoints[n3]["group"])
            annotationGroups.append(annotationGroup)
            annotationGroup.getNodesetGroup(nodes).addNode(node)
            if n3+1 in nodes_w_ver:
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1_list[n3])
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, d1_list[n3])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, nerveMarkerPoints[n3]["d1"])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, nerveMarkerPoints[n3]["d1"])

            radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nerveMarkerPoints[n3]["r"])
            nodeIdentifier = nodeIdentifier + 1

        # connectivity
        # elem_list = [
        #     [1, 2],
        #     [2, 3],
        #     [3, 4],
        #     [4, 5],
        #     [5, 6],
        #     [6, 7],
        #     [7, 8],
        #     [8, 9],
        #     [9, 10],
        #     [10, 11],
        #     [11, 12],
        #     [12, 13],
        #     [13, 14],
        #     [14, 15],
        #     [3, 16],
        #     [3, 17],
        #     [4, 18],
        #     [18, 20],
        #     [20, 22],
        #     [18, 19],
        #     [20, 21],
        #     [22, 23],
        #     [22, 25],
        #     [22, 24],
        #     [5, 27],
        #     [2, 26],
        #     [27, 29],
        #     [27, 28],
        #     [6, 31],
        #     [2, 30],
        #     [7, 33],
        #     [2, 32],
        #     [33, 34],
        #     [33, 35],
        #     [34, 36],
        #     [34, 37],
        #     [37, 38],
        #     [37, 39],
        #     [37, 40],
        #     [37, 41],
        #     [36, 42],
        #     [36, 43],
        #     [36, 44],
        #     [36, 45],
        #     [36, 46],
        #     [35, 47],
        #     [35, 48],
        #     [35, 49],
        #     [35, 50],
        #     [35, 51],
        #     [35, 52],
        #     [8, 54],
        #     [36, 53],
        #     [54, 55],
        #     [55, 56],
        #     [54, 56],
        #     [9, 57],
        #     [57, 58],
        #     [57, 59],
        #     [57, 60],
        #     [57, 62],
        #     [57, 61],
        #     [57, 65],
        #     [57, 64],
        #     [57, 63],
        #     [57, 68],
        #     [59, 69],
        #     [59, 67],
        #     [59, 70],
        #     [59, 66],
        #     [59, 73],
        #     [59, 72],
        #     [59, 71],
        #     [3, 74],
        #     [3, 75],
        #     [2, 76],
        #     [10, 77],
        #     [11, 78],
        #     [12, 79],
        #     [79, 80],
        #     [79, 81],
        #     [13, 82],
        #     [14, 83],
        #     [15, 84],
        #     [59, 85],
        # ]

        centrelines = [
            {'id': 'n_1', 'connects': ['brain_44-1', 'point_1']},
            {'id': 'n_5', 'connects': ['point_1', 'ganglion_1-1']},
            {'id': 'n_5', 'connects': ['ganglion_1-1', 'ganglion_2-1']},
            {'id': 'n_5', 'connects': ['ganglion_2-1', 'point_11']},
            {'id': 'n_5', 'connects': ['point_11', 'point_14']},
            {'id': 'n_5', 'connects': ['point_14', 'point_16']},
            {'id': 'n_5', 'connects': ['point_16', 'point_22']},
            {'id': 'n_5', 'connects': ['point_22', 'point_26']},
            {'id': 'n_5', 'connects': ['point_26', 'point_30']},
            {'id': 'n_5', 'connects': ['point_30', 'point_31']},
            {'id': 'n_5', 'connects': ['point_31', 'point_32']},
            {'id': 'n_5', 'connects': ['point_32', 'point_34']},
            {'id': 'n_5', 'connects': ['point_34', 'point_35']},
            {'id': 'n_5', 'connects': ['point_35', 'point_36']},
            {'id': 'n_68', 'connects': ['ganglion_1-1', 'point_3']},
            {'id': 'n_6', 'connects': ['ganglion_1-1', 'label_1-1']},
            {'id': 'n_7', 'connects': ['ganglion_2-1', 'point_6']},
            {'id': 'n_7', 'connects': ['point_6', 'point_8']},
            {'id': 'n_7', 'connects': ['point_8', 'point_10']},
            {'id': 'n_12', 'connects': ['point_6', 'point_7']},
            {'id': 'n_8', 'connects': ['point_8', 'point_9']},
            {'id': 'n_9', 'connects': ['point_10', 'label_3-1']},
            {'id': 'n_11', 'connects': ['point_10', 'label_4-1']},
            {'id': 'n_10', 'connects': ['point_10', 'label_2-1']},
            {'id': 'n_13', 'connects': ['point_11', 'point_13']},
            {'id': 'n_2', 'connects': ['point_1', 'brain_12-1']},
            {'id': 'n_15', 'connects': ['point_13', 'label_6-1']},
            {'id': 'n_14', 'connects': ['point_13', 'label_5-1']},
            {'id': 'n_16', 'connects': ['point_14', 'cardio_4-1']},
            {'id': 'n_3', 'connects': ['point_1', 'brain_34-1']},
            {'id': 'n_17', 'connects': ['point_16', 'point_18']},
            {'id': 'n_4', 'connects': ['point_1', 'brain_33-1']},
            {'id': 'n_18', 'connects': ['point_18', 'point_19']},
            {'id': 'n_29', 'connects': ['point_18', 'external_laryngeal_n_branching_point']},
            {'id': 'n_24', 'connects': ['point_19', 'point_21']},
            {'id': 'n_19', 'connects': ['point_19', 'point_20']},
            {'id': 'n_22', 'connects': ['point_20', 'label_7-1']},
            {'id': 'n_20', 'connects': ['point_20', 'label_8-1']},
            {'id': 'n_23', 'connects': ['point_20', 'respiratory_5-1']},
            {'id': 'n_21', 'connects': ['point_20', 'respiratory_11-1']},
            {'id': 'n_27', 'connects': ['point_21', 'label_12-1']},
            {'id': 'bolser_5', 'connects': ['point_21', 'bolser_X-5']},
            {'id': 'n_26', 'connects': ['point_21', 'label_11-1']},
            {'id': 'bolser_4', 'connects': ['point_21', 'bolser_X-4']},
            {'id': 'n_28', 'connects': ['point_21', 'point_37']},
            {'id': 'n_67', 'connects': ['external_laryngeal_n_branching_point', 'plexus_2-1']},
            {'id': 'n_68', 'connects': ['external_laryngeal_n_branching_point', 'point_41']},
            {'id': 'bolser_15', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-15']},
            {'id': 'bolser_14', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-14']},
            {'id': 'n_65', 'connects': ['external_laryngeal_n_branching_point', 'label_14-1']},
            {'id': 'n_66', 'connects': ['external_laryngeal_n_branching_point', 'label_15-1']},
            {'id': 'n_30', 'connects': ['point_22', 'point_24']},
            {'id': 'n_25', 'connects': ['point_21', 'label_10-1']},
            {'id': 'n_32', 'connects': ['point_24', 'point_25']},
            {'id': 'n_34', 'connects': ['point_25', 'plexus_3-1']},
            {'id': 'n_35', 'connects': ['point_24', 'plexus_3-1']},
            {'id': 'n_36', 'connects': ['point_26', 'point_28']},
            {'id': 'n_50', 'connects': ['point_28', 'point_38']},
            {'id': 'n_38', 'connects': ['point_28', 'point_29']},
            {'id': 'n_47', 'connects': ['point_28', 'label_13-1']},
            {'id': 'bolser_2', 'connects': ['point_28', 'bolser_X-2']},
            {'id': 'n_49', 'connects': ['point_28', 'respiratory_11-1']},
            {'id': 'n_48', 'connects': ['point_28', 'respiratory_8-1']},
            {'id': 'n_37', 'connects': ['point_28', 'label_9-1']},
            {'id': 'brain_D_1', 'connects': ['point_28', 'brain_D_1']},
            {'id': 'bolser_3', 'connects': ['point_28', 'bolser_X-3']},
            {'id': 'n_42', 'connects': ['point_29', 'label_18-1']},
            {'id': 'n_43', 'connects': ['point_29', 'label_19-1']},
            {'id': 'n_40', 'connects': ['point_29', 'label_16-1']},
            {'id': 'n_46', 'connects': ['point_29', 'label_22-1']},
            {'id': 'n_44', 'connects': ['point_29', 'label_20-1']},
            {'id': 'n_45', 'connects': ['point_29', 'label_21-1']},
            {'id': 'n_41', 'connects': ['point_29', 'label_17-1']},
            {'id': 'bolser_1', 'connects': ['ganglion_1-1', 'bolser_X_1']},
            {'id': 'ardell_3', 'connects': ['ganglion_1-1', 'ardell_3_branching_point']},
            {'id': 'ardell_4', 'connects': ['point_1', 'ardell_4_branching_point']},
            {'id': 'n_51', 'connects': ['point_30', 'plexus_4-1']},
            {'id': 'n_52', 'connects': ['point_31', 'respiratory_17-1']},
            {'id': 'n_53', 'connects': ['point_32', 'point_33']},
            {'id': 'n_54', 'connects': ['point_33', 'digestive_12-1']},
            {'id': 'n_55', 'connects': ['point_33', 'plexus_5-1']},
            {'id': 'n_56', 'connects': ['point_34', 'digestive_9-1']},
            {'id': 'n_57', 'connects': ['point_35', 'plexus_6-1']},
            {'id': 'n_58', 'connects': ['point_36', 'urinary_13-1']},
            {'id': 'n_39', 'connects': ['point_29', 'label_15-1']},
        ]

        # elem_list_marker = [
        #     ['brain_44-1', 'point_1'],
        #     ['point_1', 'ganglion_1-1'],
        #     ['ganglion_1-1', 'ganglion_2-1'],
        #     ['ganglion_2-1', 'point_11'],
        #     ['point_11', 'point_14'],
        #     ['point_14', 'point_16'],
        #     ['point_16', 'point_22'],
        #     ['point_22', 'point_26'],
        #     ['point_26', 'point_30'],
        #     ['point_30', 'point_31'],
        #     ['point_31', 'point_32'],
        #     ['point_32', 'point_34'],
        #     ['point_34', 'point_35'],
        #     ['point_35', 'point_36'],
        #     ['ganglion_1-1', 'point_3'],
        #     ['ganglion_1-1', 'label_1-1'],
        #     ['ganglion_2-1', 'point_6'],
        #     ['point_6', 'point_8'],
        #     ['point_8', 'point_10'],
        #     ['point_6', 'point_7'],
        #     ['point_8', 'point_9'],
        #     ['point_10', 'label_3-1'],
        #     ['point_10', 'label_4-1'],
        #     ['point_10', 'label_2-1'],
        #     ['point_11', 'point_13'],
        #     ['point_1', 'brain_12-1'],
        #     ['point_13', 'label_6-1'],
        #     ['point_13', 'label_5-1'],
        #     ['point_14', 'cardio_4-1'],
        #     ['point_1', 'brain_34-1'],
        #     ['point_16', 'point_18'],
        #     ['point_1', 'brain_33-1'],
        #     ['point_18', 'point_19'],
        #     ['point_18', 'external_laryngeal_n_branching_point'],
        #     ['point_19', 'point_21'],
        #     ['point_19', 'point_20'],
        #     ['point_20', 'label_7-1'],
        #     ['point_20', 'label_8-1'],
        #     ['point_20', 'respiratory_5-1'],
        #     ['point_20', 'respiratory_11-1'],
        #     ['point_21', 'label_12-1'],
        #     ['point_21', 'bolser_X-5'],
        #     ['point_21', 'label_11-1'],
        #     ['point_21', 'bolser_X-4'],
        #     ['point_21', 'point_37'],
        #     ['external_laryngeal_n_branching_point', 'plexus_2-1'],
        #     ['external_laryngeal_n_branching_point', 'point_41'],
        #     ['external_laryngeal_n_branching_point', 'bolser_X-15'],
        #     ['external_laryngeal_n_branching_point', 'bolser_X-14'],
        #     ['external_laryngeal_n_branching_point', 'label_14-1'],
        #     ['external_laryngeal_n_branching_point', 'label_15-1'],
        #     ['point_22', 'point_24'],
        #     ['point_21', 'label_10-1'],
        #     ['point_24', 'point_25'],
        #     ['point_25', 'plexus_3-1'],
        #     ['point_24', 'plexus_3-1'],
        #     ['point_26', 'point_28'],
        #     ['point_28', 'point_38'],
        #     ['point_28', 'point_29'],
        #     ['point_28', 'label_13-1'],
        #     ['point_28', 'bolser_X-2'],
        #     ['point_28', 'respiratory_11-1'],
        #     ['point_28', 'respiratory_8-1'],
        #     ['point_28', 'label_9-1'],
        #     ['point_28', 'brain_D_1'],
        #     ['point_28', 'bolser_X-3'],
        #     ['point_29', 'label_18-1'],
        #     ['point_29', 'label_19-1'],
        #     ['point_29', 'label_16-1'],
        #     ['point_29', 'label_22-1'],
        #     ['point_29', 'label_20-1'],
        #     ['point_29', 'label_21-1'],
        #     ['point_29', 'label_17-1'],
        #     ['ganglion_1-1', 'bolser_X_1'],
        #     ['ganglion_1-1', 'ardell_3_branching_point'],
        #     ['point_1', 'ardell_4_branching_point'],
        #     ['point_30', 'plexus_4-1'],
        #     ['point_31', 'respiratory_17-1'],
        #     ['point_32', 'point_33'],
        #     ['point_33', 'digestive_12-1'],
        #     ['point_33', 'plexus_5-1'],
        #     ['point_34', 'brain_44-1'],
        #     ['point_35', 'plexus_6-1'],
        #     ['point_36', 'urinary_13-1'],
        #     ['point_29', 'label_15-1'],
        # ]

        # create elements
        # vagusGroup = AnnotationGroup(region, ("n_5", ""))
        # annotationGroups = [vagusGroup]
        # vagusMeshGroup = vagusGroup.getMeshGroup(mesh)


        elementIdentifier = 1
        for e3 in range(len(centrelines)):
            elementtemplate = elementtemplate1
            eft = eft1
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result = element.setNodesByIdentifier(eft, [marker_node_dict[centrelines[e3]['connects'][0]], marker_node_dict[centrelines[e3]['connects'][1]]])
            result = element.setNodesByIdentifier(eftRadius, [marker_node_dict[centrelines[e3]['connects'][0]], marker_node_dict[centrelines[e3]['connects'][1]]])
            annotationGroup = AnnotationGroup(region, nerve_terms.get_nerve_term(centrelines[e3]['id']))
            annotationGroups.append(annotationGroup)
            annotationMeshGroup = annotationGroup.getMeshGroup(mesh)
            annotationMeshGroup.addElement(element)
            elementIdentifier = elementIdentifier + 1

        # Groups of different parts of the vagus
        # is_vagus = fm.createFieldConstant(1)
        # vagusMeshGroup = vagusGroup.getMeshGroup(mesh)
        # vagusMeshGroup.addElementsConditional(is_vagus)

        # is_all = fm.createFieldConstant(1)
        # markerGroup = findOrCreateFieldGroup(fm, 'centreline')
        # markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        # markerPoints.addNodesConditional(is_all)
        #
        # elementGroup = findOrCreateFieldGroup(fm, 'centreline')
        # meshGroup = elementGroup.createFieldElementGroup(mesh).getMeshGroup()
        # meshGroup.addElementsConditional(is_all)


        fm.endChange()
        return annotationGroups
