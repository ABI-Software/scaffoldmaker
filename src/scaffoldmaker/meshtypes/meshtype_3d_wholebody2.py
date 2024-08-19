"""
Generates a 3D body coordinates using tube network mesh.
"""

from cmlibs.zinc.element import Element
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from cmlibs.zinc.node import Node


def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames()
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            "scaffoldSettings": {
                "Structure": "1-2-3, 3-4,4.2-5,4.3-6,4.1-7,7-8,8-9.1,9.2-10,9.3-11, \
                     5-12-13-14,6-15-16-17,10-18-19-20, 11-21-22-23",
                "Define inner coordinates": True
            },
            "meshEdits": exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[0.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, -0.000]]),
                    (2, [[0.500, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, -0.000]]),
                    (3, [[1.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, 0.100, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, 0.100]]),
                    (4, [[1.500, 0.000, 0.000], [[0.643, 0.000, 0.000], [0.028, -1.246, 0.000], [0.028, 1.246, 0.000]], [[0.000, 0.600, 0.000], [0.396, 0.009, 0.000], [-0.395, 0.009, 0.000]], [[0.000, 0.100, 0.000], [-0.374, 0.001, 0.000], [0.376, 0.005, 0.000]], [[0.000, 0.000, 0.600], [-0.000, 0.000, 0.400], [0.000, -0.000, 0.400]], [[0.000, 0.000, 0.100], [0.000, 0.000, -0.270], [0.000, 0.000, -0.268]]]),
                    (5, [[1.788, -1.014, 0.000], [0.546, -0.709, 0.000], [0.121, 0.093, 0.000], [-0.189, 0.048, 0.000], [-0.000, 0.000, 0.200], [0.000, 0.000, -0.130]]),
                    (6, [[1.788, 1.014, 0.000], [0.547, 0.709, 0.000], [-0.120, 0.093, 0.000], [0.186, 0.045, 0.000], [0.000, -0.000, 0.201], [0.000, 0.000, -0.130]]),
                    (7, [[2.400, -0.000, 0.000], [0.788, 0.000, 0.000], [0.000, 0.600, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.600], [0.000, 0.000, -0.000]]),
                    (8, [[3.100, -0.000, 0.000], [0.747, 0.000, 0.000], [0.000, 0.600, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.600], [0.000, 0.000, -0.000]]),
                    (9, [[3.900, -0.000, 0.000], [[0.853, 0.000, 0.000], [0.731, -0.720, 0.000], [0.731, 0.720, 0.000]], [[0.000, 0.600, 0.000], [0.496, 0.504, 0.000], [-0.470, 0.478, 0.000]], [[-0.000, 0.000, -0.000], [-0.839, -0.173, -0.006], [0.775, -0.153, 0.006]], [[0.000, 0.000, 0.600], [-0.000, 0.000, 0.350], [0.000, -0.000, 0.350]], [[-0.000, -0.000, 0.000], [0.006, 0.004, -0.022], [0.006, -0.004, -0.022]]]),
                    (10, [[5.250, -0.600, 0.000], [1.915, -0.427, 0.000], [0.056, 0.249, 0.000], [-0.211, -0.182, 0.006], [-0.000, 0.000, 0.300], [-0.003, -0.005, -0.078]]),
                    (11, [[5.250, 0.600, 0.000], [1.915, 0.427, 0.000], [-0.056, 0.249, 0.000], [0.206, -0.165, -0.006], [0.000, -0.000, 0.300], [-0.003, 0.005, -0.078]]),
                    (12, [[2.439, -1.399, 0.000], [0.920, -0.382, 0.000], [0.057, 0.137, 0.000], [-0.009, -0.011, 0.001], [-0.000, 0.000, 0.140], [0.006, -0.002, -0.030]]),
                    (13, [[3.807, -1.755, 0.000], [1.181, -0.272, 0.078], [0.028, 0.121, -0.000], [-0.018, -0.063, 0.000], [-0.009, 0.002, 0.140], [-0.023, 0.005, -0.056]]),
                    (14, [[4.841, -1.961, 0.137], [0.882, -0.141, 0.194], [0.002, 0.015, -0.000], [-0.023, -0.149, 0.002], [-0.006, 0.001, 0.029], [0.054, -0.011, -0.158]]),
                    (15, [[2.439, 1.399, 0.000], [0.946, 0.385, 0.000], [-0.056, 0.138, 0.000], [0.012, -0.011, -0.001], [0.000, -0.000, 0.140], [0.006, 0.002, -0.031]]),
                    (16, [[3.911, 1.761, 0.000], [1.211, 0.255, 0.091], [-0.025, 0.121, 0.000], [0.018, -0.056, -0.000], [-0.010, -0.002, 0.140], [-0.028, -0.005, -0.042]]),
                    (17, [[4.932, 1.939, 0.153], [0.825, 0.101, 0.214], [-0.004, 0.029, 0.000], [0.014, -0.129, -0.002], [-0.015, -0.002, 0.058], [0.044, 0.008, -0.114]]),
                    (18, [[8.100, -0.665, 0.000], [0.775, -0.027, 0.275], [0.007, 0.250, 0.003], [0.032, -0.000, 0.006], [-0.069, -0.000, 0.194], [-0.224, -0.006, -0.056]]),
                    (19, [[8.459, -0.683, 0.298], [0.258, -0.010, 0.370], [0.004, 0.250, 0.004], [-0.008, 0.000, -0.013], [-0.181, 0.001, 0.126], [0.051, -0.005, -0.001]]),
                    (20, [[8.601, -0.685, 0.694], [0.024, 0.006, 0.399], [-0.001, 0.250, -0.004], [-0.002, -0.000, 0.005], [-0.141, -0.000, 0.009], [0.109, 0.005, -0.120]]),
                    (21, [[8.100, 0.665, 0.000], [0.775, 0.027, 0.275], [-0.007, 0.250, -0.003], [-0.032, -0.000, -0.009], [-0.069, 0.000, 0.194], [-0.224, 0.008, -0.055]]),
                    (22, [[8.459, 0.683, 0.298], [0.258, 0.010, 0.370], [0.000, 0.250, -0.007], [0.008, 0.000, 0.013], [-0.182, 0.004, 0.127], [0.049, 0.005, 0.002]]),
                    (23, [[8.601, 0.685, 0.694], [0.024, -0.006, 0.399], [0.001, 0.250, 0.004], [-0.007, 0.000, -0.005], [-0.148, 0.000, 0.009], [0.101, -0.015, -0.124]])], [
                    (1, [[0.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, 0.000, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, 0.000]]),
                    (2, [[0.500, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, -0.000]]),
                    (3, [[1.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, 0.070, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, 0.070]]),
                    (4, [[1.500, 0.000, 0.000], [[0.643, 0.000, 0.000], [0.028, -1.246, 0.000], [0.028, 1.246, 0.000]], [[0.000, 0.420, 0.000], [0.277, 0.006, 0.000], [-0.277, 0.006, 0.000]], [[0.000, 0.070, 0.000], [-0.262, 0.001, 0.000], [0.263, 0.003, 0.000]], [[0.000, 0.000, 0.420], [-0.000, 0.000, 0.280], [0.000, -0.000, 0.280]], [[0.000, 0.000, 0.070], [0.000, 0.000, -0.189], [0.000, 0.000, -0.188]]]),
                    (5, [[1.788, -1.014, 0.000], [0.546, -0.709, 0.000], [0.085, 0.065, 0.000], [-0.133, 0.033, 0.000], [-0.000, 0.000, 0.140], [0.000, 0.000, -0.091]]),
                    (6, [[1.788, 1.014, 0.000], [0.547, 0.709, 0.000], [-0.084, 0.065, 0.000], [0.130, 0.032, 0.000], [0.000, -0.000, 0.141], [0.000, 0.000, -0.091]]),
                    (7, [[2.400, -0.000, 0.000], [0.788, 0.000, 0.000], [0.000, 0.420, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.420], [0.000, 0.000, 0.000]]),
                    (8, [[3.100, -0.000, 0.000], [0.747, 0.000, 0.000], [0.000, 0.420, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.420], [0.000, 0.000, -0.000]]),
                    (9, [[3.900, -0.000, 0.000], [[0.853, 0.000, 0.000], [0.731, -0.720, 0.000], [0.731, 0.720, 0.000]], [[0.000, 0.420, 0.000], [0.347, 0.353, 0.000], [-0.329, 0.334, 0.000]], [[-0.000, 0.000, -0.000], [-0.587, -0.121, -0.004], [0.543, -0.107, 0.004]], [[0.000, 0.000, 0.420], [-0.000, 0.000, 0.245], [0.000, -0.000, 0.245]], [[-0.000, -0.000, 0.000], [0.004, 0.003, -0.015], [0.004, -0.003, -0.015]]]),
                    (10, [[5.250, -0.600, 0.000], [1.915, -0.427, 0.000], [0.039, 0.174, 0.000], [-0.148, -0.127, 0.004], [-0.000, 0.000, 0.210], [-0.002, -0.004, -0.055]]),
                    (11, [[5.250, 0.600, 0.000], [1.915, 0.427, 0.000], [-0.039, 0.174, 0.000], [0.144, -0.115, -0.004], [0.000, -0.000, 0.210], [-0.002, 0.004, -0.055]]),
                    (12, [[2.439, -1.399, 0.000], [0.920, -0.382, 0.000], [0.040, 0.096, 0.000], [-0.006, -0.008, 0.001], [-0.000, 0.000, 0.098], [0.004, -0.002, -0.021]]),
                    (13, [[3.807, -1.755, 0.000], [1.181, -0.272, 0.078], [0.019, 0.084, -0.000], [-0.013, -0.044, 0.000], [-0.006, 0.001, 0.098], [-0.016, 0.003, -0.039]]),
                    (14, [[4.841, -1.961, 0.137], [0.882, -0.141, 0.194], [0.002, 0.011, -0.000], [-0.016, -0.104, 0.001], [-0.004, 0.001, 0.021], [0.038, -0.008, -0.110]]),
                    (15, [[2.439, 1.399, 0.000], [0.946, 0.385, 0.000], [-0.039, 0.096, 0.000], [0.008, -0.007, -0.001], [0.000, -0.000, 0.098], [0.004, 0.002, -0.022]]),
                    (16, [[3.911, 1.761, 0.000], [1.211, 0.255, 0.091], [-0.018, 0.085, 0.000], [0.013, -0.039, -0.000], [-0.007, -0.002, 0.098], [-0.019, -0.004, -0.029]]),
                    (17, [[4.932, 1.939, 0.153], [0.825, 0.101, 0.214], [-0.003, 0.020, 0.000], [0.010, -0.090, -0.002], [-0.010, -0.002, 0.041], [0.031, 0.006, -0.080]]),
                    (18, [[8.100, -0.665, 0.000], [0.775, -0.027, 0.275], [0.005, 0.175, 0.002], [0.022, -0.000, 0.004], [-0.048, -0.000, 0.136], [-0.157, -0.004, -0.039]]),
                    (19, [[8.459, -0.683, 0.298], [0.258, -0.010, 0.370], [0.003, 0.175, 0.003], [-0.005, 0.000, -0.009], [-0.127, 0.001, 0.088], [0.036, -0.004, -0.000]]),
                    (20, [[8.601, -0.685, 0.694], [0.024, 0.006, 0.399], [-0.001, 0.175, -0.002], [-0.001, -0.000, 0.003], [-0.099, -0.000, 0.006], [0.077, 0.004, -0.084]]),
                    (21, [[8.100, 0.665, 0.000], [0.775, 0.027, 0.275], [-0.005, 0.175, -0.002], [-0.022, -0.000, -0.006], [-0.048, 0.000, 0.136], [-0.157, 0.006, -0.039]]),
                    (22, [[8.459, 0.683, 0.298], [0.258, 0.010, 0.370], [0.000, 0.175, -0.005], [0.005, 0.000, 0.009], [-0.127, 0.002, 0.089], [0.034, 0.004, 0.001]]),
                    (23, [[8.601, 0.685, 0.694], [0.024, -0.006, 0.399], [0.001, 0.175, 0.002], [-0.005, 0.000, -0.003], [-0.104, 0.000, 0.006], [0.071, -0.010, -0.087]])
                ]]),
            "userAnnotationGroups": [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_body_term('head')[0],
                    'ontId': get_body_term('head')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-8',
                    'name': 'torso',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4,11-13, 5,14-16',
                    'name': get_body_term('arm')[0],
                    'ontId': get_body_term('arm')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9,17-18-19, 10,20-21-22',
                    'name': get_body_term('leg')[0],
                    'ontId': get_body_term('leg')[1]
                }
            ]
        })
    elif parameterSetName in ("Human 2 Coarse", "Human 2 Medium", "Human 2 Fine"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            "scaffoldSettings": {
                "Structure": "1-2-3, 3-4,4.2-5,4.3-6,4.1-7,7-8,8-9.1,9.2-10,9.3-11, \
                     5-12-13-14,6-15-16-17,10-18-19-20, 11-21-22-23",
                "Define inner coordinates": True
            },
            "meshEdits": exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[0.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, -0.000]]),
                    (2, [[0.500, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, -0.000]]),
                    (3, [[1.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.400, 0.000], [0.000, 0.100, 0.000], [0.000, 0.000, 0.400], [0.000, 0.000, 0.100]]),
                    (4, [[1.500, 0.000, 0.000], [[0.643, 0.000, 0.000], [0.012, -1.338, 0.000], [0.054, 1.358, 0.000]], [[0.000, 0.600, 0.000], [0.277, 0.003, 0.000], [-0.277, 0.011, 0.000]], [[0.000, 0.100, 0.000], [-0.083, -0.206, 0.000], [0.082, -0.198, 0.000]], [[0.000, 0.000, 0.600], [-0.000, 0.000, 0.280], [0.000, -0.000, 0.280]], [[0.000, 0.000, 0.100], [0.000, 0.000, -0.195], [0.000, 0.000, -0.195]]]),
                    (5, [[1.806, -1.023, 0.000], [0.599, -0.599, 0.000], [0.158, 0.158, 0.000], [-0.253, 0.182, 0.000], [-0.000, 0.000, 0.150], [0.000, 0.000, -0.065]]),
                    (6, [[1.825, 1.023, 0.000], [0.592, 0.587, 0.000], [-0.156, 0.158, 0.000], [0.242, 0.173, 0.000], [0.000, -0.000, 0.150], [0.000, 0.000, -0.065]]),
                    (7, [[2.400, -0.000, 0.000], [0.788, 0.000, 0.000], [0.000, 0.600, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.600], [0.000, 0.000, -0.000]]),
                    (8, [[3.100, -0.000, 0.000], [0.747, 0.000, 0.000], [0.000, 0.600, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.600], [0.000, 0.000, -0.000]]),
                    (9, [[3.900, -0.000, 0.000], [[0.853, 0.000, 0.000], [0.731, -0.720, 0.000], [0.731, 0.720, 0.000]], [[0.000, 0.600, 0.000], [0.496, 0.504, 0.000], [-0.470, 0.478, 0.000]], [[-0.000, 0.000, -0.000], [-0.839, -0.173, -0.006], [0.776, -0.153, 0.006]], [[0.000, 0.000, 0.600], [-0.000, 0.000, 0.350], [0.000, -0.000, 0.350]], [[-0.000, -0.000, 0.000], [0.006, 0.004, -0.022], [0.006, -0.004, -0.022]]]),
                    (10, [[5.250, -0.600, 0.000], [1.915, -0.427, 0.000], [0.056, 0.249, 0.000], [-0.211, -0.182, 0.006], [-0.000, 0.000, 0.300], [-0.003, -0.005, -0.078]]),
                    (11, [[5.250, 0.600, 0.000], [1.915, 0.427, 0.000], [-0.056, 0.249, 0.000], [0.206, -0.165, -0.006], [0.000, -0.000, 0.300], [-0.003, 0.005, -0.078]]),
                    (12, [[2.450, -1.213, 0.000], [1.000, -0.145, 0.000], [0.022, 0.151, 0.000], [0.031, -0.007, 0.000], [-0.000, 0.000, 0.150], [0.000, 0.000, -0.015]]),
                    (13, [[4.327, -1.224, 0.000], [0.756, 0.003, 0.000], [-0.001, 0.200, 0.000], [-0.001, 0.075, -0.001], [0.000, -0.000, 0.120], [-0.000, 0.000, -0.025]]),
                    (14, [[4.800, -1.218, 0.000], [0.190, 0.009, 0.000], [-0.015, 0.300, -0.001], [-0.045, 0.124, -0.001], [-0.000, 0.000, 0.100], [-0.000, 0.000, -0.015]]),
                    (15, [[2.450, 1.213, 0.000], [0.979, 0.146, 0.000], [-0.023, 0.151, 0.000], [-0.029, -0.008, 0.000], [0.000, -0.000, 0.150], [0.000, 0.000, -0.015]]),
                    (16, [[4.327, 1.224, 0.000], [0.756, -0.003, 0.000], [0.001, 0.200, 0.000], [0.001, 0.075, 0.001], [-0.000, 0.000, 0.120], [-0.000, -0.000, -0.025]]),
                    (17, [[4.800, 1.218, 0.000], [0.190, -0.009, 0.000], [0.015, 0.300, 0.001], [0.045, 0.124, 0.001], [-0.000, -0.000, 0.100], [-0.000, -0.000, -0.015]]),
                    (18, [[8.100, -0.665, 0.000], [0.775, -0.027, 0.275], [0.007, 0.250, 0.003], [0.032, -0.000, 0.006], [-0.069, -0.000, 0.194], [-0.224, -0.006, -0.056]]),
                    (19, [[8.459, -0.683, 0.298], [0.258, -0.010, 0.370], [0.004, 0.250, 0.004], [-0.008, 0.000, -0.013], [-0.181, 0.001, 0.126], [0.051, -0.006, -0.001]]),
                    (20, [[8.601, -0.685, 0.694], [0.024, 0.006, 0.399], [-0.001, 0.250, -0.004], [-0.003, -0.000, 0.005], [-0.141, -0.000, 0.009], [0.111, 0.004, -0.120]]),
                    (21, [[8.100, 0.665, 0.000], [0.775, 0.027, 0.275], [-0.007, 0.250, -0.003], [-0.032, -0.000, -0.008], [-0.069, 0.000, 0.194], [-0.224, 0.008, -0.055]]),
                    (22, [[8.459, 0.683, 0.298], [0.258, 0.010, 0.370], [-0.000, 0.250, -0.007], [0.008, 0.000, 0.013], [-0.182, 0.003, 0.127], [0.049, 0.005, 0.002]]),
                    (23, [[8.601, 0.685, 0.694], [0.024, -0.006, 0.399], [0.001, 0.250, 0.004], [-0.006, 0.000, -0.005], [-0.148, 0.000, 0.009], [0.101, -0.013, -0.124]])], [
                    (1, [[0.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, 0.000, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, 0.000]]),
                    (2, [[0.500, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, -0.000]]),
                    (3, [[1.000, 0.000, 0.000], [0.500, 0.000, 0.000], [0.000, 0.280, 0.000], [0.000, 0.070, 0.000], [0.000, 0.000, 0.280], [0.000, 0.000, 0.070]]),
                    (4, [[1.500, 0.000, 0.000], [[0.643, 0.000, 0.000], [0.012, -1.338, 0.000], [0.054, 1.358, 0.000]], [[0.000, 0.420, 0.000], [0.194, 0.002, 0.000], [-0.194, 0.008, 0.000]], [[0.000, 0.070, 0.000], [-0.058, -0.144, 0.000], [0.057, -0.139, 0.000]], [[0.000, 0.000, 0.420], [-0.000, 0.000, 0.196], [0.000, -0.000, 0.196]], [[0.000, 0.000, 0.070], [0.000, 0.000, -0.137], [0.000, 0.000, -0.137]]]),
                    (5, [[1.806, -1.023, 0.000], [0.599, -0.599, 0.000], [0.111, 0.111, 0.000], [-0.177, 0.128, 0.000], [-0.000, 0.000, 0.105], [0.000, 0.000, -0.045]]),
                    (6, [[1.825, 1.023, 0.000], [0.592, 0.587, 0.000], [-0.109, 0.110, 0.000], [0.169, 0.121, 0.000], [0.000, -0.000, 0.105], [0.000, 0.000, -0.045]]),
                    (7, [[2.400, -0.000, 0.000], [0.788, 0.000, 0.000], [0.000, 0.420, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.420], [0.000, 0.000, 0.000]]),
                    (8, [[3.100, -0.000, 0.000], [0.747, 0.000, 0.000], [0.000, 0.420, 0.000], [0.000, -0.000, 0.000], [0.000, 0.000, 0.420], [0.000, 0.000, -0.000]]),
                    (9, [[3.900, -0.000, 0.000], [[0.853, 0.000, 0.000], [0.731, -0.720, 0.000], [0.731, 0.720, 0.000]], [[0.000, 0.420, 0.000], [0.347, 0.353, 0.000], [-0.329, 0.335, 0.000]], [[-0.000, 0.000, -0.000], [-0.587, -0.121, -0.004], [0.543, -0.107, 0.004]], [[0.000, 0.000, 0.420], [-0.000, 0.000, 0.245], [0.000, -0.000, 0.245]], [[-0.000, -0.000, 0.000], [0.004, 0.003, -0.015], [0.004, -0.003, -0.015]]]),
                    (10, [[5.250, -0.600, 0.000], [1.915, -0.427, 0.000], [0.039, 0.174, 0.000], [-0.147, -0.127, 0.004], [-0.000, 0.000, 0.210], [-0.002, -0.004, -0.055]]),
                    (11, [[5.250, 0.600, 0.000], [1.915, 0.427, 0.000], [-0.039, 0.174, 0.000], [0.144, -0.115, -0.004], [0.000, -0.000, 0.210], [-0.002, 0.004, -0.055]]),
                    (12, [[2.450, -1.213, 0.000], [1.000, -0.145, 0.000], [0.015, 0.106, 0.000], [0.022, -0.005, 0.000], [-0.000, 0.000, 0.105], [0.000, 0.000, -0.011]]),
                    (13, [[4.327, -1.224, 0.000], [0.756, 0.003, 0.000], [-0.000, 0.140, 0.000], [-0.000, 0.052, -0.000], [0.000, -0.000, 0.084], [-0.000, 0.000, -0.018]]),
                    (14, [[4.800, -1.218, 0.000], [0.190, 0.009, 0.000], [-0.010, 0.210, -0.001], [-0.032, 0.087, -0.001], [-0.000, 0.000, 0.070], [-0.000, 0.000, -0.011]]),
                    (15, [[2.450, 1.213, 0.000], [0.979, 0.146, 0.000], [-0.016, 0.106, 0.000], [-0.020, -0.006, 0.000], [0.000, -0.000, 0.105], [0.000, 0.000, -0.011]]),
                    (16, [[4.327, 1.224, 0.000], [0.756, -0.003, 0.000], [0.000, 0.140, 0.000], [0.001, 0.052, 0.000], [-0.000, 0.000, 0.084], [-0.000, -0.000, -0.018]]),
                    (17, [[4.800, 1.218, 0.000], [0.190, -0.009, 0.000], [0.010, 0.210, 0.001], [0.031, 0.087, 0.001], [-0.000, -0.000, 0.070], [-0.000, -0.000, -0.011]]),
                    (18, [[8.100, -0.665, 0.000], [0.775, -0.027, 0.275], [0.005, 0.175, 0.002], [0.022, -0.000, 0.004], [-0.048, -0.000, 0.136], [-0.157, -0.004, -0.039]]),
                    (19, [[8.459, -0.683, 0.298], [0.258, -0.010, 0.370], [0.003, 0.175, 0.003], [-0.005, 0.000, -0.009], [-0.127, 0.001, 0.089], [0.036, -0.004, -0.000]]),
                    (20, [[8.601, -0.685, 0.694], [0.024, 0.006, 0.399], [-0.001, 0.175, -0.002], [-0.002, -0.000, 0.003], [-0.099, -0.000, 0.006], [0.077, 0.003, -0.084]]),
                    (21, [[8.100, 0.665, 0.000], [0.775, 0.027, 0.275], [-0.005, 0.175, -0.002], [-0.022, -0.000, -0.006], [-0.048, 0.000, 0.136], [-0.157, 0.006, -0.039]]),
                    (22, [[8.459, 0.683, 0.298], [0.258, 0.010, 0.370], [-0.000, 0.175, -0.005], [0.005, 0.000, 0.009], [-0.127, 0.002, 0.089], [0.034, 0.004, 0.001]]),
                    (23, [[8.601, 0.685, 0.694], [0.024, -0.006, 0.399], [0.001, 0.175, 0.002], [-0.004, 0.000, -0.003], [-0.104, 0.000, 0.006], [0.071, -0.009, -0.087]])
                ]]),
            "userAnnotationGroups": [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_body_term('head')[0],
                    'ontId': get_body_term('head')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-8',
                    'name': 'torso',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4,11-13, 5,14-16',
                    'name': get_body_term('arm')[0],
                    'ontId': get_body_term('arm')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9,17-18-19, 10,20-21-22',
                    'name': get_body_term('leg')[0],
                    'ontId': get_body_term('leg')[1]
                }
            ]
        })


class MeshType_3d_wholebody2(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network, with linear basis through wall.
    """

    @staticmethod
    def getName():
        return "3D Whole Body 2"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1",
            "Human 2 Coarse",
            "Human 2 Medium",
            "Human 2 Fine"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):

        options = {
            'Base parameter set': parameterSetName,
            "Network layout": getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            "Number of elements around head": 12,
            "Number of elements around torso": 12,
            "Number of elements around arm": 8,
            "Number of elements around leg": 8,
            "Number of elements through wall": 1,
            "Target element density along longest segment": 5.0,
            "Show trim surfaces": False,
            "Use Core": True,
            "Number of elements across head core major": 6,
            "Number of elements across torso core major": 6,
            "Number of elements across arm core major": 4,
            "Number of elements across leg core major": 4,
            "Number of elements across core transition": 1
        }

        if "Human 2 Medium" in parameterSetName:
            options["Number of elements around head"] = 16
            options["Number of elements around torso"] = 16
            options["Number of elements around arm"] = 12
            options["Number of elements around leg"] = 12

            options["Number of elements across head core major"] = 8
            options["Number of elements across torso core major"] = 8
            options["Number of elements across arm core major"] = 6
            options["Number of elements across leg core major"] = 6

        elif "Human 2 Fine" in parameterSetName:
            options["Number of elements around head"] = 24
            options["Number of elements around torso"] = 24
            options["Number of elements around arm"] = 20
            options["Number of elements around leg"] = 20

            options["Number of elements across head core major"] = 10
            options["Number of elements across torso core major"] = 10
            options["Number of elements across arm core major"] = 8
            options["Number of elements across leg core major"] = 8

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Network layout",
            "Number of elements around head",
            "Number of elements around torso",
            "Number of elements around arm",
            "Number of elements around leg",
            "Number of elements through wall",
            "Target element density along longest segment",
            "Show trim surfaces",
            "Use Core",
            "Number of elements across head core major",
            "Number of elements across torso core major",
            "Number of elements across arm core major",
            "Number of elements across leg core major"]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Network layout":
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == "Network layout":
            return cls.getParameterSetNames()
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + ".getOptionScaffoldTypeParameterSetNames.  " + \
            "Invalid option \"" + optionName + "\" scaffold type " + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                "Invalid parameter set " + str(parameterSetName) + " for scaffold " + str(scaffoldType.getName()) + \
                " in option " + str(optionName) + " of scaffold " + cls.getName()
        if optionName == "Network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        for key in [
            "Number of elements around head",
            "Number of elements around torso",
            "Number of elements around arm",
            "Number of elements around leg"
        ]:
            if options[key] < 8:
                options[key] = 8

        for key in [
            "Number of elements across head core major",
            "Number of elements across torso core major",
            "Number of elements across arm core major",
            "Number of elements across leg core major"
        ]:
            if options[key] < 4:
                options[key] = 4

        if options["Number of elements through wall"] < 0:
            options["Number of elements through wall"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0

        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ["Default", "Human 1", "Human 2 Coarse", "Human 2 Medium", "Human 2 Fine"]

        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        annotationElementsCountsAround = []
        annotationElementsCountsAcross = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            elementsCountAround = 0
            elementsCountAcrossMajor = 0
            name = layoutAnnotationGroup.getName()
            if name in ["head", "torso", "arm", "leg"]:
                elementsCountAround = options["Number of elements around " + name]
                elementsCountAcrossMajor = options["Number of elements across " + name + " core major"]
            annotationElementsCountsAround.append(elementsCountAround)
            annotationElementsCountsAcross.append(elementsCountAcrossMajor)

        isCore = options["Use Core"]

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Number of elements around head"],
            elementsCountThroughWall=options["Number of elements through wall"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=annotationElementsCountsAround,
            defaultElementsCountAcrossMajor=options['Number of elements across head core major'],
            elementsCountTransition=options['Number of elements across core transition'],
            annotationElementsCountsAcrossMajor=annotationElementsCountsAcross,
            isCore=isCore)

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None

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

        # create 2d surface mesh groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))

        is_exterior = fm.createFieldIsExterior()
        is_on_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_skin = fm.createFieldAnd(is_exterior, is_on_face_xi3_1)

        skinMeshGroup = skinGroup.getMeshGroup(mesh2d)
        skinMeshGroup.addElementsConditional(is_skin)
