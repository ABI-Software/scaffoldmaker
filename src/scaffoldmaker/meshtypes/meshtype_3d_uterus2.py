"""
Generates a 3-D uterus mesh from a 1-D network layout, with variable
numbers of elements around, along and through wall.
"""
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, group_add_connected_elements


class UterusTubeNetworkMeshGenerateData(TubeNetworkMeshGenerateData):

    def __init__(self, region, meshDimension, isLinearThroughWall, isShowTrimSurfaces,
            coordinateFieldName="coordinates", startNodeIdentifier=1, startElementIdentifier=1):
        """
        :param isLinearThroughWall: Callers should only set if 3-D with no core.
        :param isShowTrimSurfaces: Tells junction generateMesh to make 2-D trim surfaces.
        """
        super(UterusTubeNetworkMeshGenerateData, self).__init__(
            region, meshDimension, isLinearThroughWall, isShowTrimSurfaces,
            coordinateFieldName, startNodeIdentifier, startElementIdentifier)
        self._fundusGroup = self.getOrCreateAnnotationGroup(get_uterus_term("fundus of uterus"))
        self._leftGroup = self.getOrCreateAnnotationGroup(("left uterus", "None"))
        self._rightGroup = self.getOrCreateAnnotationGroup(("right uterus", "None"))
        self._dorsalGroup = self.getOrCreateAnnotationGroup(("dorsal uterus", "None"))
        self._ventralGroup = self.getOrCreateAnnotationGroup(("ventral uterus", "None"))

    def getFundusMeshGroup(self):
        return self._fundusGroup.getMeshGroup(self._mesh)

    def getLeftMeshGroup(self):
        return self._leftGroup.getMeshGroup(self._mesh)

    def getRightMeshGroup(self):
        return self._rightGroup.getMeshGroup(self._mesh)

    def getDorsalMeshGroup(self):
        return self._dorsalGroup.getMeshGroup(self._mesh)

    def getVentralMeshGroup(self):
        return self._ventralGroup.getMeshGroup(self._mesh)

class UterusTubeNetworkMeshBuilder(TubeNetworkMeshBuilder):

    def __init__(self, networkMesh: NetworkMesh, targetElementDensityAlongLongestSegment: float,
                 defaultElementsCountAround: int, elementsCountThroughWall: int,
                 layoutAnnotationGroups: list = [], annotationElementsCountsAround: list = []):
        super(UterusTubeNetworkMeshBuilder, self).__init__(
            networkMesh, targetElementDensityAlongLongestSegment, defaultElementsCountAround,
            elementsCountThroughWall, layoutAnnotationGroups, annotationElementsCountsAround)

    def generateMesh(self, generateData):
        super(UterusTubeNetworkMeshBuilder, self).generateMesh(generateData)
        # build temporary left/right dorsal/ventral groups
        mesh = generateData.getMesh()
        fundusMeshGroup = generateData.getFundusMeshGroup()
        leftMeshGroup = generateData.getLeftMeshGroup()
        rightMeshGroup = generateData.getRightMeshGroup()
        dorsalMeshGroup = generateData.getDorsalMeshGroup()
        ventralMeshGroup = generateData.getVentralMeshGroup()
        for networkSegment in self._networkMesh.getNetworkSegments():
            segment = self._segments[networkSegment]
            elementsCountRim = segment.getElementsCountRim()
            elementsCountAlong = segment.getSampledElementsCountAlong()
            junctions = segment.getJunctions()
            elementsCountAround = segment.getElementsCountAround()
            annotationTerms = segment.getAnnotationTerms()
            preBifurcation = ("pre-bifurcation segments", "None") in annotationTerms
            allLeft = False
            allRight = False
            fundus = preBifurcation and (get_uterus_term("body of uterus") in annotationTerms)
            if preBifurcation:
                for annotationTerm in annotationTerms:
                    if "left" in annotationTerm[0]:
                        allLeft = True
                        break
                    elif "right" in annotationTerm[0]:
                        allRight = True
                        break
                if not (allLeft or allRight):
                    prevSegment = junctions[0].getSegments()[0]
                    if prevSegment is not segment:
                        prevAnnotationTerms = prevSegment.getAnnotationTerms()
                        for annotationTerm in prevAnnotationTerms:
                            if "left" in annotationTerm[0]:
                                allLeft = True
                                break
                            elif "right" in annotationTerm[0]:
                                allRight = True
                                break
            e1FundusLimit = elementsCountAround // 2
            e1RightStart = elementsCountAround // 2
            e1DVStart = e1RightStart // 2
            e1DVEnd = elementsCountAround - e1DVStart - 1
            for e1 in range(elementsCountAround):
                meshGroups = []
                if fundus and ((allLeft and (e1 < e1FundusLimit)) or (allRight and (e1 >= e1FundusLimit))):
                    meshGroups.append(fundusMeshGroup)
                if allLeft or ((not allRight) and (e1 < e1RightStart)):
                    meshGroups.append(leftMeshGroup)
                else:
                    meshGroups.append(rightMeshGroup)
                if ((preBifurcation and (e1 >= e1DVStart) and (e1 <= e1DVEnd)) or
                        ((not preBifurcation) and ((e1 < e1DVStart) or (e1 > e1DVEnd)))):
                    meshGroups.append(dorsalMeshGroup)
                else:
                    meshGroups.append(ventralMeshGroup)
                for e2 in range(elementsCountAlong):
                    for e3 in range(elementsCountRim):
                        elementIdentifier = segment.getRimElementId(e1, e2, e3)
                        if elementIdentifier is not None:
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            for meshGroup in meshGroups:
                                meshGroup.addElement(element)


def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames() # make sure parameter set is in list of parameters of parent scaffold
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3, 4-5-6, 3-7-8-11.1, 6-9-10-11.2, 11.3-12-13-14,14-15-16,16-17-18",
                "Define inner coordinates": True,
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[-9.160,10.270,-3.090], [ 0.889, 0.730,-1.149], [ 0.042,-0.228,-0.113], [ 0.083, 0.023, 0.006], [-0.202, 0.031,-0.137], [ 0.077,-0.062,-0.091]]),
                (2, [[-8.060,10.640,-4.320], [ 1.271,-0.020,-1.261], [-0.002,-0.250, 0.002], [-0.103, 0.007, 0.139], [-0.169, 0.000,-0.171], [-0.024,-0.002, 0.016]]),
                (3, [[-6.710,10.190,-5.510], [ 1.730,-0.520,-0.850], [-0.060,-0.250, 0.030], [-0.106,-0.253,-0.053], [-0.118, 0.000,-0.241], [ 0.118, 0.006,-0.372]]),
                (4, [[ 9.160,10.270,-3.090], [-0.889, 0.730,-1.149], [-0.042,-0.228,-0.113], [-0.083, 0.023, 0.006], [-0.202,-0.031, 0.137], [ 0.077, 0.062, 0.091]]),
                (5, [[ 8.060,10.640,-4.320], [-1.271,-0.020,-1.261], [ 0.002,-0.250, 0.002], [ 0.103, 0.007, 0.139], [-0.169, 0.000, 0.171], [-0.024, 0.002,-0.016]]),
                (6, [[ 6.710,10.190,-5.510], [-1.730,-0.520,-0.850], [ 0.060,-0.250, 0.030], [ 0.106,-0.253,-0.053], [-0.118, 0.000, 0.241], [ 0.118,-0.006, 0.372]]),
                (7, [[-4.710, 9.650,-5.890], [ 2.182,-0.210,-0.350], [-0.090,-0.770,-0.100], [ 0.319,-1.159,-0.237], [-0.114, 0.114,-0.779], [-0.339, 0.211,-0.805]]),
                (8, [[-2.420, 9.800,-6.200], [ 2.364, 0.210,-0.170], [ 0.188,-2.460,-0.420], [-0.027,-1.024,-0.146], [-0.169, 0.320,-1.950], [ 0.243, 0.218,-1.088]]),
                (9, [[ 4.710, 9.650,-5.890], [-2.182,-0.210,-0.350], [ 0.090,-0.770,-0.100], [-0.319,-1.159,-0.237], [-0.114,-0.114, 0.779], [-0.339,-0.211, 0.805]]),
                (10, [[ 2.420, 9.800,-6.200], [-2.364, 0.210,-0.170], [-0.188,-2.460,-0.420], [ 0.027,-1.024,-0.146], [-0.169,-0.320, 1.950], [ 0.243,-0.218, 1.088]]),
                (11, [[ 0.000,10.080,-6.230],
                      [[ 2.471, 0.340, 0.110],[-2.471, 0.340, 0.110],[ 0.000,-0.350, 2.812]],
                      [[ 0.400,-2.761,-0.461],[-0.400,-2.761,-0.461],[ 0.000, 2.929, 0.365]],
                      [[ 0.196, 0.407, 0.058],[-0.196, 0.407, 0.058],[ 0.000,-0.216,-0.046]],
                      [[0.062, 0.495, -2.911], [0.062, -0.495, 2.911], [-3.520, 0.000, 0.000]],
                      [[-0.011, 0.101,-0.853],[-0.011,-0.101, 0.853],[ 0.470, 0.000, 0.000]]]),
                (12, [[ 0.000, 9.500,-3.620], [ 0.000,-0.810, 2.390], [ 0.000, 2.369, 0.803], [ 0.000,-0.819, 0.290], [-3.080, 0.000, 0.000], [ 0.410, 0.000, 0.000]]),
                (13, [[ 0.000, 8.500,-1.480], [ 0.000,-1.162, 1.773], [ 0.000, 1.329, 0.871], [ 0.000,-0.823, 0.095], [-2.700, 0.000, 0.000], [ 0.285, 0.000, 0.000]]),
                (14, [[ 0.000, 7.270,-0.080], [ 0.000,-1.091, 0.821], [ 0.000, 0.747, 0.993], [ 0.000,-0.491, 0.093], [-2.510, 0.000, 0.000], [ 0.180, 0.000, 0.000]]),
                (15, [[ 0.000, 6.500, 0.280], [ 0.000,-0.828, 0.179], [ 0.000, 0.216, 1.001], [ 0.000,-0.302,-0.201], [-2.340, 0.000, 0.000], [ 0.220, 0.000, 0.000]]),
                (16, [[ 0.000, 5.670, 0.280], [ 0.000,-1.578,-0.050], [ 0.000,-0.026, 0.810], [ 0.000, 0.109,-0.265], [-2.070, 0.000, 0.000], [ 0.630, 0.000, 0.000]]),
                (17, [[ 0.000, 3.350, 0.140], [ 0.000,-2.850,-0.140], [ 0.000,-0.023, 0.460], [ 0.000, 0.021,-0.129], [-1.080, 0.000, 0.000], [ 0.760, 0.000, 0.000]]),
                (18, [[ 0.000,-0.030, 0.000], [ 0.000,-3.910,-0.130], [ 0.000,-0.018, 0.550], [ 0.000,-0.014, 0.309], [-0.550, 0.000, 0.000], [ 0.300, 0.000, 0.000]])], [
                (1, [[-9.160,10.270,-3.090], [ 0.889, 0.730,-1.149], [ 0.018,-0.112,-0.057], [ 0.043, 0.018, 0.011], [-0.102, 0.018,-0.067], [ 0.040,-0.036,-0.050]]),
                (2, [[-8.060,10.640,-4.320], [ 1.271,-0.020,-1.261], [-0.001,-0.120, 0.001], [-0.050, 0.000, 0.067], [-0.085, 0.000,-0.085], [-0.012, 0.000, 0.008]]),
                (3, [[-6.710,10.190,-5.510], [ 1.730,-0.520,-0.850], [-0.033,-0.129, 0.012], [-0.050,-0.127,-0.029], [-0.058, 0.004,-0.121], [ 0.059, 0.004,-0.186]]),
                (4, [[ 9.160,10.270,-3.090], [-0.889, 0.730,-1.149], [-0.018,-0.112,-0.057], [-0.043, 0.018, 0.011], [-0.102,-0.018, 0.067], [ 0.040, 0.036, 0.050]]),
                (5, [[ 8.060,10.640,-4.320], [-1.271,-0.020,-1.261], [ 0.001,-0.120, 0.001], [ 0.050, 0.000, 0.067], [-0.085, 0.000, 0.085], [-0.012, 0.000,-0.008]]),
                (6, [[ 6.710,10.190,-5.510], [-1.730,-0.520,-0.850], [ 0.033,-0.129, 0.012], [ 0.050,-0.127,-0.029], [-0.058,-0.004, 0.121], [ 0.059,-0.004, 0.186]]),
                (7, [[-4.710, 9.650,-5.890], [ 2.182,-0.210,-0.350], [-0.045,-0.381,-0.051], [ 0.157,-0.594,-0.107], [-0.057, 0.059,-0.391], [-0.172, 0.094,-0.403]]),
                (8, [[-2.420, 9.800,-6.200], [ 2.364, 0.210,-0.170], [ 0.098,-1.262,-0.191], [-0.014,-0.567,-0.074], [-0.083, 0.141,-0.976], [ 0.121, 0.102,-0.547]]),
                (9, [[ 4.710, 9.650,-5.890], [-2.182,-0.210,-0.350], [ 0.045,-0.381,-0.051], [-0.157,-0.594,-0.107], [-0.057,-0.059, 0.391], [-0.172,-0.094, 0.403]]),
                (10, [[ 2.420, 9.800,-6.200], [-2.364, 0.210,-0.170], [-0.098,-1.262,-0.191], [ 0.014,-0.567,-0.074], [-0.083,-0.141, 0.976], [ 0.121,-0.102, 0.547]]),
                (11, [[ 0.000,10.080,-6.230],
                      [[ 2.471, 0.340, 0.110],[-2.471, 0.340, 0.110],[ 0.000,-0.350, 2.812]],
                      [[ 0.215,-1.485,-0.234],[-0.215,-1.485,-0.234],[ 0.000, 1.460, 0.182]],
                      [[ 0.122, 0.114,-0.016],[-0.122, 0.114,-0.016],[ 0.000,-0.114,-0.024]],
                      [[ 0.033, 0.235,-1.463],[ 0.033,-0.235, 1.463],[-1.760, 0.000, 0.000]],
                      [[-0.005, 0.071,-0.436],[-0.005,-0.071, 0.436],[ 0.235, 0.000, 0.000]]]),
                (12, [[ 0.000, 9.500,-3.620], [ 0.000,-0.810, 2.390], [ 0.000, 1.180, 0.400], [ 0.000,-0.404, 0.146], [-1.540, 0.000, 0.000], [ 0.205, 0.000, 0.000]]),
                (13, [[ 0.000, 8.500,-1.480], [ 0.000,-1.162, 1.773], [ 0.000, 0.670, 0.439], [ 0.000,-0.393, 0.063], [-1.350, 0.000, 0.000], [ 0.590, 0.000, 0.000]]),
                (14, [[ 0.000, 7.270,-0.080], [ 0.000,-1.091, 0.821], [ 0.000, 0.399, 0.531], [ 0.000,-0.261, 0.055], [-0.360, 0.000, 0.000], [ 0.470, 0.000, 0.000]]),
                (15, [[ 0.000, 6.500, 0.280], [ 0.000,-0.828, 0.179], [ 0.000, 0.108, 0.500], [ 0.000,-0.151,-0.120], [-0.410, 0.000, 0.000], [ 0.000, 0.000, 0.000]]),
                (16, [[ 0.000, 5.670, 0.280], [ 0.000,-1.578,-0.050], [ 0.000,-0.013, 0.410], [ 0.000, 0.055,-0.132], [-0.360, 0.000, 0.000], [-0.065, 0.000, 0.000]]),
                (17, [[ 0.000, 3.350, 0.140], [ 0.000,-2.850,-0.140], [ 0.000,-0.011, 0.230], [ 0.000, 0.010,-0.064], [-0.540, 0.000, 0.000], [ 0.040, 0.000, 0.000]]),
                (18, [[ 0.000,-0.030, 0.000], [ 0.000,-3.910,-0.130], [ 0.000,-0.009, 0.280], [ 0.000,-0.007, 0.165], [-0.280, 0.000, 0.000], [ 0.480, 0.000, 0.000]])]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-17',
                    'name': get_uterus_term('uterus')[0],
                    'ontId': get_uterus_term('uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_uterus_term('left uterine tube')[0],
                    'ontId': get_uterus_term('left uterine tube')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_uterus_term('right uterine tube')[0],
                    'ontId': get_uterus_term('right uterine tube')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-13',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14-15',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '16-17',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': 'pre-bifurcation segments',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-17',
                    'name': 'post-bifurcation segments',
                    'ontId': 'None'
                }
            ]
        })
    elif "Mouse 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-10, 5-11.1, 10-11.2, 11.3-12-13-14, 14-15, 15-16",
                "Define inner coordinates": True,
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[                    
                (1, [[-24.20,54.56,0.00], [2.53,-14.75,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.42,0.41,0.00], [0.14,-0.02,0.00]]),
                (2, [[-20.88,41.34,0.00], [4.10,-11.62,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.47,0.87,0.00], [0.05,0.43,0.00]]),
                (3, [[-16.28,31.37,0.00], [5.76,-8.94,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.38,1.54,0.00], [-0.55,0.82,0.00]]),
                (4, [[-9.59,23.63,0.00], [6.74,-6.07,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.74,1.94,0.00], [-0.66,0.45,0.00]]),
                (5, [[-3.14,19.13,0.00], [4.91,-2.86,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.31,2.26,0.00], [-0.29,0.21,0.00]]),
                (6, [[24.20,54.56,0.00], [-2.53,-14.75,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.42,-0.41,0.00], [0.14,0.02,0.00]]),
                (7, [[20.88,41.34,0.00], [-4.10,-11.62,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.47,-0.87,0.00], [0.05,-0.43,0.00]]),
                (8, [[16.28,31.37,0.00], [-5.76,-8.94,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.38,-1.54,0.00], [-0.55,-0.82,0.00]]),
                (9, [[9.59,23.63,0.00], [-6.74,-6.07,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.74,-1.94,0.00], [-0.66,-0.45,0.00]]),
                (10, [[3.14,19.13,0.00], [-4.91,-2.86,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.31,-2.26,0.00], [-0.29,-0.21,0.00]]),
                (11, [[0.00,17.63,-0.00], [[1.29,-0.13,-0.00],[-1.29,-0.13,-0.00],[0.00,-3.86,-0.00]], [[0.00,0.00,-2.50],[0.00,0.00,-2.50],[0.00,0.00,2.50]], [[-0.00,-0.00,-0.00],[-0.00,0.00,-0.00],[0.00,0.00,0.00]], [[0.28,2.81,0.00],[0.28,-2.81,0.00],[-2.67,-0.00,0.00]], [[-3.19,1.05,-0.00],[-3.19,-1.05,-0.00],[0.25,0.00,0.00]]]),
                (12, [[0.00,13.95,-0.00], [0.00,-3.50,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.09,0.00,0.00]]),
                (13, [[0.00,10.63,-0.00], [0.00,-3.47,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (14, [[0.00,7.02,-0.00], [0.00,-3.31,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (15, [[0.00,4.00,-0.00], [0.00,-3.51,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (16, [[0.00,0.00,-0.00], [0.00,-4.49,-0.00], [0.00,0.00,2.50], [-0.00,-0.00,0.00], [-2.50,-0.00,0.00], [-0.00,-0.00,-0.00]])],[
                (1, [[-24.20,54.56,0.00], [2.53,-14.75,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.45,0.25,0.00], [0.08,-0.01,0.00]]),
                (2, [[-20.88,41.34,0.00], [4.10,-11.62,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.48,0.52,0.00], [0.03,0.26,0.00]]),
                (3, [[-16.28,31.37,0.00], [5.76,-8.94,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.43,0.92,0.00], [-0.33,0.50,0.00]]),
                (4, [[-9.59,23.63,0.00], [6.74,-6.07,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.05,1.16,0.00], [-0.40,0.27,0.00]]),
                (5, [[-3.14,19.13,0.00], [4.91,-2.86,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [0.79,1.35,0.00], [-0.18,0.12,0.00]]),
                (6, [[24.20,54.56,0.00], [-2.53,-14.75,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.45,-0.25,0.00], [0.08,0.01,0.00]]),
                (7, [[20.88,41.34,0.00], [-4.10,-11.62,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.48,-0.52,0.00], [0.03,-0.26,0.00]]),
                (8, [[16.28,31.37,0.00], [-5.76,-8.94,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.43,-0.92,0.00], [-0.33,-0.50,0.00]]),
                (9, [[9.59,23.63,0.00], [-6.74,-6.07,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.05,-1.16,0.00], [-0.40,-0.27,0.00]]),
                (10, [[3.14,19.13,0.00], [-4.91,-2.86,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [0.79,-1.35,0.00], [-0.18,-0.12,0.00]]),
                (11, [[0.00,17.63,-0.00], [[1.29,-0.13,-0.00],[-1.29,-0.13,-0.00],[0.00,-3.86,-0.00]], [[0.00,0.00,-1.50],[0.00,0.00,-1.50],[0.00,0.00,1.50]], [[-0.00,-0.00,-0.00],[-0.00,0.00,-0.00],[0.00,0.00,0.00]], [[0.17,1.68,0.00],[0.17,-1.68,0.00],[-1.60,-0.00,0.00]], [[-1.91,0.62,-0.00],[-1.91,-0.62,-0.00],[0.15,0.00,0.00]]]),
                (12, [[0.00,13.95,-0.00], [0.00,-3.50,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [0.05,0.00,0.00]]),
                (13, [[0.00,10.63,-0.00], [0.00,-3.47,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [-0.00,0.00,0.00]]),
                (14, [[0.00,7.02,-0.00], [0.00,-3.31,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [-0.00,0.00,0.00]]),
                (15, [[0.00,4.00,-0.00], [0.00,-3.51,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (16, [[0.00,0.00,-0.00], [0.00,-4.49,-0.00], [0.00,0.00,1.50], [-0.00,-0.00,-0.00], [-1.50,-0.00,0.00], [-0.00,-0.00,-0.00]])]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-15',
                    'name': get_uterus_term('uterus')[0],
                    'ontId': get_uterus_term('uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-13',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '15',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': 'pre-bifurcation segments',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-15',
                    'name': 'post-bifurcation segments',
                    'ontId': 'None'
                },
            ]
        })


class MeshType_3d_uterus2(Scaffold_base):
    """
    Generates a 3-D uterus mesh from a 1-D network layout with variable numbers of elements around, along and through
    wall.
    Magnitude of D2 and D3 are the radii of the uterus in the respective directions.
    """

    @classmethod
    def getName(cls):
        return '3D Uterus 2'

    @classmethod
    def getParameterSetNames(cls):
        return [
            'Default',
            'Human 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):

        options = {
            'Base parameter set': parameterSetName,
            'Network layout': getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            'Number of elements around': 12,
            'Number of elements around horns': 12,
            'Number of elements through wall': 1,
            'Target element density along longest segment': 5.0,
            'Use linear through wall': True,
            'Show trim surfaces': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Mouse' in parameterSetName:
            options['Number of elements around'] = 8
            options['Number of elements around horns'] = 8
            options['Target element density along longest segment'] = 10.0

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            'Network layout',
            'Number of elements around',
            'Number of elements around horns',
            'Number of elements through wall',
            'Target element density along longest segment',
            'Use linear through wall',
            'Show trim surfaces',
            'Refine',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return cls.getParameterSetNames()
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Network layout':
            if not parameterSetName:
                parameterSetName = "Default"
            return getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Network layout'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Network layout'):
            options['Network layout'] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        for key in [
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Number of elements around',
            'Number of elements around horns']:
            if options[key] < 4:
                options[key] = 4
            elif (options[key] % 4) > 0:
                options[key] += options[key] % 4
        if options["Number of elements through wall"] < 1:
            options["Number of elements through wall"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic hermite or bicubic hermite-linear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ("Default", "Human 1")

        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        annotationElementsCountsAround = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            elementsCountAround = 0
            if layoutAnnotationGroup.getName() == "pre-bifurcation segments":
                elementsCountAround = options['Number of elements around horns']
            elif layoutAnnotationGroup.getName() == "post-bifurcation segments":
                elementsCountAround = options['Number of elements around']
            annotationElementsCountsAround.append(elementsCountAround)

        uterusTubeNetworkMeshBuilder = UterusTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options['Number of elements around'],
            elementsCountThroughWall=options["Number of elements through wall"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=annotationElementsCountsAround)
        uterusTubeNetworkMeshBuilder.build()
        generateData = UterusTubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=options["Use linear through wall"],
            isShowTrimSurfaces=options["Show trim surfaces"])
        uterusTubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()
        nodeIdentifier, _ = generateData.getNodeElementIdentifiers()

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldmodule.findMeshByDimension(3)

        uterusGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("uterus"))
        myometriumGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_uterus_term("myometrium"))
        myometriumGroup.getMeshGroup(mesh).addElementsConditional(uterusGroup.getGroup())

        if isHuman:
            # add human specific annotations

            allMarkers = {"junction of left round ligament with uterus": {"x": [-4.1312, 9.96436, -7.11994]},
                          "junction of right round ligament with uterus": {"x": [4.13116, 9.96438, -7.11997]}}

            for key in allMarkers:
                x = allMarkers[key]["x"]
                group = findOrCreateAnnotationGroupForTerm( annotationGroups, region,
                                                            get_uterus_term(key), isMarker=True)
                markerNode = group.createMarkerNode(nodeIdentifier, coordinates, x)
                nodeIdentifier = markerNode.getIdentifier() + 1
                for group in annotationGroups:
                    group.getNodesetGroup(nodes).addNode(markerNode)

        # remove temporary annotation groups
        for annotationTerms in [("pre-bifurcation segments", "None"), ("post-bifurcation segments", "None")]:
            annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, annotationTerms)
            annotationGroups.remove(annotationGroup)

        return annotationGroups, None

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

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
        isHuman = parameterSetName in ("Default", "Human 1")
        isMouse = parameterSetName in "Mouse 1"

        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        uterusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterus"))
        cervixGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterine cervix"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("body of uterus"))
        vaginaGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("vagina"))

        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
        is_exterior_face_xi2_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))

        is_uterus = uterusGroup.getGroup()
        is_uterus_outer = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_1)
        is_uterus_inner = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_0)

        is_body = bodyGroup.getGroup()
        is_body_outer = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_inner = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_cervix = cervixGroup.getGroup()
        is_cervix_outer = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_1)
        is_cervix_inner = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_0)

        is_vagina = vaginaGroup.getGroup()
        is_vagina_outer = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_1)
        is_vagina_inner = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_0)
        is_vagina_xi2_0 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_0)
        is_vagina_xi2_1 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_1)
        is_vagina_xi2_01 = fm.createFieldXor(is_vagina_xi2_0, is_vagina_xi2_1)

        serosaOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterus"))
        serosaOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_outer)

        lumenOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterus"))
        lumenOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_inner)

        serosaOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                          get_uterus_term("serosa of body of uterus"))
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_outer)

        lumenOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_uterus_term("lumen of body of uterus"))
        lumenOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_inner)

        serosaOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterine cervix"))
        serosaOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_outer)

        lumenOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterine cervix"))
        lumenOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_inner)

        serosaOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of vagina"))
        serosaOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_outer)

        lumenOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of vagina"))
        lumenOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_inner)

        leftGroup = getAnnotationGroupForTerm(annotationGroups, ("left uterus", "None"))
        rightGroup = getAnnotationGroupForTerm(annotationGroups, ("right uterus", "None"))
        dorsalGroup = getAnnotationGroupForTerm(annotationGroups, ("dorsal uterus", "None"))
        ventralGroup = getAnnotationGroupForTerm(annotationGroups, ("ventral uterus", "None"))

        if isHuman:
            leftUterineTubeGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left uterine tube"))
            is_leftUterineTube = leftUterineTubeGroup.getGroup()
            is_leftUterineTube_outer = fm.createFieldAnd(is_leftUterineTube, is_exterior_face_xi3_1)
            is_leftUterineTube_inner = fm.createFieldAnd(is_leftUterineTube, is_exterior_face_xi3_0)

            serosaOfLeftUterineTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of left uterine tube"))
            serosaOfLeftUterineTube.getMeshGroup(mesh2d).addElementsConditional(is_leftUterineTube_outer)

            rightUterineTubeGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right uterine tube"))
            is_rightUterineTube = rightUterineTubeGroup.getGroup()
            is_rightUterineTube_outer = fm.createFieldAnd(is_rightUterineTube, is_exterior_face_xi3_1)
            is_rightUterineTube_inner = fm.createFieldAnd(is_rightUterineTube, is_exterior_face_xi3_0)

            serosaOfRightUterineTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                         get_uterus_term("serosa of right uterine tube"))
            serosaOfRightUterineTube.getMeshGroup(mesh2d).addElementsConditional(is_rightUterineTube_outer)

            lumenOfLeftUterineTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of left uterine tube"))
            lumenOfLeftUterineTube.getMeshGroup(mesh2d).addElementsConditional(is_leftUterineTube_inner)

            lumenOfRightUterineTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                        get_uterus_term("lumen of right uterine tube"))
            lumenOfRightUterineTube.getMeshGroup(mesh2d).addElementsConditional(is_rightUterineTube_inner)

            is_pubocervical = fm.createFieldAnd(is_body_outer, is_cervix_outer)
            pubocervical = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("pubocervical ligament (TA98)"))
            pubocervical.getMeshGroup(mesh1d).addElementsConditional(is_pubocervical)

            is_internal_os = fm.createFieldAnd(is_body_inner, is_cervix_inner)
            internalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_uterus_term("internal cervical os"))
            internalOs.getMeshGroup(mesh1d).addElementsConditional(is_internal_os)

            is_external_os = fm.createFieldAnd(is_vagina_inner, is_cervix_inner)
            externalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("external cervical os"))
            externalOs.getMeshGroup(mesh1d).addElementsConditional(is_external_os)

            is_vagina_orifice = fm.createFieldAnd(is_vagina_xi2_01, is_exterior_face_xi3_0)
            vaginaOrifice = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("vagina orifice"))
            vaginaOrifice.getMeshGroup(mesh1d).addElementsConditional(is_vagina_orifice)

            # ligaments
            is_dorsalVentral = fm.createFieldAnd(dorsalGroup.getGroup(), ventralGroup.getGroup())
            is_dorsalVentralSerosa = fm.createFieldAnd(is_dorsalVentral, is_exterior_face_xi3_1)
            is_leftDorsalVentralSerosa = fm.createFieldAnd(leftGroup.getGroup(), is_dorsalVentralSerosa)
            is_rightDorsalVentralSerosa = fm.createFieldAnd(rightGroup.getGroup(), is_dorsalVentralSerosa)
            fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("fundus of uterus"))
            is_bodyNotFundus = fm.createFieldAnd(bodyGroup.getGroup(), fm.createFieldNot(fundusGroup.getGroup()))

            # Broad ligament of uterus
            is_leftBroadLigament = fm.createFieldAnd(is_bodyNotFundus, is_leftDorsalVentralSerosa)
            leftBroadLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("left broad ligament of uterus"))
            leftBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_leftBroadLigament)
            # add connected edges from left uterine tube, avoiding adding dorsal-ventral edges on the superior edge
            leftBroadLigament.addSubelements()  # need current nodes in ligament for group_add_connected_elements
            tmpGroup = fm.createFieldGroup()
            tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
            tmpMeshGroup.addElementsConditional(fm.createFieldAnd(is_leftUterineTube, is_leftDorsalVentralSerosa))
            group_add_connected_elements(leftBroadLigament.getGroup(), tmpMeshGroup)
            del tmpMeshGroup
            del tmpGroup

            is_rightBroadLigament = fm.createFieldAnd(is_bodyNotFundus, is_rightDorsalVentralSerosa)
            rightBroadLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("right broad ligament of uterus"))
            rightBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_rightBroadLigament)
            # add connected edges from right uterine tube, avoiding adding dorsal-ventral edges on the superior edge
            rightBroadLigament.addSubelements()  # need current nodes in ligament for group_add_connected_elements
            tmpGroup = fm.createFieldGroup()
            tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
            tmpMeshGroup.addElementsConditional(fm.createFieldAnd(is_rightUterineTube, is_rightDorsalVentralSerosa))
            group_add_connected_elements(rightBroadLigament.getGroup(), tmpMeshGroup)
            del tmpMeshGroup
            del tmpGroup

            # Transverse cervical ligament
            is_leftTransverseCervicalLigament = fm.createFieldAnd(cervixGroup.getGroup(), is_leftDorsalVentralSerosa)
            leftTransverseCervicalLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("left transverse cervical ligament"))
            leftTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(
                is_leftTransverseCervicalLigament)

            is_rightTransverseCervicalLigament = fm.createFieldAnd(cervixGroup.getGroup(), is_rightDorsalVentralSerosa)
            rightTransverseCervicalLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("right transverse cervical ligament"))
            rightTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(
                is_rightTransverseCervicalLigament)

        if isMouse:
            rightHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right uterine horn"))
            leftHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left uterine horn"))

            is_rightHorn = rightHornGroup.getGroup()
            is_rightHorn_outer = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_1)
            is_rightHorn_inner = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_0)

            is_leftHorn = leftHornGroup.getGroup()
            is_leftHorn_outer = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_1)
            is_leftHorn_inner = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_0)

            serosaOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of right horn"))
            serosaOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_outer)

            lumenOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of right horn"))
            lumenOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_inner)

            serosaOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("serosa of left horn"))
            serosaOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_outer)

            lumenOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_uterus_term("lumen of left horn"))
            lumenOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_inner)

        # keeping these in for now
        # annotationGroups.remove(leftGroup)
        # annotationGroups.remove(rightGroup)
        # annotationGroups.remove(dorsalGroup)
        # annotationGroups.remove(ventralGroup)
