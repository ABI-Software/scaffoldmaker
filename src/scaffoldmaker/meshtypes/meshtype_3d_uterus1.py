"""
Generates a 3-D uterus mesh along the central line, with variable
numbers of elements around, along and through wall.
"""

import copy
import math

from cmlibs.maths.vectorops import add, mult
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation import get_tube_bifurcation_connection_elements_counts, get_bifurcation_triple_point
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, get_nodeset_path_field_parameters

class MeshType_3d_uterus1(Scaffold_base):
    """
    Generates a 3-D uterus mesh with variable numbers of elements around, along the central line, and through wall.
    The uterus is created using a central path as the longitudinal axis of the uterus. Magnitude of D2 and D3 are
    the radii of the uterus in the respective directions.
    """
    parameterSetStructureStrings = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-5.2, 5.3-10-11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[13.84, 2.07, 20.50], [-0.53, -0.68, -6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24], [0.08, -0.99, 0.10], [-0.02, -0.01, -0.01]]),
                    (2, [[12.56, 1.40, 14.10], [-2.03, -0.65, -6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22], [0.06, -1.00, 0.09], [-0.01, -0.00, -0.01]]),
                    (3, [[9.85, 0.78, 8.39], [-3.33, -0.59, -5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22], [0.06, -1.00, 0.07], [-0.00, -0.00, -0.03]]),
                    (4, [[5.97, 0.24, 3.53], [-4.99, -0.41, -4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20], [0.05, -1.00, 0.03], [0.01, 0.01, -0.11]]),
                    (5, [[0.00, 0.00, 0.00], [[-6.82, -0.07, -2.67], [6.82, -0.07, -2.67], [0.00, 0.00, -2.00]], [0.36, 0.18, -0.92], [-0.10, 0.09, 0.10], [0.07, -0.98, -0.16], [-0.01, 0.01, -0.10]]),
                    (6, [[-13.84, 2.07, 20.50], [0.53, -0.68, -6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40], [-0.05, -0.99, 0.10], [-0.03, -0.00, 0.04]]),
                    (7, [[-12.56, 1.40, 14.10], [2.03, -0.65, -6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22], [-0.06, -0.99, 0.09], [-0.01, -0.01, -0.01]]),
                    (8, [[-9.85, 0.78, 8.39], [3.33, -0.59, -5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22], [-0.07, -1.00, 0.07], [0.00, -0.01, -0.03]]),
                    (9, [[-5.97, 0.24, 3.53], [4.99, -0.41, -4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20], [-0.05, -1.00, 0.03], [0.03, 0.00, -0.03]]),
                    (10, [[0.00, 0.00, -2.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, -1.00, 0.00], [0.01, -0.00, -0.01]]),
                    (11, [[0.00, 0.00, -4.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (12, [[0.00, 0.00, -6.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (13, [[0.00, 0.00, -8.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (14, [[0.00, 0.00, -10.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]])]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-10',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-13',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        })
    }

    @staticmethod
    def getName():
        return '3D Uterus 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        networkLayoutOption = cls.parameterSetStructureStrings['Mouse 1']
        options = {
            'Network layout': copy.deepcopy(networkLayoutOption),
            'Target element length': 4.0,
            'Number of elements around': 8,
            'Wall thickness': 1.0,
            'Number of elements through wall': 1,
            'Use linear through wall': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Target element length',
            'Number of elements around',
            'Wall thickness',
            'Number of elements through wall',
            'Use linear through wall',
            'Use cross derivatives',
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
            return list(cls.parameterSetStructureStrings.keys())
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
                parameterSetName = list(cls.parameterSetStructureStrings.keys())[1]
            return copy.deepcopy(cls.parameterSetStructureStrings[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Network layout'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Network layout'):
            options['Network layout'] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        for key in [
            'Number of elements around',
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around'] % 2 != 0:
            options['Number of elements around'] += 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        networkLayout = options['Network layout']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        targetElementLength = options['Target element length']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        mesh = fm.findMeshByDimension(3)
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        firstNodeIdentifier = 1
        elementIdentifier = 1

        # Create annotation groups
        rightHornGroup = AnnotationGroup(region, get_uterus_term("right uterine horn"))
        leftHornGroup = AnnotationGroup(region, get_uterus_term("left uterine horn"))
        cervixGroup = AnnotationGroup(region, get_uterus_term("uterine cervix"))
        vaginaGroup = AnnotationGroup(region, get_uterus_term("vagina"))
        uterusGroup = AnnotationGroup(region, get_uterus_term("uterus"))
        annotationGroups = [cervixGroup, vaginaGroup, leftHornGroup, rightHornGroup, uterusGroup]

        rightHornMeshGroup = rightHornGroup.getMeshGroup(mesh)
        leftHornMeshGroup = leftHornGroup.getMeshGroup(mesh)
        cervixMeshGroup = cervixGroup.getMeshGroup(mesh)
        vaginaMeshGroup = vaginaGroup.getMeshGroup(mesh)
        uterusMeshGroup = uterusGroup.getMeshGroup(mesh)

        # Geometric coordinates
        geometricCentralPath = UterusCentralPath(region, networkLayout, targetElementLength)

        rightHornLength = geometricCentralPath.rightHornLength
        leftHornLength = geometricCentralPath.leftHornLength
        cervixLength = geometricCentralPath.cervixLength
        vaginaLength = geometricCentralPath.vaginaLength

        elementsCountInRightHorn = math.ceil(rightHornLength / targetElementLength)
        elementsCountInLeftHorn = math.ceil(leftHornLength / targetElementLength)
        elementsCountInCervix = math.ceil(cervixLength / targetElementLength)
        elementsCountInVagina = math.ceil(vaginaLength / targetElementLength)

        cx_right_horn_group = geometricCentralPath.cx_right_horn_group
        cx_left_horn_group = geometricCentralPath.cx_left_horn_group
        cx_cervix_group = geometricCentralPath.cx_cervix_group
        cx_vagina_group = geometricCentralPath.cx_vagina_group

        sx_right_horn_group = geometricCentralPath.sx_right_horn_group
        sx_left_horn_group = geometricCentralPath.sx_left_horn_group
        sx_cervix_group = geometricCentralPath.sx_cervix_group
        sx_vagina_group = geometricCentralPath.sx_vagina_group

        # Get right horn nodes
        rightHornCoordinates = getCoordinatesAlongTube3D(cx_right_horn_group, elementsCountAround,
                                                         elementsCountInRightHorn, elementsCountThroughWall,
                                                         wallThickness, startRadian=-math.pi/2)

        rhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(rightHornCoordinates, elementsCountAround,
                                                                    elementsCountInRightHorn, elementsCountThroughWall,
                                                                    omitStartRows=0, omitEndRows=1)

        rhLastRingNodeId, nodeCount = getTargetedRingNodesId(firstNodeIdentifier, elementsCountAround,
                                                             elementsCountInRightHorn, elementsCountThroughWall,
                                                             omitStartRows=0, omitEndRows=1)

        # Get left horn nodes
        leftHornCoordinates = getCoordinatesAlongTube3D(cx_left_horn_group, elementsCountAround, elementsCountInLeftHorn,
                                                        elementsCountThroughWall, wallThickness, startRadian=-math.pi/2)

        lhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(leftHornCoordinates, elementsCountAround,
                                                                    elementsCountInLeftHorn, elementsCountThroughWall,
                                                                    omitStartRows=0, omitEndRows=1)

        lhLastRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountInLeftHorn,
                                                             elementsCountThroughWall, omitStartRows=0, omitEndRows=1)

        # Get cervix nodes
        cervixCoordinates = getCoordinatesAlongTube3D(cx_cervix_group, elementsCountAround, elementsCountInCervix,
                                                      elementsCountThroughWall, wallThickness, startRadian=-math.pi/2)

        cFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(cervixCoordinates, elementsCountAround,
                                                                    elementsCountInCervix, elementsCountThroughWall,
                                                                    omitStartRows=1, omitEndRows=0)

        # Get vagina nodes
        vaginaCoordinates = getCoordinatesAlongTube3D(cx_vagina_group, elementsCountAround, elementsCountInVagina,
                                                      elementsCountThroughWall, wallThickness, startRadian=-math.pi/2)

        vFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(vaginaCoordinates, elementsCountAround,
                                                                    elementsCountInVagina, elementsCountThroughWall,
                                                                    omitStartRows=0, omitEndRows=0)

        # Create nodes
        # Create right horn nodes
        nodeIdentifier = generateTubeNodes(fm, firstNodeIdentifier, rightHornCoordinates, elementsCountInRightHorn,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=0,
                                           omitEndRows=1, startNodes=None)

        # Create left horn nodes
        nodeIdentifier = generateTubeNodes(fm, nodeIdentifier, leftHornCoordinates, elementsCountInLeftHorn,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=0,
                                           omitEndRows=1, startNodes=None)

        # Create bifurcation nodes
        paCentre = sx_cervix_group[0][1]
        c1Centre = sx_right_horn_group[0][-2]
        c2Centre = sx_left_horn_group[0][-2]
        paxList = cFirstRingNodeCoordinates[0]
        pad1List = cFirstRingNodeCoordinates[1]
        pad2List = cFirstRingNodeCoordinates[2]
        pad2 = cFirstRingNodeCoordinates[2]
        c1xList = rhLastRingNodeCoordinates[0]
        c1d2 = rhLastRingNodeCoordinates[2]
        c2xList = lhLastRingNodeCoordinates[0]
        c2d2 = lhLastRingNodeCoordinates[2]
        nodeIdentifier, roNodeId, coNodeId = create3dBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2,
                                                                      c1Centre, c1xList, c1d2, c2Centre, c2xList, c2d2,
                                                                      elementsCountThroughWall)

        # Create cervix nodes
        nodeCount = nodeIdentifier
        nodeIdentifier = generateTubeNodes(fm, nodeIdentifier, cervixCoordinates, elementsCountInCervix,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                           omitEndRows=0, startNodes=None)

        # Create vagina nodes
        nodeIdentifier = generateTubeNodes(fm, nodeIdentifier, vaginaCoordinates, elementsCountInVagina,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                           omitEndRows=0, startNodes=None)

        # Create elements
        # Create right horn elements
        startNodeId = firstNodeIdentifier
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInRightHorn,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=0, omitEndRows=1,
                                                 meshGroups=[rightHornMeshGroup, uterusMeshGroup])

        # Create left horn elements
        startNodeId = rhLastRingNodeId[-1][-1] + 1
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInLeftHorn,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=0, omitEndRows=1,
                                                 meshGroups=[leftHornMeshGroup, uterusMeshGroup])

        # Create bifurcation elements
        cFirstRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountInCervix,
                                                             elementsCountThroughWall, omitStartRows=1, omitEndRows=0)
        paNodeId = cFirstRingNodeId
        c1NodeId = rhLastRingNodeId
        c2NodeId = lhLastRingNodeId
        elementIdentifier = make_tube_bifurcation_elements_3d(region, coordinates, elementIdentifier,
                                                              elementsCountAround, elementsCountThroughWall, paNodeId,
                                                              c1NodeId, c2NodeId, roNodeId, coNodeId,
                                                              meshGroups=[cervixMeshGroup, rightHornMeshGroup,
                                                                          leftHornMeshGroup, uterusMeshGroup])

        # Create cervix elements
        startNodeId = paNodeId[0][0]
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInCervix,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=1, omitEndRows=0,
                                                 meshGroups=[cervixMeshGroup, uterusMeshGroup])

        # Create vagina elements
        startNodeId = paNodeId[0][0] + (elementsCountInCervix - 1) * elementsCountAround * (elementsCountThroughWall + 1)
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInVagina,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=0, omitEndRows=0,
                                                 meshGroups=[vaginaMeshGroup, uterusMeshGroup])

        fm.endChange()
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

def findDerivativeBetweenPoints(v1, v2):
    """
    Find vector difference between two points and rescale vector difference using cubic hermite arclength
    between the points to derive the derivative between the points.
    :param v1: start vector
    :param v2: end vector
    :return: derivative of between v1 and v2
    """
    d = [v2[c] - v1[c] for c in range(3)]
    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
    d = [c * arcLengthAround for c in vector.normalise(d)]

    return d


class UterusCentralPath:
    """
    Generates sampled central path for uterus scaffold.
    """
    def __init__(self, region, centralPath, targetElementLength):
        """
        :param region: Zinc region to define model in.
        :param centralPath: Central path subscaffold from MeshType_1d_network_layout1
        """

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        tmpFieldmodule = tmpRegion.getFieldmodule()
        tmpNodes = tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        tmpCoordinates = tmpFieldmodule.findFieldByName('coordinates')

        tmpGroup = tmpFieldmodule.findFieldByName('right uterine horn').castGroup()
        cx_rightHorn, cd1_rightHorn, cd2_rightHorn, cd12_rightHorn = get_nodeset_path_field_parameters(
            tmpGroup.getFieldNodeGroup(tmpNodes).getNodesetGroup(), tmpCoordinates,
            [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])

        tmpGroup = tmpFieldmodule.findFieldByName('left uterine horn').castGroup()
        cx_leftHorn, cd1_leftHorn, cd2_leftHorn, cd12_leftHorn = get_nodeset_path_field_parameters(
            tmpGroup.getFieldNodeGroup(tmpNodes).getNodesetGroup(), tmpCoordinates,
            [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])

        tmpGroup = tmpFieldmodule.findFieldByName('uterine cervix').castGroup()
        cx_cervix, cd1_cervix, cd2_cervix, cd12_cervix = get_nodeset_path_field_parameters(
            tmpGroup.getFieldNodeGroup(tmpNodes).getNodesetGroup(), tmpCoordinates,
            [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])

        tmpGroup = tmpFieldmodule.findFieldByName('vagina').castGroup()
        cx_vagina, cd1_vagina, cd2_vagina, cd12_vagina = get_nodeset_path_field_parameters(
            tmpGroup.getFieldNodeGroup(tmpNodes).getNodesetGroup(), tmpCoordinates,
            [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])


        # Find arcLength
        # Right horn
        rightHornLength = 0.0
        rightHornSegmentLengthList = []
        elementsCountInRightHorn = len(cx_rightHorn) - 1
        for e in range(elementsCountInRightHorn):
            arcLength = interp.getCubicHermiteArcLength(cx_rightHorn[e], cd1_rightHorn[e],
                                                        cx_rightHorn[e + 1], cd1_rightHorn[e + 1])
            rightHornSegmentLengthList.append(arcLength)
            rightHornLength += arcLength
        elementsAlongRightHorn = math.ceil(rightHornLength / targetElementLength)

        # Left horn
        leftHornLength = 0.0
        leftHornSegmentLengthList = []
        elementsCountInLeftHorn = len(cx_leftHorn) - 1
        for e in range(elementsCountInLeftHorn):
            arcLength = interp.getCubicHermiteArcLength(cx_leftHorn[e], cd1_leftHorn[e],
                                                        cx_leftHorn[e + 1], cd1_leftHorn[e + 1])
            leftHornSegmentLengthList.append(arcLength)
            leftHornLength += arcLength
        elementsAlongLeftHorn = math.ceil(leftHornLength / targetElementLength)

        # Cervix
        cervixLength = 0.0
        cervixSegmentLengthList = []
        elementsCountInCervix = len(cx_cervix) - 1
        for e in range(elementsCountInCervix):
            arcLength = interp.getCubicHermiteArcLength(cx_cervix[e], cd1_cervix[e],
                                                        cx_cervix[e + 1], cd1_cervix[e + 1])
            cervixSegmentLengthList.append(arcLength)
            cervixLength += arcLength
        elementsAlongCervix = math.ceil(cervixLength / targetElementLength)

        # Vagina
        vaginaLength = 0.0
        vaginaSegmentLengthList = []
        elementsCountInVagina = len(cx_vagina) - 1
        for e in range(elementsCountInVagina):
            arcLength = interp.getCubicHermiteArcLength(cx_vagina[e], cd1_vagina[e],
                                                        cx_vagina[e + 1], cd1_vagina[e + 1])
            vaginaSegmentLengthList.append(arcLength)
            vaginaLength += arcLength
        elementsAlongVagina = math.ceil(vaginaLength / targetElementLength)

        # Sample central path
        sx_rightHorn, sd1_rightHorn, se_rightHorn, sxi_rightHorn, ssf_rightHorn = interp.sampleCubicHermiteCurves(
            cx_rightHorn, cd1_rightHorn, elementsAlongRightHorn)
        sd2_rightHorn, sd12_rightHorn = interp.interpolateSampleCubicHermite(cd2_rightHorn, cd12_rightHorn,
                                                                             se_rightHorn, sxi_rightHorn, ssf_rightHorn)
        # sd3_rightHorn, sd13_rightHorn = interp.interpolateSampleCubicHermite(cd3_rightHorn, cd13_rightHorn,
        #                                                                      se_rightHorn, sxi_rightHorn, ssf_rightHorn)

        sx_leftHorn, sd1_leftHorn, se_leftHorn, sxi_leftHorn, ssf_leftHorn = interp.sampleCubicHermiteCurves(
            cx_leftHorn, cd1_leftHorn, elementsAlongLeftHorn)
        sd2_leftHorn, sd12_leftHorn = interp.interpolateSampleCubicHermite(cd2_leftHorn, cd12_leftHorn, se_leftHorn,
                                                                           sxi_leftHorn, ssf_leftHorn)

        sx_cervix, sd1_cervix, se_cervix, sxi_cervix, ssf_cervix = interp.sampleCubicHermiteCurves(
            cx_cervix, cd1_cervix, elementsAlongCervix)
        sd2_cervix, sd12_cervix = interp.interpolateSampleCubicHermite(cd2_cervix, cd12_cervix, se_cervix, sxi_cervix,
                                                                       ssf_cervix)

        sx_vagina, sd1_vagina, se_vagina, sxi_vagina, ssf_vagina = interp.sampleCubicHermiteCurves(
            cx_vagina, cd1_vagina, elementsAlongVagina)
        sd2_vagina, sd12_vagina = interp.interpolateSampleCubicHermite(cd2_vagina, cd12_vagina, se_vagina, sxi_vagina,
                                                                       ssf_vagina)

        sx_right_horn_group = [sx_rightHorn, sd1_rightHorn, sd2_rightHorn, sd12_rightHorn]
        sx_left_horn_group = [sx_leftHorn, sd1_leftHorn, sd2_leftHorn, sd12_leftHorn]
        sx_cervix_group = [sx_cervix, sd1_cervix, sd2_cervix, sd12_cervix]
        sx_vagina_group = [sx_vagina, sd1_vagina, sd2_vagina, sd12_vagina]

        del tmpGroup
        del tmpCoordinates
        del tmpNodes
        del tmpFieldmodule
        del tmpRegion

        # Find nodes of uterus cervix and horn
        self.cx_right_horn_group = [cx_rightHorn, cd1_rightHorn, cd2_rightHorn, cd12_rightHorn]
        self.cx_left_horn_group = [cx_leftHorn, cd1_leftHorn, cd2_leftHorn, cd12_leftHorn]
        self.cx_cervix_group = [cx_cervix, cd1_cervix, cd2_cervix, cd12_cervix]
        self.cx_vagina_group = [cx_vagina, cd1_vagina, cd2_vagina, cd12_vagina]
        self.rightHornLength = rightHornLength
        self.rightHornSegmentLengthList = rightHornSegmentLengthList
        self.leftHornLength = leftHornLength
        self.leftHornSegmentLengthList = leftHornSegmentLengthList
        self.cervixLength = cervixLength
        self.cervixSegmentLengthList = cervixSegmentLengthList
        self.vaginaLength = vaginaLength
        self.vaginaSegmentLengthList = vaginaSegmentLengthList

        self.sx_right_horn_group = sx_right_horn_group
        self.sx_left_horn_group = sx_left_horn_group
        self.sx_cervix_group = sx_cervix_group
        self.sx_vagina_group = sx_vagina_group


# class UterusCentralPath:
#     """
#     Generates sampled central path for uterus scaffold.
#     """
#     def __init__(self, region, centralPath, targetElementLength):
#         """
#         :param region: Zinc region to define model in.
#         :param centralPath: Central path subscaffold from MeshType_1d_network_layout1
#         """
#
#         # Central path
#         tmpRegion = region.createRegion()
#         centralPath.generate(tmpRegion)
#         cx_rightHorn, cd1_rightHorn, cd2_rightHorn, cd12_rightHorn = \
#             extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
#                                                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
#                                             groupName='right uterine horn')
#         cx_leftHorn, cd1_leftHorn, cd2_leftHorn, cd12_leftHorn = \
#             extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
#                                                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
#                                             groupName='left uterine horn')
#         # should be fixed here (find a general way that works for any segment sequences)
#         # Remove the first node of the left horn and add it to the end of the list
#         cx_leftHorn.append(cx_leftHorn[0])
#         cd1_leftHorn.append(cd1_leftHorn[0])
#         cd2_leftHorn.append(cd2_leftHorn[0])
#         cd12_leftHorn.append(cd12_leftHorn[0])
#         cx_leftHorn = cx_leftHorn[1:]
#         cd1_leftHorn = cd1_leftHorn[1:]
#         cd2_leftHorn = cd2_leftHorn[1:]
#         cd12_leftHorn = cd12_leftHorn[1:]
#
#         # [[-6.82, -0.07, -2.67], [6.82, -0.07, -2.67], [0.00, 0.00, -2.00]]
#         cd1_leftHorn = cd1_leftHorn[:-1]
#         cd1_leftHorn.append([6.82, -0.07, -2.67])
#         cd2_leftHorn = cd2_leftHorn[:-1]
#         cd2_leftHorn.append([0.85, -0.02, 0.53])
#
#         cx_cervix, cd1, cd2, cd12_cervix = \
#             extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
#                                                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
#                                             groupName='uterine cervix')
#
#         cd1_cervix = []
#         cd1_cervix.append([0.00, 0.00, -2.00])
#         cd1_cervix += cd1[1:]
#         cd2_cervix = []
#         cd2_cervix.append([1.00, 0.00, 0.00])
#         cd2_cervix += cd2[1:]
#
#         cx_vagina, cd1_vagina, cd2_vagina, cd12_vagina = \
#             extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
#                                                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
#                                             groupName='vagina')
#
#         # Find arcLength
#         # Right horn
#         rightHornLength = 0.0
#         rightHornSegmentLengthList = []
#         elementsCountInRightHorn = len(cx_rightHorn) - 1
#         for e in range(elementsCountInRightHorn):
#             arcLength = interp.getCubicHermiteArcLength(cx_rightHorn[e], cd1_rightHorn[e],
#                                                         cx_rightHorn[e + 1], cd1_rightHorn[e + 1])
#             rightHornSegmentLengthList.append(arcLength)
#             rightHornLength += arcLength
#         elementsAlongRightHorn = math.ceil(rightHornLength / targetElementLength)
#
#         # Left horn
#         leftHornLength = 0.0
#         leftHornSegmentLengthList = []
#         elementsCountInLeftHorn = len(cx_leftHorn) - 1
#         for e in range(elementsCountInLeftHorn):
#             arcLength = interp.getCubicHermiteArcLength(cx_leftHorn[e], cd1_leftHorn[e],
#                                                         cx_leftHorn[e + 1], cd1_leftHorn[e + 1])
#             leftHornSegmentLengthList.append(arcLength)
#             leftHornLength += arcLength
#         elementsAlongLeftHorn = math.ceil(leftHornLength / targetElementLength)
#
#         # Cervix
#         cervixLength = 0.0
#         cervixSegmentLengthList = []
#         elementsCountInCervix = len(cx_cervix) - 1
#         for e in range(elementsCountInCervix):
#             arcLength = interp.getCubicHermiteArcLength(cx_cervix[e], cd1_cervix[e],
#                                                         cx_cervix[e + 1], cd1_cervix[e + 1])
#             cervixSegmentLengthList.append(arcLength)
#             cervixLength += arcLength
#         elementsAlongCervix = math.ceil(cervixLength / targetElementLength)
#
#         # Vagina
#         vaginaLength = 0.0
#         vaginaSegmentLengthList = []
#         elementsCountInVagina = len(cx_vagina) - 1
#         for e in range(elementsCountInVagina):
#             arcLength = interp.getCubicHermiteArcLength(cx_vagina[e], cd1_vagina[e],
#                                                         cx_vagina[e + 1], cd1_vagina[e + 1])
#             vaginaSegmentLengthList.append(arcLength)
#             vaginaLength += arcLength
#         elementsAlongVagina = math.ceil(vaginaLength / targetElementLength)
#
#         # Sample central path
#         sx_rightHorn, sd1_rightHorn, se_rightHorn, sxi_rightHorn, ssf_rightHorn = interp.sampleCubicHermiteCurves(
#             cx_rightHorn, cd1_rightHorn, elementsAlongRightHorn)
#         sd2_rightHorn, sd12_rightHorn = interp.interpolateSampleCubicHermite(cd2_rightHorn, cd12_rightHorn,
#                                                                              se_rightHorn, sxi_rightHorn, ssf_rightHorn)
#         # sd3_rightHorn, sd13_rightHorn = interp.interpolateSampleCubicHermite(cd3_rightHorn, cd13_rightHorn,
#         #                                                                      se_rightHorn, sxi_rightHorn, ssf_rightHorn)
#
#         sx_leftHorn, sd1_leftHorn, se_leftHorn, sxi_leftHorn, ssf_leftHorn = interp.sampleCubicHermiteCurves(
#             cx_leftHorn, cd1_leftHorn, elementsAlongLeftHorn)
#         sd2_leftHorn, sd12_leftHorn = interp.interpolateSampleCubicHermite(cd2_leftHorn, cd12_leftHorn, se_leftHorn,
#                                                                            sxi_leftHorn, ssf_leftHorn)
#
#         sx_cervix, sd1_cervix, se_cervix, sxi_cervix, ssf_cervix = interp.sampleCubicHermiteCurves(
#             cx_cervix, cd1_cervix, elementsAlongCervix)
#         sd2_cervix, sd12_cervix = interp.interpolateSampleCubicHermite(cd2_cervix, cd12_cervix, se_cervix, sxi_cervix,
#                                                                        ssf_cervix)
#
#         sx_vagina, sd1_vagina, se_vagina, sxi_vagina, ssf_vagina = interp.sampleCubicHermiteCurves(
#             cx_vagina, cd1_vagina, elementsAlongVagina)
#         sd2_vagina, sd12_vagina = interp.interpolateSampleCubicHermite(cd2_vagina, cd12_vagina, se_vagina, sxi_vagina,
#                                                                        ssf_vagina)
#
#         sx_right_horn_group = [sx_rightHorn, sd1_rightHorn, sd2_rightHorn, sd12_rightHorn]
#         sx_left_horn_group = [sx_leftHorn, sd1_leftHorn, sd2_leftHorn, sd12_leftHorn]
#         sx_cervix_group = [sx_cervix, sd1_cervix, sd2_cervix, sd12_cervix]
#         sx_vagina_group = [sx_vagina, sd1_vagina, sd2_vagina, sd12_vagina]
#
#         del tmpRegion
#
#         # Find nodes of uterus cervix and horn
#         self.cx_right_horn_group = [cx_rightHorn, cd1_rightHorn, cd2_rightHorn, cd12_rightHorn]
#         self.cx_left_horn_group = [cx_leftHorn, cd1_leftHorn, cd2_leftHorn, cd12_leftHorn]
#         self.cx_cervix_group = [cx_cervix, cd1_cervix, cd2_cervix, cd12_cervix]
#         self.cx_vagina_group = [cx_vagina, cd1_vagina, cd2_vagina, cd12_vagina]
#         self.rightHornLength = rightHornLength
#         self.rightHornSegmentLengthList = rightHornSegmentLengthList
#         self.leftHornLength = leftHornLength
#         self.leftHornSegmentLengthList = leftHornSegmentLengthList
#         self.cervixLength = cervixLength
#         self.cervixSegmentLengthList = cervixSegmentLengthList
#         self.vaginaLength = vaginaLength
#         self.vaginaSegmentLengthList = vaginaSegmentLengthList
#
#         self.sx_right_horn_group = sx_right_horn_group
#         self.sx_left_horn_group = sx_left_horn_group
#         self.sx_cervix_group = sx_cervix_group
#         self.sx_vagina_group = sx_vagina_group

def getCoordinatesAlongTube3D(sx_group, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                              wallThickness, startRadian):

    d3List = []
    for n in range(len(sx_group[0])):
        d3Dir = vector.normalise(vector.crossproduct3(sx_group[1][n], sx_group[2][n]))
        d3Mag = vector.magnitude(sx_group[2][n])
        v = vector.setMagnitude(d3Dir, d3Mag)
        d3List.append(v)

    # Create ellipses along tube around the central path
    xEllipsesAlong = []
    d1EllipsesAlong = []
    for n in range(len(sx_group[0])):
        px, pd1 = createEllipsePoints(sx_group[0][n], 2 * math.pi, sx_group[2][n], d3List[n], elementsCountAround,
                                      startRadians=startRadian)
        xEllipsesAlong.append(px)
        d1EllipsesAlong.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong) - 1):
            v1 = xEllipsesAlong[n2][n1]
            v2 = xEllipsesAlong[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipsesAlong[-1][n1])
        d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2EllipsesAlong = []
    for n2 in range(len(xEllipsesAlong)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2EllipsesAlong.append(d2Around)

    # Spread out elements along tube
    xRaw = []
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong)):
            xAlong.append(xEllipsesAlong[n2][n1])
            d2Along.append(d2EllipsesAlong[n2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongTube,
                                                                        arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledTube = []
    d1SampledTube = []
    d2SampledTube = []
    for n2 in range(elementsCountAlongTube + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
            # if n2 > 0:
            #     v1 = xRaw[n1][n2]
            #     v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 2 else 0][n2]
            #     d1 = findDerivativeBetweenPoints(v1, v2)
            #     d1Around.append(d1)
            # else:
            #     d1Around.append(d2Raw[n2][0])
        # if n2 > 0:
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        # else:
        #     d1Smoothed = d1Around
        xSampledTube.append(xAround)
        d1SampledTube.append(d1Smoothed)
        d2SampledTube.append(d2Around)

    d3Tube = []
    for n2 in range(elementsCountAlongTube + 1):
        d3Around = []
        for n1 in range(elementsCountAround):
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1SampledTube[n2][n1]), vector.normalise(d2SampledTube[n2][n1]))))
        d3Tube.append(d3Around)

    xInner = []
    d1Inner = []
    d2Inner = []
    d3Inner = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            xInner.append(xSampledTube[n2][n1])
            d1Inner.append(d1SampledTube[n2][n1])
            d2Inner.append(d2SampledTube[n2][n1])
            d3Inner.append(d3Tube[n2][n1])

    transitElementList = [0] * elementsCountAround
    # if elementsCountThroughWall == 2:
    #     relativeThicknessList = [0.3, 0.7]
    # else:
    #     relativeThicknessList = []
    relativeThicknessList = []
    xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xInner, d1Inner,
        d2Inner, d3Inner, [wallThickness]*(elementsCountAlongTube+1), relativeThicknessList,
        elementsCountAround, elementsCountAlongTube, elementsCountThroughWall, transitElementList)

    coordinatesList = [xList, d1List, d2List, d3List]

    return coordinatesList

def generateTubeNodes(fm, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                      elementsCountThroughWall, omitStartRows, omitEndRows, startNodes=None):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        if n2 == 0:
            startNodeId = nodeIdentifier - 1
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeIdentifier)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        nodeIdentifier += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeIdentifier)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        nodeIdentifier += 1
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    nodeIdentifier += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    d3FirstRing.append(d3FirstRingThroughWall)

    return nodeIdentifier

def generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountAlongTube, elementsCountAround,
                         elementsCountThroughWall, useCrossDerivatives, omitStartRows, omitEndRows, startNodes=None,
                         meshGroups=None):

    cache = fm.createFieldcache()
    mesh = fm.findMeshByDimension(3)
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    # Create tube elements
    now = elementsCountAround * (elementsCountThroughWall + 1)
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlongTube - 1 if omitStartRows or omitEndRows == 1 else elementsCountAlongTube):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = e2 * now + e3 * elementsCountAround + e1 + startNodeId
                bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + startNodeId
                bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                # print('nodeIdentifiers', nodeIdentifiers)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

    # if omitStartRows == 1 or omitEndRows == 1:
    #     lastNodeId = elementsCountAlongTube * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    # else:
    #     lastNodeId = (elementsCountAlongTube + 1) * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    return elementIdentifier


def make_tube_bifurcation_elements_3d(region, coordinates, elementIdentifier, elementsCountAround,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    fm = region.getFieldmodule()
    mesh = fm.findMeshByDimension(3)
    eftfactory = eftfactory_bicubichermitelinear(mesh, None)
    eftStd = eftfactory.createEftBasic()

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStd.defineField(coordinates, -1, eftStd)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    # Right tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c1Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = coNodeId[e3][e1 - elementsCountAround // 2]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = coNodeId[e3][e1 - elementsCountAround // 2 - 1]
                bni5 = bni1 + elementsCountAround
                if e1 == c1Count - 1:
                    bni4 = roNodeId[e3][0]
                    bni6 = bni5 - elementsCountAround + 1
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 in (0, pac1Count - 1, pac1Count, c1Count - 1):
                eft = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == c1Count - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
               if meshGroups.index(meshGroup) == 1:
                    meshGroup.addElement(element)
               elif meshGroups.index(meshGroup) == 3:
                   meshGroup.addElement(element)

    # Left tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c2Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround//2:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                if e1 == 0:
                    bni3 = roNodeId[e3][e1]
                    bni4 = coNodeId[e3][-1 - e1]
                else:
                    bni3 = coNodeId[e3][-e1]
                    if e1 == elementsCountAround//2 - 1:
                        bni4 = roNodeId[e3][0] + pac1Count
                    else:
                        bni4 = bni3 - 1
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1)%c2Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 <= c1c2Count:
                eft = eftfactory.createEftBasic()
                scalefactors = [ -1.0 ]
                setEftScaleFactorIds(eft, [1], [])
                if e1 == 0:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    scaleEftNodeValueLabels(eft, [4, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 < (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 4, 7, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 == (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 7], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == c1c2Count:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            elif e1 == c2Count - 1:
                eft = eftfactory.createEftBasic()
                remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 2:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(paCount):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = roNodeId[e3][e1]
            bni2 = roNodeId[e3][(e1 + 1) % elementsCountAround]
            bni3 = paNodeId[e3][e1]
            bni4 = paNodeId[e3][(e1 + 1) % elementsCountAround]
            bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
            bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
            bni7 = bni3 + elementsCountAround
            bni8 = bni4 + elementsCountAround
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            if e1 in (0, pac1Count):
                eft = eftfactory.createEftBasic()
                if e1 in (0, pac1Count):
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)
    # lastNodeId = startNodeIndex + (elementsCountThroughWall + 1) * (len(roNodeId) + len(coNodeId)) + 1
    # lastNodeId = paNodeId[0][0]
    return elementIdentifier


def getTargetedRingNodesCoordinates(tubeCoordinates, elementsCountAround, elementsCountAlongTube,
                                    elementsCountThroughWall, omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            # lastRingNodeIdThroughWall.append(nodeCount)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        # nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            # firstRingNodeIdThroughWall.append(nodeCount)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        # nodeCount += 1
                # else:
                #     nodeCount += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    d3FirstRing.append(d3FirstRingThroughWall)

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing, d3FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing, d3LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                           omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    firstRingsNodeId = []
    for n2 in range(elementsCountAlongTube + 1):
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            firstRingNodeIdThroughWall = []
            for n1 in range(elementsCountAround):
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                else:
                    nodeCount += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)

    if omitStartRows == 1:
        targetedRingNodeId = firstRingsNodeId
    elif omitEndRows == 1:
        targetedRingNodeId = lastRingsNodeId
    else:
        targetedRingNodeId = []

    return targetedRingNodeId, nodeCount


def create3dBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList, c1d2, c2Centre, c2xList,
                             c2d2, elementsCountThroughWall):

    cache = fm.createFieldcache()
    mesh = fm.findMeshByDimension(3)
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    roxList = []
    rod1List = []
    rod2List = []
    coxList = []
    cod1List = []
    cod2List = []
    for n3 in range(elementsCountThroughWall + 1):
        roxOuter, rod1Outer, rod2Outer, coxOuter, cod1Outer, cod2Outer, paStartIndex, c1StartIndex, c2StartIndex = \
            make_tube_bifurcation_points_converging(paCentre, paxList[n3], pad2[n3], c1Centre, c1xList[n3], c1d2[n3],
                                                    c2Centre, c2xList[n3], c2d2[n3])
        roxList.append(roxOuter)
        rod1List.append(rod1Outer)
        rod2List.append(rod2Outer)
        coxList.append(coxOuter)
        cod1List.append(cod1Outer)
        cod2List.append(cod2Outer)

    # Create bifurcation nodes
    roNodeId = []
    coNodeId = []
    for n3 in range(elementsCountThroughWall + 1):
        coNodeIdThroughWall = []
        for n in range(len(coxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, coxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2List[n3][n])
            coNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        coNodeId.append(coNodeIdThroughWall)
        roNodeIdThroughWall = []
        for n in range(len(roxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, roxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2List[n3][n])
            roNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        roNodeId.append(roNodeIdThroughWall)

    return nodeIdentifier, roNodeId, coNodeId


def make_tube_bifurcation_points_converging(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    '''
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2
    '''
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
    # convert to number of nodes, includes both 6-way points
    pac1NodeCount = pac1Count + 1
    pac2NodeCount = pac2Count + 1
    c1c2NodeCount = c1c2Count + 1
    paStartIndex = 0
    c1StartIndex = 0
    c2StartIndex = 0
    pac1x  = [ None ]*pac1NodeCount
    pac1d1 = [ None ]*pac1NodeCount
    pac1d2 = [ None ]*pac1NodeCount
    for n in range(pac1NodeCount):
        pan = (paStartIndex + n) % paCount
        c1n = (c1StartIndex + n) % c1Count
        # x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac1x [n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac1d1[n] = [ 0.0, 0.0, 0.0 ]
        pac1d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    paStartIndex2 = paStartIndex + pac1Count
    c1StartIndex2 = c1StartIndex + pac1Count
    c2StartIndex2 = c2StartIndex + c1c2Count
    pac2x  = [ None ]*pac2NodeCount
    pac2d1 = [ None ]*pac2NodeCount
    pac2d2 = [ None ]*pac2NodeCount
    for n in range(pac2NodeCount):
        pan = (paStartIndex2 + n) % paCount
        c2n = (c2StartIndex2 + n) % c2Count
        # x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c2x[c2n], mult(c2d2[c2n], 2.0)
        x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac2x [n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac2d1[n] = [ 0.0, 0.0, 0.0 ]
        pac2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    c1c2x  = [ None ]*c1c2NodeCount
    c1c2d1 = [ None ]*c1c2NodeCount
    c1c2d2 = [ None ]*c1c2NodeCount
    for n in range(c1c2NodeCount):
        c1n = (c1StartIndex2 + n) % c1Count
        c2n = (c2StartIndex2 - n) % c2Count  # note: reversed
        # x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], -2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), c2x[c2n], mult(c2d2[c2n], -2.0)
        c1c2x [n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        c1c2d1[n] = [ 0.0, 0.0, 0.0 ]
        c1c2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    # get hex triple points
    # hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
    #     pax[paStartIndex], mult(pad2[paStartIndex], -1.0),
    #     c1x[c1StartIndex], c1d2[c1StartIndex],
    #     c2x[c1StartIndex], c2d2[c2StartIndex])
    # hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
    #     pax[paStartIndex2], mult(pad2[paStartIndex2], -1.0),
    #     c2x[c2StartIndex2], c2d2[c2StartIndex2],
    #     c1x[c1StartIndex2], c1d2[c1StartIndex2])
    hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
        c2x[c1StartIndex], mult(c2d2[c2StartIndex], -1.0),
        c1x[c1StartIndex], mult(c1d2[c1StartIndex], -1.0),
        pax[paStartIndex], pad2[paStartIndex])
    hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
        c1x[c1StartIndex2], mult(c1d2[c1StartIndex2], -1.0),
        c2x[c2StartIndex2], mult(c2d2[c2StartIndex2], -1.0),
        pax[paStartIndex2], pad2[paStartIndex2])
    # smooth around loops through hex points to get d1
    loop1x  = [ hex2 ] + pac2x[1:-1] + [ hex1 ]
    loop1d1 = [ [ -d for d in hex2d2 ] ] + pac2d1[1:-1] + [ hex1d1 ]
    loop2x  = [ hex1 ] + pac1x[1:-1] + [ hex2 ]
    loop2d1 = [ [ -d for d in hex1d2 ] ] + pac1d1[1:-1] + [ hex2d1 ]
    loop1d1 = interp.smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    loop2d1 = interp.smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # smooth over "crotch" between c1 and c2
    crotchx = [ hex2 ] + c1c2x[1:-1] + [ hex1 ]
    crotchd1 = [ add(hex2d1, hex2d2) ] + c1c2d1[1:-1] + [ [ (-hex1d1[c] - hex1d2[c]) for c in range(3) ] ]
    crotchd1 = interp.smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True,
                                                        fixEndDerivative=True,
                                                        magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    rox  = [ hex1 ] + pac1x[1:-1] + [ hex2 ] + pac2x[1:-1]
    rod1 = [ loop1d1[-1] ] + loop2d1[1:] + loop1d1[1:-1]
    rod2 = [ [ -d for d in loop2d1[ 0] ] ] + pac1d2[1:-1] + [ [ -d for d in loop1d1[0] ] ] + pac2d2[1:-1]
    cox  = crotchx [1:-1]
    cod1 = crotchd1[1:-1]
    cod2 = c1c2d2[1:-1]
    return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex