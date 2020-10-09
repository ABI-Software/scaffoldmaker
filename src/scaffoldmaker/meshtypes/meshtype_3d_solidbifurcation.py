"""
Generates a solid cylinder using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor and length directions.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, ConeBaseProgression, Tapered, \
    CylinderEnds, CylinderCentralPath, Ellipse2D, EllipseShape
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from opencmiss.zinc.node import Node
from opencmiss.zinc.field import Field
from scaffoldmaker.utils.shieldmesh import ShieldMesh, ShieldShape, ShieldRimDerivativeMode
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite,\
    smoothCubicHermiteDerivativesLine, getCubicHermiteArcLength,DerivativeScalingMode, sampleParameterAlongLine
from scaffoldmaker.utils import centralpath
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element


class MeshType_3d_solidbifurcation(Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor and length directions.
    """
    centralPathDefaultScaffoldPackages = {
        'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.5]],
                    [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.5]]])
        })
    }

    @staticmethod
    def getName():
        return '3D Solid Bifurcation 1'


    @classmethod
    def getDefaultOptions(cls,parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Use central path': False,
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements along': 1,
            'Lower half': False,
            'Length': 1.0,
            'Major radius': 1.0,
            'Major radius geometric progression change': True,
            'Major radius end ratio': 0.92,
            'Minor radius': 1.0,
            'Minor radius geometric progression change': True,
            'Minor radius end ratio': 0.92,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Use central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements along',
            'Lower half',
            'Length',
            'Major radius',
            'Major radius geometric progression change',
            'Major radius end ratio',
            'Minor radius',
            'Minor radius geometric progression change',
            'Minor radius end ratio',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls,options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1
        if options['Number of elements along'] < 1:
            options['Number of elements along'] = 1

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        centralPath = options['Central path']
        useCentralPath = options['Use central path']
        full = not options['Lower half']
        length = options['Length']
        majorRadius = options['Major radius']
        majorGeometric = options['Major radius geometric progression change']
        minorGeometric = options['Minor radius geometric progression change']
        minorRadius = options['Minor radius']
        majorRadiusEndRatio = options['Major radius end ratio']
        minorRadiusEndRatio = options['Minor radius end ratio']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = options['Use cross derivatives']

        if useCentralPath:
            cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)
        else:
            cylinderCentralPath = None

            # segmentLength = length / segmentCount
            # elementAlongLength = length / elementsCountAlong
            # # print('Length = ', length)

            # Sample central path
            # sx, sd1, se, sxi, ssf = sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment * segmentCount)
            # sd2, sd12 = interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

            # # Generate variation of radius & tc width along length
            # lengthList = [0.0, duodenumLength, duodenumLength + jejunumLength, length]
            # innerRadiusList = [duodenumInnerRadius, duodenumJejunumInnerRadius, jejunumIleumInnerRadius, ileumInnerRadius]
            # innerRadiusSegmentList, dInnerRadiusSegmentList = sampleParameterAlongLine(lengthList, innerRadiusList,
            #                                                                                   segmentCount)

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fm.findMeshByDimension(3)

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]

        ellipseBRL = createBranchEllipseLeft()

        ellipseBRU = createBranchEllipseUp()

        ellipseBL = createBaseEllipseLeft()

        ellipseDL = createDaughterEllipseLeft()
        elementsCountAlong = 1
        shieldDL = copyEllipsesNodesToShieldNodes(elementsCountAlong, ellipseBRU, ellipseBRL, ellipseDL)

        # now pd2 derivatives.
        for n3 in range(shieldDL.elementsCountAlong + 1):
            calculateD2Derivatives(shieldDL, n3, 1)
        smoothd2Derivatives(shieldDL)
        skipNodes = createNodeId(shieldDL, nodes)
        generateNodes(shieldDL, nodes, fm, coordinates, skipNodes)
        generateElements(shieldDL, mesh, fm, coordinates)

        ellipseBR = createBaseEllipseRight()
        ellipseBRR = createBranchEllipseRight()
        ellipseBRUR = createBranchEllipseUpRight()
        ellipseDR = createDaughterEllipseRight()
        shieldDR = copyEllipsesNodesToShieldNodes(elementsCountAlong, ellipseBRUR, ellipseBRR, ellipseDR)
        for n3 in range(shieldDR.elementsCountAlong + 1):
            calculateD2Derivatives(shieldDR, n3, 1)
        smoothd2Derivatives(shieldDR)
        skipNodes = createNodeId(shieldDR, nodes)
        generateNodes(shieldDR, nodes, fm, coordinates, skipNodes)
        generateElements(shieldDR, mesh, fm, coordinates)

        shieldInB = createCylinderInBetween(ellipseBRU, ellipseBRUR)
        skipNodes = setNodeIdForInBetween(shieldInB,shieldDL,shieldDR, nodes)
        # generateNodes(shield, nodes, fm, coordinates, skipNodes)
        generateElementsInBetween(shieldInB, mesh, fm, coordinates)

        ellipseBASER = createBaseEllipseRight()
        shieldBR = createCylinderRight(ellipseBASER, ellipseBRR)
        skipNodes = createNodeIdForCylRight(shieldBR,shieldDR,nodes)
        # skipNodes = createNodeId(shieldBR, nodes)
        generateNodes(shieldBR, nodes, fm, coordinates, skipNodes)
        # generateElementsBaseRight(shieldBR, mesh, fm, coordinates)

        # ellipseBLF = createBaseFullLeft()
        # t=1
        # shieldBL = createBaseCylinderLeft(elementsCountAlong, ellipseBL, ellipseBRL)
        # skipNodes = createNodeId(shieldBL, nodes)
        # generateNodes(shieldBL, nodes, fm, coordinates, skipNodes)
        # generateElements(shieldBL, mesh, fm, coordinates)




        # if not useCentralPath:
        #     majorRatio, majorProgression = radiusChange(majorRadius, majorRadiusEndRatio, elementsCountAlong, geometric=majorGeometric)
        #     minorRatio, minorProgression = radiusChange(minorRadius, minorRadiusEndRatio, elementsCountAlong, geometric=minorGeometric)
        #     radiusChanges = Tapered(majorRatio, majorProgression, minorRatio, minorProgression)
        # else:
        #     radiusChanges = []

        # cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF
        #
        # base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, [0.0, 0.0, 0.0],
        #                     vector.setMagnitude(axis3, length), vector.setMagnitude(axis1, majorRadius), minorRadius)
        # cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
        #                          cylinderShape=cylinderShape, tapered=radiusChanges,
        #                          cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # shieldD1, skipNodes = createDaughterCylinder(nodes, axis1, axis2, axis3, 45, 90, 45, -45)

        # generateParentCylinder

        #####################
        # nodeIdentifier = 1
        # startNodeIdentifier = nodeIdentifier
        # generateNodes(shieldD1, nodes, fm, coordinates, skipNodes)
        # generateElements(shieldD1, mesh, fm, coordinates)
        #
        # shieldD2, skipNodes = createDaughterCylinder(nodes, axis1, axis2, axis3, 135, 90, 135, -135)
        # modifyNodeId(shieldD1,shieldD2,range(1), range(shieldD2.elementsCountUp,shieldD2.elementsCountUpFull+1),
        #              range(shieldD2.elementsCountAcross+1))
        # #####################
        # # nodeIdentifier = 1
        # # startNodeIdentifier = nodeIdentifier
        # # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # # mesh = fm.findMeshByDimension(3)
        # generateNodes(shieldD2, nodes, fm, coordinates, skipNodes)
        # generateElements(shieldD2, mesh, fm, coordinates)

        annotationGroup = []
        return annotationGroup

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)


def createNodeIdForCylRight(shieldBR, shieldDR, nodes):
    nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
    skipNodes = [[] for _ in range(shieldBR.elementsCountAlong + 1)]
    for n3 in range(shieldBR.elementsCountAlong + 1):
        for n2 in range(shieldBR.elementsCountUpFull + 1):
            for p in [skipNodes[n3]]:
                p.append([None] * (shieldBR.elementsCountAcross + 1))
    for n2 in range(shieldBR.elementsCountUpFull + 1):
        for n1 in range(shieldBR.elementsCountAcross + 1):
            shieldBR.nodeId[1][n2][n1] = shieldDR.nodeId[0][4-n2][n1]
            skipNodes[1][n2][n1] = 1
            if shieldBR.px[0][n2][n1]:
                shieldBR.nodeId[0][n2][n1] = nodeIdentifier
                nodeIdentifier += 1



    # for n2 in n2range:
    #     for n3 in n3range:
    #         for n1 in n1range:
    #             if shield2.px[n3][n2][n1]:
    #                 shield2.nodeId[n3][n2][n1] = shield1.nodeId[n3][n2][n1]
    #                 skipNodes[n3][n2][n1] = 1
    return skipNodes


def createCylinderRight(ellipseBASER, ellipseBRR):
    elementsCountRim = 0
    elementsCountAlong = 1
    shieldMode = ShieldShape.SHIELD_SHAPE_LOWER_HALF
    shield = ShieldMesh(ellipseBASER.elementsCountAcrossMinor, ellipseBASER.elementsCountAcrossMajor, elementsCountRim,
                        None, elementsCountAlong, shieldMode,
                        shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

    shield.px[0] = ellipseBASER.px
    shield.pd1[0] = ellipseBASER.pd1
    shield.pd3[0] = ellipseBASER.pd3

    shield.px[1] = ellipseBRR.px
    shield.pd1[1] = ellipseBRR.pd1
    shield.pd3[1] = ellipseBRR.pd3
    for n3 in range(shield.elementsCountAlong + 1):
        calculateD2Derivatives(shield, n3, 1)
    smoothd2Derivatives(shield)

    return shield


def createCylinderInBetween(ellipseL,ellipseR):
    elementsCountRim = 0
    elementsCountAlong = 1
    shieldMode = ShieldShape.SHIELD_SHAPE_LOWER_HALF
    shield = ShieldMesh(ellipseL.elementsCountAcrossMinor, ellipseL.elementsCountAcrossMajor, elementsCountRim,
                        None, elementsCountAlong, shieldMode,
                        shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

    shield.px[0] = ellipseL.px
    shield.pd1[0] = ellipseL.pd1
    shield.pd3[0] = ellipseL.pd3

    shield.px[1] = ellipseR.px
    shield.pd1[1] = ellipseR.pd1
    shield.pd3[1] = ellipseR.pd3
    for n3 in range(shield.elementsCountAlong + 1):
        calculateD2Derivatives(shield, n3, 1)
    smoothd2Derivatives(shield)

    return shield

def createBaseCylinderLeft(elementsCountAlong, ellipseB, ellipseBR):
    # now pd2 derivatives.
    elementsCountRim = 0
    shieldMode = ShieldShape.SHIELD_SHAPE_LOWER_HALF
    shield = ShieldMesh(ellipseB.elementsCountAcrossMinor, ellipseB.elementsCountAcrossMajor, elementsCountRim,
                        None, elementsCountAlong, shieldMode,
                        shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

    shield.px[0] = ellipseB.px
    shield.pd1[0] = ellipseB.pd1
    shield.pd3[0] = ellipseB.pd3

    shield.px[1] = ellipseBR.px
    shield.pd1[1] = ellipseBR.pd1
    shield.pd3[1] = ellipseBR.pd3
    for n3 in range(shield.elementsCountAlong + 1):
        calculateD2Derivatives(shield, n3, 1)
    smoothd2Derivatives(shield)

    return shield



def generateEllipseNodes(ellipse, fieldmodule, coordinates, startNodeIdentifier):
    """
    Create shield nodes from coordinates.
    :param fieldmodule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
    :param coordinates: Coordinate field to define.
    :param startNodeIdentifier: First node identifier to use.
    :param mirrorPlane: mirror plane ax+by+cz=d in form of [a,b,c,d]
    :return: next nodeIdentifier.
     """
    nodeIdentifier = startNodeIdentifier
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    cache = fieldmodule.createFieldcache()


    for n2 in range(ellipse.elementsCountAcrossMajor + 1):
        for n1 in range(ellipse.elementsCountAcrossMinor + 1):
            if ellipse.px[n2][n1]:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                ellipse.nodeId[n2][n1] = nodeIdentifier
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ellipse.px[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ellipse.pd1[n2][n1])
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ellipse.pd2[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, ellipse.pd3[n2][n1])
                nodeIdentifier += 1

    return nodeIdentifier


def createBaseFullLeft():
    centre = [-0.25, 0.0, -3.0]
    majorAxis = [1.0, 0.0, 0.0]
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 4
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

    return ellipse


def createBaseEllipseLeft():
    centre = [-0.25, 0.0, -3.0]
    majorAxis = [1.0, 0.0, 0.0]
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_UPPER_HALF)

    return ellipse


def createBaseEllipseRight():
    centre = [0.25, 0.0, -3.0]
    majorAxis = [-1.0, 0.0, 0.0]
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_UPPER_HALF)

    return ellipse


def createBranchEllipseLeft():
    centre = [-0.25, 0.0, 0.0]
    axis1 = [1.0, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    theta = -45
    majorAxis = rotateVector(axis1, axis2, theta)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_UPPER_HALF)

    return ellipse


def createBranchEllipseRight():
    centre = [0.25, 0.0, 0.0]
    axis1 = [-1.0, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    theta = 45
    majorAxis = rotateVector(axis1, axis2, theta)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_UPPER_HALF)

    return ellipse


def createBranchEllipseUp():
    centre = [-0.25, 0.0, 0.0]
    axis1 = [1.0, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    theta = -90
    majorAxis = rotateVector(axis1, axis2, theta)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_LOWER_HALF)

    return ellipse


def createBranchEllipseUpRight():
    centre = [0.25, 0.0, 0.0]
    axis1 = [1.0, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    theta = -90
    majorAxis = rotateVector(axis1, axis2, theta)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_LOWER_HALF)

    return ellipse


def createDaughterEllipseLeft():
    dc = 1.0
    # theta = 45
    centre = [-0.25, 0.0, 0.0]
    axis1 = [0.5, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    displacement = vector.setMagnitude(rotateVector(axis1, axis2, -157.5), dc)
    centre = vector.addVectors(centre, displacement)
    majorAxis = rotateVector(axis1, axis2, -50.0)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 4
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                         ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)
    return ellipse


def createDaughterEllipseRight():
    dc = 1.0
    # theta = 45
    centre = [0.25, 0.0, 0.0]
    axis1 = [-0.5, 0.0, 0.0]
    axis2 = [0.0, 1.0, 0.0]
    displacement = vector.setMagnitude(rotateVector(axis1, axis2, 150), dc)
    centre = vector.addVectors(centre, displacement)
    majorAxis = rotateVector(axis1, axis2, 22.5)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 4
    elementsCountAcrossMinor = 4
    ellipse = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                         ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)
    return ellipse


def radiusChange(radius, radiusEndRatio,elementsCountAlong,geometric=True):
    '''
    returns common ratio ro common difference for radius change.
    :param radius: cylinder base radius
    :param geometric: if True the radius change as r_n+1=r_n*ratio otherwise r_n+1 = r_n + ratio
    :return: common ratio (difference) and type of progression (either geometric or arithmetic).
    '''
    if geometric:
        ratio = math.pow(radiusEndRatio, 1.0 / elementsCountAlong)
        progression = ConeBaseProgression.GEOMETRIC_PROGRESSION
    else:
        ratio = (radiusEndRatio * radius - radius) / elementsCountAlong
        progression = ConeBaseProgression.ARITHMETIC_PROGRESSION

    return ratio, progression


def createDaughterCylinder(nodes, axis1, axis2, axis3, theta1, theta2, theta3, theta32):
    # theta = 45
    centre = [0.0, 0.0, 0.0]
    majorAxis = rotateVector(axis1, axis2, theta1)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse1p1 = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_LOWER_HALF)

    # create part 2 of first ellipse
    centre = [0.0, 0.0, 0.0]
    majorAxis = rotateVector(axis1, axis2, theta2)
    minorAxis = [0.0, 1.0, 0.0]
    elementsCountAcrossMajor = 2
    elementsCountAcrossMinor = 4
    ellipse1p2 = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                           ellipseShape=EllipseShape.Ellipse_SHAPE_UPPER_HALF)

    # Create another the daughter1 ellipse.
    dc = 2.0
    # theta = 45
    displacement = vector.setMagnitude(rotateVector(axis1, axis2, theta32), dc)
    centre = vector.addVectors(centre, displacement)
    majorAxis = rotateVector(axis1, axis2, theta3)
    minorAxis = minorAxis
    elementsCountAcrossMajor = 4
    elementsCountAcrossMinor = 4
    ellipse2 = Ellipse2D(centre, majorAxis, minorAxis, elementsCountAcrossMajor, elementsCountAcrossMinor,
                         ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

    # make a shield using three ellipses. Fill in px and pd1, pd3.
    elementsCountAlong = 1
    shield = copyEllipsesNodesToShieldNodes(elementsCountAlong, ellipse1p1, ellipse1p2, ellipse2)

    # now pd2 derivatives.
    for n3 in range(shield.elementsCountAlong + 1):
        calculateD2Derivatives(shield, n3, 1)
    smoothd2Derivatives(shield)

    skipNodes = createNodeId(shield, nodes)

    return shield, skipNodes


def rotateVector(v, axis, theta):
    a3 = vector.normalise(vector.crossproduct3(axis, v))
    v = vector.addVectors(v, a3, math.cos(math.radians(theta)), math.sin(math.radians(theta)))
    return v


def copyEllipsesNodesToShieldNodes(elementsCountAlong, ellipse1p1, ellipse1p2, ellipse2):
    """
    Copy coordinates and derivatives of ellipses to shield.
    :param n3: the index number of ellipse along the central path.
    """
    elementsCountRim = 0
    shieldMode = ShieldShape.SHIELD_SHAPE_FULL
    shield = ShieldMesh(ellipse2.elementsCountAcrossMinor, ellipse2.elementsCountAcrossMajor, elementsCountRim,
                        None, elementsCountAlong, shieldMode,
                        shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

    shield.px[1] = ellipse2.px
    shield.pd1[1] = ellipse2.pd1
    shield.pd3[1] = ellipse2.pd3

    for n2 in range(shield.elementsCountUp):
        shield.px[0][n2] = ellipse1p1.px[n2]
        shield.pd1[0][n2] = ellipse1p1.pd1[n2]
        shield.pd3[0][n2] = ellipse1p1.pd3[n2]
    for n2 in range(shield.elementsCountUp, shield.elementsCountUpFull+1):
        shield.px[0][n2] = ellipse1p2.px[shield.elementsCountUpFull-n2]
        shield.pd1[0][n2] = ellipse1p2.pd1[shield.elementsCountUpFull-n2]
        shield.pd3[0][n2] = ellipse1p2.pd3[shield.elementsCountUpFull-n2]

    return shield


def calculateD2Derivatives(shield, n3, n3Count):
    """
    calculate d2 derivatives.
    :param n3: Index of along cylinder axis coordinates to use
    :param n3Count: number of bases to create coordinates for.
    """
    btx = shield.px
    btd1 = shield.pd1
    btd2 = shield.pd2
    btd3 = shield.pd3
    # get ellipse d2 and next ellipse x, d2
    for n2 in range(shield.elementsCountUpFull + 1):
        for n1 in range(shield.elementsCountAcross + 1):
            if btd1[n3][n2][n1]:
                n3n = n3 if (n3 < n3Count) else n3 - 1
                btd2[n3][n2][n1] = [(btx[n3n + 1][n2][n1][c] - btx[n3n][n2][n1][c]) for c in range(3)]


def smoothd2Derivatives(shield):
    """
    smooth d2 derivatives using initial values calculated by calculateD2Derivatives
    """
    btx = shield.px
    btd1 = shield.pd1
    btd2 = shield.pd2
    btd3 = shield.pd3
    for n2 in range(shield.elementsCountUpFull + 1):
        for n1 in range(shield.elementsCountAcross + 1):
            td2 = []
            tx = []
            if btx[0][n2][n1]:
                for n3 in range(shield.elementsCountAlong + 1):
                    tx.append(btx[n3][n2][n1])
                    td2.append(btd2[n3][n2][n1])
                td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection=True)
                for n3 in range(shield.elementsCountAlong + 1):
                    btd2[n3][n2][n1] = td2[n3]


def createNodeId(shield, nodes):
    nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
    for n2 in range(shield.elementsCountUpFull + 1):
        for n3 in range(shield.elementsCountAlong + 1):
            for n1 in range(shield.elementsCountAcross + 1):
                if shield.px[n3][n2][n1]:
                    shield.nodeId[n3][n2][n1] = nodeIdentifier
                    nodeIdentifier += 1

    skipNodes = modifyNodeId(shield, shield, [], [], [])
    return skipNodes


def setNodeIdForInBetween(shieldInB,shieldL,shieldR, nodes):
    for n2 in range(shieldInB.elementsCountUpFull + 1):
        for n1 in range(shieldInB.elementsCountAcross + 1):
            if shieldL.px[0][n2][n1]:
                shieldInB.nodeId[0][n2][n1] = shieldL.nodeId[0][n2][n1]
            if shieldR.px[0][n2][n1]:
                shieldInB.nodeId[1][n2][n1] = shieldR.nodeId[0][n2][n1]

    skipNodes = modifyNodeId(shieldL, shieldL, [], [], [])
    return skipNodes


def modifyNodeId(shield1, shield2, n3range, n2range, n1range):
    skipNodes = [ [] for _ in range(shield2.elementsCountAlong+1) ]
    for n3 in range(shield2.elementsCountAlong + 1):
        for n2 in range(shield2.elementsCountUpFull + 1):
            for p in [skipNodes[n3]]:
                p.append([None] * (shield2.elementsCountAcross + 1))
    for n2 in n2range:
        for n3 in n3range:
            for n1 in n1range:
                if shield2.px[n3][n2][n1]:
                    shield2.nodeId[n3][n2][n1] = shield1.nodeId[n3][n2][n1]
                    skipNodes[n3][n2][n1] = 1
    return skipNodes


def generateNodes(shield, nodes, fieldmodule, coordinates, skipNodes):
    """
    Create cylinder nodes from coordinates.
    :param nodes: nodes from coordinates.
    :param fieldmodule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
    :param coordinates: Coordinate field to define.
    """

    nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
    # nodeIdentifier = shield.generateNodes(fieldModule, coordinates, nodeIdentifier)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    cache = fieldmodule.createFieldcache()

    for n2 in range(shield.elementsCountUpFull + 1):
        for n3 in range(shield.elementsCountAlong + 1):
            for n1 in range(shield.elementsCountAcross + 1):
                if shield.px[n3][n2][n1]:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    if not skipNodes[n3][n2][n1]:
                        # shield.nodeId[n3][n2][n1] = nodeIdentifier
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, shield.px[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, shield.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, shield.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, shield.pd3[n3][n2][n1])
                        nodeIdentifier += 1

    return nodeIdentifier


def generateElements(shield, mesh, fieldmodule, coordinates):
    """
    Create cylinder elements from nodes.
    :param mesh:
    :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
    :param coordinates: Coordinate field to define.
    """
    elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)


    """
      Create shield elements from nodes.
      :param fieldmodule: Zinc fieldmodule to create elements in.
      :param coordinates: Coordinate field to define.
      :param startElementIdentifier: First element identifier to use.
      :param meshGroups: Zinc mesh groups to add elements to.
      :return: next elementIdentifier.
       """
    useCrossDerivatives = False
    mesh = fieldmodule.findMeshByDimension(3)

    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    eft = tricubichermite.createEftNoCrossDerivatives()
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    isEven = (shield.elementsCountAcross % 2) == 0
    e1a = shield.elementsCountRim
    e1b = e1a + 1
    e1z = shield.elementsCountAcross - 1 - shield.elementsCountRim
    e1y = e1z - 1
    e2a = shield.elementsCountRim
    e2b = shield.elementsCountRim + 1
    e2c = shield.elementsCountRim + 2
    e2d = 2 * shield.elementsCountUp - 1

    for e3 in range(shield.elementsCountAlong):
        for e2 in range(shield.elementsCountUpFull):
            for e1 in range(shield.elementsCountAcross):
                eft1 = eft
                scalefactors = None
                if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                    nids = [shield.nodeId[e3][e2][e1], shield.nodeId[e3][e2 + 1][e1], shield.nodeId[e3 + 1][e2][e1],
                            shield.nodeId[e3 + 1][e2 + 1][e1],
                            shield.nodeId[e3][e2][e1 + 1], shield.nodeId[e3][e2 + 1][e1 + 1], shield.nodeId[e3 + 1][e2][e1 + 1],
                            shield.nodeId[e3 + 1][e2 + 1][e1 + 1]]
                elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                    nids = [shield.nodeId[0][e2][e1], shield.nodeId[0][e2][e1 + 1], shield.nodeId[0][e2 + 1][e1],
                            shield.nodeId[0][e2 + 1][e1 + 1],
                            shield.nodeId[1][e2][e1], shield.nodeId[1][e2][e1 + 1], shield.nodeId[1][e2 + 1][e1],
                            shield.nodeId[1][e2 + 1][e1 + 1]]
                if (e2 < e2b) or (e2 == e2d):
                    if (e1 < e1b) or (e1 > e1y):
                        continue  # no element due to triple point closure
                    if (e2 == e2a) or (e2 == e2d):
                        # bottom and top row elements
                        if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                            elif e2 == e2d:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                        elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                            if (e1 == e1b) or (e1 == e1y):
                                # map bottom triple point element
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e1 == e1b:
                                    remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])

                elif (e2 == e2b) or (e2 == e2d - e2b):
                    if (e1 <= e1a) or (e1 >= e1z):
                        # map top 2 triple point elements
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e1 < e1a:
                            e2r = e1
                            nids[0] = shield.nodeId[0][e2r][e1b]
                            nids[1] = shield.nodeId[0][e2r + 1][e1b]
                            nids[4] = shield.nodeId[1][e2r][e1b]
                            nids[5] = shield.nodeId[1][e2r + 1][e1b]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1a:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[0] = shield.nodeId[e3][e2a][e1b]
                                    nids[2] = shield.nodeId[e3 + 1][e2a][e1b]
                                    tripleN = [5, 7]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[1] = shield.nodeId[e3][e2d + 1][e1b]
                                    nids[3] = shield.nodeId[e3 + 1][e2d + 1][e1b]
                                    tripleN = [6, 8]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[0] = shield.nodeId[0][e2a][e1b]
                                nids[4] = shield.nodeId[1][e2a][e1b]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1z:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[4] = shield.nodeId[e3][e2a][e1z]
                                    nids[6] = shield.nodeId[e3 + 1][e2a][e1z]
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[5] = shield.nodeId[e3][e2d + 1][e1z]
                                    nids[7] = shield.nodeId[e3 + 1][e2d + 1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[1] = shield.nodeId[0][e2a][e1z]
                                nids[5] = shield.nodeId[1][e2a][e1z]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        elif e1 > e1z:
                            e2r = shield.elementsCountAcross - e1
                            nids[0] = shield.nodeId[0][e2r][e1z]
                            nids[1] = shield.nodeId[0][e2r - 1][e1z]
                            nids[4] = shield.nodeId[1][e2r][e1z]
                            nids[5] = shield.nodeId[1][e2r - 1][e1z]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                else:
                    if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                        if (e1 <= e1a):
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])

                if eft1 is not eft:
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                else:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                # print('create element shield', elementIdentifier, result2, result3, nids)
                shield.elementId[e2][e1] = elementIdentifier
                elementIdentifier += 1

                # for meshGroup in meshGroups:
                #     meshGroup.addElement(element)

    return elementIdentifier


def generateElementsInBetween(shield, mesh, fieldmodule, coordinates):
    """
    Create cylinder elements from nodes.
    :param mesh:
    :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
    :param coordinates: Coordinate field to define.
    """
    elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)


    """
      Create shield elements from nodes.
      :param fieldmodule: Zinc fieldmodule to create elements in.
      :param coordinates: Coordinate field to define.
      :param startElementIdentifier: First element identifier to use.
      :param meshGroups: Zinc mesh groups to add elements to.
      :return: next elementIdentifier.
       """
    useCrossDerivatives = False
    mesh = fieldmodule.findMeshByDimension(3)

    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    eft = tricubichermite.createEftNoCrossDerivatives()
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    isEven = (shield.elementsCountAcross % 2) == 0
    e1a = shield.elementsCountRim
    e1b = e1a + 1
    e1z = shield.elementsCountAcross - 1 - shield.elementsCountRim
    e1y = e1z - 1
    e2a = shield.elementsCountRim
    e2b = shield.elementsCountRim + 1
    e2c = shield.elementsCountRim + 2
    e2d = 2 * shield.elementsCountUp - 1

    for e3 in range(shield.elementsCountAlong):
        for e2 in range(shield.elementsCountUpFull):
            for e1 in range(shield.elementsCountAcross):
                eft1 = eft
                scalefactors = None
                if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                    nids = [shield.nodeId[e3][e2][e1], shield.nodeId[e3][e2 + 1][e1], shield.nodeId[e3 + 1][e2][e1],
                            shield.nodeId[e3 + 1][e2 + 1][e1],
                            shield.nodeId[e3][e2][e1 + 1], shield.nodeId[e3][e2 + 1][e1 + 1], shield.nodeId[e3 + 1][e2][e1 + 1],
                            shield.nodeId[e3 + 1][e2 + 1][e1 + 1]]
                elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                    nids = [shield.nodeId[0][e2][e1], shield.nodeId[0][e2][e1 + 1], shield.nodeId[0][e2 + 1][e1],
                            shield.nodeId[0][e2 + 1][e1 + 1],
                            shield.nodeId[1][e2][e1], shield.nodeId[1][e2][e1 + 1], shield.nodeId[1][e2 + 1][e1],
                            shield.nodeId[1][e2 + 1][e1 + 1]]
                if (e2 < e2b) or (e2 == e2d):
                    if (e1 < e1b) or (e1 > e1y):
                        continue  # no element due to triple point closure
                    if (e2 == e2a) or (e2 == e2d):
                        # bottom and top row elements
                        if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                            elif e2 == e2d:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                        elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                            if (e1 == e1b) or (e1 == e1y):
                                # map bottom triple point element
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e1 == e1b:
                                    remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])

                elif (e2 == e2b) or (e2 == e2d - e2b):
                    if (e1 <= e1a) or (e1 >= e1z):
                        # map top 2 triple point elements
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e1 < e1a:
                            e2r = e1
                            nids[0] = shield.nodeId[0][e2r][e1b]
                            nids[1] = shield.nodeId[0][e2r + 1][e1b]
                            nids[4] = shield.nodeId[1][e2r][e1b]
                            nids[5] = shield.nodeId[1][e2r + 1][e1b]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1a:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[0] = shield.nodeId[e3][e2a][e1b]
                                    nids[2] = shield.nodeId[e3 + 1][e2a][e1b]
                                    tripleN = [5, 7]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[1] = shield.nodeId[e3][e2d + 1][e1b]
                                    nids[3] = shield.nodeId[e3 + 1][e2d + 1][e1b]
                                    tripleN = [6, 8]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[0] = shield.nodeId[0][e2a][e1b]
                                nids[4] = shield.nodeId[1][e2a][e1b]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1z:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[4] = shield.nodeId[e3][e2a][e1z]
                                    nids[6] = shield.nodeId[e3 + 1][e2a][e1z]
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[5] = shield.nodeId[e3][e2d + 1][e1z]
                                    nids[7] = shield.nodeId[e3 + 1][e2d + 1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[1] = shield.nodeId[0][e2a][e1z]
                                nids[5] = shield.nodeId[1][e2a][e1z]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        elif e1 > e1z:
                            e2r = shield.elementsCountAcross - e1
                            nids[0] = shield.nodeId[0][e2r][e1z]
                            nids[1] = shield.nodeId[0][e2r - 1][e1z]
                            nids[4] = shield.nodeId[1][e2r][e1z]
                            nids[5] = shield.nodeId[1][e2r - 1][e1z]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                else:
                    if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                        if (e1 <= e1a):
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])

                # if e3 == 0:
                #
                #     eft1 = tricubichermite.createEftNoCrossDerivatives()
                #     setEftScaleFactorIds(eft1, [1], [])
                #     scalefactors = [-1.0]
                #     remapEftNodeValueLabel(eft1, [1,2,5,6], Node.VALUE_LABEL_D_DS2,
                #                            [(Node.VALUE_LABEL_D_DS2, [1])])
                    # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                    #                        [(Node.VALUE_LABEL_D2_DS1DS2, [])])

                if eft1 is not eft:
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                else:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                # print('create element shield', elementIdentifier, result2, result3, nids)
                shield.elementId[e2][e1] = elementIdentifier
                elementIdentifier += 1

                # for meshGroup in meshGroups:
                #     meshGroup.addElement(element)

    return elementIdentifier


def generateElementsBaseRight(shield, mesh, fieldmodule, coordinates):
    """
    Create cylinder elements from nodes.
    :param mesh:
    :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
    :param coordinates: Coordinate field to define.
    """
    elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)


    """
      Create shield elements from nodes.
      :param fieldmodule: Zinc fieldmodule to create elements in.
      :param coordinates: Coordinate field to define.
      :param startElementIdentifier: First element identifier to use.
      :param meshGroups: Zinc mesh groups to add elements to.
      :return: next elementIdentifier.
       """
    useCrossDerivatives = False
    mesh = fieldmodule.findMeshByDimension(3)

    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    eft = tricubichermite.createEftNoCrossDerivatives()
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    isEven = (shield.elementsCountAcross % 2) == 0
    e1a = shield.elementsCountRim
    e1b = e1a + 1
    e1z = shield.elementsCountAcross - 1 - shield.elementsCountRim
    e1y = e1z - 1
    e2a = shield.elementsCountRim
    e2b = shield.elementsCountRim + 1
    e2c = shield.elementsCountRim + 2
    e2d = 2 * shield.elementsCountUp - 1

    for e3 in range(shield.elementsCountAlong):
        for e2 in range(shield.elementsCountUpFull):
            for e1 in range(shield.elementsCountAcross):
                eft1 = eft
                scalefactors = None
                if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                    nids = [shield.nodeId[e3][e2][e1], shield.nodeId[e3][e2 + 1][e1], shield.nodeId[e3 + 1][e2][e1],
                            shield.nodeId[e3 + 1][e2 + 1][e1],
                            shield.nodeId[e3][e2][e1 + 1], shield.nodeId[e3][e2 + 1][e1 + 1], shield.nodeId[e3 + 1][e2][e1 + 1],
                            shield.nodeId[e3 + 1][e2 + 1][e1 + 1]]
                elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                    nids = [shield.nodeId[0][e2][e1], shield.nodeId[0][e2][e1 + 1], shield.nodeId[0][e2 + 1][e1],
                            shield.nodeId[0][e2 + 1][e1 + 1],
                            shield.nodeId[1][e2][e1], shield.nodeId[1][e2][e1 + 1], shield.nodeId[1][e2 + 1][e1],
                            shield.nodeId[1][e2 + 1][e1 + 1]]
                if (e2 < e2b) or (e2 == e2d):
                    if (e1 < e1b) or (e1 > e1y):
                        continue  # no element due to triple point closure
                    if (e2 == e2a) or (e2 == e2d):
                        # bottom and top row elements
                        if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                            elif e2 == e2d:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                        elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                            if (e1 == e1b) or (e1 == e1y):
                                # map bottom triple point element
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e1 == e1b:
                                    remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])

                elif (e2 == e2b) or (e2 == e2d - e2b):
                    if (e1 <= e1a) or (e1 >= e1z):
                        # map top 2 triple point elements
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e1 < e1a:
                            e2r = e1
                            nids[0] = shield.nodeId[0][e2r][e1b]
                            nids[1] = shield.nodeId[0][e2r + 1][e1b]
                            nids[4] = shield.nodeId[1][e2r][e1b]
                            nids[5] = shield.nodeId[1][e2r + 1][e1b]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1a:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[0] = shield.nodeId[e3][e2a][e1b]
                                    nids[2] = shield.nodeId[e3 + 1][e2a][e1b]
                                    tripleN = [5, 7]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[1] = shield.nodeId[e3][e2d + 1][e1b]
                                    nids[3] = shield.nodeId[e3 + 1][e2d + 1][e1b]
                                    tripleN = [6, 8]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[0] = shield.nodeId[0][e2a][e1b]
                                nids[4] = shield.nodeId[1][e2a][e1b]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == e1z:
                            if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                                if e2 == e2b:
                                    nids[4] = shield.nodeId[e3][e2a][e1z]
                                    nids[6] = shield.nodeId[e3 + 1][e2a][e1z]
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2d - e2b:
                                    nids[5] = shield.nodeId[e3][e2d + 1][e1z]
                                    nids[7] = shield.nodeId[e3 + 1][e2d + 1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_REGULAR:
                                nids[1] = shield.nodeId[0][e2a][e1z]
                                nids[5] = shield.nodeId[1][e2a][e1z]
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        elif e1 > e1z:
                            e2r = shield.elementsCountAcross - e1
                            nids[0] = shield.nodeId[0][e2r][e1z]
                            nids[1] = shield.nodeId[0][e2r - 1][e1z]
                            nids[4] = shield.nodeId[1][e2r][e1z]
                            nids[5] = shield.nodeId[1][e2r - 1][e1z]
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                else:
                    if shield._type == ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND:
                        if (e1 <= e1a):
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])

                # if e3 == 0:
                #
                #     eft1 = tricubichermite.createEftNoCrossDerivatives()
                #     setEftScaleFactorIds(eft1, [1], [])
                #     scalefactors = [-1.0]
                #     remapEftNodeValueLabel(eft1, [1,2,5,6], Node.VALUE_LABEL_D_DS2,
                #                            [(Node.VALUE_LABEL_D_DS2, [1])])
                    # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                    #                        [(Node.VALUE_LABEL_D2_DS1DS2, [])])

                if eft1 is not eft:
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                else:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                # print('create element shield', elementIdentifier, result2, result3, nids)
                shield.elementId[e2][e1] = elementIdentifier
                elementIdentifier += 1

                # for meshGroup in meshGroups:
                #     meshGroup.addElement(element)

    return elementIdentifier