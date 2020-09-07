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
    CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite,\
    smoothCubicHermiteDerivativesLine, getCubicHermiteArcLength,DerivativeScalingMode, sampleParameterAlongLine
from scaffoldmaker.utils import centralpath

class MeshType_3d_solidcylinder1(Scaffold_base):
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
        return '3D Solid Cylinder 1'


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

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]


        if not useCentralPath:
            majorRatio, majorProgression = radiusChange(majorRadius, majorRadiusEndRatio, elementsCountAlong, geometric=majorGeometric)
            minorRatio, minorProgression = radiusChange(minorRadius, minorRadiusEndRatio, elementsCountAlong, geometric=minorGeometric)
            radiusChanges = Tapered(majorRatio, majorProgression, minorRatio, minorProgression)
        else:
            radiusChanges = []

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, [0.0, 0.0, 0.0],
                            vector.setMagnitude(axis3, length), vector.setMagnitude(axis1, majorRadius), minorRadius)
        cylinder1 = CylinderMesh(fm, coordinates, base, elementsCountAlong,
                                 cylinderShape=cylinderShape, tapered=radiusChanges,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

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

def radiusChange(radius,radiusEndRatio,elementsCountAlong,geometric=True):
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
