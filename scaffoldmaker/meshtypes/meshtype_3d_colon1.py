"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_haustra1 import MeshType_3d_haustra1, getColonHaustraSegmentInnerPoints
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.matrix import *
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh import *
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a haustra
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 2D 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 2,
                'Length' : 1.0,
                'Number of elements' : 5
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ], [
                [ [  0.0,  0.0 ], [  0.0, 10.0 ] ],
                [ [  0.0, 10.0 ], [  5.0,  5.0 ] ],
                [ [  5.0,  9.0 ], [  5.0,  0.0 ] ],
                [ [ 10.0, 10.0 ], [  5.0, -5.0 ] ],
                [ [ 10.0, -2.0 ], [ -3.0, -5.0 ] ],
                [ [  7.0, -4.0 ], [ -3.0,  0.0 ] ] ] )
            } ),
        'Human 3D 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 6
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ], [
                [ [  0.0,  0.0,  0.0 ], [  0.0, 10.0,  3.0 ] ],
                [ [  0.0, 10.0,  3.0 ], [  5.0,  5.0,  0.0 ] ],
                [ [  5.0,  9.0,  0.0 ], [  5.0,  0.0,  0.0 ] ],
                [ [ 10.0, 10.0,  2.0 ], [ 10.0, -5.0,  0.0 ] ],
                [ [ 15.0, 15.0,  7.0 ], [ 12.0, 12.0,  0.0 ] ],
                [ [ 20.0, -2.0,  0.0 ], [  5.0,-12.0, -5.0 ] ],
                [ [ 10.0, -4.0, -0.0 ], [ -8.0,  0.0,  0.0 ] ] ] )
            } ),
        'Test Line' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ], [
                [ [ -4.0,  1.0,  3.0 ], [  5.0,  1.0, -3.0 ] ],
                [ [  1.0,  2.0,  0.0 ], [  5.0,  1.0, -3.0 ] ] ] )
            } )
        }

    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = MeshType_3d_haustra1.getDefaultOptions(parameterSetName)
        options['Number of elements around'] = 15
        options['Number of elements along haustrum'] = 4
        options['Inner radius'] = 1.0
        options['Haustrum length mid derivative factor'] = 2.0
        options['Wall thickness'] = 0.05
        if 'Human 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 3D 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2D 1']
        optionsColon = {
            'Central path' : copy.deepcopy(centralPathOption),
            'Number of haustra segments': 30
            }
        options.update(optionsColon)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_haustra1.getOrderedOptionNames()
        optionNames.remove('Haustrum length')
        for optionName in [
            'Number of haustra segments',
            'Central path']:
            optionNames.insert(0, optionName)
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        return Scaffold_base.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType)

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if options['Number of haustra segments'] < 1:
            options['Number of haustra segments'] = 1

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        centralPath = options['Central path']
        elementsCountAround = options['Number of elements around']
        elementsCountAlongHaustrum = options['Number of elements along haustrum']
        elementsCountThroughWall = options['Number of elements through wall']
        haustraSegmentCount = options['Number of haustra segments']
        radius = options['Inner radius']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        haustrumLengthEndDerivativeFactor = options['Haustrum length end derivative factor']
        haustrumLengthMidDerivativeFactor = options['Haustrum length mid derivative factor']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongHaustrum*haustraSegmentCount)

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections = True,
            magnitudeScalingMode = DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            length += arcLength
        haustrumLength = length / haustraSegmentCount

        # Generate inner surface of a haustra segment
        xHaustraInner, d1HaustraInner, d2HaustraInner, haustraSegmentAxis = getColonHaustraSegmentInnerPoints(elementsCountAround, elementsCountAlongHaustrum, radius, cornerInnerRadiusFactor,
            haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier = generatetubemesh(region, elementsCountAround, elementsCountAlongHaustrum, elementsCountThroughWall, haustraSegmentCount,
            cx, cd1, xHaustraInner, d1HaustraInner, d2HaustraInner, wallThickness, haustraSegmentAxis, haustrumLength, useCrossDerivatives, useCubicHermiteThroughWall)

        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along haustrum']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()
