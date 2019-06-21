"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_haustra1 import MeshType_3d_haustra1, getColonHaustraSegmentInnerPoints, getTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.node import Node

class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a haustra
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [  0.0,  0.0,  0.0 ], [ -1.6,  5.8,  0.0 ], [ -2.4, -0.6, -1.2 ], [ -1.4, -0.1, -1.2 ] ],
                [ [ -1.6,  6.1,  0.0 ], [ -0.8,  6.2,  0.0 ], [ -2.2, -0.4, -0.8 ], [ -0.4,  1.9,  2.2 ] ],
                [ [ -0.3, 12.3,  0.0 ], [  4.8,  1.2,  0.0 ], [ -1.0,  2.0,  0.8 ], [ -0.6,  0.0,  5.1 ] ],
                [ [  3.0, 12.2,  0.0 ], [  4.0, -0.3,  0.0 ], [ -0.5,  0.4,  2.9 ], [  0.0,  0.1,  2.4 ] ],
                [ [  7.3, 12.3,  0.0 ], [  4.6,  1.1,  0.0 ], [ -0.2,  1.0,  2.2 ], [  0.5,  2.5, -2.0 ] ],
                [ [ 12.2, 12.9,  0.0 ], [  5.1, -3.6,  0.0 ], [  1.0,  1.7,  0.6 ], [  0.1, -0.6, -3.5 ] ],
                [ [ 13.9,  6.4,  0.0 ], [  0.2, -4.0,  0.0 ], [  2.0,  0.0, -2.0 ], [  1.5, -0.1, -1.0 ] ],
                [ [ 12.9, -0.3,  0.0 ], [ -4.1, -3.0,  0.0 ], [  0.6, -0.9, -1.4 ], [  0.8, -1.1, -1.3 ] ],
                [ [  7.7,  0.3,  0.0 ], [ -3.0,  1.7,  0.0 ], [  0.1, -1.1, -1.8 ], [  0.4, -1.2, -1.2 ] ] ] )
            } ),
        'Human 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  0.0,  0.0,  0.0 ], [ -1.6,  5.8,  0.0 ], [ -2.4, -0.6, -1.2 ], [ -1.4, -0.1, -1.2 ] ],
                [ [ -1.6,  6.1,  0.0 ], [ -0.8,  6.2,  0.0 ], [ -2.2, -0.4, -0.8 ], [ -0.4,  1.9,  2.2 ] ],
                [ [ -0.3, 12.3,  0.0 ], [  3.8, -0.4,  2.6 ], [ -1.0,  2.0,  0.8 ], [ -0.6,  0.0,  5.1 ] ],
                [ [  3.6, 10.4,  2.3 ], [  4.4,  3.9,  1.4 ], [ -0.5,  0.4,  2.9 ], [  0.0,  0.1,  2.4 ] ],
                [ [  9.5, 11.7,  2.4 ], [  5.9, -0.6,  2.0 ], [ -0.2,  1.0,  2.2 ], [  0.5,  2.5, -2.0 ] ],
                [ [ 12.9, 17.3, -1.1 ], [  1.8, -1.5, -5.5 ], [  1.0,  1.7,  0.6 ], [  0.1, -0.6, -3.5 ] ],
                [ [ 14.6,  7.8, -1.5 ], [ -0.2, -3.8, -0.1 ], [  2.0,  0.0, -2.0 ], [  1.5, -0.1, -1.0 ] ],
                [ [ 12.9, -0.3,  0.0 ], [ -4.1, -3.0,  0.0 ], [  0.6, -0.9, -1.4 ], [  0.8, -1.1, -1.3 ] ],
                [ [  7.7,  0.3,  0.0 ], [ -3.0,  1.7,  0.0 ], [  0.1, -1.1, -1.8 ], [  0.4, -1.2, -1.2 ] ] ] )
            } ),
        'Test Line' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [ -4.0, 1.0, 3.0 ], [ 5.0, 1.0, -3.0 ], [ 0.8, -4.0, 0.0 ], [0.0, 0.0, 0.0 ] ],
                [ [  1.0, 2.0, 0.0 ], [ 5.0, 1.0, -3.0 ], [ 0.8, -4.0, 0.0 ], [0.0, 0.0, 0.0 ] ] ])
            } )
        }

    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Human 2']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = MeshType_3d_haustra1.getDefaultOptions(parameterSetName)
        options['Number of elements along haustrum'] = 4
        options['Inner radius'] = 1.0
        options['Haustrum length mid derivative factor'] = 2.0
        options['Wall thickness'] = 0.02
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
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
        MeshType_3d_haustra1.checkOptions(options)
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
        elementsCountAroundTC = options['Number of elements around tenia coli']
        elementsCountAroundHaustrum = options['Number of elements around haustrum']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*3
        elementsCountAlongHaustrum = options['Number of elements along haustrum']
        elementsCountThroughWall = options['Number of elements through wall']
        haustraSegmentCount = options['Number of haustra segments']
        radius = options['Inner radius']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        haustrumLengthEndDerivativeFactor = options['Haustrum length end derivative factor']
        haustrumLengthMidDerivativeFactor = options['Haustrum length mid derivative factor']
        widthTC = options['Tenia coli width']
        TCThickness = options['Tenia coli thickness']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongHaustrum*haustraSegmentCount)

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections = True,
            magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            length += arcLength
        haustrumLength = length / haustraSegmentCount

        # Generate inner surface of a haustra segment
        xHaustraInner, d1HaustraInner, d2HaustraInner, haustraSegmentAxis = getColonHaustraSegmentInnerPoints(elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
            haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier, xList, d1List, d2List, d3List, sx, curvatureAlong, factorList = tubemesh.generatetubemesh(region, elementsCountAround, elementsCountAlongHaustrum, elementsCountThroughWall, haustraSegmentCount,
            cx, cd1, cd2, cd12, xHaustraInner, d1HaustraInner, d2HaustraInner, wallThickness, haustraSegmentAxis, haustrumLength, useCrossDerivatives, useCubicHermiteThroughWall)

        # Generate tenia coli
        annotationGroupsTC, nextNodeIdentifier, nextElementIdentifier = getTeniaColi(region, nextNodeIdentifier, nextElementIdentifier, useCrossDerivatives, useCubicHermiteThroughWall,
            xList, d1List, d2List, d3List, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall, widthTC, TCThickness, sx, curvatureAlong, factorList)

        annotationGroups += annotationGroupsTC

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
            return cls.generateBaseMesh(region, options)

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along haustrum']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()
