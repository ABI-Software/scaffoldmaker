"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

from scaffoldmaker.meshtypes.meshtype_3d_haustra1 import MeshType_3d_haustra1, getColonHaustraSegmentInnerPoints
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.matrix import *
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh import *

class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a haustra
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''
    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = MeshType_3d_haustra1.getDefaultOptions(parameterSetName)
        options['Number of elements around'] = 15
        options['Number of elements along haustrum'] = 4
        options['Inner radius'] = 1.0
        options['Haustrum length mid derivative factor'] = 2.0
        options['Wall thickness'] = 0.05
        optionsColon = {
            'Number of haustra segments': 30,
            'Tube type': 2
            }
        options.update(optionsColon)
        if 'Human 1' in parameterSetName:
            options['Tube type'] = 3
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_haustra1.getOrderedOptionNames()
        optionNames.remove('Haustrum length')
        for optionName in [
            'Number of haustra segments',
            'Tube type']:
            optionNames.insert(3, optionName)
        return optionNames

    def checkOptions(options):
        MeshType_3d_haustra1.checkOptions(options)
        if options['Number of haustra segments'] < 1:
            options['Number of haustra segments'] = 1
        if options['Tube type'] < 1:
            options['Tube type'] = 1
        if options['Tube type'] > 3:
            options['Tube type'] = 3

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
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
        tubeType = options['Tube type']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongHaustrum*haustraSegmentCount)

        if tubeType == 1: # Straight tube
            cx = [[-4.0, 1.0, 3.0], [ 1.0, 2.0, 0.0 ] ]
            cd1 = [[ 5.0, 1.0, -3.0 ], [ 5.0, 1.0, -3.0 ]]
        elif tubeType == 2: # Human colon in x-y plane
            cx = [ [ 0.0, 0.0, 0.0], [0.0, 10.0, 0.0], [5.0, 9.0, 0.0], [ 10.0, 10.0, 0.0 ], [ 10.0, -2.0, 0.0], [ 7.0, -4.0, 0.0] ]
            cd1 = [ [ 0.0, 10.0, 0.0 ], [ 5.0, 5.0, 0.0 ], [5.0, 0.0, 0.0], [ 5.0, -5.0, 0.0 ], [ -3.0, -5.0, 0.0 ], [ -3.0, 0.0, 0.0 ]]
        elif tubeType == 3: # Human colon in 3D
            cx = [ [ 0.0, 0.0, 0.0], [0.0, 10.0, 3.0], [5.0, 9.0, 0.0], [ 10.0, 10.0, 2.0 ], [15.0, 15.0, 7.0], [ 20.0, -2.0, 0.0], [ 10.0, -4.0, -0.0] ]
            cd1 = [ [ 0.0, 10.0, 3.0 ], [ 5.0, 5.0, 0.0 ], [5.0, 0.0, 0.0], [ 10.0, -5.0, 0.0 ], [12.0, 12.0, 0.0], [ 5.0, -12.0, -5.0 ], [ -8.0, 0.0, 0.0 ]]

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
