"""
Generates a 3-D cecum mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
import math
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import ColonSegmentTubeMeshInnerPoints, getTeniaColi, createFlatAndTextureCoordinatesTeniaColi, createNodesAndElementsTeniaColi, createHalfSetInterHaustralSegment, createHalfSetIntraHaustralSegment, getFullProfileFromHalfHaustrum, getXiListFromOuterLengthProfile
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates # KM
from opencmiss.zinc.field import Field #KM

class MeshType_3d_cecum1(Scaffold_base):
    '''
    Generates a 3-D cecum mesh with variable numbers
    of elements around, along the central line, and through wall.
    The cecum is created by a function that generates a cecum
    segment and uses tubemesh to map the segment along a central
    line profile. The proximal end of the cecum is closed up with
    an apex plate. An ostium is included to generate the
    ileo-ceco-colic junction.
    '''

    centralPathDefaultScaffoldPackages = {
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 120.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    # [[0.0, 0.0, 0.0], [40.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    # [[40.0, 0.0, 0.0], [40.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    # [[80.0, 0.0, 0.0], [40.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    # [[120.0, 0.0, 0.0], [40.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]])
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 40.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 40.0], [0.0, 0.0, 40.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 80.0], [0.0, 0.0, 40.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 120.0], [0.0, 0.0, 40.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]])

    } ),
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [   0.0,   0.0, 0.0 ], [ -50.7,  178.2, 0.0 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ] ],
                [ [ -47.4, 188.6, 0.0 ], [ -19.3,  177.1, 0.0 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ] ],
                [ [  -4.4, 396.5, 0.0 ], [ 206.0,   40.1, 0.0 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ] ],
                [ [ 130.0, 384.1, 0.0 ], [ 130.8,  -40.5, 0.0 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ] ],
                [ [ 279.4, 383.0, 0.0 ], [ 118.0,   48.7, 0.0 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ] ],
                [ [ 443.9, 390.8, 0.0 ], [ 111.3,  -97.0, 0.0 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ] ],
                [ [ 475.2, 168.0, 0.0 ], [  -0.8, -112.4, 0.0 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ] ],
                [ [ 432.6, -32.3, 0.0 ], [ -90.5,  -59.0, 0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ] ],
                [ [ 272.4,   7.5, 0.0 ], [ -79.0,   47.4, 0.0 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ] ] ] )
            } ),
        }

    @staticmethod
    def getName():
        return '3D Cecum 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']

        # if 'Human 2' in parameterSetName:
        #     centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        # elif 'Mouse 1' in parameterSetName:
        #     centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        # elif 'Mouse 2' in parameterSetName:
        #     centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 2']
        # elif 'Pig 1' in parameterSetName:
        #     centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        # else:
        #     centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        # if 'Mouse' in parameterSetName:
        #     segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Mouse 1')
        # elif 'Pig' in parameterSetName:
        #     segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Pig 1')
        # else:
        #     segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Human 1')

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of segments': 4,
            'Number of elements around tenia coli': 2,
            'Number of elements around haustrum': 8,
            'Number of elements along segment': 8,
            'Number of elements through wall': 1,
            'Start inner radius': 20.0,
            'Start inner radius derivative': 0.0,
            'End inner radius': 20.0,
            'End inner radius derivative': 0.0,
            'Corner inner radius factor': 0.5,
            'Haustrum inner radius factor': 0.6,
            'Segment length end derivative factor': 0.5,
            'Segment length mid derivative factor': 4.0,
            'Number of tenia coli': 3,
            'Start tenia coli width': 5.0,
            'Start tenia coli width derivative': 0.0,
            'End tenia coli width': 5.0,
            'End tenia coli width derivative': 0.0,
            'Tenia coli thickness': 0.5,
            'Wall thickness': 0.5, #2.0,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along segment': 1,
            'Refine number of elements through wall': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of segments',
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Start inner radius',
            'Start inner radius derivative',
            'End inner radius',
            'End inner radius derivative',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        # if optionName == 'Segment profile':
        #     return [ MeshType_3d_colonsegment1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
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
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        # if optionName == 'Segment profile':
        #     if not parameterSetName:
        #         parameterSetName = scaffoldType.getParameterSetNames()[0]
        #     return ScaffoldPackage(scaffoldType, defaultParameterSetName = parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        # if not options['Segment profile'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Segment profile'):
        #     options['Segment profile'] = cls.getOptionScaffoldPackage('Segment profile', MeshType_3d_colonsegmentteniacoli1)
        for key in [
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Start inner radius',
            'Start inner radius derivative',
            'End inner radius',
            'End inner radius derivative',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        centralPath = options['Central path']
        segmentCount = options['Number of segments']
        startPhase = 0.0
        elementsCountAroundTC = options['Number of elements around tenia coli']
        elementsCountAroundHaustrum = options['Number of elements around haustrum']
        elementsCountAlongSegment = options['Number of elements along segment']
        elementsCountThroughWall = options['Number of elements through wall']
        startInnerRadius = options['Start inner radius']
        startInnerRadiusDerivative = options['Start inner radius derivative']
        endInnerRadius = options['End inner radius']
        endInnerRadiusDerivative = options['End inner radius derivative']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = options['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = options['Segment length mid derivative factor']
        tcCount = options['Number of tenia coli']
        startTCWidth = options['Start tenia coli width']
        startTCWidthDerivative = options['Start tenia coli width derivative']
        endTCWidth = options['End tenia coli width']
        endTCWidthDerivative = options['End tenia coli width derivative']
        tcThickness = options['Tenia coli thickness']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment*segmentCount)
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*tcCount

        ##################################################################################
        # zero = [0.0, 0.0, 0.0]
        # fm = region.getFieldmodule()
        # fm.beginChange()
        # cache = fm.createFieldcache()
        # nodeIdentifier = 1
        #
        # # Coordinates field
        # coordinates = findOrCreateFieldCoordinates(fm)
        # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(coordinates)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        # if useCubicHermiteThroughWall:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        #     if useCrossDerivatives:
        #         nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        #         nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        #         nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        #####################################################################################

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], '],')
        del tmpRegion

        # find arclength of cecum
        cecumLength = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            # print(e+1, arcLength)
            cecumLength += arcLength

        # Sample central path
        smoothd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        sxCecum, sd1Cecum, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, smoothd1, elementsCountAlong)
        sd2Cecum, sd12Cecum = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Calculate segment length
        segmentLength = cecumLength / segmentCount

        # Generate variation of radius & tc width along length
        innerRadiusAlongCecum = []
        dInnerRadiusAlongCecum = []
        tcWidthAlongCecum = []

        for n2 in range(elementsCountAlongSegment*segmentCount + 1):
            xi = 1/(elementsCountAlongSegment*segmentCount) * n2

            radius = interp.interpolateCubicHermite([startInnerRadius], [startInnerRadiusDerivative],
                                                    [endInnerRadius], [endInnerRadiusDerivative], xi)[0]
            innerRadiusAlongCecum.append(radius)
            dRadius = interp.interpolateCubicHermiteDerivative([startInnerRadius], [startInnerRadiusDerivative],
                                                               [endInnerRadius], [endInnerRadiusDerivative], xi)[0]
            dInnerRadiusAlongCecum.append(dRadius)
            tcWidth = interp.interpolateCubicHermite([startTCWidth], [startTCWidthDerivative],
                                                     [endTCWidth], [endTCWidthDerivative], xi)[0]
            tcWidthAlongCecum.append(tcWidth)

        xToSample = []
        d1ToSample = []
        d2ToSample = []
        zFirstNodeAlongList = []

        elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum)*0.5)

        # Create object
        colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactor,
            innerRadiusAlongCecum, dInnerRadiusAlongCecum, tcWidthAlongCecum, startPhase)

        for nSegment in range(segmentCount):
            # Make regular segments
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroups, annotationArray \
                = colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

            # Add apex to second half of first segment and sample along
            if nSegment == 0:
                xFirstSegmentSampled, d1FirstSegmentSampled, d2FirstSegmentSampled, d1FirstDirectionVector = \
                    getApexSegmentForCecum(xInner, d1Inner, d2Inner, elementsCountAroundHalfHaustrum,
                                           elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment,
                                           tcCount)

                xToSample += xFirstSegmentSampled
                d1ToSample += d1FirstSegmentSampled
                d2ToSample += d2FirstSegmentSampled
            else:
                xInnerExtrude = []
                for n in range(len(xInner)):
                    xInnerExtrude.append([xInner[n][0], xInner[n][1], xInner[n][2] + segmentLength*nSegment])
                xToSample += xInnerExtrude[elementsCountAround:]
                d1ToSample += d1Inner[elementsCountAround:]
                d2ToSample += d2Inner[elementsCountAround:]

        # Sample along length
        xToWarp, d1ToWarp, d2ToWarp = sampleCecumAlongLength(xToSample, d1ToSample, d2ToSample, d1FirstDirectionVector,
                                                             elementsCountAroundHalfHaustrum, elementsCountAroundTC,
                                                             elementsCountAround, elementsCountAlong, tcCount)

        # Project reference point for warping onto central path
        sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
            tubemesh.getPlaneProjectionOnCentralPath(xToWarp, elementsCountAround, elementsCountAlong,
                                                     cecumLength, sxCecum, sd1Cecum, sd2Cecum, sd12Cecum)

        # Warp points
        xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
            tubemesh.warpSegmentPoints(xToWarp, d1ToWarp, d2ToWarp, segmentAxis, sxRefList, sd1RefList,
                                       sd2ProjectedListRef, elementsCountAround, elementsCountAlong,
                                       zRefList, innerRadiusAlongCecum, closedProximalEnd=True)

        # xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
        #     xToWarp, d1ToWarp, d2ToWarp, segmentAxis,
        #     sxCentroidList, sd1CentroidList, sd2ProjectedListCentroid,
        #     elementsCountAround, elementsCountAlong, 0, zFirstNodeAlongList, innerRadiusAlongCecum, #closedProximalEnd=True)
        #     region, useCubicHermiteThroughWall, useCrossDerivatives, nodeIdentifier, closedProximalEnd=True)
        # ##############################################################################################################
        # Create nodes
        # Coordinates field
        # for n in range(len(xToWarped)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xToWarp[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     if useCrossDerivatives:
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        #     #print('NodeIdentifier = ', nodeIdentifier, xWarpedList[n], d1WarpedList[n], d2WarpedList[n])
        #     nodeIdentifier = nodeIdentifier + 1
        # fm.endChange()
        # ###############################################################################################################

        # Calculate unit d3
        # # xWarpedList = xToWarp
        # # d1WarpedList = d1ToWarp
        # # d2WarpedList = d2ToWarp
        d3UnitWarpedList = []
        for n in range(len(xWarpedList)):
            d3Unit = vector.normalise(
                vector.crossproduct3(vector.normalise(d1WarpedList[n]), vector.normalise(d2WarpedList[n])))
            d3UnitWarpedList.append(d3Unit)

        # Create coordinates and derivatives
        wallThicknessList = [wallThickness] * (elementsCountAlong + 1)

        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList,
            d2WarpedList, d3UnitWarpedList, wallThicknessList,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        # Deal with multiple nodes at end point for closed proximal end
        xApexInner = xList[0]
        mag = interp.getCubicHermiteArcLength(xList[0], d2List[0], xList[elementsCountAround*2], d2List[elementsCountAround*2]) # arclength between apex point and corresponding point on next face
        d2ApexInner = vector.setMagnitude(sd2Cecum[0], mag)
        d1ApexInner = vector.crossproduct3(sd1Cecum[0], d2ApexInner)
        d1ApexInner = vector.setMagnitude(d1ApexInner, mag)
        d3ApexUnit = vector.normalise(vector.crossproduct3(vector.normalise(d1ApexInner), vector.normalise(d2ApexInner)))
        d3ApexInner = [d3ApexUnit[c] * wallThickness/elementsCountThroughWall for c in range(3)]

        xApexOuter = [xApexInner[c] + d3ApexUnit[c] * wallThickness for c in range(3)]
        d1ApexOuter = d1ApexInner # Probably need to scale to curvature
        d2ApexOuter = d2ApexInner
        d3ApexOuter = d3ApexInner

        xCecum = [xApexInner] + [xApexOuter] + xList[elementsCountAround*2:]
        d1Cecum = [d1ApexInner] + [d1ApexOuter] + d1List[elementsCountAround * 2:]
        d2Cecum = [d2ApexInner] + [d2ApexOuter] + d2List[elementsCountAround * 2:]
        d3Cecum = [d3ApexInner] + [d3ApexOuter] + d3List[elementsCountAround * 2:]

        xFlat = d1Flat = d2Flat = d3Flat = []
        xTexture = d1Texture = d2Texture = d3Texture = []

        # Create nodes and elements
        nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)

        # ################################################################################################################
        # nodeIdentifier = nextNodeIdentifier
        # for n in range(len(sxCentroidList)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sxCentroidList[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1CentroidList[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [30.0* sd2ProjectedListCentroid[n][c] for c in range(3)])
        #     if useCrossDerivatives:
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        #     # print('NodeIdentifier = ', nodeIdentifier, xWarpedList[n], d1WarpedList[n], d2WarpedList[n])
        #     nodeIdentifier = nodeIdentifier + 1
        # fm.endChange()
        ###############################################################################################################

        # if tcThickness > 0:
        #     tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
        #
        #     xList, d1List, d2List, d3List, annotationGroups, annotationArray = getTeniaColi(
        #         region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
        #         elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
        #         tubeTCWidthList, tcThickness, sxCecum, annotationGroups, annotationArray, closedProximalEnd = True)
        #
        #     # Create nodes and elements
        #     nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
        #         region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
        #         elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
        #         tcCount, annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
        #         useCubicHermiteThroughWall, useCrossDerivatives)
        #
        # else:
        #     # Create nodes and elements
        #     nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
        #         region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
        #         elementsCountAround, elementsCountAlong, elementsCountThroughWall,
        #         annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
        #         useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)

        ##############################################################################################################
        # # Create nodes
        # # Coordinates field
        # for n in range(len(xSampled)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampled[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Sampled[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Sampled[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     if useCrossDerivatives:
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        #     #print('NodeIdentifier = ', nodeIdentifier, xWarpedList[n], d1WarpedList[n], d2WarpedList[n])
        #     nodeIdentifier = nodeIdentifier + 1
        # fm.endChange()
        ###############################################################################################################

        # if tcThickness > 0:
        #     tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
        #     xList, d1List, d2List, d3List, annotationGroups, annotationArray = getTeniaColi(
        #         region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
        #         elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
        #         tubeTCWidthList, tcThickness, sxCecum, annotationGroups, annotationArray)
        #
        #     # Create flat and texture coordinates
        #     xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = createFlatAndTextureCoordinatesTeniaColi(
        #         xiList, relaxedLengthList, cecumLength, wallThickness, tcCount, tcThickness,
        #         elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
        #         elementsCountThroughWall, transitElementList)
        #
        #     # Create nodes and elements
        #     nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
        #         region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
        #         elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
        #         tcCount, annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
        #         useCubicHermiteThroughWall, useCrossDerivatives)
        #
        # else:
        #     # Create flat and texture coordinates
        #     xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = tubemesh.createFlatAndTextureCoordinates(
        #         xiList, relaxedLengthList, length, wallThickness, elementsCountAround,
        #         elementsCountAlong, elementsCountThroughWall, transitElementList)
        #
        #     # Create nodes and elements
        #     nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
        #         region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
        #         elementsCountAround, elementsCountAlong, elementsCountThroughWall,
        #         annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
        #         useCubicHermiteThroughWall, useCrossDerivatives)
        #
        return annotationGroups
        # ########################################################################################################
        # nodeIdentifier = nextNodeIdentifier
        # for n in range(len(sxCecum)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sxCecum[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, 0.0])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, sd1Cecum[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, [50*sd2Cecum[n][c] for c in range(3)])
        #     # print('nodeIdentifier = ', nodeIdentifier, 'sd2Cecum = ', sd2Cecum[n])
        #     nodeIdentifier = nodeIdentifier + 1
        #
        # fm.endChange()
        # ######################################################################################################
        #
        # return

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
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()

def getApexSegmentForCecum(xInner, d1Inner, d2Inner, elementsCountAroundHalfHaustrum,
                           elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment, tcCount):
    """
    Generates the inner coordinates and derivatives for a cecum segment on the closed end.
    The closed end is a single node and segment is created by sampling curves between the
    point on closed end with nodes on the length along second half of a colon segment.
    :param xInner: coordinates of a colon segment.
    :param d1Inner: derivative around colon segment.
    :param d2Inner: derivative along colon segment.
    :param elementsCountAroundHalfHaustrum: half of total number of elements in haustrum and tenia coli.
    :param elementsCountAroundTC: number of elements around a tenia coli.
    :param elementsCountAround: number of elements around a cecum.
    :param elementsCountAlongSegment: number of elements along a segment of cecum.
    :param tcCount: number of tenia coli.
    :return: coordinates and derivatives around and along closed segment of cecum, and directional derivative for node
    on middle of tenia coli.
    """

    # Make apex cap - apex points like multiple points
    xFirstSegment = [[0.0, 0.0, 0.0] for c in range(elementsCountAround)]

    # Compile nodes and d2 for sampling
    xFirstSegment += xInner[elementsCountAround * int(elementsCountAlongSegment * 0.5):] # second half of first regular segment
    d1FirstDirectionVector = vector.normalise(d1Inner[elementsCountAround]) # Store direction vector of first d1 intra-haustral for later
    d2Vector = xInner[elementsCountAround*int(elementsCountAlongSegment*0.5):elementsCountAround*(int(elementsCountAlongSegment*0.5)+1)] # half face of segment - apex
    d2FirstSegment = []
    for c in range(elementsCountAround):
        d2 = [d2Vector[c][0], d2Vector[c][1], 0.0 ]  # project onto x-y plane to get d2 pointing vertically
        d2FirstSegment.append(d2)
    d2FirstSegment += d2Inner[elementsCountAround*int(elementsCountAlongSegment*0.5):]

    # Sample along first segment
    xFirstSegmentSampledRaw = []
    d2FirstSegmentSampledRaw = []
    xFirstSegmentSampled = []
    d1FirstSegmentSampled = []
    d2FirstSegmentSampled = []

    for n1 in range(elementsCountAround):
        xForSamplingAlong = []
        d2ForSamplingAlong = []
        for n2 in range(1 + elementsCountAlongSegment - int(elementsCountAlongSegment*0.5) + 1):
            idx = elementsCountAround*n2 + n1
            xForSamplingAlong.append(xFirstSegment[idx])
            d2ForSamplingAlong.append(d2FirstSegment[idx])
        xResampled, d1Resampled, se, sxi, _ = interp.sampleCubicHermiteCurves(xForSamplingAlong, d2ForSamplingAlong,
                                                                              elementsCountAlongSegment,
                                                                              arcLengthDerivatives = True)

        xFirstSegmentSampledRaw.append(xResampled)
        d2FirstSegmentSampledRaw.append(d1Resampled)

    # Re-arrange sample order
    for n2 in range(elementsCountAlongSegment + 1):
        xAround = []
        for n1 in range(elementsCountAround):
            x = xFirstSegmentSampledRaw[n1][n2]
            d2 = d2FirstSegmentSampledRaw[n1][n2]
            if n1 < elementsCountAroundHalfHaustrum + 1:
                xAround.append(x)
            xFirstSegmentSampled.append(x)
            d2FirstSegmentSampled.append(d2)
            if n2 == 0:
                d1 = matrix.rotateAboutZAxis(d2, math.pi*0.5)
                d1FirstSegmentSampled.append(d1)

        if n2 > 0:
            d1Around = []
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c*arcLengthAround for c in vector.normalise(d1)]
                d1Around.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = vector.normalise(
                [xAround[elementsCountAroundHalfHaustrum][c] - xAround[elementsCountAroundHalfHaustrum - 1][c]
                 for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1Around.append(dx_ds1)

            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1Around, fixStartDerivative=True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC * 0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC * 0.5 + 2):]
            d1Full = getD1ForFullProfileFromHalfHaustrum(d1Corrected, tcCount)
            d1FirstSegmentSampled += d1Full

    return xFirstSegmentSampled, d1FirstSegmentSampled, d2FirstSegmentSampled, d1FirstDirectionVector

def getD1ForFullProfileFromHalfHaustrum(d1HaustrumHalfSet, tcCount):
    """
    Get full profile from half haustrum
    :param d1HaustrumHalfSet:
    :param tcCount:
    :return:
    """
    d1HaustrumHalfSet2 = []
    d1Haustra = []

    rotAng = 2 * math.pi / tcCount
    for n in range(1, len(d1HaustrumHalfSet)):
        idx = -n + len(d1HaustrumHalfSet) - 1
        d1 = d1HaustrumHalfSet[idx]
        d1Reflect = [d1[0], -d1[1], d1[2]]
        d1Rot = [-(d1Reflect[0] * math.cos(rotAng) - d1Reflect[1] * math.sin(rotAng)),
                 -(d1Reflect[0] * math.sin(rotAng) + d1Reflect[1] * math.cos(rotAng)),
                 -d1Reflect[2]]
        d1HaustrumHalfSet2.append(d1Rot)

    d1Haustrum = d1HaustrumHalfSet + d1HaustrumHalfSet2

    # Rotate to get all 3 sectors
    d1Haustra = d1Haustra + d1Haustrum[:-1]
    ang = [2 / 3 * math.pi, -2 / 3 * math.pi] if tcCount == 3 else [math.pi]
    for i in range(tcCount - 1):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(d1Haustrum) - 1):
            d1 = d1Haustrum[n]
            dx_ds1 = [d1[0] * cosRotAng - d1[1] * sinRotAng, d1[0] * sinRotAng + d1[1] * cosRotAng, d1[2]]
            d1Haustra.append(dx_ds1)

    return d1Haustra

def sampleCecumAlongLength(xToSample, d1ToSample, d2ToSample, d1FirstDirectionVector, elementsCountAroundHalfHaustrum,
                           elementsCountAroundTC, elementsCountAround, elementsCountAlong, tcCount):
    """
    Get systematically spaced points and derivatives over cubic Hermite interpolated curves along the
    length of the cecum.
    :param xToSample: coordinates of nodes.
    :param d1ToSample: derivative around elements.
    :param d2ToSample: derivative along cecum length.
    :param d1FirstDirectionVector: directional vector of derivative around for node on the middle of the tenia coli.
    :param elementsCountAroundHalfHaustrum:half the total number of elements around tenia coli and haustrum.
    :param elementsCountAroundTC: number of elements around tenia coli.
    :param elementsCountAround: number of elements around cecum.
    :param elementsCountAlong: number of elements along cecum length.
    :param tcCount: number of tenia coli.
    :return: nodes and derivatives for equally spaced points.
    """

    xInnerRaw = []
    d2InnerRaw = []
    xSampledAlongLength = []
    d1SampledAlongLength = []
    d2SampledAlongLength = []

    for n1 in range(elementsCountAroundHalfHaustrum + 1):
        xForSamplingAlong = []
        d2ForSamplingAlong = []
        for n2 in range(elementsCountAlong + 1):
            idx = n2 * elementsCountAround + n1
            xForSamplingAlong.append(xToSample[idx])
            d2ForSamplingAlong.append(d2ToSample[idx])
        xSampled, d2Sampled, se, sxi, _ = interp.sampleCubicHermiteCurves(xForSamplingAlong, d2ForSamplingAlong,
                                                                          elementsCountAlong,
                                                                          arcLengthDerivatives=True)
        xInnerRaw.append(xSampled)
        d2InnerRaw.append(d2Sampled)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    for n2 in range(elementsCountAlong + 1):
        xAround = []
        d2Around = []

        for n1 in range(elementsCountAroundHalfHaustrum + 1):
            x = xInnerRaw[n1][n2]
            d2 = d2InnerRaw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)

        d1InnerAroundList = []
        if n2 == 0:
            d1Corrected = d1ToSample[:elementsCountAroundHalfHaustrum + 1]

        else:
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c * arcLengthAround for c in vector.normalise(d1)]
                d1InnerAroundList.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = vector.normalise([xAround[elementsCountAroundHalfHaustrum][c] -
                                   xAround[elementsCountAroundHalfHaustrum - 1][c] for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1InnerAroundList.append(dx_ds1)

        if d1InnerAroundList:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1InnerAroundList, fixStartDerivative=True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC * 0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC * 0.5 + 2):]

        xAlongList, d1AlongList, d2AlongList = getFullProfileFromHalfHaustrum(xAround, d1Corrected, d2Around, tcCount)

        xSampledAlongLength += xAlongList
        d1SampledAlongLength += d1AlongList
        d2SampledAlongLength += d2AlongList

    return xSampledAlongLength, d1SampledAlongLength, d2SampledAlongLength
