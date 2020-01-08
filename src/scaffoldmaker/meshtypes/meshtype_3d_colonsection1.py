"""
Generates a 3-D partial colon mesh along the central line,
with variable numbers of elements around, along and through
wall, with variable radius and thickness along.
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, ColonSegmentTubeMeshInnerPoints, getTeniaColi, createFlatAndTextureCoordinatesTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.node import Node
from opencmiss.zinc.field import Field # KM

class MeshType_3d_colonsection1(Scaffold_base):
    '''
    Generates a 3-D partial colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a colon
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 90.0,
                'Number of elements': 3
            },
            'meshEdits': zinc_utils.exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [30.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[30.0, 0.0, 0.0], [30.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[60.0, 0.0, 0.0], [30.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[90.0, 0.0, 0.0], [30.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]])
            })
        }

    @staticmethod
    def getName():
        return '3D Colon Section 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):

        centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Pig 1')

        options = {
            'Central path' : copy.deepcopy(centralPathOption),
            'Segment profile' : segmentProfileOption,
            'Number of segments': 3,
            'Start phase': 0.0,
            'Proximal length': 30.0,
            'Transverse length': 30.0,
            'Distal length': 30.0,
            'Proximal inner radius': 16.0,
            'Proximal tenia coli width': 5.0,
            'Proximal-transverse inner radius': 16.0,
            'Proximal-transverse tenia coli width': 5.0,
            'Transverse-distal inner radius': 16.0,
            'Transverse-distal tenia coli width': 5.0,
            'Distal inner radius': 16.0,
            'Distal tenia coli width': 5.0,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Segment profile',
            'Number of segments',
            'Start phase',
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Segment profile':
            return [ MeshType_3d_colonsegment1 ]
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
        if optionName == 'Segment profile':
            if not parameterSetName:
                parameterSetName = scaffoldType.getParameterSetNames()[0]
            return ScaffoldPackage(scaffoldType, defaultParameterSetName = parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Segment profile'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Segment profile'):
            options['Segment profile'] = cls.getOptionScaffoldPackage('Segment profile', MeshType_3d_colonsegmentteniacoli1)
        for key in [
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width']:
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
        segmentProfile = options['Segment profile']
        segmentCount = options['Number of segments']
        startPhase = options['Start phase'] % 360.0
        proximalLength = options['Proximal length']
        transverseLength = options['Transverse length']
        distalLength = options['Distal length']
        proximalInnerRadius = options['Proximal inner radius']
        proximalTCWidth = options['Proximal tenia coli width']
        proximalTransverseInnerRadius = options['Proximal-transverse inner radius']
        proximalTransverseTCWidth = options['Proximal-transverse tenia coli width']
        transverseDistalInnerRadius = options['Transverse-distal inner radius']
        transverseDistalTCWidth = options['Transverse-distal tenia coli width']
        distalInnerRadius = options['Distal inner radius']
        distalTCWidth = options['Distal tenia coli width']
        segmentSettings = segmentProfile.getScaffoldSettings()

        elementsCountAroundTC = segmentSettings['Number of elements around tenia coli']
        elementsCountAroundHaustrum = segmentSettings['Number of elements around haustrum']
        cornerInnerRadiusFactor = segmentSettings['Corner inner radius factor']
        haustrumInnerRadiusFactor = segmentSettings['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = segmentSettings['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = segmentSettings['Segment length mid derivative factor']
        tcCount = segmentSettings['Number of tenia coli']
        tcThickness = segmentSettings['Tenia coli thickness']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*tcCount

        elementsCountAlongSegment = segmentSettings['Number of elements along segment']
        elementsCountThroughWall = segmentSettings['Number of elements through wall']
        wallThickness = segmentSettings['Wall thickness']
        useCrossDerivatives = segmentSettings['Use cross derivatives']
        useCubicHermiteThroughWall = not(segmentSettings['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment*segmentCount)

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i],',', cd12[i], '],')
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections = True,
            magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            # print(e+1, arcLength)
            length += arcLength
        segmentLength = length / segmentCount
        print('Length = ', length)
        # print('segmentLength = ', segmentLength)

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)[0]

        # Create length list where inter-haustral creases sit
        lengthAlongSegmentElement = length / elementsCountAlong
        segmentCounter = segmentCount + 1 if startPhase > 0.0 else segmentCount
        faceIdx = int(startPhase / 360.0 * elementsCountAlongSegment)  # should be phaseCount
        # print(faceIdx)
        elementsCountAlongFirstSegment = elementsCountAlongSegment - faceIdx
        elementsCountAlongLastSegment = faceIdx
        outLengthList = [ 0.0, elementsCountAlongFirstSegment*lengthAlongSegmentElement] # Segment 1
        outFaceIdxList = [ 0, elementsCountAlongFirstSegment ]
        for nSegment in range(segmentCounter - 2):
            outLengthList += [elementsCountAlongFirstSegment*lengthAlongSegmentElement + (nSegment + 1)*segmentLength]
            outFaceIdxList += [elementsCountAlongFirstSegment + (nSegment + 1)*elementsCountAlongSegment]
        outLengthList += [outLengthList[-1] + (segmentLength if elementsCountAlongLastSegment == 0 else
                                               elementsCountAlongLastSegment *lengthAlongSegmentElement)] # Last segment
        outFaceIdxList += [outFaceIdxList[-1] + (elementsCountAlongSegment if elementsCountAlongLastSegment == 0 else
                                               elementsCountAlongLastSegment)]  # Last segment
        # print('outLengthList = ', outLengthList)
        # print('outFaceIdxList = ', outFaceIdxList)

        # Generate variation of radius & tc width along length
        lengthList = [0.0, proximalLength, proximalLength + transverseLength, length]
        innerRadiusList = [proximalInnerRadius, proximalTransverseInnerRadius, transverseDistalInnerRadius, distalInnerRadius]
        innerRadiusSegmentList, dInnerRadiusSegmentList = interp.getParameterAlongLine(lengthList, innerRadiusList, outLengthList)

        tcWidthList = [proximalTCWidth, proximalTransverseTCWidth, transverseDistalTCWidth, distalTCWidth]
        tcWidthSegmentList, dTCWidthSegmentList = interp.getParameterAlongLine(lengthList, tcWidthList, outLengthList)

        xExtrude = []
        d1Extrude = []
        d2Extrude = []
        d3UnitExtrude = []

        # Create object
        colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactor,
            innerRadiusSegmentList, dInnerRadiusSegmentList, tcWidthSegmentList, dTCWidthSegmentList)

        # nodeIdentifier = firstNodeIdentifier
        # elementIdentifier = firstElementIdentifier
        # zero = [0.0, 0.0, 0.0]
        #
        # fm = region.getFieldmodule()
        # fm.beginChange()
        # cache = fm.createFieldcache()
        #
        # # Coordinates field
        # coordinates = zinc_utils.getOrCreateCoordinateField(fm)
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

        for nSegment in range(segmentCounter):
            # Create inner points
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroups, annotationArray = \
                colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

            faceMidPointsZ = []

            # Create first segment separately as it can start from start or in a different phase
            if nSegment == 0 and startPhase > 0.0:
                xInnerPhase = xInner[faceIdx*elementsCountAround:]
                d1InnerPhase = d1Inner[faceIdx * elementsCountAround:]
                d2InnerPhase = d2Inner[faceIdx * elementsCountAround:]
                for n2 in range(elementsCountAlongFirstSegment + 1):
                    faceMidPointsZ += [faceIdx * lengthAlongSegmentElement + n2*lengthAlongSegmentElement]
            elif nSegment > segmentCounter - 1 and startPhase > 0.0:
                xInnerPhase = xInner[: (faceIdx + 1) * elementsCountAround]
                d1InnerPhase = d1Inner[: (faceIdx + 1) * elementsCountAround]
                d2InnerPhase = d2Inner[: (faceIdx + 1) * elementsCountAround]
                for n2 in range(elementsCountAlongLastSegment + 1):
                    faceMidPointsZ += [n2*lengthAlongSegmentElement]
            else:
                xInnerPhase = xInner
                d1InnerPhase = d1Inner
                d2InnerPhase = d2Inner
                for n2 in range(elementsCountAlongSegment + 1):
                    faceMidPointsZ += [n2*lengthAlongSegmentElement]

            # Warp segment points
            nSegmentStart = outFaceIdxList[nSegment]
            nSegmentEnd = outFaceIdxList[nSegment + 1]
            xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
                xInnerPhase, d1InnerPhase, d2InnerPhase, segmentAxis, segmentLength, sx, sd1, sd2,
                elementsCountAround, elementsCountAlongSegment, nSegmentStart, nSegmentEnd, faceMidPointsZ)

            # # Create nodes
            # # Coordinates field
            # for n in range(len(xWarpedList)):
            #     node = nodes.createNode(nodeIdentifier, nodetemplate)
            #     cache.setNode(node)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xWarpedList[n])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1WarpedList[n])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2WarpedList[n])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            #     if useCrossDerivatives:
            #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            #     # print('NodeIdentifier = ', nodeIdentifier, xList[n])
            #     nodeIdentifier = nodeIdentifier + 1

            # Store points along length
            xExtrude = xExtrude + (xWarpedList if nSegment == 0 else xWarpedList[elementsCountAround:])
            d1Extrude = d1Extrude + (d1WarpedList if nSegment == 0 else d1WarpedList[elementsCountAround:])

            # Smooth d2 for nodes between segments and recalculate d3
            if nSegment == 0:
                d2Extrude = d2Extrude + (d2WarpedList[:-elementsCountAround])
                d3UnitExtrude = d3UnitExtrude + (d3WarpedUnitList[:-elementsCountAround])
            else:
                xSecondFace = xWarpedList[elementsCountAround:elementsCountAround*2]
                d2SecondFace = d2WarpedList[elementsCountAround:elementsCountAround*2]
                for n1 in range(elementsCountAround):
                    nx = [xLastTwoFaces[n1], xLastTwoFaces[n1 + elementsCountAround], xSecondFace[n1]]
                    nd2 = [d2LastTwoFaces[n1], d2LastTwoFaces[n1 + elementsCountAround], d2SecondFace[n1]]
                    d2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)[1]
                    d2Extrude.append(d2)
                    d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1LastTwoFaces[n1 + elementsCountAround]), vector.normalise(d2)))
                    d3UnitExtrude.append(d3Unit)
                d2Extrude = d2Extrude + (d2WarpedList[elementsCountAround:-elementsCountAround] if nSegment < segmentCounter - 1 else d2WarpedList[elementsCountAround:])
                d3UnitExtrude = d3UnitExtrude + (d3WarpedUnitList[elementsCountAround:-elementsCountAround] if nSegment < segmentCounter - 1 else d3WarpedUnitList[elementsCountAround:])
            xLastTwoFaces = xWarpedList[-elementsCountAround*2:]
            d1LastTwoFaces = d1WarpedList[-elementsCountAround*2:]
            d2LastTwoFaces = d2WarpedList[-elementsCountAround*2:]

        contractedWallThicknessList = colonSegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xExtrude, d1Extrude,
            d2Extrude, d3UnitExtrude, contractedWallThicknessList,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        relaxedLengthList, xiList = colonSegmentTubeMeshInnerPoints.getRelaxedLengthAndXiList()

        if tcThickness > 0:
            tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
            xList, d1List, d2List, d3List, annotationGroups, annotationArray = getTeniaColi(
                region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
                elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tubeTCWidthList, tcThickness, sx, annotationGroups, annotationArray)

            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = createFlatAndTextureCoordinatesTeniaColi(
                xiList, relaxedLengthList, length, wallThickness, tcCount, tcThickness,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
                elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tcCount, annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
                useCubicHermiteThroughWall, useCrossDerivatives)

        else:
            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = tubemesh.createFlatAndTextureCoordinates(
                xiList, relaxedLengthList, length, wallThickness, elementsCountAround,
                elementsCountAlong, elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
                useCubicHermiteThroughWall, useCrossDerivatives)

        # fm.endChange() # Km

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
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()
