"""
Generates a single 3-D colon segment mesh along a central
line, with variable numbers of elements around, along and
through wall and number of tenia coli, with variable radius
and thickness along.
"""

import copy
import math

from cmlibs.maths.vectorops import normalize, set_magnitude, magnitude, dot, cross, rotate_about_z_axis
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups, \
    findOrCreateAnnotationGroupForTerm, findAnnotationGroupByName, getAnnotationGroupForTerm
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.annotation.cecum_terms import get_cecum_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import createCirclePoints


class MeshType_3d_colonsegment1(Scaffold_base):
    """
    Generates a single 3-D colon segment mesh with variable
    numbers of tenia coli, elements around, along the central
    line, and through wall. The cross-section profile of the colon
    segment varies with species and is dependent on the number
    of tenia coli.
    Mouse: One "tenia coli", circular profile
    Pig: Two tenia coli, bow tie profile
    Human (Default): Three tenia coli, triangular profile with
    rounded corners at the inter-haustral septa, and a clover
    profile in the intra-haustral region.
    """

    @staticmethod
    def getName():
        return '3D Colon Segment 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cattle 1',
            'Human 1',
            'Human 2',
            'Mouse 1',
            'Pig 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {
            'Base parameter set': parameterSetName,
            'Number of elements around tenia coli': 2,
            'Number of elements around haustrum': 8,
            'Number of elements along segment': 4,
            'Number of elements through wall': 4,
            'Start phase': 0.0,
            'Start outer radius': 45.4,
            'Start outer radius derivative': 0.0,
            'End outer radius': 35.3,
            'End outer radius derivative': 0.0,
            'Corner outer radius factor': 0.536,
            'Haustrum outer radius factor': 0.464,
            'Segment length end derivative factor': 0.5,
            'Segment length mid derivative factor': 3.0,
            'Segment length': 50.0,
            'Number of tenia coli': 3,
            'Start tenia coli width': 11.1,
            'Start tenia coli width derivative': 0.4,
            'End tenia coli width': 11.5,
            'End tenia coli width derivative': 0.4,
            'Tenia coli thickness': 0.6,
            'Wall thickness': 1.6,
            'Mucosa relative thickness': 0.18,
            'Submucosa relative thickness': 0.25,
            'Circular muscle layer relative thickness': 0.52,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along segment': 1,
            'Refine number of elements through wall': 1
        }
        if 'Cattle' in parameterSetName:
            options['Start outer radius'] = 13.52
            options['End outer radius'] = 13.52
            options['Corner outer radius factor'] = 0.0
            options['Haustrum outer radius factor'] = 0.0
            options['Segment length end derivative factor'] = 0.0
            options['Segment length mid derivative factor'] = 0.0
            options['Number of tenia coli'] = 1
            options['Start tenia coli width'] = 3.85
            options['End tenia coli width'] = 3.85
            options['Tenia coli thickness'] = 0.0
            options['Wall thickness'] = 3.02
        elif 'Human 2' in parameterSetName:
            options['Number of elements through wall'] = 1
            options['Haustrum outer radius factor'] = 0.364
            options['Tenia coli thickness'] = 1.6
        elif 'Mouse' in parameterSetName:
            options['Start outer radius'] = 1.49
            options['End outer radius'] = 1.49
            options['Corner outer radius factor'] = 0.0
            options['Haustrum outer radius factor'] = 0.0
            options['Segment length end derivative factor'] = 0.0
            options['Segment length mid derivative factor'] = 0.0
            options['Number of tenia coli'] = 1
            options['Start tenia coli width'] = 1.29
            options['End tenia coli width'] = 1.29
            options['Tenia coli thickness'] = 0.0
            options['Wall thickness'] = 0.55
            options['Mucosa relative thickness'] = 0.4
            options['Submucosa relative thickness'] = 0.1
            options['Circular muscle layer relative thickness'] = 0.3
            options['Longitudinal muscle layer relative thickness'] = 0.2
        elif 'Pig' in parameterSetName:
            options['Start outer radius'] = 22.23
            options['End outer radius'] = 22.23
            options['Corner outer radius factor'] = 0.0
            options['Haustrum outer radius factor'] = 0.17
            options['Segment length end derivative factor'] = 0.8
            options['Segment length mid derivative factor'] = 2.0
            options['Segment length'] = 25.0
            options['Number of tenia coli'] = 2
            options['Start tenia coli width'] = 5.5
            options['End tenia coli width'] = 5.5
            options['Tenia coli thickness'] = 0.5
            options['Wall thickness'] = 2.0
            options['Mucosa relative thickness'] = 0.34
            options['Submucosa relative thickness'] = 0.25
            options['Circular muscle layer relative thickness'] = 0.25
            options['Longitudinal muscle layer relative thickness'] = 0.16

        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Start phase',
            'Start outer radius',
            'Start outer radius derivative',
            'End outer radius',
            'End outer radius derivative',
            'Corner outer radius factor',
            'Haustrum outer radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Segment length',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements along segment']:
            if options[key] < 2:
                options[key] = 2
        if options['Number of elements around haustrum'] < 4:
            options['Number of elements around haustrum'] = 4
        if options['Number of elements through wall'] != (1 or 4):
            options['Number of elements through wall'] = 4
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum']:
            if options[key] % 2 > 0:
                options[key] = options[key] + 1
        for key in [
            'Start outer radius',
            'End outer radius',
            'Haustrum outer radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Segment length',
            'Tenia coli thickness',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner outer radius factor'] < 0.1:
            options['Corner outer radius factor'] = 0.1
        for key in [
            'Corner outer radius factor',
            'Segment length end derivative factor']:
            if options[key] > 1.0:
                options[key] = 1.0
        if options['Number of tenia coli'] < 1:
            options['Number of tenia coli'] = 1
        elif options['Number of tenia coli'] > 3:
            options['Number of tenia coli'] = 3
        for key in [
            'Start tenia coli width',
            'End tenia coli width']:
            if options[key] < 0.2 * min(options['Start outer radius'], options['End outer radius']):
                options[key] = round(0.2 * min(options['Start outer radius'], options['End outer radius']), 2)
        if options['Start tenia coli width'] > round(math.sqrt(3) * 0.5 * options['Start outer radius'], 2):
            options['Start tenia coli width'] = round(math.sqrt(3) * 0.5 * options['Start outer radius'], 2)
        if options['End tenia coli width'] > round(math.sqrt(3) * 0.5 * options['End outer radius'], 2):
            options['End tenia coli width'] = round(math.sqrt(3) * 0.5 * options['End outer radius'], 2)

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        parameterSetName = options['Base parameter set']
        isHuman = 'Human' in parameterSetName
        elementsCountAroundTC = options['Number of elements around tenia coli']
        elementsCountAroundHaustrum = options['Number of elements around haustrum']
        elementsCountAlongSegment = options['Number of elements along segment']
        elementsCountThroughWall = options['Number of elements through wall']
        startPhase = options['Start phase'] % 360.0
        startRadius = options['Start outer radius']
        startRadiusDerivative = options['Start outer radius derivative']
        endRadius = options['End outer radius']
        endRadiusDerivative = options['End outer radius derivative']
        cornerOuterRadiusFactor = options['Corner outer radius factor']
        haustrumOuterRadiusFactor = options['Haustrum outer radius factor']
        segmentLengthEndDerivativeFactor = options['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = options['Segment length mid derivative factor']
        segmentLength = options['Segment length']
        tcCount = options['Number of tenia coli']
        startTCWidth = options['Start tenia coli width']
        startTCWidthDerivative = options['Start tenia coli width derivative']
        endTCWidth = options['End tenia coli width']
        endTCWidthDerivative = options['End tenia coli width derivative']
        tcThickness = options['Tenia coli thickness']
        wallThickness = options['Wall thickness']
        mucosaRelThickness = options['Mucosa relative thickness']
        submucosaRelThickness = options['Submucosa relative thickness']
        circularRelThickness = options['Circular muscle layer relative thickness']
        longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not (options['Use linear through wall'])
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount
        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Network layout
        cx = [[0.0, 0.0, 0.0], [segmentLength, 0.0, 0.0]]
        cd1 = [[segmentLength, 0.0, 0.0], [segmentLength, 0.0, 0.0]]
        cd2 = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        cd12 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        # Sample network layout
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Find parameter variation along elementsCountAlongSegment
        radiusAlongSegment = []
        dRadiusAlongSegment = []
        tcWidthAlongSegment = []

        for n2 in range(elementsCountAlongSegment + 1):
            xi = 1 / elementsCountAlongSegment * n2
            radius = interp.interpolateCubicHermite([startRadius], [startRadiusDerivative],
                                                    [endRadius], [endRadiusDerivative], xi)[0]
            radiusAlongSegment.append(radius)
            dRadius = interp.interpolateCubicHermiteDerivative([startRadius], [startRadiusDerivative],
                                                               [endRadius], [endRadiusDerivative], xi)[0]
            dRadiusAlongSegment.append(dRadius)
            tcWidth = interp.interpolateCubicHermite([startTCWidth], [startTCWidthDerivative],
                                                     [endTCWidth], [endTCWidthDerivative], xi)[0]
            tcWidthAlongSegment.append(tcWidth)

        haustrumOuterRadiusFactorAlongSegment = [haustrumOuterRadiusFactor] * (elementsCountAlongSegment + 1)

        colonSegmentTubeMeshOuterPoints = ColonSegmentTubeMeshOuterPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerOuterRadiusFactor, haustrumOuterRadiusFactorAlongSegment,
            radiusAlongSegment, dRadiusAlongSegment, tcWidthAlongSegment, startPhase)

        # Create annotation
        colonGroup = AnnotationGroup(region, get_colon_term("colon"))
        annotationGroupsAlong = []
        for i in range(elementsCountAlongSegment):
            annotationGroupsAlong.append([colonGroup])

        # Create outer points
        nSegment = 0
        closedProximalEnd = False

        xOuter, d1Outer, d2Outer, transitElementList, segmentAxis, annotationGroupsAround = \
            colonSegmentTubeMeshOuterPoints.getColonSegmentTubeMeshOuterPoints(nSegment)

        # Project reference point for warping onto network layout
        sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
            tubemesh.getPlaneProjectionOnCentralPath(xOuter, elementsCountAround, elementsCountAlongSegment,
                                                     segmentLength, sx, sd1, sd2, sd12)

        # Warp segment points
        xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
            xOuter, d1Outer, d2Outer, segmentAxis, sxRefList, sd1RefList, sd2ProjectedListRef,
            elementsCountAround, elementsCountAlongSegment, zRefList)

        contractedWallThicknessList = colonSegmentTubeMeshOuterPoints.getContractedWallThicknessList()

        if elementsCountThroughWall == 1:
            relativeThicknessList = [1.0]
            annotationGroupsThroughWall = [[]]
        else:
            relativeThicknessList = [mucosaRelThickness, submucosaRelThickness,
                                     circularRelThickness, longitudinalRelThickness]
            mucosaGroup = AnnotationGroup(region, get_colon_term("colonic mucosa"))
            submucosaGroup = AnnotationGroup(region, get_colon_term("submucosa of colon"))
            circularMuscleGroup = AnnotationGroup(region, get_colon_term("circular muscle layer of colon"))
            longitudinalMuscleGroup = AnnotationGroup(region, get_colon_term("longitudinal muscle layer of colon"))
            annotationGroupsThroughWall = [[mucosaGroup], [submucosaGroup],
                                           [circularMuscleGroup], [longitudinalMuscleGroup]]

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = \
            tubemesh.extrudeSurfaceCoordinates(xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList,
                                             contractedWallThicknessList, relativeThicknessList, elementsCountAround,
                                             elementsCountAlongSegment, elementsCountThroughWall, transitElementList,
                                               outward=False)[0:5]

        xColonSegment = d1ColonSegment = d2ColonSegment = []

        relaxedLengthList, xiList = colonSegmentTubeMeshOuterPoints.getRelaxedLengthAndXiList()

        if tcThickness > 0:
            tubeTCWidthList = colonSegmentTubeMeshOuterPoints.getTubeTCWidthList()
            xList, d1List, d2List, d3List, annotationGroupsAround = getTeniaColi(
                region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
                elementsCountAroundHaustrum, elementsCountAlongSegment, elementsCountThroughWall,
                tubeTCWidthList, tcThickness, annotationGroupsAround, closedProximalEnd, isHuman)[0:5]

            # Create flat coordinates
            xFlat, d1Flat, d2Flat = createFlatCoordinatesTeniaColi(
                xiList, relaxedLengthList, segmentLength, wallThickness, relativeThicknessList, tcCount, tcThickness,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
                elementsCountThroughWall, transitElementList, closedProximalEnd)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xColonSegment, d1ColonSegment,
                d2ColonSegment, None, elementsCountAroundTC, elementsCountAroundHaustrum,
                elementsCountAlongSegment, elementsCountThroughWall, tcCount, annotationGroupsAround,
                annotationGroupsAlong, annotationGroupsThroughWall, firstNodeIdentifier, firstElementIdentifier,
                useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd)[0:3]
        else:
            # Create flat coordinates
            xFlat, d1Flat, d2Flat = tubemesh.createFlatCoordinates(
                xiList, relaxedLengthList, segmentLength, wallThickness, relativeThicknessList, elementsCountAround,
                elementsCountAlongSegment, elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xColonSegment, d1ColonSegment,
                d2ColonSegment, None, elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall,
                annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
                firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
                closedProximalEnd)[0:3]

        return annotationGroups, None

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along segment']
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
        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        colonGroup = getAnnotationGroupForTerm(annotationGroups, get_colon_term("colon"))
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_colon = colonGroup.getGroup()
        is_serosa = fm.createFieldAnd(is_colon, is_exterior_face_xi3_1)
        serosa = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_colon_term("serosa of colon"))
        serosa.getMeshGroup(mesh2d).addElementsConditional(is_serosa)
        is_mucosaInnerSurface = fm.createFieldAnd(is_colon, is_exterior_face_xi3_0)
        mucosaInnerSurface = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_colon_term("luminal surface of the colonic mucosa"))
        mucosaInnerSurface.getMeshGroup(mesh2d).addElementsConditional(is_mucosaInnerSurface)


class ColonSegmentTubeMeshOuterPoints:
    """
    Generates outer profile of a colon segment for use by tubemesh.
    """

    def __init__(self, region, elementsCountAroundTC, elementsCountAroundHaustrum,
                 elementsCountAlongSegment, tcCount, segmentLengthEndDerivativeFactor,
                 segmentLengthMidDerivativeFactor, segmentLength, wallThickness,
                 cornerOuterRadiusFactor, haustrumOuterRadiusFactorAlongElementList, outerRadiusAlongElementList,
                 dOuterRadiusAlongElementList, tcWidthAlongElementList, startPhase):
        self._region = region
        self._elementsCountAroundTC = elementsCountAroundTC
        self._elementsCountAroundHaustrum = elementsCountAroundHaustrum
        self._elementsCountAlongSegment = elementsCountAlongSegment
        self._tcCount = tcCount
        self._segmentLengthEndDerivativeFactor = segmentLengthEndDerivativeFactor
        self._segmentLengthMidDerivativeFactor = segmentLengthMidDerivativeFactor
        self._segmentLength = segmentLength
        self._wallThickness = wallThickness
        self._cornerOuterRadiusFactor = cornerOuterRadiusFactor
        self._haustrumOuterRadiusFactorAlongElementList = haustrumOuterRadiusFactorAlongElementList
        self._outerRadiusAlongElementList = outerRadiusAlongElementList
        self._dOuterRadiusAlongElementList = dOuterRadiusAlongElementList
        self._tcWidthAlongElementList = tcWidthAlongElementList
        self._tubeTCWidthList = []
        self._xiList = []
        self._relaxedLengthList = []
        self._contractedWallThicknessList = []
        self._startPhase = startPhase

    def getColonSegmentTubeMeshOuterPoints(self, nSegment):
        # Unpack parameter variation along elements
        radiusSegmentList = self._outerRadiusAlongElementList[nSegment * self._elementsCountAlongSegment:
                                                              (nSegment + 1) * self._elementsCountAlongSegment + 1]
        dRadiusSegmentList = self._dOuterRadiusAlongElementList[nSegment * self._elementsCountAlongSegment:
                                                                (nSegment + 1) * self._elementsCountAlongSegment + 1]
        tcWidthSegmentList = self._tcWidthAlongElementList[nSegment * self._elementsCountAlongSegment:
                                                           (nSegment + 1) * self._elementsCountAlongSegment + 1]

        haustrumOuterRadiusFactorSegmentList = self._haustrumOuterRadiusFactorAlongElementList[
                                               nSegment * self._elementsCountAlongSegment:
                                               (nSegment + 1) * self._elementsCountAlongSegment + 1]

        xOuter, d1Outer, d2Outer, transitElementList, xiSegment, relaxedLengthSegment, contractedWallThicknessSegment, \
        segmentAxis, annotationGroupsAround = \
            getColonSegmentOuterPoints(self._region,
                                       self._elementsCountAroundTC, self._elementsCountAroundHaustrum,
                                       self._elementsCountAlongSegment, self._tcCount,
                                       self._segmentLengthEndDerivativeFactor, self._segmentLengthMidDerivativeFactor,
                                       self._segmentLength, self._wallThickness,
                                       self._cornerOuterRadiusFactor, haustrumOuterRadiusFactorSegmentList,
                                       radiusSegmentList, dRadiusSegmentList, tcWidthSegmentList,
                                       self._startPhase)

        startIdx = 0 if nSegment == 0 else 1
        for i in range(startIdx, self._elementsCountAlongSegment + 1):
            self._tubeTCWidthList.append(tcWidthSegmentList[i])

        xi = xiSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._xiList += xi

        relaxedLength = relaxedLengthSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._relaxedLengthList += relaxedLength

        contractedWallThickness = contractedWallThicknessSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._contractedWallThicknessList += contractedWallThickness

        return xOuter, d1Outer, d2Outer, transitElementList, segmentAxis, annotationGroupsAround

    def getTubeTCWidthList(self):
        return self._tubeTCWidthList

    def getRelaxedLengthAndXiList(self):
        return self._relaxedLengthList, self._xiList

    def getContractedWallThicknessList(self):
        return self._contractedWallThicknessList


def getColonSegmentOuterPoints(region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
                               tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
                               segmentLength, wallThickness, cornerOuterRadiusFactor,
                               haustrumOuterRadiusFactorSegmentList, radiusSegmentList, dRadiusSegmentList,
                               tcWidthSegmentList, startPhase):
    """
    Generates a 3-D colon segment mesh with variable numbers of tenia coli,
    numbers of elements around, along the network layout, and through wall.
    Colon segment with one "tenia coli" (mouse) has a circular profile along
    its length. Colon segment with two tenia coli (pig) has a circular profile
    at the inter-haustral septa, and a bowtie profile in the intra-haustral
    region. Colon segment with three tenia coli (human) has a triangular profile
    with rounded corners at the inter-haustral septa, and a clover profile
    in the intra-haustral region.
    :param elementsCountAroundTC: Number of elements around each tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlongSegment: Number of elements along colon segment.
    :param tcCount: Number of tenia coli.
    :param segmentLengthEndDerivativeFactor: Factor is multiplied by segment
    length to scale derivative along the end of a segment length.
    :param segmentLengthMidDerivativeFactor: Factor is multiplied by segment
    length to scale derivative along the mid length of the segment.
    :param segmentLength: Length of a colon segment.
    :param wallThickness: Thickness of wall.
    :param cornerOuterRadiusFactor: Roundness of triangular corners of
    inter-haustral septa. Factor is multiplied by outer radius
    to get a radius of curvature at the corners. Only applicable for three tenia
    coli. Set to zero for two tenia coli.
    :param haustrumOuterRadiusFactorSegmentList: Factor is multiplied by inner
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra segment.
    :param radiusSegmentList: List of inner radius defined from center of triangular
    profile to vertex of the triangle at proximal end of the colon segment for each
    element along.
    :param dRadiusSegmentList: List of rate of change of outer radius at proximal end
    for each element along.
    :param tcWidthSegmentList: List of tenia coli width at proximal end of the colon segment
    for each element along.
    :param startPhase: Phase at start.
    :return coordinates, derivatives on outer surface of a colon segment.
    :return transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :return xiList: List of xi for each node around. xi refers to node position
    along the relaxed length when colon segment is opened into a flat preparation,
    nominally in [0.0, 1.0].
    :return relaxedLengthList: List of relaxed length around elements for each element
    along colon segment when the segment is opened into a flat preparation.
    :return contractedWallThicknessList: List of wall thickness for each element
    along colon segment. Assume incompressiblity and a shortened length around will
    result in a thicker wall and vice-versa.
    :return segmentAxis: Axis of segment.
    :return annotationGroupsAround: annotation groups for elements around.
    """

    transitElementListHaustrum = ([0] * int(elementsCountAroundTC * 0.5) + [1] +
                                  [0] * int(elementsCountAroundHaustrum - 2) + [1] +
                                  [0] * int(elementsCountAroundTC * 0.5))
    transitElementList = transitElementListHaustrum * tcCount

    # create nodes
    sampleElementOut = 20
    segmentAxis = [0.0, 0.0, 1.0]

    d2HalfSet = []
    d2Raw = []
    xOuterRaw = []
    dx_ds2OuterRaw = []
    xFinal = []
    d1Final = []
    d2Final = []
    xiList = []
    relaxedLengthList = []
    contractedWallThicknessList = []
    xForSampling = []
    d1ForSampling = []

    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount

    if tcCount == 1:
        for n2 in range(elementsCountAlongSegment + 1):
            radius = radiusSegmentList[n2]
            tcWidth = tcWidthSegmentList[n2]

            xHalfSet, d1HalfSet = \
                createHalfSetInterHaustralSegment(elementsCountAroundTC, elementsCountAroundHaustrum, tcCount, tcWidth,
                                                  radius, cornerOuterRadiusFactor, sampleElementOut)

            for i in range(len(xHalfSet)):
                d2HalfSet.append([0.0, 0.0, 0.0])
            x, d1, _ = getFullProfileFromHalfHaustrum(xHalfSet, d1HalfSet, d2HalfSet, tcCount)
            z = segmentLength / elementsCountAlongSegment * n2
            d1Final = d1Final + d1
            xFace = []
            for n1 in range(elementsCountAroundTC + elementsCountAroundHaustrum):
                xFinal.append([x[n1][0], x[n1][1], z])
                xFace.append([x[n1][0], x[n1][1], z])
            xiFace, lengthAroundFace = getXiListFromOuterLengthProfile(xFace, d1)
            xiList.append(xiFace)
            relaxedLengthList.append(lengthAroundFace)
            contractedWallThicknessList.append(wallThickness)

        for n1 in range(elementsCountAround):
            xUp = []
            d2Up = []
            for n2 in range(elementsCountAlongSegment + 1):
                n = elementsCountAround * n2 + n1
                xUp.append(xFinal[n])
                d2 = [xFinal[n + elementsCountAround][i] - xFinal[n][i] if n2 < elementsCountAlongSegment else
                      xFinal[n][i] - xFinal[n - elementsCountAround][i] for i in range(3)]
                d2Up.append(d2)
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xUp, d2Up)
            d2Raw.append(d2Smoothed)

        # Re-arrange d2Raw
        for n2 in range(elementsCountAlongSegment + 1):
            for n1 in range(elementsCountAround):
                d2Final.append(d2Raw[n1][n2])

        # Create annotation groups for mouse colon
        mzGroup = AnnotationGroup(region, get_colon_term("mesenteric zone"))
        nonmzGroup = AnnotationGroup(region, get_colon_term("non-mesenteric zone"))
        elementsCountAroundGroups = [int(elementsCountAroundTC * 0.5),
                                     elementsCountAroundHaustrum,
                                     int(elementsCountAroundTC * 0.5)]

        annotationGroupAround = [[mzGroup], [nonmzGroup], [mzGroup]]

        annotationGroupsAround = []
        for i in range(len(elementsCountAroundGroups)):
            elementsCount = elementsCountAroundGroups[i]
            for n in range(elementsCount):
                annotationGroupsAround.append(annotationGroupAround[i])

    else:
        elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum) * 0.5)
        d1AtStartOfEachMidFace = []
        d1Corrected = []

        for n2 in range(elementsCountAlongSegment + 1):
            radius = radiusSegmentList[n2]
            sdRadius = dRadiusSegmentList[n2]
            tcWidth = tcWidthSegmentList[n2]
            haustrumOuterRadiusFactor = haustrumOuterRadiusFactorSegmentList[n2]

            # Create segment of inner radius
            # Calculate x and d1 at the start, mid, and end faces
            xHalfSetStart, d1HalfSetStart = createHalfSetInterHaustralSegment(
                elementsCountAroundTC, elementsCountAroundHaustrum, tcCount, tcWidth, radius,
                cornerOuterRadiusFactor, sampleElementOut)

            if startPhase == 0.0:
                if n2 == 0:
                    d1Phase0FirstFace = copy.deepcopy(d1HalfSetStart)
                elif n2 > elementsCountAlongSegment - 1:
                    d1Phase0LastFace = copy.deepcopy(d1HalfSetStart)

            if startPhase == 180.0:
                if elementsCountAlongSegment % 2 == 0 and n2 == int(elementsCountAlongSegment * 0.5):
                    d1180MidFace = copy.deepcopy(d1HalfSetStart)

            xHalfSetMid, d1HalfSetMid = createHalfSetIntraHaustralSegment(
                elementsCountAroundTC, elementsCountAroundHaustrum, tcCount, tcWidth, radius,
                cornerOuterRadiusFactor, sampleElementOut, haustrumOuterRadiusFactor)

            d1AtStartOfEachMidFace.append(d1HalfSetMid[0])

            if startPhase == 0.0 and elementsCountAlongSegment % 2 == 0 and n2 == int(elementsCountAlongSegment * 0.5):
                d1Phase0MidFace = copy.deepcopy(d1HalfSetMid)

            if startPhase == 180.0:
                if n2 == 0:
                    d1180FirstFace = copy.deepcopy(d1HalfSetMid)
                if n2 > elementsCountAlongSegment - 1:
                    d1180LastFace = copy.deepcopy(d1HalfSetMid)

            # Sample arclength of haustra segment
            xFace = []
            d1Face = []
            for n1 in range(elementsCountAroundHalfHaustrum + 1):
                radiansAround = (math.pi * 2 / (tcCount * 2)) / elementsCountAroundHalfHaustrum * n1
                d2Start = [sdRadius * 0.5, 0.0, segmentLength * 0.5]
                d2Mid = [sdRadius * 0.5, 0.0, segmentLength * 0.5]
                d2End = [sdRadius * 0.5, 0.0, segmentLength * 0.5]

                if n1 > elementsCountAroundTC * 0.5:
                    # Non tenia coli
                    # Rotate about segment axis (z-axis)
                    d2StartRot = rotate_about_z_axis(d2Start, radiansAround)
                    d2MidRot = rotate_about_z_axis(d2Mid, radiansAround)
                    d2EndRot = rotate_about_z_axis(d2End, radiansAround)

                    d1 = [c * segmentLengthEndDerivativeFactor for c in d2StartRot]
                    d2 = [c * segmentLengthMidDerivativeFactor for c in d2MidRot]
                    d3 = [c * segmentLengthEndDerivativeFactor for c in d2EndRot]

                else:
                    # Tenia coli do not have haustra so not subjected to derivative scaling along segment
                    d1 = d2Start
                    d2 = d2Mid
                    d3 = d2End

                v1 = [xHalfSetStart[n1][0], xHalfSetStart[n1][1], 0.0]
                v2 = [xHalfSetMid[n1][0], xHalfSetMid[n1][1], segmentLength / 2]
                v3 = [xHalfSetStart[n1][0], xHalfSetStart[n1][1], segmentLength]

                nx = [v1, v2, v3]
                nd1 = [d1, d2, d3]

                sx, sd1, se, sxi, _ = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAlongSegment)

                # Find start phase
                faceStartPhase = (360.0 / elementsCountAlongSegment * n2 + startPhase) % 360.0
                offsetLength = (360.0 / elementsCountAlongSegment * n2 + startPhase) // 360.0

                # Find start point
                phasePerElementAlongSegment = 360.0 / elementsCountAlongSegment
                downstreamStartFaceIdx = int(math.floor(faceStartPhase / phasePerElementAlongSegment))
                xiStartPhase = (faceStartPhase - downstreamStartFaceIdx * phasePerElementAlongSegment) / \
                               phasePerElementAlongSegment
                xStart = interp.interpolateCubicHermite(sx[downstreamStartFaceIdx], sd1[downstreamStartFaceIdx],
                                                        sx[downstreamStartFaceIdx + 1], sd1[downstreamStartFaceIdx + 1],
                                                        xiStartPhase)
                xStart[2] = xStart[2] + offsetLength * segmentLength
                d1Start = interp.interpolateCubicHermiteDerivative(sx[downstreamStartFaceIdx],
                                                                   sd1[downstreamStartFaceIdx],
                                                                   sx[downstreamStartFaceIdx + 1],
                                                                   sd1[downstreamStartFaceIdx + 1], xiStartPhase)
                xFace.append(xStart)
                d1Face.append(d1Start)

            # Collate points for re-sampling
            xForSampling.append(xFace)
            d1ForSampling.append(d1Face)

        # Re-sample to have points spread out evenly
        for n1 in range(elementsCountAroundHalfHaustrum + 1):
            xForSamplingAlong = []
            d1ForSamplingAlong = []
            for n2 in range(elementsCountAlongSegment + 1):
                xForSamplingAlong.append(xForSampling[n2][n1])
                d1ForSamplingAlong.append(d1ForSampling[n2][n1])
            xResampled, d1Resampled, se, sxi, _ = interp.sampleCubicHermiteCurves(xForSamplingAlong, d1ForSamplingAlong,
                                                                                  elementsCountAlongSegment,
                                                                                  arcLengthDerivatives=True)
            xOuterRaw.append(xResampled)
            dx_ds2OuterRaw.append(d1Resampled)

        # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
        for n2 in range(elementsCountAlongSegment + 1):
            lengthToFirstPhase = startPhase / 360.0 * segmentLength

            xAround = []
            d2Around = []

            for n1 in range(elementsCountAroundHalfHaustrum + 1):
                x = xOuterRaw[n1][n2]
                # Bring first face back to origin
                x = [x[0], x[1], x[2] - lengthToFirstPhase]
                dx_ds2 = dx_ds2OuterRaw[n1][n2]
                xAround.append(x)
                d2Around.append(dx_ds2)

            dx_ds1OuterAroundList = []
            if startPhase == 0.0 and n2 == 0:
                d1Corrected = d1Phase0FirstFace
            elif startPhase == 0.0 and elementsCountAlongSegment % 2 == 0 and \
                    n2 == int(elementsCountAlongSegment * 0.5):
                dx_ds1OuterAroundList = dx_ds1OuterAroundList + d1Phase0MidFace
            elif startPhase == 0.0 and n2 > elementsCountAlongSegment - 1:
                d1Corrected = d1Phase0LastFace

            elif startPhase == 180.0 and n2 == 0:
                dx_ds1OuterAroundList = dx_ds1OuterAroundList + d1180FirstFace
            elif startPhase == 180.0 and n2 > elementsCountAlongSegment - 1:
                dx_ds1OuterAroundList = dx_ds1OuterAroundList + d1180LastFace
            elif startPhase == 180.0 and elementsCountAlongSegment % 2 == 0 and \
                    n2 == int(elementsCountAlongSegment * 0.5):
                d1Corrected = d1180MidFace

            else:
                for n1 in range(elementsCountAroundHalfHaustrum):
                    v1 = xAround[n1]
                    v2 = xAround[n1 + 1]
                    d1 = d1AtStartOfEachMidFace[n2] if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                    d2 = [v2[c] - v1[c] for c in range(3)]
                    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                    dx_ds1 = [c * arcLengthAround for c in normalize(d1)]
                    dx_ds1OuterAroundList.append(dx_ds1)
                # Account for d1 of node sitting on half haustrum
                d1 = normalize([xAround[elementsCountAroundHalfHaustrum][c] -
                                       xAround[elementsCountAroundHalfHaustrum - 1][c] for c in range(3)])
                dx_ds1 = [c * arcLengthAround for c in d1]
                dx_ds1OuterAroundList.append(dx_ds1)

            if dx_ds1OuterAroundList:
                d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, dx_ds1OuterAroundList,
                                                                      fixStartDerivative=True)
                d1TCEdge = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                               magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
                d1Transition = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                                   magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
                d1Corrected = []
                d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC * 0.5)]
                d1Corrected.append(d1TCEdge)
                d1Corrected.append(d1Transition)
                d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC * 0.5 + 2):]

            xAlongList, d1AlongList, d2AlongList = getFullProfileFromHalfHaustrum(xAround, d1Corrected, d2Around,
                                                                                  tcCount)

            # Calculate xiList for relaxed state
            xHalfSetRelaxed, d1HalfSetRelaxed = \
                createHalfSetIntraHaustralSegment(elementsCountAroundTC, elementsCountAroundHaustrum, tcCount,
                                                  tcWidthSegmentList[n2], radiusSegmentList[n2],
                                                  cornerOuterRadiusFactor, sampleElementOut, haustrumOuterRadiusFactor)
            xRelaxed, d1Relaxed, _ = getFullProfileFromHalfHaustrum(xHalfSetRelaxed, d1HalfSetRelaxed, d2Around,
                                                                    tcCount)
            xiFace, relaxedLengthAroundFace = getXiListFromOuterLengthProfile(xRelaxed, d1Relaxed)
            xiList.append(xiFace)
            relaxedLengthList.append(relaxedLengthAroundFace)

            contractedLengthAroundFace = getXiListFromOuterLengthProfile(xAlongList, d1AlongList)[1]
            contractedWallThickness = relaxedLengthAroundFace * wallThickness / contractedLengthAroundFace
            contractedWallThicknessList.append(contractedWallThickness)

            xFinal = xFinal + xAlongList
            d1Final = d1Final + d1AlongList
            d2Final = d2Final + d2AlongList

        # Create annotation groups
        annotationGroupsAround = []
        for i in range(elementsCountAround):
            annotationGroupsAround.append([])

    return xFinal, d1Final, d2Final, transitElementList, xiList, relaxedLengthList, contractedWallThicknessList, \
           segmentAxis, annotationGroupsAround


def createHalfSetInterHaustralSegment(elementsCountAroundTC, elementsCountAroundHaustrum,
                                      tcCount, tcWidth, radius, cornerOuterRadiusFactor, sampleElementOut):
    """
    Find locations and derivative of nodes in half of an
    inter-haustral segment. Circular profile for segment
    with two or less tenia coli and triangular profile
    with round corners for segment with three tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param tcCount: Number of tenia coli.
    :param tcWidth: Width of tenia coli.
    :param radius: Outer radius of circular profile with two tenia coli,
    radius of circle enclosing triangle for profile with three tenia coli.
    :param cornerOuterRadiusFactor: Roundness of triangular corners of
    inter-haustral septa. Factor is multiplied by outer radius
    to get a radius of curvature at the corners. Only applicable for three tenia
    coli. Set to zero for two tenia coli.
    :param sampleElementOut: Number of sample points used to set up profile
    :return: Node location and derivative on half of a haustrum segment.
    """

    xAround = []
    d1Around = []
    xHalfSetInterHaustra = []
    d1HalfSetInterHaustra = []

    # Set up profile
    if tcCount < 3:  # Circular profile
        xLoop, d1Loop = createCirclePoints([0.0, 0.0, 0.0], [radius, 0.0, 0.0], [0.0, radius, 0.0],
                                           sampleElementOut, startRadians=0.0)

    else:  # tcCount == 3, Triangular profile
        cornerRC = cornerOuterRadiusFactor * radius
        radiansRangeRC = [7 * math.pi / 4, 0.0, math.pi / 4]

        for n1 in range(3):
            radiansAround = n1 * 2.0 * math.pi / 3.0
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            xc = [(radius - cornerRC) * cosRadiansAround, (radius - cornerRC) * sinRadiansAround, 0.0]

            for n in range(3):
                radiansRC = radiansAround + radiansRangeRC[n]
                cosRadiansRC = math.cos(radiansRC)
                sinRadiansRC = math.sin(radiansRC)
                x = [xc[0] + cornerRC * cosRadiansRC, xc[1] + cornerRC * sinRadiansRC, 0.0]
                xAround.append(x)
                d1 = [cornerRC * math.pi / 4.0 * -sinRadiansRC, cornerRC * math.pi / 4.0 * cosRadiansRC, 0.0]
                d1Around.append(d1)

        xSample = xAround[1:9] + [xAround[0], xAround[1]]
        d1Sample = d1Around[1:9] + [d1Around[0], d1Around[1]]
        sx, sd1, _, _, _ = interp.sampleCubicHermiteCurves(xSample, d1Sample, sampleElementOut)
        xLoop = sx[:-1]
        d1Loop = interp.smoothCubicHermiteDerivativesLoop(sx[:-1], sd1[:-1])

    # Calculate arc length
    arcDistance = interp.getCubicHermiteArcLength(xLoop[0], d1Loop[0], xLoop[1], d1Loop[1])
    arcLength = arcDistance * sampleElementOut
    arcStart = 0.0
    numberOfHalves = tcCount * 2
    arcEnd = arcLength / numberOfHalves

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(xLoop, d1Loop, tcWidth, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC, _ = sampleTeniaColi(xLoop, d1Loop, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum, _, _ = sampleHaustrum(xLoop, d1Loop, xTC[-1], d1TC[-1],
                                                 arcLength / numberOfHalves, arcDistanceTCEdge,
                                                 elementsCountAroundHaustrum)

    xHalfSetInterHaustra = xHalfSetInterHaustra + xTC + xHaustrum[1:]
    d1HalfSetInterHaustra = d1HalfSetInterHaustra + d1TC + d1Haustrum[1:]

    return xHalfSetInterHaustra, d1HalfSetInterHaustra


def createHalfSetIntraHaustralSegment(elementsCountAroundTC, elementsCountAroundHaustrum,
                                      tcCount, tcWidth, radius, cornerOuterRadiusFactor, sampleElementOut,
                                      haustrumOuterRadiusFactor):
    """
    Find locations and derivative of nodes in half of an intra-haustral
    segment. Bow-tie profile for segment with two tenia coli and
    clover profile for segment with three tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param tcCount: Number of tenia coli.
    :param tcWidth: Width of tenia coli.
    :param radius: Outer radius of circular inter-haustral profile with two
    tenia coli, radius of circle enclosing triangle for inter-haustral profile
    with three tenia coli.
    :param cornerOuterRadiusFactor: Roundness of triangular corners of
    inter-haustral septa. Factor is multiplied by outer radius
    to get a radius of curvature at the corners. Only applicable for three tenia
    coli. Set to zero for two tenia coli.
    :param sampleElementOut: Number of sample points used to set up profile
    :param haustrumOuterRadiusFactor: Factor is multiplied by outer
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra segment.
    :return: Node location and derivative on half of a haustrum segment.
    """

    nxHaustrum = []
    nd1Haustrum = []
    xHalfSetIntraHaustra = []
    d1HalfSetIntraHaustra = []

    # Set up profile
    cornerRC = cornerOuterRadiusFactor * radius
    haustrumRadius = (haustrumOuterRadiusFactor + 1) * radius
    if tcCount == 2:  # Bow-tie profile
        originRC = (radius * radius - haustrumRadius * haustrumRadius) / (-2.0 * haustrumRadius)
        RC = haustrumRadius - originRC
        thetaStart = math.asin(-originRC / RC)
        thetaEnd = math.pi - thetaStart
        rotOriginRC = [0.0, -originRC, 0.0]

    else:  # tcCount == 3, Clover profile
        xc = [(radius - cornerRC) * math.cos(0.0), (radius - cornerRC) * math.sin(0.0), 0.0]
        pt1 = [xc[0] + cornerRC * math.cos(0.0), xc[1] + cornerRC * math.sin(0.0), 0.0]
        xTC2 = radius * math.cos(2.0 * math.pi / 3.0)
        yTC2 = radius * math.sin(2.0 * math.pi / 3.0)
        originRC = (xTC2 * xTC2 + yTC2 * yTC2 - haustrumRadius * haustrumRadius) / (2 * (-xTC2 - haustrumRadius))
        RC = haustrumRadius - originRC

        # Rotate to find originRC of 1st haustrum
        yTC1 = pt1[1]
        rotOriginRC = [originRC * math.cos(-2.0 / 3.0 * math.pi), originRC * math.sin(-2.0 / 3.0 * math.pi), 0.0]

        thetaStart = math.asin((yTC1 + rotOriginRC[1]) / RC)
        thetaEnd = math.pi - math.asin((yTC2 + rotOriginRC[1]) / RC)

    thetaHalfHaustrum = (thetaEnd + thetaStart) * 0.5
    concaveFactor = 0.15
    thetaConcave = (thetaEnd - thetaStart) * concaveFactor + thetaStart
    thetaCircularStart = (thetaEnd - thetaStart) * (concaveFactor + 0.5) * 0.5 + thetaStart

    thetaSet = [thetaStart, thetaConcave]
    xConcave, _ = getCircleXandD1FromRadians(thetaSet, RC, rotOriginRC)
    d1 = [0.0, xConcave[1][1], 0.0]
    nxHaustrum = nxHaustrum + xConcave
    nd1Haustrum = nd1Haustrum + [d1, d1]
    thetaCircularSet = [thetaCircularStart, thetaHalfHaustrum]
    xCircular, d1Circular = getCircleXandD1FromRadians(thetaCircularSet, RC, rotOriginRC)
    nxHaustrum = nxHaustrum + xCircular
    nd1Haustrum = nd1Haustrum + d1Circular
    smoothd1 = interp.smoothCubicHermiteDerivativesLine(nxHaustrum, nd1Haustrum, fixStartDirection=True,
                                                        fixEndDirection=True)

    # Find max x along path
    sxHaustrum, sd1Haustrum, _, _, _ = interp.sampleCubicHermiteCurves(nxHaustrum, smoothd1, sampleElementOut)
    xVal = []
    for i in range(len(sxHaustrum)):
        xVal.append(sxHaustrum[i][0])
    maxValue = max(xVal)
    maxIndex = xVal.index(maxValue)
    yAtMaxValue = sxHaustrum[maxIndex][1]

    arcLength = 0.0
    for e in range(len(nxHaustrum) - 1):
        arcDistance = interp.getCubicHermiteArcLength(nxHaustrum[e], smoothd1[e], nxHaustrum[e + 1], smoothd1[e + 1])
        arcLength = arcLength + arcDistance

    arcDistanceAtMaxValue = arcLength / sampleElementOut * maxIndex
    arcStart = 0.0 if yAtMaxValue > tcWidth * 0.5 else arcDistanceAtMaxValue
    arcEnd = arcDistanceAtMaxValue if yAtMaxValue > tcWidth * 0.5 else arcLength

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(nxHaustrum, smoothd1, tcWidth, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC, arcLengthPerTC = sampleTeniaColi(nxHaustrum, smoothd1, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum, arcLengthPerHaustrum, arcLengthPerTransition = \
        sampleHaustrum(nxHaustrum, smoothd1, xTC[-1], d1TC[-1], arcLength, arcDistanceTCEdge,
                       elementsCountAroundHaustrum)

    xHalfSetIntraHaustra = xHalfSetIntraHaustra + xTC + xHaustrum[1:]
    d1HalfSetIntraHaustra = d1HalfSetIntraHaustra + d1TC + d1Haustrum[1:]

    return xHalfSetIntraHaustra, d1HalfSetIntraHaustra


def findEdgeOfTeniaColi(nx, nd1, tcWidth, arcStart, arcEnd):
    """
    Locate edge of tenia coli on a cubic hermite interpolated
    curve defined by nodes nx and derivatives nd1.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param tcWidth: Width of tenia coli.
    :param arcStart: Lower limit of arc distance to search for
    edge of tenia coli.
    :param arcEnd: Upper limit of arc distance to search for
    edge of tenia coli.
    :return: arc distance covered by tenia coli.
    """
    xTol = 1.0E-6
    for iter in range(100):
        arcDistance = (arcStart + arcEnd) * 0.5
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        diff = x[1] - tcWidth * 0.5
        if abs(diff) > xTol:
            if diff < 0.0:
                arcStart = arcDistance
            else:
                arcEnd = arcDistance
        else:
            arcDistanceTCEdge = arcDistance
            break
    if iter > 99:
        print('Search for TC boundary - Max iters reached:', iter)

    return arcDistanceTCEdge


def sampleTeniaColi(nx, nd1, arcDistanceTCEdge, elementsCountAroundTC):
    """
    Get systematically spaced points and derivatives over the width of
    tenia coli over a cubic hermite interpolated curves defined by nodes
    nx and derivatives nd1 lying along half of a haustrum.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param arcDistanceTCEdge: Arc distance covered by tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :return: coordinates, derivatives of tenia coli, and arclength of
    each tenia coli element.
    """
    xTC = []
    d1TC = []
    arcDistancePerElementTC = arcDistanceTCEdge / (elementsCountAroundTC * 0.5)
    for e in range(int(elementsCountAroundTC * 0.5) + 1):
        arcDistance = arcDistancePerElementTC * e
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = set_magnitude(d1, arcDistancePerElementTC)
        xTC.append(x)
        d1TC.append(d1Scaled)

    return xTC, d1TC, arcDistancePerElementTC


def sampleHaustrum(nx, nd1, xTCLast, d1TCLast, arcLength, arcDistanceTCEdge,
                   elementsCountAroundHaustrum):
    """
    Get systematically spaced points and derivatives over the width of
    haustrum over a cubic hermite interpolated curves defined by nodes
    nx and derivatives nd1 lying along half of a haustrum.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param xTCLast: Coordinates of node lying on edge of tenia coli.
    :param d1TCLast: Derivative of node lying on edge of tenia coli.
    :param arcLength: Arc length of curve.
    :param arcDistanceTCEdge: Arc distance covered by tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :return: coordinates, derivatives on haustrum.
    :return: arclength of haustrum element and transition elements.
    """
    xHaustrum = []
    d1Haustrum = []
    elementLengths = []
    length = arcLength - arcDistanceTCEdge
    elementsCountOut = int(elementsCountAroundHaustrum * 0.5)
    addLengthStart = 0.5 * magnitude(d1TCLast)
    lengthFractionStart = 0.5
    proportionStart = 1.0
    elementLengthMid = (length - addLengthStart) / (elementsCountOut - 1.0 + proportionStart * lengthFractionStart)
    elementLengthProportionStart = proportionStart * lengthFractionStart * elementLengthMid
    elementLengthProportionEnd = elementLengthMid

    for e in range(elementsCountOut):
        elementLengths.append(elementLengthMid)
    elementLengths[0] = addLengthStart + elementLengthProportionStart
    elementLengths[-1] = 0.0 + elementLengthProportionEnd
    haustrumElementLength = elementLengthMid
    transitionElementLength = elementLengths[0]

    arcDistance = arcDistanceTCEdge
    xHaustrum.append(xTCLast)
    d1Scaled = set_magnitude(d1TCLast, elementLengths[0])
    d1Haustrum.append(d1Scaled)

    for e in range(elementsCountOut):
        arcDistance = arcDistance + elementLengths[e]
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = set_magnitude(d1, elementLengths[e] if e > 0 else elementLengths[e + 1])
        xHaustrum.append(x)
        d1Haustrum.append(d1Scaled)

    return xHaustrum, d1Haustrum, haustrumElementLength, transitionElementLength


def getCircleXandD1FromRadians(thetaSet, radius, origin):
    """
    Gets the coordinates and derivatives along a circular path
    based on an angular range.
    :param thetaSet: Lower and upper limit of theta.
    :param radius: Radius of circle.
    :param origin: Origin of circle.
    :return: coordinates, derivatives on lower and upper
    limit of thetaSet.
    """
    nx = []
    nd1 = []
    dTheta = thetaSet[1] - thetaSet[0]
    for n in range(2):
        theta = thetaSet[n]
        x = [radius * math.cos(theta) - origin[0],
             radius * math.sin(theta) - origin[1],
             0.0]
        d1 = [-radius * math.sin(theta) * dTheta,
              radius * math.cos(theta) * dTheta,
              0.0]
        nx.append(x)
        nd1.append(d1)

    return nx, nd1


def getFullProfileFromHalfHaustrum(xHaustrumHalfSet, d1HaustrumHalfSet,
                                   d2HaustrumHalfSet, tcCount):
    """
    Gets the coordinates and derivatives of the entire profile
    using points from first half of the first sector. The first
    sector starts from the x-axis. The first half set of points
    are reflected across the x-axis followed by rotation to get
    the points in the second half of the first sector. The full
    set of points in the first sector are then rotated to obtain
    points in the other two sectors.
    :param xHaustrumHalfSet: Coordinates of points in first
    half of the first sector.
    :param d1HaustrumHalfSet: Derivatives of points around first
    half of the first sector.
    :param d2HaustrumHalfSet: Derivatives of points along first
    half of the first sector.
    :param tcCount: Number of tenia coli.
    :return: coordinates, derivatives of points over entire profile.
    """
    xHaustrumHalfSet2 = []
    d1HaustrumHalfSet2 = []
    d2HaustrumHalfSet2 = []
    xHaustra = []
    d1Haustra = []
    d2Haustra = []
    rotAng = 2 * math.pi / tcCount

    for n in range(1, len(xHaustrumHalfSet)):
        idx = -n + len(xHaustrumHalfSet) - 1
        x = xHaustrumHalfSet[idx]
        d1 = d1HaustrumHalfSet[idx]
        xReflect = [x[0], -x[1], x[2]]
        d1Reflect = [d1[0], -d1[1], d1[2]]
        xRot = [xReflect[0] * math.cos(rotAng) - xReflect[1] * math.sin(rotAng),
                xReflect[0] * math.sin(rotAng) + xReflect[1] * math.cos(rotAng),
                xReflect[2]]
        d1Rot = [-(d1Reflect[0] * math.cos(rotAng) - d1Reflect[1] * math.sin(rotAng)),
                 -(d1Reflect[0] * math.sin(rotAng) + d1Reflect[1] * math.cos(rotAng)),
                 -d1Reflect[2]]
        d2 = d2HaustrumHalfSet[idx]
        d2Reflect = [d2[0], -d2[1], d2[2]]
        d2Rot = [(d2Reflect[0] * math.cos(rotAng) - d2Reflect[1] * math.sin(rotAng)),
                 (d2Reflect[0] * math.sin(rotAng) + d2Reflect[1] * math.cos(rotAng)),
                 d2Reflect[2]]
        xHaustrumHalfSet2.append(xRot)
        d1HaustrumHalfSet2.append(d1Rot)
        d2HaustrumHalfSet2.append(d2Rot)

    xHaustrum = xHaustrumHalfSet + xHaustrumHalfSet2
    d1Haustrum = d1HaustrumHalfSet + d1HaustrumHalfSet2
    d2Haustrum = d2HaustrumHalfSet + d2HaustrumHalfSet2

    # Rotate to get all 3 sectors
    xHaustra = xHaustra + xHaustrum[:-1]
    d1Haustra = d1Haustra + d1Haustrum[:-1]
    d2Haustra = d2Haustra + d2Haustrum[:-1]

    ang = [2 / 3 * math.pi, -2 / 3 * math.pi] if tcCount == 3 else [math.pi]
    for i in range(tcCount - 1):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(xHaustrum) - 1):
            x = xHaustrum[n]
            d1 = d1Haustrum[n]
            d2 = d2Haustrum[n]
            x = [x[0] * cosRotAng - x[1] * sinRotAng, x[0] * sinRotAng + x[1] * cosRotAng, x[2]]
            xHaustra.append(x)
            dx_ds1 = [d1[0] * cosRotAng - d1[1] * sinRotAng, d1[0] * sinRotAng + d1[1] * cosRotAng, d1[2]]
            d1Haustra.append(dx_ds1)
            dx_ds2 = [d2[0] * cosRotAng - d2[1] * sinRotAng, d2[0] * sinRotAng + d2[1] * cosRotAng, d2[2]]
            d2Haustra.append(dx_ds2)

    return xHaustra, d1Haustra, d2Haustra


def getXiListFromOuterLengthProfile(xOuter, d1Outer):
    """
    Gets a list of xi for flat coordinates calculated
    from outer arclength of elements around a segment (most relaxed state).
    :param xOuter: Coordinates of points on outer surface around segment.
    :param d1Outer: Derivatives of points on outer surface around segment.
    :param segmentAxis: Axis of segment.
    :param wallThickness: Thickness of wall.
    :param transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :return xiList: List containing xi for each point on outer surface extruded by
    wall thickness from inner points.
    :return totalArcLengthOuter: Total arclength around outer surface of elements.
    """

    arcLengthList = []
    for n1 in range(len(xOuter)):
        arcLengthPerElement = \
            interp.getCubicHermiteArcLength(xOuter[n1], d1Outer[n1],
                                            xOuter[(n1 + 1) % len(xOuter)], d1Outer[(n1 + 1) % len(xOuter)])
        arcLengthList.append(arcLengthPerElement)

    # Total arcLength
    totalArcLengthOuter = 0.0
    for n1 in range(len(arcLengthList)):
        totalArcLengthOuter += arcLengthList[n1]

    xiList = [0.0]
    arcDistance = 0
    for n in range(len(arcLengthList)):
        arcDistance = arcDistance + arcLengthList[n]
        xi = arcDistance / totalArcLengthOuter
        xiList.append(xi)

    return xiList, totalArcLengthOuter


def getTeniaColi(region, xList, d1List, d2List, d3List, curvatureList,
                 tcCount, elementsCountAroundTC, elementsCountAroundHaustrum,
                 elementsCountAlong, elementsCountThroughWall,
                 tubeTCWidthList, tcThickness, annotationGroupsAround, closedProximalEnd, isHuman,
                 xProximal=[], d1Proximal=[], d2Proximal=[], d3Proximal=[]):
    """
    Create equally spaced points for tenia coli over the outer
    surface of the colon. Points are sampled from a cubic
    hermite curve running through the left edge of tenia coli
    boundary on the outer surface of the haustra, a midpoint
    lying at tenia coli thickness normal to the midpoint tenia coli
    boundary on the outer surface, and the right edge of tenia coli
    boundary. Function incorporates annotation groups for tenia coli,
    and arranges coordinates, derivatives of tenia coli to follow
    after points on outer surface of the haustra as the mesh extends
    along its length.
    :param xList, d1List, d2List, d3List: Coordinates and derivatives of nodes
    on haustra.
    :param curvatureList: Curvature of points along length.
    :param tcCount: Number of tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlong: Number of elements along colon.
    :param elementsCountThroughWall: Number of elements through wall.
    :param tubeTCWidthList: List of tenia coli width along tube length.
    :param tcThickness: Thickness of tenia coli at its thickest part.
    :param annotationGroupsAround: annotation groups for elements around tube.
    :param closedProximalEnd: True when proximal end of tube is closed.
    :param isHuman: True if making scaffold for human.
    :param xProximal, d1Proximal, d2Proximal, d3Proximal: Optional for passing coordinates and derivatives of points on
    proximal end of colon segment.
    :return: coordinates, derivatives, annotationGroupsAround for colon with tenia coli and local index, coordinates and
    derivatives for nodes around distal end of colon.
    """

    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount
    TCEdgeFactor = 1.5
    xTCArranged = []
    d1TCArranged = []
    d2TCArranged = []
    d3TCArranged = []

    # Create points for tenia coli
    for n2 in range(1 if closedProximalEnd else 0, elementsCountAlong + 1):
        xTCRaw = []
        d1TCRaw = []
        d2TCRaw = []
        d3TCRaw = []
        tcWidth = tubeTCWidthList[n2]
        for N in range(tcCount):
            idxTCMid = elementsCountThroughWall + 1 + (n2 - 1) * \
                       elementsCountAround * (elementsCountThroughWall + 1) + \
                       elementsCountAround * elementsCountThroughWall + \
                       N * (elementsCountAroundTC + elementsCountAroundHaustrum) if closedProximalEnd \
                else n2 * elementsCountAround * (elementsCountThroughWall + 1) + \
                     elementsCountAround * elementsCountThroughWall + \
                     N * (elementsCountAroundTC + elementsCountAroundHaustrum)
            unitNorm = normalize(d3List[idxTCMid])
            xMid = [xList[idxTCMid][i] + unitNorm[i] * tcThickness for i in range(3)]
            d1Mid = d1List[idxTCMid]
            TCStartIdx = idxTCMid - int(elementsCountAroundTC * 0.5) if N > 0 else \
                idxTCMid + tcCount * (elementsCountAroundTC + elementsCountAroundHaustrum) - \
                int(elementsCountAroundTC * 0.5)
            TCEndIdx = idxTCMid + int(elementsCountAroundTC * 0.5)
            if closedProximalEnd:
                tcWidth = magnitude([xList[TCStartIdx][c] - xList[TCEndIdx][c] for c in range(3)])

            v1 = xList[TCStartIdx]
            v2 = xMid
            d1MidScaled = [c * tcWidth * TCEdgeFactor for c in normalize(d1Mid)]
            v3 = xList[TCEndIdx]
            nx = [v1, v2, v3]
            nd1 = [d1List[TCStartIdx], d1MidScaled, d1List[TCEndIdx]]
            sxTC, sd1TC, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)
            xTCRaw = xTCRaw + sxTC[1:-1]
            if elementsCountAroundTC == 2:
                p = [v2[i] - v1[i] for i in range(3)]
                A = dot(unitNorm, p)  # A<0 if v2 is higher than v1
                d1 = [c * tcWidth * 0.5 for c in normalize(d1Mid)] if A < 0 else d1MidScaled
            d1TCRaw = d1TCRaw + sd1TC[1:-1] if elementsCountAroundTC > 2 else d1TCRaw + [d1]
            xTCInnerSet = list(range(TCStartIdx + 1, TCEndIdx)) if N > 0 else \
                list(range(TCStartIdx + 1, TCStartIdx + int(elementsCountAroundTC * 0.5))) + \
                list(range(idxTCMid, idxTCMid + int(elementsCountAroundTC * 0.5)))

            for n in range(elementsCountAroundTC - 1):
                xTCInner = xList[xTCInnerSet[n]]
                xTCOuter = sxTC[n + 1]
                d3 = [xTCOuter[i] - xTCInner[i] for i in range(3)]
                d3TCRaw.append(d3)
                d3List[xTCInnerSet[n]] = d3

                curvature = curvatureList[xTCInnerSet[n]]
                distanceToInnerIdx = magnitude([xTCOuter[i] - xTCInner[i] for i in range(3)])
                factor = 1.0 - curvature * distanceToInnerIdx
                d2 = [factor * c for c in d2List[xTCInnerSet[n]]]
                d2TCRaw.append(d2)

        xTCArranged = xTCArranged + xTCRaw[int(elementsCountAroundTC * 0.5 - 1):] + \
                      xTCRaw[:int(elementsCountAroundTC * 0.5 - 1)]
        d1TCArranged = d1TCArranged + d1TCRaw[int(elementsCountAroundTC * 0.5 - 1):] + \
                       d1TCRaw[:int(elementsCountAroundTC * 0.5 - 1)]
        d2TCArranged = d2TCArranged + d2TCRaw[int(elementsCountAroundTC * 0.5 - 1):] + \
                       d2TCRaw[:int(elementsCountAroundTC * 0.5 - 1)]
        d3TCArranged = d3TCArranged + d3TCRaw[int(elementsCountAroundTC * 0.5 - 1):] + \
                       d3TCRaw[:int(elementsCountAroundTC * 0.5 - 1)]

    x, d1, d2, d3, localIdxDistal, xDistal, d1Distal, d2Distal, d3Distal = \
        combineTeniaColiWithColon(xList, d1List, d2List, d3List, xTCArranged, d1TCArranged,
                                  d2TCArranged, d3TCArranged, (elementsCountAroundTC - 1) * tcCount,
                                  elementsCountAround, elementsCountAlong, elementsCountThroughWall, closedProximalEnd,
                                  xProximal, d1Proximal, d2Proximal, d3Proximal)

    # Update annotation groups
    if tcCount == 2 or not isHuman:
        tcGroup = AnnotationGroup(region, get_colon_term("taenia coli"))
        for i in range(elementsCountAroundTC * tcCount):
            annotationGroupsAround.append([tcGroup])

    elif tcCount == 3 and isHuman:
        tlGroup = AnnotationGroup(region, get_colon_term("taenia libera"))
        tmGroup = AnnotationGroup(region, get_colon_term("taenia mesocolica"))
        toGroup = AnnotationGroup(region, get_colon_term("taenia omentalis"))
        annotationGroupAround = [[toGroup], [tlGroup], [tmGroup], [toGroup]]
        elementsCountAroundGroups = [elementsCountAroundTC // 2, elementsCountAroundTC, elementsCountAroundTC,
                                     elementsCountAroundTC // 2]

        for i in range(len(elementsCountAroundGroups)):
            elementsCount = elementsCountAroundGroups[i]
            for n in range(elementsCount):
                annotationGroupsAround.append(annotationGroupAround[i])

    return x, d1, d2, d3, annotationGroupsAround, localIdxDistal, xDistal, d1Distal, d2Distal, d3Distal


def combineTeniaColiWithColon(xList, d1List, d2List, d3List, xTC, d1TC, d2TC, d3TC, nodesCountAroundTC,
                              elementsCountAround, elementsCountAlong, elementsCountThroughWall, closedProximalEnd,
                              xProximal=[], d1Proximal=[], d2Proximal=[], d3Proximal=[]):
    """
    Arranges coordinates and derivatives around inner surface to outer surface, followed by tenia coli points before
    extending along length of colon.
    :param xList, d1List, d2List, d3List: coordinates and derivatives of colon.
    :param xTC, d1TC, d2TC, d3TC: coordinates and derivatives of tenia coli.
    :param nodesCountAroundTC: Number of nodes around tenia coli.
    :param elementsCountAround: Number of elements around colon.
    :param elementsCountAlong: Number of elements along colon.
    :param elementsCountThroughWall: Number of elements through wall.
    :param closedProximalEnd: True when proximal end of tube is closed.
    :param xProximal, d1Proximal, d2Proximal, d3Proximal: Optional for passing coordinates and derivatives of points on
    proximal end of colon segment.
    :return: reordered coordinates and derivatives and local index, coordinates and derivatives of nodes on distal end.
    """
    x = []
    d1 = []
    d2 = []
    d3 = []
    count = 0
    localIdxDistal = []
    xDistal = []
    d1Distal = []
    d2Distal = []
    d3Distal = []

    # Add tenia coli points to coordinates list
    for n2 in range(elementsCountAlong + 1):
        for n3 in range(elementsCountThroughWall + 1):
            xDistalAround = []
            d1DistalAround = []
            d2DistalAround = []
            d3DistalAround = []
            localIdxDistalAround = []

            if closedProximalEnd and n2 == 0:
                x.append(xList[n3])
                d1.append(d1List[n3])
                d2.append(d2List[n3])
                if d3List:
                    d3.append(d3List[n3])
                count += 1
            else:
                for n1 in range(elementsCountAround):
                    # Append colon wall coordinates from inside to outside wall
                    n = (elementsCountThroughWall + 1) + \
                        (n2 - 1) * elementsCountAround * (elementsCountThroughWall + 1) + \
                        n3 * elementsCountAround + n1 if closedProximalEnd \
                        else n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1

                    if xProximal and n2 == 0 and not closedProximalEnd:
                        x.append(xProximal[n3][n1])
                        d1.append(d1Proximal[n3][n1])
                        d2.append(d2Proximal[n3][n1])
                        if d3List:
                            d3.append(d3Proximal[n3][n1])
                    else:
                        x.append(xList[n])
                        d1.append(d1List[n])
                        d2.append(d2List[n])
                        if d3List:
                            d3.append(d3List[n])

                        if n2 == elementsCountAlong:
                            localIdxDistalAround.append(count)
                            xDistalAround.append(xList[n])
                            d1DistalAround.append(d1List[n])
                            d2DistalAround.append(d2List[n])
                            if d3List:
                                d3DistalAround.append(d3List[n])
                        count += 1

                if n2 == elementsCountAlong:
                    xDistal.append(xDistalAround)
                    d1Distal.append(d1DistalAround)
                    d2Distal.append(d2DistalAround)
                    d3Distal.append(d3DistalAround)
                    localIdxDistal.append(localIdxDistalAround)

        # Append tenia coli coordinates
        if not (closedProximalEnd and n2 == 0):
            for nTC in range(nodesCountAroundTC):
                if xProximal and n2 == 0:
                    x.append(xProximal[-1][elementsCountAround + nTC])
                    d1.append(d1Proximal[-1][elementsCountAround + nTC])
                    d2.append(d2Proximal[-1][elementsCountAround + nTC])
                    if d3List:
                        d3.append(d3Proximal[-1][elementsCountAround + nTC])
                else:
                    nTCCount = (n2 - 1 if closedProximalEnd else n2) * nodesCountAroundTC + nTC
                    x.append(xTC[nTCCount])
                    d1.append(d1TC[nTCCount])
                    d2.append(d2TC[nTCCount])
                    if d3TC:
                        d3.append(d3TC[nTCCount])
                    if n2 == elementsCountAlong:
                        localIdxDistal[-1].append(count)
                        xDistal[-1].append(xTC[nTCCount])
                        d1Distal[-1].append(d1TC[nTCCount])
                        d2Distal[-1].append(d2TC[nTCCount])
                        if d3TC:
                            d3Distal[-1].append(d3TC[nTCCount])
                    count += 1

    return x, d1, d2, d3, localIdxDistal, xDistal, d1Distal, d2Distal, d3Distal


def createFlatCoordinatesTeniaColi(xiList, relaxedLengthList,
                                   totalLengthAlong, wallThickness, relativeThicknessList, tcCount, tcThickness,
                                   elementsCountAroundTC, elementsCountAroundHaustrum,
                                   elementsCountAlong, elementsCountThroughWall, transitElementList, closedProximalEnd):
    """
    Calculates flat coordinates for a colon scaffold with tenia coli when it is opened
    up into a flat preparation.
    :param xiList: List containing xi for each point around the outer surface of
    colon in its most relaxed state.
    :param relaxedLengthList: List of total arclength around the outer surface in
    its most relaxed state for each element along.
    :param totalLengthAlong: Total length along colon.
    :param wallThickness: Thickness of wall.
    :param relativeThicknessList: Relative thickness of each element through wall.
    :param tcCount: Number of tenia coli.
    :param tcThickness: Thickness of tenia coli at its thickest region.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlong: Number of elements along colon.
    :param elementsCountThroughWall: Number of elements through wall.
    :param transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :param closedProximalEnd: True when proximal end of tube is closed.
    :return: coordinates and derivatives of flat coordinates fields.
    """

    # Calculate flat coordinates
    factor = 3.0 if tcCount == 3 else 2.0
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount

    # Find flat coordinates for colon
    xFlatColon, d1FlatColon, d2FlatColon = tubemesh.createFlatCoordinates(xiList, relaxedLengthList, totalLengthAlong,
                                                                          wallThickness, relativeThicknessList,
                                                                          elementsCountAround, elementsCountAlong,
                                                                          elementsCountThroughWall, transitElementList)

    # Find flat coordinates for tenia coli
    xFlatListTC = []
    d1FlatListTC = []
    d2FlatListTC = []

    for n2 in range(elementsCountAlong + 1):
        xiFace = xiList[n2]
        relaxedLength = relaxedLengthList[n2]
        xPad = (relaxedLengthList[0] - relaxedLength) * 0.5
        for N in range(tcCount + 1):
            idxTCMid = N * (elementsCountAroundTC + elementsCountAroundHaustrum)
            TCStartIdx = idxTCMid - int(elementsCountAroundTC * 0.5)
            TCEndIdx = idxTCMid + int(elementsCountAroundTC * 0.5)
            dTC = (xiFace[idxTCMid] - xiFace[TCStartIdx]) * relaxedLength if N > 0 else \
                (xiFace[TCEndIdx] - xiFace[idxTCMid]) * relaxedLength
            v1 = [xiFace[TCStartIdx] * relaxedLength, 0.0, wallThickness] if N > 0 else \
                [-dTC, 0.0, wallThickness]
            v2 = [xiFace[idxTCMid] * relaxedLength, 0.0, wallThickness + tcThickness]
            v3 = [xiFace[TCEndIdx] * relaxedLength, 0.0, wallThickness] if N < tcCount else \
                [relaxedLength + dTC, 0.0, wallThickness]
            d1 = d3 = [dTC, 0.0, 0.0]
            d2 = [c * factor for c in [dTC, 0.0, 0.0]]
            nx = [v1, v2, v3]
            nd1 = [d1, d2, d3]
            sx, sd1, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)
            if 0 < N < tcCount:
                w = sx[1:-1]
                dw = sd1[1:-1] if elementsCountAroundTC > 2 else nd1[1:-1]
            elif N == 0:
                w = sx[int(elementsCountAroundTC * 0.5):-1]
                dw = sd1[int(elementsCountAroundTC * 0.5):-1] if elementsCountAroundTC > 2 else \
                    nd1[int(elementsCountAroundTC * 0.5):-1]
            else:
                w = sx[1:int(elementsCountAroundTC * 0.5) + 1]
                dw = sd1[1:int(elementsCountAroundTC * 0.5) + 1] if elementsCountAroundTC > 2 else \
                    nd1[1:int(elementsCountAroundTC * 0.5) + 1]
            for n in range(len(w)):
                x = [xPad + w[n][0],
                     totalLengthAlong / elementsCountAlong * n2,
                     w[n][2]]
                xFlatListTC.append(x)
                d1FlatListTC.append(dw[n])

    for n2 in range(elementsCountAlong):
        for n1 in range((elementsCountAroundTC - 1) * tcCount + 1):
            nodeIdx = n2 * ((elementsCountAroundTC - 1) * tcCount + 1) + n1
            nodeNextElementAlong = nodeIdx + ((elementsCountAroundTC - 1) * tcCount + 1)
            v1 = xFlatListTC[nodeNextElementAlong]
            v2 = xFlatListTC[nodeIdx]
            d1 = d2 = [v1[i] - v2[i] for i in range(3)]
            arclength = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
            d2Flat = set_magnitude(d1, arclength)
            d2FlatListTC.append(d2Flat)
    d2FlatListTC = d2FlatListTC + d2FlatListTC[-((elementsCountAroundTC - 1) * tcCount + 1):]

    xFlat, d1Flat, d2Flat = combineTeniaColiWithColon(xFlatColon, d1FlatColon, d2FlatColon, [],
                                                         xFlatListTC, d1FlatListTC, d2FlatListTC, [],
                                                         (elementsCountAroundTC - 1) * tcCount + 1,
                                                         elementsCountAround + 1, elementsCountAlong,
                                                         elementsCountThroughWall, closedProximalEnd)[0:3]

    return xFlat, d1Flat, d2Flat


def createColonCoordinatesTeniaColi(xiList, relativeThicknessList, lengthToDiameterRatio, wallThicknessToDiameterRatio,
                                    teniaColiThicknessToDiameterRatio, tcCount, elementsCountAroundTC,
                                    elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                                    transitElementList, closedProximalEnd):
    """
    Calculates organ coordinates for a colon scaffold with tenia coli. Organ coordinate takes the form of a cylinder
    with unit inner diameter, length of lengthToDiameterRatio, wall thickness of wallThicknessToDiameterRatio, with
    tenia coli of teniaColiThicknessToDiameterRatio running along its length.
    :param xiList: List containing xi for each point around the outer surface of colon in its most relaxed state.
    :param relativeThicknessList: Relative thickness for each element through wall for colon coordinates.
    :param lengthToDiameterRatio: Ratio of total length along organ to inner diameter of organ
    :param wallThicknessToDiameterRatio: Ratio of wall thickness to inner diameter of organ.
    :param teniaColiThicknessToDiameterRatio: Ratio of tenia coli thickness to inner diameter of organ.
    :param tcCount: Number of tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlong: Number of elements along colon.
    :param elementsCountThroughWall: Number of elements through wall.
    :param transitElementList: stores true if element around is an element that transits from tenia coli to haustrum.
    :param closedProximalEnd: True when proximal end of tube is closed.
    :return: coordinates and derivatives of colon coordinates field.
    """
    # Calculate organ coordinates
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount

    # Find organ coordinates for colon
    xColon, d1Colon, d2Colon = tubemesh.createOrganCoordinates(xiList, relativeThicknessList, lengthToDiameterRatio,
                                                               wallThicknessToDiameterRatio, elementsCountAround,
                                                               elementsCountAlong, elementsCountThroughWall,
                                                               transitElementList)

    # Find organ coordinates for tenia coli
    xTC = []
    d1TC = []
    d2 = [0.0, lengthToDiameterRatio / elementsCountAlong, 0.0]

    for n2 in range(elementsCountAlong + 1):
        faceStartIdx = elementsCountAround * (elementsCountThroughWall + 1) * n2 + \
                       elementsCountAround * elementsCountThroughWall
        xTCRaw = []
        d1TCRaw = []
        for N in range(tcCount):
            TCMidIdx = faceStartIdx + N * (elementsCountAroundTC + elementsCountAroundHaustrum)
            TCStartIdx = TCMidIdx - int(elementsCountAroundTC * 0.5) if N > 0 else \
                TCMidIdx + tcCount * (elementsCountAroundTC + elementsCountAroundHaustrum) - int(
                    elementsCountAroundTC * 0.5)
            TCEndIdx = TCMidIdx + int(elementsCountAroundTC * 0.5)
            v1 = xColon[TCStartIdx]
            norm = set_magnitude(
                cross(normalize(d1Colon[TCMidIdx]), normalize(d2Colon[TCMidIdx])),
                teniaColiThicknessToDiameterRatio)
            v2 = [xColon[TCMidIdx][c] + norm[c] for c in range(3)]
            v3 = xColon[TCEndIdx]
            nx = [v1, v2, v3]
            nd1 = [d1Colon[TCStartIdx],
                   set_magnitude(d1Colon[TCMidIdx], 2.0 * magnitude(d1Colon[TCMidIdx])),
                   d1Colon[TCEndIdx]]
            sx, sd1 = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)[0:2]
            xTCRaw += sx[1:-1]
            d1TCRaw += sd1[1:-1]

        xTC += xTCRaw[int(elementsCountAroundTC * 0.5 - 1):] + xTCRaw[:int(elementsCountAroundTC * 0.5 - 1)]
        d1TC += d1TCRaw[int(elementsCountAroundTC * 0.5 - 1):] + d1TCRaw[:int(elementsCountAroundTC * 0.5 - 1)]

    d2TC = [d2 for n in range(len(xTC))]

    xColon, d1Colon, d2Colon = combineTeniaColiWithColon(xColon, d1Colon, d2Colon, [], xTC, d1TC, d2TC, [],
                                                            (elementsCountAroundTC - 1) * tcCount,
                                                            elementsCountAround, elementsCountAlong,
                                                            elementsCountThroughWall, closedProximalEnd)[0:3]

    return xColon, d1Colon, d2Colon


def createNodesAndElementsTeniaColi(region,
                                    x, d1, d2, d3,
                                    xFlat, d1Flat, d2Flat,
                                    xOrgan, d1Organ, d2Organ, organCoordinateFieldName,
                                    elementsCountAroundTC, elementsCountAroundHaustrum,
                                    elementsCountAlong, elementsCountThroughWall, tcCount,
                                    annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
                                    firstNodeIdentifier, firstElementIdentifier,
                                    useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd,
                                    localIdxDistal=[], nodeIdProximal=[]):
    """
    Create nodes and elements for the coordinates and flat coordinates fields.
    Note that flat coordinates not implemented for closedProximalEnd yet.
    :param x, d1, d2, d3: coordinates and derivatives of coordinates field.
    :param xFlat, d1Flat, d2Flat: coordinates and derivatives of
    flat coordinates field.
    :param xOrgan, d1Organ, d2Organ, d3Organ, organCoordinateFieldName: coordinates,
    derivatives and name of organ coordinates field.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlong: Number of elements along colon.
    :param elementsCountThroughWall: Number of elements through wall.
    :param tcCount: Number of tenia coli.
    :param annotationGroupsAround: Annotation groups of elements around.
    :param annotationGroupsAlong: Annotation groups of elements along.
    :param annotationGroupsThroughWall: Annotation groups of elements through wall.
    :param firstNodeIdentifier, firstElementIdentifier: first node and element
    identifier to use.
    :param useCubicHermiteThroughWall: use linear when false.
    :param useCrossDerivatives: use cross derivatives when true.
    :param closedProximalEnd: True when proximal end of tube is closed.
    :param localIdxDistal: Local identifiers derived for distal nodes of the tube.
    :param nodeIdProximal: Node identifiers to use as proximal end of tube.
    :return nodeIdentifier, elementIdentifier, allAnnotationGroups, nodes on distal end of scaffold.
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount
    startNode = firstNodeIdentifier

    # Create coordinates field
    zero = [0.0, 0.0, 0.0]
    nodesDistal = []

    for i in range(len(localIdxDistal)):
        nodesDistal.append([firstNodeIdentifier + c for c in localIdxDistal[i]])

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    if useCubicHermiteThroughWall:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

    mesh = fm.findMeshByDimension(3)

    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    # Tenia coli edge elements
    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    eft1 = eftfactory.createEftWedgeXi1One()
    elementtemplate1.defineField(coordinates, -1, eft1)

    elementtemplate2 = mesh.createElementtemplate()
    elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    eft2 = eftfactory.createEftWedgeXi1Zero()
    elementtemplate2.defineField(coordinates, -1, eft2)

    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftFlat3 = bicubichermitelinear.createEftBasic()
    eftFlat4 = bicubichermitelinear.createEftOpenTube()
    eftFlat5 = bicubichermitelinear.createEftWedgeXi1One()
    eftFlat6 = bicubichermitelinear.createEftWedgeXi1Zero()
    eftFlat7 = bicubichermitelinear.createEftWedgeXi1ZeroOpenTube()

    if xFlat:
        # Create flat coordinates field
        flatCoordinates = findOrCreateFieldCoordinates(fm, name="flat coordinates")
        flatNodetemplate1 = nodes.createNodetemplate()
        flatNodetemplate1.defineField(flatCoordinates)
        flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        flatNodetemplate2 = nodes.createNodetemplate()
        flatNodetemplate2.defineField(flatCoordinates)
        flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
        flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
        flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
        if useCrossDerivatives:
            flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)

        flatElementtemplate1 = mesh.createElementtemplate()
        flatElementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate1.defineField(flatCoordinates, -1, eftFlat3)

        flatElementtemplate2 = mesh.createElementtemplate()
        flatElementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate2.defineField(flatCoordinates, -1, eftFlat4)

        flatElementtemplate3 = mesh.createElementtemplate()
        flatElementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate3.defineField(flatCoordinates, -1, eftFlat5)

        flatElementtemplate4 = mesh.createElementtemplate()
        flatElementtemplate4.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate4.defineField(flatCoordinates, -1, eftFlat6)

        flatElementtemplate5 = mesh.createElementtemplate()
        flatElementtemplate5.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate5.defineField(flatCoordinates, -1, eftFlat7)

    if xOrgan:
        # Organ coordinates field
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eftOrgan = bicubichermitelinear.createEftBasic()

        organCoordinates = findOrCreateFieldCoordinates(fm, name=organCoordinateFieldName)
        organNodetemplate = nodes.createNodetemplate()
        organNodetemplate.defineField(organCoordinates)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        organElementtemplate = mesh.createElementtemplate()
        organElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        organElementtemplate.defineField(organCoordinates, -1, eftOrgan)

        # Tenia coli edge elements
        organElementtemplate1 = mesh.createElementtemplate()
        organElementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        eftOrgan1 = bicubichermitelinear.createEftWedgeXi1One()
        organElementtemplate1.defineField(organCoordinates, -1, eftOrgan1)

        organElementtemplate2 = mesh.createElementtemplate()
        organElementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        eftOrgan2 = bicubichermitelinear.createEftWedgeXi1Zero()
        organElementtemplate2.defineField(organCoordinates, -1, eftOrgan2)

    # create nodes for coordinates field
    if nodeIdProximal:
        proximalNodesOffset = 0
        for m in range(len(nodeIdProximal)):
            proximalNodesOffset += len(nodeIdProximal[m])
        # proximalNodesOffset = len(nodeIdProximal) * len(nodeIdProximal[0])
        nodeList = []
        newNodeList = []

    if nodeIdProximal:
        for n3 in range(len(nodeIdProximal)):
            for n1 in range(len(nodeIdProximal[n3])):
                nodeList.append(nodeIdentifier)
                newNodeList.append(nodeIdProximal[n3][n1])
                nodeIdentifier = nodeIdentifier + 1

    for n in range(proximalNodesOffset if nodeIdProximal else 0, len(x)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes for flat coordinates field
    if xFlat:
        nodeIdentifier = firstNodeIdentifier
        for n2 in range(elementsCountAlong + 1):
            for n3 in range(elementsCountThroughWall + 1):
                for n1 in range(elementsCountAround):
                    i = n2 * (elementsCountAround + 1) * (elementsCountThroughWall + 1) + \
                        (elementsCountAround + 1) * n3 + n1 + n2 * ((elementsCountAroundTC - 1) * tcCount + 1)
                    node = nodes.findNodeByIdentifier(nodeIdentifier)
                    node.merge(flatNodetemplate2 if n1 == 0 else flatNodetemplate1)
                    cache.setNode(node)
                    # print('NodeIdentifier', nodeIdentifier, 'version 1', xFlatList[i])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlat[i])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Flat[i])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Flat[i])
                    if useCrossDerivatives:
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    if n1 == 0:
                        # print('NodeIdentifier', nodeIdentifier, 'version 2', xFlatList[i+elementsCountAround])
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2,
                                                          xFlat[i + elementsCountAround])
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2,
                                                          d1Flat[i + elementsCountAround])
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2,
                                                          d2Flat[i + elementsCountAround])
                        if useCrossDerivatives:
                            flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                    nodeIdentifier = nodeIdentifier + 1

            # Create flat coordinates nodes for tenia coli
            for nTC in range((elementsCountAroundTC - 1) * tcCount):
                j = i + 2 + nTC
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                node.merge(flatNodetemplate2 if nTC == 0 else flatNodetemplate1)
                cache.setNode(node)
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlat[j])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Flat[j])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Flat[j])
                if useCrossDerivatives:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                if nTC == 0:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2,
                                                      xFlat[j + (elementsCountAroundTC - 1) * tcCount])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2,
                                                      d1Flat[j + (elementsCountAroundTC - 1) * tcCount])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2,
                                                      d2Flat[j + (elementsCountAroundTC - 1) * tcCount])
                    if useCrossDerivatives:
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                nodeIdentifier = nodeIdentifier + 1

    # Create nodes for organ coordinates field
    if xOrgan:
        nodeIdentifier = firstNodeIdentifier
        for n in range(len(xOrgan)):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(organNodetemplate)
            cache.setNode(node)
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOrgan[n])
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Organ[n])
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Organ[n])
            if useCrossDerivatives:
                organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

    # Create elements
    allAnnotationGroups = []

    for group in annotationGroupsThroughWall:
        longitudinalMuscle = findAnnotationGroupByName(group, ("longitudinal muscle layer of cecum" if closedProximalEnd
                                                               else "longitudinal muscle layer of colon"))

    if longitudinalMuscle:
        longitudinalMuscleGroup = \
        AnnotationGroup(region, (get_cecum_term("longitudinal muscle layer of cecum") if closedProximalEnd
                                 else get_colon_term("longitudinal muscle layer of colon")))

    if closedProximalEnd:
        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        radiansPerElementAround = math.pi * 2.0 / elementsCountAround

        elementtemplate4 = mesh.createElementtemplate()
        elementtemplate4.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementtemplate5 = mesh.createElementtemplate()
        elementtemplate5.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementtemplate6 = mesh.createElementtemplate()
        elementtemplate6.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Create apex
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1) % elementsCountAround
                eft3 = eftfactory.createEftShellPoleBottom(va * 100, vb * 100)
                elementtemplate3.defineField(coordinates, -1, eft3)
                element = mesh.createElement(elementIdentifier, elementtemplate3)
                bni1 = e3 + 1
                bni2 = elementsCountThroughWall + 1 + elementsCountAround * e3 + e1 + 1
                bni3 = elementsCountThroughWall + 1 + elementsCountAround * e3 + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni1 + 1, bni2 + elementsCountAround, bni3 + elementsCountAround]
                nodeIdentifiers = [c + startNode - 1 for c in nodeIdentifiers]
                element.setNodesByIdentifier(eft3, nodeIdentifiers)
                # set general linear map coefficients
                radiansAround = e1 * radiansPerElementAround
                radiansAroundNext = ((e1 + 1) % elementsCountAround) * radiansPerElementAround
                scalefactors = [
                    -1.0,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                ]
                element.setScaleFactors(eft3, scalefactors)
                elementIdentifier = elementIdentifier + 1
                annotationGroups = annotationGroupsAround[e1] + annotationGroupsAlong[0] + \
                                   annotationGroupsThroughWall[e3]
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

        # Create tenia coli
        for eTC in range(int(elementsCountAroundTC * 0.5)):
            va = eTC
            vb = eTC + 1
            # set general linear map coefficients
            radiansAround = va * radiansPerElementAround
            radiansAroundNext = vb * radiansPerElementAround
            scalefactors = [-1.0,
                            math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                            math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]

            bni21 = elementsCountThroughWall + 1 + startNode - 1
            bni22 = bni21 + elementsCountAround * elementsCountThroughWall + 1 + eTC
            bni23 = bni22 + 1
            bni31 = bni22 + elementsCountAround

            tetrahedralElement = (eTC == int(elementsCountAroundTC * 0.5) - 1)
            if tetrahedralElement:
                eft4 = eftfactory.createEftTetrahedronXi1One(va * 10000, vb * 10000)
                elementtemplate4.defineField(coordinates, -1, eft4)
                nodeIdentifiers = [bni21, bni22, bni23, bni31]
            else:
                eft6 = eftfactory.createEftPyramidBottomSimple(va * 10000, vb * 10000)
                elementtemplate6.defineField(coordinates, -1, eft6)
                nodeIdentifiers = [bni21, bni22, bni23, bni31, bni31 + 1]

            element = \
                mesh.createElement(elementIdentifier, elementtemplate4 if tetrahedralElement else elementtemplate6)
            element.setNodesByIdentifier(eft4 if tetrahedralElement else eft6, nodeIdentifiers)
            element.setScaleFactors(eft4 if tetrahedralElement else eft6, scalefactors)
            elementIdentifier = elementIdentifier + 1
            annotationGroups = annotationGroupsAround[elementsCountAround + eTC] + annotationGroupsAlong[0] + \
                               ([longitudinalMuscleGroup] if longitudinalMuscle else [])
            if annotationGroups:
                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                for annotationGroup in annotationGroups:
                    meshGroup = annotationGroup.getMeshGroup(mesh)
                    meshGroup.addElement(element)

        for N in range(tcCount - 1):
            for eTC in range(elementsCountAroundTC):
                va = int(elementsCountAroundTC * 0.5) + elementsCountAroundHaustrum + \
                     N * (elementsCountAroundHaustrum + elementsCountAroundTC) + eTC
                vb = va + 1
                # set general linear map coefficients
                radiansAround = va * radiansPerElementAround
                radiansAroundNext = vb * radiansPerElementAround
                scalefactors = [
                    -1.0,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]

                bni21 = elementsCountThroughWall + 1 + startNode - 1
                bni22 = bni21 + elementsCountAround * elementsCountThroughWall + int(elementsCountAroundTC * 0.5) + \
                        elementsCountAroundHaustrum + N * (elementsCountAroundTC + elementsCountAroundHaustrum) + \
                        eTC + 1
                bni23 = bni22 + 1
                bni31 = bni21 + elementsCountAround * (elementsCountThroughWall + 1) + \
                        int(elementsCountAroundTC * 0.5) + N * (elementsCountAroundTC - 1) + eTC + 1

                if eTC == 0:
                    eft5 = eftfactory.createEftTetrahedronXi1Zero(va * 10000, vb * 10000)
                    elementtemplate5.defineField(coordinates, -1, eft5)
                    nodeIdentifiers = [bni21, bni22, bni23, bni31]
                    element = mesh.createElement(elementIdentifier, elementtemplate5)
                    element.setNodesByIdentifier(eft5, nodeIdentifiers)
                    element.setScaleFactors(eft5, scalefactors)

                elif 0 < eTC < elementsCountAroundTC - 1:
                    eft6 = eftfactory.createEftPyramidBottomSimple(va * 10000, vb * 10000)
                    elementtemplate6.defineField(coordinates, -1, eft6)
                    nodeIdentifiers = [bni21, bni22, bni23, bni31 - 1, bni31]
                    element = mesh.createElement(elementIdentifier, elementtemplate6)
                    element.setNodesByIdentifier(eft6, nodeIdentifiers)
                    element.setScaleFactors(eft6, scalefactors)

                else:
                    eft4 = eftfactory.createEftTetrahedronXi1One(va * 10000, vb * 10000)
                    elementtemplate4.defineField(coordinates, -1, eft4)
                    nodeIdentifiers = [bni21, bni22, bni23, bni31 - 1]
                    element = mesh.createElement(elementIdentifier, elementtemplate4)
                    element.setNodesByIdentifier(eft4, nodeIdentifiers)
                    element.setScaleFactors(eft4, scalefactors)

                elementIdentifier = elementIdentifier + 1
                annotationGroups = \
                    annotationGroupsAround[elementsCountAround + int(elementsCountAroundTC * 0.5) +
                                           N * elementsCountAroundTC + eTC] + annotationGroupsAlong[0] + \
                    ([longitudinalMuscleGroup] if longitudinalMuscle else [])
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

        for eTC in range(int(elementsCountAroundTC * 0.5)):
            va = int(elementsCountAround - elementsCountAroundTC * 0.5) + eTC
            vb = (va + 1) % elementsCountAround
            # set general linear map coefficients
            radiansAround = va * radiansPerElementAround
            radiansAroundNext = vb * radiansPerElementAround
            scalefactors = [-1.0,
                            math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                            math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]

            bni21 = elementsCountThroughWall + 1 + startNode - 1
            bni22 = bni21 + elementsCountAround * (elementsCountThroughWall + 1) - \
                    int(elementsCountAroundTC * 0.5) + eTC + 1
            if elementsCountAroundTC > 2:
                bni23 = bni21 + elementsCountAround * elementsCountThroughWall + 1 \
                    if eTC == int(elementsCountAroundTC * 0.5 - 1) \
                    else bni22 + 1
                bni31 = bni21 + elementsCountAround * (elementsCountThroughWall + 1) + \
                        int(tcCount - 1) * (elementsCountAroundTC - 1) + \
                        int(elementsCountAroundTC * 0.5) + eTC + (0 if eTC > 0 else 1)
            else:
                bni23 = bni21 + elementsCountAround * elementsCountThroughWall + eTC + 1
                bni31 = bni21 + elementsCountAround * (elementsCountThroughWall + 1) + 1
            bni32 = bni21 + elementsCountAround * (elementsCountThroughWall + 1) + 1 \
                if eTC == int(elementsCountAroundTC * 0.5 - 1) else bni31 + 1

            tetrahedralElement = (eTC == 0)
            if tetrahedralElement:
                eft5 = eftfactory.createEftTetrahedronXi1Zero(va * 10000, vb * 10000)
                elementtemplate5.defineField(coordinates, -1, eft5)
                nodeIdentifiers = [bni21, bni22, bni23, bni31]
            else:
                eft6 = eftfactory.createEftPyramidBottomSimple(va * 10000, vb * 10000)
                elementtemplate6.defineField(coordinates, -1, eft6)
                nodeIdentifiers = [bni21, bni22, bni23, bni31, bni32]

            element = \
                mesh.createElement(elementIdentifier, elementtemplate5 if tetrahedralElement else elementtemplate6)
            element.setNodesByIdentifier(eft5 if tetrahedralElement else eft6, nodeIdentifiers)
            element.setScaleFactors(eft5 if tetrahedralElement else eft6, scalefactors)
            elementIdentifier = elementIdentifier + 1
            annotationGroups = annotationGroupsAround[elementsCountAround + eTC] + annotationGroupsAlong[0] + \
                               ([longitudinalMuscleGroup] if longitudinalMuscle else [])
            if annotationGroups:
                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                for annotationGroup in annotationGroups:
                    meshGroup = annotationGroup.getMeshGroup(mesh)
                    meshGroup.addElement(element)

    # create regular elements
    now = elementsCountAround * (elementsCountThroughWall + 1)
    tcOffset = (elementsCountAroundTC - 1) * tcCount
    for e2 in range(1 if closedProximalEnd else 0, elementsCountAlong):
        tcOffset1 = (e2 - 1 if closedProximalEnd else e2) * (elementsCountAroundTC - 1) * tcCount
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                if closedProximalEnd:
                    bni11 = (e2 - 1) * now + e3 * elementsCountAround + e1 + 1 + (elementsCountThroughWall + 1) + \
                            tcOffset1
                    bni12 = (e2 - 1) * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1) + tcOffset1
                    bni21 = (e2 - 1) * now + (e3 + 1) * elementsCountAround + e1 + 1 + \
                            (elementsCountThroughWall + 1) + tcOffset1
                    bni22 = (e2 - 1) * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1) + tcOffset1
                else:
                    bni11 = e2 * now + e3 * elementsCountAround + e1 + 1 + tcOffset1
                    bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + tcOffset1
                    bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + 1 + tcOffset1
                    bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + tcOffset1
                nodeIdentifiers = [bni11, bni12, bni11 + now + tcOffset, bni12 + now + tcOffset,
                                   bni21, bni22, bni21 + now + tcOffset, bni22 + now + tcOffset]
                nodeIdentifiers = [c + startNode - 1 for c in nodeIdentifiers]
                if e2 == 0 and nodeIdProximal:
                    for m in range(len(nodeIdentifiers)):
                        if nodeIdentifiers[m] in nodeList:
                            idx = nodeList.index(nodeIdentifiers[m])
                            nodeIdentifiers[m] = newNodeList[idx]
                onOpening = e1 > elementsCountAround - 2
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if xFlat:
                    element.merge(flatElementtemplate2 if onOpening else flatElementtemplate1)
                    element.setNodesByIdentifier(eftFlat4 if onOpening else eftFlat3, nodeIdentifiers)
                if xOrgan:
                    element.merge(organElementtemplate)
                    element.setNodesByIdentifier(eftOrgan, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1
                annotationGroups = annotationGroupsAround[e1] + annotationGroupsAlong[e2] + \
                                   annotationGroupsThroughWall[e3]
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

        # Add elements for tenia coli
        for eTC in range(int(elementsCountAroundTC * 0.5)):
            if closedProximalEnd:
                bni21 = elementsCountThroughWall + 1 + (e2 - 1) * now + \
                        elementsCountThroughWall * elementsCountAround + eTC + 1 + tcOffset1
                bni22 = elementsCountThroughWall + 1 + (e2 - 1) * now + \
                        elementsCountThroughWall * elementsCountAround + eTC + 2 + tcOffset1
                bni31 = elementsCountThroughWall + 1 + e2 * now + eTC + 1 + tcOffset1
                bni32 = elementsCountThroughWall + 1 + e2 * now + eTC + 2 + tcOffset1
            else:
                bni21 = e2 * now + elementsCountThroughWall * elementsCountAround + eTC + 1 + tcOffset1
                bni22 = e2 * now + elementsCountThroughWall * elementsCountAround + eTC + 2 + tcOffset1
                bni31 = (e2 + 1) * now + eTC + 1 + tcOffset1
                bni32 = (e2 + 1) * now + eTC + 2 + tcOffset1

            if eTC < int(elementsCountAroundTC * 0.5) - 1:
                nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset, bni22 + now + tcOffset,
                                   bni31, bni32, bni31 + now + tcOffset, bni32 + now + tcOffset]
            else:
                nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset,
                                   bni22 + now + tcOffset, bni31, bni31 + now + tcOffset]

            nodeIdentifiers = [c + startNode - 1 for c in nodeIdentifiers]
            if e2 == 0 and nodeIdProximal:
                for m in range(len(nodeIdentifiers)):
                    if nodeIdentifiers[m] in nodeList:
                        idx = nodeList.index(nodeIdentifiers[m])
                        nodeIdentifiers[m] = newNodeList[idx]

            element = \
                mesh.createElement(elementIdentifier,
                                   elementtemplate if eTC < int(elementsCountAroundTC * 0.5) - 1 else elementtemplate1)
            element.setNodesByIdentifier(eft if eTC < int(elementsCountAroundTC * 0.5) - 1 else eft1, nodeIdentifiers)
            if xFlat:
                element.merge(flatElementtemplate1 if eTC < int(elementsCountAroundTC * 0.5) - 1 else
                              flatElementtemplate3)
                element.setNodesByIdentifier(eftFlat3 if eTC < int(elementsCountAroundTC * 0.5) - 1 else
                                             eftFlat5, nodeIdentifiers)
            if xOrgan:
                element.merge(organElementtemplate if eTC < int(elementsCountAroundTC * 0.5) - 1 else
                              organElementtemplate1)
                element.setNodesByIdentifier(eftOrgan if eTC < int(elementsCountAroundTC * 0.5) - 1 else
                                             eftOrgan1, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1
            annotationGroups = annotationGroupsAround[elementsCountAround + eTC] + annotationGroupsAlong[e2] + \
                               ([longitudinalMuscleGroup] if longitudinalMuscle else [])
            if annotationGroups:
                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                for annotationGroup in annotationGroups:
                    meshGroup = annotationGroup.getMeshGroup(mesh)
                    meshGroup.addElement(element)

        for N in range(tcCount - 1):
            for eTC in range(elementsCountAroundTC):
                if closedProximalEnd:
                    bni21 = elementsCountThroughWall + 1 + (
                            e2 - 1) * now + elementsCountThroughWall * elementsCountAround \
                            + eTC + 1 + tcOffset1 + int(elementsCountAroundTC * 0.5) + \
                            (N + 1) * elementsCountAroundHaustrum + N * elementsCountAroundTC + startNode - 1
                    bni22 = elementsCountThroughWall + 1 + (
                            e2 - 1) * now + elementsCountThroughWall * elementsCountAround + \
                            eTC + 2 + tcOffset1 + int(elementsCountAroundTC * 0.5) + \
                            (N + 1) * elementsCountAroundHaustrum + N * elementsCountAroundTC + startNode - 1
                    bni31 = elementsCountThroughWall + 1 + e2 * now + eTC + 1 + tcOffset1 + \
                            int(elementsCountAroundTC * 0.5) - 1 + N * (elementsCountAroundTC - 1) + startNode - 1
                    bni32 = elementsCountThroughWall + 1 + e2 * now + eTC + 2 + tcOffset1 + \
                            int(elementsCountAroundTC * 0.5) - 1 + N * (elementsCountAroundTC - 1) + startNode - 1
                else:
                    bni21 = e2 * now + elementsCountThroughWall * elementsCountAround + eTC + 1 + tcOffset1 + \
                            int(elementsCountAroundTC * 0.5) + (
                                    N + 1) * elementsCountAroundHaustrum + N * elementsCountAroundTC + startNode - 1
                    bni22 = e2 * now + elementsCountThroughWall * elementsCountAround + eTC + 2 + tcOffset1 + \
                            int(elementsCountAroundTC * 0.5) + (N + 1) * elementsCountAroundHaustrum + \
                            N * elementsCountAroundTC + startNode - 1
                    bni31 = (e2 + 1) * now + eTC + 1 + tcOffset1 + int(elementsCountAroundTC * 0.5) - 1 + \
                            N * (elementsCountAroundTC - 1) + startNode - 1
                    bni32 = (e2 + 1) * now + eTC + 2 + tcOffset1 + int(elementsCountAroundTC * 0.5) - 1 + \
                            N * (elementsCountAroundTC - 1) + startNode - 1
                if eTC == 0:
                    nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset,
                                       bni22 + now + tcOffset, bni32, bni32 + now + tcOffset]
                    if e2 == 0 and nodeIdProximal:
                        for m in range(len(nodeIdentifiers)):
                            if nodeIdentifiers[m] in nodeList:
                                idx = nodeList.index(nodeIdentifiers[m])
                                nodeIdentifiers[m] = newNodeList[idx]
                    element = mesh.createElement(elementIdentifier, elementtemplate2)
                    element.setNodesByIdentifier(eft2, nodeIdentifiers)
                    if xFlat:
                        element.merge(flatElementtemplate4)
                        element.setNodesByIdentifier(eftFlat6, nodeIdentifiers)
                    if xOrgan:
                        element.merge(organElementtemplate2)
                        element.setNodesByIdentifier(eftOrgan2, nodeIdentifiers)
                elif 0 < eTC < elementsCountAroundTC - 1:
                    nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset, bni22 + now + tcOffset,
                                       bni31, bni32, bni31 + now + tcOffset, bni32 + now + tcOffset]
                    if e2 == 0 and nodeIdProximal:
                        for m in range(len(nodeIdentifiers)):
                            if nodeIdentifiers[m] in nodeList:
                                idx = nodeList.index(nodeIdentifiers[m])
                                nodeIdentifiers[m] = newNodeList[idx]
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, nodeIdentifiers)
                    if xFlat:
                        element.merge(flatElementtemplate1)
                        element.setNodesByIdentifier(eftFlat3, nodeIdentifiers)
                    if xOrgan:
                        element.merge(organElementtemplate)
                        element.setNodesByIdentifier(eftOrgan, nodeIdentifiers)
                else:
                    nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset,
                                       bni22 + now + tcOffset, bni31, bni31 + now + tcOffset]
                    if e2 == 0 and nodeIdProximal:
                        for m in range(len(nodeIdentifiers)):
                            if nodeIdentifiers[m] in nodeList:
                                idx = nodeList.index(nodeIdentifiers[m])
                                nodeIdentifiers[m] = newNodeList[idx]
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if xFlat:
                        element.merge(flatElementtemplate3)
                        element.setNodesByIdentifier(eftFlat5, nodeIdentifiers)
                    if xOrgan:
                        element.merge(organElementtemplate1)
                        element.setNodesByIdentifier(eftOrgan1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1
                annotationGroups = \
                    annotationGroupsAround[elementsCountAround + int(elementsCountAroundTC * 0.5) +
                                           N * elementsCountAroundTC + eTC] + annotationGroupsAlong[e2] + \
                    ([longitudinalMuscleGroup] if longitudinalMuscle else [])
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

        for eTC in range(int(elementsCountAroundTC * 0.5)):
            if closedProximalEnd:
                bni21 = elementsCountThroughWall + 1 + (e2 - 1) * now + \
                        elementsCountThroughWall * elementsCountAround + eTC + 1 + tcOffset1 + \
                        int(elementsCountAroundTC * 0.5) + tcCount * elementsCountAroundHaustrum + \
                        (tcCount - 1) * elementsCountAroundTC + startNode - 1
                bni22 = + startNode - 1 + elementsCountThroughWall + 1 + (e2 - 1) * now + \
                        elementsCountThroughWall * elementsCountAround + 1 + \
                        tcOffset1 if eTC == int(elementsCountAroundTC * 0.5 - 1) else bni21 + 1
                bni31 = elementsCountThroughWall + 1 + e2 * now + eTC + 1 + tcOffset1 + \
                        int(elementsCountAroundTC * 0.5) - 1 + (tcCount - 1) * (elementsCountAroundTC - 1) + startNode - 1
                bni32 = startNode - 1 + elementsCountThroughWall + 1 + e2 * now + 1 + tcOffset1 if eTC == int(
                    elementsCountAroundTC * 0.5 - 1) \
                    else bni31 + 1
            else:
                bni21 = e2 * now + elementsCountThroughWall * elementsCountAround + eTC + 1 + tcOffset1 + \
                        int(elementsCountAroundTC * 0.5) + tcCount * elementsCountAroundHaustrum + \
                        (tcCount - 1) * elementsCountAroundTC + startNode - 1
                bni22 = e2 * now + elementsCountThroughWall * elementsCountAround + 1 + tcOffset1 + startNode - 1 if eTC == int(
                    elementsCountAroundTC * 0.5 - 1) else bni21 + 1
                bni31 = (e2 + 1) * now + eTC + 1 + tcOffset1 + int(elementsCountAroundTC * 0.5) - 1 + \
                        (tcCount - 1) * (elementsCountAroundTC - 1) + startNode - 1
                bni32 = (e2 + 1) * now + 1 + tcOffset1 + startNode - 1 if eTC == int(elementsCountAroundTC * 0.5 - 1) else bni31 + 1
            if eTC > 0:
                nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset, bni22 + now + tcOffset,
                                   bni31, bni32, bni31 + now + tcOffset, bni32 + now + tcOffset]
            else:
                nodeIdentifiers = [bni21, bni22, bni21 + now + tcOffset,
                                   bni22 + now + tcOffset, bni32, bni32 + now + tcOffset]

            if e2 == 0 and nodeIdProximal:
                for m in range(len(nodeIdentifiers)):
                    if nodeIdentifiers[m] in nodeList:
                        idx = nodeList.index(nodeIdentifiers[m])
                        nodeIdentifiers[m] = newNodeList[idx]

            onOpening = (eTC == int(elementsCountAroundTC * 0.5 - 1))
            element = mesh.createElement(elementIdentifier, elementtemplate if eTC > 0 else elementtemplate2)
            element.setNodesByIdentifier(eft if eTC > 0 else eft2, nodeIdentifiers)
            if xFlat:
                if eTC > 0:
                    element.merge(flatElementtemplate2 if onOpening else flatElementtemplate1)
                    element.setNodesByIdentifier(eftFlat4 if onOpening else eftFlat3, nodeIdentifiers)
                else:
                    element.merge(flatElementtemplate5 if onOpening else flatElementtemplate4)
                    element.setNodesByIdentifier(eftFlat7 if onOpening else eftFlat6, nodeIdentifiers)
            if xOrgan:
                element.merge(organElementtemplate if eTC > 0 else organElementtemplate2)
                element.setNodesByIdentifier(eftOrgan if eTC > 0 else eftOrgan2, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1
            annotationGroups = \
                annotationGroupsAround[elementsCountAround + int(elementsCountAroundTC * (tcCount - 0.5)) + eTC] + \
                annotationGroupsAlong[e2] + ([longitudinalMuscleGroup] if longitudinalMuscle else [])
            if annotationGroups:
                allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                for annotationGroup in annotationGroups:
                    meshGroup = annotationGroup.getMeshGroup(mesh)
                    meshGroup.addElement(element)

    fm.endChange()

    return nodeIdentifier, elementIdentifier, allAnnotationGroups, nodesDistal
