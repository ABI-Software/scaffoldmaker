"""
Generates a single or double/common ostium, where one or more vessels enters a chamber.
"""

from __future__ import division

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import createCirclePoints, getCircleProjectionAxes
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes


class MeshType_3d_ostium1(Scaffold_base):
    """
    Generates a 3-D single or double/common ostium inlet or outlet.
    """
    @staticmethod
    def getName():
        return '3D Ostium 1'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        return {
            'Number of vessels': 2,
            'Number of elements across common': 2,
            'Number of elements around ostium': 10,
            'Number of elements along': 1,
            'Number of elements through wall': 1,
            'Unit scale': 1.0,
            'Outlet': False,
            'Ostium diameter': 1.0,
            'Ostium length': 0.4,
            'Ostium wall thickness': 0.08,
            'Ostium wall relative thicknesses': [1.0],
            'Ostium inter-vessel distance': 0.8,
            'Ostium inter-vessel height': 0.0,
            'Use linear through ostium wall': False,
            'Vessel end length factor': 1.0,
            'Vessel inner diameter': 0.6,
            'Vessel wall thickness': 0.04,
            'Vessel wall relative thicknesses': [1.0],
            'Vessel angle 1 degrees': 0.0,
            'Vessel angle 1 spread degrees': 0.0,
            'Vessel angle 2 degrees': 0.0,
            'Use linear through vessel wall': True,
            #  'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements along': 4,
            'Refine number of elements through wall': 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of vessels',
            'Number of elements across common',
            'Number of elements around ostium',
            'Number of elements along',
            'Number of elements through wall',
            'Unit scale',
            'Outlet',
            'Ostium diameter',
            'Ostium length',
            'Ostium wall thickness',
            'Ostium wall relative thicknesses',
            'Ostium inter-vessel distance',
            'Ostium inter-vessel height',
            'Use linear through ostium wall',
            'Vessel end length factor',
            'Vessel inner diameter',
            'Vessel wall thickness',
            'Vessel wall relative thicknesses',
            'Vessel angle 1 degrees',
            'Vessel angle 1 spread degrees',
            'Vessel angle 2 degrees',
            'Use linear through vessel wall',
            #  'Use cross derivatives',  # not implemented
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        vesselsCount = options['Number of vessels']
        if vesselsCount < 1:
            vesselsCount = options['Number of vessels'] = 1
        elif vesselsCount > 3:
            vesselsCount = options['Number of vessels'] = 3
        for key in [
            'Number of elements across common',
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around ostium'] < 2 * vesselsCount:
            options['Number of elements around ostium'] = 2 * vesselsCount
            dependentChanges = True  # because can happen by changing number of vessels
        # currently must have even number around ostium if multiple vessels:
        if (vesselsCount > 1) and (options['Number of elements around ostium'] % 2):
            options['Number of elements around ostium'] += 1
            dependentChanges = True  # because can happen by changing number of vessels
        for key in [
            'Unit scale',
            'Ostium length',
            'Ostium wall thickness',
            'Ostium inter-vessel distance',
            'Vessel inner diameter',
            'Vessel wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Ostium diameter'] <= 0.0:
            options['Ostium diameter'] = 0.000001  # avoid division by zero
        elementsThroughWall = options['Number of elements through wall']
        ostiumThicknessProportionsCountKey = 'Ostium wall relative thicknesses'
        vesselThicknessProportionsCountKey = 'Vessel wall relative thicknesses'
        ostiumWallCount = len(options[ostiumThicknessProportionsCountKey])
        vesselWallCount = len(options[vesselThicknessProportionsCountKey])
        if elementsThroughWall == 1:
            options[ostiumThicknessProportionsCountKey] = options[vesselThicknessProportionsCountKey] = [1.0]
        if ostiumWallCount < elementsThroughWall:
            options[ostiumThicknessProportionsCountKey] += \
                [options[ostiumThicknessProportionsCountKey][-1] for i in range(elementsThroughWall - ostiumWallCount)]
            dependentChanges = True
        elif ostiumWallCount > elementsThroughWall:
            options[ostiumThicknessProportionsCountKey] = \
                options[ostiumThicknessProportionsCountKey][:elementsThroughWall]
            dependentChanges = True
        if vesselWallCount < elementsThroughWall:
            options[vesselThicknessProportionsCountKey] += \
                [options[vesselThicknessProportionsCountKey][-1] for i in range(elementsThroughWall - vesselWallCount)]
            dependentChanges = True
        elif vesselWallCount > elementsThroughWall:
            options[vesselThicknessProportionsCountKey] = \
                options[vesselThicknessProportionsCountKey][:elementsThroughWall]
            dependentChanges = True
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic/bicubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        unitScale = options['Unit scale']
        ostiumRadius = 0.5 * unitScale * options['Ostium diameter']
        interVesselDistance = unitScale * options['Ostium inter-vessel distance']

        scale = 1.1 * (ostiumRadius * 2.0 + interVesselDistance)
        nx = [[-scale, -scale, 0.0], [scale, -scale, 0.0], [-scale, scale, 0.0], [scale, scale, 0.0]]
        nd1 = [[2.0 * scale, 0.0, 0.0]] * 4
        nd2 = [[0.0, 2.0 * scale, 0.0]] * 4
        # #  curve track surface:
        # zf = 2.0
        # nd1 = [ [ 2.0 * scale, 0.0, zf * scale ], [ 2.0 * scale, 0.0, -zf * scale ],
        #         [ 2.0 * scale, 0.0, zf * scale ], [ 2.0 * scale, 0.0, -zf * scale ] ]
        # nd2 = [ [ 0.0, 2.0 * scale, zf * scale ], [ 0.0, 2.0 * scale, zf * scale ],
        #         [ 0.0, 2.0 * scale, -zf * scale ], [ 0.0, 2.0 * scale, -zf * scale ] ]
        trackSurface = TrackSurface(1, 1, nx, nd1, nd2)
        centrePosition = TrackSurfacePosition(0, 0, 0.5, 0.5)
        axis1 = [1.0, 0.0, 0.0]
        generateOstiumMesh(region, options, trackSurface, centrePosition, axis1)
        return []  # no annotation groups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)


def getOstiumElementsCountsAroundVessels(elementsCountAroundOstium, elementsCountAcross, vesselsCount):
    """
    Determine numbers of elements around each vessel to fit numbers of elements around ostium.
    :param elementsCountAroundOstium: Number of elements around outside of ostium.
    Minimum of 2 + 2 * vesselsCount. Must be even if more than 1 vessels.
    :param elementsCountAcross: Number of elements across ostium between multiple vessels.
    Unused if vesselsCount is 1.
    :param vesselsCount: Number of vessels from 1 to 3.
    :return: List of numbers of elements around each vessel to fit numbers around ostium,
    number of elements mid side (non-zero only for 3+ vessels).
    """
    assert 1 <= vesselsCount <= 3
    assert elementsCountAroundOstium >= 2 * vesselsCount
    assert (1 == vesselsCount) or (elementsCountAroundOstium % 2 == 0)
    assert elementsCountAcross > 0
    if vesselsCount == 1:
        return [elementsCountAroundOstium], 0
    if vesselsCount == 2:
        count = elementsCountAroundOstium // 2 + elementsCountAcross
        return [count, count], 0
    # vesselsCount == 3:
    # Around oinc   Total:Vessels 1  Total:Vessels 2  
    # 6      1       8: 3-4-3        10: 4-6-4        
    # 8      1      10: 4-4-4        12: 5-6-5        
    # 10     1      12: 5-4-5        14: 6-6-6        
    # 12     2      14: 5-6-5        16: 6-8-6        
    # 14     2      16: 6-6-6        18: 7-8-7        
    # 16     2      18: 7-6-7        20: 8-8-8        
    elementsCountAroundMid = ((elementsCountAroundOstium + 1) // 6)
    countInner = 2 * (elementsCountAroundMid + elementsCountAcross)
    countOuter = (elementsCountAroundOstium - 2 * elementsCountAroundMid) // 2 + elementsCountAcross
    return [countOuter, countInner, countOuter], elementsCountAroundMid


def generateOstiumMesh(region, options, trackSurface, centrePosition, axis1, startNodeIdentifier=1,
                       startElementIdentifier=1, vesselMeshGroups=None, ostiumMeshGroups=None,
                       wallAnnotationGroups=None, coordinates=None):
    """
    :param vesselMeshGroups: List (over number of vessels) of list of mesh groups to add vessel elements to.
    :param ostiumMeshGroups: List of mesh groups to add only row of elements at ostium end to.
    :param wallAnnotationGroups: list of annotation groups to add to wall elements.
    :return: nextNodeIdentifier, nextElementIdentifier, Ostium points tuple
    (ox[n3][n1][c], od1[n3][n1][c], od2[n3][n1][c], od3[n3][n1][c], oNodeId[n3][n1], oPositions).
    """
    vesselsCount = options['Number of vessels']
    elementsCountAroundOstium = options['Number of elements around ostium']
    elementsCountAcross = options['Number of elements across common']
    elementsCountsAroundVessels, elementsCountAroundMid = \
        getOstiumElementsCountsAroundVessels(elementsCountAroundOstium, elementsCountAcross, vesselsCount)
    elementsCountAroundEnd = (elementsCountAroundOstium - 2 * elementsCountAroundMid) // 2
    # print('\nvesselsCount', vesselsCount, 'elementsCountsAroundOstium', elementsCountAroundOstium,
    #       'elementsCountAcross', elementsCountAcross)
    # print('--> elementsCountsAroundVessels', elementsCountsAroundVessels,
    #       'elementsCountAroundMid', elementsCountAroundMid)
    elementsCountAlong = options['Number of elements along']
    elementsCountThroughWall = options['Number of elements through wall']
    unitScale = options['Unit scale']

    isOutlet = options['Outlet']
    ostiumRadius = 0.5 * unitScale * options['Ostium diameter']
    ostiumLength = unitScale * options['Ostium length']
    ostiumWallThickness = unitScale * options['Ostium wall thickness']
    ostiumWallThicknessProportionsUI = copy.deepcopy(options['Ostium wall relative thicknesses'])
    interVesselHeight = unitScale * options['Ostium inter-vessel height']
    interVesselDistance = unitScale * options['Ostium inter-vessel distance'] if (vesselsCount > 1) else 0.0
    halfInterVesselDistance = 0.5 * interVesselDistance
    useCubicHermiteThroughOstiumWall = not(options['Use linear through ostium wall'])
    vesselEndDerivative = ostiumLength * options['Vessel end length factor'] / elementsCountAlong
    vesselInnerRadius = 0.5 * unitScale * options['Vessel inner diameter']
    vesselWallThickness = unitScale * options['Vessel wall thickness']
    vesselWallThicknessProportionsUI = copy.deepcopy(options['Vessel wall relative thicknesses'])
    vesselOuterRadius = vesselInnerRadius + vesselWallThickness
    vesselAngle1Radians = math.radians(options['Vessel angle 1 degrees'])
    vesselAngle1SpreadRadians = math.radians(options['Vessel angle 1 spread degrees'])
    vesselAngle2Radians = math.radians(options['Vessel angle 2 degrees'])
    useCubicHermiteThroughVesselWall = not(options['Use linear through vessel wall'])
    # useCrossDerivatives = False  # options['Use cross derivatives']  # not implemented

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    if not coordinates:
        coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodeIdentifier = startNodeIdentifier

    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    nodetemplateLinearS3 = nodes.createNodetemplate()
    nodetemplateLinearS3.defineField(coordinates)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

    mesh = fm.findMeshByDimension(3)
    elementIdentifier = startElementIdentifier

    # track points in shape of ostium

    # get directions in plane of surface at centre:
    cx, cd1, cd2 = trackSurface.evaluateCoordinates(centrePosition, True)
    trackDirection1, trackDirection2, centreNormal = calculate_surface_axes(cd1, cd2, axis1)
    trackDirection2reverse = [-d for d in trackDirection2]

    halfCircumference = math.pi * ostiumRadius
    circumference = 2.0 * halfCircumference
    distance = 0.0
    elementLengthAroundOstiumMid = 0.0
    vesselsSpanAll = interVesselDistance * (vesselsCount - 1)
    vesselsSpanMid = interVesselDistance * (vesselsCount - 2)
    if vesselsCount == 1:
        elementLengthAroundOstiumEnd = circumference / elementsCountAroundOstium
        vesselOstiumPositions = [centrePosition]
        ocx = [cx]
        ocd1 = [trackDirection1]
        ocd2 = [trackDirection2]
        ocd3 = [centreNormal]
    else:
        elementLengthAroundOstiumEnd = (circumference + 2.0 * interVesselDistance) / \
                                       (elementsCountAroundOstium - 2 * elementsCountAroundMid)
        if elementsCountAroundMid > 0:
            elementLengthAroundOstiumMid = interVesselDistance * (vesselsCount - 2) / elementsCountAroundMid
        vesselOstiumPositions = []
        ocx = []
        ocd1 = []
        ocd2 = []
        ocd3 = []
        for v in range(vesselsCount):
            vesselOstiumPositions.append(trackSurface.trackVector(centrePosition, trackDirection1,
                                                                  (v / (vesselsCount - 1) - 0.5) * vesselsSpanAll))
            x, d1, d2 = trackSurface.evaluateCoordinates(vesselOstiumPositions[-1], -1)
            d1, d2, d3 = calculate_surface_axes(d1, d2, trackDirection1)
            ocx .append(x)
            ocd1.append(d1)
            ocd2.append(d2)
            ocd3.append(d3)

    # coordinates around ostium
    ox = [[] for n3 in range(elementsCountThroughWall + 1)]
    od1 = [[] for n3 in range(elementsCountThroughWall + 1)]
    od2 = [[] for n3 in range(elementsCountThroughWall + 1)]
    od3 = [[] for n3 in range(elementsCountThroughWall + 1)]
    oPositions = []
    ostiumWallThicknessProportions = [ostiumWallThicknessProportion / sum(ostiumWallThicknessProportionsUI)
                                      for ostiumWallThicknessProportion in ostiumWallThicknessProportionsUI]
    vesselWallThicknessProportions = [vesselWallThicknessProportion / sum(vesselWallThicknessProportionsUI)
                                      for vesselWallThicknessProportion in vesselWallThicknessProportionsUI]

    ostiumWallThicknessProportions.append(ostiumWallThicknessProportions[-1])
    vesselWallThicknessProportions.append(vesselWallThicknessProportions[-1])

    for n1 in range(elementsCountAroundOstium):
        elementLength = elementLengthAroundOstiumEnd
        if distance <= (vesselsSpanMid + halfInterVesselDistance):
            position = trackSurface.trackVector(centrePosition, trackDirection1, 0.5 * vesselsSpanMid - distance)
            sideDirection = trackDirection2reverse
            if n1 < elementsCountAroundMid:
                elementLength = elementLengthAroundOstiumMid
        elif distance < (vesselsSpanMid + halfInterVesselDistance + halfCircumference):
            position = vesselOstiumPositions[0]
            angleRadians = (distance - (vesselsSpanMid + halfInterVesselDistance)) / ostiumRadius
            w1 = -math.sin(angleRadians)
            w2 = -math.cos(angleRadians)
            sideDirection = [(w1 * trackDirection1[c] + w2 * trackDirection2[c]) for c in range(3)]
        elif distance < (2.0 * vesselsSpanMid + halfInterVesselDistance + halfCircumference + interVesselDistance):
            position = trackSurface.trackVector(centrePosition, trackDirection1, distance -
                                                (1.5 * vesselsSpanMid + interVesselDistance + halfCircumference))
            sideDirection = trackDirection2
            if 0 <= (n1 - elementsCountAroundEnd - elementsCountAroundMid) < elementsCountAroundMid:
                elementLength = elementLengthAroundOstiumMid
        elif distance < (2.0 * vesselsSpanMid + halfInterVesselDistance + circumference + interVesselDistance):
            position = vesselOstiumPositions[-1]
            angleRadians = (distance - (2.0 * vesselsSpanMid + halfInterVesselDistance + halfCircumference +
                                        interVesselDistance)) / ostiumRadius
            w1 = math.sin(angleRadians)
            w2 = math.cos(angleRadians)
            sideDirection = [(w1 * trackDirection1[c] + w2 * trackDirection2[c]) for c in range(3)]
        else:
            position = \
                trackSurface.trackVector(centrePosition, trackDirection1,
                                         0.5 * vesselsSpanMid + (circumference + 2.0 *
                                                                 (vesselsSpanMid + interVesselDistance)) - distance)
            sideDirection = trackDirection2reverse
        position = trackSurface.trackVector(position, sideDirection, ostiumRadius)
        oPositions.append(position)
        px, d1, d2 = trackSurface.evaluateCoordinates(position, True)
        pd2, pd1, pd3 = calculate_surface_axes(d1, d2, sideDirection)
        # get outer coordinates
        opx = px
        opd1 = vector.setMagnitude([-d for d in pd1], elementLengthAroundOstiumEnd)
        opd2 = vector.setMagnitude(pd2, elementLengthAroundOstiumEnd)  # smoothed later
        opd3 = vector.setMagnitude(pd3, ostiumWallThickness)

        ostiumTotalXi3 = 0.0
        vesselTotalXi3 = 0.0
        ostiumWallThicknessXi3List = [0.0]
        vesselWallThicknessXi3List = [0.0]

        for n3 in range(elementsCountThroughWall):
            ostiumTotalXi3 += ostiumWallThicknessProportions[n3]
            vesselTotalXi3 += vesselWallThicknessProportions[n3]
            ostiumWallThicknessXi3List.append(ostiumTotalXi3)
            vesselWallThicknessXi3List.append(vesselTotalXi3)

        # set coordinates through wall (use copy to avoid references to same list later)
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1 - ostiumWallThicknessXi3List[n3]
            ox[n3].append([(opx[c] - opd3[c] * xi3) for c in range(3)])
            od1[n3].append(copy.copy(opd1))
            od2[n3].append(copy.copy(opd2))
            if useCubicHermiteThroughOstiumWall:
                od3[n3].append([opd3[c] * ostiumWallThicknessProportions[n3] for c in range(3)])
        distance += elementLength
    for n3 in range(elementsCountThroughWall + 1):
        od1[n3] = interp.smoothCubicHermiteDerivativesLoop(ox[n3], od1[n3], fixAllDirections=True)

    xx = []
    xd1 = []
    xd2 = []
    xd3 = []
    if (vesselWallThickness > 0.0) and (ostiumWallThickness > 0.0):
        commonOstiumWallThickness = 2.0 / (1.0 / vesselWallThickness + 1.0 / ostiumWallThickness)
        commonOstiumWallThicknessProportions = []
        commonOstiumWallTotalXi3 = 0.0
        commonOstiumWallThicknessXi3List = [0.0]
        for c in range(elementsCountThroughWall + 1):
            commonOstiumWallThicknessProportions.append(
                (vesselWallThicknessProportions[c] + ostiumWallThicknessProportions[c]) * 0.5)
        for n3 in range(elementsCountThroughWall):
            commonOstiumWallTotalXi3 += commonOstiumWallThicknessProportions[n3]
            commonOstiumWallThicknessXi3List.append(commonOstiumWallTotalXi3)
    else:
        commonOstiumWallThickness = vesselWallThickness
        commonOstiumWallThicknessProportions = vesselWallThicknessProportions
        commonOstiumWallThicknessXi3List = vesselWallThicknessXi3List

    # coordinates across common ostium, between vessels
    nodesCountFreeEnd = elementsCountsAroundVessels[0] + 1 - elementsCountAcross
    oinc = 0 if (vesselsCount <= 2) else elementsCountAroundMid // (vesselsCount - 2)
    for iv in range(vesselsCount - 1):
        xx .append([None for n3 in range(elementsCountThroughWall + 1)])
        xd1.append([None for n3 in range(elementsCountThroughWall + 1)])
        xd2.append([None for n3 in range(elementsCountThroughWall + 1)])
        xd3.append([None for n3 in range(elementsCountThroughWall + 1)])
        oa = elementsCountAroundMid - iv * oinc
        ob = elementsCountAroundMid + nodesCountFreeEnd - 1 + iv * oinc
        nx = [ox[elementsCountThroughWall][oa], ox[elementsCountThroughWall][ob]]
        nd1 = [[-d for d in od1[elementsCountThroughWall][oa]], od1[elementsCountThroughWall][ob]]
        nd2 = [[-d for d in od2[elementsCountThroughWall][oa]], od2[elementsCountThroughWall][ob]]

        if elementsCountAcross > 1:
            # add centre point, displaced by interVesselHeight
            if vesselsCount == 2:
                position = centrePosition
            else:
                position = trackSurface.trackVector(centrePosition, trackDirection1,
                                                    (iv / (vesselsCount - 2) - 0.5) * vesselsSpanMid)
            mx, d1, d2 = trackSurface.evaluateCoordinates(position, derivatives=True)
            md1, md2, md3 = calculate_surface_axes(d1, d2, trackDirection1)
            nx .insert(1, [(mx[c] + interVesselHeight * md3[c]) for c in range(3)])
            nd1.insert(1, vector.setMagnitude(md1, elementLengthAroundOstiumMid if (0 < iv < (vesselsCount - 2)) else
                                                   elementLengthAroundOstiumEnd))
            nd2.insert(1, vector.setMagnitude(md2, ostiumRadius))
        nd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True)
        px, pd2, pe, pxi = interp.sampleCubicHermiteCurves(nx, nd2, elementsCountAcross)[0:4]
        pd1 = interp.interpolateSampleLinear(nd1, pe, pxi)
        pd3 = [vector.setMagnitude(vector.crossproduct3(pd1[n2], pd2[n2]), commonOstiumWallThickness)
               for n2 in range(elementsCountAcross + 1)]
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1 - commonOstiumWallThicknessXi3List[n3]
            lx = [([(px[n2][c] - xi3 * pd3[n2][c]) for c in range(3)]) for n2 in range(elementsCountAcross + 1)]
            ld2 = interp.smoothCubicHermiteDerivativesLine(lx, pd2, fixAllDirections=True)
            xx[iv][n3] = lx[1:elementsCountAcross]
            xd1[iv][n3] = copy.deepcopy(pd1[1:elementsCountAcross])  # to be smoothed later
            xd2[iv][n3] = ld2[1:elementsCountAcross]
            # set smoothed d2 on ostium circumference
            od2[n3][oa] = [-d for d in ld2[0]]
            od2[n3][ob] = ld2[-1]
            if useCubicHermiteThroughOstiumWall:
                pd3Element = [vector.setMagnitude(vector.crossproduct3(pd1[n2], pd2[n2]),
                                                  commonOstiumWallThickness * commonOstiumWallThicknessProportions[n3])
                              for n2 in range(elementsCountAcross + 1)]
                xd3[iv][n3] = copy.deepcopy(pd3Element[1:elementsCountAcross])

    # get positions of vessel end centres and rings
    vcx = []
    vcd1 = []
    vcd2 = []
    vcd3 = []
    vox = []
    vod1 = []
    vod2 = []
    vod3 = []
    for v in range(vesselsCount):
        elementsCountAroundVessel = elementsCountsAroundVessels[v]
        radiansPerElementVessel = 2.0 * math.pi / elementsCountAroundVessel
        useVesselAngleRadians = vesselAngle1Radians
        if vesselsCount > 1:
            useVesselAngleRadians += (v / (vesselsCount - 1) - 0.5) * vesselAngle1SpreadRadians
        vx, vd1, vd2, vd3 = getCircleProjectionAxes(ocx[v], ocd1[v], ocd2[v], ocd3[v], ostiumLength,
                                                    useVesselAngleRadians, vesselAngle2Radians)
        vd1 = [vesselOuterRadius * d for d in vd1]
        vd2 = [-vesselOuterRadius * d for d in vd2]
        vd3 = [-vesselEndDerivative * d for d in vd3]
        vcx.append(vx)
        vcd1.append(vd1)
        vcd2.append(vd2)
        vcd3.append(vd3)
        vox.append([])
        vod1.append([])
        vod2.append([])
        vod3.append([])
        for n3 in range(elementsCountThroughWall + 1):
            radius = vesselInnerRadius + vesselWallThicknessXi3List[n3] * vesselWallThickness
            vAxis1 = vector.setMagnitude(vd1, radius)
            vAxis2 = vector.setMagnitude(vd2, radius)
            if vesselsCount == 1:
                startRadians = 0.5 * math.pi
            else:
                startRadians = 0.5 * radiansPerElementVessel * elementsCountAcross
                if v == (vesselsCount - 1):
                    startRadians -= math.pi
            px, pd1 = createCirclePoints(vx, vAxis1, vAxis2, elementsCountAroundVessel, startRadians)
            vox[-1].append(px)
            vod1[-1].append(pd1)
            vod2[-1].append([vd3] * elementsCountAroundVessel)
            if useCubicHermiteThroughVesselWall:
                vod3[-1].append([vector.setMagnitude(vector.crossproduct3(d1, vd3),
                                                     vesselWallThickness * vesselWallThicknessProportions[n3])
                                 for d1 in pd1])

    # calculate common ostium vessel node derivatives map
    mvPointsx = [None] * vesselsCount
    mvPointsd1 = [None] * vesselsCount
    mvPointsd2 = [None] * vesselsCount
    mvPointsd3 = [None] * vesselsCount
    mvDerivativesMap = [None] * vesselsCount
    mvMeanCount = [None] * vesselsCount  # stores 1 if first reference to common point between vessels, 2 if second. Otherwise 0.
    for v in range(vesselsCount):
        if vesselsCount == 1:
            mvPointsx[v], mvPointsd1[v], mvPointsd2[v], mvPointsd3[v], mvDerivativesMap[v] = \
                ox, od1, od2, od3 if useCubicHermiteThroughOstiumWall else None, None
            mvMeanCount[v] = [0] * elementsCountsAroundVessels[v]
        else:
            iv = max(0, v - 1)
            oa = elementsCountAroundMid - iv * oinc
            ob = elementsCountAroundMid + nodesCountFreeEnd - 1 + iv * oinc
            mvPointsx[v] = []
            mvPointsd1[v] = []
            mvPointsd2[v] = []
            mvPointsd3[v] = [] if useCubicHermiteThroughOstiumWall else None
            mvDerivativesMap[v] = []
            for n3 in range(elementsCountThroughWall + 1):
                mvPointsd1[v].append([])
                mvPointsd2[v].append([])
                mvPointsx[v].append([])
                if useCubicHermiteThroughOstiumWall:
                    mvPointsd3[v].append([])
                mvDerivativesMap[v].append([])
                if v == 0:  # first end vessel
                    mvPointsd1[v][n3] += od1[n3][oa:ob + 1]
                    mvPointsd2[v][n3] += od2[n3][oa:ob + 1]
                    mvPointsx[v][n3] += ox[n3][oa:ob + 1]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += od3[n3][oa:ob + 1]
                    mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 1, 0), None, (1, 0, 0)))
                    for i in range(nodesCountFreeEnd - 2):
                        mvDerivativesMap[v][n3].append((None, None, None))
                    mvDerivativesMap[v][n3].append(((1, 0, 0), (1, 1, 0), None, (0, -1, 0)))
                    mvPointsx[v][n3] += reversed(xx[iv][n3])
                    mvPointsd1[v][n3] += reversed(xd1[iv][n3])
                    mvPointsd2[v][n3] += reversed(xd2[iv][n3])
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += reversed(xd3[iv][n3])
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append(((0, -1, 0), (1, 0, 0), None))
                    if n3 == 0:
                        mvMeanCount[v] = [1] + [0] * (nodesCountFreeEnd - 2) + [1] * elementsCountAcross
                elif v < (vesselsCount - 1):  # middle vessels
                    # left:
                    mvPointsx[v][n3] += ox[n3][oa - oinc:oa + 1]
                    mvPointsd1[v][n3] += od1[n3][oa - oinc:oa + 1]
                    mvPointsd2[v][n3] += od2[n3][oa - oinc:oa + 1]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += od3[n3][oa - oinc:oa + 1]
                    mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 1, 0), None, (1, 0, 0)))
                    for i in range(oinc - 1):
                        mvDerivativesMap[v][n3].append((None, None, None))
                    mvDerivativesMap[v][n3].append(((1, 0, 0), (1, 1, 0), None, (0, -1, 0)))
                    # across
                    mvPointsx[v][n3] += xx[iv][n3]
                    mvPointsd1[v][n3] += xd1[iv][n3]
                    mvPointsd2[v][n3] += xd2[iv][n3]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += xd3[iv][n3]
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 0, 0), None))
                    # right
                    mvPointsx[v][n3] += ox[n3][ob:ob + oinc + 1]
                    mvPointsd1[v][n3] += od1[n3][ob:ob + oinc + 1]
                    mvPointsd2[v][n3] += od2[n3][ob:ob + oinc + 1]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += od3[n3][ob:ob + oinc + 1]
                    mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 1, 0), None, (1, 0, 0)))
                    for i in range(oinc - 1):
                        mvDerivativesMap[v][n3].append((None, None, None))
                    mvDerivativesMap[v][n3].append(((1, 0, 0), (1, 1, 0), None, (0, -1, 0)))
                    # across reverse
                    mvPointsx[v][n3] += reversed(xx[iv + 1][n3])
                    mvPointsd1[v][n3] += reversed(xd1[iv + 1][n3])
                    mvPointsd2[v][n3] += reversed(xd2[iv + 1][n3])
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += reversed(xd3[iv + 1][n3])
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append(((0, -1, 0), (1, 0, 0), None))
                    if n3 == 0:
                        mvMeanCount[v] = [1] + [0] * (oinc - 1) + [2] * (elementsCountAcross + 1) + [0] * (oinc - 1) + \
                                         [1] * elementsCountAcross
                else:  # last end vessel
                    mvPointsx[v][n3] += ox[n3][ob:] + [ox[n3][0]]
                    mvPointsd1[v][n3] += od1[n3][ob:] + [od1[n3][0]]
                    mvPointsd2[v][n3] += od2[n3][ob:] + [od2[n3][0]]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += od3[n3][ob:] + [od3[n3][0]]
                    mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 1, 0), None, (1, 0, 0)))
                    for i in range(nodesCountFreeEnd - 2):
                        mvDerivativesMap[v][n3].append((None, None, None))
                    mvDerivativesMap[v][n3].append(((1, 0, 0), (1, 1, 0), None, (0, -1, 0)))
                    mvPointsx[v][n3] += xx[iv][n3]
                    mvPointsd1[v][n3] += xd1[iv][n3]
                    mvPointsd2[v][n3] += xd2[iv][n3]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += xd3[iv][n3]
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append(((0, 1, 0), (-1, 0, 0), None))
                    if n3 == 0:
                        mvMeanCount[v] = [2] + [0] * (nodesCountFreeEnd - 2) + [2] * elementsCountAcross

    # calculate derivative 2 around free sides of inlets to fit vessel derivatives
    for v in range(vesselsCount):
        for n3 in [elementsCountThroughWall]:  # was range(2), now using curvature for inside:
            # print('v',v,'n3',n3,'elementsAround',elementsCountsAroundVessels[v])
            # print('mvPointsx [v][n3]', mvPointsx [v][n3])
            # print('mvPointsd1[v][n3]', mvPointsd1[v][n3])
            # print('mvPointsd2[v][n3]', mvPointsd2[v][n3])
            # print('mvDerivativesMap[v][n3]', mvDerivativesMap[v][n3])
            for n1 in range(elementsCountsAroundVessels[v]):
                d2Map = mvDerivativesMap[v][n3][n1][1] if (mvDerivativesMap[v] and mvDerivativesMap[v][n3][n1]) \
                    else None
                sf1 = d2Map[0] if d2Map else 0.0
                sf2 = d2Map[1] if d2Map else 1.0
                nx = [vox[v][n3][n1], mvPointsx[v][n3][n1]]
                nd2 = [[d * elementsCountAlong for d in vod2[v][n3][n1]],
                       [(sf1 * mvPointsd1[v][n3][n1][c] + sf2 * mvPointsd2[v][n3][n1][c]) for c in range(3)]]
                nd2f = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative=True, fixEndDirection=True)
                ndf = [d / elementsCountAlong for d in nd2f[1]]
                # assign components to set original values:
                if sf1 == 0:
                    for c in range(3):
                        mvPointsd2[v][n3][n1][c] = sf2 * ndf[c]
                elif sf2 == 0:
                    if mvMeanCount[v][n1] < 2:
                        for c in range(3):
                            mvPointsd1[v][n3][n1][c] = sf1 * ndf[c]
                    else:
                        # take mean of values from this and last vessel
                        for c in range(3):
                            mvPointsd1[v][n3][n1][c] = 0.5 * (mvPointsd1[v][n3][n1][c] + sf1 * ndf[c])
                else:
                    # print('v', v, 'n3', n3, 'n1', n1, ':', vector.magnitude(ndf), 'vs.', vector.magnitude(nd2[1]),
                    #       'd2Map', d2Map)
                    pass

    # calculate inner d2 derivatives around ostium from outer using track surface curvature
    factor = 1.0
    for n1 in range(elementsCountAroundOstium):
        trackDirection = vector.normalise(od2[elementsCountThroughWall][n1])
        trackDistance = factor * vector.magnitude(od2[elementsCountThroughWall][n1])
        tx = [None, ox[elementsCountThroughWall][n1], None]
        td1 = [None, vector.setMagnitude(od2[elementsCountThroughWall][n1], trackDistance), None]
        td2 = [None, vector.setMagnitude(od1[elementsCountThroughWall][n1], -trackDistance), None]
        positionBackward = trackSurface.trackVector(oPositions[n1], trackDirection, -trackDistance)
        tx[0], d1, d2 = trackSurface.evaluateCoordinates(positionBackward, derivatives=True)
        sd1, sd2, sd3 = calculate_surface_axes(d1, d2, trackDirection)
        td1[0] = vector.setMagnitude(sd1, trackDistance)
        td2[0] = vector.setMagnitude(sd2, trackDistance)
        positionForward = trackSurface.trackVector(oPositions[n1], trackDirection, trackDistance)
        tx[2], d1, d2 = trackSurface.evaluateCoordinates(positionForward, derivatives=True)
        sd1, sd2, sd3 = calculate_surface_axes(d1, d2, trackDirection)
        td1[2] = vector.setMagnitude(sd1, trackDistance)
        td2[2] = vector.setMagnitude(sd2, trackDistance)
        for n3 in range(elementsCountThroughWall):
            xi3 = 1 - ostiumWallThicknessXi3List[n3]
            newd2 = interp.projectHermiteCurvesThroughWall(tx, td1, td2, 1, -ostiumWallThickness * xi3)[1]
            # assign components to set in all lists:
            for c in range(3):
                od2[n3][n1][c] = newd2[c] / factor
        # for n in [ 0, 1, 2 ]:
        #    node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
        #    cache.setNode(node)
        #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, tx [n])
        #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, td1[n])
        #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, td2[n])
        #    nodeIdentifier += 1

    if isOutlet:
        # reverse directions of d1 and d2 on vessels and ostium base
        for c in range(3):
            for n3 in range(elementsCountThroughWall + 1):
                for n1 in range(elementsCountAroundOstium):
                    od1[n3][n1][c] = -od1[n3][n1][c]
                    od2[n3][n1][c] = -od2[n3][n1][c]
                for iv in range(vesselsCount - 1):
                    for n1 in range(elementsCountAcross - 1):
                        xd1[iv][n3][n1][c] = -xd1[iv][n3][n1][c]
                        xd2[iv][n3][n1][c] = -xd2[iv][n3][n1][c]
                for v in range(vesselsCount):
                    for n1 in range(elementsCountsAroundVessels[v]):
                        vod1[v][n3][n1][c] = -vod1[v][n3][n1][c]
            # d2 is referenced all around, so only change once per vessel
            for v in range(vesselsCount):
                vod2[v][0][0][c] = -vod2[v][0][0][c]

    ##############
    # Create nodes
    ##############

    oNodeId = []
    for n3 in range(elementsCountThroughWall + 1):
        oNodeId.append([])
        for n1 in range(elementsCountAroundOstium):
            node = nodes.createNode(nodeIdentifier,
                                    nodetemplate if useCubicHermiteThroughOstiumWall else nodetemplateLinearS3)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ox[n3][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, od1[n3][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, od2[n3][n1])
            if useCubicHermiteThroughOstiumWall:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, od3[n3][n1])
            oNodeId[n3].append(nodeIdentifier)
            nodeIdentifier += 1

    xNodeId = []
    for iv in range(vesselsCount - 1):
        xNodeId.append([])
        for n3 in range(elementsCountThroughWall + 1):
            xNodeId[iv].append([])
            for n2 in range(elementsCountAcross - 1):
                node = nodes.createNode(nodeIdentifier,
                                        nodetemplate if useCubicHermiteThroughOstiumWall else nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xx[iv][n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xd1[iv][n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xd2[iv][n3][n2])
                if useCubicHermiteThroughOstiumWall:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xd3[iv][n3][n2])
                xNodeId[iv][n3].append(nodeIdentifier)
                nodeIdentifier += 1

    # for v in range(vesselsCount):
    #    node = nodes.createNode(nodeIdentifier, nodetemplate)
    #    cache.setNode(node)
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vcx [v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vcd1[v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vcd2[v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vcd3[v])
    #    nodeIdentifier += 1
    #    for n3 in range(elementsCountThroughWall + 1):
    #        for n1 in range(elementsCountsAroundVessels[v]):
    #            node = nodes.createNode(nodeIdentifier,
    #                                    nodetemplate if useCubicHermiteThroughVesselWall else nodetemplateLinearS3)
    #            cache.setNode(node)
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vox [v][n3][n1])
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vod1[v][n3][n1])
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vod2[v][n3][n1])
    #            if useCubicHermiteThroughVesselWall:
    #                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vod3[v][n3][n1])
    #            #vNodeId.append(nodeIdentifier)
    #            nodeIdentifier += 1

    # get identifiers of nodes around each vessel at ostium end
    mvNodeId = [None] * vesselsCount
    for v in range(vesselsCount):
        if vesselsCount == 1:
            mvNodeId[v] = oNodeId
        else:
            iv = max(0, v - 1)
            mvNodeId[v] = [None for n3 in range(elementsCountThroughWall + 1)]
            oa = elementsCountAroundMid - iv * oinc
            ob = elementsCountAroundMid + nodesCountFreeEnd - 1 + iv * oinc
            for n3 in range(elementsCountThroughWall + 1):
                if v == 0:  # first end vessel
                    mvNodeId[v][n3] = oNodeId[n3][oa:ob + 1] + \
                                      (list(reversed(xNodeId[iv][n3])) if (v == 0) else xNodeId[iv][n3])
                elif v == (vesselsCount - 1):  # last end vessels
                    mvNodeId[v][n3] = oNodeId[n3][ob:] + [oNodeId[n3][0]] + \
                                      (list(reversed(xNodeId[iv][n3])) if (v == 0) else xNodeId[iv][n3])
                else:  # mid vessels
                    mvNodeId[v][n3] = oNodeId[n3][oa - oinc:oa + 1] + xNodeId[iv][n3] + \
                                      oNodeId[n3][ob:ob + oinc + 1] + list(reversed(xNodeId[iv + 1][n3]))

    #################
    # Create elements
    #################

    # tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    # tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    #
    # eft = tricubichermite.createEftBasic()
    # elementtemplate = mesh.createElementtemplate()
    # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    # elementtemplate.defineField(coordinates, -1, eft)
    #
    # elementtemplateX = mesh.createElementtemplate()
    # elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    for v in range(vesselsCount):
        if vesselMeshGroups or ostiumMeshGroups:
            rowMeshGroups = []
            for i in range(elementsCountAlong):
                rowMeshGroups.append(copy.copy(vesselMeshGroups[v]) if vesselMeshGroups else [])
        else:
            rowMeshGroups = None
        if isOutlet:
            startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap = \
                mvPointsx[v], mvPointsd1[v], mvPointsd2[v], mvPointsd3[v], mvNodeId[v], mvDerivativesMap[v]
            endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap = \
                vox[v], vod1[v], vod2[v], vod3[v] if useCubicHermiteThroughVesselWall else None, None, None
            # reverse order of nodes around:
            for px in [startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
                        endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap]:
                if px:
                    for n3 in range(elementsCountThroughWall + 1):
                        px[n3] = [px[n3][0]] + px[n3][len(px[n3]) - 1:0:-1]
            if vesselsCount > 1:
                # must switch in and out xi1 maps around corners in startDerivativesMap
                for n3 in range(elementsCountThroughWall + 1):
                    for n1 in range(elementsCountsAroundVessels[v]):
                        derivativesMap = startDerivativesMap[n3][n1]
                        if len(derivativesMap) == 4:
                            startDerivativesMap[n3][n1] = derivativesMap[3], derivativesMap[1], \
                                                          derivativesMap[2], derivativesMap[0]
            if ostiumMeshGroups:
                rowMeshGroups[0] += ostiumMeshGroups
        else:
            startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap = \
                vox[v], vod1[v], vod2[v], vod3[v] if useCubicHermiteThroughVesselWall else None, None, None
            endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap = \
                mvPointsx[v], mvPointsd1[v], mvPointsd2[v], mvPointsd3[v], mvNodeId[v], mvDerivativesMap[v]
            if ostiumMeshGroups:
                rowMeshGroups[-1] += ostiumMeshGroups

        # print('endPointsx ', endPointsx )
        # print('endPointsd1', endPointsd1)
        # print('endPointsd2', endPointsd2)
        # print('endPointsd3', endPointsd3)
        # print('endNodeId', endNodeId)
        # print('endDerivativesMap', endDerivativesMap)

        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
            endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
            forceMidLinearXi3=not useCubicHermiteThroughVesselWall,
            maxStartThickness=vesselWallThickness, maxEndThickness=vesselWallThickness,
            elementsCountRadial=elementsCountAlong, meshGroups=rowMeshGroups, wallAnnotationGroups=wallAnnotationGroups,
            coordinates=coordinates)

    fm.endChange()
    return nodeIdentifier, elementIdentifier, (ox, od1, od2, od3, oNodeId, oPositions)
