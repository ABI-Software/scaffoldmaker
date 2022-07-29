"""
Utility functions for generating annulus mesh between start and end loops of points.
"""
from __future__ import division

import copy
from collections.abc import Sequence

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.node import Node
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite


def derivativeSignsToExpressionTerms(valueLabels, signs, scaleFactorIdx=None):
    """
    Return remap expression terms for summing derivative[i] * sign[i] * scaleFactor
    :param valueLabels: List of node value labels to possibly include.
    :param signs: List of 1 (no scaling), -1 (scale by scale factor 1) or 0 (no term).
    :param scaleFactorIdx: Optional index of local scale factor to scale all non-zero terms. Default None means no
    extra scaling.
    """
    expressionTerms = []
    for i in range(len(valueLabels)):
        if signs[i] == 1:
            expressionTerms.append((valueLabels[i], ([scaleFactorIdx] if scaleFactorIdx else [])))
        elif signs[i] == -1:
            expressionTerms.append((valueLabels[i], ([1, scaleFactorIdx] if scaleFactorIdx else [1])))
    return expressionTerms


def getMappedD1D2(gds, derivativesMaps):
    """
    Get vector combinations of d1In, d2In, d3In indicated by derivativesMap.
    :param gds: List of global d1, d2 and optionally d3.
    :param derivativesMaps: List over d1, d2, d3, and optionally d1b (for
    different d1 exiting global node) of list of 3 weights of gds,
    each limited to  -1.0, 0.0, or -1.0.
    :return: Effective d1, d2. Where d1 is around, d2 is radial.
    """
    dslimit = len(gds)
    if not (derivativesMaps and derivativesMaps[0]):
        d1 = gds[0]
    else:
        derivativesMap = derivativesMaps[0]
        d1 = [0.0, 0.0, 0.0]
        for ds in range(dslimit):
            if derivativesMap[ds] != 0.0:
                for c in range(3):
                    d1[c] += derivativesMap[ds] * gds[ds][c]
        if len(derivativesMaps) > 3:
            # average with d1 map for other side
            derivativesMap = derivativesMaps[3]
            d1 = [0.5 * d for d in d1]
            if not derivativesMap:
                for c in range(3):
                    d1[c] += 0.5 * gds[0][c]
            else:
                for ds in range(dslimit):
                    if derivativesMap[ds] != 0.0:
                        for c in range(3):
                            d1[c] += 0.5 * derivativesMap[ds] * gds[ds][c]
    if not (derivativesMaps and derivativesMaps[1]):
        d2 = gds[1]
    else:
        derivativesMap = derivativesMaps[1]
        d2 = [0.0, 0.0, 0.0]
        for ds in range(dslimit):
            if derivativesMap[ds] != 0.0:
                for c in range(3):
                    d2[c] += derivativesMap[ds] * gds[ds][c]
    return d1, d2


def createAnnulusMesh3d(nodes, mesh, nextNodeIdentifier, nextElementIdentifier, startPointsx, startPointsd1,
                        startPointsd2, startPointsd3, startNodeId, startDerivativesMap, endPointsx, endPointsd1,
                        endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
                        forceStartLinearXi3=False, forceMidLinearXi3=False, forceEndLinearXi3=False,
                        maxStartThickness=None, maxEndThickness=None, useCrossDerivatives=False,
                        elementsCountRadial=1, meshGroups=None, wallAnnotationGroups=None,
                        tracksurface=None, startProportions=None, endProportions=None,
                        rescaleStartDerivatives=False, rescaleEndDerivatives=False, sampleBlend=0.0,
                        fixMinimumStart=False, fixMinimumEnd=False, coordinates=None):
    """
    Create an annulus mesh from a loop of start points/nodes with specified derivative mappings to
    a loop of end points/nodes with specified derivative mappings.
    Derivative d3 is through the wall. Currently limited to single element layer through wall.
    Points/nodes order cycles fastest around the annulus, then through the wall.
    Note doesn't support cross derivatives.
    Arrays are indexed by n3 (node through wall, size 2), n2 (node along/radial), n1 (node around, variable size)
    and coordinate component c.
    :param nodes: The nodeset to create nodes in.
    :param mesh: The mesh to create elements in.
    :param nextNodeIdentifier, nextElementIdentifier: Next identifiers to use and increment.
    :param startPointsx, startPointsd1, startPointsd2, startPointsd3, endPointsx, endPointsd1, endPointsd2, endPointsd3:
        List array[n3][n1][c] or start/point coordinates and derivatives. To linearise through the wall, pass None to
        d3. If both ends are linear through the wall, interior points are linear through the wall.
    :param startNodeId, endNodeId: List array [n3][n1] of existing node identifiers to use at start/end. Pass None for
        argument if no nodes are specified at end. These arguments are 'all or nothing'.
    :param startDerivativesMap, endDerivativesMap: List array[n3][n1] of mappings for d/dxi1, d/dxi2, d/dxi3 at
        start/end of form:
        ( (1, -1, 0), (1, 0, 0), None ) where the first tuple means d/dxi1 = d/ds1 - d/ds2. Only 0, 1 and -1 may be
        used.
        None means use default e.g. d/dxi2 = d/ds2.
        Pass None for the entire argument to use the defaults d/dxi1 = d/ds1, d/dxi2 = d/ds2, d/dxi3 = d/ds3.
        Pass a 4th mapping to apply to d/dxi1 on other side of node; if not supplied first mapping applies both sides.
    :param forceStartLinearXi3, forceMidLinearXi3, forceEndLinearXi3: Force start, middle or
        end elements to be linear through the wall, even if d3 is supplied at either end.
        Can only use forceMidLinearXi3 only if at least one end is linear in d3.
    :param maxStartThickness, maxEndThickness: Optional maximum override on start/end thicknesses.
    :param useCrossDerivatives: May only be True if no derivatives maps are in use.
    :param elementsCountRadial: Optional number of elements in radial direction between start and end.
    :param meshGroups:  Optional sequence of Zinc MeshGroup for adding all new elements to, or a sequence of
    length elementsCountRadial containing sequences of mesh groups to add rows of radial elements to
    from start to end.
    :param wallAnnotationGroups: Annotation groups for adding all new elements to a sequence
    of groups to add to elements through wall.
    :param tracksurface: Description for outer surface representation used for creating annulus mesh. Provides
    information for creating radial nodes on annulus that sit on tracksurface. Need startProportions and endProportions
    to work.
    :param startProportions: Proportion around and along of startPoints on tracksurface. These vary with nodes
    around as for startPoints. Values only given for tracksurface for outer layer (xi3 == 1).
    :param endProportions: Proportion around and along of endPoints on track surface. These vary with nodes
    around as for endPoints. Values only given for tracksurface for outer layer (xi3 == 1).
    :param rescaleStartDerivatives, rescaleEndDerivatives: Optional flags to compute and multiply additional scale
    factors on start, end or both radial derivatives to fit arc length, needed if derivatives are of the wrong scale
    for the radial distances and the chosen elementsCountRadial. If either is True, derivatives and sampled radial
    nodes are spaced for a gradual change of derivative from that at the other end. If both are True, scaling is set to
    give even sampling and arclength derivatives.
    :param sampleBlend: Real value varying from 0.0 to 1.0 controlling weighting of start and end
    derivatives when interpolating extra points in-between, where 0.0 = sample with equal end derivatives,
    and 1.0 = proportional to current magnitudes, interpolated in between.
    :param fixMinimumStart, fixMinimumEnd: If either set (but not both) the shortest arc length is used to fix
    the interpolation derivative at that end, reducing stretching on long diagonals. Sample blend must be 0.0.
    :return: Final values of nextNodeIdentifier, nextElementIdentifier
    """
    assert (elementsCountRadial >= 1), 'createAnnulusMesh3d:  Invalid number of radial elements'
    startLinearXi3 = (not startPointsd3) or forceStartLinearXi3
    endLinearXi3 = (not endPointsd3) or forceEndLinearXi3
    midLinearXi3 = (startLinearXi3 and endLinearXi3) or ((startLinearXi3 or endLinearXi3) and forceMidLinearXi3)
    # get list whether each row of nodes in elements is linear in Xi3
    # this is for element use; start/end nodes may have d3 even if element is linear
    rowLinearXi3 = [startLinearXi3] + [midLinearXi3] * (elementsCountRadial - 1) + [endLinearXi3]
    assert (not useCrossDerivatives) or ((not startDerivativesMap) and (not endDerivativesMap)), \
        'createAnnulusMesh3d:  Cannot use cross derivatives with derivatives map'
    nodesCountWall = len(startPointsx)
    assert (len(startPointsd1) == nodesCountWall) and (len(startPointsd2) == nodesCountWall) and \
           (startLinearXi3 or (len(startPointsd3) == nodesCountWall)) and \
           (len(endPointsx) == nodesCountWall) and (len(endPointsd1) == nodesCountWall) and \
           (len(endPointsd2) == nodesCountWall) and (endLinearXi3 or (len(endPointsd3) == nodesCountWall)) and \
           ((startNodeId is None) or (len(startNodeId) == nodesCountWall)) and \
           ((endNodeId is None) or (len(endNodeId) == nodesCountWall)) and \
           ((startDerivativesMap is None) or (len(startDerivativesMap) == nodesCountWall)) and \
           ((endDerivativesMap is None) or (len(endDerivativesMap) == nodesCountWall)),\
        'createAnnulusMesh3d:  Mismatch in number of layers through wall'
    elementsCountAround = nodesCountAround = len(startPointsx[0])
    assert (nodesCountAround > 1), 'createAnnulusMesh3d:  Invalid number of points/nodes around annulus'
    for n3 in range(nodesCountWall):
        assert (len(startPointsx[n3]) == nodesCountAround) and (len(startPointsd1[n3]) == nodesCountAround) and \
               (len(startPointsd2[n3]) == nodesCountAround) and \
               (startLinearXi3 or (len(startPointsd3[n3]) == nodesCountAround)) and\
               (len(endPointsx[n3]) == nodesCountAround) and (len(endPointsd1[n3]) == nodesCountAround) and \
               (len(endPointsd2[n3]) == nodesCountAround) and \
               (endLinearXi3 or (len(endPointsd3[n3]) == nodesCountAround)) and \
               ((startNodeId is None) or (len(startNodeId[n3]) == nodesCountAround)) and\
               ((endNodeId is None) or (len(endNodeId[n3]) == nodesCountAround)) and \
               ((startDerivativesMap is None) or (len(startDerivativesMap[n3]) == nodesCountAround)) and \
               ((endDerivativesMap is None) or (len(endDerivativesMap[n3]) == nodesCountAround)), \
            'createAnnulusMesh3d:  Mismatch in number of points/nodes in layers through wall'
    assert not (fixMinimumStart and fixMinimumEnd and (sampleBlend > 0.0))
    rowMeshGroups = meshGroups
    if meshGroups:
        assert isinstance(meshGroups, Sequence), 'createAnnulusMesh3d:  Mesh groups is not a sequence'
        if (len(meshGroups) == 0) or (not isinstance(meshGroups[0], Sequence)):
            rowMeshGroups = [meshGroups] * elementsCountRadial
        else:
            assert len(meshGroups) == elementsCountRadial, \
                'createAnnulusMesh3d:  Length of meshGroups sequence does not equal elementsCountRadial'
    if wallAnnotationGroups:
        assert len(wallAnnotationGroups) == nodesCountWall - 1, \
            'createAnnulusMesh3d:  Length of wallAnnotationGroups sequence does not equal elementsCountThroughWall'
    if tracksurface:
        assert startProportions and endProportions, \
            'createAnnulusMesh3d: Missing start and/or end proportions for use with tracksurface'
        assert len(startProportions) == nodesCountAround, \
            'createAnnulusMesh3d: Length of startProportions does not equal nodesCountAround'
        assert len(endProportions) == nodesCountAround, \
            'createAnnulusMesh3d: Length of endProportions does not equal nodesCountAround'

    fm = mesh.getFieldmodule()
    if not coordinates:
        coordinates = findOrCreateFieldCoordinates(fm)
    fm.beginChange()
    cache = fm.createFieldcache()

    # Build arrays of points from start to end
    px = [[] for n3 in range(nodesCountWall)]
    pd1 = [[] for n3 in range(nodesCountWall)]
    pd2 = [[] for n3 in range(nodesCountWall)]
    pd3 = [[] for n3 in range(nodesCountWall)]

    # Find total wall thickness
    thicknessProportions = []
    thicknesses = []
    thicknesses.append([vector.magnitude([(startPointsx[nodesCountWall - 1][n1][c] - startPointsx[0][n1][c])
                                          for c in range(3)]) for n1 in range(nodesCountAround)])
    for n2 in range(1, elementsCountRadial):
        thicknesses.append([None] * nodesCountAround)
    thicknesses.append([vector.magnitude([(endPointsx[nodesCountWall - 1][n1][c] - endPointsx[0][n1][c])
                                          for c in range(3)]) for n1 in range(nodesCountAround)])

    for n3 in range(nodesCountWall):
        px[n3] = [startPointsx[n3], endPointsx[n3]]
        pd1[n3] = [startPointsd1[n3], endPointsd1[n3]]
        pd2[n3] = [startPointsd2[n3], endPointsd2[n3]]
        pd3[n3] = [startPointsd3[n3] if (startPointsd3 is not None) else None,
                   endPointsd3[n3] if (endPointsd3 is not None) else None]

        startThicknessList = \
            [vector.magnitude([(startPointsx[n3][n1][c] - startPointsx[n3 - (1 if n3 > 0 else 0)][n1][c])
                               for c in range(3)]) for n1 in range(len(startPointsx[n3]))]
        endThicknessList = \
            [vector.magnitude([(endPointsx[n3][n1][c] - endPointsx[n3 - (1 if n3 > 0 else 0)][n1][c])
                               for c in range(3)]) for n1 in range(len(endPointsx[n3]))]
        thicknessList = [startThicknessList, endThicknessList]  # thickness of each layer

        startThicknessProportions = [thicknessList[0][c] / thicknesses[0][c] for c in range(nodesCountAround)]
        endThicknessProportions = [thicknessList[1][c] / thicknesses[-1][c] for c in range(nodesCountAround)]
        thicknessProportions.append([startThicknessProportions, endThicknessProportions])

    if rescaleStartDerivatives:
        scaleFactorMapStart = [[] for n3 in range(nodesCountWall)]
    if rescaleEndDerivatives:
        scaleFactorMapEnd = [[] for n3 in range(nodesCountWall)]

    fixStartMagnitude = None
    fixEndMagnitude = None
    if fixMinimumStart or fixMinimumEnd:
        n3 = nodesCountWall - 1
        # get minimum derivative for shortest arc length to bias interpolation
        for n1 in range(nodesCountAround):
            ax = startPointsx[n3][n1]
            ad1, ad2 = getMappedD1D2([startPointsd1[n3][n1], startPointsd2[n3][n1]] +
                                     ([startPointsd3[n3][n1]] if startPointsd3 else []),
                                     startDerivativesMap[n3][n1] if startDerivativesMap else None)
            bx = endPointsx[n3][n1]
            bd1, bd2 = getMappedD1D2([endPointsd1[n3][n1], endPointsd2[n3][n1]] +
                                     ([endPointsd3[n3][n1]] if endPointsd3 else []),
                                     endDerivativesMap[n3][n1] if endDerivativesMap else None)
            bd2 = vector.setMagnitude(bd2, vector.magnitude(ad2))
            scaling = interp.computeCubicHermiteDerivativeScaling(ax, ad2, bx, bd2)
            if fixMinimumStart:
                mag = vector.magnitude(ad2) * scaling
                if (fixStartMagnitude is None) or (mag < fixStartMagnitude):
                    fixStartMagnitude = mag
            else:  # fixMinimumEnd:
                mag = vector.magnitude(bd2) * scaling
                if (fixEndMagnitude is None) or (mag < fixEndMagnitude):
                    fixEndMagnitude = mag

    # following code adds in-between points, but also handles rescaling for 1 radial element
    for n3 in range(nodesCountWall):
        for n2 in range(1, elementsCountRadial):
            px[n3].insert(n2, [None] * nodesCountAround)
            pd1[n3].insert(n2, [None] * nodesCountAround)
            pd2[n3].insert(n2, [None] * nodesCountAround)
            pd3[n3].insert(n2, None if midLinearXi3 else [None] * nodesCountAround)
            thicknessProportions[n3].insert(n2, [None] * nodesCountAround)

    if maxStartThickness:
        for n1 in range(nodesCountAround):
            thicknesses[0][n1] = min(thicknesses[0][n1], maxStartThickness)
    if maxEndThickness:
        for n1 in range(nodesCountAround):
            thicknesses[-1][n1] = min(thicknesses[-1][n1], maxEndThickness)
    n3 = nodesCountWall - 1
    for n1 in range(nodesCountAround):
        ax = startPointsx[n3][n1]
        ad1, ad2 = getMappedD1D2([startPointsd1[n3][n1], startPointsd2[n3][n1]] +
                                 ([startPointsd3[n3][n1]] if startPointsd3 else []),
                                 startDerivativesMap[n3][n1] if startDerivativesMap else None)
        bx = endPointsx[n3][n1]
        bd1, bd2 = getMappedD1D2([endPointsd1[n3][n1], endPointsd2[n3][n1]] +
                                 ([endPointsd3[n3][n1]] if endPointsd3 else []),
                                 endDerivativesMap[n3][n1] if endDerivativesMap else None)

        # sample between start and end points and derivatives
        # scaling end derivatives to arc length gives even curvature along the curve
        aMag = vector.magnitude(ad2)
        bMag = vector.magnitude(bd2)
        ad2mag = 0.5 * ((1.0 + sampleBlend) * aMag + (1.0 - sampleBlend) * bMag)
        ad2Scaled = vector.setMagnitude(ad2, ad2mag) if (aMag > 0.0) else [0.0, 0.0, 0.0]
        bd2mag = 0.5 * ((1.0 + sampleBlend) * bMag + (1.0 - sampleBlend) * aMag)
        bd2Scaled = vector.setMagnitude(bd2, bd2mag) if (bMag > 0.0) else [0.0, 0.0, 0.0]
        scaling = interp.computeCubicHermiteDerivativeScaling(ax, ad2Scaled, bx, bd2Scaled)
        ad2Scaled = [d * scaling for d in ad2Scaled]
        bd2Scaled = [d * scaling for d in bd2Scaled]
        derivativeMagnitudeStart = None if rescaleStartDerivatives else vector.magnitude(ad2)
        derivativeMagnitudeEnd = None if rescaleEndDerivatives else vector.magnitude(bd2)
        if fixStartMagnitude:
            derivativeMagnitudeStart = fixStartMagnitude / elementsCountRadial
        elif fixEndMagnitude:
            derivativeMagnitudeEnd = fixEndMagnitude / elementsCountRadial

        if tracksurface:
            mx, md2, md1, md3, mProportions = \
                tracksurface.createHermiteCurvePoints(startProportions[n1][0], startProportions[n1][1],
                                                      endProportions[n1][0], endProportions[n1][1], elementsCountRadial,
                                                      derivativeStart=[d / elementsCountRadial for d in ad2Scaled],
                                                      derivativeEnd=[d / elementsCountRadial for d in bd2Scaled])
            mx, md2, md1 = \
                tracksurface.resampleHermiteCurvePointsSmooth(mx, md2, md1, md3, mProportions,
                                                              derivativeMagnitudeStart, derivativeMagnitudeEnd)[0:3]
            # interpolate thicknesses using xi calculated from radial arclength distances to points
            arcLengthInsideToRadialPoint = \
                [0.0] + [interp.getCubicHermiteArcLength(mx[n2], md2[n2], mx[n2 + 1], md2[n2 + 1])
                         for n2 in range(elementsCountRadial)]
            arclengthInsideToOutside = sum(arcLengthInsideToRadialPoint)
            thi = []
            for n2 in range(elementsCountRadial + 1):
                xi2 = arcLengthInsideToRadialPoint[n2 - 1] / arclengthInsideToOutside
                thi.append(thicknesses[-1][n1] * xi2 + thicknesses[0][n1] * (1.0 - xi2))
            thiProportion = []
            for m3 in range(nodesCountWall):
                thiProportionRadial = []
                for n2 in range(elementsCountRadial + 1):
                    xi2 = arcLengthInsideToRadialPoint[n2 - 1] / arclengthInsideToOutside
                    thiProportionRadial.append(thicknessProportions[m3][-1][n1] * xi2 +
                                               thicknessProportions[m3][0][n1] * (1.0 - xi2))
                thiProportion.append(thiProportionRadial)
        else:
            mx, md2, me, mxi = interp.sampleCubicHermiteCurvesSmooth([ax, bx], [ad2Scaled, bd2Scaled],
                                                                     elementsCountRadial, derivativeMagnitudeStart,
                                                                     derivativeMagnitudeEnd)[0:4]
            md1 = interp.interpolateSampleLinear([ad1, bd1], me, mxi)
            thi = interp.interpolateSampleLinear([thicknesses[0][n1], thicknesses[-1][n1]], me, mxi)
            thiProportion = []
            for m3 in range(nodesCountWall):
                thiProportion.append(interp.interpolateSampleLinear([thicknessProportions[m3][0][n1],
                                                                     thicknessProportions[m3][-1][n1]], me, mxi))

        # set scalefactors if rescaling, make same on inside for now
        if rescaleStartDerivatives:
            scaleFactor = vector.magnitude(md2[0]) / vector.magnitude(ad2)
            scaleFactorMapStart[n3].append(scaleFactor)
        if rescaleEndDerivatives:
            scaleFactor = vector.magnitude(md2[-1]) / vector.magnitude(bd2)
            scaleFactorMapEnd[n3].append(scaleFactor)

        for n2 in range(1, elementsCountRadial):
            px[n3][n2][n1] = mx[n2]
            pd1[n3][n2][n1] = md1[n2]
            pd2[n3][n2][n1] = md2[n2]
            thicknesses[n2][n1] = thi[n2]
            for m3 in range(nodesCountWall):
                thicknessProportions[m3][n2][n1] = thiProportion[m3][n2]

    xi3List = [[[[] for n1 in range(nodesCountAround)] for n2 in range(elementsCountRadial + 1)] for n3 in
               range(nodesCountWall)]
    for n1 in range(nodesCountAround):
        for n2 in range(elementsCountRadial + 1):
            xi3 = 0.0
            for n3 in range(nodesCountWall):
                xi3 += thicknessProportions[n3][n2][n1]
                xi3List[n3][n2][n1] = xi3

    # now get inner positions from normal and thickness, derivatives from curvature
    for n2 in range(1, elementsCountRadial):
        # first smooth derivative 1 around outer loop
        pd1[-1][n2] = \
            interp.smoothCubicHermiteDerivativesLoop(px[-1][n2], pd1[-1][n2],
                                                     magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)

        for n3 in range(0, nodesCountWall - 1):
            for n1 in range(nodesCountAround):
                xi3 = 1 - xi3List[n3][n2][n1]
                normal = vector.normalise(vector.crossproduct3(pd1[-1][n2][n1], pd2[-1][n2][n1]))
                thickness = thicknesses[n2][n1] * xi3
                d3 = [d * thickness for d in normal]
                px[n3][n2][n1] = [(px[-1][n2][n1][c] - d3[c]) for c in range(3)]
                # calculate inner d1 from curvature around
                n1m = n1 - 1
                n1p = (n1 + 1) % nodesCountAround
                curvature = 0.5 * (
                        interp.getCubicHermiteCurvature(px[-1][n2][n1m], pd1[-1][n2][n1m],
                                                        px[-1][n2][n1], pd1[-1][n2][n1], normal, 1.0) +
                        interp.getCubicHermiteCurvature(px[-1][n2][n1], pd1[-1][n2][n1], px[-1][n2][n1p],
                                                        pd1[-1][n2][n1p], normal, 0.0))
                factor = 1.0 + curvature * thickness
                pd1[n3][n2][n1] = [factor * d for d in pd1[-1][n2][n1]]
                # calculate inner d2 from curvature radially
                n2m = n2 - 1
                n2p = n2 + 1
                curvature = 0.5 * (
                        interp.getCubicHermiteCurvature(px[-1][n2m][n1], pd2[-1][n2m][n1],
                                                        px[-1][n2][n1], pd2[-1][n2][n1], normal, 1.0) +
                        interp.getCubicHermiteCurvature(px[-1][n2][n1], pd2[-1][n2][n1],
                                                        px[-1][n2p][n1], pd2[-1][n2p][n1], normal, 0.0))
                factor = 1.0 + curvature * thickness
                pd2[n3][n2][n1] = [factor * d for d in pd2[-1][n2][n1]]
                d2Scaled = [factor * d for d in pd2[-1][n2][n1]]
                if vector.dotproduct(vector.normalise(pd2[-1][n2][n1]), vector.normalise(d2Scaled)) == -1:
                    pd2[n3][n2][n1] = [-factor * d for d in pd2[-1][n2][n1]]
                if not midLinearXi3:
                    pd3[n3][n2][n1] = pd3[-1][n2][n1] = \
                        [d * thicknesses[n2][n1] * thicknessProportions[n3 + 1][n2][n1] for d in normal]

            # smooth derivative 1 around inner loop
            pd1[n3][n2] = interp.smoothCubicHermiteDerivativesLoop(px[n3][n2], pd1[n3][n2],
                                                                   magnitudeScalingMode=interp.DerivativeScalingMode.
                                                                   HARMONIC_MEAN)

    for n3 in range(0, nodesCountWall):
        # smooth derivative 2 radially/along annulus
        for n1 in range(nodesCountAround):
            mx = [px[n3][n2][n1] for n2 in range(elementsCountRadial + 1)]
            md2 = [pd2[n3][n2][n1] for n2 in range(elementsCountRadial + 1)]
            # replace mapped start/end d2
            md2[0] = getMappedD1D2([startPointsd1[n3][n1], startPointsd2[n3][n1]] +
                                   ([startPointsd3[n3][n1]] if startPointsd3 else []),
                                   startDerivativesMap[n3][n1] if startDerivativesMap else None)[1]
            md2[-1] = getMappedD1D2([endPointsd1[n3][n1], endPointsd2[n3][n1]] +
                                    ([endPointsd3[n3][n1]] if endPointsd3 else []),
                                    endDerivativesMap[n3][n1] if endDerivativesMap else None)[1]

            sd2 = interp.smoothCubicHermiteDerivativesLine(mx, md2, fixAllDirections=True,
                                                           fixStartDerivative=not rescaleStartDerivatives,
                                                           fixStartDirection=rescaleStartDerivatives,
                                                           fixEndDerivative=not rescaleEndDerivatives,
                                                           fixEndDirection=rescaleEndDerivatives,
                                                           magnitudeScalingMode=interp.DerivativeScalingMode.
                                                           HARMONIC_MEAN)
            if rescaleStartDerivatives:
                scaleFactor = vector.magnitude(sd2[0]) / vector.magnitude(md2[0])
                scaleFactorMapStart[n3].append(scaleFactor)
            if rescaleEndDerivatives:
                scaleFactor = vector.magnitude(sd2[-1]) / vector.magnitude(md2[-1])
                scaleFactorMapEnd[n3].append(scaleFactor)

            for n2 in range(1, elementsCountRadial):
                pd2[n3][n2][n1] = sd2[n2]

    ##############
    # Create nodes
    ##############

    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    if useCrossDerivatives:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
    nodetemplateLinearS3 = nodes.createNodetemplate()
    nodetemplateLinearS3.defineField(coordinates)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    nodeIdentifier = nextNodeIdentifier
    nodeId = [[] for n3 in range(nodesCountWall)]
    for n2 in range(elementsCountRadial + 1):
        for n3 in range(nodesCountWall):
            if (n2 == 0) and (startNodeId is not None):
                rowNodeId = copy.deepcopy(startNodeId[n3])
            elif (n2 == elementsCountRadial) and (endNodeId is not None):
                rowNodeId = copy.deepcopy(endNodeId[n3])
            else:
                rowNodeId = []
                nodetemplate1 = nodetemplate if pd3[n3][n2] else nodetemplateLinearS3
                for n1 in range(nodesCountAround):
                    node = nodes.createNode(nodeIdentifier, nodetemplate1)
                    rowNodeId.append(nodeIdentifier)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, px[n3][n2][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n3][n2][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pd2[n3][n2][n1])
                    if pd3[n3][n2]:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, pd3[n3][n2][n1])
                    nodeIdentifier = nodeIdentifier + 1
            nodeId[n3].append(rowNodeId)

    #################
    # Create elements
    #################

    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)

    elementIdentifier = nextElementIdentifier

    elementtemplateStandard = mesh.createElementtemplate()
    elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateX = mesh.createElementtemplate()
    elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementsCountWall = nodesCountWall - 1

    for e2 in range(elementsCountRadial):
        nonlinearXi3 = (not rowLinearXi3[e2]) or (not rowLinearXi3[e2 + 1])
        eftFactory = tricubichermite if nonlinearXi3 else bicubichermitelinear
        eftStandard = eftFactory.createEftBasic()
        elementtemplateStandard.defineField(coordinates, -1, eftStandard)
        mapStartDerivatives = (e2 == 0) and (startDerivativesMap or rescaleStartDerivatives)
        mapStartLinearDerivativeXi3 = nonlinearXi3 and rowLinearXi3[e2]
        mapEndDerivatives = (e2 == (elementsCountRadial - 1)) and (endDerivativesMap or rescaleEndDerivatives)
        mapEndLinearDerivativeXi3 = nonlinearXi3 and rowLinearXi3[e2 + 1]
        mapDerivatives = mapStartDerivatives or mapStartLinearDerivativeXi3 or \
                         mapEndDerivatives or mapEndLinearDerivativeXi3
        for e3 in range(elementsCountWall):
            for e1 in range(elementsCountAround):
                en = (e1 + 1) % elementsCountAround
                nids = [nodeId[e3][e2][e1], nodeId[e3][e2][en], nodeId[e3][e2 + 1][e1], nodeId[e3][e2 + 1][en],
                        nodeId[e3 + 1][e2][e1], nodeId[e3 + 1][e2][en], nodeId[e3 + 1][e2 + 1][e1],
                        nodeId[e3 + 1][e2 + 1][en]]
                scaleFactors = []

                if mapDerivatives:
                    eft1 = eftFactory.createEftNoCrossDerivatives()
                    # work out if scaling by global -1
                    scaleMinus1 = mapStartLinearDerivativeXi3 or mapEndLinearDerivativeXi3
                    if (not scaleMinus1) and mapStartDerivatives and startDerivativesMap:
                        for n3 in range(2):
                            n3Idx = n3 + e3
                            # need to handle 3 or 4 maps (e1 uses last 3, en uses first 3)
                            for map in startDerivativesMap[n3Idx][e1][-3:]:
                                if map and (-1 in map):
                                    scaleMinus1 = True
                                    break
                            for map in startDerivativesMap[n3Idx][en][:3]:
                                if map and (-1 in map):
                                    scaleMinus1 = True
                                    break
                    if (not scaleMinus1) and mapEndDerivatives and endDerivativesMap:
                        for n3 in range(2):
                            n3Idx = n3 + e3
                            # need to handle 3 or 4 maps (e1 uses last 3, en uses first 3)
                            for map in endDerivativesMap[n3Idx][e1][-3:]:
                                if map and (-1 in map):
                                    scaleMinus1 = True
                                    break
                            for map in endDerivativesMap[n3Idx][en][:3]:
                                if map and (-1 in map):
                                    scaleMinus1 = True
                                    break
                    # make node scale factors vary fastest by local node varying across lower xi
                    nodeScaleFactorIds = []
                    for n3 in range(2):
                        n3Idx = n3 + e3
                        if mapStartDerivatives and rescaleStartDerivatives:
                            for i in range(2):
                                derivativesMap = (startDerivativesMap[n3Idx][e1][1] if (i == 0) else
                                                  startDerivativesMap[n3Idx][en][1]) if startDerivativesMap else None
                                nodeScaleFactorIds.append(getQuadrantID(derivativesMap if derivativesMap else
                                                                        (0, 1, 0)))
                        if mapEndDerivatives and rescaleEndDerivatives:
                            for i in range(2):
                                derivativesMap = (endDerivativesMap[n3Idx][e1][1] if (i == 0) else
                                                  endDerivativesMap[n3Idx][en][1]) if endDerivativesMap else None
                                nodeScaleFactorIds.append(getQuadrantID(derivativesMap if derivativesMap else
                                                                        (0, 1, 0)))
                    setEftScaleFactorIds(eft1, [1] if scaleMinus1 else [], nodeScaleFactorIds)
                    firstNodeScaleFactorIndex = 2 if scaleMinus1 else 1
                    firstStartNodeScaleFactorIndex = \
                        firstNodeScaleFactorIndex if (mapStartDerivatives and rescaleStartDerivatives) else None
                    firstEndNodeScaleFactorIndex = \
                        (firstNodeScaleFactorIndex + (2 if firstStartNodeScaleFactorIndex else 0)) \
                            if (mapEndDerivatives and rescaleEndDerivatives) else None
                    layerNodeScaleFactorIndexOffset = \
                        4 if (firstStartNodeScaleFactorIndex and firstEndNodeScaleFactorIndex) else 2
                    if scaleMinus1:
                        scaleFactors.append(-1.0)
                    for n3 in range(2):
                        n3Idx = n3 + e3
                        if firstStartNodeScaleFactorIndex:
                            scaleFactors.append(scaleFactorMapStart[n3Idx][e1])
                            scaleFactors.append(scaleFactorMapStart[n3Idx][en])
                        if firstEndNodeScaleFactorIndex:
                            scaleFactors.append(scaleFactorMapEnd[n3Idx][e1])
                            scaleFactors.append(scaleFactorMapEnd[n3Idx][en])

                    if mapStartLinearDerivativeXi3:
                        eftFactory.setEftLinearDerivative2(eft1, [1, 5, 2, 6], Node.VALUE_LABEL_D_DS3,
                                                           [Node.VALUE_LABEL_D2_DS1DS3])
                    if mapStartDerivatives:
                        for i in range(2):
                            lns = [1, 5] if (i == 0) else [2, 6]
                            for n3 in range(2):
                                n3Idx = n3 + e3
                                derivativesMap = \
                                    (startDerivativesMap[n3Idx][e1] if (i == 0) else startDerivativesMap[n3Idx][en]) \
                                        if startDerivativesMap else (None, None, None)
                                # handle different d1 on each side of node
                                d1Map = \
                                    derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else derivativesMap[3]
                                d2Map = derivativesMap[1] if derivativesMap[1] else (0, 1, 0)
                                d3Map = derivativesMap[2]
                                # use temporary to safely swap DS1 and DS2:
                                ln = [lns[n3]]
                                if d1Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                if d3Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS2DS3, [])])
                                if d2Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d2Map,
                                                                                            (firstStartNodeScaleFactorIndex + i + n3 * layerNodeScaleFactorIndexOffset) if rescaleStartDerivatives else None))
                                if d1Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d1Map))
                                if d3Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d3Map))
                    if mapEndLinearDerivativeXi3:
                        eftFactory.setEftLinearDerivative2(eft1, [3, 7, 4, 8], Node.VALUE_LABEL_D_DS3,
                                                           [Node.VALUE_LABEL_D2_DS1DS3])
                    if mapEndDerivatives:
                        for i in range(2):
                            lns = [3, 7] if (i == 0) else [4, 8]
                            for n3 in range(2):
                                n3Idx = n3 + e3
                                derivativesMap = \
                                    (endDerivativesMap[n3Idx][e1] if (i == 0) else endDerivativesMap[n3Idx][en]) \
                                        if endDerivativesMap else (None, None, None)
                                # handle different d1 on each side of node
                                d1Map = derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else \
                                    derivativesMap[3]
                                d2Map = derivativesMap[1] if derivativesMap[1] else (0, 1, 0)
                                d3Map = derivativesMap[2]

                                # use temporary to safely swap DS1 and DS2:
                                ln = [lns[n3]]
                                if d1Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                if d3Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS2DS3, [])])
                                if d2Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d2Map,
                                                                                            (firstEndNodeScaleFactorIndex + i + n3 * layerNodeScaleFactorIndexOffset) if rescaleEndDerivatives else None))
                                if d1Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d1Map))
                                if d3Map:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3,
                                                           derivativeSignsToExpressionTerms((Node.VALUE_LABEL_D_DS1,
                                                                                             Node.VALUE_LABEL_D_DS2,
                                                                                             Node.VALUE_LABEL_D_DS3),
                                                                                            d3Map))

                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                else:
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scaleFactors:
                    result3 = element.setScaleFactors(eft1, scaleFactors)
                # print('create element annulus', element.isValid(), elementIdentifier, eft1.validate(),
                #       result2, result3 if scaleFactors else None, nids)
                elementIdentifier += 1

                if rowMeshGroups:
                    for meshGroup in rowMeshGroups[e2]:
                        meshGroup.addElement(element)

                if wallAnnotationGroups:
                    for annotationGroup in wallAnnotationGroups[e3]:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

    fm.endChange()
    return nodeIdentifier, elementIdentifier


def getQuadrantID(d):
    """
    Returns a scale factor ID based on direction of the derivative. Index starts
    at 1 in direction (1, 0, 0) and progresses in anti-clockwise direction.
    :param d: directional derivative
    :return: scale factor ID
    """
    maps = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0), (-1, 0, 0), (-1, -1, 0), (0, -1, 0), (1, -1, 0),  # d3 = 0
            (1, 0, 1), (1, 1, 1), (0, 1, 1), (-1, 1, 1), (-1, 0, 1), (-1, -1, 1), (0, -1, 1), (1, -1, 1), (0, 0, 1),  # d3 = 1
            (1, 0, -1), (1, 1, -1), (0, 1, -1), (-1, 1, -1), (-1, 0, -1), (-1, -1, -1), (0, -1, -1), (1, -1, -1),
            (0, 0, -1)]  # d3 = -1
    for i in range(len(maps)):
        if d == maps[i]:
            return i + 1
