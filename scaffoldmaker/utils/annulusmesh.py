'''
Utility functions for generating annulus mesh between start and end loops of points.
'''
from __future__ import division
import copy
import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

def derivativeSignsToExpressionTerms(valueLabels, signs):
    '''
    Return remap expression terms for summing derivative[i]*sign[i]
    :param valueLabels: List of node value labels to possibly include.
    :param signs: List of 1 (no scaling), -1 (scale by scale factor 1) or 0 (no term).
    '''
    expressionTerms = []
    for i in range(len(valueLabels)):
        if signs[i] is 1:
            expressionTerms.append( ( valueLabels[i], [] ) )
        elif signs[i] is -1:
            expressionTerms.append( ( valueLabels[i], [1] ) )
    return expressionTerms

def createAnnulusMesh3d(nodes, mesh, nextNodeIdentifier, nextElementIdentifier,
    startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
    endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
    forceStartLinearXi3 = False, forceMidLinearXi3 = False, forceEndLinearXi3 = False,
    maxStartThickness = None, maxEndThickness = None, useCrossDerivatives = False,
    elementsCountRadial = 1, meshGroups = []):
    '''
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
        List array[n3][n1][c] or start/point coordinates and derivatives. To linearise through the wall, pass None to d3.
        If both ends are linear through the wall, interior points are linear through the wall.
    :param startNodeId, endNodeId: List array [n3][n1] of existing node identifiers to use at start/end. Pass None for
        argument if no nodes are specified at end. These arguments are 'all or nothing'.
    :param startDerivativesMap, endDerivativesMap: List array[n3][n1] of mappings for d/dxi1, d/dxi2, d/dxi3 at start/end of form:
        ( (1, -1, 0), (1, 0, 0), None ) where the first tuple means d/dxi1 = d/ds1 - d/ds2. Only 0, 1 and -1 may be used.
        None means use default e.g. d/dxi2 = d/ds2.
        Pass None for the entire argument to use the defaults d/dxi1 = d/ds1, d/dxi2 = d/ds2, d/dxi3 = d/ds3.
        Pass a 4th mapping to apply to d/dxi1 on other side of node; if not supplied first mapping applies both sides.
    :param nodetemplate: Full tricubic Hermite node template, can omit cross derivatives.
    :param forceStartLinearXi3, forceMidLinearXi3, forceEndLinearXi3: Force start, middle or
        end elements to be linear through the wall, even if d3 is supplied at either end.
        Can only use forceMidLinearXi3 only if at least one end is linear in d3.
    :param maxStartThickness, maxEndThickness: Optional maximum override on start/end thicknesses.
    :param useCrossDerivatives: May only be True if no derivatives maps are in use.
    :param elementsCountRadial: Optional number of elements in radial direction between start and end.
    :param meshGroups:  Optional list of Zinc MeshGroup for adding new elements to.
    :return: Final values of nextNodeIdentifier, nextElementIdentifier
    '''
    assert (elementsCountRadial >= 1), 'createAnnulusMesh3d:  Invalid number of radial elements'
    startLinearXi3 = (not startPointsd3) or forceStartLinearXi3
    endLinearXi3 = (not endPointsd3) or forceEndLinearXi3
    midLinearXi3 = (startLinearXi3 and endLinearXi3) or ((startLinearXi3 or endLinearXi3) and forceMidLinearXi3)
    # get list whether each row of nodes in elements is linear in Xi3
    # this is for element use; start/end nodes may have d3 even if element is linear
    rowLinearXi3 = [ startLinearXi3 ] + [ midLinearXi3 ]*(elementsCountRadial - 1) + [ endLinearXi3 ]
    assert (not useCrossDerivatives) or ((not startDerivativesMap) and (not endDerivativesMap)), \
        'createAnnulusMesh3d:  Cannot use cross derivatives with derivatives map'
    elementsCountWall = 1
    nodesCountWall = elementsCountWall + 1
    assert (len(startPointsx) == nodesCountWall) and (len(startPointsd1) == nodesCountWall) and (len(startPointsd2) == nodesCountWall) and \
        (startLinearXi3 or (len(startPointsd3) == nodesCountWall)) and \
        (len(endPointsx) == nodesCountWall) and (len(endPointsd1) == nodesCountWall) and (len(endPointsd2) == nodesCountWall) and \
        (endLinearXi3 or (len(endPointsd3) == nodesCountWall)) and \
        ((startNodeId is None) or (len(startNodeId) == nodesCountWall)) and \
        ((endNodeId is None) or (len(endNodeId) == nodesCountWall)) and \
        ((startDerivativesMap is None) or (len(startDerivativesMap) == nodesCountWall)) and \
        ((endDerivativesMap is None) or (len(endDerivativesMap) == nodesCountWall)), \
        'createAnnulusMesh3d:  Mismatch in number of layers through wall'
    elementsCountAround = nodesCountAround = len(startPointsx[0])
    assert (nodesCountAround > 1), 'createAnnulusMesh3d:  Invalid number of points/nodes around annulus'
    for n3 in range(nodesCountWall):
        assert (len(startPointsx[n3]) == nodesCountAround) and (len(startPointsd1[n3]) == nodesCountAround) and (len(startPointsd2[n3]) == nodesCountAround) and \
            (startLinearXi3 or (len(startPointsd3[n3]) == nodesCountAround)) and \
            (len(endPointsx[n3]) == nodesCountAround) and (len(endPointsd1[n3]) == nodesCountAround) and (len(endPointsd2[n3]) == nodesCountAround) and \
            (endLinearXi3 or (len(endPointsd3[n3]) == nodesCountAround)) and \
            ((startNodeId is None) or (len(startNodeId[n3]) == nodesCountAround)) and \
            ((endNodeId is None) or (len(endNodeId[n3]) == nodesCountAround)) and \
            ((startDerivativesMap is None) or (len(startDerivativesMap[n3]) == nodesCountAround)) and \
            ((endDerivativesMap is None) or (len(endDerivativesMap[n3]) == nodesCountAround)), \
            'createAnnulusMesh3d:  Mismatch in number of points/nodes in layers through wall'

    fm = mesh.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    coordinates = zinc_utils.getOrCreateCoordinateField(fm)

    # Build arrays of points from start to end
    px  = [ [], [] ]
    pd1 = [ [], [] ]
    pd2 = [ [], [] ]
    pd3 = [ [], [] ]
    for n3 in range(2):
        px [n3] = [ startPointsx [n3], endPointsx [n3] ]
        pd1[n3] = [ startPointsd1[n3], endPointsd1[n3] ]
        pd2[n3] = [ startPointsd2[n3], endPointsd2[n3] ]
        pd3[n3] = [ startPointsd3[n3] if (startPointsd3 is not None) else None, \
                    endPointsd3[n3] if (endPointsd3 is not None) else None ]
    if elementsCountRadial > 1:
        # add in-between points
        startPointsd = [ startPointsd1, startPointsd2, startPointsd3 ]
        startPointsdslimit = 2 if (startPointsd3 is None) else 3
        endPointsd = [ endPointsd1, endPointsd2, endPointsd3 ]
        endPointsdslimit = 2 if (endPointsd3 is None) else 3
        for n3 in range(2):
            for n2 in range(1, elementsCountRadial):
                px [n3].insert(n2, [ None ]*nodesCountAround)
                pd1[n3].insert(n2, [ None ]*nodesCountAround)
                pd2[n3].insert(n2, [ None ]*nodesCountAround)
                pd3[n3].insert(n2, None if midLinearXi3 else [ None ]*nodesCountAround)
        # compute on outside / n3 = 1, then map to inside using thickness
        thicknesses = []
        thicknesses.append([ vector.magnitude([ (startPointsx[1][n1][c] - startPointsx[0][n1][c]) for c in range(3) ]) for n1 in range(nodesCountAround) ])
        if maxStartThickness:
            for n1 in range(nodesCountAround):
                thicknesses[0][n1] = min(thicknesses[0][n1], maxStartThickness)
        for n2 in range(1, elementsCountRadial):
            thicknesses.append([ None ]*nodesCountAround)
        thicknesses.append([ vector.magnitude([ (endPointsx[1][n1][c] - endPointsx[0][n1][c]) for c in range(3) ]) for n1 in range(nodesCountAround) ])
        if maxEndThickness:
            for n1 in range(nodesCountAround):
                thicknesses[-1][n1] = min(thicknesses[-1][n1], maxEndThickness)
        n3 == 1
        for n1 in range(nodesCountAround):
            ax  = startPointsx [n3][n1]
            if (startDerivativesMap is None) or (startDerivativesMap[n3][n1][0] is None):
                ad1 = startPointsd1[n3][n1]
            else:
                derivativesMap = startDerivativesMap[n3][n1][0]
                ad1 = [ 0.0, 0.0, 0.0 ]
                for ds in range(startPointsdslimit):
                    if derivativesMap[ds] != 0.0:
                        for c in range(3):
                            ad1[c] += derivativesMap[ds]*startPointsd[ds][n3][n1][c]
                if len(startDerivativesMap[n3][n1]) > 3:
                    # average with d1 map for other side
                    derivativesMap = startDerivativesMap[n3][n1][3]
                    ad1 = [ 0.5*d for d in ad1 ]
                    if not derivativesMap:
                        for c in range(3):
                            ad1[c] += 0.5*startPointsd[0][n3][n1][c]
                    else:
                        for ds in range(startPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    ad1[c] += 0.5*derivativesMap[ds]*startPointsd[ds][n3][n1][c]
            if (startDerivativesMap is None) or (startDerivativesMap[n3][n1][1] is None):
                ad2 = startPointsd2[n3][n1]
            else:
                derivativesMap = startDerivativesMap[n3][n1][1]
                ad2 = [ 0.0, 0.0, 0.0 ]
                for ds in range(startPointsdslimit):
                    if derivativesMap[ds] != 0.0:
                        for c in range(3):
                            ad2[c] += derivativesMap[ds]*startPointsd[ds][n3][n1][c]

            bx  = endPointsx [n3][n1]
            if (endDerivativesMap is None) or (endDerivativesMap[n3][n1][0] is None):
                bd1 = endPointsd1[n3][n1]
            else:
                derivativesMap = endDerivativesMap[n3][n1][0]
                bd1 = [ 0.0, 0.0, 0.0 ]
                for ds in range(endPointsdslimit):
                    if derivativesMap[ds] != 0.0:
                        for c in range(3):
                            bd1[c] += derivativesMap[ds]*endPointsd[ds][n3][n1][c]
                if len(endDerivativesMap[n3][n1]) > 3:
                    # average with d1 map for other side
                    derivativesMap = endDerivativesMap[n3][n1][3]
                    bd1 = [ 0.5*d for d in bd1 ]
                    if not derivativesMap:
                        for c in range(3):
                            bd1[c] += 0.5*endPointsd[0][n3][n1][c]
                    else:
                        for ds in range(endPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    bd1[c] += 0.5*derivativesMap[ds]*endPointsd[ds][n3][n1][c]
            if (endDerivativesMap is None) or (endDerivativesMap[n3][n1][1] is None):
                bd2 = endPointsd2[n3][n1]
            else:
                derivativesMap = endDerivativesMap[n3][n1][1]
                bd2 = [ 0.0, 0.0, 0.0 ]
                for ds in range(endPointsdslimit):
                    if derivativesMap[ds] != 0.0:
                        for c in range(3):
                            bd2[c] += derivativesMap[ds]*endPointsd[ds][n3][n1][c]

            # scaling end derivatives to arc length gives even curvature along the curve
            arcLength = interp.computeCubicHermiteArcLength(ax, ad2, bx, bd2, rescaleDerivatives = False)
            scaledDerivatives = [ vector.setMagnitude(d2, arcLength) for d2 in [ ad2, bd2 ]]
            mx, md2, me, mxi = interp.sampleCubicHermiteCurvesSmooth([ ax, bx ], scaledDerivatives, elementsCountRadial,
                derivativeMagnitudeStart = vector.magnitude(ad2),
                derivativeMagnitudeEnd = vector.magnitude(bd2))[0:4]
            md1 = interp.interpolateSampleLinear([ ad1, bd1 ], me, mxi)
            thi = interp.interpolateSampleLinear([ thicknesses[0][n1], thicknesses[-1][n1] ], me, mxi)
            #md2 = interp.smoothCubicHermiteDerivativesLine(mx, md2, fixStartDerivative = True, fixEndDerivative = True)
            for n2 in range(1, elementsCountRadial):
                px [n3][n2][n1] = mx [n2]
                pd1[n3][n2][n1] = md1[n2]
                pd2[n3][n2][n1] = md2[n2]
                thicknesses[n2][n1] = thi[n2]

        # now get inner positions from normal and thickness, derivatives from curvature
        for n2 in range(1, elementsCountRadial):
            # first smooth derivative 1 around outer loop
            pd1[1][n2] = interp.smoothCubicHermiteDerivativesLoop(px[1][n2], pd1[1][n2], magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)

            for n1 in range(nodesCountAround):
                normal = vector.normalise(vector.crossproduct3(pd1[1][n2][n1], pd2[1][n2][n1]))
                thickness = thicknesses[n2][n1]
                d3 = [ d*thickness for d in normal ]
                px [0][n2][n1] = [ (px [1][n2][n1][c] - d3[c]) for c in range(3) ]
                # calculate inner d1 from curvature around
                n1m = n1 - 1
                n1p = (n1 + 1)%nodesCountAround
                curvature = 0.5*(
                    interp.getCubicHermiteCurvature(px[1][n2][n1m], pd1[1][n2][n1m], px[1][n2][n1 ], pd1[1][n2][n1 ], normal, 1.0) +
                    interp.getCubicHermiteCurvature(px[1][n2][n1 ], pd1[1][n2][n1 ], px[1][n2][n1p], pd1[1][n2][n1p], normal, 0.0))
                factor = 1.0 + curvature*thickness
                pd1[0][n2][n1] = [ factor*d for d in pd1[1][n2][n1] ]
                # calculate inner d2 from curvature radially
                n2m = n2 - 1
                n2p = n2 + 1
                curvature = 0.5*(
                    interp.getCubicHermiteCurvature(px[1][n2m][n1], pd2[1][n2m][n1], px[1][n2 ][n1], pd2[1][n2 ][n1], normal, 1.0) +
                    interp.getCubicHermiteCurvature(px[1][n2 ][n1], pd2[1][n2 ][n1], px[1][n2p][n1], pd2[1][n2p][n1], normal, 0.0))
                factor = 1.0 + curvature*thickness
                pd2[0][n2][n1] = [ factor*d for d in pd2[1][n2][n1] ]
                if not midLinearXi3:
                    pd3[0][n2][n1] = pd3[1][n2][n1] = d3

            # smooth derivative 1 around inner loop
            pd1[0][n2] = interp.smoothCubicHermiteDerivativesLoop(px[0][n2], pd1[0][n2], magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)

        for n3 in range(0, 1):  # was (0, nodesCountWall)
            # smooth derivative 2 radially/along annulus
            for n1 in range(nodesCountAround):
                sd2 = interp.smoothCubicHermiteDerivativesLine(
                    [ px [n3][n2][n1] for n2 in range(elementsCountRadial + 1) ],
                    [ pd2[n3][n2][n1] for n2 in range(elementsCountRadial + 1) ],
                    fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True,
                    magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
                for n2 in range(elementsCountRadial + 1):
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
    nodeId = [ [], [] ]
    for n3 in range(2):
        for n2 in range(elementsCountRadial + 1):
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

    for e2 in range(elementsCountRadial):
        nonlinearXi3 = (not rowLinearXi3[e2]) or (not rowLinearXi3[e2 + 1])
        eftFactory = tricubichermite if nonlinearXi3 else bicubichermitelinear
        eftStandard = eftFactory.createEftBasic()
        elementtemplateStandard.defineField(coordinates, -1, eftStandard)
        mapStartDerivatives = (e2 == 0) and (startDerivativesMap is not None)
        mapStartLinearDerivativeXi3 = nonlinearXi3 and rowLinearXi3[e2]
        mapEndDerivatives = (e2 == (elementsCountRadial - 1)) and (endDerivativesMap is not None)
        mapEndLinearDerivativeXi3 = nonlinearXi3 and rowLinearXi3[e2 + 1]
        mapDerivatives = mapStartDerivatives or mapStartLinearDerivativeXi3 or mapEndDerivatives or mapEndLinearDerivativeXi3
        for e1 in range(elementsCountAround):
            en = (e1 + 1)%elementsCountAround
            nids = [ nodeId[0][e2][e1], nodeId[0][e2][en], nodeId[0][e2 + 1][e1], nodeId[0][e2 + 1][en],
                     nodeId[1][e2][e1], nodeId[1][e2][en], nodeId[1][e2 + 1][e1], nodeId[1][e2 + 1][en] ]

            if mapDerivatives:
                eft1 = eftFactory.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                if mapStartLinearDerivativeXi3:
                    eftFactory.setEftLinearDerivative(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS3, 1, 5, 1)
                    eftFactory.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                if mapStartDerivatives:
                    for i in range(2):
                        lns = [ 1, 5 ] if (i == 0) else [ 2, 6 ]
                        for n3 in range(2):
                            derivativesMap = startDerivativesMap[n3][e1] if (i == 0) else startDerivativesMap[n3][en]
                            # handle different d1 on each side of node
                            d1Map = derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else derivativesMap[3]
                            d2Map = derivativesMap[1]
                            d3Map = derivativesMap[2]
                            # use temporary to safely swap DS1 and DS2:
                            ln = [ lns[n3] ]
                            if d1Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])
                            if d3Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS2DS3, [] ) ])
                            if d2Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d2Map))
                            if d1Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d1Map))
                            if d3Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d3Map))
                if mapEndLinearDerivativeXi3:
                    eftFactory.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                    eftFactory.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                if mapEndDerivatives:
                    for i in range(2):
                        lns = [ 3, 7 ] if (i == 0) else [ 4, 8 ]
                        for n3 in range(2):
                            derivativesMap = endDerivativesMap[n3][e1] if (i == 0) else endDerivativesMap[n3][en]
                            # handle different d1 on each side of node
                            d1Map = derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else derivativesMap[3]
                            d2Map = derivativesMap[1]
                            d3Map = derivativesMap[2]
                            # use temporary to safely swap DS1 and DS2:
                            ln = [ lns[n3] ]
                            if d1Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])
                            if d3Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS2DS3, [] ) ])
                            if d2Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d2Map))
                            if d1Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d1Map))
                            if d3Map is not None:
                                remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3, \
                                    derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d3Map))
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateX
            else:
                eft1 = eftStandard
                elementtemplate1 = elementtemplateStandard

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if mapDerivatives:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            #else:
            #    result3 = '-'
            #print('create element annulus', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

    fm.endChange()

    return nodeIdentifier, elementIdentifier
