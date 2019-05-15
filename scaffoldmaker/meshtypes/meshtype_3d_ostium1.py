"""
Generates a single or double/common ostium, where one or more vessels enters a chamber.
"""

from __future__ import division
import copy
import math
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import createCirclePoints, getCircleProjectionAxes
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_ostium1(Scaffold_base):
    '''
    Generates a 3-D single or double/common ostium inlet or outlet.
    '''
    @staticmethod
    def getName():
        return '3D Ostium 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of vessels' : 2,
            'Number of elements around vessel' : 8,
            'Number of elements across multiple' : 2,
            'Number of elements along' : 1,
            'Number of elements through wall' : 1,  # not implemented for > 1
            'Unit scale' : 1.0,
            'Inlet' : True,
            'Ostium diameter' : 0.5,
            'Ostium length' : 0.2,
            'Ostium vessel spacing' : 0.4,
            'Ostium wall thickness' : 0.05,
            'Use linear through ostium wall' : False,
            'Vessel end length factor' : 1.0,
            'Vessel inner diameter' : 0.25,
            'Vessel wall thickness' : 0.025,
            'Vessel angle 1 degrees' : 0.0,
            'Vessel angle 1 spread degrees' : 30.0,
            'Vessel angle 2 degrees' : -20.0,
            'Use linear through vessel wall' : True,
            #'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of vessels',
            'Number of elements around vessel',
            'Number of elements across multiple',
            'Number of elements along',
            'Number of elements through wall',  # not implemented for > 1
            'Unit scale',
            'Inlet',  # GRC not implemented for outlet
            'Ostium diameter',
            'Ostium length',
            'Ostium vessel spacing',
            'Ostium wall thickness',
            'Use linear through ostium wall',
            'Vessel end length factor',
            'Vessel inner diameter',
            'Vessel wall thickness',
            'Vessel angle 1 degrees',
            'Vessel angle 1 spread degrees',
            'Vessel angle 2 degrees',
            'Use linear through vessel wall',
            #'Use cross derivatives',  # not implemented
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        dependentChanges = False
        for key in ['Number of vessels']:
            if options[key] < 1:
                options[key] = 1
            elif options[key] > 2:
                options[key] = 2
        for key in [
            'Number of elements across multiple',
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements through wall'] > 1:
            options['Number of elements through wall'] = 1
        if (options['Number of elements around vessel'] < 3) :
            options['Number of elements around vessel'] = 3
        if (options['Number of elements around vessel'] - options['Number of elements across multiple']) < 2:
            options['Number of elements around vessel'] = options['Number of elements across multiple'] + 2
            dependentChanges = True
        for key in [
            'Unit scale',
            'Ostium length',
            'Ostium vessel spacing',
            'Ostium wall thickness',
            'Vessel inner diameter',
            'Vessel wall thickness']:
           if options[key] < 0.0:
                options[key] = 0.0
        if options['Ostium diameter'] <= 0.0:
            options['Ostium diameter'] = 0.000001  # avoid division by zero
        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic/bicubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        unitScale = options['Unit scale']
        ostiumRadius = 0.5*unitScale*options['Ostium diameter']
        ostiumVesselSpacing = unitScale*options['Ostium vessel spacing']

        scale = unitScale*(ostiumRadius*2.0 + ostiumVesselSpacing)
        nx = [ [ -scale, -scale, 0.0 ], [ scale, -scale, 0.0 ], [ -scale, scale, 0.0 ], [ scale, scale, 0.0 ] ]
        nd1 = [ [ 2.0*scale, 0.0, 0.0 ] ]*4
        nd2 = [ [ 0.0, 2.0*scale, 0.0 ] ]*4
        trackSurface = TrackSurface(1, 1, nx, nd1, nd2)
        centrePosition = TrackSurfacePosition(0, 0, 0.5, 0.5)
        axis1 = [ 1.0, 0.0, 0.0 ]

        generateOstiumMesh(region, options, trackSurface, centrePosition, axis1)

        return

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
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


def generateOstiumMesh(region, options, trackSurface, centrePosition, axis1, startNodeIdentifier = 1, startElementIdentifier = 1,
        vesselMeshGroups = None, spacingMeshGroups = None):
    '''
    :param vesselMeshGroups: List (over number of vessels) of list of mesh groups to add vessel elements to.
    :param spacingMeshGroups: List of mesh groups to add spacing elements to, if any.
    :return: nextNodeIdentifier, nextElementIdentifier, Ostium points tuple
    (ox[n3][n1][c], od1[n3][n1][c], od2[n3][n1][c], od3[n3][n1][c], oNodeId[n3][n1]).
    '''
    vesselsCount = options['Number of vessels']
    elementsCountAroundVessel = options['Number of elements around vessel']
    elementsCountAcross = options['Number of elements across multiple']
    elementsCountAlong = options['Number of elements along']
    elementsCountThroughWall = options['Number of elements through wall']
    unitScale = options['Unit scale']

    isInlet = options['Inlet']
    ostiumRadius = 0.5*unitScale*options['Ostium diameter']
    ostiumLength = unitScale*options['Ostium length']
    ostiumVesselSpacing = unitScale*options['Ostium vessel spacing'] if (vesselsCount > 1) else 0.0
    halfOstiumVesselSpacing = 0.5*ostiumVesselSpacing
    ostiumWallThickness = unitScale*options['Ostium wall thickness']
    useCubicHermiteThroughOstiumWall = not(options['Use linear through ostium wall'])
    vesselEndDerivative = ostiumLength*options['Vessel end length factor']/elementsCountAlong
    vesselInnerRadius = 0.5*unitScale*options['Vessel inner diameter']
    vesselWallThickness = unitScale*options['Vessel wall thickness']
    vesselOuterRadius = vesselInnerRadius + vesselWallThickness
    vesselAngle1Radians = math.radians(options['Vessel angle 1 degrees'])
    vesselAngle1SpreadRadians = math.radians(options['Vessel angle 1 spread degrees'])
    vesselAngle2Radians = math.radians(options['Vessel angle 2 degrees'])
    useCubicHermiteThroughVesselWall = not(options['Use linear through vessel wall'])
    useCrossDerivatives = False  # options['Use cross derivatives']  # not implemented

    fm = region.getFieldmodule()
    fm.beginChange()
    coordinates = zinc_utils.getOrCreateCoordinateField(fm)
    cache = fm.createFieldcache()

    # track points in shape of ostium

    # get directions in plane of surface at centre:
    cx, cd1, cd2 = trackSurface.evaluateCoordinates(centrePosition, True)
    trackDirection1, trackDirection2, centreNormal = calculate_surface_axes(cd1, cd2, axis1)
    trackDirection1reverse = [ -d for d in trackDirection1 ]
    trackDirection2reverse = [ -d for d in trackDirection2 ]

    halfCircumference = math.pi*ostiumRadius
    circumference = 2.0*halfCircumference
    distance = 0.0
    if vesselsCount == 1:
        elementsCountAroundOstium = elementsCountAroundVessel
        distanceAround = circumference
        elementLengthAroundOstium = distanceAround/elementsCountAroundOstium
        if elementsCountAroundVessel % 2:
            distance += 0.5*elementLengthAroundOstium
        ostiumCentre1 = ostiumCentre2 = centrePosition
        ocx = [ cx ]
    else:
        elementsCountAroundOstium = 2*(elementsCountAroundVessel - elementsCountAcross)
        distanceAround = circumference + 2.0*ostiumVesselSpacing
        elementLengthAroundOstium = distanceAround/elementsCountAroundOstium
        ostiumCentre1 = trackSurface.trackVector(centrePosition, trackDirection1reverse, halfOstiumVesselSpacing)
        ostiumCentre2 = trackSurface.trackVector(centrePosition, trackDirection1, halfOstiumVesselSpacing)
        ocx = [ trackSurface.evaluateCoordinates(ostiumCentre1), trackSurface.evaluateCoordinates(ostiumCentre2) ]

    # coordinates around ostium
    ox = [ [], [] ]
    od1 = [ [], [] ]
    od2 = [ [], [] ]
    od3 = [ [], [] ]
    for n1 in range(elementsCountAroundOstium):
        if distance <= halfOstiumVesselSpacing:
            position = trackSurface.trackVector(centrePosition, trackDirection1reverse, distance)
            sideDirection = trackDirection2reverse
        elif distance < (halfOstiumVesselSpacing + halfCircumference):
            position = ostiumCentre1
            angleRadians = (distance - halfOstiumVesselSpacing)/ostiumRadius
            w1 = -math.sin(angleRadians)
            w2 = -math.cos(angleRadians)
            sideDirection = [ (w1*trackDirection1[c] + w2*trackDirection2[c]) for c in range(3) ]
        elif distance < (halfOstiumVesselSpacing + halfCircumference + ostiumVesselSpacing):
            position = trackSurface.trackVector(centrePosition, trackDirection1, distance - (ostiumVesselSpacing + halfCircumference))
            sideDirection = trackDirection2
        elif distance < (halfOstiumVesselSpacing + circumference + ostiumVesselSpacing):
            position = ostiumCentre2
            angleRadians = (distance - (halfOstiumVesselSpacing + halfCircumference + ostiumVesselSpacing))/ostiumRadius
            w1 = math.sin(angleRadians)
            w2 = math.cos(angleRadians)
            sideDirection = [ (w1*trackDirection1[c] + w2*trackDirection2[c]) for c in range(3) ]
        else:
            position = trackSurface.trackVector(centrePosition, trackDirection1reverse, distance - (circumference + 2.0*ostiumVesselSpacing))
            sideDirection = trackDirection2reverse
        position = trackSurface.trackVector(position, sideDirection, ostiumRadius)
        px, d1, d2 = trackSurface.evaluateCoordinates(position, True)
        pd2, pd1, pd3 = calculate_surface_axes(d1, d2, sideDirection)
        # get outer coordinates
        opx = px
        opd1 = vector.setMagnitude([ -d for d in pd1 ], elementLengthAroundOstium)
        opd2 = vector.setMagnitude(pd2, elementLengthAroundOstium)  # smoothed later
        opd3 = vector.setMagnitude(pd3, ostiumWallThickness)
        # set inner and outer coordinates (use copy to avoid references to same list later)
        ox [0].append([ (opx[c] - opd3[c]) for c in range(3) ])
        od1[0].append(copy.copy(opd1))
        od2[0].append(copy.copy(opd2))
        ox [1].append(opx)
        od1[1].append(opd1)
        od2[1].append(opd2)
        if useCubicHermiteThroughOstiumWall:
            od3[0].append(copy.copy(opd3))
            od3[1].append(opd3)
        distance += elementLengthAroundOstium
    for n3 in range(2):
        od1[n3] = interp.smoothCubicHermiteDerivativesLoop(ox[n3], od1[n3], fixAllDirections = True)

    if vesselsCount == 2:
        # coordinates across ostium, between multiple vessels
        xx = [ [], [] ]
        xd1 = [ [], [] ]
        xd2 = [ [], [] ]
        xd3 = [ [], [] ]
        elementLengthAcrossOstium = 2.0*ostiumRadius/elementsCountAcross
        for n2 in range(elementsCountAcross - 1):
            position2 = trackSurface.trackVector(centrePosition, trackDirection2, (n2 - 0.5*(elementsCountAcross - 2))*elementLengthAcrossOstium)
            px, d1, d2 = trackSurface.evaluateCoordinates(position2, True)
            pd2, pd1, pd3 = calculate_surface_axes(d1, d2, trackDirection2)
            # get outer coordinates
            opx = px
            opd1 = vector.setMagnitude([ -d for d in pd1 ], elementLengthAroundOstium)  # smoothed later
            opd2 = vector.setMagnitude(pd2, elementLengthAcrossOstium)  # smoothed later
            opd3 = vector.setMagnitude(pd3, ostiumWallThickness)
            # set inner and outer coordinates (use copy to avoid references to same list later)
            xx [0].append([ (opx[c] - opd3[c]) for c in range(3) ])
            xd1[0].append(copy.copy(opd1))
            xd2[0].append(copy.copy(opd2))
            xx [1].append(px)
            xd1[1].append(opd1)
            xd2[1].append(opd2)
            if useCubicHermiteThroughOstiumWall:
                xd3[0].append(copy.copy(opd3))
                xd3[1].append(opd3)

        # smooth d2 across
        oa = 0
        ob = elementsCountAroundOstium//2
        for n3 in range(2):
            nx = [ ox[n3][oa] ] + xx[n3] + [ ox[n3][ob] ]
            nd2 = [ [ -d for d in od2[n3][oa] ] ] + xd2[n3] + [ od2[n3][ob] ]
            nd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections = True)
            od2[n3][oa] = [ -d for d in nd2[0] ]
            xd2[n3] = nd2[1:elementsCountAcross]
            od2[n3][ob] = nd2[-1]

    # get positions of vessel end centres and rings
    cx, trackDirection1, trackDirection2, centreNormal
    vcx = []
    vcd1 = []
    vcd2 = []
    vcd3 = []
    vox = []
    vod1 = []
    vod2 = []
    vod3 = []
    radiansPerElementVessel = 2.0*math.pi/elementsCountAroundVessel
    for v in range(vesselsCount):
        useVesselAngleRadians = vesselAngle1Radians
        if vesselsCount > 1:
            useVesselAngleRadians += (v - 0.5)*vesselAngle1SpreadRadians
        vx, vd1, vd2, vd3 = getCircleProjectionAxes(ocx[v], trackDirection1, trackDirection2, centreNormal, useVesselAngleRadians, vesselAngle2Radians, ostiumLength)
        vd1 = [    vesselOuterRadius*d for d in vd1 ]
        vd2 = [   -vesselOuterRadius*d for d in vd2 ]
        vd3 = [ -vesselEndDerivative*d for d in vd3 ]
        vcx.append(vx)
        vcd1.append(vd1)
        vcd2.append(vd2)
        vcd3.append(vd3)

        vox.append([])
        vod1.append([])
        vod2.append([])
        vod3.append([])
        for n3 in range(2):
            radius = vesselInnerRadius if (n3 == 0) else vesselOuterRadius
            vAxis1 = vector.setMagnitude(vd1, radius)
            vAxis2 = vector.setMagnitude(vd2, radius)
            if vesselsCount == 1:
                startRadians = 0.5*math.pi
                if elementsCountAroundVessel % 2:
                    startRadians += 0.5*radiansPerElementVessel
            else:
                startRadians = 0.5*radiansPerElementVessel*elementsCountAcross
                if v == 1:
                    startRadians -= math.pi
            px, pd1 = createCirclePoints(vx, vAxis1, vAxis2, elementsCountAroundVessel, startRadians)
            vox [-1].append(px)
            vod1[-1].append(pd1)
            vod2[-1].append([ vd3 ]*elementsCountAroundVessel)
            if useCubicHermiteThroughVesselWall:
                vod3[-1].append([ vector.setMagnitude(vector.crossproduct3(d1, vd3), vesselWallThickness) for d1 in pd1 ])

    # calculate multiple vessel ostium node derivatives map
    nodesCountFree = elementsCountAroundVessel + 1 - elementsCountAcross
    mvPointsx = [ None ]*vesselsCount
    mvPointsd1 = [ None ]*vesselsCount
    mvPointsd2 = [ None ]*vesselsCount
    mvPointsd3 = [ None ]*vesselsCount
    mvDerivativesMap = [ None ]*vesselsCount
    for v in range(vesselsCount):
        if vesselsCount == 1:
            mvPointsx[v], mvPointsd1[v], mvPointsd2[v], mvPointsd3[v], mvDerivativesMap[v] = \
                ox, od1, od2, od3 if useCubicHermiteThroughOstiumWall else None, None
        else:
            os = 0 if (v == 0) else (nodesCountFree - 1)
            of = os + nodesCountFree
            mvPointsx [v] = []
            mvPointsd1[v] = []
            mvPointsd2[v] = []
            mvPointsd3[v] = [] if useCubicHermiteThroughOstiumWall else None
            mvDerivativesMap[v] = []
            for n3 in range(2):
                mvPointsx [v].append((ox [n3]*2)[os:of])
                mvPointsd1[v].append((od1[n3]*2)[os:of])
                mvPointsd2[v].append((od2[n3]*2)[os:of])
                if useCubicHermiteThroughOstiumWall:
                    mvPointsd3[v].append((od3[n3]*2)[os:of])
                mvDerivativesMap[v].append([])
                mvDerivativesMap[v][n3].append( ( (0, 1, 0), (-1, 1, 0), None, (1, 0, 0) ) )
                for i in range(nodesCountFree - 2):
                    mvDerivativesMap[v][n3].append( ( None, None, None ) )
                mvDerivativesMap[v][n3].append( ( (1, 0, 0), (1, 1, 0), None, (0, -1, 0) ) )
                if v == 0:
                    mvPointsx [v][n3] += reversed(xx [n3])
                    mvPointsd1[v][n3] += reversed(xd1[n3])
                    mvPointsd2[v][n3] += reversed(xd2[n3])
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += reversed(xd3[n3])
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append( ( (0, -1, 0), (1, 0, 0), None ) )
                else:
                    mvPointsx [v][n3] += xx [n3]
                    mvPointsd1[v][n3] += xd1[n3]
                    mvPointsd2[v][n3] += xd2[n3]
                    if useCubicHermiteThroughOstiumWall:
                        mvPointsd3[v][n3] += xd3[n3]
                    for i in range(elementsCountAcross - 1):
                        mvDerivativesMap[v][n3].append( ( (0, 1, 0), (-1, 0, 0), None ) )

    # calculate derivative 2 around free sides of inlets to fit vessel derivatives
    for v in range(vesselsCount):
        for n3 in range(2):
            for n1 in range(elementsCountAroundVessel):
                d2Map = mvDerivativesMap[v][n3][n1][1] if (mvDerivativesMap[v][n3][n1]) else None
                sf1 = d2Map[0] if d2Map else 0.0
                sf2 = d2Map[1] if d2Map else 1.0
                nx = [ vox[v][n3][n1], mvPointsx[v][n3][n1] ]
                nd2 = [ [ d*elementsCountAlong for d in vod2[v][n3][n1] ], [ (sf1*mvPointsd1[v][n3][n1][c] + sf2*mvPointsd2[v][n3][n1][c]) for c in range(3) ] ]
                nd2f = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDirection = True)[1]
                ndf = [ d/elementsCountAlong for d in nd2f ]
                # assign components to set original values:
                if sf1 == 0:
                    for c in range(3):
                        mvPointsd2[v][n3][n1][c] = sf2*ndf[c]
                elif sf2 == 0:
                    for c in range(3):
                        mvPointsd1[v][n3][n1][c] = sf1*ndf[c]
                else:
                    #print('v', v, 'n3', n3, ':', vector.magnitude(ndf), 'vs.', vector.magnitude(nd2[1]), 'd2Map', d2Map)
                    pass

    ##############
    # Create nodes
    ##############

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

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

    nodeIdentifier = startNodeIdentifier

    oNodeId = []
    for n3 in range(2):
        oNodeId.append([])
        for n1 in range(elementsCountAroundOstium):
            node = nodes.createNode(nodeIdentifier, nodetemplate if useCubicHermiteThroughOstiumWall else nodetemplateLinearS3)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ox [n3][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, od1[n3][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, od2[n3][n1])
            if useCubicHermiteThroughOstiumWall:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, od3[n3][n1])
            oNodeId[n3].append(nodeIdentifier)
            nodeIdentifier += 1

    if vesselsCount > 1:
        xNodeId = []
        for n3 in range(2):
            xNodeId.append([])
            for n2 in range(elementsCountAcross - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate if useCubicHermiteThroughOstiumWall else nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xx [n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xd1[n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xd2[n3][n2])
                if useCubicHermiteThroughOstiumWall:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xd3[n3][n2])
                xNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

    #for v in range(vesselsCount):
    #    node = nodes.createNode(nodeIdentifier, nodetemplate)
    #    cache.setNode(node)
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vcx [v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vcd1[v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vcd2[v])
    #    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vcd3[v])
    #    nodeIdentifier += 1
    #    for n3 in range(2):
    #        for n1 in range(elementsCountAroundVessel):
    #            node = nodes.createNode(nodeIdentifier, nodetemplate if useCubicHermiteThroughVesselWall else nodetemplateLinearS3)
    #            cache.setNode(node)
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vox [v][n3][n1])
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vod1[v][n3][n1])
    #            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vod2[v][n3][n1])
    #            if useCubicHermiteThroughVesselWall:
    #                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vod3[v][n3][n1])
    #            #vNodeId.append(nodeIdentifier)
    #            nodeIdentifier += 1

    # get identifiers of nodes around each vessel at ostium end
    mvNodeId = [ None ]*vesselsCount
    for v in range(vesselsCount):
        if vesselsCount == 1:
            mvNodeId[v] = oNodeId
        else:
            os = 0 if (v == 0) else (nodesCountFree - 1)
            of = os + nodesCountFree
            mvNodeId[v] = []
            for n3 in range(2):
                mvNodeId[v].append((oNodeId[n3]*2)[os:of])
                mvNodeId[v][n3] += reversed(xNodeId[n3]) if (v == 0) else xNodeId[n3]

    #################
    # Create elements
    #################

    mesh = fm.findMeshByDimension(3)
    elementIdentifier = startElementIdentifier

    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    #tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

    #eft = tricubichermite.createEftBasic()
    #elementtemplate = mesh.createElementtemplate()
    #elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    #elementtemplate.defineField(coordinates, -1, eft)

    #elementtemplateX = mesh.createElementtemplate()
    #elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    for v in range(vesselsCount):
        startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap = \
            vox[v], vod1[v], vod2[v], vod3[v] if useCubicHermiteThroughVesselWall else None, None, None
        endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap = \
            mvPointsx[v], mvPointsd1[v], mvPointsd2[v], mvPointsd3[v], mvNodeId[v], mvDerivativesMap[v]
        nodeIdentifier, elementIdentifier = tricubichermite.createAnnulusMesh3d(
            startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
            endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
            nodetemplate, nodetemplateLinearS3, nodeIdentifier, elementIdentifier,
            elementsCountRadial = elementsCountAlong,
            meshGroups = vesselMeshGroups[v] if vesselMeshGroups else [])

    oNodeId = []

    fm.endChange()
    return nodeIdentifier, elementIdentifier, (ox, od1, od2, od3, oNodeId)
