"""
Generates 3-D mesh of left and right ventricles below base plane.
Variant using collapsed/wedge elements at septum junction.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import createEllipsoidPoints, getApproximateEllipsePerimeter, getEllipseArcLength, getEllipseRadiansToX
from scaffoldmaker.utils.interpolation import computeCubicHermiteDerivativeScaling, getCubicHermiteArcLength, interpolateSampleCubicHermite, \
    sampleCubicHermiteCurves, sampleCubicHermiteCurvesSmooth, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.shieldmesh import ShieldMesh
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes


class MeshType_3d_heartventricles3(Scaffold_base):
    '''
    Generates 3-D mesh of left and right ventricles below base plane.
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles 3'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Unit Human 1',
            'Unit Mouse 1',
            'Unit Pig 1',
            'Unit Rat 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        isHuman = 'Human' in parameterSetName
        isMouse = 'Mouse' in parameterSetName
        isPig = 'Pig' in parameterSetName
        isRat = 'Rat' in parameterSetName
        unitScale = 'Unit' in parameterSetName
        options = {}
        # first 3 are all even or all odd
        options['Number of elements across septum'] = 6  # S
        options['Number of elements around LV free wall'] = 8  # L
        options['Number of elements around RV free wall'] = 8  # R
        #options['Number of elements through LV wall'] = 1
        options['Number of elements up LV free wall'] = 5  # U
        options['Number of elements up LV apex'] = 0  # A each contributes 2 to number around LV free wall
        # N = L - 2A = 8
        # M = 2U - N = 6 (require 
        options['Unit scale'] = 1.0
        options['Interventricular septum thickness'] = 0.12
        options['Interventricular sulcus apex transition length'] = 0.2
        options['Interventricular sulcus base transition length'] = 0.15
        options['LV apex thickness'] = 0.08
        options['LV outer height axis length'] = 0.8
        options['LV outer height'] = 1.0
        options['LV outer diameter'] = 1.0
        options['LV free wall thickness'] = 0.15
        options['RV apex cusp angle degrees'] = 10.0
        options['RV apex length factor'] = 0.35
        options['RV proportion up LV'] = 0.85
        options['RV free wall thickness'] = 0.05
        options['RV inlet angle degrees'] = 80.0
        options['RV inlet cusp angle degrees'] = 80.0
        options['RV inlet position around LV'] = 0.35
        options['RV outlet angle degrees'] = 45.0
        options['RV outlet cusp angle degrees'] = 70.0
        options['RV outlet position around LV'] = 0.7
        options['RV width'] = 0.3
        options['Use cross derivatives'] = False
        options['Refine'] = False
        options['Refine number of elements surface'] = 4
        options['Refine number of elements through wall'] = 1

        if isHuman:
            if not unitScale:
                options['Unit scale'] = 80.0
        elif isMouse or isRat:
            if not unitScale:
                options['Unit scale'] = 5.0 if isMouse else 12.0
        elif isPig:
            options['Number of elements up LV apex'] = 1
            if not unitScale:
                options['Unit scale'] = 80.0
            options['LV apex thickness'] = 0.07
            options['LV free wall thickness'] = 0.17
            options['Interventricular septum thickness'] = 0.14
            options['RV apex cusp angle degrees'] = 45.0
            options['RV apex length factor'] = 0.7
            options['RV proportion up LV'] = 0.6
            options['RV free wall thickness'] = 0.06
            options['RV outlet angle degrees'] = 60.0
            options['RV width'] = 0.35
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements across septum',
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            #'Number of elements through LV wall',
            'Number of elements up LV free wall',
            'Number of elements up LV apex',
            'Unit scale',
            'Interventricular septum thickness',
            'Interventricular sulcus apex transition length',
            'Interventricular sulcus base transition length',
            'LV apex thickness',
            'LV outer height axis length',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
            'RV apex cusp angle degrees',
            'RV apex length factor',
            'RV proportion up LV',
            'RV free wall thickness',
            'RV inlet angle degrees',
            'RV inlet cusp angle degrees',
            'RV inlet position around LV',
            'RV outlet angle degrees',
            'RV outlet cusp angle degrees',
            'RV outlet position around LV',
            'RV width',
            #'Use cross derivatives',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall'
        ]

    @classmethod
    def getElementsCounts(cls, options):
        elementsCountAcrossIVSeptum = options['Number of elements across septum']
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        #elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountUpLVFreeWall = options['Number of elements up LV free wall']
        elementsCountUpLVApex = options['Number of elements up LV apex']
        nas = elementsCountAcrossIVSeptum//2
        nal = elementsCountAroundLVFreeWall//2 - elementsCountUpLVApex
        nar = elementsCountAroundRVFreeWall//2
        nul = elementsCountUpLVFreeWall
        # number around half = must be matched on RV and septum
        nah = nal + elementsCountUpLVFreeWall - 2
        elementsCountAroundFull = 2*nah + (options['Number of elements across septum'] % 2)
        # number up rv, min 2
        elementsCountUpRVFreeWall = nah + 2 - nar
        elementsCountUpIVSeptum = nah + 2 - elementsCountAcrossIVSeptum//2
        return elementsCountAroundFull, elementsCountUpLVApex, \
            elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
            elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
            elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall

    @classmethod
    def checkOptions(cls, options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        if options['Number of elements up LV apex'] < 0:
            options['Number of elements up LV apex'] = 0
        if options['Number of elements up LV free wall'] < (options['Number of elements up LV apex'] + 2):
            options['Number of elements up LV free wall'] = options['Number of elements up LV apex'] + 2
            dependentChanges = True
        for key in [
            'Number of elements across septum',
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            ]:
            if options[key] < 4:
                options[key] = 4
        if options['Number of elements across septum'] % 1:
            if 0 == (options['Number of elements around LV free wall'] % 1):
                options['Number of elements around LV free wall'] += 1
                dependentChanges = True
            if 0 == (options['Number of elements around RV free wall'] % 1):
                options['Number of elements around RV free wall'] += 1
                dependentChanges = True
        nas = options['Number of elements across septum']//2
        nal = options['Number of elements around LV free wall']//2 - options['Number of elements up LV apex']
        nar = options['Number of elements around RV free wall']//2
        nul = options['Number of elements up LV free wall']
        # number around half = must be matched on RV and septum
        nah = nal + nul - 2
        # number up rv, min 2
        nur = nah + 2 - nar
        if nur < 2:
            options['Number of elements around RV free wall'] += 2*(2 - nur)
            dependentChanges = True
        # number up septum, min 2
        nus = nah + 2 - nas
        if nus < 2:
            options['Number of elements across septum'] += 2*(2 - nus)
            dependentChanges = True
        for key in [
            #'Number of elements through LV wall',
            'Refine number of elements surface',
            'Refine number of elements through wall',
            ]:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Unit scale',
            'Interventricular septum thickness',
            'Interventricular sulcus apex transition length',
            'Interventricular sulcus base transition length',
            'LV apex thickness',
            'LV outer height axis length',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
            'RV apex cusp angle degrees',
            'RV free wall thickness',
            'RV inlet cusp angle degrees',
            'RV outlet cusp angle degrees',
            'RV width'
            ]:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'RV apex length factor',
            'RV inlet position around LV',
            'RV outlet position around LV',
            'RV proportion up LV'
            ]:
            if options[key] < 0.01:
                options[key] = 0.01
            elif options[key] > 0.99:
                options[key] = 0.99
        return dependentChanges

    @staticmethod
    def getRVEdgePoints(lvTrackSurface, position, cuspDirection, cuspAngleRadians, rvFreeWallThickness):
        """
        Given position on LV outer track surface, cusp direction, cusp angle, wall thickness, get position on outer wall thickness.
        :return: x, d1, d2, d3
        """
        x, sd1, sd2 = lvTrackSurface.evaluateCoordinates(position, derivatives = True)
        td1, td2, td3 = calculate_surface_axes(sd1, sd2, vector.normalise(cuspDirection))
        cosCuspAngle = math.cos(cuspAngleRadians)
        sinCuspAngle = math.sin(cuspAngleRadians)
        pd1 = [ (cosCuspAngle*td1[c] - sinCuspAngle*td3[c]) for c in range(3) ]
        pd2 = td2
        pd3 = [ rvFreeWallThickness*(cosCuspAngle*td3[c] + sinCuspAngle*td1[c]) for c in range(3) ]
        px = [ (x[c] + pd3[c]) for c in range(3) ]
        return [x, px], [ pd1, pd1 ], [ pd2, pd2 ], [ pd3, pd3 ]

    @staticmethod
    def getTransitionRVtoLV(lvTrackSurface, elementsCount, cuspProportions, cuspDerivative, endProportions, ivSulcusTransitionLength, rvx, rvd1):
        '''
        Transition from RV to a position on lvTrackSurface.
        :return: tx, td1, td2, td3, tProportions
        '''
        tx, td1, td2, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(
            cuspProportions[0], cuspProportions[1], endProportions[0], endProportions[1], elementsCount, derivativeStart=cuspDerivative)
        if elementsCount > 1:
            uniformLength = -ivSulcusTransitionLength
            for n in range(elementsCount):
                uniformLength += getCubicHermiteArcLength(tx[n], td1[n], tx[n + 1], td1[n + 1])
            uniformLength /= (elementsCount - 1)
            tx, td1 = sampleCubicHermiteCurves(tx, td1, elementsCount, addLengthStart=(ivSulcusTransitionLength - uniformLength))[0:2]
        tx [0] = rvx
        td1[0] = rvd1
        tx, td1 = sampleCubicHermiteCurves(tx, td1, elementsCount, addLengthStart=0.5*vector.magnitude(td1[0]), lengthFractionStart=0.5, arcLengthDerivatives=True)[0:2]
        p1 = lvTrackSurface.createPositionProportion(*cuspProportions)
        for i in range(0, elementsCount + 1):
            p1 = lvTrackSurface.findNearestPosition(tx[i], p1)
            x, sd1, sd2 = lvTrackSurface.evaluateCoordinates(p1, derivatives = True)
            d1, d2, d3 = calculate_surface_axes(sd1, sd2, vector.normalise(td1[i]))
            magd2 = vector.magnitude(d2)
            if magd2 > 0.0:
                scale = vector.magnitude(td1[i])/magd2
                td2[i] = [ d*scale for d in d2 ]
            else:
                td2[i] = [ 0.0, 0.0, 0.0 ]
            td3[i] = d3
            tProportions[i] = lvTrackSurface.getProportion(p1)
        return tx, td1, td2, td3, tProportions

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundFull, elementsCountUpLVApex, \
            elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
            elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
            elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall = cls.getElementsCounts(options)
        unitScale = options['Unit scale']
        ivSeptumThickness = unitScale*options['Interventricular septum thickness']
        ivSulcusApexTransitionLength = options['Interventricular sulcus apex transition length']
        ivSulcusBaseTransitionLength = options['Interventricular sulcus base transition length']
        lvApexThickness = unitScale*options['LV apex thickness']
        lvOuterHeightAxisLength = unitScale*options['LV outer height axis length']
        lvOuterHeight = unitScale*options['LV outer height']
        lvOuterRadius = unitScale*0.5*options['LV outer diameter']
        lvFreeWallThickness = unitScale*options['LV free wall thickness']
        rvApexCuspAngleRadians = math.radians(options['RV apex cusp angle degrees'])
        rvApexLengthFactor = options['RV apex length factor']
        rvProportionUpLV = options['RV proportion up LV']
        rvFreeWallThickness = unitScale*options['RV free wall thickness']
        rvInletAngleRadians = math.radians(options['RV inlet angle degrees'])
        rvInletCuspAngleRadians = math.radians(options['RV inlet cusp angle degrees'])
        rvInletPositionAroundLV = options['RV inlet position around LV']
        rvOutletAngleRadians = math.radians(options['RV outlet angle degrees'])
        rvOutletCuspAngleRadians = math.radians(options['RV outlet cusp angle degrees'])
        rvOutletPositionAroundLV = options['RV outlet position around LV']
        rvWidth = unitScale*options['RV width']
        useCrossDerivatives = options['Use cross derivatives']

        #print("elementsCountAroundFull", elementsCountAroundFull)
        #print("elementsCountUpRVFreeWall", elementsCountUpRVFreeWall)

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        cache = fieldmodule.createFieldcache()

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldmodule.findMeshByDimension(3)

        lvGroup = AnnotationGroup(region, "left ventricle myocardium", FMANumber = 9558, lyphID = 'Lyph ID unknown')
        rvGroup = AnnotationGroup(region, "right ventricle myocardium", FMANumber = 9535, lyphID = 'Lyph ID unknown')
        vSeptumGroup = AnnotationGroup(region, "interventricular septum", FMANumber = 7133, lyphID = 'Lyph ID unknown')
        annotationGroups = [ lvGroup, rvGroup, vSeptumGroup ]

        centre = [ 0.0, 0.0, 0.0 ]
        poleAxis = vector.setMagnitude([ 0.0, 0.0, -1.0 ], lvOuterHeightAxisLength)
        sideAxis = vector.setMagnitude([ -1.0, 0.0, 0.0 ], lvOuterRadius)

        elementsCountAroundLVTrackSurface = 16
        elementsCountUpLVTrackSurface = 8

        useHeight = min(max(0.0, lvOuterHeight), 2.0*lvOuterHeightAxisLength)
        baseRadiansUp = getEllipseRadiansToX(lvOuterHeightAxisLength, 0.0, lvOuterHeightAxisLength - useHeight, initialTheta = 0.5*math.pi*useHeight/lvOuterHeightAxisLength)
        baseProportionUp = 2.0*getEllipseArcLength(lvOuterHeightAxisLength, lvOuterRadius, 0.0, baseRadiansUp)/getApproximateEllipsePerimeter(lvOuterHeightAxisLength, lvOuterRadius)

        nx, nd1, nd2 = createEllipsoidPoints(centre, poleAxis, sideAxis, elementsCountAroundLVTrackSurface, elementsCountUpLVTrackSurface, 2.0*lvOuterHeightAxisLength)

        lvTrackSurface = TrackSurface(elementsCountAroundLVTrackSurface, elementsCountUpLVTrackSurface, nx, nd1, nd2, loop1 = True)

        rvApexProportion1 = 0.5
        rvApexProportion2 = baseProportionUp*(1.0 - rvProportionUpLV)
        rvApexPosition = lvTrackSurface.createPositionProportion(rvApexProportion1, rvApexProportion2)
        rvApexCuspCoordinates, sd1, sd2 = lvTrackSurface.evaluateCoordinates(rvApexPosition, derivatives = True)
        rvApexLengthDerivative = vector.setMagnitude(sd1, rvApexLengthFactor)
        rvApexCuspDirection = vector.setMagnitude(sd2, -1.0)
        rvax, rvad1, rvad2, rvad3 = cls.getRVEdgePoints(lvTrackSurface, rvApexPosition, rvApexCuspDirection, rvApexCuspAngleRadians, rvFreeWallThickness)

        elementLengthUpLV = vector.magnitude(lvTrackSurface.createHermiteCurvePoints(0.0, 0.0, 0.0, baseProportionUp, elementsCountUpLVFreeWall)[2][-1])

        # Get RV inlet, outlet
        # sample curves around RV cusp
        elementsCountAroundHalf = elementsCountAroundFull//2
        elementsCountUpRVRegular = elementsCountUpRVFreeWall - 2
        elementsCountAroundRVFreeWallHalf = elementsCountAroundRVFreeWall//2

        # get RV cusp curves from inlet -> apex -> outlet

        inletPosition = lvTrackSurface.createPositionProportion(rvInletPositionAroundLV, baseProportionUp)
        rvInletCuspCoordinates, sd1, sd2 = lvTrackSurface.evaluateCoordinates(inletPosition, derivatives = True)
        td1, td2, rvInletSurfaceNormal = calculate_surface_axes(sd1, sd2, vector.normalise(sd1))
        cosAngle = math.cos(rvInletAngleRadians)
        sinAngle = math.sin(rvInletAngleRadians)
        rvBaseLengthFactor = (1.0 - rvApexLengthFactor)
        rvInletDerivative = [ rvBaseLengthFactor*(cosAngle*td1[c] - sinAngle*td2[c]) for c in range(3) ]

        outletPosition = lvTrackSurface.createPositionProportion(rvOutletPositionAroundLV, baseProportionUp)
        rvOutletCuspCoordinates, sd1, sd2 = lvTrackSurface.evaluateCoordinates(outletPosition, derivatives = True)
        td1, td2, rvOutletSurfaceNormal = calculate_surface_axes(sd1, sd2, vector.normalise(sd1))
        cosAngle = math.cos(rvOutletAngleRadians)
        sinAngle = math.sin(rvOutletAngleRadians)
        rvOutletDerivative = [ rvBaseLengthFactor*(cosAngle*td1[c] + sinAngle*td2[c]) for c in range(3) ]

        scaling = computeCubicHermiteDerivativeScaling(rvApexCuspCoordinates, rvApexLengthDerivative, rvOutletCuspCoordinates, rvOutletDerivative)/elementsCountAroundHalf
        rvApexLengthDerivative = [ scaling*d for d in rvApexLengthDerivative ]
        rvInletDerivative = [ scaling*d for d in rvInletDerivative ]
        rvOutletDerivative = [ scaling*d for d in rvOutletDerivative ]

        rvox, rvod1, rvod2, rvod3, rvoProportions = lvTrackSurface.createHermiteCurvePoints(rvApexProportion1, rvApexProportion2, rvOutletPositionAroundLV, baseProportionUp,
            elementsCountAroundHalf, derivativeStart = rvApexLengthDerivative, derivativeEnd = rvOutletDerivative, curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
        rvox, rvod1, rvod2, rvod3, rvoProportions = lvTrackSurface.resampleHermiteCurvePointsSmooth(rvox, rvod1, rvod2, rvod3, rvoProportions,
            derivativeMagnitudeEnd=elementLengthUpLV)
        rvix, rvid1, rvid2, rvid3, rviProportions = lvTrackSurface.createHermiteCurvePoints(rvInletPositionAroundLV, baseProportionUp, rvApexProportion1, rvApexProportion2,
            elementsCountAroundHalf, derivativeStart = rvInletDerivative, derivativeEnd = rvApexLengthDerivative, curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
        rvix, rvid1, rvid2, rvid3, rviProportions = lvTrackSurface.resampleHermiteCurvePointsSmooth(rvix, rvid1, rvid2, rvid3, rviProportions,
            derivativeMagnitudeStart=elementLengthUpLV, derivativeMagnitudeEnd=vector.magnitude(rvod1[0]))
        rvix  += rvox [1:]
        rvid1 += rvod1[1:]
        rvid2 += rvod2[1:]
        rvid3 += rvod3[1:]
        rviProportions += rvoProportions[1:]

        #print('elementsCountAroundHalf',elementsCountAroundHalf)
        #print('elementsCountUpRVRegular',elementsCountUpRVRegular)
        #print('elementsCountAroundRVFreeWallHalf',elementsCountAroundRVFreeWallHalf)
        rvox  = []
        rvod1 = []
        rvod2 = []
        rvod3 = []
        for n in range(elementsCountAroundFull + 1):
            proportion1, proportion2 = rviProportions[n]
            cuspDirection = vector.setMagnitude(rvid2[n], -1.0)
            nRadians = n*math.pi/elementsCountAroundFull
            cuspAngleRadians = math.sin(nRadians)*math.sin(nRadians)*rvApexCuspAngleRadians + \
                               math.cos(nRadians)*math.cos(nRadians)*(rvInletCuspAngleRadians if (n <= elementsCountAroundHalf) else rvOutletCuspAngleRadians)
            rx, rd1, rd2, rd3 =  cls.getRVEdgePoints(lvTrackSurface, lvTrackSurface.createPositionProportion(proportion1, proportion2), 
                                                     cuspDirection, cuspAngleRadians, rvFreeWallThickness)
            rvix [n] = rx [0]
            rvid1[n] = rd1[0]
            rvid2[n] = rd2[0]
            rvid3[n] = rd3[0]
            rvox .append(rx [1])
            rvod1.append(rd1[1])
            rvod2.append(rd2[1])
            rvod3.append(rd3[1])
        length = 0.0
        for i in range(1, len(rvod2)):
            for c in range(3):
                length += (rvox[i][c] - rvox[i - 1][c])*(rvox[i][c] - rvox[i - 1][c])
        length /= (len(rvod2) - 1)
        rvid2 = smoothCubicHermiteDerivativesLine(rvox, [ vector.setMagnitude(rvid2[i], length) for i in range(len(rvid2)) ], fixAllDirections = True)
        rvod2 = smoothCubicHermiteDerivativesLine(rvox, [ vector.setMagnitude(rvod2[i], length) for i in range(len(rvod2)) ], fixAllDirections = True)

        # get position of centre of RV freewall
        rvInletLateralDirection = vector.normalise(vector.crossproduct3(rvInletSurfaceNormal, rvInletDerivative))
        rvOutletLateralDirection = vector.normalise(vector.crossproduct3(rvOutletDerivative, rvOutletSurfaceNormal))
        td1 = smoothCubicHermiteDerivativesLine([ rvInletCuspCoordinates, rvOutletCuspCoordinates ], [ rvInletLateralDirection, rvOutletLateralDirection ], fixAllDirections = True)
        px, pd1, pd2, pd3, pProportions = lvTrackSurface.createHermiteCurvePoints(rvInletPositionAroundLV, baseProportionUp, rvOutletPositionAroundLV, baseProportionUp, elementsCount = 2,
            derivativeStart = [ 0.5*d for d in td1[0] ], derivativeEnd = [ 0.5*d for d in td1[1] ], curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
        # get regular spacing up centre of RV
        rbcx = [ (px[1][c] + (rvWidth + rvFreeWallThickness)*pd3[1][c]) for c in range(3) ]
        rcx = [ rvax[1], rbcx ]
        rcd1 = [ vector.setMagnitude(rvad1[1], -1.0), pd2[1] ]
        #print("rcx ", rcx)
        #print("rcd1", rcd1)
        rcd1 = smoothCubicHermiteDerivativesLine(rcx, rcd1, fixAllDirections = True)
        #rscx, rscd1, rscpe, rscpxi, rscpsf = sampleCubicHermiteCurves(rcx, rcd1, elementsCountUpRVFreeWall, arcLengthDerivatives = True)
        rscx, rscd2 = sampleCubicHermiteCurves(rcx, rcd1, elementsCountUpRVFreeWall, lengthFractionStart=1.0, arcLengthDerivatives = True)[0:2]  # GRC fudge factor
        # get d1, d3
        rscd1 = []
        rscd3 = []
        for n in range(len(rscx)):
            a = [ (rscx[n][c] - centre[c]) for c in range(3) ]
            d1 = vector.normalise(vector.crossproduct3(rscd2[n], a))
            d3 = vector.normalise(vector.crossproduct3(d1, rscd2[n]))
            rscd1.append(d1)
            rscd3.append(d3)

        # transition from RV apex to LV apex:
        lax, lad2, lad1, lad3, laProportions = cls.getTransitionRVtoLV(lvTrackSurface, elementsCountUpLVApex + 1,
            [ rvApexProportion1, rvApexProportion2 ], vector.setMagnitude(rvApexCuspDirection, ivSulcusApexTransitionLength),
            [ 0.5, 0.0 ],
            ivSulcusApexTransitionLength, rvx=rscx[0], rvd1=[ -s for s in rscd2[0] ])
        # fix apex
        lad3[-1] = vector.normalise(poleAxis)
        lad1[-1] = vector.crossproduct3(lad3[-1], lad2[-1])
        # reverse d1
        lad1 = [ [-s for s in d ] for d in lad1 ]
        lad3 = [ vector.setMagnitude(d, lvFreeWallThickness) for d in lad3 ]

        rvShield = ShieldMesh(elementsCountUpRVFreeWall, elementsCountAroundRVFreeWall, 0)
        rvx  = rvShield.px
        rvd1 = rvShield.pd1
        rvd2 = rvShield.pd2
        rvd3 = rvShield.pd3
        for n in range(elementsCountAroundFull + 1):
            if n <= elementsCountUpRVRegular:
                n1 = 0
                n2 = elementsCountUpRVFreeWall - n
            elif n < (elementsCountAroundFull - elementsCountUpRVRegular):
                n1 = n - elementsCountAroundHalf + elementsCountAroundRVFreeWallHalf
                n2 = 0
            else:
                n1 = -1
                n2  = n - elementsCountAroundFull + elementsCountUpRVFreeWall
            #print('n', n, 'n1', n1, 'n2', n2)
            rvx [0][n2][n1] = rvix[n]
            rvx [1][n2][n1] = rvox[n]
            mag = -1 if (n <= elementsCountAroundHalf) else 1.0
            rvd1[0][n2][n1] = vector.setMagnitude(rvid1[n], mag)
            rvd1[1][n2][n1] = vector.setMagnitude(rvod1[n], mag)
            if n < elementsCountAroundHalf:
                rvd2[0][n2][n1] = [ -rvid2[n][c] for c in range(3) ]
                rvd3[0][n2][n1] = [ -rvid3[n][c] for c in range(3) ]
                rvd2[1][n2][n1] = [ -rvod2[n][c] for c in range(3) ]
                rvd3[1][n2][n1] = [ -rvod3[n][c] for c in range(3) ]
            else:
                rvd2[0][n2][n1] = rvid2[n]
                rvd3[0][n2][n1] = rvid3[n]
                rvd2[1][n2][n1] = rvod2[n]
                rvd3[1][n2][n1] = rvod3[n]

        # around regular rows of RV
        for n2 in range(2, elementsCountUpRVFreeWall + 1):
            rvx[1][n2], rvd1[1][n2], pe, pxi, psf = sampleCubicHermiteCurves(
                [ rvx[1][n2][0], rscx[n2], rvx[1][n2][-1] ], [ rvd1[1][n2][0], rscd1[n2], rvd1[1][n2][-1] ], elementsCountAroundRVFreeWall,
                lengthFractionStart=1.0, lengthFractionEnd=1.0, arcLengthDerivatives = True)  # GRC fudge factor
            rvd2[1][n2] = interpolateSampleCubicHermite([ rvd2[1][n2][0], rscd2[n2], rvd2[1][n2][-1] ], [ [ 0.0, 0.0, 0.0 ] ]*3, pe, pxi, psf)[0]

        # up regular columns of RV
        for n1 in range(2, elementsCountAroundRVFreeWall - 1):
            left = n1 < elementsCountAroundRVFreeWallHalf
            right = n1 > (elementsCountAroundRVFreeWall - elementsCountAroundRVFreeWallHalf)
            startd2 = rvd1[1][0][n1] if (not right) else [ -d for d in rvd1[1][0][n1] ]
            tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(
                [ rvx[1][0][n1], rvx[1][2][n1] ], [ startd2, rvd2[1][2][n1] ], 2, lengthFractionStart=1.0, arcLengthDerivatives = True)  # GRC fudge factor
            for n2 in range(3, elementsCountUpRVFreeWall + 1):
                tx .append(rvx [1][n2][n1])
                td2.append(rvd2[1][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection = True)
            startd1 = [ -d for d in rvd2[1][0][n1] ] if left else rvd2[1][0][n1]
            td1 = interpolateSampleCubicHermite([ startd1, rvd1[1][2][n1] ], [ [ 0.0, 0.0, 0.0 ] ]*2, pe, pxi, psf)[0]
            rvx [1][0][n1] = tx [0]
            if left:
                rvd1[1][0][n1] = td2[0]
                rvd2[1][0][n1] = [ -d for d in td1[0] ]
            elif right:
                rvd1[1][0][n1] = [ -d for d in td2[0] ]
                rvd2[1][0][n1] = td1[0]
            else:  # middle
                rvd1[1][0][n1] = td1[0]
                rvd2[1][0][n1] = td2[0]
            rvx [1][1][n1] = tx [1]
            rvd1[1][1][n1] = td1[1]
            for n2 in range(1, elementsCountUpRVFreeWall + 1):
                rvd2[1][n2][n1] = td2[n2]

        rvShield.getTriplePoints(n3=1)
        m1a = elementsCountAroundRVFreeWall
        m1b = m1a - 1
        m1c = m1a - 2

        # smooth RV freewall row 1
        n1b = 1
        n2b = 1
        rvd1[1][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(rvx[1][n2b][n1b:m1a], rvd1[1][n2b][n1b:m1a])

        # smooth RV columns 1, -2
        for n1 in [ 1, -2 ]:
            tx = []
            td2 = []
            for n2 in range(1, elementsCountUpRVFreeWall + 1):
                tx .append(rvx [1][n2][n1])
                td2.append(rvd2[1][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2)
            for n in range(elementsCountUpRVFreeWall):
                rvd2[1][n + 1][n1] = td2[n]

        # fix outer derivatives leading to triple points
        rvd1[1][0][1] = smoothCubicHermiteDerivativesLine([ rvx[1][0][1], rvx[1][1][1] ], [ rvd1[1][0][1], [ (rvd1[1][1][1][c] + rvd2[1][1][1][c]) for c in range(3) ] ],
                                                          fixEndDerivative = True, fixStartDirection = True)[0]
        rvd1[1][0][m1b] = smoothCubicHermiteDerivativesLine([ rvx[1][1][m1b], rvx[1][0][m1b] ], [ [ (rvd1[1][1][m1b][c] - rvd2[1][1][m1b][c]) for c in range(3) ], rvd1[1][0][m1b] ],
                                                          fixStartDerivative = True, fixEndDirection = True)[1]

        # get outer d3 and inner x, d3
        for n2 in range(elementsCountUpRVFreeWall + 1):
            for n1 in range(elementsCountAroundRVFreeWall + 1):
                if rvd1[1][n2][n1]:
                    rvd3[0][n2][n1] = rvd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(rvd1[1][n2][n1], rvd2[1][n2][n1]), rvFreeWallThickness)
                    rvx [0][n2][n1] = [ (rvx [1][n2][n1][c] - rvd3[1][n2][n1][c]) for c in range(3) ]

        # get inner d1, d2
        # row 1
        rvd1[0][1][1:-1] = smoothCubicHermiteDerivativesLine(rvx[0][1][1:-1], rvd1[1][1][1:-1], fixAllDirections = True)
        # regular rows 2+
        for n2 in range(2, elementsCountUpRVFreeWall + 1):
            rvd1[0][n2] = smoothCubicHermiteDerivativesLine(rvx[0][n2], rvd1[1][n2], fixAllDirections = True)
        # column 1 and -2
        # regular columns
        for n1 in range(1, elementsCountAroundRVFreeWall):
            startn2 = 1 if (n1 in [1, elementsCountAroundRVFreeWall - 1]) else 0
            tx  = []
            td2 = []
            if startn2 == 0:
                tx .append(rvx [0][0][n1])
                left = n1 < elementsCountAroundRVFreeWallHalf
                right = n1 > (elementsCountAroundRVFreeWall - elementsCountAroundRVFreeWallHalf)
                if left:
                    td2.append(rvd1[1][0][n1])
                elif right:
                    td2.append([ -d for d in rvd1[1][0][n1] ])
                else:  # middle:
                    td2.append(rvd2[1][0][n1])
            for n in range(elementsCountUpRVFreeWall):
                tx .append(rvx [0][n + 1][n1])
                td2.append(rvd2[1][n + 1][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixAllDirections = True)
            if startn2 == 0:
                if left:
                    rvd1[0][0][n1] = td2[0]
                elif right:
                    rvd1[0][0][n1] = [ -d for d in td2[0] ]
                else:  # middle
                    rvd1[0][0][n1] = rvd2[0][0][n1]
                    rvd2[0][0][n1] = td2[0]
            for n2 in range(1, elementsCountUpRVFreeWall + 1):
                rvd2[0][n2][n1] = td2[n2 - startn2]

        # fix inner derivatives leading to triple points
        rvd1[0][0][1] = smoothCubicHermiteDerivativesLine([ rvx[0][0][1], rvx[0][1][1] ], [ rvd1[0][0][1], [ (rvd1[0][1][1][c] + rvd2[0][1][1][c]) for c in range(3) ] ],
                                                          fixEndDerivative = True, fixStartDirection = True)[0]
        rvd1[0][0][m1b] = smoothCubicHermiteDerivativesLine([ rvx[0][1][m1b], rvx[0][0][m1b] ], [ [ (rvd1[0][1][m1b][c] - rvd2[0][1][m1b][c]) for c in range(3) ], rvd1[0][0][m1b] ],
                                                          fixStartDerivative = True, fixEndDirection = True)[1]

        # LV free wall
        elementsCountUpLV = elementsCountUpLVFreeWall + elementsCountUpLVApex
        lvShield = ShieldMesh(elementsCountUpLV, elementsCountAroundLVFreeWall, elementsCountUpLVApex, lvTrackSurface)
        lvx  = lvShield.px
        lvd1 = lvShield.pd1
        lvd2 = lvShield.pd2
        lvd3 = lvShield.pd3
        lvProportions = lvShield.pProportions

        # sample around top of LV free wall to get centre
        #td1 = smoothCubicHermiteDerivativesLine([ rvOutletCuspCoordinates, rvInletCuspCoordinates ], [ rvOutletLateralDirection, rvInletLateralDirection ], fixAllDirections = True)
        # GRC fudge factor
        sulcusDerivativeFactor = 2.0
        td1 = [ vector.setMagnitude(d, sulcusDerivativeFactor*ivSulcusBaseTransitionLength) for d in [ rvOutletLateralDirection, rvInletLateralDirection ] ]
        px, pd1, pd2, pd3, pProportions = lvTrackSurface.createHermiteCurvePoints(rvOutletPositionAroundLV, baseProportionUp, rvInletPositionAroundLV + 1.0, baseProportionUp,
            elementsCount = 2, derivativeStart = td1[0], derivativeEnd = td1[1], curveMode = TrackSurface.HermiteCurveMode.SMOOTH)
        # sample up centre of LV free wall
        tx  = [ nx[0], px[1] ]
        td2 = [ nd2[0], pd2[1] ]
        td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixAllDirections = True)
        for i in range(2):
            td2[i] = [ d/elementsCountUpLVFreeWall for d in td2[i] ]
        #print('tx ', tx)
        #print('td2', td2)
        #print('pProportions', 1.0, 0.0, pProportions[1][0], pProportions[1][1])
        lscx, lscd2, lscd1, lscd3, lscProportions = lvTrackSurface.createHermiteCurvePoints(1.0, 0.0, pProportions[1][0], pProportions[1][1],
            elementsCount = elementsCountUpLVFreeWall, derivativeStart = lad2[-1], derivativeEnd = td2[1], curveMode = TrackSurface.HermiteCurveMode.TRANSITION_START)
        lscd3[0] = vector.normalise(vector.crossproduct3(nd2[0], nd2[elementsCountAroundLVTrackSurface*3//4]))
        lscd1 = [ vector.crossproduct3(lscd2[i], lscd3[i]) for i in range(len(lscd2)) ]

        # around regular rows of LV free wall, transition from RV
        for n2 in range(2 + elementsCountUpLVApex, elementsCountUpLV + 1):
            ix = elementsCountUpLV - n2
            ox = -1 - ix
            #print('n2',n2,'ix',ix,'ox',ox)
            # GRC fix:
            td1 = smoothCubicHermiteDerivativesLine([ lscx[ox], rvix[ix] ],
                [ lscd1[ox], vector.setMagnitude(rvid1[ix], -sulcusDerivativeFactor*ivSulcusBaseTransitionLength) ], fixStartDirection=True, fixEndDerivative=True)
            td1 = [ [ 0.25*s for s in d ] for d in td1 ]
            #print('pProp', lscProportions[ox][0], lscProportions[ox][1], rviProportions[ix][0] + 1.0, rviProportions[ix][1])
            px, pd1, pd2, pd3, pProportions = lvTrackSurface.createHermiteCurvePoints(
                lscProportions[ox][0], lscProportions[ox][1],
                rviProportions[ix][0] + 1.0, rviProportions[ix][1],
                elementsCount = 4, derivativeStart = td1[0], derivativeEnd = td1[1],
                curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
            #pex, ped1 = getCubicHermiteCurvesPointAtArcDistance(reversed(px), [ vector.setMagnitude(d, -1.0) for d in reversed(pd1) ], ivSulcusBaseTransitionLength)[0:2]
            #lvTrackSurface.findNearestPosition(x, lvTrackSurface.createPositionProportion(*rviProportions[ix]))
            #print('qProp', rviProportions[ox][0], rviProportions[ox][1], lscProportions[ox][0], lscProportions[ox][1])
            qx, qd1, qd2, qd3, qProportions = lvTrackSurface.createHermiteCurvePoints(
                rviProportions[ox][0], rviProportions[ox][1],
                lscProportions[ox][0], lscProportions[ox][1],
                elementsCount = 4, derivativeStart = vector.setMagnitude(rvid1[ox], 0.25*sulcusDerivativeFactor*ivSulcusBaseTransitionLength), derivativeEnd = td1[0],
                curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
            qx  += px [1:]
            qd1 += pd1[1:]
            qd2 += pd2[1:]
            qd3 += pd3[1:]
            qProportions += pProportions[1:]

            uniformLength = -2.0*ivSulcusBaseTransitionLength
            for n in range(8):
                uniformLength += getCubicHermiteArcLength(qx[n], qd1[n], qx[n + 1], qd1[n + 1])
            uniformLength /= 6
            addLength = ivSulcusBaseTransitionLength - uniformLength
            # T = 2S + 6U
            # U = (T - 2S)/6
            # A = S - U = (8S - T)/6
            rx, rd1 = sampleCubicHermiteCurves(qx, qd1, elementsCountOut=8, addLengthStart=addLength, addLengthEnd=addLength)[0:2]
            hx = 1 - elementsCountUpRVFreeWall - ox
            #print('hx', hx) GRC remove hx
            if hx < 1:
                rx [0] = rvx [1][ox][-1]
                rd1[0] = rvd1[1][ox][-1]
                rx [-1] = rvx[1][ox][0]
                rd1[-1] = rvd1[1][ox][0]
            else:
                rx [0] = rvx [1][0][-1 - hx]
                rd1[0] = rvd1[1][0][-1 - hx]
                rx [-1] = rvx[1][0][hx]
                rd1[-1] = rvd1[1][0][hx]
            tx, td1 = sampleCubicHermiteCurves(rx, rd1, elementsCountAroundLVFreeWall + 2,
                addLengthStart=0.5*vector.magnitude(rd1[0]), lengthFractionStart=0.5,
                addLengthEnd=0.5*vector.magnitude(rd1[-1]), lengthFractionEnd=0.5,
                arcLengthDerivatives=True)[0:2]
            td2 = [ None ]*(elementsCountAroundLVFreeWall + 3)
            td3 = [ None ]*(elementsCountAroundLVFreeWall + 3)
            tProportions = [ None ]*(elementsCountAroundLVFreeWall + 3)
            p1 = lvTrackSurface.createPositionProportion(*qProportions[1])
            for i in range(1, elementsCountAroundLVFreeWall + 2):
                p1 = lvTrackSurface.findNearestPosition(tx[i], p1)
                x, sd1, sd2 = lvTrackSurface.evaluateCoordinates(p1, derivatives = True)
                d1, d2, d3 = calculate_surface_axes(sd1, sd2, vector.normalise(td1[i]))
                td2[i] = vector.setMagnitude(d2, vector.magnitude(td1[i]))
                td3[i] = vector.setMagnitude(d3, lvFreeWallThickness)
                tProportions[i] = lvTrackSurface.getProportion(p1)
            lvx [1][ox] = tx [1:-1]
            lvd1[1][ox] = td1[1:-1]
            lvd2[1][ox] = td2[1:-1]
            lvd3[1][ox] = td3[1:-1]
            lvProportions[ox] = tProportions[1:-1]

        # smooth d2 up regular columns of LV
        n2reg = 2 + elementsCountUpLVApex
        for n1 in range(elementsCountAroundLVFreeWall + 1):
            tx = []
            td2 = []
            for n2 in range(2 + elementsCountUpLVApex, elementsCountUpLV + 1):
                tx .append(lvx [1][n2][n1])
                td2.append(lvd2[1][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2)  # GRC fix to make on lvTrackSurface, where possible?
            for n2 in range(n2reg, elementsCountUpLV + 1):
                lvd2[1][n2][n1] = td2[n2 - n2reg]

        # transition to LV apex, anterior and posterior
        elementsCountAroundLVFreeWallHalf = elementsCountAroundLVFreeWall//2
        #print('elementsCountAroundLVFreeWallHalf', elementsCountAroundLVFreeWallHalf)
        elementsCountRemaining = elementsCountAroundHalf - (elementsCountUpLVFreeWall - 2)
        for n2 in range(0, elementsCountUpLVApex + 1):
            n1reg = n2
            lan = n2 + 1  # GRC elementsCountUpLVApex - n2
            startProportion = laProportions[lan]
            if n2 == elementsCountUpLVApex:
                startProportion = [ 0.75, 0.0 ]
            #print('ant n2', n2, '/', elementsCountUpLVApex, 'start', startProportion, 'end', lvProportions[n2reg][n1reg], 'derivativeStart', lad1[lan])
            tx, td2, td1, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(startProportion[0], startProportion[1],
                lvProportions[n2reg][n1reg][0], lvProportions[n2reg][n1reg][1], elementsCountRemaining,
                derivativeStart=[ -d for d in lad1[lan] ],
                derivativeEnd=lvd2[1][n2reg][n1reg])
            tx, td2, td1, td3, tProportions = lvTrackSurface.resampleHermiteCurvePointsSmooth(tx, td2, td1, td3, tProportions,
                derivativeMagnitudeStart=vector.magnitude(lad1[lan]),
                derivativeMagnitudeEnd=vector.magnitude(lvd2[1][n2reg][n1reg]))
            for n in range(elementsCountRemaining):
                n1 = elementsCountAroundLVFreeWallHalf - n
                lvx [1][n2][n1] = tx [n]
                lvd1[1][n2][n1] = [ -d for d in td1[n] ] if n else [ -d for d in td2[n] ]
                lvd2[1][n2][n1] = td2[n] if n else lad2[lan]
                lvd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(lvd1[1][n2][n1], lvd2[1][n2][n1]), lvFreeWallThickness)
                #print('ant n2', n2, 'n1', n1)
                lvProportions[n2][n1] = tProportions[n]
                if False:
                    # transition from RV
                    if n:
                        ix = elementsCountAroundHalf + n
                        n1r, n2r = rvShield.convertRimIndex(ix)
                        #print('ant ix', ix, '-->', n1r, n2r)
                        elementsCount = 1 + elementsCountUpLVApex
                        qx, qd1, qd2, qd3, qProportions = cls.getTransitionRVtoLV(lvTrackSurface, elementsCount,
                            rviProportions[ix], vector.setMagnitude(rvid1[ix], ivSulcusBaseTransitionLength),  # GRC vary transition length from Apex to Base?
                            tProportions[n],
                            ivSulcusBaseTransitionLength, rvx[1][n2r][n1r], rvd1[1][n2r][n1r])
                        for i in range(1, elementsCount + 1):
                            n2l = i - 1
                            lvd1[1][n2l][n1] = qd1[i]
                            if i < elementsCount:
                                lvx [1][n2l][n1] = qx[i]
                                lvd2[1][n2l][n1] = vector.setMagnitude(qd2[i], vector.magnitude(qd1[i]))
                                lvd3[1][n2l][n1] = vector.setMagnitude(qd3[i], lvFreeWallThickness)
                                lvProportions[n2l][n1] = qProportions[i]
            n1reg = -1 - n2
            if n2 == elementsCountUpLVApex:
                startProportion = [ 1.25, 0.0 ]
            else:
                startProportion = [ laProportions[lan][0] + 1.0, laProportions[lan][1] ]
            #print('pos n2', n2, '/', elementsCountUpLVApex, 'start', startProportion, 'end', lvProportions[n2reg][n1reg])
            tx, td2, td1, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(startProportion[0], startProportion[1],
                lvProportions[n2reg][n1reg][0], lvProportions[n2reg][n1reg][1], elementsCountRemaining,
                #derivativeStart=lvd2[1][n2][elementsCountAroundLVFreeWallHalf],
                derivativeStart=lad1[lan],
                derivativeEnd=lvd2[1][n2reg][n1reg])
            tx, td2, td1, td3, tProportions = lvTrackSurface.resampleHermiteCurvePointsSmooth(tx, td2, td1, td3, tProportions,
                derivativeMagnitudeStart=vector.magnitude(lad1[lan]),
                derivativeMagnitudeEnd=vector.magnitude(lvd2[1][n2reg][n1reg]))
            for n in range(1, elementsCountRemaining):
                n1 = elementsCountAroundLVFreeWallHalf + n
                #print('pos n2', n2, 'n1', n1)
                lvx [1][n2][n1] = tx [n]
                lvd1[1][n2][n1] = [ -d for d in td1[n] ]
                lvd2[1][n2][n1] = td2[n]
                lvd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(lvd1[1][n2][n1], lvd2[1][n2][n1]), lvFreeWallThickness)
                lvProportions[n2][n1] = tProportions[n]
                if False:
                    # transition from RV
                    ix = elementsCountAroundHalf - n
                    n1r, n2r = rvShield.convertRimIndex(ix)
                    #print('pos ix', ix, '-->', n1r, n2r)
                    elementsCount = 1 + elementsCountUpLVApex
                    qx, qd1, qd2, qd3, qProportions = cls.getTransitionRVtoLV(lvTrackSurface, elementsCount,
                        [ rviProportions[ix][0] + 1.0, rviProportions[ix][1] ], vector.setMagnitude(rvid1[ix], -ivSulcusBaseTransitionLength),  # GRC vary transition length from Apex to Base?
                        tProportions[n],
                        ivSulcusBaseTransitionLength, rvx[1][n2r][n1r], [ -s for s in rvd1[1][n2r][n1r] ])
                    for i in range(1, elementsCount + 1):
                        n2l = i - 1
                        lvd1[1][n2l][n1] = [ -s for s in qd1[i] ]
                        if i < elementsCount:
                            lvx [1][n2l][n1] = qx[i]
                            lvd2[1][n2l][n1] = vector.setMagnitude(qd2[i], -vector.magnitude(qd1[i]))
                            lvd3[1][n2l][n1] = vector.setMagnitude(qd3[i], lvFreeWallThickness)
                            lvProportions[n2l][n1] = qProportions[i]

        # up regular columns of LV
        for n1 in range(2 + elementsCountUpLVApex, elementsCountAroundLVFreeWall - 1 - elementsCountUpLVApex):
            n2apex = elementsCountUpLVApex
            left = n1 < elementsCountAroundLVFreeWallHalf
            right = n1 > (elementsCountAroundLVFreeWall - elementsCountAroundLVFreeWallHalf)
            aProportions = lvProportions[n2apex][n1] if (left or right) else [ 1.0, 0.0 ]
            tx, td2, td1, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(
                aProportions[0], aProportions[1],
                lvProportions[n2apex + 2][n1][0], lvProportions[n2apex + 2][n1][1],
                elementsCount=2,
                derivativeStart=lvd1[1][n2apex][n1] if left else [ -d for d in lvd1[1][n2apex][n1] ] if right else lvd2[1][n2apex][n1],
                derivativeEnd = lvd2[1][n2apex + 2][n1])
            lvx [1][n2apex + 1][n1] = tx [1]
            lvd1[1][n2apex + 1][n1] = [-d for d in td1[1] ]
            lvd2[1][n2apex + 1][n1] = td2[1]
            lvd3[1][n2apex + 1][n1] = td3[1]
            lvProportions[n2apex + 1][n1] = tProportions[1]

        lvShield.getTriplePoints(n3=1)
        n1a = elementsCountUpLVApex
        n1b = n1a + 1
        n1c = n1a + 2
        m1a = elementsCountAroundLVFreeWall - elementsCountUpLVApex
        m1b = m1a - 1
        m1c = m1a - 2
        n2a = elementsCountUpLVApex
        n2b = n2a + 1
        n2c = n2a + 2

        # smooth LV freewall row 1
        lvd1[1][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(lvx[1][n2b][n1b:m1a], lvd1[1][n2b][n1b:m1a])

        # smooth LV columns 1, -2
        for n1 in [ n1b, -1 - n1b ]:
            tx = []
            td2 = []
            for n2 in range(n2b, elementsCountUpLV + 1):
                tx .append(lvx [1][n2][n1])
                td2.append(lvd2[1][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2)
            for n2 in range(n2b, elementsCountUpLV + 1):
                lvd2[1][n2][n1] = td2[n2 - n2b]

        # fix outer derivatives leading to triple points
        lvd1[1][n2a][n1b] = smoothCubicHermiteDerivativesLine([ lvx[1][n2a][n1b], lvx[1][n2b][n1b] ], [ lvd1[1][n2a][n1b], [ (lvd1[1][n2b][n1b][c] + lvd2[1][n2b][n1b][c]) for c in range(3) ] ],
                                                              fixEndDerivative = True, fixStartDirection = True)[0]
        lvd1[1][n2a][m1b] = smoothCubicHermiteDerivativesLine([ lvx[1][n2b][m1b], lvx[1][n2a][m1b] ], [ [ (lvd1[1][n2b][m1b][c] - lvd2[1][n2b][m1b][c]) for c in range(3) ], lvd1[1][n2a][m1b] ],
                                                              fixStartDerivative = True, fixEndDirection = True)[1]

        # get outer d3 and inner x, d3
        for n2 in range(elementsCountUpLV + 1):
            for n1 in range(elementsCountAroundLVFreeWall + 1):
                if lvd1[1][n2][n1]:
                    if lvd3[1][n2][n1]:
                        d3 = lvd3[1][n2][n1]
                    else:
                        d3 = vector.crossproduct3(lvd1[1][n2][n1], lvd2[1][n2][n1])
                    lvd3[0][n2][n1] = lvd3[1][n2][n1] = vector.setMagnitude(d3, lvFreeWallThickness)
                    lvx [0][n2][n1] = [ (lvx [1][n2][n1][c] - lvd3[1][n2][n1][c]) for c in range(3) ]

        # get inner d2 (d1 apex) around full rows up LV apex
        for r in [ elementsCountUpLVApex ]:  # grc for now
            tx = []
            td2 = []
            n1 = r
            for n2 in range(elementsCountUpLV, n2b, -1):
                tx .append(lvx [0][n2][n1])
                td2.append([ -d for d in lvd2[1][n2][n1] ])
                #print('  down', n2, n1, td2[-1])
            n2 = n2a
            for n1 in range(n1b, m1a):
                tx .append(lvx [0][n2][n1])
                if n1 < elementsCountAroundLVFreeWallHalf:
                    td2.append([ -d for d in lvd2[1][n2][n1] ])
                elif n1 > (elementsCountAroundLVFreeWall - elementsCountAroundLVFreeWallHalf):
                    td2.append(lvd2[1][n2][n1])
                else:
                    td2.append(lvd1[1][n2][n1])
                #print('across', n2, n1, td2[-1])
            n1 = elementsCountAroundLVFreeWall - r
            for n2 in range(n2c, elementsCountUpLV + 1):
                tx .append(lvx [0][n2][n1])
                td2.append(lvd2[1][n2][n1])
                #print('    up', n2, n1, td2[-1])
            #print(len(td2), td2)
            td2 = smoothCubicHermiteDerivativesLine(tx, td2)
            n1 = r
            for n2 in range(elementsCountUpLV, n2b, -1):
                lvd2[0][n2][n1] = [ -d for d in td2.pop(0) ]
            n2 = n2a
            for n1 in range(n1b, m1a):
                if n1 < elementsCountAroundLVFreeWallHalf:
                    lvd2[0][n2][n1] = [ -d for d in td2.pop(0) ]
                elif n1 > (elementsCountAroundLVFreeWall - elementsCountAroundLVFreeWallHalf):
                    lvd2[0][n2][n1] = td2.pop(0)
                else:
                    #print('n1',n1,'middle', td2[0])
                    lvd1[0][n2][n1] = td2.pop(0)
            n1 = elementsCountAroundLVFreeWall - r
            for n2 in range(n2c, elementsCountUpLV + 1):
                lvd2[0][n2][n1] = td2.pop(0)

        # get inner d1, d2
        # row 1
        lvd1[0][n2b][n1b:-n1b] = smoothCubicHermiteDerivativesLine(lvx[0][n2b][n1b:-n1b], lvd1[1][n2b][n1b:-n1b], fixAllDirections = True)
        # regular rows 2+
        for n2 in range(n2c, elementsCountUpLV + 1):
            lvd1[0][n2] = smoothCubicHermiteDerivativesLine(lvx[0][n2], lvd1[1][n2], fixAllDirections = True)
        # column 1 and -2
        # regular columns
        for n1 in range(n1b, elementsCountAroundLVFreeWall + 1 - n1b):
            startn2 = n2b if (n1 in [n1b, elementsCountAroundLVFreeWall - n1b]) else n2a
            tx  = []
            td2 = []
            if startn2 == n2a:
                tx .append(lvx [0][n2a][n1])
                left = n1 < elementsCountAroundLVFreeWallHalf
                right = n1 > (elementsCountAroundLVFreeWall - elementsCountAroundLVFreeWallHalf)
                if left:
                    td2.append(lvd1[1][n2a][n1])
                elif right:
                    td2.append([ -d for d in lvd1[1][n2a][n1] ])
                else:  # middle:
                    td2.append(lvd2[1][n2a][n1])
            for n2 in range(n2b, elementsCountUpLV + 1):
                tx .append(lvx [0][n2][n1])
                td2.append(lvd2[1][n2][n1])
            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixAllDirections = True)
            if startn2 == n2a:
                if left:
                    lvd1[0][n2a][n1] = td2[0]
                elif right:
                    lvd1[0][n2a][n1] = [ -d for d in td2[0] ]
                else:  # middle
                    #lvd1[0][n2a][n1] = lvd2[0][n2a][n1]
                    lvd2[0][n2a][n1] = td2[0]
            for n2 in range(n2b, elementsCountUpLV + 1):
                lvd2[0][n2][n1] = td2[n2 - startn2]

        # fix inner derivatives leading to triple points
        lvd1[0][n2a][n1b] = smoothCubicHermiteDerivativesLine([ lvx[0][n2a][n1b], lvx[0][n2b][n1b] ], [ lvd1[1][n2a][n1b], [ (lvd1[0][n2b][n1b][c] + lvd2[0][n2b][n1b][c]) for c in range(3) ] ],
                                                              fixEndDerivative = True, fixStartDirection = True)[0]
        lvd1[0][n2a][m1b] = smoothCubicHermiteDerivativesLine([ lvx[0][n2b][m1b], lvx[0][n2a][m1b] ], [ [ (lvd1[0][n2b][m1b][c] - lvd2[0][n2b][m1b][c]) for c in range(3) ], lvd1[1][n2a][m1b] ],
                                                              fixStartDerivative = True, fixEndDirection = True)[1]

        #################
        # Create nodes
        #################

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        nodeIdentifier = startNodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        if False:
            for n in range((elementsCountUpLVTrackSurface + 1)*elementsCountAroundLVTrackSurface):
                node = nodes.createNode(nodeIdentifier, nodetemplate12)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n])
                nodeIdentifier += 1

        if False:
            for n in range(len(lscx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lscx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lscd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lscd2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lscd3[n])
                nodeIdentifier += 1

        if False:
            for n in range(len(px)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, px [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pd2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, pd3[n])
                nodeIdentifier += 1

        if False:
            for n in range(2):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3[n])
                nodeIdentifier += 1

        if False:
            for n in range(len(rvix)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rvix [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rvid1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rvid2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rvid3[n])
                nodeIdentifier += 1

        if False:
            for n in range(len(rvox)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rvox [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rvod1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rvod2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rvod3[n])
                nodeIdentifier += 1

        if False:
            for n in range(elementsCountUpRVFreeWall + 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rscx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rscd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rscd2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rscd3[n])
                nodeIdentifier += 1

        if True:
            for n in range(1, elementsCountUpLVApex + 2):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lax [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lad1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lad2[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lad3[n])
                nodeIdentifier += 1

        #print("LV nodes")
        nodeIdentifier = lvShield.generateNodes(fieldmodule, coordinates, nodeIdentifier)
        #print("RV nodes")
        nodeIdentifier = rvShield.generateNodes(fieldmodule, coordinates, nodeIdentifier)

        #################
        # Create elements
        #################

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        #print("LV elements")
        elementIdentifier = lvShield.generateElements(fieldmodule, coordinates, elementIdentifier, [ lvMeshGroup ])
        #print("RV elements")
        elementIdentifier = rvShield.generateElements(fieldmodule, coordinates, elementIdentifier, [ rvMeshGroup ])

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        Stops at end of ventricles, hence can be called from ventriclesbase.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAroundFull, elementsCountUpLVApex, \
            elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
            elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
            elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall = cls.getElementsCounts(options)
        refineElementsCountSurface = options['Refine number of elements surface']
        #refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        lvGroup = findAnnotationGroupByName(meshrefinement._sourceAnnotationGroups, "left ventricle myocardium")
        rvGroup = findAnnotationGroupByName(meshrefinement._sourceAnnotationGroups, "right ventricle myocardium")
        vSeptumGroup = findAnnotationGroupByName(meshrefinement._sourceAnnotationGroups, "interventricular septum")
        lvMeshGroup = lvGroup.getMeshGroup(meshrefinement._sourceMesh)
        rvMeshGroup = rvGroup.getMeshGroup(meshrefinement._sourceMesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(meshrefinement._sourceMesh)
        elementsCountVentricles = refineElementsCountThroughWall*(
            (elementsCountUpLVFreeWall - elementsCountUpLVApex - 2)*elementsCountAroundLVFreeWall + (elementsCountUpLVApex + 2)*(elementsCountAroundLVFreeWall - 2) +
            (elementsCountUpRVFreeWall - 2)*elementsCountAroundRVFreeWall + 2*(elementsCountAroundRVFreeWall - 2) +
            (elementsCountUpIVSeptum - 2)*elementsCountAcrossIVSeptum + 2*(elementsCountAcrossIVSeptum - 2) +
            elementsCountAroundFull)
        element = meshrefinement._sourceElementiterator.next()
        lastVentriclesElementIdentifier = element.getIdentifier() + elementsCountVentricles - 1
        while element.isValid():
            elementIdentifier = element.getIdentifier()
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            numberInXi3 = refineElementsCountThroughWall
            #if lvMeshGroup.containsElement(element):
            #    numberInXi3 = refineElementsCountThroughLVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementIdentifier == lastVentriclesElementIdentifier:
                return  # finish on last so can continue in ventriclesbase
            element = meshrefinement._sourceElementiterator.next()

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
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()
