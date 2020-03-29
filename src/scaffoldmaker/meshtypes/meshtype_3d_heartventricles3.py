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
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import createEllipsoidPoints, getApproximateEllipsePerimeter, getEllipseArcLength, getEllipseRadiansToX
from scaffoldmaker.utils.interpolation import computeCubicHermiteDerivativeScaling, getCubicHermiteArcLength, interpolateSampleCubicHermite, sampleCubicHermiteCurves, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
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

    @staticmethod
    def checkOptions(options):
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

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAcrossIVSeptum = options['Number of elements across septum']
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        #elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountUpLVFreeWall = options['Number of elements up LV free wall']
        elementsCountUpLVApex = options['Number of elements up LV apex']
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

        nas = options['Number of elements across septum']//2
        nal = options['Number of elements around LV free wall']//2 - options['Number of elements up LV apex']
        nar = options['Number of elements around RV free wall']//2
        nul = options['Number of elements up LV free wall']
        # number around half = must be matched on RV and septum
        nah = nal + nul - 2
        elementsCountAroundFull = 2*nah + (options['Number of elements across septum'] % 2)
        # number up rv, min 2
        elementsCountUpRVFreeWall = nah + 2 - nar
        print("elementsCountAroundFull", elementsCountAroundFull)
        print("elementsCountUpRVFreeWall", elementsCountUpRVFreeWall)

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

        annotationGroups = []
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

        rvix, rvid1, rvid2, rvid3, rviProportions = lvTrackSurface.createHermiteCurvePoints(rvInletPositionAroundLV, baseProportionUp, rvApexProportion1, rvApexProportion2,
            elementsCountAroundHalf, derivativeStart = rvInletDerivative, derivativeEnd = rvApexLengthDerivative, curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
        rvox, rvod1, rvod2, rvod3, rvoProportions = lvTrackSurface.createHermiteCurvePoints(rvApexProportion1, rvApexProportion2, rvOutletPositionAroundLV, baseProportionUp,
            elementsCountAroundHalf, derivativeStart = rvApexLengthDerivative, derivativeEnd = rvOutletDerivative, curveMode = TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
        rvix  += rvox [1:]
        rvid1 += rvod1[1:]
        rvid2 += rvod2[1:]
        rvid3 += rvod3[1:]
        rviProportions += rvoProportions[1:]

        print('elementsCountAroundHalf',elementsCountAroundHalf)
        print('elementsCountUpRVRegular',elementsCountUpRVRegular)
        print('elementsCountAroundRVFreeWallHalf',elementsCountAroundRVFreeWallHalf)
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
        rscx, rscd2 = sampleCubicHermiteCurves(rcx, rcd1, elementsCountUpRVFreeWall, arcLengthDerivatives = True)[0:2]
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
        p2 = lvTrackSurface.trackVector(rvApexPosition, rvApexCuspDirection, ivSulcusApexTransitionLength)
        proportion2 = lvTrackSurface.getProportion(p2)[1]
        x, sd1, sd2 = lvTrackSurface.evaluateCoordinates(p2, derivatives = True)
        lax = [ rscx[0] ]
        lad2 = [ [ -s for s in rscd2[0] ] ]
        if (elementsCountUpLVApex > 0) and (proportion2 > 0.0):
            lax.append(x)
            lad2.append(vector.setMagnitude(sd2, -ivSulcusApexTransitionLength))
        lax.append(nx[0])
        lad2.append(nd2[0])

        if elementsCountUpLVApex > 0:
            lax, lad2 = sampleCubicHermiteCurves(lax, lad2, elementsCountUpLVApex + 1,
                addLengthStart = vector.magnitude(rscd2[0])*0.5, lengthFractionStart = 0.5, arcLengthDerivatives = True)[0:2]
        else:
            #print('lax', lax)
            #print('lad2', lad2)
            lad2 = smoothCubicHermiteDerivativesLine(lax, lad2, fixStartDerivative = True, fixEndDirection = True)
        #print('mag rv apex', vector.magnitude(rscd2[0]), vector.magnitude(lad2[0]))
        #print('mag lv apex', vector.magnitude(lad2[-1]), vector.magnitude(pd1[0]))

        rvx  = [ [], [] ]
        rvd1 = [ [], [] ]
        rvd2 = [ [], [] ]
        rvd3 = [ [], [] ]
        for n3 in range(2):
            for n2 in range(elementsCountUpRVFreeWall + 1):
                for rv in [ rvx[n3], rvd1[n3], rvd2[n3], rvd3[n3] ]:
                    rv.append([ None ]*(elementsCountAroundRVFreeWall + 1))
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
                                     arcLengthDerivatives = True)
            rvd2[1][n2] = interpolateSampleCubicHermite([ rvd2[1][n2][0], rscd2[n2], rvd2[1][n2][-1] ], [ [ 0.0, 0.0, 0.0 ] ]*3, pe, pxi, psf)[0]

        # up regular columns of RV
        for n1 in range(2, elementsCountAroundRVFreeWall - 1):
            left = n1 < elementsCountAroundRVFreeWallHalf
            right = n1 > (elementsCountAroundRVFreeWall - elementsCountAroundRVFreeWallHalf)
            startd2 = rvd1[1][0][n1] if (not right) else [ -d for d in rvd1[1][0][n1] ]
            tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(
                [ rvx[1][0][n1], rvx[1][2][n1] ], [ startd2, rvd2[1][2][n1] ], 2, arcLengthDerivatives = True)
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

        # RV triple points
        ltx = []
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][0][2], rvx[1][2][1] ], [ [ (rvd1[1][0][2][c] + rvd2[1][0][2][c]) for c in range(3) ], rvd2[1][2][1] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][0][1], rvx[1][2][2] ], [ rvd1[1][0][1], [ (rvd1[1][2][2][c] + rvd2[1][2][2][c]) for c in range(3) ] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][2][0], rvx[1][1][2] ], [ [ (rvd1[1][2][0][c] - rvd2[1][2][0][c]) for c in range(3) ], rvd1[1][1][2] ], 2, arcLengthDerivatives = True)[0:2]
        ltx.append(tx[1])
        ltx.append(rvx[1][2][0])
        rtx = []
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][0][-3], rvx[1][2][-2] ], [ [ (-rvd1[1][0][-3][c] + rvd2[1][0][-3][c]) for c in range(3) ], rvd2[1][2][-2] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][0][-2], rvx[1][2][-3] ], [ [ -d for d in rvd1[1][0][-2] ], [ (-rvd1[1][2][-3][c] + rvd2[1][2][-3][c]) for c in range(3) ] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        tx, td1 = sampleCubicHermiteCurves(
            [ rvx[1][2][-1], rvx[1][1][-3] ], [ [ (-rvd1[1][2][-1][c] - rvd2[1][2][-1][c]) for c in range(3) ], [ -d for d in rvd1[1][1][-3] ] ], 2, arcLengthDerivatives = True)[0:2]
        rtx.append(tx[1])
        rtx.append(rvx[1][2][0])

        rvx [1][1][1] = [ (ltx[0][c] + ltx[1][c] + ltx[2][c])/3.0 for c in range(3) ]
        rvd1[1][1][1] = [ (rvx[1][1][2][c] - rvx[1][1][1][c]) for c in range(3) ]
        rvd2[1][1][1] = [ (rvx[1][2][1][c] - rvx[1][1][1][c]) for c in range(3) ]
        rvx [1][1][-2] = [ (rtx[0][c] + rtx[1][c] + rtx[2][c])/3.0 for c in range(3) ]
        rvd1[1][1][-2] = [ (rvx[1][1][-2][c] - rvx[1][1][-3][c]) for c in range(3) ]
        rvd2[1][1][-2] = [ (rvx[1][2][-2][c] - rvx[1][1][-2][c]) for c in range(3) ]

        # smooth RV row 1
        tx = []
        td1 = []
        for n1 in range(1, elementsCountAroundRVFreeWall):
            tx .append(rvx [1][1][n1])
            td1.append(rvd1[1][1][n1])
        td1 = smoothCubicHermiteDerivativesLine(tx, td1)
        for n in range(elementsCountAroundRVFreeWall - 1):
            rvd1[1][1][n + 1] = td1[n]

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
        rvd1[1][0][-2] = smoothCubicHermiteDerivativesLine([ rvx[1][1][-2], rvx[1][0][-2] ], [ [ (rvd1[1][1][-2][c] - rvd2[1][1][-2][c]) for c in range(3) ], rvd1[1][0][-2] ],
                                                          fixStartDerivative = True, fixEndDirection = True)[1]

        # get outer d3 and inner x, d3
        for n2 in range(elementsCountUpRVFreeWall + 1):
            for n1 in range(elementsCountAroundLVFreeWall + 1):
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
        rvd1[0][0][-2] = smoothCubicHermiteDerivativesLine([ rvx[0][1][-2], rvx[0][0][-2] ], [ [ (rvd1[0][1][-2][c] - rvd2[0][1][-2][c]) for c in range(3) ], rvd1[0][0][-2] ],
                                                          fixStartDerivative = True, fixEndDirection = True)[1]

        # LV free wall
        elementsCountUpLV = elementsCountUpLVFreeWall + elementsCountUpLVApex
        lvx  = [ [], [] ]
        lvd1 = [ [], [] ]
        lvd2 = [ [], [] ]
        lvd3 = [ [], [] ]
        lvProportions = []
        for n2 in range(elementsCountUpLV + 1):
            for n3 in range(2):
                for lv in [ lvx[n3], lvd1[n3], lvd2[n3], lvd3[n3] ]:
                    lv.append([ None ]*(elementsCountAroundLVFreeWall + 1))
            lvProportions.append([ None ]*(elementsCountAroundLVFreeWall + 1))

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
            print('n2',n2,'ix',ix,'ox',ox)
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
            #print('hx', hx)
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
                addLengthEnd=0.5*vector.magnitude(rd1[0]), lengthFractionEnd=0.5,
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
        n1 = elementsCountUpLVApex
        print('derivativeStart', vector.setMagnitude(nd2[elementsCountAroundLVTrackSurface*3//4], -vector.magnitude(lad2[-1])))
        print('derivativeEnd', lvd2[1][n2reg][n1])
        print('n2reg', n2reg)
        print('lvProportions[n2reg]', lvProportions[n2reg])
        elementsCountRemaining = elementsCountAroundHalf - (elementsCountUpLVFreeWall - 2)
        print('elementsCountRemaining',elementsCountRemaining)
        elementsCountAroundLVFreeWallHalf = elementsCountAroundLVFreeWall//2
        tx, td2, td1, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(0.75, 0.0, lvProportions[n2reg][n1][0], lvProportions[n2reg][n1][1], elementsCountRemaining,
            derivativeStart=vector.setMagnitude(nd2[elementsCountAroundLVTrackSurface*3//4], vector.magnitude(lad2[-1])), derivativeEnd=lvd2[1][n2reg][n1])
        n2 = elementsCountUpLVApex
        for n in range(elementsCountRemaining):
            n1 = elementsCountAroundLVFreeWallHalf - n
            lvx [1][n2][n1] = tx [n]
            lvd1[1][n2][n1] = [ -d for d in td1[n] ] if n else [ -d for d in td2[n] ]
            lvd2[1][n2][n1] = td2[n] if n else lad2[-1]
            lvd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(lvd1[1][n2][n1], lvd2[1][n2][n1]), lvFreeWallThickness)
            lvProportions[n2][n1] = tProportions[n]
        n1 = -1 - elementsCountUpLVApex
        tx, td2, td1, td3, tProportions = lvTrackSurface.createHermiteCurvePoints(1.25, 0.0, lvProportions[n2reg][n1][0], lvProportions[n2reg][n1][1], elementsCountRemaining,
            derivativeStart=lvd2[1][n2][elementsCountAroundLVFreeWallHalf], derivativeEnd=lvd2[1][n2reg][n1])
        for n in range(1, elementsCountRemaining):
            n1 = elementsCountAroundLVFreeWallHalf + n
            lvx [1][n2][n1] = tx [n]
            lvd1[1][n2][n1] = [ -d for d in td1[n] ]
            lvd2[1][n2][n1] = td2[n]
            lvd3[1][n2][n1] = vector.setMagnitude(vector.crossproduct3(lvd1[1][n2][n1], lvd2[1][n2][n1]), lvFreeWallThickness)
            lvProportions[n2][n1] = tProportions[n]

        #################
        # Create nodes
        #################

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # LV/RV outlet elements are linear through the wall, hence their nodes do not have D_DS3 parameters
        nodetemplate12 = nodes.createNodetemplate()
        nodetemplate12.defineField(coordinates)
        nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = startNodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

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
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lad1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lad2[n])
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lad3[n])
                nodeIdentifier += 1

        lvNodeId = [ [], [] ]
        if True:
            for n2 in range(elementsCountUpLV + 1):
                for n3 in range(2):
                    lvNodeId[n3].append([ None ]*(elementsCountAroundLVFreeWall + 1))
                    for n1 in range(elementsCountAroundLVFreeWall + 1):
                        if lvx[n3][n2][n1]:
                            node = nodes.createNode(nodeIdentifier, nodetemplate)
                            lvNodeId[n3][n2][n1] = nodeIdentifier
                            cache.setNode(node)
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lvx [n3][n2][n1])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lvd1[n3][n2][n1])
                            if lvd2[n3][n2][n1]:  # GRC temp
                                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lvd2[n3][n2][n1])
                                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lvd3[n3][n2][n1])
                            nodeIdentifier += 1

        rvNodeId = [ [], [] ]
        if True:
            for n2 in range(elementsCountUpRVFreeWall + 1):
                for n3 in range(2):
                    rvNodeId[n3].append([ None ]*(elementsCountAroundRVFreeWall + 1))
                    for n1 in range(elementsCountAroundLVFreeWall + 1):
                        if rvx[n3][n2][n1]:
                            node = nodes.createNode(nodeIdentifier, nodetemplate)
                            rvNodeId[n3][n2][n1] = nodeIdentifier
                            cache.setNode(node)
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rvx [n3][n2][n1])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rvd1[n3][n2][n1])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rvd2[n3][n2][n1])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rvd3[n3][n2][n1])
                            nodeIdentifier += 1

        if False:
            # RV triple points
            for n in range(len(ltx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ltx[n])
                nodeIdentifier += 1
            for n in range(len(rtx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rtx[n])
                nodeIdentifier += 1

        #################
        # Create elements
        #################

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        elementIdentifier = getMaximumElementIdentifier(mesh) + 1

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        isEven = (elementsCountAroundRVFreeWall % 2) == 0
        elementsCountAroundRVFreeWallPlusOneHalf = (elementsCountAroundRVFreeWall + 1)//2
        for e2 in range(elementsCountUpRVFreeWall):

            meshGroups = [ rvMeshGroup ]
            for e1 in range(elementsCountAroundRVFreeWall):
                eft1 = eft
                scalefactors = None
                nids = [ rvNodeId[0][e2][e1], rvNodeId[0][e2][e1 + 1], rvNodeId[0][e2 + 1][e1], rvNodeId[0][e2 + 1][e1 + 1], 
                         rvNodeId[1][e2][e1], rvNodeId[1][e2][e1 + 1], rvNodeId[1][e2 + 1][e1], rvNodeId[1][e2 + 1][e1 + 1] ]
                if e2 == 0:
                    if (e1 == 0) or (e1 == (elementsCountAroundRVFreeWall - 1)):
                        continue
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    if e1 < elementsCountAroundRVFreeWallPlusOneHalf:
                        if e1 < (elementsCountAroundRVFreeWallPlusOneHalf - 1):
                            if e1 == 1:
                                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    if e1 >= elementsCountAroundRVFreeWallPlusOneHalf:
                        if e1 > elementsCountAroundRVFreeWallPlusOneHalf:
                            if e1 == (elementsCountAroundRVFreeWall - 2):
                                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                elif e2 == 1:
                    if e1 == 0:
                        nids[0] = rvNodeId[0][0][1]
                        nids[4] = rvNodeId[1][0][1]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    elif e1 == (elementsCountAroundRVFreeWall - 1):
                        nids[1] = rvNodeId[0][0][-2]
                        nids[5] = rvNodeId[1][0][-2]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])

                if eft1 is not eft:
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                else:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                #print('create element rv', elementIdentifier, result2, result3, nids)
                elementIdentifier = elementIdentifier + 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

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
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through wall']

        startPostApexElementIdentifier = elementsCountAroundLV*elementsCountUpLVApex + 1
        elementsCountPostApexLayer = elementsCountAroundLV + 2 + elementsCountAroundVSeptum
        lastVentriclesElementIdentifier = startPostApexElementIdentifier + elementsCountAroundVSeptum + elementsCountUpRV*elementsCountPostApexLayer - 1

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            numberInXi3 = refineElementsCountThroughLVWall
            elementIdentifier = element.getIdentifier()
            if elementIdentifier >= startPostApexElementIdentifier:
                n1 = (elementIdentifier - startPostApexElementIdentifier - elementsCountAroundVSeptum)
                if n1 < 0:
                    # collapsed elements at RV apex
                    numberInXi2 = refineElementsCountThroughLVWall
                else:
                    n1 = n1 % elementsCountPostApexLayer
                    if (n1 == 0) or (n1 == (elementsCountAroundLVFreeWall + 1)):
                        # collapsed elements on posterior or anterior interventricular sulcus:
                        numberInXi1 = refineElementsCountThroughLVWall
                        #print(n1,'refine collapse element', elementIdentifier, numberInXi1, numberInXi2, numberInXi3)
                    elif (n1 > (elementsCountAroundLVFreeWall + 1)) and (n1 <= (elementsCountAroundLV + 1)):
                        numberInXi3 = refineElementsCountThroughRVWall
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
