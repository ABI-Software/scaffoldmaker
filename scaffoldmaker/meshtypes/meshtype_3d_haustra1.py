"""
Generates a single 3-D haustra segment mesh along a central
line, with variable numbers of elements around, along and
through wall, with variable radius and thickness along.
"""

import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear # Remove
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite # Remove
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils # Remove
from opencmiss.zinc.element import Element, Elementbasis # Remove
from opencmiss.zinc.field import Field # Remove
from opencmiss.zinc.node import Node # Remove

class MeshType_3d_haustra1(Scaffold_base):
    '''
    Generates a single 3-D haustra segment mesh with variable
    numbers of elements around, along the central line, and
    through wall. The haustra segment has a triangular profile
    with rounded corners at the inter-haustral septa, and a
    clover profile in the intra-haustral region.
    '''
    @staticmethod
    def getName():
        return '3D Haustra 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around tenia coli' : 2,
            'Number of elements around haustrum' : 8,
            'Number of elements along haustrum' : 3,
            'Number of elements through wall' : 1,
            'Inner radius': 0.5,
            'Corner inner radius factor': 0.5,
            'Haustrum inner radius factor': 0.5,
            'Haustrum length end derivative factor': 0.5,
            'Haustrum length mid derivative factor': 1.0,
            'Haustrum length': 1.0,
            'Tenia coli width': 0.2,
            'Wall thickness': 0.01,
            'Use cross derivatives' : False,
            'Use linear through wall' : True,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along haustrum' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along haustrum',
            'Number of elements through wall',
            'Inner radius',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Haustrum length end derivative factor',
            'Haustrum length mid derivative factor',
            'Haustrum length',
            'Tenia coli width',
            'Wall thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along haustrum',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along haustrum',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements along haustrum']:
            if options[key] < 2:
                options[key] = 2
        if options['Number of elements around haustrum'] < 4:
            options['Number of elements around haustrum'] = 4
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum']:
            if options[key] % 2 > 0:
                options[key] = options[key] + 1
        for key in [
            'Inner radius',
            'Haustrum inner radius factor',
            'Haustrum length end derivative factor',
            'Haustrum length mid derivative factor',
            'Haustrum length',
            'Wall thickness'
            ]:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner inner radius factor'] < 0.1:
            options['Corner inner radius factor'] = 0.1
        for key in [
            'Corner inner radius factor',
            'Haustrum length end derivative factor']:
            if options[key] > 1.0:
                options[key] = 1.0
        if options['Tenia coli width'] < 0.2*options['Inner radius']:
            options['Tenia coli width'] = 0.2*options['Inner radius']
        if options['Tenia coli width'] > round(math.sqrt(3)*0.5*options['Inner radius'],2):
            options['Tenia coli width'] = round(math.sqrt(3)*0.5*options['Inner radius'],2)

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAroundTC = options['Number of elements around tenia coli']
        elementsCountAroundHaustrum = options['Number of elements around haustrum']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*3
        elementsCountAlongHaustrum = options['Number of elements along haustrum']
        elementsCountThroughWall = options['Number of elements through wall']
        radius = options['Inner radius']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        haustrumLengthEndDerivativeFactor = options['Haustrum length end derivative factor']
        haustrumLengthMidDerivativeFactor = options['Haustrum length mid derivative factor']
        haustrumLength = options['Haustrum length']
        widthTC = options['Tenia coli width']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = elementsCountAlongHaustrum
        haustraSegmentCount = 1

        cx = [ [ 0.0, 0.0, 0.0 ], [ haustrumLength, 0.0, 0.0 ] ]
        cd1 = [ [ haustrumLength, 0.0, 0.0 ], [ haustrumLength, 0.0, 0.0 ] ]

        # Generate inner surface of a haustra segment
        xHaustraInner, d1HaustraInner, d2HaustraInner, haustraSegmentAxis = getColonHaustraSegmentInnerPoints(elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
            haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength)
        #haustraSegmentAxis = getColonHaustraSegmentInnerPoints(region, useCubicHermiteThroughWall, useCrossDerivatives, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
        #    haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier = tubemesh.generatetubemesh(region, elementsCountAround, elementsCountAlongHaustrum, elementsCountThroughWall, haustraSegmentCount,
            cx, cd1, xHaustraInner, d1HaustraInner, d2HaustraInner, wallThickness, haustraSegmentAxis, haustrumLength, useCrossDerivatives, useCubicHermiteThroughWall)

        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()

def getColonHaustraSegmentInnerPoints(elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
    haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength):
#def getColonHaustraSegmentInnerPoints(region, useCubicHermiteThroughWall, useCrossDerivatives, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
#    haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength):
    """
    Generates a 3-D haustra segment mesh with variable numbers
    of elements around, along the central line, and through wall.
    The haustra segment has a triangular profile with rounded corners
    at the inter-haustral septa, and a clover profile in the intra-haustral
    region.
    :param elementsCountAroundTC: Number of elements around each tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlongHaustrum: Number of elements along haustrum.
    :param widthTC: Width of tenia coli.
    :param radius: Inner radius defined from center of triangular
    profile to vertex of the triangle.
    :param cornerInnerRadiusFactor: Roundness of triangular corners of
    interhaustral septa. Factor is multiplied by inner radius
    to get a radius of curvature at the corners.
    :param haustrumInnerRadiusFactor: Factor is multiplied by inner
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra segment.
    :param haustrumLengthEndDerivativeFactor: Factor is multiplied by haustrum
    length to scale derivative along the end of a haustrum length.
    :param haustrumLengthMidDerivativeFactor: Factor is multiplied by haustrum
    length to scale derivative along the mid length of the haustrum.
    :param haustrumLength: Length of a haustrum.
    :return: coordinates, derivatives on inner surface of haustra segment.
    """

########################### Delete later ###########################################
    # fm = region.getFieldmodule()
    # fm.beginChange()
    # cache = fm.createFieldcache()
    # coordinates = zinc_utils.getOrCreateCoordinateField(fm)

    # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    # nodetemplate = nodes.createNodetemplate()
    # nodetemplate.defineField(coordinates)
    # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    # if useCrossDerivatives:
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    # if useCubicHermiteThroughWall:
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # if useCrossDerivatives:
            # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

    # mesh = fm.findMeshByDimension(3)

    # if useCubicHermiteThroughWall:
        # eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    # else:
        # eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    # eft = eftfactory.createEftBasic()

    # elementtemplate = mesh.createElementtemplate()
    # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    # result = elementtemplate.defineField(coordinates, -1, eft)

    # nodeIdentifier = 1
##################################################################################

    # create nodes
    zero = [0.0, 0.0, 0.0] # Delete later if not in use
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    cornerRC = cornerInnerRadiusFactor*radius
    radiansRangeRC = [7*math.pi/4, 0.0, math.pi/4]
    sampleElementOut = 20
    unitZ = [0.0, 0.0, 1.0]
    haustrumRadius = (haustrumInnerRadiusFactor + 1)*radius

    xHalfSetInterHaustra = []
    d1HalfSetInterHaustra = []
    nxHaustrum = []
    nd1Haustrum = []
    xAround = []
    d1Around = []
    # dx_ds1InnerList = []
    xInnerRaw = []
    dx_ds2InnerRaw = []
    xInnerList = []
    dx_ds2InnerList = []
    xFinal = []
    d1Final = []
    d2Final = []

    # Inter-haustral segment
    # Set up profile
    for n1 in range(3):
        radiansAround = n1*2*math.pi / 3
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        xc = [(radius - cornerRC) * cosRadiansAround, (radius - cornerRC) * sinRadiansAround, 0.0]

        for n in range(3):
            radiansRC = radiansAround + radiansRangeRC[n]
            cosRadiansRC = math.cos(radiansRC)
            sinRadiansRC = math.sin(radiansRC)
            x = [xc[0] + cornerRC*cosRadiansRC, xc[1] + cornerRC*sinRadiansRC, 0.0]
            xAround.append(x)
            d1 = [ cornerRC*math.pi/4 * -sinRadiansRC, cornerRC*math.pi/4 * cosRadiansRC, 0.0]
            d1Around.append(d1)

    xSample = xAround[1:9] +[xAround[0], xAround[1]]
    d1Sample = d1Around[1:9] +[d1Around[0], d1Around[1]]
    sx, sd1, se, sxi, _= interp.sampleCubicHermiteCurves(xSample, d1Sample, sampleElementOut)
    xLoop = sx[0:-1]
    d1Loop = interp.smoothCubicHermiteDerivativesLoop(sx[0:-1], sd1[0:-1])

    # Calculate arc length
    arcLength = 0.0
    arcDistance = interp.getCubicHermiteArcLength(xLoop[0], d1Loop[0], xLoop[1], d1Loop[1])
    arcLength = arcDistance * sampleElementOut
    arcStart = 0.0
    arcEnd = arcLength/6.0

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(xLoop, d1Loop, widthTC, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC = sampleTeniaColi(xLoop, d1Loop, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum = sampleHaustrum(xLoop, d1Loop, xTC[-1], d1TC[-1], arcLength/6.0, arcDistanceTCEdge, elementsCountAroundHaustrum)

    xHalfSetInterHaustra = xHalfSetInterHaustra + xTC + xHaustrum[1:]
    d1HalfSetInterHaustra = d1HalfSetInterHaustra + d1TC + d1Haustrum[1:]

    # Haustra segment
    # Re-initialise for haustra segment
    xHalfSetIntraHaustra = []
    d1HalfSetIntraHaustra = []

    # Set up profile
    xc = [(radius - cornerRC)* math.cos(0.0), (radius - cornerRC)*math.sin(0.0), 0.0]
    pt1 = [xc[0] + cornerRC*math.cos(0.0), xc[1] + cornerRC*math.sin(0.0), 0.0]
    xTC2 = radius* math.cos(2*math.pi/3)
    yTC2 = radius* math.sin(2*math.pi/3)
    originRC = (xTC2*xTC2 + yTC2*yTC2 - haustrumRadius*haustrumRadius) / (2*(-xTC2 - haustrumRadius))
    RC = haustrumRadius - originRC

    # Rotate to find originRC of 1st haustrum
    yTC1 = pt1[1]
    rotOriginRC = [ originRC*math.cos(-2/3*math.pi), originRC*math.sin(-2/3*math.pi), 0.0]

    thetaStart = math.asin((yTC1 + rotOriginRC[1]) / RC)
    thetaEnd = math.pi - math.asin((yTC2 + rotOriginRC[1])/ RC)
    thetaHalfHaustrum = (thetaEnd + thetaStart)*0.5
    concaveFactor = 0.15
    thetaConcave = (thetaEnd - thetaStart)*concaveFactor + thetaStart
    thetaCircularStart = (thetaEnd - thetaStart)*(concaveFactor + 0.5)*0.5 + thetaStart

    thetaSet = [thetaStart, thetaConcave]
    xConcave, _ = getCircleXandD1FromRadians(thetaSet, RC, rotOriginRC)
    d1 = [0.0, xConcave[1][1], 0.0]
    nxHaustrum = nxHaustrum + xConcave
    nd1Haustrum = nd1Haustrum + [d1, d1]
    thetaCircularSet = [thetaCircularStart, thetaHalfHaustrum]
    xCircular, d1Circular = getCircleXandD1FromRadians(thetaCircularSet, RC, rotOriginRC)
    nxHaustrum = nxHaustrum + xCircular
    nd1Haustrum = nd1Haustrum + d1Circular
    smoothd1 = interp.smoothCubicHermiteDerivativesLine(nxHaustrum, nd1Haustrum, fixStartDirection = True, fixEndDirection = True)

    # Find max x along path
    sxHaustrum, sd1Haustrum, _, _ , _ = interp.sampleCubicHermiteCurves(nxHaustrum, smoothd1, sampleElementOut)
    xVal = []
    for i in range(len(sxHaustrum)):
        xVal.append(sxHaustrum[i][0])
    maxValue = max(xVal)
    maxIndex = xVal.index(maxValue)
    yAtMaxValue = sxHaustrum[maxIndex][1]

    arcLength = 0.0
    for e in range(len(nxHaustrum)-1):
        arcDistance = interp.getCubicHermiteArcLength(nxHaustrum[e], smoothd1[e], nxHaustrum[e+1], smoothd1[e+1])
        arcLength = arcLength + arcDistance

    arcDistanceAtMaxValue = arcLength / sampleElementOut * maxIndex
    arcStart = 0.0 if yAtMaxValue > widthTC*0.5 else arcDistanceAtMaxValue
    arcEnd = arcDistanceAtMaxValue if yAtMaxValue > widthTC*0.5 else arcLength

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(nxHaustrum, smoothd1, widthTC, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC = sampleTeniaColi(nxHaustrum, smoothd1, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum = sampleHaustrum(nxHaustrum, smoothd1, xTC[-1], d1TC[-1], arcLength, arcDistanceTCEdge, elementsCountAroundHaustrum)

    xHalfSetIntraHaustra = xHalfSetIntraHaustra + xTC + xHaustrum[1:]
    d1HalfSetIntraHaustra = d1HalfSetIntraHaustra + d1TC + d1Haustrum[1:]

    # Sample arclength of haustra segment
    elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum)*0.5)

    for n1 in range(elementsCountAroundHalfHaustrum + 1):
        v1 = [xHalfSetInterHaustra[n1][0], xHalfSetInterHaustra[n1][1], 0.0]
        startArcLength = haustrumLengthEndDerivativeFactor * haustrumLength
        d1 = [ c*startArcLength for c in unitZ]
        v2 = [xHalfSetIntraHaustra[n1][0], xHalfSetIntraHaustra[n1][1], haustrumLength/2]
        midArcLength = haustrumLengthMidDerivativeFactor * haustrumLength
        d2 = [c*midArcLength for c in unitZ]
        v3 = [xHalfSetInterHaustra[n1][0], xHalfSetInterHaustra[n1][1], haustrumLength]
        d3 = [ c*startArcLength for c in unitZ]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, sd1, se, sxi, _  = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAlongHaustrum)
        xInnerRaw.append(sx)
        dx_ds2InnerRaw.append(sd1)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    for n2 in range(elementsCountAlongHaustrum + 1):
        xAround = []
        unitdx_ds1Around = []
        d2Around = []
        for n1 in range(elementsCountAroundHalfHaustrum+1):
            x = xInnerRaw[n1][n2]
            xInnerList.append(x)
            dx_ds2 = dx_ds2InnerRaw[n1][n2]
            dx_ds2InnerList.append(dx_ds2)
            unitTangent = vector.normalise(dx_ds2)
            # Intra-Haustra
            if n1 == 0:
                unitdx_ds1 = vector.normalise(d1HalfSetIntraHaustra[n1])
            else: # points on clover
                if n2 <= int(elementsCountAlongHaustrum/2): # first half of haustrumLength
                    axisRot = vector.crossproduct3(unitZ, unitTangent)
                elif n2 > int(elementsCountAlongHaustrum/2): # second half of haustrumLength
                    axisRot = vector.crossproduct3(unitTangent, unitZ)
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, math.pi/2)
                rotNormal = [rotFrame[j][0]*unitTangent[0] + rotFrame[j][1]*unitTangent[1] + rotFrame[j][2]*unitTangent[2] for j in range(3)]
                unitdx_ds3 = vector.normalise(rotNormal)
                unitdx_ds1 = vector.crossproduct3(unitTangent, unitdx_ds3)
            xAround.append(x)
            d2Around.append(dx_ds2)
            unitdx_ds1Around.append(unitdx_ds1)

        if n2 > 0 and n2 < elementsCountAlongHaustrum:
            dx_ds1InnerAroundList = []
            if elementsCountAlongHaustrum%2 == 0 and n2 == int(elementsCountAlongHaustrum*0.5):
                dx_ds1InnerAroundList = dx_ds1InnerAroundList + d1HalfSetIntraHaustra
            else:
                for n1 in range(elementsCountAroundHalfHaustrum):
                    v1 = xAround[n1]
                    d1 = unitdx_ds1Around[n1]
                    v2 = xAround[n1+1]
                    d2 = unitdx_ds1Around[n1+1]
                    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                    dx_ds1 = [c*arcLengthAround for c in d1]
                    dx_ds1InnerAroundList.append(dx_ds1)
                # Account for d1 of node on half haustrum
                dx_ds1 = [c*arcLengthAround for c in unitdx_ds1Around[elementsCountAroundHalfHaustrum]]
                dx_ds1InnerAroundList.append(dx_ds1)
            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, dx_ds1InnerAroundList, fixStartDerivative = True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC*0.5)], vector.magnitude(d1Smoothed[int(elementsCountAroundTC*0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC*0.5 + 1)], vector.magnitude(d1Smoothed[int(elementsCountAroundTC*0.5)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC*0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC*0.5 + 2):]
            # dx_ds1InnerList = dx_ds1InnerList + d1Corrected
            # d1Corrected = d1Smoothed # dx_ds1InnerAroundList
        else:
            d1Corrected = d1HalfSetInterHaustra
            #dx_ds1InnerList = dx_ds1InnerList + d1HalfSetInterHaustra

        xAlongList, d1AlongList, d2AlongList = getFullProfileFromHalfHaustrum(xAround, d1Corrected, d2Around)
        xFinal = xFinal + xAlongList
        d1Final = d1Final + d1AlongList
        d2Final = d2Final + d2AlongList

    # ############# Delete later ################################################
    # for n in range(len(xInnerList)):
        # node = nodes.createNode(nodeIdentifier, nodetemplate)
        # cache.setNode(node)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xInnerList[n])
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1InnerList[n])
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        # if useCrossDerivatives:
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        # # print('NodeIdentifier = ', nodeIdentifier, xFinal[n])
        # nodeIdentifier = nodeIdentifier + 1

    # fm.endChange()
# # ###########################################################################

    return xFinal, d1Final, d2Final, unitZ
    #return unitZ

def findEdgeOfTeniaColi(nx, nd1, widthTC, arcStart, arcEnd):
    """
    Locate edge of tenia coli on a cubic hermite interpolated
    curve defined by nodes nx and derivatives nd1.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param widthTC: Width of tenia coli.
    :param arcStart: Lower limit of arc distance to search for
    edge of tenia coli.
    :param arcEnd: Upper limit of arc distance to search for
    edge of tenia coli.
    :return: arc distance covered by tenia coli.
    """
    xTol = 1.0E-6
    for iter in range(100):
        arcDistance = (arcStart + arcEnd)*0.5
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        diff = x[1] - widthTC*0.5
        if abs(diff) > xTol:
            if diff < 0.0:
                arcStart = arcDistance
            else:
                arcEnd = arcDistance
        else:
            # xTCEdge = x
            # d1TCEdge = d1
            arcDistanceTCEdge = arcDistance
            break
    if iter > 99:
        print('Search for TC boundary - Max iters reached:',iter)

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
    :return: coordinates, derivatives on tenia coli.
    """
    xTC = []
    d1TC = []
    arcDistancePerElementTC = arcDistanceTCEdge / (elementsCountAroundTC*0.5)
    for e in range(int(elementsCountAroundTC*0.5)+1):
        arcDistance = arcDistancePerElementTC * e
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = vector.setMagnitude(d1, arcDistancePerElementTC)
        xTC.append(x)
        d1TC.append(d1Scaled)

    return xTC, d1TC

def sampleHaustrum(nx, nd1, xTCLast, d1TCLast, arcLength, arcDistanceTCEdge, elementsCountAroundHaustrum):
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
    """
    xHaustrum = []
    d1Haustrum = []
    elementLengths = []
    length = arcLength - arcDistanceTCEdge
    elementsCountOut = int(elementsCountAroundHaustrum * 0.5)
    addLengthStart = 0.5 * vector.magnitude(d1TCLast)
    lengthFractionStart = 0.5
    proportionEnd = 1.0
    proportionStart = 1.0
    elementLengthMid = (length - addLengthStart) / (elementsCountOut - 1.0 + proportionStart*lengthFractionStart )
    elementLengthProportionStart = proportionStart*lengthFractionStart*elementLengthMid
    elementLengthProportionEnd = elementLengthMid

    for e in range(elementsCountOut):
        xi = e/(elementsCountOut - 1)
        elementLengths.append(elementLengthMid)
    elementLengths[0] = addLengthStart + elementLengthProportionStart
    elementLengths[-1] = 0.0 + elementLengthProportionEnd

    arcDistance = arcDistanceTCEdge
    xHaustrum.append(xTCLast)
    d1Scaled = vector.setMagnitude(d1TCLast, elementLengths[0])
    d1Haustrum.append(d1Scaled)

    for e in range(elementsCountOut):
        arcDistance = arcDistance + elementLengths[e]
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = vector.setMagnitude(d1, elementLengths[e])
        xHaustrum.append(x)
        d1Haustrum.append(d1Scaled)

    return xHaustrum, d1Haustrum

def getFullProfileFromHalfHaustrum(xHaustrumHalfSet, d1HaustrumHalfSet, d2HaustrumHalfSet):
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
    :return: coordinates, derivatives of points over entire profile.
    """
    xHaustrumHalfSet2 = []
    d1HaustrumHalfSet2 = []
    d2HaustrumHalfSet2 = []
    xHaustrum = []
    d1Haustrum = []
    d2Haustrum = []
    xHaustra = []
    d1Haustra = []
    d2Haustra = []

    for n in range(1,len(xHaustrumHalfSet)):
        idx =  -n + len(xHaustrumHalfSet) - 1
        x = xHaustrumHalfSet[idx]
        d1 = d1HaustrumHalfSet[idx]
        xReflect = [x[0], -x[1], x[2]]
        d1Reflect = [d1[0], -d1[1], d1[2]]
        xRot = [xReflect[0]*math.cos(2/3*math.pi) - xReflect[1]*math.sin(2/3*math.pi),
                xReflect[0]*math.sin(2/3*math.pi) + xReflect[1]*math.cos(2/3*math.pi),
                xReflect[2]]
        d1Rot = [-(d1Reflect[0]*math.cos(2/3*math.pi) - d1Reflect[1]*math.sin(2/3*math.pi)),
                -(d1Reflect[0]*math.sin(2/3*math.pi) + d1Reflect[1]*math.cos(2/3*math.pi)),
                d1Reflect[2]]
        d2 = d2HaustrumHalfSet[idx]
        d2Reflect = [d2[0], -d2[1], d2[2]]
        d2Rot = [(d2Reflect[0]*math.cos(2/3*math.pi) - d2Reflect[1]*math.sin(2/3*math.pi)),
                (d2Reflect[0]*math.sin(2/3*math.pi) + d2Reflect[1]*math.cos(2/3*math.pi)),
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

    ang = [ 2/3*math.pi, -2/3*math.pi]
    for i in range(2):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(xHaustrum)- 1):
            x = xHaustrum[n]
            d1 = d1Haustrum[n]
            d2 = d2Haustrum[n]
            x = [ x[0]*cosRotAng - x[1]*sinRotAng, x[0]*sinRotAng + x[1]*cosRotAng, x[2]]
            xHaustra.append(x)
            dx_ds1 = [ d1[0]*cosRotAng - d1[1]*sinRotAng, d1[0]*sinRotAng + d1[1]*cosRotAng, d1[2]]
            d1Haustra.append(dx_ds1)
            dx_ds2 = [ d2[0]*cosRotAng - d2[1]*sinRotAng, d2[0]*sinRotAng + d2[1]*cosRotAng, d2[2]]
            d2Haustra.append(dx_ds2)

    return xHaustra, d1Haustra, d2Haustra

def getCircleXandD1FromRadians(thetaSet, radius, origin):
    """
    Gets the coordinates and derivatives along a circular path
    based on the angular range.
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
        x = [radius*math.cos(theta) - origin[0],
            radius*math.sin(theta) - origin[1],
            0.0]
        d1 = [-radius*math.sin(theta)*dTheta,
            radius*math.cos(theta)*dTheta,
            0.0]
        nx.append(x)
        nd1.append(d1)

    return nx, nd1
