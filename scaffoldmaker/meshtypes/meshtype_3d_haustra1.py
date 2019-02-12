"""
Generates a 3-D haustra mesh along a central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh import *

class MeshType_3d_haustra1(object):
    '''
    Generates a 3-D haustra mesh with variable numbers
    of elements around, along the central line, and through wall.
    The haustra has a triangular profile with rounded corners at the
    interhaustral septa, and a clover profile outside the septa.
    '''
    @staticmethod
    def getName():
        return '3D Haustra 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around' : 9,
            'Number of elements along haustra' : 3,
            'Number of elements through wall' : 1,
            'Inner radius': 0.5,
            'Corner radius fraction': 0.5,
            'Wall thickness': 0.01,
            'Number of haustra': 1,
            'Haustra radius fraction': 0.5,
            'Haustra length': 1.0,
            'Interhaustra fold factor': 0.5,
            'Haustra curvature factor': 1.0,
            'Use cross derivatives' : False,
            'Use linear through wall' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements along haustra',
            'Number of elements through wall',
            'Inner radius',
            'Corner radius fraction',
            'Wall thickness',
            'Number of haustra',
            'Haustra radius fraction',
            'Haustra length',
            'Interhaustra fold factor',
            'Haustra curvature factor',
           'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements through wall',
            'Number of haustra',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements around'] < 9) :
            options['Number of elements around'] = 9
        if (options['Number of elements around'] % 3 != 0) :
            options['Number of elements around'] = (options['Number of elements around']//3)*3
        for key in [
            'Inner radius',
            'Corner radius fraction',
            'Wall thickness',
            'Haustra radius fraction',
            'Haustra length',
            'Interhaustra fold factor',
            'Haustra curvature factor']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner radius fraction'] < 0.1:
            options['Corner radius fraction'] = 0.1
        for key in [
            'Corner radius fraction',
            'Interhaustra fold factor']:
            if options[key] > 1.0:
                options[key] = 1.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountAlongHaustra = options['Number of elements along haustra']
        elementsCountThroughWall = options['Number of elements through wall']
        radius = options['Inner radius']
        cornerRadiusFraction = options['Corner radius fraction']
        wallThickness = options['Wall thickness']
        haustraCount = options['Number of haustra']
        haustraRadiusFraction = options['Haustra radius fraction']
        haustraLength = options['Haustra length']
        foldFactor = options['Interhaustra fold factor']
        haustraCurvatureFactor = options['Haustra curvature factor']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongHaustra*haustraCount)

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        coordinates = getOrCreateCoordinateField(fm)

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
        result = elementtemplate.defineField(coordinates, -1, eft)

        # create nodes
        nodeIdentifier = 1 #nextNodeIdentifier
        zero = [ 0.0, 0.0, 0.0 ]
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        radiansRangeRC = [7*math.pi/4, 0.0, math.pi/4]
        cornerRC = cornerRadiusFraction*radius
        unitZ = [0.0, 0.0, 1.0]
        haustraRadius = (haustraRadiusFraction + 1)*radius
        xAround = []
        d1Around = []
        xInner = []
        dx_ds1InnerList = []
        xHaustraSide = []
        xHaustraInner = []
        d1InnerHaustraRaw = []
        xInnerRaw = []
        dx_ds2InnerRaw = []
        xInnerList = []
        dx_ds2InnerList = []
        dx_ds3InnerUnitList = []
        xList = []
        dx_ds1List = []
        dx_ds2List = []
        dx_ds3List = []

        # Pre-calculate node locations and derivatives on inner triangle
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

        xSample = xAround[1:9]
        xSample.append(xAround[0])
        xSample.append(xAround[1])
        d1Sample = d1Around[1:9]
        d1Sample.append(d1Around[0])
        d1Sample.append(d1Around[1])
        sx, sd1, se, sxi, _= sampleCubicHermiteCurves(xSample, d1Sample, elementsCountAround)
        xInner = xInner + sx[0:-1]
        d1Inner = smoothCubicHermiteDerivativesLoop(sx[0:-1], sd1[0:-1])

        # Pre-calculate node locations and derivatives on haustra inner cross-section
        elementsCountAroundSide = int(elementsCountAround/3)
        Ax = xInner[elementsCountAroundSide][0]
        Ay = xInner[elementsCountAroundSide][1]
        originRC = (Ax*Ax + Ay*Ay - haustraRadius*haustraRadius) / (2*(-Ax - haustraRadius))
        RC = haustraRadius - originRC

        if originRC > -Ax:
            startTheta = math.asin(Ay/RC)
            thetaRC = (math.pi - startTheta)*2
        else:
            startTheta = math.pi - math.asin(Ay/RC)
            thetaRC = math.asin(Ay/RC)*2
        thetaPerElementAround = thetaRC/(elementsCountAround/3)

        for n in range(elementsCountAroundSide + 1):
            theta = startTheta + thetaPerElementAround * n
            x = [RC*math.cos(theta) - originRC,
                RC*math.sin(theta),
                0.0]
            xHaustraSide.append(x)

        ang = [-2/3*math.pi, 0.0, 2/3*math.pi]

        for i in range(3):
            rotAng = ang[i]
            cosRotAng = math.cos(rotAng)
            sinRotAng = math.sin(rotAng)
            for n in range(elementsCountAroundSide):
                theta = startTheta + thetaPerElementAround * n
                cosTheta = math.cos(theta)
                sinTheta = math.sin(theta)
                x = [ (RC*cosTheta - originRC)*cosRotAng - RC*sinTheta*sinRotAng,
                      (RC*cosTheta - originRC)*sinRotAng + RC*sinTheta*cosRotAng,
                      0.0]
                xHaustraInner.append(x)
                dx_ds1 = [(-RC*sinTheta*cosRotAng - RC*cosTheta*sinRotAng)*thetaPerElementAround,
                      (-RC*sinTheta*sinRotAng + RC*cosTheta*cosRotAng)*thetaPerElementAround,
                      0.0]
                d1InnerHaustraRaw.append(dx_ds1)
        d1InnerHaustra = smoothCubicHermiteDerivativesLoop(xHaustraInner, d1InnerHaustraRaw)

        # Sample arclength of haustra
        for n1 in range(elementsCountAround):
            if n1%(elementsCountAround/3) > 0.0:
                v1 = [xInner[n1][0], xInner[n1][1], 0.0]
                startArcLength = foldFactor * haustraLength
                d1 = [ c*startArcLength for c in unitZ]
                v2 = [xHaustraInner[n1][0], xHaustraInner[n1][1], haustraLength/2]
                midArcLength = haustraCurvatureFactor * haustraLength
                d2 = [ c*midArcLength for c in unitZ]
                v3 = [xInner[n1][0], xInner[n1][1], haustraLength]
                d3 = [ c*startArcLength for c in unitZ]
            else:
                v1 = [xInner[n1][0], xInner[n1][1], 0.0]
                v2 = [xInner[n1][0], xInner[n1][1], haustraLength/2]
                v3 = [xInner[n1][0], xInner[n1][1], haustraLength]
                d3 = d2 = d1 = [c* haustraLength/3 for c in unitZ]
            nx = [v1, v2, v3]
            nd1 = [d1, d2, d3]
            sx, sd1, se, sxi, _  = sampleCubicHermiteCurves(nx, nd1, elementsCountAlongHaustra)
            xInnerRaw.append(sx)
            dx_ds2InnerRaw.append(sd1)
            # if n1 == 28: # troubleshoot issue when samples suddenly go back to start
                # print('nx = ', nx)
                # print(sx)

        # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
        dx_ds1InnerList = dx_ds1InnerList + d1Inner
        for n2 in range(elementsCountAlongHaustra + 1):
            xAround = []
            unitdx_ds1Around = []
            for n1 in range(elementsCountAround):
                x = xInnerRaw[n1][n2]
                xInnerList.append(x)
                dx_ds2 = dx_ds2InnerRaw[n1][n2]
                dx_ds2InnerList.append(dx_ds2)

                unitTangent = normalise(dx_ds2)
                # Interhaustra
                if n2 == 0 or n2 > elementsCountAlongHaustra - 1:
                    dx_ds1 = d1Inner[n1]
                    unitdx_ds3 = crossproduct3(normalise(dx_ds1), unitTangent)
                else:
                    # Haustra
                    if elementsCountAlongHaustra == 2:
                        unitdx_ds1 = normalise(d1InnerHaustra[n1])
                    else:
                        if n1%(elementsCountAround/3) == 0: # intersection points
                            unitdx_ds1 = normalise(d1InnerHaustra[n1])
                        else: # points on clover
                            if elementsCountAlongHaustra > 3: 
                                if n2 < int(elementsCountAlongHaustra/2): # first half of haustraLength
                                    axisRot = crossproduct3(unitZ, unitTangent)
                                elif n2 > int(elementsCountAlongHaustra/2): # second half of haustraLength
                                    axisRot = crossproduct3(unitTangent, unitZ)
                            elif elementsCountAlongHaustra == 3: # 3 elementsAlongHaustra
                                axisRot = crossproduct3(unitTangent, unitZ)

                            rotFrame = rotationMatrixAboutAxis(axisRot, math.pi/2)
                            rotNormal = [rotFrame[j][0]*unitTangent[0] + rotFrame[j][1]*unitTangent[1] + rotFrame[j][2]*unitTangent[2] for j in range(3)]
                            unitdx_ds3 = normalise(rotNormal)
                            unitdx_ds1 = crossproduct3(unitTangent, unitdx_ds3)
                    xAround.append(x)
                    unitdx_ds1Around.append(unitdx_ds1)

            if n2 > 0 and n2 < elementsCountAlongHaustra:
                dx_ds1InnerAroundList = []
                for n1 in range(elementsCountAround):
                    v1 = xAround[n1]
                    d1 = unitdx_ds1Around[n1]
                    v2 = xAround[(n1+1)%elementsCountAround]
                    d2 = unitdx_ds1Around[(n1+1)%elementsCountAround]
                    arcLengthAround = computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                    dx_ds1 = [c*arcLengthAround for c in d1]
                    dx_ds1InnerAroundList.append(dx_ds1)
                d1Smoothed = smoothCubicHermiteDerivativesLoop(xAround, dx_ds1InnerAroundList)
                dx_ds1InnerList = dx_ds1InnerList + d1Smoothed

        dx_ds1InnerList = dx_ds1InnerList + d1Inner

        for n in range(len(xInnerList)):
            dx_ds3 = crossproduct3(normalise(dx_ds1InnerList[n]), normalise( dx_ds2InnerList[n]))
            unitdx_ds3 = normalise(dx_ds3)
            dx_ds3InnerUnitList.append(unitdx_ds3)

        # Pre-calculate node locations and derivatives on outer boundary
        xOuterList, d1Outer, curvatureInner = getOuterCoordinatesAndDerivativesFromInner(xInnerList, dx_ds1InnerList, dx_ds3InnerUnitList, wallThickness, elementsCountAlongHaustra, elementsCountAround)

        # Interpolate to get nodes through wall
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for nH in range(haustraCount):
                if nH == 0:
                    x, dx_ds1, dx_ds2, dx_ds3 = interpolatefromInnerAndOuter( xInnerList, xOuterList, wallThickness, xi3, curvatureInner, dx_ds1InnerList, dx_ds2InnerList, dx_ds3InnerUnitList, elementsCountAround, elementsCountAlongHaustra, elementsCountThroughWall)
                    xList = xList + x
                    dx_ds1List = dx_ds1List + dx_ds1
                    dx_ds2List = dx_ds2List + dx_ds2
                    dx_ds3List = dx_ds3List + dx_ds3
                else:
                    xNext = []
                    for n in range(elementsCountAround, len(x)):
                        newX = [x[n][0], x[n][1], x[n][2] + haustraLength*nH]
                        xNext.append(newX)
                    xList = xList + xNext
                    dx_ds1List = dx_ds1List + dx_ds1[elementsCountAround:]
                    dx_ds2List = dx_ds2List + dx_ds2[elementsCountAround:]
                    dx_ds3List = dx_ds3List + dx_ds3[elementsCountAround:]

        for n in range(len(xList)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3List[n])
            # print('NodeIdentifier = ', nodeIdentifier, xList[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1 # nextElementIdentifier
        now = (elementsCountAlong + 1)*elementsCountAround
        for e3 in range(elementsCountThroughWall):
            for e2 in range(elementsCountAlong):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni11 = e3*now + e2*elementsCountAround + e1 + 1
                    bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                    bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

        fm.endChange()

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
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)

def getOuterCoordinatesAndDerivativesFromInner(xInner, d1Inner, d3Inner, wallThickness, elementsCountAlongHaustra, elementsCountAround):
    """
    Generates coordinates and derivatives of outer surface from 
    coordinates and derivatives of inner surface using wall thickness
    and normals.
    param xInner: Coordinates on inner surface
    param d1Inner: Derivatives on inner surface around colon
    param d3Inner: Derivatives on inner surface through wall
    param wallThickness: Thickness of wall
    param elementsCountAlongHaustra: Number of elements along haustra
    param elementsCountAround: Number of elements around colon
    return xOuter: Coordinates on outer surface
    return nd1: Derivatives on outer surface around colon
    return curvatureInner: Curvature of coordinates on inner surface
    """
    xOuter = []
    nd1 = []
    dWall = []
    curvatureInner = []
    for n2 in range(elementsCountAlongHaustra + 1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            x = [xInner[n][i] + d3Inner[n][i]*wallThickness for i in range(3)]
            prevIdx = n-1 if (n1 != 0) else (n2+1)*elementsCountAround - 1
            nextIdx = n+1 if (n1 < elementsCountAround-1) else n2*elementsCountAround
            curvatureAround = 0.5*(
                getCubicHermiteCurvature(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], d3Inner[n], 1.0) +
                getCubicHermiteCurvature(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], d3Inner[n], 0.0))
            factor = 1.0 - curvatureAround*wallThickness
            nd1Outer = [ factor*c for c in d1Inner[n]]
            xOuter.append(x)
            nd1.append(nd1Outer)
            curvatureInner.append(curvatureAround)

    return xOuter, nd1, curvatureInner
    
def interpolatefromInnerAndOuter( xInner, xOuter, thickness, xi3, curvatureInner, d1Inner, d2Inner, d3InnerUnit, elementsCountAround, elementsCountAlongHaustra, elementsCountThroughWall):
    """
    Generate coordinates and derivatives at xi3 by interpolating with 
    inner and outer coordinates and derivatives.
    param xInner: Coordinates on inner surface
    param xOuter: Coordinates on outer surface
    param thickness: Thickness of wall
    param curvatureInner: Curvature of coordinates on inner surface
    param d1Inner: Derivatives on inner surface around colon
    param d2Inner: Derivatives on inner surface along haustra
    param d3InnerUnit: Unit derivatives on inner surface through wall
    param elementsCountAround: Number of elements around colon
    param elementsCountAlongHaustra: Number of elements along haustra
    param elementsCountThroughWall: Number of elements through wall
    return xList, dx_ds1List, dx_ds2List, dx_ds3List: Coordinates and derivatives on xi3
    """
    xList = []
    dx_ds1List = []
    dx_ds2List = []
    dx_ds3List =[]

    for n2 in range(elementsCountAlongHaustra+1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            # x
            innerx = xInner[n]
            outerx = xOuter[n]
            dWall = [thickness*c for c in d3InnerUnit[n]] 
            x = interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xList.append(x)
            # dx_ds1
            factor = 1.0 - curvatureInner[n]*thickness*xi3
            dx_ds1 = [ factor*c for c in d1Inner[n]]
            dx_ds1List.append(dx_ds1)
            # dx_ds2
            norm = d3InnerUnit[n]
            if n2 > 0 and n2 < elementsCountAlongHaustra:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
            elif n2 == 0:
                prevIdx = (elementsCountAlongHaustra-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
            else:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = elementsCountAround + n1
            curvatureAround = 0.5*(
                getCubicHermiteCurvature(xInner[prevIdx], d2Inner[prevIdx], xInner[n], d2Inner[n], norm, 1.0)+
                getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[nextIdx], d2Inner[nextIdx], norm, 0.0))
            factor = 1.0 - curvatureAround*thickness*xi3
            dx_ds2 = [ factor*c for c in d2Inner[n]]
            dx_ds2List.append(dx_ds2)
            #dx_ds3
            dx_ds3 = [c * thickness/elementsCountThroughWall for c in d3InnerUnit[n]]
            dx_ds3List.append(dx_ds3)

    return xList, dx_ds1List, dx_ds2List, dx_ds3List

def rotationMatrixAboutAxis(rotAxis, theta):
    """
    Generate the rotation matrix for rotation about an axis.
    :param rotAxis: axis of rotation
    :param theta: angle of rotation
    :return: rotation matrix
    """
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    C = 1 - cosTheta
    rotMatrix = ([[rotAxis[0]*rotAxis[0]*C + cosTheta, rotAxis[0]*rotAxis[1]*C - rotAxis[2]*sinTheta, rotAxis[0]*rotAxis[2]*C + rotAxis[1]*sinTheta],
        [rotAxis[1]*rotAxis[0]*C + rotAxis[2]*sinTheta, rotAxis[1]*rotAxis[1]*C + cosTheta, rotAxis[1]*rotAxis[2]*C - rotAxis[0]*sinTheta],
        [rotAxis[2]*rotAxis[0]*C - rotAxis[1]*sinTheta, rotAxis[2]*rotAxis[1]*C + rotAxis[0]*sinTheta, rotAxis[2]*rotAxis[2]*C + cosTheta]])
    return rotMatrix