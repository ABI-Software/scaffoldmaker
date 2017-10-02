"""
Generates a 3-D unit box mesh with variable numbers of elements in 3 directions,
and a hole through one axis. The mesh is designed to fit in with regular square
meshes on two sides, and tube meshes at the hole; derivatives and element xi
directions go around the hole.
"""

import math
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

def interpolateCubicHermite(v1, d1, v2, d2, xi):
    """
    Return cubic Hermite interpolated value of tuples v1, d1 (end 1) to v2, d2 (end 2) for xi in [0,1]
    :return: tuple containing result
    """
    xi2 = xi*xi
    xi3 = xi2*xi
    f1 = 1.0 - 3.0*xi2 + 2.0*xi3
    f2 = xi - 2.0*xi2 + xi3
    f3 = 3.0*xi2 - 2.0*xi3
    f4 = -xi2 + xi3
    result = []
    for i in range(len(v1)):
        result.append(f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i])
    return tuple(result)

def interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi):
    """
    Return cubic Hermite interpolated derivatives of tuples v1, d1 (end 1) to v2, d2 (end 2) for xi in [0,1]
    :return: tuple containing result
    """
    xi2 = xi*xi
    f1 = -6.0*xi + 6.0*xi2
    f2 = 1.0 - 4.0*xi + 3.0*xi2
    f3 = 6.0*xi - 6.0*xi2
    f4 = -2.0*xi + 3.0*xi2
    result = []
    for i in range(len(v1)):
        result.append(f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i])
    return tuple(result)

class MeshType_3d_boxhole1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Box Hole 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements 1' : 1,
            'Number of elements 2' : 1,
            'Number of elements 3' : 1,
            'Number of elements through wall' : 1,
            'Hole diameter' : 0.5,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Number of elements through wall',
            'Hole diameter',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements 1'] < 1) :
            options['Number of elements 1'] = 1
        if (options['Number of elements 2'] < 1) :
            options['Number of elements 2'] = 1
        if (options['Number of elements 3'] < 1) :
            options['Number of elements 3'] = 1
        if (options['Number of elements through wall'] < 1) :
            options['Number of elements through wall'] = 1
        if (options['Hole diameter'] < 0.0) :
            options['Hole diameter'] = 0.0
        elif (options['Hole diameter'] > 1.0) :
            options['Hole diameter'] = 1.0

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCount1 = options['Number of elements 1']
        elementsCount2 = options['Number of elements 2']
        elementsCount3 = options['Number of elements 3']
        elementsCountThroughWall = options['Number of elements through wall']
        elementsCountAround = 2*(elementsCount1 + elementsCount2)
        holeRadius = options['Hole diameter']*0.5
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = fm.createFieldFiniteElement(3)
        coordinates.setName('coordinates')
        coordinates.setManaged(True)
        coordinates.setTypeCoordinate(True)
        coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        coordinates.setComponentName(1, 'x')
        coordinates.setComponentName(2, 'y')
        coordinates.setComponentName(3, 'z')

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
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

        mesh = fm.findMeshByDimension(3)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eft.setNumberOfLocalScaleFactors(1)
        eft.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        # GRC: allow scale factor identifier for global -1 to be prescribed
        eft.setScaleFactorIdentifier(1, 1)
        if not useCrossDerivatives:
            for n in range(8):
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        # map dxi2 to -ds3, dxi3 -> ds2
        # xi2 derivatives use negative s3 derivative parameter:
        for n in range(8):
            ln = n + 1
            eft.setTermNodeParameter(n*8 + 3, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 3, 1, [1])
            eft.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS2, 1)
            if useCrossDerivatives:
                eft.setTermNodeParameter(n*8 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS3, 1)
                eft.setTermScaling(n*8 + 4, 1, [1])
                eft.setTermNodeParameter(n*8 + 6, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 1)
                eft.setTermNodeParameter(n*8 + 7, 1, ln, Node.VALUE_LABEL_D2_DS2DS3, 1)
                eft.setTermScaling(n*8 + 7, 1, [1])
                eft.setTermNodeParameter(n*8 + 8, 1, ln, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                eft.setTermScaling(n*8 + 8, 1, [1])

        # note I'm cheating here by using element-based scale factors
        # i.e. not shared by neighbouring elements, and I'm also using the same scale
        # factors for the bottom and top of each element
        eftOuter = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftOuter.setNumberOfLocalScaleFactors(9)
        eft.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        # GRC: allow scale factor identifier for global -1 to be prescribed
        eft.setScaleFactorIdentifier(1, 1)
        nonCrossDerivativesNodes = [0, 1, 4, 5] if useCrossDerivatives else range(8)
        for n in nonCrossDerivativesNodes:
            eftOuter.setFunctionNumberOfTerms(n*8 + 4, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 6, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 7, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 8, 0)
        # nodes 3,4,7,8: map dxi2 to -ds3, dxi3 -> ds2
        # xi2 derivatives use negative s3 derivative parameter:
        for n in [2, 3, 6, 7]:
            ln = n + 1
            eftOuter.setTermNodeParameter(n*8 + 3, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftOuter.setTermScaling(n*8 + 3, 1, [1])
            eftOuter.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS2, 1)
            if useCrossDerivatives:
                eftOuter.setTermNodeParameter(n*8 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS3, 1)
                eftOuter.setTermScaling(n*8 + 4, 1, [1])
                eftOuter.setTermNodeParameter(n*8 + 6, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 1)
                eftOuter.setTermNodeParameter(n*8 + 7, 1, ln, Node.VALUE_LABEL_D2_DS2DS3, 1)
                eftOuter.setTermScaling(n*8 + 7, 1, [1])
                eftOuter.setTermNodeParameter(n*8 + 8, 1, ln, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                eftOuter.setTermScaling(n*8 + 8, 1, [1])
        # nodes 1,2,5,6: general map dxi1, dsxi2 from ds1, ds2
        for n in [0, 1, 4, 5]:
            ln = n + 1
            sfo = (n % 2)*4 + 1
            eftOuter.setFunctionNumberOfTerms(n*8 + 2, 2)
            eftOuter.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftOuter.setTermScaling(n*8 + 2, 1, [1 + sfo])
            eftOuter.setTermNodeParameter(n*8 + 2, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftOuter.setTermScaling(n*8 + 2, 2, [2 + sfo])
            eftOuter.setFunctionNumberOfTerms(n*8 + 3, 2)
            eftOuter.setTermNodeParameter(n*8 + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftOuter.setTermScaling(n*8 + 3, 1, [3 + sfo])
            eftOuter.setTermNodeParameter(n*8 + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftOuter.setTermScaling(n*8 + 3, 2, [4 + sfo])

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateOuter = mesh.createElementtemplate()
        elementtemplateOuter.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuter.defineField(coordinates, -1, eftOuter)

        cache = fm.createFieldcache()

        # create nodes
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        wallThicknessMin = (0.5 - holeRadius)
        wallThicknessMax = math.sqrt(0.5) - holeRadius
        wallThicknessPerElement = wallThicknessMin/elementsCountThroughWall
        radius = holeRadius
        inner_x = []
        inner_d1 = []
        inner_d2 = []
        startRadians = math.pi*(-0.5 - elementsCount1/elementsCountAround)
        for n1 in range(elementsCountAround):
            radiansAround = startRadians + n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            inner_x.append((radius*cosRadiansAround, radius*sinRadiansAround))
            inner_d1.append((radiansPerElementAround*radius*-sinRadiansAround, radiansPerElementAround*radius*cosRadiansAround))
            #inner_d2.append((wallThicknessPerElement*cosRadiansAround, wallThicknessPerElement*sinRadiansAround))
            inner_d2.append((wallThicknessMin*cosRadiansAround, wallThicknessMin*sinRadiansAround))
        outer_x = []
        outer_d1 = []
        outer_d2 = []
        mag = math.sqrt(elementsCount1*elementsCount1 + elementsCount2*elementsCount2)
        outer_dx1_ds1 = 1.0 / elementsCount1
        outer_dx2_ds2 = 1.0 / elementsCount2
        cornerScale1 = 1.0 / max(elementsCount1 + 1, elementsCount2 + 1)
        for n in range(elementsCount1):
            x = -0.5 + n/elementsCount1
            outer_x.append((x, -0.5))
            mag = math.sqrt(x*x + 0.25)
            wallThicknessHere = mag - holeRadius
            scale2 = 2.0*wallThicknessHere - wallThicknessMin
            rx = x/mag
            ry = -0.5/mag
            r = 2.0*math.fabs(x)
            if n == 0:
                outer_d1.append((-ry*cornerScale1, rx*cornerScale1))
            else:
                outer_d1.append((outer_dx1_ds1, 0.0))
            outer_d2.append((rx*scale2, ry*scale2))
        for n in range(elementsCount2):
            y =  -0.5 + n/elementsCount2
            outer_x.append((0.5, y))
            mag = math.sqrt(0.25 + y*y)
            wallThicknessHere = mag - holeRadius
            scale2 = 2.0*wallThicknessHere - wallThicknessMin
            rx = 0.5/mag
            ry = y/mag
            r = 2.0*math.fabs(y)
            if n == 0:
                outer_d1.append((-ry*cornerScale1, rx*cornerScale1))
            else:
                outer_d1.append((0.0, outer_dx2_ds2))
            outer_d2.append((rx*scale2, ry*scale2))
        for n in range(elementsCount1):
            x = 0.5 - n/elementsCount1
            outer_x.append((x, 0.5))
            mag = math.sqrt(x*x + 0.25)
            wallThicknessHere = mag - holeRadius
            scale2 = 2.0*wallThicknessHere - wallThicknessMin
            rx = x/mag
            ry = 0.5/mag
            if n == 0:
                outer_d1.append((-ry*cornerScale1, rx*cornerScale1))
            else:
                outer_d1.append((-outer_dx1_ds1, 0.0))
            outer_d2.append((rx*scale2, ry*scale2))
        for n in range(elementsCount2):
            y =  0.5 - n/elementsCount2
            outer_x.append((-0.5, y))
            mag = math.sqrt(0.25 + y*y)
            wallThicknessHere = mag - holeRadius
            scale2 = 2.0*wallThicknessHere - wallThicknessMin
            rx = -0.5/mag
            ry = y/mag
            if n == 0:
                outer_d1.append((-ry*cornerScale1, rx*cornerScale1))
            else:
                outer_d1.append((0.0, -outer_dx2_ds2))
            outer_d2.append((rx*scale2, ry*scale2))

        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCount3 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        outer_dx_ds1 = [ 1.0 / elementsCount1, 0.0, 0.0 ]
        outer_dx_ds2 = [ 0.0, 1.0 / elementsCount2, 0.0 ]
        outer_dx_ds3 = [ 0.0, 0.0, 1.0 / elementsCount3 ]
        zero = [ 0.0, 0.0, 0.0 ]
        for n3 in range(elementsCount3 + 1):
            x[2] = n3 / elementsCount3
            # outer nodes
            for n1 in range(elementsCountAround):
                x[0] = outer_x[n1][0]
                x[1] = outer_x[n1][1]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, outer_dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, outer_dx_ds2)
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ outer_d1[n1][0], outer_d1[n1][1], 0.0 ])
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ outer_d2[n1][0], outer_d2[n1][1], 0.0 ])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, outer_dx_ds3)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier = nodeIdentifier + 1
            # inner nodes
            for n2 in range(elementsCountThroughWall):
                xir = (n2 + 1)/elementsCountThroughWall
                xi = 1.0 - xir
                for n1 in range(elementsCountAround):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    v = interpolateCubicHermite(inner_x[n1], inner_d2[n1], outer_x[n1], outer_d2[n1], xi)
                    x[0] = v[0]
                    x[1] = v[1]
                    dx_ds1[0] = xir*inner_d1[n1][0] + xi*outer_d1[n1][0]
                    dx_ds1[1] = xir*inner_d1[n1][1] + xi*outer_d1[n1][1]
                    #dx_ds3[0] = xir*inner_d2[n1][0] + xi*outer_d2[n1][0]
                    #dx_ds3[1] = xir*inner_d2[n1][1] + xi*outer_d2[n1][1]
                    d2 = interpolateCubicHermiteDerivative(inner_x[n1], inner_d2[n1], outer_x[n1], outer_d2[n1], xi)
                    dx_ds3[0] = d2[0]/elementsCountThroughWall  # *wallThicknessPerElement
                    dx_ds3[1] = d2[1]/elementsCountThroughWall  # *wallThicknessPerElement
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        no3 = (elementsCountThroughWall + 1)*elementsCountAround
        for e3 in range(elementsCount3):
            # first row general maps ds1, ds2 to dxi1, dxi2
            for e1 in range(elementsCountAround):
                en = (e1 + 1)%elementsCountAround
                element = mesh.createElement(elementIdentifier, elementtemplateOuter)
                bni11 = e3*no3 + e1 + 1
                bni12 = e3*no3 + en + 1
                bni21 = e3*no3 + elementsCountAround + e1 + 1
                bni22 = e3*no3 + elementsCountAround + en + 1
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + no3, bni12 + no3, bni21 + no3, bni22 + no3 ]
                result = element.setNodesByIdentifier(eftOuter, nodeIdentifiers)
                onX = (e1 % (elementsCount1 + elementsCount2)) < elementsCount1
                rev = e1 >= (elementsCount1 + elementsCount2)
                one = -1.0 if rev else 1.0
                vx = one if onX else 0.0
                vy = 0.0 if onX else one
                scaleFactors = [
                    -1.0,
                    vx, vy,
                    -outer_d2[e1][0]/outer_dx_ds1[0]/elementsCountThroughWall, -outer_d2[e1][1]/outer_dx_ds2[1]/elementsCountThroughWall,
                    vx, vy,
                    -outer_d2[en][0]/outer_dx_ds1[0]/elementsCountThroughWall, -outer_d2[en][1]/outer_dx_ds2[1]/elementsCountThroughWall ]
                element.setScaleFactors(eftOuter, scaleFactors)
                elementIdentifier = elementIdentifier + 1

            # remaining rows
            for e2 in range(1, elementsCountThroughWall):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni11 = e3*no3 + e2*elementsCountAround + e1 + 1
                    bni12 = e3*no3 + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    bni21 = e3*no3 + (e2 + 1)*elementsCountAround + e1 + 1
                    bni22 = e3*no3 + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + no3, bni12 + no3, bni21 + no3, bni22 + no3 ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    element.setScaleFactors(eft, [-1.0])
                    elementIdentifier = elementIdentifier + 1

        fm.endChange()

