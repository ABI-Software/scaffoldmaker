"""
Generates a 3-D unit solid sphere mesh with variable numbers of elements
around, up the central axis, and radially.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_solidsphere1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Solid Sphere 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 4,
            'Number of elements up' : 4,
            'Number of elements radial' : 1,
            'Diameter' : 1.0,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements up' : 1,
            'Refine number of elements radial' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements up',
            'Number of elements radial',
            'Diameter',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements radial'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements radial',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements radial']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 2:
            options['Number of elements around'] = 2
        if options['Diameter'] < 0.0:
            options['Diameter'] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountRadial = options['Number of elements radial']
        useCrossDerivatives = options['Use cross derivatives']
        radius = 0.5*options['Diameter']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        else:
            nodetemplate = nodetemplateApex

        cache = fm.createFieldcache()

        #################
        # Create nodes
        #################

        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementUp = math.pi/elementsCountUp

        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]

        cubicArcLengthList = [0.0]*(elementsCountUp+1)

        # Pre-calculate cubicArcLength along elementsCountUp
        for n2 in range(1,elementsCountUp + 1):
            radiansUp = n2*radiansPerElementUp
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)

            # Calculate cubic hermite arclength linking point on axis to surface on sphere
            v1 = [0.0, 0.0, -radius+n2*2.0*radius/elementsCountUp]
            d1 = [0.0, 1.0, 0.0]
            v2 = [
                 radius*math.cos(math.pi/2.0)*sinRadiansUp,
                 radius*math.sin(math.pi/2.0)*sinRadiansUp,
                 -radius*cosRadiansUp
            ]
            d2 = [math.cos(math.pi/2.0)*sinRadiansUp,math.sin(math.pi/2.0)*sinRadiansUp,-cosRadiansUp]
            cubicArcLengthList[n2] = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)

        # Create node for bottom pole
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -radius ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ 0.0, radius*radiansPerElementUp, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, -radius*2.0/elementsCountUp])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

        # Create nodes along axis between top and bottom poles
        for n2 in range(1,elementsCountUp):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -radius+n2*2.0*radius/elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ cubicArcLengthList[n2]/elementsCountRadial, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ 0.0, 0.0, radius*2.0/elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, cubicArcLengthList[n2]/elementsCountRadial, 0.0 ])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        # Create nodes for top pole 
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, radius])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ 0.0, radius*radiansPerElementUp, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, radius*2.0/elementsCountUp ])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

        # Create other nodes
        for n3 in range(1,elementsCountRadial+1):
            xi = 1/elementsCountRadial*n3
            radiansUpArcOriginList = [0.0]*(elementsCountUp)

            # Pre-calculate RC for points on vertical arc running between top and bottom poles
            pt = [0.0, radius*xi, 0.0]
            arcOrigin = (radius*radius - pt[2]*pt[2] - pt[1]*pt[1])/(-2.0*pt[1])
            RC = math.sqrt(arcOrigin*arcOrigin + radius*radius)

            radiansUpArcOriginList[0] = math.acos(-radius/RC)

            # Identify nodes on the vertical arc using radiansAround = pi/2
            for n2 in range(1,elementsCountUp):
                radiansUp = n2*radiansPerElementUp
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)

                # Calculate node coordinates on arc using cubic hermite interpolation
                cubicArcLength = cubicArcLengthList[n2]
                v1 = [0.0, 0.0, -radius+n2*2.0*radius/elementsCountUp]
                d1 = [math.cos(math.pi/2.0), math.sin(math.pi/2.0), 0.0]
                d1 = vector.normalise(d1)
                d1 = [d*cubicArcLength for d in d1]
                v2 = [
                     radius*math.cos(math.pi/2.0)*sinRadiansUp,
                     radius*math.sin(math.pi/2.0)*sinRadiansUp,
                     -radius*cosRadiansUp
                    ]     
                d2 = [math.cos(math.pi/2.0)*sinRadiansUp,math.sin(math.pi/2.0)*sinRadiansUp,-cosRadiansUp]
                d2 = vector.normalise(d2)
                d2 = [d*cubicArcLength for d in d2]
                x = interp.interpolateCubicHermite(v1, d1, v2, d2, xi)

                # Calculate radiansUp for each point wrt arcOrigin
                radiansUpArcOriginList[n2] = math.acos(x[2]/RC)

            for n2 in range(1,elementsCountUp):
                radiansUp = n2*radiansPerElementUp
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)

                for n1 in range(elementsCountAround):
                    radiansAround = n1*radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    cubicArcLength = cubicArcLengthList[n2]

                    # Calculate node coordinates on arc using cubic hermite interpolation
                    v1 = [0.0, 0.0, -radius+n2*2.0*radius/elementsCountUp]
                    d1 = [cosRadiansAround, sinRadiansAround, 0.0]
                    d1 = vector.normalise(d1)
                    d1 = [d*cubicArcLength for d in d1]
                    v2 = [
                         radius*cosRadiansAround*sinRadiansUp,
                         radius*sinRadiansAround*sinRadiansUp,
                        -radius*cosRadiansUp
                    ]
                    d2 = [cosRadiansAround*sinRadiansUp,sinRadiansAround*sinRadiansUp,-cosRadiansUp]
                    d2 = vector.normalise(d2)
                    d2 = [d*cubicArcLength for d in d2]
                    x = interp.interpolateCubicHermite(v1, d1, v2, d2, xi)

                    # For dx_ds1 - Calculate radius wrt origin where interpolated points lie on
                    orthoRadius = vector.magnitude(x)
                    orthoRadiansUp = math.pi - math.acos(x[2]/orthoRadius)
                    sinOrthoRadiansUp = math.sin(orthoRadiansUp)
                    cosOrthoRadiansUp = math.cos(orthoRadiansUp)

                    # For dx_ds2 - Assign radiansUp from radiansUpArcOriginList and calculate diff between radiansUp as we move up
                    radiansUpArcOrigin =  radiansUpArcOriginList[n2]
                    sinRadiansUpArcOrigin = math.sin(radiansUpArcOrigin)
                    cosRadiansUpArcOrigin = math.cos(radiansUpArcOrigin)
                    radiansPerElementUpArcOrigin = radiansUpArcOriginList[n2]-radiansUpArcOriginList[n2-1]

                    dx_ds1 = [
                        orthoRadius*-sinRadiansAround*sinOrthoRadiansUp*radiansPerElementAround,
                        orthoRadius*cosRadiansAround*sinOrthoRadiansUp*radiansPerElementAround,
                        0.0
                    ]

                    dx_ds2 = [
                        RC*cosRadiansAround*cosRadiansUpArcOrigin*radiansPerElementUpArcOrigin,
                        RC*sinRadiansAround*cosRadiansUpArcOrigin*radiansPerElementUpArcOrigin,
                        -RC*sinRadiansUpArcOrigin*radiansPerElementUpArcOrigin
                    ]

                    dx_ds3 = interp.interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
                    dx_ds3 = vector.normalise(dx_ds3)
                    dx_ds3 = [d*cubicArcLength/elementsCountRadial for d in dx_ds3]

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
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

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        # Regular elements
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        # Bottom tetrahedon elements
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Axial elements
        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Top tetrahedron elements
        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Bottom pyramid elements
        elementtemplate4 = mesh.createElementtemplate()
        elementtemplate4.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Top pyramid elements
        elementtemplate5 = mesh.createElementtemplate()
        elementtemplate5.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementIdentifier = 1

        no2 = elementsCountAround
        no3 = elementsCountAround*(elementsCountUp - 1)
        rni = (1 + elementsCountUp) - no3 - no2 + 1 # regular node identifier
        radiansPerElementAround = 2.0*math.pi/elementsCountAround

        for e3 in range(elementsCountRadial):

            # Create elements on bottom pole
            radiansIncline = math.pi*0.5*e3/elementsCountRadial
            radiansInclineNext = math.pi*0.5*(e3 + 1)/elementsCountRadial
            
            if e3 == 0:
                # Create tetrahedron elements on the bottom pole
                bni1 = elementsCountUp + 2

                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft1 = tricubichermite.createEftTetrahedronBottom(va*100, vb*100, 10000)
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    nodeIdentifiers = [ 1, 2, bni1 + va, bni1 + vb ]
                    result1 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = va*radiansPerElementAround
                    radiansAroundNext = vb*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansIncline), math.sin(radiansIncline),
                        math.cos(radiansInclineNext), math.sin(radiansInclineNext)
                    ]
                    result2 = element.setScaleFactors(eft1, scalefactors)
                    # print('Tetrahedron Bottom element', elementIdentifier, result1, result2, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            else:
                # Create pyramid elements on the bottom pole
                bni4 = elementsCountUp + 1 + (e3-1)*no3 + 1

                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft4 = tricubichermite.createEftPyramidBottom(va*100, vb*100, 100000 + e3*2)
                    elementtemplate4.defineField(coordinates, -1, eft4)
                    element = mesh.createElement(elementIdentifier, elementtemplate4)
                    nodeIdentifiers = [ 1, bni4 + va, bni4 + vb, bni4 + no3 + va, bni4 + no3 + vb ]
                    result1 = element.setNodesByIdentifier(eft4, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = va*radiansPerElementAround
                    radiansAroundNext = vb*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansIncline), math.sin(radiansIncline),
                        math.cos(radiansInclineNext), math.sin(radiansInclineNext)
                    ]
                    result2 = element.setScaleFactors(eft4, scalefactors)
                    # print('pyramid bottom element', elementIdentifier, result1, result2, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            # create regular radial elements
            for e2 in range(1, elementsCountUp-1):
                if e3 == 0:
                    for e1 in range(elementsCountAround):
                        # create central radial elements: 6 node wedges
                        va = e1
                        vb = (e1 + 1)%elementsCountAround
                        eft2 = tricubichermite.createEftWedgeRadial(va*100, vb*100)
                        elementtemplate2.defineField(coordinates, -1, eft2)
                        element = mesh.createElement(elementIdentifier, elementtemplate2)
                        bni2 = elementsCountUp + 1 + (e2-1) * no2 + 1
                        nodeIdentifiers = [ e3 + e2 + 1, e3 + e2 + 2, bni2 + va, bni2 + vb, bni2 + va + elementsCountAround, bni2 + vb + elementsCountAround ]
                        result1 = element.setNodesByIdentifier(eft2, nodeIdentifiers)
                        # set general linear map coefficients
                        radiansAround = va*radiansPerElementAround
                        radiansAroundNext = vb*radiansPerElementAround
                        scalefactors = [
                            -1.0,
                            math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                            math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                            math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                            math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround
                        ]
                        result2 = element.setScaleFactors(eft2, scalefactors)
                        # print('axis element', elementIdentifier, result1, result2, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
                else:
                    # Regular elements 
                    bni = rni + e3*no3 + e2*no2

                    for e1 in range(elementsCountAround):
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        na = e1
                        nb = (e1 + 1)%elementsCountAround
                        nodeIdentifiers = [
                            bni       + na, bni       + nb, bni       + no2 + na, bni       + no2 + nb,
                            bni + no3 + na, bni + no3 + nb, bni + no3 + no2 + na, bni + no3 + no2 + nb
                        ]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        # print('regular element', elementIdentifier, result, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1

            # Create elements on top pole
            radiansIncline = math.pi*0.5*e3/elementsCountRadial
            radiansInclineNext = math.pi*0.5*(e3 + 1)/elementsCountRadial

            if e3 == 0:
                # # Create tetrahedron elements on the top pole
                bni3 = elementsCountUp + 1 + (elementsCountUp-2) * no2 + 1

                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft3 = tricubichermite.createEftTetrahedronTop(va*100, vb*100, 100000)
                    elementtemplate3.defineField(coordinates, -1, eft3)
                    element = mesh.createElement(elementIdentifier, elementtemplate3)
                    nodeIdentifiers = [ elementsCountUp, elementsCountUp + 1, bni3 + va, bni3 + vb ]
                    result1 = element.setNodesByIdentifier(eft3, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = va*radiansPerElementAround
                    radiansAroundNext = vb*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansIncline), math.sin(radiansIncline),
                        math.cos(radiansInclineNext), math.sin(radiansInclineNext)
                    ]
                    result2 = element.setScaleFactors(eft3, scalefactors)
                    # print('Tetrahedron top element', elementIdentifier, result1, result2, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            else:
                # Create pyramid elements on the top pole
                bni5 = elementsCountUp + 1 + (e3-1)*no3 + (elementsCountUp-2) * no2 + 1

                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft5 = tricubichermite.createEftPyramidTop(va*100, vb*100, 100000 + e3*2)

                    elementtemplate5.defineField(coordinates, -1, eft5)
                    element = mesh.createElement(elementIdentifier, elementtemplate5)
                    nodeIdentifiers = [ bni5 + va, bni5 + vb, elementsCountUp + 1, bni5 + no3 + va, bni5 + no3 + vb ]
                    result1 = element.setNodesByIdentifier(eft5, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = va*radiansPerElementAround
                    radiansAroundNext = vb*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        math.cos(radiansAround), math.sin(radiansAround), radiansPerElementAround,
                        math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAround,
                        math.cos(radiansIncline), math.sin(radiansIncline),
                        math.cos(radiansInclineNext), math.sin(radiansInclineNext)
                    ]
                    result2 = element.setScaleFactors(eft5, scalefactors)
                    # print('pyramid top element', elementIdentifier, result1, result2, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountRadial = options['Refine number of elements radial']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountRadial)
