'''
Utility functions for generating tubular mesh.
Created on Oct 9, 2018

@author: Mabelle Lin
'''
from __future__ import division
import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.geometry import *
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

def generatetubemesh(region,
    elementsCountAlong,
    elementsCountAround,
    elementsCountThroughWall,
    a, b, wallThickness,
    useCrossDerivatives,
    useCubicHermiteThroughWall, # or Zinc Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE etc.
    nextNodeIdentifier = 1, nextElementIdentifier = 1
    ):
    '''
    Generates a 3-D tubular mesh with variable numbers of elements
    around, along the central axis, and radially through wall.
    The tubular mesh follows the profile of a central line and 
    has a elliptical cross-section.
    :param elementsCountAlong: number of elements along tube
    :param elementsCountAround: number of elements around tube
    :param elementsCountThroughWall: number of elements through wall thickness
    :param a: major axis length of inner ellipse
    :param b: minor axis length of inner ellipse
    :param wallThickness: Thickness of wall 
    :param useCubicHermiteThroughWall: use linear when false
    :param nextNodeIdentifier, nextElementIdentifier: Next identifiers to use and increment.
    :return: Final values of nextNodeIdentifier, nextElementIdentifier
    '''

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
    nodeIdentifier = nextNodeIdentifier
    wallThicknessPerElement = wallThickness/elementsCountThroughWall
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    zero = [ 0.0, 0.0, 0.0 ]

    for n3 in range(elementsCountThroughWall + 1):
        aThroughWallElement = a + wallThickness*(n3/elementsCountThroughWall)
        bThroughWallElement = b + wallThickness*(n3/elementsCountThroughWall)
        perimeterAroundWallElement = getApproximateEllipsePerimeter(aThroughWallElement, bThroughWallElement)
        arcLengthPerElementAround = perimeterAroundWallElement / elementsCountAround
        for n2 in range(elementsCountAlong + 1):
            x[2] = n2 / elementsCountAlong
            prevRadiansAround = -1*updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0, arcLengthPerElementAround)
            for n1 in range(elementsCountAround):
                arcLengthAround = n1*arcLengthPerElementAround
                radiansAround = updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0, arcLengthAround)
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = aThroughWallElement*cosRadiansAround
                x[1] = bThroughWallElement*sinRadiansAround
                dx_ds1[0] = (radiansAround - prevRadiansAround)*aThroughWallElement*-sinRadiansAround
                dx_ds1[1] = (radiansAround - prevRadiansAround)*bThroughWallElement*cosRadiansAround
                dx_ds3[0] = wallThicknessPerElement*cosRadiansAround
                dx_ds3[1] = wallThicknessPerElement*sinRadiansAround
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if useCubicHermiteThroughWall:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier = nodeIdentifier + 1
                prevRadiansAround = radiansAround

    # create elements
    elementIdentifier = nextElementIdentifier
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
    
    return nodeIdentifier, elementIdentifier