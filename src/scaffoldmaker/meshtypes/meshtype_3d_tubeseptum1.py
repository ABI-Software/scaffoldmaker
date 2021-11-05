"""
Generates a 3-D 'tube septum' mesh crossing from one side of a tube to the
other, including elements for the outside, suitable for joining one half of a
tube on one side, and the other half on the other.
It is the middle line in (|).
The number of elements along the tube and across the septum can be varied.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite


class MeshType_3d_tubeseptum1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Tube Septum 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 1,
            'Number of elements across' : 2,
            'Wall thickness left' : 0.25,
            'Wall thickness right' : 0.25,
            'Flange length' : 0.125,
            'Bulge radius' : 0.0,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements across',
            'Wall thickness left',
            'Wall thickness right',
            'Flange length',
            'Bulge radius',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements along'] < 1) :
            options['Number of elements along'] = 1
        if (options['Number of elements across'] < 2) :
            options['Number of elements across'] = 2
        if (options['Wall thickness left'] < 0.0) :
            options['Wall thickness left'] = 0.0
        elif (options['Wall thickness left'] > 0.5) :
            options['Wall thickness left'] = 0.5
        if (options['Wall thickness right'] < 0.0) :
            options['Wall thickness right'] = 0.0
        elif (options['Wall thickness right'] > 0.5) :
            options['Wall thickness right'] = 0.5
        if (options['Flange length'] < 0.0) :
            options['Flange length'] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountAlong = options['Number of elements along']
        elementsCountAcross = options['Number of elements across']
        wallThickness = [ options['Wall thickness left'], options['Wall thickness right'] ]
        flangeLength = options['Flange length']
        bulgeRadius = options['Bulge radius']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

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

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()
        eftOuter = tricubichermite.createEftTubeSeptumOuter()
        eftInner1 = tricubichermite.createEftTubeSeptumInner1()
        eftInner2 = tricubichermite.createEftTubeSeptumInner2()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)
        elementtemplateOuter = mesh.createElementtemplate()
        elementtemplateOuter.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuter.defineField(coordinates, -1, eftOuter)
        elementtemplateInner1 = mesh.createElementtemplate()
        elementtemplateInner1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInner1.defineField(coordinates, -1, eftInner1)
        elementtemplateInner2 = mesh.createElementtemplate()
        elementtemplateInner2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInner2.defineField(coordinates, -1, eftInner2)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        #radiansPerElementAcross = math.pi/elementsCountAcross
        bevelAngle = math.pi/2
        sinBevelAngle = math.sin(bevelAngle)
        cosBevelAngle = math.cos(bevelAngle)

        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        wallThicknessSeptum = wallThickness[0]
        #wallThicknessMin = min(wallThickness)
        wallThicknessMax = max(wallThickness)
        for n3 in range(2):
            sign = -1.0 if (n3 == 0) else 1.0
            radiusY = 0.5*wallThicknessSeptum
            radiusX = 0.5 - wallThicknessMax
            for n2 in range(elementsCountAlong + 1):
                x[2] = n2 / elementsCountAlong
                for n1 in range(elementsCountAcross + 3):
                    if (n1 == 0) or (n1 == (elementsCountAcross + 2)):
                        flip = -1.0 if (n1 == 0) else 1.0
                        x[0] = radiusX + wallThickness[n3]
                        if n1 == 0:
                            x[0] = -x[0]
                        x[1] = -sign*(radiusY + flangeLength)
                        #if wallThickness[0] > wallThickness[1]:
                        #    x[1] -= flangeLength
                        #elif wallThickness[0] < wallThickness[1]:
                        #    x[1] += flangeLength
                        dx_ds1[0] = 0.0
                        dx_ds1[1] = 2.0*flip*(radiusY + flangeLength)
                        dx_ds3[0] = flip*wallThickness[n3]
                        dx_ds3[1] = 0.0
                    elif (n1 == 1) or (n1 == (elementsCountAcross + 1)):
                        flip = -1.0 if (n1 == 1) else 1.0
                        radius = radiusX
                        x[0] = flip*radius
                        x[1] = -sign*(radiusY + flangeLength)
                        #dx_ds1[0] = -sign*2.0*radiusX/elementsCountAcross*cosBevelAngle
                        #dx_ds1[1] = 2.0*radiusX/elementsCountAcross*sinBevelAngle
                        mag1x = radiusX/elementsCountAcross
                        mag1y = flangeLength
                        #mag1min = radiusX/elementsCountAcross + flangeLength + math.sqrt(mag1x*mag1x + mag1y*mag1y)
                        mag1min = 2.0*math.sqrt(mag1x*mag1x + mag1y*mag1y)
                        mag1max = 2.0*(radiusY + flangeLength)
                        mag1 = mag1min if (mag1min < mag1max) else mag1max
                        dx_ds1[0] = -sign*mag1*cosBevelAngle
                        dx_ds1[1] = mag1*sinBevelAngle
                        #if ((n3 == 1) and (wallThickness[0] > wallThickness[1])) or \
                        #    ((n3 == 0) and (wallThickness[0] < wallThickness[1])):
                        dx_ds3[0] = wallThickness[n3]
                        dx_ds3[1] = 0.0
                        #else:
                        #    dx_ds3[0] = wallThickness[n3]*sinBevelAngle
                        #    dx_ds3[1] = sign*wallThickness[n3]*cosBevelAngle
                        if n1 == 1:
                            dx_ds1[1] = -dx_ds1[1]
                            dx_ds3[0] = -dx_ds3[0]
                    else:
                        f1 = (n1 - 1)/elementsCountAcross
                        f2 = 1.0 - f1
                        x[0] = (f1 - f2)*radiusX
                        x[1] = -sign*radiusY
                        dx_ds1[0] = 2.0*radiusX/elementsCountAcross
                        dx_ds1[1] = 0.0
                        dx_ds3[0] = 0.0
                        dx_ds3[1] = -wallThicknessSeptum
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

        rimCrossAngle = math.pi/4
        sinRimCrossAngle = math.sin(rimCrossAngle)
        cosRimCrossAngle = math.cos(rimCrossAngle)
        scaleFactorsOuter = [
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsInner1 = [ -1.0, 4.0,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, sinRimCrossAngle,
            cosRimCrossAngle, -sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsInner2 = [ -1.0, 4.0,
            cosRimCrossAngle, -sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, sinRimCrossAngle
        ]

        # create elements
        elementIdentifier = 1
        rno = elementsCountAcross + 3
        wno = (elementsCountAlong + 1)*rno
        for e2 in range(elementsCountAlong):
            bn = e2*rno
            element = mesh.createElement(elementIdentifier, elementtemplateOuter)
            bni11 = bn + 2
            bni12 = bn + wno + 2
            bni21 = bn + rno + 2
            bni22 = bn + wno + rno + 2
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 - 1, bni12 - 1, bni21 - 1, bni22 - 1 ]
            result = element.setNodesByIdentifier(eftOuter, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, nodeIdentifiers)
            element.setScaleFactors(eftOuter, scaleFactorsOuter)
            elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateInner1)
            bni11 = bn + 2
            bni12 = bn + 3
            bni21 = bn + rno + 2
            bni22 = bn + rno + 3
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
            result = element.setNodesByIdentifier(eftInner1, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
            element.setScaleFactors(eftInner1, scaleFactorsInner1)
            #print(result, 'element', elementIdentifier, 'scale factors', scaleFactorsInner1)
            elementIdentifier = elementIdentifier + 1

            for e1 in range(elementsCountAcross - 2):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = bn + e1 + 3
                bni12 = bn + e1 + 4
                bni21 = bn + rno + e1 + 3
                bni22 = bn + rno + e1 + 4
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                #print(result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateInner2)
            bni11 = bn + rno - 2
            bni12 = bn + rno - 1
            bni21 = bn + 2*rno -2
            bni22 = bn + 2*rno - 1
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
            result = element.setNodesByIdentifier(eftInner2, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
            element.setScaleFactors(eftInner2, scaleFactorsInner2)
            #print(result, 'element', elementIdentifier, 'scale factors', scaleFactorsInner1)
            elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateOuter)
            bni11 = bn + wno + rno - 1
            bni12 = bn + rno - 1
            bni21 = bn + wno + 2*rno - 1
            bni22 = bn + 2*rno - 1
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + 1, bni12 + 1, bni21 + 1, bni22 + 1 ]
            result = element.setNodesByIdentifier(eftOuter, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, nodeIdentifiers)
            element.setScaleFactors(eftOuter, scaleFactorsOuter)
            elementIdentifier = elementIdentifier + 1

        if bulgeRadius != 0.0:
            # cylindrical polar coordinates:
            # r = y - bulgeRadius
            # theta = -x / bulgeRadius
            # z = z
            yxzCoordinates = fm.createFieldComponent(coordinates, [2, 1, 3])
            scale = fm.createFieldConstant([1.0, -1.0/bulgeRadius, 1.0 ])
            scaleCoordinates = fm.createFieldMultiply(yxzCoordinates, scale)
            offset = fm.createFieldConstant([-bulgeRadius, 0.0, 0.0 ])
            polarCoordinates = fm.createFieldAdd(scaleCoordinates, offset)
            polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
            rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
            rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
            newyxzCoordinates = fm.createFieldSubtract(rcCoordinates, offset)
            newCoordinates = fm.createFieldComponent(newyxzCoordinates, [2, 1, 3])
            fieldassignment = coordinates.createFieldassignment(newCoordinates)
            result = fieldassignment.assign()

        fm.endChange()
        return []
