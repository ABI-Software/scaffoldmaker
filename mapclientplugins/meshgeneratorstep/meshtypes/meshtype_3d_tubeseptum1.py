"""
Generates a 3-D 'tube septum' mesh crossing from one side of a tube to the
other, including elements for the outside, suitable for joining one half of a
tube on one side, and the other half on the other.
It is the middle line in (|).
The number of elements along the tube and across the septum can be varied.
"""

import math
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_tubeseptum1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Tube Septum 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements along' : 1,
            'Number of elements across' : 2,
            'Wall thickness left' : 0.25,
            'Wall thickness right' : 0.25,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements across',
            'Wall thickness left',
            'Wall thickness right',
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

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAlong = options['Number of elements along']
        elementsCountAcross = options['Number of elements across']
        wallThickness = [ options['Wall thickness left'], options['Wall thickness right'] ]
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
        if not useCrossDerivatives:
            for n in range(8):
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        eftOuter = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        # general linear map at 4 nodes for one derivative
        eftOuter.setNumberOfLocalScaleFactors(8)
        for s in range(8):
            eftOuter.setScaleFactorType(s + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftOuter.setScaleFactorIdentifier(s + 1, (s % 2) + 1)
        if useCrossDerivatives:
            noCrossRange = range(4)
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eftOuter.setFunctionNumberOfTerms(n*8 + 4, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 6, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 7, 0)
            eftOuter.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi1 from ds1 + ds3 for first 4 nodes
        for n in range(4):
            ln = n + 1
            eftOuter.setFunctionNumberOfTerms(n*8 + 2, 2)
            eftOuter.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftOuter.setTermScaling(n*8 + 2, 1, [n*2 + 1])
            eftOuter.setTermNodeParameter(n*8 + 2, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftOuter.setTermScaling(n*8 + 2, 2, [n*2 + 2])

        eftInner1 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        # negate dxi1 plus general linear map at 4 nodes for one derivative
        eftInner1.setNumberOfLocalScaleFactors(9)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftInner1.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftInner1.setScaleFactorIdentifier(1, 1)
        for s in range(8):
            eftInner1.setScaleFactorType(s + 2, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftInner1.setScaleFactorIdentifier(s + 2, (s % 2) + 1)
        if useCrossDerivatives:
            noCrossRange = [ 0, 2, 4, 6 ]
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eftInner1.setFunctionNumberOfTerms(n*8 + 4, 0)
            eftInner1.setFunctionNumberOfTerms(n*8 + 6, 0)
            eftInner1.setFunctionNumberOfTerms(n*8 + 7, 0)
            eftInner1.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi3 from ds1 + ds3 for 4 odd nodes
        s = 0
        for n in [ 0, 2, 4, 6 ]:
            ln = n + 1
            eftInner1.setFunctionNumberOfTerms(n*8 + 5, 2)
            eftInner1.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftInner1.setTermScaling(n*8 + 5, 1, [s*2 + 2])
            eftInner1.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftInner1.setTermScaling(n*8 + 5, 2, [s*2 + 3])
            s += 1
        # negate d/dxi1 at 2 nodes
        for n in [4, 6]:
            result = eftInner1.setTermScaling(n*8 + 2, 1, [1])

        eftInner2 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        # negate dxi1 plus general linear map at 4 nodes for one derivative
        eftInner2.setNumberOfLocalScaleFactors(9)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftInner2.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftInner2.setScaleFactorIdentifier(1, 1)
        for s in range(8):
            eftInner2.setScaleFactorType(s + 2, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftInner2.setScaleFactorIdentifier(s + 2, (s % 2) + 1)
        if useCrossDerivatives:
            noCrossRange = [ 1, 3, 5, 7 ]
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eftInner2.setFunctionNumberOfTerms(n*8 + 4, 0)
            eftInner2.setFunctionNumberOfTerms(n*8 + 6, 0)
            eftInner2.setFunctionNumberOfTerms(n*8 + 7, 0)
            eftInner2.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi3 from ds1 + ds3 for 4 even nodes
        s = 0
        for n in [ 1, 3, 5, 7 ]:
            ln = n + 1
            eftInner2.setFunctionNumberOfTerms(n*8 + 5, 2)
            eftInner2.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftInner2.setTermScaling(n*8 + 5, 1, [1, s*2 + 2])
            eftInner2.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftInner2.setTermScaling(n*8 + 5, 2, [1, s*2 + 3])
            s += 1
        # negate d/dxi1 at 2 nodes
        for n in [5, 7]:
            eftInner2.setTermScaling(n*8 + 2, 1, [1])

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
        radiansPerElementAcross = math.pi/elementsCountAcross
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        wallThicknessSeptum = wallThickness[0]
        for n3 in range(2):
            sign = -1.0 if (n3 == 0) else 1.0
            baseY = wallThicknessSeptum
            radiusY = 0.5*wallThicknessSeptum
            radiusX = 0.5 - wallThickness[n3]
            for n2 in range(elementsCountAlong + 1):
                x[2] = n2 / elementsCountAlong
                for n1 in range(elementsCountAcross + 3):
                    if (n1 == 0) or (n1 == (elementsCountAcross + 2)):
                        x[0] = -0.5 if (n1 == 0) else 0.5
                        x[1] = -sign*baseY*0.5
                        dx_ds1[0] = 0.0
                        dx_ds1[1] = -1.0*baseY*(1.0 if (n1 == 0) else -1.0)
                        dx_ds3[0] = wallThickness[n3]*(-1.0 if (n1 == 0) else 1.0)
                        dx_ds3[1] = 0.0
                    elif (n1 == 1) or (n1 == (elementsCountAcross + 1)):
                        x[0] = -radiusX if (n1 == 1) else radiusX
                        x[1] = sign*radiusY*(-1.0 - 4.0*radiusX/(elementsCountAcross + 1))
                        angle_radians = math.pi/4
                        dx_ds1[0] = -sign*2.0*radiusX/(elementsCountAcross + 1)*math.cos(angle_radians)
                        dx_ds1[1] = 2.0*radiusX/(elementsCountAcross + 1)*math.sin(angle_radians)
                        dx_ds3[0] = wallThickness[n3]*math.sin(angle_radians)
                        dx_ds3[1] = sign*wallThickness[n3]*math.cos(angle_radians)
                        if n1 == 1:
                            dx_ds1[1] = -dx_ds1[1]
                            dx_ds3[0] = -dx_ds3[0]
                    else:
                        f1 = (n1 - 1)/elementsCountAcross
                        f2 = 1.0 - f1
                        x[0] = (f1 - f2)*radiusX
                        x[1] = -sign*radiusY
                        dx_ds1[0] = 2.0*radiusX/(elementsCountAcross + 1)
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

        scaleFactorsOuter = [
            math.cos(0.25*math.pi), math.cos(0.25*math.pi), math.cos(0.25*math.pi), -math.cos(0.25*math.pi),
            math.cos(0.25*math.pi), math.cos(0.25*math.pi), math.cos(0.25*math.pi), -math.cos(0.25*math.pi)
        ]
        scaleFactorsInner1 = [-1,
            math.cos(0.25*math.pi), math.cos(0.25*math.pi), math.cos(0.25*math.pi), math.cos(0.25*math.pi),
            math.cos(0.25*math.pi), -math.cos(0.25*math.pi), math.cos(0.25*math.pi), -math.cos(0.25*math.pi)
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
            result = element.setScaleFactor(eftInner1, 1, -1.0)
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
            result = element.setScaleFactor(eftInner2, 1, -1.0)
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

        fm.endChange()

