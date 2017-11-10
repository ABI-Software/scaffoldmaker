"""
Generates a 3-D unit sphere shell mesh with variable numbers of elements
around, up and through the thickness.
"""

import math
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_sphereshell1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Sphere Shell 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements up' : 4,
            'Number of elements around' : 4,
            'Number of elements through wall' : 1,
            'Wall thickness' : 0.25,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements around',
            'Number of elements through wall',
            'Wall thickness',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements up'] < 2) :
            options['Number of elements up'] = 2
        if (options['Number of elements around'] < 2) :
            options['Number of elements around'] = 2
        if (options['Number of elements through wall'] < 1) :
            options['Number of elements through wall'] = 1
        if (options['Wall thickness'] < 0.0) :
            options['Wall thickness'] = 0.0
        elif (options['Wall thickness'] > 0.5) :
            options['Wall thickness'] = 0.5

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountUp = options['Number of elements up']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
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

        mesh = fm.findMeshByDimension(3)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(8):
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # Apex1: collapsed on xi1 = 0
        eftApex1 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftApex1.setNumberOfLocalNodes(6)
        eftApex1.setNumberOfLocalScaleFactors(13)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftApex1.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftApex1.setScaleFactorIdentifier(1, 1)
        for layer in range(2):
            so = layer*6 + 1
            no = layer*3
            fo = layer*32
            for s in range(6):
                si = so + s + 1
                # 3 scale factors per node: cos(theta), sin(theta), arc angle radians
                sid = (s // 3)*100 + s + 1  # add 100 for different 'version'
                eftApex1.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                eftApex1.setScaleFactorIdentifier(si, sid)

            # basis node 1 -> local node 1
            ln = no + 1
            eftApex1.setTermNodeParameter(fo + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
            # 0 terms = zero parameter for d/dxi1 basis
            eftApex1.setFunctionNumberOfTerms(fo + 2, 0)
            # 2 terms for d/dxi2 via general linear map:
            eftApex1.setFunctionNumberOfTerms(fo + 3, 2)
            eftApex1.setTermNodeParameter(fo + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex1.setTermScaling(fo + 3, 1, [so + 1])
            eftApex1.setTermNodeParameter(fo + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex1.setTermScaling(fo + 3, 2, [so + 2])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            eftApex1.setFunctionNumberOfTerms(fo + 4, 2)
            eftApex1.setTermNodeParameter(fo + 4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex1.setTermScaling(fo + 4, 1, [so + 2, so + 3])
            eftApex1.setTermNodeParameter(fo + 4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex1.setTermScaling(fo + 4, 2, [1, so + 1, so + 3])
            eftApex1.setTermNodeParameter(fo + 5, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 0 terms = zero parameter for cross derivative 1 3, 2 3 and 1 2 3
            eftApex1.setFunctionNumberOfTerms(fo + 6, 0)
            eftApex1.setFunctionNumberOfTerms(fo + 7, 0)
            eftApex1.setFunctionNumberOfTerms(fo + 8, 0)

            # basis node 2 -> local node 1
            eftApex1.setTermNodeParameter(fo + 9, 1, ln, Node.VALUE_LABEL_VALUE, 1)
            # 0 terms = zero parameter for d/dxi1 basis
            eftApex1.setFunctionNumberOfTerms(fo + 10, 0)
            # 2 terms for d/dxi2 via general linear map:
            eftApex1.setFunctionNumberOfTerms(fo + 11, 2)
            eftApex1.setTermNodeParameter(fo + 11, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex1.setTermScaling(fo + 11, 1, [so + 4])
            eftApex1.setTermNodeParameter(fo + 11, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex1.setTermScaling(fo + 11, 2, [so + 5])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            eftApex1.setFunctionNumberOfTerms(fo + 12, 2)
            eftApex1.setTermNodeParameter(fo + 12, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex1.setTermScaling(fo + 12, 1, [so + 5, so + 6])
            eftApex1.setTermNodeParameter(fo + 12, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex1.setTermScaling(fo + 12, 2, [1, so + 4, so + 6])
            eftApex1.setTermNodeParameter(fo + 13, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 0 terms = zero parameter for cross derivative 1 3, 2 3 and 1 2 3
            eftApex1.setFunctionNumberOfTerms(fo + 14, 0)
            eftApex1.setFunctionNumberOfTerms(fo + 15, 0)
            eftApex1.setFunctionNumberOfTerms(fo + 16, 0)

            # basis nodes 3, 4 -> regular local nodes 2, 3
            for bn in range(2,4):
                fo2 = fo + bn*8
                ni = no + bn
                eftApex1.setTermNodeParameter(fo2 + 1, 1, ni, Node.VALUE_LABEL_VALUE, 1)
                eftApex1.setTermNodeParameter(fo2 + 2, 1, ni, Node.VALUE_LABEL_D_DS1, 1)
                eftApex1.setTermNodeParameter(fo2 + 3, 1, ni, Node.VALUE_LABEL_D_DS2, 1)
                eftApex1.setTermNodeParameter(fo2 + 5, 1, ni, Node.VALUE_LABEL_D_DS3, 1)
                if useCrossDerivatives:
                    eftApex1.setTermNodeParameter(fo2 + 4, 1, ni, Node.VALUE_LABEL_D2_DS1DS2, 1)
                    eftApex1.setTermNodeParameter(fo2 + 6, 1, ni, Node.VALUE_LABEL_D2_DS1DS3, 1)
                    eftApex1.setTermNodeParameter(fo2 + 7, 1, ni, Node.VALUE_LABEL_D2_DS2DS3, 1)
                    eftApex1.setTermNodeParameter(fo2 + 8, 1, ni, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                else:
                     eftApex1.setFunctionNumberOfTerms(fo2 + 4, 0)
                     eftApex1.setFunctionNumberOfTerms(fo2 + 6, 0)
                     eftApex1.setFunctionNumberOfTerms(fo2 + 7, 0)
                     eftApex1.setFunctionNumberOfTerms(fo2 + 8, 0)

        # Apex2: collapsed on xi1 = 1
        eftApex2 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftApex2.setNumberOfLocalNodes(6)
        eftApex2.setNumberOfLocalScaleFactors(8)
        eftApex2.setNumberOfLocalScaleFactors(13)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftApex2.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftApex2.setScaleFactorIdentifier(1, 1)
        for layer in range(2):
            so = layer*6 + 1
            no = layer*3
            fo = layer*32
            for s in range(6):
                si = so + s + 1
                # 3 scale factors per node: cos(theta), sin(theta), arc angle radians
                sid = (s // 3)*100 + s + 1  # add 100 for different 'version'
                eftApex2.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                eftApex2.setScaleFactorIdentifier(si, sid)

            # basis nodes 1, 2 -> regular local nodes 1, 2 (for each layer)
            for bn in range(2):
                fo2 = fo + bn*8
                ni = no + bn + 1
                eftApex2.setTermNodeParameter(fo2 + 1, 1, ni, Node.VALUE_LABEL_VALUE, 1)
                eftApex2.setTermNodeParameter(fo2 + 2, 1, ni, Node.VALUE_LABEL_D_DS1, 1)
                eftApex2.setTermNodeParameter(fo2 + 3, 1, ni, Node.VALUE_LABEL_D_DS2, 1)
                eftApex2.setTermNodeParameter(fo2 + 5, 1, ni, Node.VALUE_LABEL_D_DS3, 1)
                if useCrossDerivatives:
                    eftApex2.setTermNodeParameter(fo2 + 4, 1, ni, Node.VALUE_LABEL_D2_DS1DS2, 1)
                    eftApex2.setTermNodeParameter(fo2 + 6, 1, ni, Node.VALUE_LABEL_D2_DS1DS3, 1)
                    eftApex2.setTermNodeParameter(fo2 + 7, 1, ni, Node.VALUE_LABEL_D2_DS2DS3, 1)
                    eftApex2.setTermNodeParameter(fo2 + 8, 1, ni, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                else:
                     eftApex2.setFunctionNumberOfTerms(fo2 + 4, 0)
                     eftApex2.setFunctionNumberOfTerms(fo2 + 6, 0)
                     eftApex2.setFunctionNumberOfTerms(fo2 + 7, 0)
                     eftApex2.setFunctionNumberOfTerms(fo2 + 8, 0)

            # basis node 3 -> local node 3
            ln = no + 3
            fo3 = fo + 16
            eftApex2.setTermNodeParameter(fo3 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
            # 0 terms = zero parameter for d/dxi1 basis
            eftApex2.setFunctionNumberOfTerms(fo3 + 2, 0)
            # 2 terms for d/dxi2 via general linear map:
            eftApex2.setFunctionNumberOfTerms(fo3 + 3, 2)
            eftApex2.setTermNodeParameter(fo3 + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex2.setTermScaling(fo3 + 3, 1, [so + 1])
            eftApex2.setTermNodeParameter(fo3 + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex2.setTermScaling(fo3 + 3, 2, [so + 2])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            eftApex2.setFunctionNumberOfTerms(fo3 + 4, 2)
            eftApex2.setTermNodeParameter(fo3 + 4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex2.setTermScaling(fo3 + 4, 1, [1, so + 2, so + 3])
            eftApex2.setTermNodeParameter(fo3 + 4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex2.setTermScaling(fo3 + 4, 2, [so + 1, so + 3])
            eftApex2.setTermNodeParameter(fo3 + 5, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 0 terms = zero parameter for cross derivative 1 3, 2 3 and 1 2 3
            eftApex2.setFunctionNumberOfTerms(fo3 + 6, 0)
            eftApex2.setFunctionNumberOfTerms(fo3 + 7, 0)
            eftApex2.setFunctionNumberOfTerms(fo3 + 8, 0)

            # basis node 4 -> local node 3
            eftApex2.setTermNodeParameter(fo3 + 9, 1, ln, Node.VALUE_LABEL_VALUE, 1)
            # 0 terms = zero parameter for d/dxi1 basis
            eftApex2.setFunctionNumberOfTerms(fo3 + 10, 0)
            # 2 terms for d/dxi2 via general linear map:
            eftApex2.setFunctionNumberOfTerms(fo3 + 11, 2)
            eftApex2.setTermNodeParameter(fo3 + 11, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex2.setTermScaling(fo3 + 11, 1, [so + 4])
            eftApex2.setTermNodeParameter(fo3 + 11, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex2.setTermScaling(fo3 + 11, 2, [so + 5])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            eftApex2.setFunctionNumberOfTerms(fo3 + 12, 2)
            eftApex2.setTermNodeParameter(fo3 + 12, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftApex2.setTermScaling(fo3 + 12, 1, [1, so + 5, so + 6])
            eftApex2.setTermNodeParameter(fo3 + 12, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftApex2.setTermScaling(fo3 + 12, 2, [so + 4, so + 6])
            eftApex2.setTermNodeParameter(fo3 + 13, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 0 terms = zero parameter for cross derivative 1 3, 2 3 and 1 2 3
            eftApex2.setFunctionNumberOfTerms(fo3 + 14, 0)
            eftApex2.setFunctionNumberOfTerms(fo3 + 15, 0)
            eftApex2.setFunctionNumberOfTerms(fo3 + 16, 0)

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)
        elementtemplateApex1 = mesh.createElementtemplate()
        elementtemplateApex1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateApex1.defineField(coordinates, -1, eftApex1)
        elementtemplateApex2 = mesh.createElementtemplate()
        elementtemplateApex2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateApex2.defineField(coordinates, -1, eftApex2)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementUp = math.pi/elementsCountUp
        wallThicknessPerElement = wallThickness/elementsCountThroughWall
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        for n3 in range(elementsCountThroughWall + 1):
            radius = 0.5 + wallThickness*(n3/elementsCountThroughWall - 1.0)

            # create apex1 node
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -radius ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, radius*radiansPerElementUp, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, -wallThicknessPerElement ])
            nodeIdentifier = nodeIdentifier + 1

            # create regular rows between apexes
            for n2 in range(1, elementsCountUp):
                radiansUp = n2*radiansPerElementUp
                cosRadiansUp = math.cos(radiansUp);
                sinRadiansUp = math.sin(radiansUp);
                for n1 in range(elementsCountAround):
                    radiansAround = n1*radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        radius*cosRadiansAround*sinRadiansUp,
                        radius*sinRadiansAround*sinRadiansUp,
                        -radius*cosRadiansUp
                    ]
                    dx_ds1 = [
                        radius*-sinRadiansAround*sinRadiansUp*radiansPerElementAround,
                        radius*cosRadiansAround*sinRadiansUp*radiansPerElementAround,
                        0.0
                    ]
                    dx_ds2 = [
                        radius*cosRadiansAround*cosRadiansUp*radiansPerElementUp,
                        radius*sinRadiansAround*cosRadiansUp*radiansPerElementUp,
                        radius*sinRadiansUp*radiansPerElementUp
                    ]
                    dx_ds3 = [
                        wallThicknessPerElement*cosRadiansAround*sinRadiansUp,
                        wallThicknessPerElement*sinRadiansAround*sinRadiansUp,
                        -wallThicknessPerElement*cosRadiansUp
                    ]
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

            # create apex2 node
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, radius ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, -radius*radiansPerElementUp, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, wallThicknessPerElement ])
            nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        now = 2 + (elementsCountUp - 1)*elementsCountAround
        for e3 in range(elementsCountThroughWall):
            no = e3*now

            # create Apex1 elements, editing eft scale factor identifiers around apex
            # scale factor identifiers follow convention of offsetting by 100 for each 'version'
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1)%elementsCountAround
                eftApex1.setScaleFactorIdentifier(2, va*100 + 1)
                eftApex1.setScaleFactorIdentifier(3, va*100 + 2)
                eftApex1.setScaleFactorIdentifier(4, va*100 + 3)
                eftApex1.setScaleFactorIdentifier(5, vb*100 + 1)
                eftApex1.setScaleFactorIdentifier(6, vb*100 + 2)
                eftApex1.setScaleFactorIdentifier(7, vb*100 + 3)
                eftApex1.setScaleFactorIdentifier(8, va*100 + 1)
                eftApex1.setScaleFactorIdentifier(9, va*100 + 2)
                eftApex1.setScaleFactorIdentifier(10, va*100 + 3)
                eftApex1.setScaleFactorIdentifier(11, vb*100 + 1)
                eftApex1.setScaleFactorIdentifier(12, vb*100 + 2)
                eftApex1.setScaleFactorIdentifier(13, vb*100 + 3)
                # redefine field in template for changes to eftApex1:
                elementtemplateApex1.defineField(coordinates, -1, eftApex1)
                element = mesh.createElement(elementIdentifier, elementtemplateApex1)
                bni1 = no + 1
                bni2 = no + e1 + 2
                bni3 = no + (e1 + 1)%elementsCountAround + 2
                nodeIdentifiers = [ bni1, bni2, bni3, bni1 + now, bni2 + now, bni3 + now ]
                element.setNodesByIdentifier(eftApex1, nodeIdentifiers)
                # set general linear map coefficients
                radiansAround = e1*radiansPerElementAround
                radiansAroundNext = ((e1 + 1)%elementsCountAround)*radiansPerElementAround
                scalefactors = [
                    -1.0,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                ]
                result = element.setScaleFactors(eftApex1, scalefactors)
                elementIdentifier = elementIdentifier + 1

            # create regular rows between apexes
            for e2 in range(elementsCountUp - 2):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni11 = no + e2*elementsCountAround + e1 + 2
                    bni12 = no + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 2
                    bni21 = no + (e2 + 1)*elementsCountAround + e1 + 2
                    bni22 = no + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 2
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            # create Apex2 elements, editing eft scale factor identifiers around apex
            # scale factor identifiers follow convention of offsetting by 100 for each 'version'
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1)%elementsCountAround
                eftApex2.setScaleFactorIdentifier(2, va*100 + 1)
                eftApex2.setScaleFactorIdentifier(3, va*100 + 2)
                eftApex2.setScaleFactorIdentifier(4, va*100 + 3)
                eftApex2.setScaleFactorIdentifier(5, vb*100 + 1)
                eftApex2.setScaleFactorIdentifier(6, vb*100 + 2)
                eftApex2.setScaleFactorIdentifier(7, vb*100 + 3)
                eftApex2.setScaleFactorIdentifier(8, va*100 + 1)
                eftApex2.setScaleFactorIdentifier(9, va*100 + 2)
                eftApex2.setScaleFactorIdentifier(10, va*100 + 3)
                eftApex2.setScaleFactorIdentifier(11, vb*100 + 1)
                eftApex2.setScaleFactorIdentifier(12, vb*100 + 2)
                eftApex2.setScaleFactorIdentifier(13, vb*100 + 3)
                # redefine field in template for changes to eftApex2:
                elementtemplateApex1.defineField(coordinates, -1, eftApex2)
                element = mesh.createElement(elementIdentifier, elementtemplateApex1)
                bni3 = no + now
                bni1 = bni3 - elementsCountAround + e1
                bni2 = bni3 - elementsCountAround + (e1 + 1)%elementsCountAround
                nodeIdentifiers = [ bni1, bni2, bni3, bni1 + now, bni2 + now, bni3 + now ]
                element.setNodesByIdentifier(eftApex2, nodeIdentifiers)
                # set general linear map coefficients
                radiansAround = math.pi + e1*radiansPerElementAround
                radiansAroundNext = math.pi + ((e1 + 1)%elementsCountAround)*radiansPerElementAround
                scalefactors = [
                    -1.0,
                    -math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    -math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                    -math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    -math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                ]
                result = element.setScaleFactors(eftApex2, scalefactors)
                elementIdentifier = elementIdentifier + 1

        fm.endChange()
