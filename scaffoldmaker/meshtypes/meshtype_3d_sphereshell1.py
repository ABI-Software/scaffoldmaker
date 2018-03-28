"""
Generates a 3-D unit sphere shell mesh with variable numbers of elements
around, up and through the thickness.
"""

from __future__ import division
import math
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_sphereshell1:
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
            'Wall thickness ratio apex' : 1.0,
            'Exclude bottom rows' : 0,
            'Exclude top rows' : 0,
            'Length ratio' : 1.0,
            'Element length ratio equator/apex' : 1.0,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements around',
            'Number of elements through wall',
            'Exclude bottom rows',
            'Exclude top rows',
            'Wall thickness',
            'Wall thickness ratio apex',
            'Length ratio',
            'Element length ratio equator/apex',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 2:
            options['Number of elements around'] = 2
        if options['Number of elements through wall'] < 1:
            options['Number of elements through wall'] = 1
        if options['Exclude bottom rows'] < 0:
            options['Exclude bottom rows'] = 0
        if options['Exclude top rows'] < 0:
            options['Exclude top rows'] = 0
        if options['Wall thickness'] < 0.0:
            options['Wall thickness'] = 0.0
        elif options['Wall thickness'] > 0.5:
            options['Wall thickness'] = 0.5
        if options['Wall thickness ratio apex'] < 0.0:
            options['Wall thickness ratio apex'] = 0.0
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6

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
        useCrossDerivatives = options['Use cross derivatives']
        excludeBottomRows = options['Exclude bottom rows']
        excludeTopRows = options['Exclude top rows']
        wallThickness = options['Wall thickness']
        wallThicknessRatioApex = options['Wall thickness ratio apex']
        lengthRatio = options['Length ratio']
        elementLengthRatioEquatorApex = options['Element length ratio equator/apex']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)

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

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        # Apex1: collapsed on xi2 = 0
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
                sid = (s // 3)*100 + (s % 3) + 1  # add 100 for different 'version'
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

        # Apex2: collapsed on xi2 = 1
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
                sid = (s // 3)*100 + (s % 3) + 1  # add 100 for different 'version'
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
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]

        # pre-calculate positions and tangent/normal vectors up (elementsCountUp + 1) node layers
        outerWidth = 0.5
        outerLength = outerWidth*lengthRatio
        bOuter = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aOuter = 1.0 - bOuter

        innerWidth = outerWidth - wallThickness
        innerLength = outerLength - wallThickness*wallThicknessRatioApex
        lengthRatioInner = innerLength/innerWidth if (innerWidth > 0.0) else lengthRatio
        bInner = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aInner = 1.0 - bInner

        positionOuterArray = [(0,0)]*(elementsCountUp + 1)
        positionInnerArray = [(0,0)]*(elementsCountUp + 1)
        radiansUpOuterArray = [0]*(elementsCountUp + 1)
        radiansUpInnerArray = [0]*(elementsCountUp + 1)
        vector2OuterArray = [(0,0)]*(elementsCountUp + 1)
        vector2InnerArray = [(0,0)]*(elementsCountUp + 1)
        for n2 in range(elementsCountUp + 1):
            if n2*2 <= elementsCountUp:
                xi = n2*2 / elementsCountUp
            else:
                xi = 2.0 - (n2*2 / elementsCountUp)

            nxiOuter = aOuter*xi*xi + bOuter*xi
            dnxiOuter = 2.0*aOuter*xi + bOuter
            radiansUpOuterArray[n2] = radiansUpOuter = nxiOuter*math.pi*0.5 if (n2*2 <= elementsCountUp) else (math.pi - nxiOuter*math.pi*0.5)
            dRadiansUpOuter = dnxiOuter*math.pi/elementsCountUp
            cosRadiansUpOuter = math.cos(radiansUpOuter);
            sinRadiansUpOuter = math.sin(radiansUpOuter);
            positionOuterArray[n2] = positionOuter = ( outerWidth*sinRadiansUpOuter, -outerLength*cosRadiansUpOuter )
            vector2OuterArray[n2] = ( outerWidth*cosRadiansUpOuter*dRadiansUpOuter, outerLength*sinRadiansUpOuter*dRadiansUpOuter )

            nxiInner = aInner*xi*xi + bInner*xi
            dnxiInner = 2.0*aInner*xi + bInner
            radiansUpInnerArray[n2] = radiansUpInner = nxiInner*math.pi*0.5 if (n2*2 <= elementsCountUp) else (math.pi - nxiInner*math.pi*0.5)
            dRadiansUpInner = dnxiInner*math.pi/elementsCountUp
            cosRadiansUpInner = math.cos(radiansUpInner);
            sinRadiansUpInner = math.sin(radiansUpInner);
            positionInnerArray[n2] = positionInner = ( innerWidth*sinRadiansUpInner, -innerLength*cosRadiansUpInner )
            vector2InnerArray[n2] = ( innerWidth*cosRadiansUpInner*dRadiansUpInner, innerLength*sinRadiansUpInner*dRadiansUpInner )

        # now create the nodes
        for n3 in range(elementsCountThroughWall + 1):

            n3_fraction = n3 / elementsCountThroughWall

            for n2 in range(excludeBottomRows, elementsCountUp + 1 - excludeTopRows):

                positionOuter = positionOuterArray[n2]
                positionInner = positionInnerArray[n2]
                position = (
                    positionOuter[0]*n3_fraction + positionInner[0]*(1.0 - n3_fraction),
                    positionOuter[1]*n3_fraction + positionInner[1]*(1.0 - n3_fraction)
                )

                radiansUpOuter = radiansUpOuterArray[n2]
                sinRadiansUpOuter = math.sin(radiansUpOuter)
                radiansUpInner = radiansUpInnerArray[n2]
                sinRadiansUpInner = math.sin(radiansUpInner)

                vector2Outer = vector2OuterArray[n2]
                vector2Inner = vector2InnerArray[n2]
                vector2 = (
                    vector2Outer[0]*n3_fraction + vector2Inner[0]*(1.0 - n3_fraction),
                    vector2Outer[1]*n3_fraction + vector2Inner[1]*(1.0 - n3_fraction)
                )

                vector3 = (
                    (positionOuter[0] - positionInner[0])/elementsCountThroughWall,
                    (positionOuter[1] - positionInner[1])/elementsCountThroughWall
                )

                if n2 == 0:
                    # create apex1 node
                    node = nodes.createNode(nodeIdentifier, nodetemplateApex)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, position[1] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, vector2[0], 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ vector2[0], 0.0, 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, vector3[1] ])
                    nodeIdentifier = nodeIdentifier + 1

                elif n2 < elementsCountUp:
                    # create regular rows between apexes
                    for n1 in range(elementsCountAround):
                        radiansAround = n1*radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            position[0]*cosRadiansAround,
                            position[0]*sinRadiansAround,
                            position[1]
                        ]
                        dx_ds1 = [
                            position[0]*-sinRadiansAround*radiansPerElementAround,
                            position[0]*cosRadiansAround*radiansPerElementAround,
                            0.0
                        ]
                        dx_ds2 = [
                            vector2[0]*cosRadiansAround,
                            vector2[0]*sinRadiansAround,
                            vector2[1]
                        ]
                        dx_ds3 = [
                            vector3[0]*cosRadiansAround,
                            vector3[0]*sinRadiansAround,
                            vector3[1]
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

                else:
                    # create apex2 node
                    node = nodes.createNode(nodeIdentifier, nodetemplateApex)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, position[1] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, vector2[0], 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ -vector2[0], 0.0, 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, vector3[1] ])
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        # now (node offset through wall) varies with number of excluded rows
        now = 0;
        if excludeBottomRows == 0:
            now += 1  # bottom apex node
        if excludeTopRows == 0:
            now += 1  # top apex node
        fullNodeRows = elementsCountUp - 1
        if excludeBottomRows > 1:
            fullNodeRows -= (excludeBottomRows - 1)
        if excludeTopRows > 1:
            fullNodeRows -= (excludeTopRows - 1)
        if fullNodeRows > 0:
            now += fullNodeRows*elementsCountAround
        row2NodeOffset = 2 if (excludeBottomRows == 0) else 1
        for e3 in range(elementsCountThroughWall):
            no = e3*now

            if excludeBottomRows == 0:
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
            rowLimit = (elementsCountUp - 2)
            if excludeBottomRows > 1:
                rowLimit -= (excludeBottomRows - 1)
            if excludeTopRows > 1:
                rowLimit -= (excludeTopRows - 1)
            for e2 in range(0, rowLimit):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni11 = no + e2*elementsCountAround + e1 + row2NodeOffset
                    bni12 = no + e2*elementsCountAround + (e1 + 1)%elementsCountAround + row2NodeOffset
                    bni21 = no + (e2 + 1)*elementsCountAround + e1 + row2NodeOffset
                    bni22 = no + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + row2NodeOffset
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            if excludeTopRows == 0:
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

        if False:  # (wallThickness < 0.5):  # and (wallThicknessRatioApex != 1.0):
            r = fm.createFieldMagnitude(coordinates)
            const05 = fm.createFieldConstant([0.5])
            d = fm.createFieldSubtract(const05, r)
            d_r = fm.createFieldDivide(d, r)
            rRatio = 1.0 - wallThicknessRatioApex
            rScale = fm.createFieldConstant([rRatio])
            zScale = fm.createFieldMultiply(d_r, rScale)
            one = fm.createFieldConstant([1.0])
            one_plus_zScale = fm.createFieldAdd(one, zScale)
            scale = fm.createFieldConcatenate([one, one, one_plus_zScale])
            newCoordinates = fm.createFieldMultiply(coordinates, scale)
            fieldassignment = coordinates.createFieldassignment(newCoordinates)
            fieldassignment.assign()

        fm.endChange()
