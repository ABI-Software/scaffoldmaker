"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valves regions.
"""

import math
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_heartventricles1 import MeshType_3d_heartventricles1
from mapclientplugins.meshgeneratorstep.utils.eft_utils import *
from mapclientplugins.meshgeneratorstep.utils.zinc_utils import *
from mapclientplugins.meshgeneratorstep.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventriclesbase1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Ventricles with Base 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'LV wall thickness' : 0.15,
            'LV wall thickness ratio apex' : 0.5,
            'LV wall thickness ratio base': 0.5,
            'RV free wall thickness' : 0.05,
            'RV width' : 0.15,
            'Length ratio' : 2.0,
            'Element length ratio equator/apex' : 1.0,
            'Base height' : 0.1,
            'Base thickness' : 0.05,
            'LV outlet inner diameter' : 0.25,
            'LV outlet wall thickness' : 0.02,
            'RV outlet inner diameter' : 0.25,
            'RV outlet wall thickness' : 0.02,
            'Outlet element length' : 0.1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'LV wall thickness ratio base',
            'RV free wall thickness',
            'RV width',
            'Length ratio',
            'Element length ratio equator/apex',
            'Base height',
            'Base thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Outlet element length'
        ]

    @staticmethod
    def checkOptions(options):
        if options['LV wall thickness'] < 0.0:
            options['LV wall thickness'] = 0.0
        elif options['LV wall thickness'] > 0.5:
            options['LV wall thickness'] = 0.5
        if options['LV wall thickness ratio apex'] < 0.0:
            options['LV wall thickness ratio apex'] = 0.0
        if options['LV wall thickness ratio base'] < 0.0:
            options['LV wall thickness ratio base'] = 0.0
        if options['RV free wall thickness'] < 0.0:
            options['RV free wall thickness'] = 0.0
        if options['RV width'] < 0.0:
            options['RV width'] = 0.0
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6
        for key in options:
            if key in [
                'Base height',
                'Base thickness',
                'LV outlet inner diameter',
                'LV outlet wall thickness',
                'RV outlet inner diameter',
                'RV outlet wall thickness',
                'Outlet element length']:
                if options[key] < 0.0:
                    options[key] = 0.0

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        lvWallThickness = options['LV wall thickness']
        lvWallThicknessRatioBase = options['LV wall thickness ratio base']
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        lvOutletRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        rvOutletRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        outletElementLength = options['Outlet element length']
        useCrossDerivatives = False

        # generate default heart ventricles model to add base plane to
        heartVentriclesOptions = MeshType_3d_heartventricles1.getDefaultOptions()
        for key in [
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'LV wall thickness ratio base',
            'RV free wall thickness',
            'RV width',
            'Length ratio',
            'Element length ratio equator/apex']:
            heartVentriclesOptions[key] = options[key]
        MeshType_3d_heartventricles1.generateMesh(region, heartVentriclesOptions)

        fm = region.getFieldmodule()
        fm.beginChange()
        # find the coordinates field created for the sphere shell
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateFull = nodes.createNodetemplate()
        nodetemplateFull.defineField(coordinates)
        nodetemplateFull.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateFull.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateFull.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateFull.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # nodes used only in bicubic-linear elements do not have D_DS3 parameters
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = startNodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

        lvDisplacementRadians = math.pi*2.0/3.0
        coslvDisplacementRadians = math.cos(lvDisplacementRadians)
        sinlvDisplacementRadians = math.sin(lvDisplacementRadians)

        outletJoinRadians = math.pi  # *0.8
        cosOutletJoinRadians = math.cos(outletJoinRadians)
        sinOutletJoinRadians = math.sin(outletJoinRadians)
        lvOutletOffset = 0.5 - lvWallThickness*(1.0 - lvWallThicknessRatioBase) - lvOutletRadius - lvOutletWallThickness
        #lvOutletCentreX = cosOutletJoinRadians*lvOutletOffset
        #lvOutletCentreY = sinOutletJoinRadians*lvOutletOffset
        lvOutletCentreX = lvOutletOffset*math.cos(lvDisplacementRadians)
        lvOutletCentreY = lvOutletOffset*math.sin(lvDisplacementRadians)

        # LV outlet - for bicubic-linear tube connection
        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        outletScale3 = 0.15
        for n3 in range(2):
            radius = lvOutletRadius + lvOutletWallThickness*n3
            for n1 in range(elementsCountAround):
                radiansAround = outletJoinRadians - math.pi + n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = lvOutletCentreX + radius*cosRadiansAround
                x[1] = lvOutletCentreY + radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                nodetemplate = nodetemplateLinearS3 if ((n3 == 0) or (n1 == 3)) else nodetemplateFull
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if nodetemplate is nodetemplateFull:
                    dx_ds3[0] = outletScale3*cosRadiansAround
                    dx_ds3[1] = outletScale3*sinRadiansAround
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # RV outlet - for bicubic-linear tube connection
        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        rvOutletOffset = lvOutletRadius + lvOutletWallThickness + rvOutletWallThickness + rvOutletRadius
        rvOutletCentreX = lvOutletCentreX + cosOutletJoinRadians*rvOutletOffset
        rvOutletCentreY = lvOutletCentreY + sinOutletJoinRadians*rvOutletOffset

        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        for n3 in range(2):
            radius = rvOutletRadius + rvOutletWallThickness*n3
            for n1 in range(elementsCountAround):
                if (n3 == 1) and (n1 == 0):
                    continue  # node is common with LV outlet
                radiansAround = outletJoinRadians - math.pi + n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = rvOutletCentreX + radius*cosRadiansAround
                x[1] = rvOutletCentreY + radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                nodeIdentifier += 1


        # create nodes on top and bottom of supraventricular crest centre where it meets
        # between nodes 141 and 130
        for i in range(2):
            radiansAround = outletJoinRadians - math.pi + (1.5 + 0.5*i)*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            scale1 = 0.15
            distance = 0.15
            if i == 0:
                scale2 = 0.25
            else:
                scale2 = 0.1
            x = [
                lvOutletCentreX + cosRadiansAround*(lvOutletRadius + lvOutletWallThickness + distance),
                lvOutletCentreY + sinRadiansAround*(lvOutletRadius + lvOutletWallThickness + distance),
                baseHeight
            ]
            dx_ds1 = [ -sinRadiansAround*scale1, cosRadiansAround*scale1, 0.0 ]
            dx_ds2 = [ -cosRadiansAround*scale2, -sinRadiansAround*scale2, 0.0 ]
            dx_ds3 = [ 0.0, 0.0, baseThickness ]
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            if i == 1:
                crest_nid1 = nodeIdentifier
            else:
                ra_nid1 = nodeIdentifier
            nodeIdentifier += 1

            x[2] = baseHeight + baseThickness
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            if i == 1:
                crest_nid2 = nodeIdentifier
            else:
                ra_nid2 = nodeIdentifier
            nodeIdentifier += 1

        # place nodes below RV outlet nodes 143, 144
        for nid in [143, 144]:
            node1 = nodes.findNodeByIdentifier(nid)
            cache.setNode(node1)
            result, x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            x[2] -= baseThickness
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            nodeIdentifier += 1
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # create extra nodes on LV 'crest' and edge of RA
        for i in range(2):
            radiansAround = math.pi*(1.4 - 0.2*i)
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            scale1 = 0.25
            distance = 0.17 + 0.08*i
            scale2 = 0.25
            x = [ cosRadiansAround*distance, sinRadiansAround*distance, baseHeight ]
            dx_ds1 = [ -sinRadiansAround*scale1, cosRadiansAround*scale1, 0.0 ]
            dx_ds2 = [ -cosRadiansAround*scale2, -sinRadiansAround*scale2, 0.0 ]
            dx_ds3 = [ 0.0, 0.0, baseThickness ]
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier += 1

            x[2] = baseHeight + baseThickness
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier += 1


        # create elements
        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        # crux of the heart - collapsed on 1 side, and partially on the other
        eft1 = tricubichermite.createEftBasic()
        ln_map = [ 1, 2, 1, 3, 4, 5, 1, 6 ]
        remapEftLocalNodes(eft1, 6, ln_map)
        setEftScaleFactorIds(eft1, [1], [])
        for n in [0, 2, 4]:
            ln = ln_map[n]
            eft1.setFunctionNumberOfTerms(n*8 + 3, 0)
        for n in [2]:
            ln = ln_map[n]
            eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        for n in [4, 6]:
            ln = ln_map[n]
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
        for n in [2, 6]:
            ln = ln_map[n]
            eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 3, 6, 1)
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 39, 40, 133, 88, 89, 139 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element crux', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1

        # Mitral - Aorta separator, collapsed on one end and one side
        eft1 = tricubichermite.createEftBasic()
        ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
        remapEftLocalNodes(eft1, 6, ln_map)
        setEftScaleFactorIds(eft1, [1], [])
        for n in [3, 7]:
            ln = ln_map[n]
            #eft1.setFunctionNumberOfTerms(n*8 + 3, 0)
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
        for n in [0, 1, 4, 5]:
            ln = ln_map[n]
            eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        for n in [1, 5]:
            ln = ln_map[n]
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
        for n in [0]:
            ln = ln_map[n]
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
        tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 5, 1)
        tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 6, 1)
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 161, 39, 138, 133, 144, 139 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element Mitral - Aorta separator', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1



        eftLinearOutlet = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftLinearOutlet, [1], [])
        # transition to no derivatives at outlet
        tricubichermite.setEftLinearDerivativeXi3(eftLinearOutlet, 3, 7, 3, 7, 1)
        tricubichermite.setEftLinearDerivativeXi3(eftLinearOutlet, 4, 8, 4, 8, 1)
        elementtemplateLinearOutlet = mesh.createElementtemplate()
        elementtemplateLinearOutlet.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateLinearOutlet.defineField(coordinates, -1, eftLinearOutlet)

        if False:
            # triangle wedges on top of ventricular septum
            # collapsed at xi2 = 1, angling in xi3 = 0 face
            eftSeptumWedge = tricubichermite.createEftBasic()
            eftSeptumWedge.setNumberOfLocalNodes(6)
            setEftScaleFactorIds(eftSeptumWedge, [1], [])
            for ln in [5, 6]:
                n = ln - 1
                # general map d/ds2 to get angle:
                mapEftFunction1Node2Terms(eftSeptumWedge, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [], Node.VALUE_LABEL_D_DS3, 1, [1])
            for ln in [3, 4]:
                n = ln + 3
                # shape value and d/ds1
                eftSeptumWedge.setTermNodeParameter(n*8 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                eftSeptumWedge.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                # general map d/ds2 to get angle:
                mapEftFunction1Node2Terms(eftSeptumWedge, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [], Node.VALUE_LABEL_D_DS3, 1, [1])
                # zero d/dxi3
                eftSeptumWedge.setFunctionNumberOfTerms(n*8 + 5, 0)
            print('eftSeptumWedge.validate()', eftSeptumWedge.validate())
            elementtemplateSeptumWedge = mesh.createElementtemplate()
            elementtemplateSeptumWedge.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplateSeptumWedge.defineField(coordinates, -1, eftSeptumWedge)

        septumNids = [
            [ 40, 41, 133, 134, 89, 90, 139, 140 ],
            [ 41, 42, 134, 135, 90, 91, 140, 141 ],
            [ 42, 43, 135, 136, 91, 92, 141, 142 ],
            [ 43, 44, 136, 137, 92, 93, 142, 143 ]
        ]
        for i in range(4):
            element = mesh.createElement(elementIdentifier, elementtemplateLinearOutlet)
            result2 = element.setNodesByIdentifier(eftLinearOutlet, septumNids[i])
            result3 = element.setScaleFactors(eftLinearOutlet, [-1])
            print('create septum element', elementIdentifier, result2, result3, septumNids[i])
            elementIdentifier += 1

        # supraventricular crest at RA by RV free wall
        eft1 = tricubichermite.createEftBasic()
        eft1.setTermNodeParameter(0*8 + 3, 1, 1, Node.VALUE_LABEL_D_DS1, 1)
        eft1.setTermNodeParameter(4*8 + 3, 1, 5, Node.VALUE_LABEL_D_DS1, 1)
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 111, 112, 156, 158, 129, 130, 157, 159 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        print('create element supra RA 1', elementIdentifier, result2, nodeIdentifiers1 )
        elementIdentifier += 1

        # supraventricular crest at PO by RV free wall
        eft1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eft1, [1], [])
        tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
        for ln in [4,8]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            #mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 112, 113, 158, 147, 130, 131, 159, 152 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element supra PO 1', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1

        # supraventricular crest at RA by RV septum
        eft1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eft1, [1], [])
        for ln in [3, 4]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        for ln in [7]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        for ln in [8]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 156, 158, 90, 91, 157, 159, 140, 141 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element supra RA 2', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1

        # supraventricular crest at PO by RV septum
        eft1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eft1, [1], [])
        tricubichermite.setEftLinearDerivativeXi3(eft1, 2, 6, 2, 6, 1)
        tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
        for ln in [3]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        for ln in [7]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        for ln in [2,6]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        for ln in [4,8]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 158, 147, 91, 146, 159, 152, 141, 151 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element supra PO 2', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1

        # supraventricular crest at PO-AO joiner
        # This element has some penetraton when very thin. Possible adding a rake angle to septal elements will help
        eft1 = tricubichermite.createEftBasic()
        eft1.setNumberOfLocalNodes(7)
        setEftScaleFactorIds(eft1, [1], [])
        tricubichermite.setEftLinearDerivativeXi3(eft1, 2, 6, 2, 6, 1)
        #tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 7, 4, 7, 1)
        tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 7, 1)
        for ln in [1]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        #print('eft1.validate() 1', eft1.validate())
        for ln in [3]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            #eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        #print('eft1.validate() 2', eft1.validate())
        for ln in [2, 4, 6]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        for ln in [5]:
            n = ln - 1
            mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
        #print('eft1.validate() 3', eft1.validate())
        for n in [6, 7]:
            ln = 7  # collapse on one edge
            mapEftFunction1Node1Term(eft1, n*8 + 1, ln, Node.VALUE_LABEL_VALUE, 1, [])
            eft1.setFunctionNumberOfTerms(n*8 + 2, 0)
            #eft1.setFunctionNumberOfTerms(n*8 + 3, 0)
            mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            #eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        for n in [6]:
            ln = 7  # collapse on one edge
            mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            #eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
        print('eft1.validate() 4', eft1.validate())
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        element = mesh.createElement(elementIdentifier, elementtemplate1)
        nodeIdentifiers1 = [ 91, 146, 92, 145, 141, 151, 142 ]
        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
        result3 = element.setScaleFactors(eft1, [-1.0])
        print('create element supra PO AO 3', elementIdentifier, result2, result3, nodeIdentifiers1 )
        elementIdentifier += 1


        if False:
            eft1 = tricubichermite.createEftBasic()
            setEftScaleFactorIds(eft1, [1], [])
            mapEftFunction1Node1Term(eft1, 0*8 + 2, 1, Node.VALUE_LABEL_D_DS2, 1, [1])
            mapEftFunction1Node1Term(eft1, 0*8 + 3, 1, Node.VALUE_LABEL_D_DS1, 1, [])
            mapEftFunction1Node1Term(eft1, 1*8 + 2, 2, Node.VALUE_LABEL_D_DS2, 1, [1])
            mapEftFunction1Node1Term(eft1, 4*8 + 2, 5, Node.VALUE_LABEL_D_DS2, 1, [1])
            mapEftFunction1Node1Term(eft1, 4*8 + 3, 5, Node.VALUE_LABEL_D_DS1, 1, [])
            mapEftFunction1Node1Term(eft1, 5*8 + 2, 6, Node.VALUE_LABEL_D_DS2, 1, [1])
            # transition to no derivatives at outlet
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 156, 112, 146, 147, 157, 130, 151, 152 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            rvNids = [
                [147, 148, 152, 153 ],
                [148, 149, 153, 154 ]
            ]

            i = 0
            for eid in [ 67, 68 ]:
                origElement = mesh.findElementByIdentifier(eid)
                origEft = origElement.getElementfieldtemplate(coordinates, -1)
                origNodeIdentifiers = getElementNodeIdentifiers(origElement, origEft)
                eft1 = origEft
                setEftScaleFactorIds(eft1, [1], [])
                # transition to no derivatives at outlet
                tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                nids = rvNids[i]
                nodeIdentifiers1 = [ origNodeIdentifiers[2], origNodeIdentifiers[3], nids[0], nids[1], origNodeIdentifiers[6], origNodeIdentifiers[7], nids[2], nids[3] ]
                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
                result3 = element.setScaleFactors(eft1, [-1.0])
                print('create element', elementIdentifier, result2, result3, nodeIdentifiers1 )
                elementIdentifier += 1
                i += 1


        rvNids = [
            [147, 148, 152, 153 ],
            [148, 149, 153, 154 ]
        ]

        scalefactors5hanging = [ -1.0, 0.5, 0.25, 0.125, 0.75 ]

        i = -1
        for eid in [ 68 ]:
            i += 1
            origElement = mesh.findElementByIdentifier(eid)
            origEft = origElement.getElementfieldtemplate(coordinates, -1)
            origNodeIdentifiers = getElementNodeIdentifiers(origElement, origEft)

            eft1 = tricubichermite.createEftBasic()
            eft2 = tricubichermite.createEftBasic()

            # general scale factors 1 -> 1, 102 -> 1/2, 104 -> 1/4, 108 -> 1/8, 304 -> 3/4
            setEftScaleFactorIds(eft1, [1, 102, 104, 108, 304], [])
            setEftScaleFactorIds(eft2, [1, 102, 104, 108, 304], [])

            tricubichermite.setEftMidsideXi1HangingNode(eft1, 2, 1, 1, 2, [1, 2, 3, 4, 5])
            tricubichermite.setEftMidsideXi1HangingNode(eft1, 6, 5, 5, 6, [1, 2, 3, 4, 5])
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)

            tricubichermite.setEftMidsideXi1HangingNode(eft2, 1, 2, 1, 2, [1, 2, 3, 4, 5])
            tricubichermite.setEftMidsideXi1HangingNode(eft2, 5, 6, 5, 6, [1, 2, 3, 4, 5])
            tricubichermite.setEftLinearDerivativeXi3(eft2, 3, 7, 3, 7, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft2, 4, 8, 4, 8, 1)

            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nids = rvNids[i*2]
            nodeIdentifiers1 = [ origNodeIdentifiers[2], origNodeIdentifiers[3], nids[0], nids[1], origNodeIdentifiers[6], origNodeIdentifiers[7], nids[2], nids[3] ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, scalefactors5hanging)
            print('create hanging element 1', elementIdentifier, result2, result3, nodeIdentifiers1, scalefactors5hanging)
            elementIdentifier += 1

            elementtemplate2 = mesh.createElementtemplate()
            elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate2.defineField(coordinates, -1, eft2)
            element = mesh.createElement(elementIdentifier, elementtemplate2)
            nids = rvNids[i*2 + 1]
            nodeIdentifiers2 = [ origNodeIdentifiers[2], origNodeIdentifiers[3], nids[0], nids[1], origNodeIdentifiers[6], origNodeIdentifiers[7], nids[2], nids[3] ]
            result2 = element.setNodesByIdentifier(eft2, nodeIdentifiers2)
            result3 = element.setScaleFactors(eft2, scalefactors5hanging)
            print('create hanging element 2', elementIdentifier, result2, result3, nodeIdentifiers2, scalefactors5hanging)
            elementIdentifier += 1


        fm.endChange()

