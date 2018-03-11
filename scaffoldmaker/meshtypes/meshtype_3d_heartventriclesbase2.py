"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valves regions.
"""

from __future__ import division
import math
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_heartventricles2 import MeshType_3d_heartventricles2
from mapclientplugins.meshgeneratorstep.utils.eft_utils import *
from mapclientplugins.meshgeneratorstep.utils.zinc_utils import *
from mapclientplugins.meshgeneratorstep.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventriclesbase2(object):
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
            'RV width' : 0.2,
            'Length ratio' : 2.0,
            'Element length ratio equator/apex' : 1.0,
            'Atria length' : 0.35,
            'Atria width' : 0.25,
            'Base height' : 0.1,
            'Base thickness' : 0.1,
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
            'Atria length',
            'Atria width',
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
                'Atria length',
                'Atria width',
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
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        lvOutletRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        rvOutletRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        outletElementLength = options['Outlet element length']
        useCrossDerivatives = False

        # generate default heart ventricles model to add base plane to
        heartVentriclesOptions = MeshType_3d_heartventricles2.getDefaultOptions()
        for key in [
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'LV wall thickness ratio base',
            'RV free wall thickness',
            'RV width',
            'Length ratio',
            'Element length ratio equator/apex']:
            heartVentriclesOptions[key] = options[key]
        MeshType_3d_heartventricles2.generateMesh(region, heartVentriclesOptions)

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

        outletJoinRadians = math.pi  # *0.8
        cosOutletJoinRadians = math.cos(outletJoinRadians)
        sinOutletJoinRadians = math.sin(outletJoinRadians)
        lvOutletOffset = 0.5 - lvWallThickness - lvOutletRadius
        #lvOutletCentreX = cosOutletJoinRadians*lvOutletOffset
        #lvOutletCentreY = sinOutletJoinRadians*lvOutletOffset
        lvOutletCentreX = 0.0
        lvOutletCentreY = lvOutletOffset

        # LV outlet - for bicubic-linear tube connection
        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
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
                nodetemplate = nodetemplateLinearS3 if (n3 == 0) else nodetemplateFull
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if nodetemplate is nodetemplateFull:
                    dx_ds3[0] = lvOutletWallThickness*cosRadiansAround
                    dx_ds3[1] = lvOutletWallThickness*sinRadiansAround
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

        # create elements
        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        rvNids = [
            [104, 105, 110, 111 ],
            [105, 106, 111, 112 ],
            [106, 107, 112, 113 ],
            [107, 108, 113, 114 ],
            [116, 117, 113, 122 ],
            [117, 118, 122, 123 ],
            [118, 119, 123, 124 ],
            [119, 120, 124, 125 ]
        ]

        scalefactors5hanging = [ -1.0, 0.5, 0.25, 0.125, 0.75 ]

        i = -1
        for eid in range(55, 59):
            i += 1
            if i < 2:
                continue
            origElement = mesh.findElementByIdentifier(eid)
            origEft = origElement.getElementfieldtemplate(coordinates, -1)
            origNodeIdentifiers = getElementNodeIdentifiers(origElement, origEft)

            # first and last elements have scaling to use and fix (GRC to do)
            if False:  #eid == 55:
                eft1 = origEft
            else:
                eft1 = tricubichermite.createEftBasic()
            if False:  #eid == 58:
                eft2 = origEft
            else:
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

