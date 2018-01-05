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
            'Atria length' : 0.35,
            'Atria width' : 0.25,
            'Base height' : 0.2,
            'Base thickness' : 0.1,
            'LV outlet inner diameter' : 0.25,
            'LV outlet wall thickness' : 0.02,
            'LV outlet element length' : 0.1,
            'RV outlet inner diameter' : 0.25,
            'RV outlet wall thickness' : 0.02,
            'RV outlet element length' : 0.1,
            'Outlet spacing' : 0.05,
            'Atrial septum thickness' : 0.075,
            'Ventricular septum thickness' : 0.1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Atria length',
            'Atria width',
            'Base height',
            'Base thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'LV outlet element length',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'RV outlet element length',
            'Outlet spacing',
            'Atrial septum thickness',
            'Ventricular septum thickness'
        ]

    @staticmethod
    def checkOptions(options):
        for key in options:
            if isinstance(options[key], float):
                if options[key] < 0.0:
                    options[key] = 0.0

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        lvOutletRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        lvOutletElementLength = options['LV outlet element length']
        rvOutletRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        rvOutletElementLength = options['RV outlet element length']
        outletSpacing = options['Outlet spacing']
        useCrossDerivatives = False


        # generate heart ventricles model to add base pla
        heartVentriclesOptions = MeshType_3d_heartventricles1.getDefaultOptions()
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
        startNodeIdentifier = nodeIdentifier = getMaximumNodeIdentifier(nodes) + 1
        # LV outflow - for bicubic-linear tube connection
        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, lvOutletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        for n3 in range(2):
            radius = lvOutletRadius + lvOutletWallThickness*n3
            for n1 in range(elementsCountAround):
                radiansAround = n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = radius*cosRadiansAround
                x[1] = radius*sinRadiansAround
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

        # RV outflow - for bicubic-linear tube connection
        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        xCentre = -(lvOutletRadius + lvOutletWallThickness + outletSpacing + rvOutletRadius + rvOutletWallThickness)
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, rvOutletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        for n3 in range(2):
            radius = lvOutletRadius + lvOutletWallThickness*n3
            for n1 in range(elementsCountAround):
                radiansAround = n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = xCentre + radius*cosRadiansAround
                x[1] = radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                nodeIdentifier += 1

        # create elements
        elementIdentifier = getMaximumElementIdentifier(mesh) + 1
        if False:
            no2 = (elementsCount1 + 1)
            no3 = (elementsCount2 + 1)*no2
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1):
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni = e3*no3 + e2*no2 + e1 + 1
                        nodeIdentifiers = [ bni, bni + 1, bni + no2, bni + no2 + 1, bni + no3, bni + no3 + 1, bni + no2 + no3, bni + no2 + no3 + 1 ]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier += 1

        fm.endChange()

