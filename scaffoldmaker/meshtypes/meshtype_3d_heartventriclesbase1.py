"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valves regions.
"""

from __future__ import division
import math
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles1 import MeshType_3d_heartventricles1
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
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
        options = MeshType_3d_heartventricles1.getDefaultOptions()
        # only works with particular numbers of elements
        options['Number of elements up'] = 4
        options['Number of elements around'] = 11
        options['Number of elements across septum'] = 5
        options['Number of elements below septum'] = 2
        options['Number of elements through LV wall'] = 1
        # additional options
        options['Base height'] = 0.1
        options['Base thickness'] = 0.05
        options['LV outlet inner diameter'] = 0.25
        options['LV outlet wall thickness'] = 0.02
        options['RV outlet inner diameter'] = 0.25
        options['RV outlet wall thickness'] = 0.02
        options['Outlet element length'] = 0.1
        options['Outlet rotation degrees'] = 0.0
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventricles1.getOrderedOptionNames()
        # do not allow numbers of elements to be edited
        for hiddenOptionName in [
            'Number of elements up',
            'Number of elements around',
            'Number of elements across septum',
            'Number of elements below septum',
            'Number of elements through LV wall']:
            optionNames.remove(hiddenOptionName)
        optionNames += [
            'Base height',
            'Base thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Outlet element length',
            'Outlet rotation degrees'
        ]
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_heartventricles1.checkOptions(options)
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
        elementsCountUp = options['Number of elements up']
        elementsCountAround = options['Number of elements around']
        print('elementsCountAround', elementsCountAround)
        elementsCountAcrossSeptum = options['Number of elements across septum']
        lvWallThickness = options['LV wall thickness']
        lvWallThicknessRatioBase = options['LV wall thickness ratio base']
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        outletRotationRadians = options['Outlet rotation degrees']*math.pi/180.0
        lvOutletRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        rvOutletRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        outletElementLength = options['Outlet element length']
        LVBaseFlattenRatio = options['LV base flatten ratio']
        LVBaseFlattenAngleRadians = options['LV base flatten angle degrees']*math.pi/180.0
        useCrossDerivatives = False

        # generate default heart ventricles model to add base plane to
        MeshType_3d_heartventricles1.generateMesh(region, options)

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

        # node offsets for each row, and wall in LV, plus first LV node on inside top
        norl = elementsCountAround
        nowl = 1 + elementsCountUp*norl
        nidl = nowl - norl + 1
        # nodes offsets through RV wall, plus first RV node on inside top
        nowr = 18
        nidr = nowl*2 + nowr - 7

        radiansPerElementOrig = 2.0*math.pi/elementsCountAround
        septumCentreRadians = (1.0 + elementsCountAcrossSeptum/2.0)*radiansPerElementOrig
        lvDisplacementRadians = septumCentreRadians + LVBaseFlattenAngleRadians

        # GRC TODO check:
        lvOutletOffset = 0.5 - lvWallThickness*lvWallThicknessRatioBase - lvOutletRadius
        lvOutletCentreX = lvOutletOffset*math.cos(lvDisplacementRadians)
        lvOutletCentreY = lvOutletOffset*math.sin(lvDisplacementRadians)

        # LV outlet - for bicubic-linear tube connection

        elementsCountAroundOutlet = 6
        radiansPerElementAroundOutlet = 2.0*math.pi/elementsCountAroundOutlet
        lvOutletNodeId = [ [-1]*elementsCountAroundOutlet, [-1]*elementsCountAroundOutlet ]
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        # GRC Set properly:
        outletScale3 = 0.15
        for n3 in range(2):
            radius = lvOutletRadius + lvOutletWallThickness*n3
            for n1 in range(elementsCountAroundOutlet):
                radiansAround = outletRotationRadians + n1*radiansPerElementAroundOutlet
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = lvOutletCentreX + radius*cosRadiansAround
                x[1] = lvOutletCentreY + radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAroundOutlet*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAroundOutlet*radius*cosRadiansAround
                nodetemplate = nodetemplateLinearS3 if ((n3 == 0) or (n1 == 3)) else nodetemplateFull
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                lvOutletNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if nodetemplate is nodetemplateFull:
                    lvOutletInclineRadians = math.pi/6.0 if (n1 == 4) else 0.0
                    dx_ds3 = [
                        outletScale3*cosRadiansAround*math.cos(lvOutletInclineRadians),
                        outletScale3*sinRadiansAround*math.cos(lvOutletInclineRadians),
                        -outletScale3*math.sin(lvOutletInclineRadians)
                    ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                if n3 == 1:
                    if n1 == 0:
                        cruxCentreNodeId = nodeIdentifier
                        cruxCentre = [ x[0], x[1], x[2] ]
                    elif n1 == 1:
                        cruxRightNodeId = nodeIdentifier
                        cruxRight = [ x[0], x[1], x[2] ]
                    elif n1 == (elementsCountAroundOutlet - 1):
                        cruxLeftNodeId = nodeIdentifier
                        cruxLeft = [ x[0], x[1], x[2] ]
                nodeIdentifier += 1

        # RV outlet - for bicubic-linear tube connection
        rvOutletNodeId = [ [-1]*elementsCountAroundOutlet, [-1]*elementsCountAroundOutlet ]
        rvOutletNodeId[1][0] = lvOutletNodeId[1][elementsCountAroundOutlet // 2]
        rvOutletOffset = lvOutletRadius + lvOutletWallThickness + rvOutletWallThickness + rvOutletRadius
        rvOutletCentreX = lvOutletCentreX + math.cos(outletRotationRadians + math.pi)*rvOutletOffset
        rvOutletCentreY = lvOutletCentreY + math.sin(outletRotationRadians + math.pi)*rvOutletOffset

        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        for n3 in range(2):
            radius = rvOutletRadius + rvOutletWallThickness*n3
            for n1 in range(elementsCountAroundOutlet):
                if (n3 == 1) and (n1 == 0):
                    continue  # node is common with LV outlet
                radiansAround = outletRotationRadians + n1*radiansPerElementAroundOutlet
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = rvOutletCentreX + radius*cosRadiansAround
                x[1] = rvOutletCentreY + radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAroundOutlet*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAroundOutlet*radius*cosRadiansAround
                need_dx_ds3 = (n3 == 1) and ((n1 == 1) or (n1 == (elementsCountAroundOutlet - 1)))
                nodetemplate = nodetemplateFull if need_dx_ds3 else nodetemplateLinearS3
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                rvOutletNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if need_dx_ds3:
                    rvOutletInclineRadians = math.pi/6.0
                    dx_ds3 = [
                        outletScale3*cosRadiansAround*math.cos(rvOutletInclineRadians),
                        outletScale3*sinRadiansAround*math.cos(rvOutletInclineRadians),
                        -outletScale3*math.sin(rvOutletInclineRadians)
                    ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # create node mid septum
        nid1 = nowl*2 - norl + 5
        nid2 = lvOutletNodeId[1][2]
        for n in range(1):
            node = nodes.findNodeByIdentifier(nid1 + n)
            cache.setNode(node)
            result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            mag2 = math.sqrt(dx_ds2[0]*dx_ds2[0] + dx_ds2[1]*dx_ds2[1] + dx_ds2[2]*dx_ds2[2])
            node = nodes.findNodeByIdentifier(nid2 + n)
            cache.setNode(node)
            result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            midseptum_nid = nodeIdentifier
            cache.setNode(node)
            r1 = 0.5*(baseHeight + baseThickness)/mag2
            r2 = 0.25
            r1d = r1
            r2d = 0.5
            x_s = [
                x1[0] + r1*dx_ds2[0] + r2*dx_ds3[0],
                x1[1] + r1*dx_ds2[1] + r2*dx_ds3[1],
                x1[2] + r1*dx_ds2[2] + r2*dx_ds3[2]
            ]
            dx_ds1_s  = dx_ds1
            dx_ds2_s = [
                r1d*dx_ds2[0] + r2d*dx_ds3[0],
                r1d*dx_ds2[1] + r2d*dx_ds3[1],
                r1d*dx_ds2[2] + r2d*dx_ds3[2]
            ]
            dx_ds3_s = [ x_s[0] - x2[0], x_s[1] - x2[1], x_s[2] - x2[2] ]
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x_s)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1_s)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2_s)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3_s)
            nodeIdentifier += 1

        # Atria nodes
        #############
        elementsCountAroundAtria = 6

        # left atrium interface nodes

        atriumInletSlopeRadians = math.pi/6.0
        # GRC change from baseThickness?
        atriumInletSlopeLength = baseThickness*math.cos(atriumInletSlopeRadians)
        atriumInletSlopeHeight = baseThickness*math.sin(atriumInletSlopeRadians)

        # given: lvDisplacementRadians = short axis angle of LA
        laCentreOffset1 = lvOutletOffset - lvOutletRadius - lvOutletWallThickness - atriumInletSlopeLength
        # GRC fudge factor:
        laCentreOffset2 = 0.9*(-0.5 + lvWallThickness)
        laCenterOffset = 0.5*(laCentreOffset1 + laCentreOffset2)

        laUnitMinorX = math.cos(lvDisplacementRadians)
        laUnitMinorY = math.sin(lvDisplacementRadians)

        laCentreX = laCenterOffset*laUnitMinorX
        laCentreY = laCenterOffset*laUnitMinorY

        atriaInnerMinorMag = 0.5*(laCentreOffset1 - laCentreOffset2)
        laInnerMinorX = laUnitMinorX*atriaInnerMinorMag
        laInnerMinorY = laUnitMinorY*atriaInnerMinorMag
        atriaOuterMinorMag = atriaInnerMinorMag + atriumInletSlopeLength
        laOuterMinorX = laUnitMinorX*atriaOuterMinorMag
        laOuterMinorY = laUnitMinorY*atriaOuterMinorMag
        lvOutletCentreToAtriaCentre = lvOutletRadius + lvOutletWallThickness + atriaOuterMinorMag

        # GRC fudge factor:
        atriaInnerMajorMag = 0.9*(0.5 - lvWallThickness)*LVBaseFlattenRatio
        laInnerMajorX = atriaInnerMajorMag*laUnitMinorY
        laInnerMajorY = atriaInnerMajorMag*-laUnitMinorX
        atriaOuterMajorMag = atriaInnerMajorMag + atriumInletSlopeLength
        laUnitMajorX = laUnitMinorY
        laUnitMajorY = -laUnitMinorX
        laOuterMajorX = atriaOuterMajorMag*laUnitMajorX
        laOuterMajorY = atriaOuterMajorMag*laUnitMajorY

        # GRC revisit:
        atrialSeptumThickness = lvWallThickness*lvWallThicknessRatioBase

        # calculate radians to atrial septum, cruxLeft
        dx = laCentreX - lvOutletCentreX
        dy = laCentreY - lvOutletCentreY
        laCentreOutletRadians = math.atan2(dy, dx)

        # following is radians around LA
        laSeptumRadians = math.asin(atriaInnerMinorMag/lvOutletCentreToAtriaCentre)
        ox = cruxLeft[0] - laCentreX
        oy = cruxLeft[1] - laCentreY
        qx = ox*laUnitMajorX + oy*laUnitMajorY
        qy = ox*laUnitMinorX + oy*laUnitMinorY
        cruxLeftRadians = math.atan2(atriaOuterMajorMag*qy, atriaOuterMinorMag*qx)
        print('qx', qx, ' qy', qy, ' cruxLeftDegrees', cruxLeftRadians*180.0/math.pi)

        cosRadiansAround = math.cos(laSeptumRadians)
        sinRadiansAround = math.sin(laSeptumRadians)
        laSeptumX = laCentreX + cosRadiansAround*laInnerMajorX + sinRadiansAround*laInnerMinorX
        laSeptumY = laCentreY + cosRadiansAround*laInnerMajorY + sinRadiansAround*laInnerMinorY
        dx = laSeptumX - lvOutletCentreX
        dy = laSeptumY - lvOutletCentreY
        lvOutletCentreToLaSeptum = math.sqrt(dx*dx + dy*dy)
        # get radians (around centre of LV outlet) across half atrial septum
        halfAtrialSeptumOutletRadians = math.asin(0.5*atrialSeptumThickness/lvOutletCentreToLaSeptum)
        laSeptumOutletRadians = math.atan2(dy, dx)
        raSeptumOutletRadians = laSeptumOutletRadians + 2.0*halfAtrialSeptumOutletRadians
        raSeptumX = lvOutletCentreX + lvOutletCentreToLaSeptum*math.cos(raSeptumOutletRadians)
        raSeptumY = lvOutletCentreY + lvOutletCentreToLaSeptum*math.sin(raSeptumOutletRadians)

        laRadians = [0.0]*elementsCountAroundAtria
        laRadians[0] = laSeptumRadians
        laRadians[1] = cruxLeftRadians
        deltaRadiansL1 = laRadians[1] - laRadians[0]
        deltaRadiansLN = (2.0*math.pi - 2.0*deltaRadiansL1)/(elementsCountAroundAtria - 2.0)
        laDeltaRadians = [ deltaRadiansL1 ]*elementsCountAroundAtria
        laRadians[2] = laRadians[1] + 0.5*(deltaRadiansL1 + deltaRadiansLN)
        laDeltaRadians[2] = deltaRadiansLN
        for i in range(3, elementsCountAroundAtria):
            laRadians[i] = laRadians[i - 1] + deltaRadiansLN
            laDeltaRadians[i] = deltaRadiansLN

        laNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        laNodeId[1][1] = cruxLeftNodeId

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = laRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    laCentreX + cosRadiansAround*laInnerMajorX + sinRadiansAround*laInnerMinorX,
                    laCentreY + cosRadiansAround*laInnerMajorY + sinRadiansAround*laInnerMinorY,
                    cruxCentre[2] - atriumInletSlopeHeight ]
                outer = [
                    laCentreX + cosRadiansAround*laOuterMajorX + sinRadiansAround*laOuterMinorX,
                    laCentreY + cosRadiansAround*laOuterMajorY + sinRadiansAround*laOuterMinorY,
                    cruxCentre[2] ]
                if (n3 == 1) and (n1 < 2):
                    continue  # already have a node from crux or will get right atrial septrum
                node = nodes.createNode(nodeIdentifier, nodetemplateFull)
                laNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        laDeltaRadians[n1]*(-sinRadiansAround*laInnerMajorX + cosRadiansAround*laInnerMinorX),
                        laDeltaRadians[n1]*(-sinRadiansAround*laInnerMajorY + cosRadiansAround*laInnerMinorY),
                        0.0 ]
                else:
                    dx_ds1 = [
                        laDeltaRadians[n1]*(-sinRadiansAround*laOuterMajorX + cosRadiansAround*laOuterMinorX),
                        laDeltaRadians[n1]*(-sinRadiansAround*laOuterMajorY + cosRadiansAround*laOuterMinorY),
                        0.0 ]
                dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                dx_ds2 = [
                    dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                    dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                    dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                mag = math.sqrt(dx_ds2[0]*dx_ds2[0] + dx_ds2[1]*dx_ds2[1] + dx_ds2[2]*dx_ds2[2])
                for i in range(3):
                    # GRC check scaling here:
                    dx_ds2[i] *= inner[2]/mag
                if n1 == 0:
                    # GRC check scaling here:
                    dx_ds2 = [ 0.0, 0.0, inner[2] ]
                    dx_ds3 = [ raSeptumX - laSeptumX, raSeptumY - laSeptumY, 0.0 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if True:
            # show axes of left atrium
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, cruxCentre[2] - atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laInnerMajorX, laInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laInnerMinorX, laInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, cruxCentre[2] ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laOuterMajorX, laOuterMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laOuterMinorX, laOuterMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        # right atrium interface nodes

        raSeptumRadians = math.pi*2.0 - laSeptumRadians
        raCentreOutletRadians = raSeptumOutletRadians + (laSeptumOutletRadians - laCentreOutletRadians)
        raUnitMinorX = math.cos(raCentreOutletRadians)
        raUnitMinorY = math.sin(raCentreOutletRadians)
        raCentreX = lvOutletCentreX + lvOutletCentreToAtriaCentre*raUnitMinorX
        raCentreY = lvOutletCentreY + lvOutletCentreToAtriaCentre*raUnitMinorY
        raInnerMinorX = raUnitMinorX*atriaInnerMinorMag
        raInnerMinorY = raUnitMinorY*atriaInnerMinorMag
        raOuterMinorX = raUnitMinorX*atriaOuterMinorMag
        raOuterMinorY = raUnitMinorY*atriaOuterMinorMag
        raUnitMajorX = raUnitMinorY
        raUnitMajorY = -raUnitMinorX
        raInnerMajorX = raUnitMajorX*atriaInnerMajorMag
        raInnerMajorY = raUnitMajorY*atriaInnerMajorMag
        raOuterMajorX = raUnitMajorX*atriaOuterMajorMag
        raOuterMajorY = raUnitMajorY*atriaOuterMajorMag

        ox = cruxRight[0] - raCentreX
        oy = cruxRight[1] - raCentreY
        qx = ox*raUnitMajorX + oy*raUnitMajorY
        qy = ox*raUnitMinorX + oy*raUnitMinorY
        cruxRightRadians = math.atan2(atriaOuterMajorMag*qy, atriaOuterMinorMag*qx)
        if cruxRightRadians < math.pi:
            cruxRightRadians += 2.0*math.pi
        print('qx', qx, ' qy', qy, ' cruxRightDegrees', cruxRightRadians*180.0/math.pi)

        raRadians = [0.0]*elementsCountAroundAtria
        raRadians[0] = raSeptumRadians
        raRadians[-1] = 0.5*(raSeptumRadians + cruxRightRadians)
        raRadians[-2] = cruxRightRadians
        deltaRadiansR1 = 0.5*(raSeptumRadians - cruxRightRadians)
        deltaRadiansRN = (2.0*math.pi - 3.0*deltaRadiansR1)/(elementsCountAroundAtria - 3.0)
        raDeltaRadians = [ deltaRadiansR1 ]*elementsCountAroundAtria
        raRadians[1] = raSeptumRadians - math.pi*2.0 + 0.5*(deltaRadiansR1 + deltaRadiansRN)
        raDeltaRadians[1] = deltaRadiansRN
        for i in range(2, elementsCountAroundAtria - 2):
            raRadians[i] = raRadians[i - 1] + deltaRadiansRN
            raDeltaRadians[i] = deltaRadiansRN

        raNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        raNodeId[1][0] = laNodeId[0][0]
        raNodeId[1][-2] = cruxRightNodeId
        raNodeId[1][-1] = cruxCentreNodeId

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = raRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    raCentreX + cosRadiansAround*raInnerMajorX + sinRadiansAround*raInnerMinorX,
                    raCentreY + cosRadiansAround*raInnerMajorY + sinRadiansAround*raInnerMinorY,
                    cruxRight[2] - atriumInletSlopeHeight ]
                outer = [
                    raCentreX + cosRadiansAround*raOuterMajorX + sinRadiansAround*raOuterMinorX,
                    raCentreY + cosRadiansAround*raOuterMajorY + sinRadiansAround*raOuterMinorY,
                    cruxRight[2] ]
                if raNodeId[n3][n1] >= 0:
                    continue  # already have a node from crux or left atrium interface
                node = nodes.createNode(nodeIdentifier, nodetemplateFull)
                raNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        raDeltaRadians[n1]*(-sinRadiansAround*raInnerMajorX + cosRadiansAround*raInnerMinorX),
                        raDeltaRadians[n1]*(-sinRadiansAround*raInnerMajorY + cosRadiansAround*raInnerMinorY),
                         0.0 ]
                else:
                    dx_ds1 = [
                        raDeltaRadians[n1]*(-sinRadiansAround*raOuterMajorX + cosRadiansAround*raOuterMinorX),
                        raDeltaRadians[n1]*(-sinRadiansAround*raOuterMajorY + cosRadiansAround*raOuterMinorY),
                         0.0 ]
                dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                dx_ds2 = [
                    dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                    dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                    dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                mag = math.sqrt(dx_ds2[0]*dx_ds2[0] + dx_ds2[1]*dx_ds2[1] + dx_ds2[2]*dx_ds2[2])
                for i in range(3):
                    # GRC check scaling here:
                    dx_ds2[i] *= inner[2]/mag
                if n1 == 0:
                    # GRC check scaling here:
                    dx_ds1 = [ -x for x in dx_ds1 ]
                    dx_ds2 = [ 0.0, 0.0, inner[2] ]
                    dx_ds3 = [ raSeptumX - laSeptumX, raSeptumY - laSeptumY, 0.0 ]
                elif n1 > (elementsCountAroundAtria - 3):
                    # negate two derivatives to align with left side across ventricular septum
                    dx_ds1 = [ -x for x in dx_ds1 ]
                    dx_ds3 = [ -x for x in dx_ds3 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if True:
            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, cruxRight[2] - atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raInnerMajorX, raInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raInnerMinorX, raInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, cruxRight[2] ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raOuterMajorX, raOuterMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raOuterMinorX, raOuterMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        # transfer right atrium septum node to outside of left:
        laNodeId[1][0] = raNodeId[0][0]

        # fix dx_ds3 on crux nodes to be difference between node coordinates across atrium wall
        for nids in [ ( laNodeId[0][1], laNodeId[1][1] ), (raNodeId[0][4], raNodeId[1][4]) ]:
            node1 = nodes.findNodeByIdentifier(nids[0])
            cache.setNode(node1)
            result, x_i = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            node2 = nodes.findNodeByIdentifier(nids[1])
            cache.setNode(node2)
            result, x_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            dx_ds3 = [ x_o[i] - x_i[i] for i in range(3) ]
            if nids[0] is raNodeId[0][4]:
                dx_ds3 = [ -v for v in dx_ds3 ]
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            cache.setNode(node1)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # create nodes on top and bottom of supraventricular crest
        crestRadians = outletRotationRadians + math.pi - radiansPerElementAroundOutlet
        cosCrestRadians = math.cos(crestRadians)
        sinCrestRadians = math.sin(crestRadians)
        scale1 = 0.15
        distance = 0.20
        scale2 = 0.1
        x = [
            lvOutletCentreX + cosCrestRadians*(lvOutletRadius + lvOutletWallThickness + distance),
            lvOutletCentreY + sinCrestRadians*(lvOutletRadius + lvOutletWallThickness + distance),
            baseHeight
        ]
        dx_ds1 = [ -sinCrestRadians*scale1, cosCrestRadians*scale1, 0.0 ]
        innerCrestScaling = 0.5
        dx_ds2 = [ -cosCrestRadians*scale2*innerCrestScaling, -sinCrestRadians*scale2*innerCrestScaling, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, baseThickness ]
        node = nodes.createNode(nodeIdentifier, nodetemplateFull)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        crest_nid1 = nodeIdentifier
        nodeIdentifier += 1

        x[2] = baseHeight + baseThickness
        dx_ds2 = [ -cosCrestRadians*scale2, -sinCrestRadians*scale2, 0.0 ]
        node = nodes.createNode(nodeIdentifier, nodetemplateFull)
        cache.setNode(node)
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        crest_nid2 = nodeIdentifier
        nodeIdentifier += 1

        # create node on LV freewall base to help close

        # left crux coordinates and derivative
        node = nodes.findNodeByIdentifier(lvOutletNodeId[1][-1])
        cache.setNode(node)
        result, x_c = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, dx_ds1_c = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        # LV freewall node inside coordinates for computing dx_ds3
        node = nodes.findNodeByIdentifier(nowl - 3)
        cache.setNode(node)
        result, x_i = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        # LV freewall node outside coordinates and derivatives
        node = nodes.findNodeByIdentifier(2*nowl - 3)
        cache.setNode(node)
        result, x_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, dx_ds1_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, dx_ds2_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        # interpolate
        # GRC controllable xi location
        xi = 0.4
        xr = 1.0 - xi
        # GRC controllable derivative scaling
        sc = 2.0
        d_o = [ sc*v for v in dx_ds2_o ]
        d_c = [ sc*v for v in dx_ds1_c ]
        x = [ v for v in interpolateCubicHermite(x_o, d_o, x_c, d_c, xi) ]
        dx_ds1 = [ (xr*dx_ds1_o[i] + xi*dx_ds1_c[i]) for i in range(3) ]
        dx_ds2 = [ 0.5*v for v in interpolateCubicHermiteDerivative(x_o, d_o, x_c, d_c, xi) ]
        dx_ds3 = [ (x[i] - x_i[i]) for i in range(3) ]
        node = nodes.createNode(nodeIdentifier, nodetemplateFull)
        cache.setNode(node)
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        lv_nid1 = nodeIdentifier
        nodeIdentifier += 1


        #################
        # Create elements
        #################

        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        # LV Base elements

        nids = [
            [ nidl +  0, nidl +  1,        laNodeId[0][5],        laNodeId[0][0], nidl + nowl +  0, nidl + nowl +  1,       laNodeId[1][5],        laNodeId[1][0] ],
            [ nidl +  1, nidl +  2,        laNodeId[0][0],        raNodeId[1][5], nidl + nowl +  1, nidl + nowl +  2,       laNodeId[1][0],        raNodeId[0][5] ],
            [ nidl +  2, nidl +  3,        raNodeId[1][5],        raNodeId[1][4], nidl + nowl +  2, nidl + nowl +  3,       raNodeId[0][5],        raNodeId[0][4] ],
            [ nidl +  3, nidl +  4,  lvOutletNodeId[1][1],  lvOutletNodeId[1][2], nidl + nowl +  3, nidl + nowl +  4,       raNodeId[0][4],         midseptum_nid ],
            [ nidl +  4, nidl +  5,  lvOutletNodeId[1][2],  lvOutletNodeId[1][3], nidl + nowl +  4, nidl + nowl +  5,        midseptum_nid,  rvOutletNodeId[0][0] ],
            [ nidl +  5, nidl +  6,  rvOutletNodeId[1][0], rvOutletNodeId[1][-1], nidl + nowl +  5, nidl + nowl +  6, rvOutletNodeId[0][0], rvOutletNodeId[0][-1] ],
            # 6-node collapsed element on inside of RV-LV join. xi2=1 face collapsed to a line by merging edges adjacent to corner node 7
            [ nidl +  6, nidl +  7, rvOutletNodeId[1][-1],       nidl + nowl + 6, nidl + nowl +  7, rvOutletNodeId[0][-1] ],
            # 4 node collapsed / tetrahedral element
            [ nidl +  7, rvOutletNodeId[1][-1], nidl + nowl + 7, lv_nid1 ],
            # 7-node collapsed element
            [ nidl +  7, nidl +  8,        laNodeId[0][2], nidl + nowl +  7, nidl + nowl +  8,              lv_nid1,        laNodeId[1][2] ],
            [ nidl +  8, nidl +  9,        laNodeId[0][2],        laNodeId[0][3], nidl + nowl +  8, nidl + nowl +  9,       laNodeId[1][2],        laNodeId[1][3] ],
            [ nidl +  9, nidl + 10,        laNodeId[0][3],        laNodeId[0][4], nidl + nowl +  9, nidl + nowl + 10,       laNodeId[1][3],        laNodeId[1][4] ],
            [ nidl + 10, nidl +  0,        laNodeId[0][4],        laNodeId[0][5], nidl + nowl + 10, nidl + nowl +  0,       laNodeId[1][4],        laNodeId[1][5] ],
            # 5-node collapsed element at cusp of LV and RV outlets, on LV freewall end of septum
            [ nidl +  5, nidl +  6,  lvOutletNodeId[1][3],  lvOutletNodeId[1][4], rvOutletNodeId[1][-1] ],
            # 5-node square pyramid
            [ rvOutletNodeId[1][-1], nidl +  6, nidl +  7, lvOutletNodeId[1][4], lv_nid1 ],
            [ nidl +  6, nidl +  7,  laNodeId[0][1], laNodeId[0][2],lvOutletNodeId[1][-2], lv_nid1, laNodeId[1][1],  laNodeId[1][2] ]
        ]

        for e in range(len(nids)):
            eft1 = eft
            if e == 0:
                eft1 = tricubichermite.createEftSplitXi1RightStraight()
            elif e == 1:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == 3:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # must remap derivatives before remapping local nodes
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
            elif e == 4:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # must remap derivatives before remapping local nodes
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
            elif e == 5:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
                scaleEftNodeValueLabels(eft1, [ 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
            elif e == 6:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3, 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, [ ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
                tricubichermite.setEftLinearDerivativeXi1(eft1, 7, 8, 7, 8, 1)
                ln_map = [ 1, 2, 3, 3, 4, 5, 6, 3 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 7:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ])
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4, 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1,  [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 1, 1, 1, 2, 3, 2, 4 ]
                remapEftLocalNodes(eft1, 4, ln_map)
            elif e == 8:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS2, [ ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                ln_map = [ 1, 2, 1, 3,  4, 5, 6, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 12:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []) ])
                remapEftNodeValueLabel(eft1, [ 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [ ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
                ln_map = [ 1, 2, 3, 4, 2, 2, 5, 5 ]
                remapEftLocalNodes(eft1, 5, ln_map)
            elif e == 13:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ])
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                ln_map = [ 1, 1, 2, 3, 1, 1, 4, 5 ]
                remapEftLocalNodes(eft1, 5, ln_map)
            elif e == 14:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [ 1 ]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                # following uses D2_DS1DS3 as a temporary swap variable
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D2_DS1DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, [ 1 ]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                scaleEftNodeValueLabels(eft1, [ 7 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])

            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids[e])
            if eft1.getNumberOfLocalScaleFactors == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 1
            print('create element lv', elementIdentifier, result, result2, result3, nids[e])
            elementIdentifier += 1

        # RV base elements

        scalefactors5hanging = [ -1.0, 0.5, 0.25, 0.125, 0.75 ]
        nids = [
            [ nidl + nowl + 1, nidr + 0, raNodeId[0][0], raNodeId[0][1], nidl + nowl + 0, nidr + nowr + 0, laNodeId[1][5] , raNodeId[1][1] ],
            [ nidr        + 0, nidr + 1, raNodeId[0][1], raNodeId[0][2], nidr + nowr + 0, nidr + nowr + 1, raNodeId[1][1], raNodeId[1][2] ],
            [ nidr        + 1, nidr + 2, raNodeId[0][2], raNodeId[0][3], nidr + nowr + 1, nidr + nowr + 2, raNodeId[1][2], raNodeId[1][3] ],
            [ nidr        + 2, nidr + 3, raNodeId[0][3],              crest_nid1, nidr + nowr + 2, nidr + nowr + 3, raNodeId[1][3],              crest_nid2 ],
            # refined pair of elements with hanging nodes across xi1=0.5 at bottom
            [ nidr        + 3, nidr + 4,              crest_nid1,    rvOutletNodeId[0][1], nidr + nowr + 3, nidr + nowr + 4,           crest_nid2, rvOutletNodeId[1][1] ],
            [ nidr        + 3, nidr + 4,    rvOutletNodeId[0][1],    rvOutletNodeId[0][2], nidr + nowr + 3, nidr + nowr + 4, rvOutletNodeId[1][1], rvOutletNodeId[1][2] ],
            # refined pair of elements with hanging nodes across xi1=0.5 at bottom
            [ nidr        + 4, nidr + 5,    rvOutletNodeId[0][2],    rvOutletNodeId[0][3], nidr + nowr + 4, nidr + nowr + 5, rvOutletNodeId[1][2], rvOutletNodeId[1][3] ],
            [ nidr        + 4, nidr + 5,    rvOutletNodeId[0][3],    rvOutletNodeId[0][4], nidr + nowr + 4, nidr + nowr + 5, rvOutletNodeId[1][3], rvOutletNodeId[1][4] ],
            # return junction of RV with LV
            [ nidr        + 5, nidl + nowl + 6,    rvOutletNodeId[0][-2],   rvOutletNodeId[0][-1], nidr + nowr + 5, nidl + nowl + 7, rvOutletNodeId[1][-2], rvOutletNodeId[1][-1] ],
            # top of supraventricular crest
            [ raNodeId[0][3],           crest_nid1, raNodeId[0][4],        midseptum_nid, raNodeId[1][3],           crest_nid2, lvOutletNodeId[1][1], lvOutletNodeId[1][2] ],
            [     crest_nid1, rvOutletNodeId[0][1],  midseptum_nid, rvOutletNodeId[0][0],  crest_nid2, rvOutletNodeId[1][1], lvOutletNodeId[1][2], rvOutletNodeId[1][0] ]
        ]
        for e in range(len(nids)):
            eft1 = eft
            if e == 0:
                eft1 = tricubichermite.createEftSplitXi1RightOut()
            elif e == 3:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
            elif e in [ 4, 6 ]:
                eft1 = tricubichermite.createEftBasic()
                # general scale factors 1 -> 1, 102 -> 1/2, 104 -> 1/4, 108 -> 1/8, 304 -> 3/4
                setEftScaleFactorIds(eft1, [1, 102, 104, 108, 304], [])
                tricubichermite.setEftMidsideXi1HangingNode(eft1, 2, 1, 1, 2, [1, 2, 3, 4, 5])
                tricubichermite.setEftMidsideXi1HangingNode(eft1, 6, 5, 5, 6, [1, 2, 3, 4, 5])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                if e == 4:
                    remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                    # must do following after above otherwise rescaled
                    remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
                    # must do following after above otherwise rescaled
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                else:
                    tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            elif e in [ 5, 7 ]:
                eft1 = tricubichermite.createEftBasic()
                # general scale factors 1 -> 1, 102 -> 1/2, 104 -> 1/4, 108 -> 1/8, 304 -> 3/4
                setEftScaleFactorIds(eft1, [1, 102, 104, 108, 304], [])
                tricubichermite.setEftMidsideXi1HangingNode(eft1, 1, 2, 1, 2, [1, 2, 3, 4, 5])
                tricubichermite.setEftMidsideXi1HangingNode(eft1, 5, 6, 5, 6, [1, 2, 3, 4, 5])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                if e == 5:
                    remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == 8:
                # return junction of RV with LV
                # Same as tricubichermite.createEftSplitXi1RightIn() on bottom layer only
                eft1 = tricubichermite.createEftBasic()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == 9:
                # supraventricular crest 1, by RA-LV outlet junction
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                scaleEftNodeValueLabels(eft1, [ 3, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [1])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == 10:
                # supraventricular crest 2, by LV-RV outlet junction)
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                # must do following after above otherwise rescaled
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 2, 6, 2, 6, 1)
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])

            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids[e])
            if e in [ 4, 5, 6, 7 ]:
                result3 = element.setScaleFactors(eft1, scalefactors5hanging)
            elif eft1.getNumberOfLocalScaleFactors == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 1
            print('create element rv', elementIdentifier, result, result2, result3, nids[e])
            elementIdentifier += 1

        fm.endChange()

