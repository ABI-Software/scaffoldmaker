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
        options['RV outlet wall thickness'] = 0.25
        options['RV outlet wall thickness'] = 0.02
        options['Outlet element length'] = 0.1
        options['Atria axis angle degrees'] = 30.0
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
            'Atria axis angle degrees'
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
        lvOutletRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        rvOutletRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        outletElementLength = options['Outlet element length']
        LVBaseFlattenRatio = options['LV base flatten ratio']
        LVBaseFlattenAngleRadians = options['LV base flatten angle degrees']*math.pi/180.0
        atriaAxisRadians = options['Atria axis angle degrees']*math.pi/180.0
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

        # node offsets for each row, and wall in LV
        norl = elementsCountAround
        nowl = 1 + elementsCountUp*norl

        radiansPerElementOrig = 2.0*math.pi/elementsCountAround
        septumCentreRadians = (1.0 + elementsCountAcrossSeptum/2.0)*radiansPerElementOrig
        lvDisplacementRadians = septumCentreRadians + LVBaseFlattenAngleRadians
        #lvDisplacementRadians = math.pi*2.0/3.0
        coslvDisplacementRadians = math.cos(lvDisplacementRadians)
        sinlvDisplacementRadians = math.sin(lvDisplacementRadians)

        outletJoinRadians = math.pi  # *0.8
        cosOutletJoinRadians = math.cos(outletJoinRadians)
        sinOutletJoinRadians = math.sin(outletJoinRadians)
        #lvOutletOffset = 0.5 - lvWallThickness*(1.0 - lvWallThicknessRatioBase) - lvOutletRadius - lvOutletWallThickness
        lvOutletOffset = 0.5 - lvWallThickness - lvOutletRadius
        #lvOutletCentreX = cosOutletJoinRadians*lvOutletOffset
        #lvOutletCentreY = sinOutletJoinRadians*lvOutletOffset
        lvOutletCentreX = lvOutletOffset*math.cos(lvDisplacementRadians)
        lvOutletCentreY = lvOutletOffset*math.sin(lvDisplacementRadians)

        # LV outlet - for bicubic-linear tube connection

        elementsCountAroundOutlet = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAroundOutlet
        lvOutletNodeId = [ [-1]*elementsCountAroundOutlet, [-1]*elementsCountAroundOutlet ]
        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        outletScale3 = 0.15
        for n3 in range(2):
            radius = lvOutletRadius + lvOutletWallThickness*n3
            for n1 in range(elementsCountAroundOutlet):
                radiansAround = outletJoinRadians - math.pi + n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = lvOutletCentreX + radius*cosRadiansAround
                x[1] = lvOutletCentreY + radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                nodetemplate = nodetemplateLinearS3 if ((n3 == 0) or (n1 == 3)) else nodetemplateFull
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                lvOutletNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if nodetemplate is nodetemplateFull:
                    dx_ds3[0] = outletScale3*cosRadiansAround
                    dx_ds3[1] = outletScale3*sinRadiansAround
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
        rvOutletCentreX = lvOutletCentreX + cosOutletJoinRadians*rvOutletOffset
        rvOutletCentreY = lvOutletCentreY + sinOutletJoinRadians*rvOutletOffset

        x = [ 0.0, 0.0, baseHeight + baseThickness ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, outletElementLength ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        for n3 in range(2):
            radius = rvOutletRadius + rvOutletWallThickness*n3
            for n1 in range(elementsCountAroundOutlet):
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
                rvOutletNodeId[n3][n1] = nodeIdentifier
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
            dx_ds2 = [ 0.0, 0.0, baseThickness ]
            result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            for c in range(3):
                dx_ds3[c] *= lvOutletWallThickness/outletScale3
            x[2] -= baseThickness
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            nodeIdentifier += 1
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        if False:
            # create extra nodes on LV 'crest' and edge of RA
            for i in range(2):
                scale1 = 0.08
                scale2 = 0.2
                x = [ -0.06 - scale1*i, -0.16, baseHeight ]
                dx_ds1 = [ scale1, 0.0, 0.0 ]
                dx_ds2 = [ 0.0, scale2, 0.0 ]
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

        # Septum nodes
        ##############

        if False:
            nid1 = nowl*2 - norl + 5
            nid2 = lvOutletNodeId[1][2]
            for n in range(3):
                node = nodes.findNodeByIdentifier(nid1 + n)
                cache.setNode(node)
                result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                result, dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                result, dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                node = nodes.findNodeByIdentifier(nid2 + n)
                cache.setNode(node)
                result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                dx_ds3 = [ x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2] ]
                cache.setNode(nodes.findNodeByIdentifier(nodePair[1]))
                if nodePair[0] is atriumLeftNodeId[0][1]:
                    dx_ds3= [ -x for x in dx_ds3 ]
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)


        # Atria nodes
        #############
        elementsCountAroundAtria = 6

        # left atrium interface nodes

        # GRC revisit:
        atrialSeptumThickness = lvWallThickness*lvWallThicknessRatioBase
        atrialSeptumLeftY = cruxCentre[1] - 0.5*atrialSeptumThickness
        atrialSeptumRightY = cruxCentre[1] + 0.5*atrialSeptumThickness

        atriumInletSlopeRadians = math.pi/6.0
        # GRC change from baseThickness?
        atriumInletSlopeLength = baseThickness*math.cos(atriumInletSlopeRadians)
        atriumInletSlopeHeight = baseThickness*math.sin(atriumInletSlopeRadians)

        if True:
            unitMinorAxisLx = -math.cos(atriaAxisRadians)
            unitMinorAxisLy = math.sin(atriaAxisRadians)
        elif False:
            mag = lvOutletRadius + lvOutletWallThickness
            unitMinorAxisLx = (lvOutletCentreX - cruxLeft[0])/mag
            unitMinorAxisLy = (lvOutletCentreY - cruxLeft[1])/mag
        else:
            radiansPerElementOrig = 2.0*math.pi/elementsCountAround
            septumCentreRadians = (1.0 + elementsCountAcrossSeptum/2.0)*radiansPerElementOrig
            minorAxisRadians = septumCentreRadians + LVBaseFlattenAngleRadians
            print('radians', septumCentreRadians, LVBaseFlattenAngleRadians, minorAxisRadians)
            unitMinorAxisLx = math.cos(minorAxisRadians)
            unitMinorAxisLy = math.sin(minorAxisRadians)
        # GRC minor length = lvOutletRadius
        minorAxisInnerMag = lvOutletRadius
        minorAxisInnerLx = unitMinorAxisLx*minorAxisInnerMag
        minorAxisInnerLy = unitMinorAxisLy*minorAxisInnerMag
        minorAxisOuterMag = minorAxisInnerMag + atriumInletSlopeLength
        minorAxisOuterLx = unitMinorAxisLx*minorAxisOuterMag
        minorAxisOuterLy = unitMinorAxisLy*minorAxisOuterMag
        print('minorAxisInner = ', minorAxisInnerLx, minorAxisInnerLy)

        atriumLeftCentre = [
            cruxLeft[0] - minorAxisOuterLx,
            cruxLeft[1] - minorAxisOuterLy,
            cruxLeft[2] - atriumInletSlopeHeight ]

        print('atriumLeftCentre = ', atriumLeftCentre)
        if True:
            ratio = minorAxisInnerLy/(atrialSeptumLeftY - atriumLeftCentre[1])
            atrialSeptumLeftRadians = math.asin(ratio)
            innerMajorMinorAxisRatio = -minorAxisInnerLy*math.cos(atrialSeptumLeftRadians)/(minorAxisInnerLx*math.sin(atrialSeptumLeftRadians))
            print('septum radians = ', atrialSeptumLeftRadians)
            majorAxisInnerLx = minorAxisInnerLy*innerMajorMinorAxisRatio
            majorAxisInnerLy = -minorAxisInnerLx*innerMajorMinorAxisRatio
            print('majorAxisInner = ', majorAxisInnerLx, majorAxisInnerLy, 'ratio = ', innerMajorMinorAxisRatio)
        else:
            # GRC 0.95 fudge factor for not being in centre of ellipse:
            majorAxisInnerMag = 0.95*(0.5 - lvWallThickness)*LVBaseFlattenRatio
            majorAxisInnerLx = unitMinorAxisLy*majorAxisInnerMag
            majorAxisInnerLy = -unitMinorAxisLx*majorAxisInnerMag
            innerMajorMinorAxisRatio = majorAxisInnerMag / minorAxisInnerMag
            # GRC temp:
            atrialSeptumLeftRadians = math.pi/4.0

        outerMajorMinorAxisRatio = (innerMajorMinorAxisRatio*minorAxisInnerMag + atriumInletSlopeLength)/minorAxisOuterMag
        majorAxisOuterLx = minorAxisOuterLy*outerMajorMinorAxisRatio
        majorAxisOuterLy = -minorAxisOuterLx*outerMajorMinorAxisRatio
        print('majorAxisOuter = ', majorAxisOuterLx, majorAxisOuterLy, 'ratio = ', outerMajorMinorAxisRatio)

        atriumLeftRadians = [0.0]*elementsCountAroundAtria
        atriumLeftRadians[0] = atrialSeptumLeftRadians
        atriumLeftRadians[1] = 0.5*math.pi
        deltaRadiansL1 = atriumLeftRadians[1] - atriumLeftRadians[0]
        deltaRadiansLN = (2.0*math.pi - 2.0*deltaRadiansL1)/(elementsCountAroundAtria - 2.0)
        atriumLeftDeltaRadians = [ deltaRadiansL1 ]*elementsCountAroundAtria
        atriumLeftRadians[2] = atriumLeftRadians[1] + 0.5*(deltaRadiansL1 + deltaRadiansLN)
        atriumLeftDeltaRadians[2] = deltaRadiansLN
        for i in range(3, elementsCountAroundAtria):
            atriumLeftRadians[i] = atriumLeftRadians[i - 1] + deltaRadiansLN
            atriumLeftDeltaRadians[i] = deltaRadiansLN

        atriumLeftNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        atriumLeftNodeId[1][1] = cruxLeftNodeId

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = atriumLeftRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    atriumLeftCentre[0] + cosRadiansAround*majorAxisInnerLx + sinRadiansAround*minorAxisInnerLx,
                    atriumLeftCentre[1] + cosRadiansAround*majorAxisInnerLy + sinRadiansAround*minorAxisInnerLy,
                    atriumLeftCentre[2] ]
                outer = [
                    atriumLeftCentre[0] + cosRadiansAround*majorAxisOuterLx + sinRadiansAround*minorAxisOuterLx,
                    atriumLeftCentre[1] + cosRadiansAround*majorAxisOuterLy + sinRadiansAround*minorAxisOuterLy,
                    atriumLeftCentre[2] + atriumInletSlopeHeight ]
                if (n3 == 1) and (n1 < 2):
                    continue  # already have a node from crux
                node = nodes.createNode(nodeIdentifier, nodetemplateFull)
                atriumLeftNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        atriumLeftDeltaRadians[n1]*(-sinRadiansAround*majorAxisInnerLx + cosRadiansAround*minorAxisInnerLx),
                        atriumLeftDeltaRadians[n1]*(-sinRadiansAround*majorAxisInnerLy + cosRadiansAround*minorAxisInnerLy),
                        0.0 ]
                else:
                    dx_ds1 = [
                        atriumLeftDeltaRadians[n1]*(-sinRadiansAround*majorAxisOuterLx + cosRadiansAround*minorAxisOuterLx),
                        atriumLeftDeltaRadians[n1]*(-sinRadiansAround*majorAxisOuterLy + cosRadiansAround*minorAxisOuterLy),
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
                    dx_ds3 = [ 0.0, atrialSeptumThickness, 0.0 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if True:
            # show axes of left atrium
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, atriumLeftCentre)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ majorAxisInnerLx, majorAxisInnerLy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ minorAxisInnerLx, minorAxisInnerLy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ atriumLeftCentre[0], atriumLeftCentre[1], atriumLeftCentre[2] + atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ majorAxisOuterLx, majorAxisOuterLy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ minorAxisOuterLx, minorAxisOuterLy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        # right atrium interface nodes

        minorAxisInnerRx = -minorAxisInnerLx
        minorAxisInnerRy = minorAxisInnerLy
        minorAxisOuterRx = -minorAxisOuterLx
        minorAxisOuterRy = minorAxisOuterLy
        atriumRightCentre = [
            cruxRight[0] + minorAxisOuterRx,
            cruxRight[1] + minorAxisOuterRy,
            cruxRight[2] - atriumInletSlopeHeight ]
        atrialSeptumRightRadians = math.pi*2.0 - atrialSeptumLeftRadians
        majorAxisInnerRx = majorAxisInnerLx
        majorAxisInnerRy = -majorAxisInnerLy
        majorAxisOuterRx = majorAxisOuterLx
        majorAxisOuterRy = -majorAxisOuterLy

        atriumRightRadians = [0.0]*elementsCountAroundAtria
        atriumRightRadians[0] = atrialSeptumRightRadians
        atriumRightRadians[elementsCountAroundAtria - 1] = 0.5*atrialSeptumRightRadians + 0.75*math.pi
        atriumRightRadians[elementsCountAroundAtria - 2] = 1.5*math.pi
        deltaRadiansR1 = 0.5*(atriumRightRadians[0] - atriumRightRadians[elementsCountAroundAtria - 2])
        deltaRadiansRN = (2.0*math.pi - 3.0*deltaRadiansR1)/(elementsCountAroundAtria - 3.0)
        atriumRightDeltaRadians = [ deltaRadiansR1 ]*elementsCountAroundAtria
        atriumRightRadians[1] = atrialSeptumRightRadians - math.pi*2.0 + 0.5*(deltaRadiansR1 + deltaRadiansRN)
        atriumRightDeltaRadians[1] = deltaRadiansRN
        for i in range(2, elementsCountAroundAtria - 2):
            atriumRightRadians[i] = atriumRightRadians[i - 1] + deltaRadiansRN
            atriumRightDeltaRadians[i] = deltaRadiansRN

        atriumRightNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        atriumRightNodeId[1][0] = atriumLeftNodeId[0][0]
        atriumRightNodeId[1][elementsCountAroundAtria - 2] = cruxRightNodeId
        atriumRightNodeId[1][elementsCountAroundAtria - 1] = cruxCentreNodeId

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = atriumRightRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    atriumRightCentre[0] + cosRadiansAround*majorAxisInnerRx + sinRadiansAround*minorAxisInnerRx,
                    atriumRightCentre[1] + cosRadiansAround*majorAxisInnerRy + sinRadiansAround*minorAxisInnerRy,
                    atriumRightCentre[2] ]
                outer = [
                    atriumRightCentre[0] + cosRadiansAround*majorAxisOuterRx + sinRadiansAround*minorAxisOuterRx,
                    atriumRightCentre[1] + cosRadiansAround*majorAxisOuterRy + sinRadiansAround*minorAxisOuterRy,
                    atriumRightCentre[2] + atriumInletSlopeHeight ]
                if atriumRightNodeId[n3][n1] >= 0:
                    continue  # already have a node from crux or left atrium interface
                node = nodes.createNode(nodeIdentifier, nodetemplateFull)
                atriumRightNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        atriumRightDeltaRadians[n1]*(-sinRadiansAround*majorAxisInnerRx + cosRadiansAround*minorAxisInnerRx),
                        atriumRightDeltaRadians[n1]*(-sinRadiansAround*majorAxisInnerRy + cosRadiansAround*minorAxisInnerRy),
                         0.0 ]
                else:
                    dx_ds1 = [
                        atriumRightDeltaRadians[n1]*(-sinRadiansAround*majorAxisOuterRx + cosRadiansAround*minorAxisOuterRx),
                        atriumRightDeltaRadians[n1]*(-sinRadiansAround*majorAxisOuterRy + cosRadiansAround*minorAxisOuterRy),
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
                    dx_ds3 = [ 0.0, atrialSeptumThickness, 0.0 ]
                elif n1 > (elementsCountAroundAtria - 3):
                    # negate two derivatives to align with left side across ventricular septum
                    dx_ds1 = [ -x for x in dx_ds1 ]
                    dx_ds3 = [ -x for x in dx_ds3 ]
                elif n1 == (elementsCountAroundAtria - 3):
                    dx_ds2 = dx_ds1
                    if n3 == 0:
                        mag = 0.25/math.sqrt(dx_ds3[0]*dx_ds3[0] + dx_ds3[1]*dx_ds3[1] + dx_ds3[2]*dx_ds3[2])
                        dx_ds1 = [ mag*dx_ds3[0], mag*dx_ds3[1], -mag*dx_ds3[2] ]
                    else:
                        mag = 0.25/math.sqrt(dx_ds3[0]*dx_ds3[0] + dx_ds3[1]*dx_ds3[1])
                        dx_ds1 = [ mag*dx_ds3[0], mag*dx_ds3[1], 0.0 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if True:
            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, atriumRightCentre)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ majorAxisInnerRx, majorAxisInnerRy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ minorAxisInnerRx, minorAxisInnerRy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplateFull)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ atriumRightCentre[0], atriumRightCentre[1], atriumRightCentre[2] + atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ majorAxisOuterRx, majorAxisOuterRy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ minorAxisOuterRx, minorAxisOuterRy, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        # transfer right atrium septum node to outside of left:
        atriumLeftNodeId[1][0] = atriumRightNodeId[0][0]

        # fix dx_ds3 on left and right crux nodes to be reverse or equal to neighbour on atrium interface
        for nodePair in [ (atriumLeftNodeId[0][1], atriumLeftNodeId[1][1]), (atriumRightNodeId[0][4], atriumRightNodeId[1][4]) ]:
            cache.setNode(nodes.findNodeByIdentifier(nodePair[0]))
            result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
            cache.setNode(nodes.findNodeByIdentifier(nodePair[1]))
            if nodePair[0] is atriumLeftNodeId[0][1]:
                dx_ds3= [ -x for x in dx_ds3 ]
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        #################
        # Create elements
        #################

        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        # LV Base

        nid1 = nowl - norl + 1
        nid3 = 43
        nids = [
            [ nid1 + 0, nid1 + 1,  atriumLeftNodeId[0][5],  atriumLeftNodeId[0][0], nid1 + nowl + 0, nid1 + nowl + 1,  atriumLeftNodeId[1][5],  atriumLeftNodeId[1][0] ],
            [ nid1 + 1, nid1 + 2,  atriumLeftNodeId[0][0], atriumRightNodeId[1][5], nid1 + nowl + 1, nid1 + nowl + 2,  atriumLeftNodeId[1][0], atriumRightNodeId[0][5] ],
            [ nid1 + 2, nid1 + 3, atriumRightNodeId[1][5], atriumRightNodeId[1][4], nid1 + nowl + 2, nid1 + nowl + 3, atriumRightNodeId[0][5], atriumRightNodeId[0][4] ],
            [ nid1 + 3, nid1 + 4,    lvOutletNodeId[1][1],    lvOutletNodeId[1][2], nid1 + nowl + 3, nid1 + nowl + 4, atriumRightNodeId[0][4] ],
            [ nid1 + 4, nid1 + 5,    lvOutletNodeId[1][2],    lvOutletNodeId[1][3], nid1 + nowl + 4, nid1 + nowl + 5, rvOutletNodeId[0][0] ],
            [ nid1 + 5, nid1 + 6,    lvOutletNodeId[1][3],    lvOutletNodeId[1][4], nid1 + nowl + 5, nid1 + nowl + 6, rvOutletNodeId[0][0] ],
            [ nid3 + 0, nid3 + 1,  atriumLeftNodeId[0][2],  atriumLeftNodeId[0][3], nid3 + nowl + 0, nid3 + nowl + 1,  atriumLeftNodeId[1][2],  atriumLeftNodeId[1][3] ],
            [ nid3 + 1, nid3 + 2,  atriumLeftNodeId[0][3],  atriumLeftNodeId[0][4], nid3 + nowl + 1, nid3 + nowl + 2,  atriumLeftNodeId[1][3],  atriumLeftNodeId[1][4] ],
            [ nid3 + 2, nid1 + 0,  atriumLeftNodeId[0][4],  atriumLeftNodeId[0][5], nid3 + nowl + 2, nid1 + nowl + 0,  atriumLeftNodeId[1][4],  atriumLeftNodeId[1][5] ]
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
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, [])
                #remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                # collapsed in xi2 on corner xi1 = 1, xi3 = 1
                ln_map = [ 1, 2, 3, 4, 5, 6, 7, 6 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e in [4, 5]:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                # add 9th node to get 'dx/ds3' by difference of x at two nodes
                eft1.setNumberOfLocalNodes(9)
                # must remap derivatives before remapping local nodes
                diffln1 = 3 if (e == 4) else 4
                remapEftNodeValueLabel(eft1, [ diffln1 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                diffln2 = 4 if (e == 4) else 3
                remapEftNodeValueLabelWithNodes(eft1, diffln2, Node.VALUE_LABEL_D_DS3, [ (diffln2, Node.VALUE_LABEL_D_DS2, [1]), (diffln2, Node.VALUE_LABEL_VALUE, [1]), (9, Node.VALUE_LABEL_VALUE, []) ])
                remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                # collapsed in xi2 on face xi3 = 1
                ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)

            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids[e])
            print('create element new', elementIdentifier, result, result2, nids[e])
            elementIdentifier += 1

        # RV
        nidr = 101
        nowr = 18
        nids = [
            [ nid1 + nowl + 1, nidr + 0, atriumRightNodeId[0][0], atriumRightNodeId[0][1], nid1 + nowl + 0, nidr + nowr + 0, atriumLeftNodeId[1][5] , atriumRightNodeId[1][1] ],
            [ nidr        + 0, nidr + 1, atriumRightNodeId[0][1], atriumRightNodeId[0][2], nidr + nowr + 0, nidr + nowr + 1, atriumRightNodeId[1][1], atriumRightNodeId[1][2] ],
            [ nidr        + 1, nidr + 2, atriumRightNodeId[0][2], atriumRightNodeId[0][3], nidr + nowr + 1, nidr + nowr + 2, atriumRightNodeId[1][2], atriumRightNodeId[1][3] ],
            [ nidr        + 2, nidr + 3, atriumRightNodeId[0][3],              crest_nid1, nidr + nowr + 2, nidr + nowr + 3, atriumRightNodeId[1][3],              crest_nid2 ],
            [ atriumRightNodeId[0][3], crest_nid1, atriumRightNodeId[0][4], 2*nowl - norl + 5, atriumRightNodeId[1][3], crest_nid2, lvOutletNodeId[1][1], lvOutletNodeId[1][2] ]
        ]
        for e in range(len(nids)):
            eft1 = eft
            if e == 0:
                eft1 = tricubichermite.createEftSplitXi1RightOut()
            elif e == 2:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
            elif e == 3:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
            elif e == 4:
                # supraventricular crest 1
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                scaleEftNodeValueLabels(eft1, [ 3, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [1])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids[e])
            print('create element new rv', elementIdentifier, result, result2, nids[e])
            elementIdentifier += 1


        if False:
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
                eft1 = tricubichermite.createEftBasic()
                eft1.setNumberOfLocalNodes(6)
                setEftScaleFactorIds(eft1, [1], [])
                for ln in [5, 6]:
                    n = ln - 1
                    # general map d/ds2 to get angle:
                    mapEftFunction1Node2Terms(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [], Node.VALUE_LABEL_D_DS3, 1, [1])
                for ln in [3, 4]:
                    n = ln + 3
                    # shape value and d/ds1
                    eft1.setTermNodeParameter(n*8 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                    eft1.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    # general map d/ds2 to get angle:
                    mapEftFunction1Node2Terms(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [], Node.VALUE_LABEL_D_DS3, 1, [1])
                    # zero d/dxi3
                    eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
                print('eft1.validate()', eft1.validate())
                elementtemplateSeptumWedge = mesh.createElementtemplate()
                elementtemplateSeptumWedge.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplateSeptumWedge.defineField(coordinates, -1, eft1)

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

            # PO segment on RV free wall h element
            eft1 = tricubichermite.createEftBasic()
            setEftScaleFactorIds(eft1, [1], [])
            for ln in [2]:
                n = ln - 1
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            for ln in [ 6]:
                n = ln - 1
                #mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 8, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 114, 93, 149, 150, 132, 94, 154, 155 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element PO 5', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            # Last PO segment on LV/AO side
            eft1 = tricubichermite.createEftBasic()
            ln_map = [ 1, 2, 3, 4, 5, 6, 7, 6 ]
            remapEftLocalNodes(eft1, 7, ln_map)
            setEftScaleFactorIds(eft1, [1], [])
            for n in [0, 1, 4, 5]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            for n in [0, 1, 4, 5, 7]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
            for n in [4]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
            for n in [5, 7]:
                ln = ln_map[n]
                eft1.setFunctionNumberOfTerms(n*8 + 3, 0)
            for n in [5]:
                ln = ln_map[n]
                eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 6, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 93, 92, 150, 145, 143, 142, 155 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element PO 6', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            # LV/RV junction 1: element collapsed on 3 edges
            eft1 = tricubichermite.createEftBasic()
            ln_map = [ 1, 1, 2, 2, 3, 4, 5, 5 ]
            remapEftLocalNodes(eft1, 5, ln_map)
            setEftScaleFactorIds(eft1, [1, 2], [])
            for n in [0, 1, 2, 3, 6, 7]:
                ln = ln_map[n]
                eft1.setFunctionNumberOfTerms(n*8 + 2, 0)
            for n in [0]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            for n in [1]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            for n in [4]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [])
                #mapEftFunction1Node2Terms(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [], Node.VALUE_LABEL_D_DS2, 1, [])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
                # cross derivative to reduce bulge:
                mapEftFunction1Node2Terms(eft1, n*8 + 4, ln, Node.VALUE_LABEL_D_DS1, 1, [2], Node.VALUE_LABEL_D_DS2, 1, [1, 2])
            for n in [5]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 2, 5, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 2, 5, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 93, 150, 94, 143, 155 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0, 4.0])
            print('create element LV/RV junction 1', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            # LV/RV junction 2: element collapsed on 1 edge
            # This needs work - non positive jacobian
            eft1 = tricubichermite.createEftBasic()
            ln_map = [ 1, 2, 3, 4, 5, 6, 7, 7 ]
            remapEftLocalNodes(eft1, 7, ln_map)
            setEftScaleFactorIds(eft1, [1], [])
            for n in [2]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
            for n in [3]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [])
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            for n in [6, 7]:
                ln = ln_map[n]
                eft1.setFunctionNumberOfTerms(n*8 + 2, 0)
            for n in [6]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            for n in [7]:
                ln = ln_map[n]
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 7, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 44, 45, 137, 160, 93, 94, 143 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element LV/RV junction 2', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            # AO 5
            eft1 = tricubichermite.createEftBasic()
            ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
            remapEftLocalNodes(eft1, 6, ln_map)
            setEftScaleFactorIds(eft1, [1], [])
            for n in [0, 1]:
                ln = ln_map[n]
                eft1.setFunctionNumberOfTerms(n*8 + 5, 0)
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            tricubichermite.setEftLinearDerivativeXi3(eft1, 3, 7, 3, 5, 1)
            tricubichermite.setEftLinearDerivativeXi3(eft1, 4, 8, 4, 6, 1)
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 160, 161, 137, 138, 143, 144 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element LV/RV junction 2', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1

            if False:
                # LV supraventricular crest at LA by LV free wall
                eft1 = tricubichermite.createEftBasic()
                setEftScaleFactorIds(eft1, [1], [])
                for n in [1, 5]:
                    ln = n + 1
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                nodeIdentifiers1 = [ 46, 47, 164, 162, 95, 96, 165, 163 ]
                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
                result3 = element.setScaleFactors(eft1, [-1.0])
                print('create element supra LA 1', elementIdentifier, result2, result3, nodeIdentifiers1 )
                elementIdentifier += 1

                # LV supraventricular crest at LA by LV septum 2
                eft1 = tricubichermite.createEftBasic()
                setEftScaleFactorIds(eft1, [1], [])
                for n in [2, 3, 6, 7]:
                    ln = n + 1
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                    mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
                for n in [2, 6]:
                    ln = n + 1
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                for n in [3, 7]:
                    ln = n + 1
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                nodeIdentifiers1 = [ 164, 162, 160, 161, 165, 163, 143, 144 ]
                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
                result3 = element.setScaleFactors(eft1, [-1.0])
                print('create element supra LA 2', elementIdentifier, result2, result3, nodeIdentifiers1 )
                elementIdentifier += 1

                # LV supraventricular crest at LA by LV septum 3
                eft1 = tricubichermite.createEftBasic()
                setEftScaleFactorIds(eft1, [1], [])
                for n in [2, 6]:
                    ln = n + 1
                    #mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                    mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                    mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
                for n in [3, 7]:
                    ln = n + 1
                    mapEftFunction1Node1Term(eft1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                nodeIdentifiers1 = [ 45, 46, 160, 164, 94, 95, 143, 165 ]
                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
                result3 = element.setScaleFactors(eft1, [-1.0])
                print('create element supra LA 3', elementIdentifier, result2, result3, nodeIdentifiers1 )
                elementIdentifier += 1

            # LV free wall to AO
            eft1 = tricubichermite.createEftBasic()
            setEftScaleFactorIds(eft1, [1], [])
            for n in [1, 5]:
                ln = n + 1
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
            for n in [2, 3, 6, 7]:
                ln = n + 1
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [])
            for n in [3, 7]:
                ln = n + 1
                mapEftFunction1Node1Term(eft1, n*8 + 3, ln, Node.VALUE_LABEL_D_DS1, 1, [])
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            nodeIdentifiers1 = [ 45, 46, 160, 161, 94, 95, 143, 144 ]
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers1)
            result3 = element.setScaleFactors(eft1, [-1.0])
            print('create element LV free wall to AO', elementIdentifier, result2, result3, nodeIdentifiers1 )
            elementIdentifier += 1



        fm.endChange()

