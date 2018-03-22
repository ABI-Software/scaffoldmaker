"""
Generates a 3-D heart atria model, suitable for attachment to the
3-D Heart Ventricles with Base 1.
"""

from __future__ import division
import math
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartatria1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Atria 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements up' : 4,
            'Number of elements around' : 6,
            'Number of extra elements along septum' : 0,
            'Total arc up degrees' : 120.0,
            'Length ratio' : 1.0,
            'Septum thickness' : 0.1,
            'Free wall thickness' : 0.03,
            'Major axis rotation degrees' : 30.0,
            'Base inner major axis length' : 0.5,
            'Base inner minor axis length' : 0.35,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements around',
            'Number of extra elements along septum',
            'Total arc up degrees',
            'Length ratio',
            'Septum thickness',
            'Free wall thickness',
            'Major axis rotation degrees',
            'Base inner major axis length',
            'Base inner minor axis length',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        if options['Number of extra elements along septum'] > (options['Number of elements around'] - 4):
            options['Number of extra elements along septum'] = options['Number of elements around'] - 4
        if options['Number of extra elements along septum'] < 0:
            options['Number of extra elements along septum'] = 0
        for key in options:
            if key in [
                'Length ratio',
                'Septum thickness',
                'Free wall thickness',
                'Base inner major axis length',
                'Base inner minor axis length']:
                if options[key] < 1.0E-6:
                    options[key] = 1.0E-6
        if options['Total arc up degrees'] < 10.0:
            options['Total arc up degrees'] = 10.0
        elif options['Total arc up degrees'] > 170.0:
            options['Total arc up degrees'] = 170.0

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountUp = options['Number of elements up']
        elementsCountAround = options['Number of elements around']
        #print('elementsCountAround', elementsCountAround)
        elementsCountAlongSeptum = options['Number of extra elements along septum']
        atriaAxisRadians = math.radians(options['Major axis rotation degrees'])
        totalArcUpRadians = math.radians(options['Total arc up degrees'])
        lengthRatio = options['Length ratio']
        septumThickness = options['Septum thickness']
        freeWallThickness = options['Free wall thickness']
        majorAxisRadians = math.radians(options['Major axis rotation degrees'])
        baseToEquatorRatio = 1.0/math.sin(totalArcUpRadians)
        innerMajorMag = 0.5*options['Base inner major axis length']
        innerMinorMag = 0.5*options['Base inner minor axis length']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)

        cache = fm.createFieldcache()

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

        ##############
        # Create nodes
        ##############

        nodeIdentifier = 1

        outerMajorMag = innerMajorMag + freeWallThickness
        outerMinorMag = innerMinorMag + freeWallThickness

        laInnerMajorX = innerMajorMag*math.sin(majorAxisRadians)
        laInnerMajorY = innerMajorMag*math.cos(majorAxisRadians)
        laInnerMinorX = -innerMinorMag*math.cos(majorAxisRadians)
        laInnerMinorY = innerMinorMag*math.sin(majorAxisRadians)

        laOuterMajorX = outerMajorMag*math.sin(majorAxisRadians)
        laOuterMajorY = outerMajorMag*math.cos(majorAxisRadians)
        laOuterMinorX = -outerMinorMag*math.cos(majorAxisRadians)
        laOuterMinorY = outerMinorMag*math.sin(majorAxisRadians)

        # following is radians around LA
        laSeptumRadians = math.atan2(laInnerMinorY, laInnerMajorY)
        laSeptumX = 0.0
        laSeptumY = -0.5*septumThickness
        laCentreX = laSeptumX - laInnerMajorX*math.cos(laSeptumRadians) - laInnerMinorX*math.sin(laSeptumRadians)
        laCentreY = laSeptumY - laInnerMajorY*math.cos(laSeptumRadians) - laInnerMinorY*math.sin(laSeptumRadians)

        raSeptumX = 0.0
        raSeptumY = 0.5*septumThickness

        raSeptumRadians = math.pi*2.0 - laSeptumRadians
        raCentreX = laCentreX
        raCentreY = -laCentreY

        raInnerMajorX = innerMajorMag*math.sin(majorAxisRadians)
        raInnerMajorY = -innerMajorMag*math.cos(majorAxisRadians)
        raInnerMinorX = innerMinorMag*math.cos(majorAxisRadians)
        raInnerMinorY = innerMinorMag*math.sin(majorAxisRadians)

        raOuterMajorX = outerMajorMag*math.sin(majorAxisRadians)
        raOuterMajorY = -outerMajorMag*math.cos(majorAxisRadians)
        raOuterMinorX = outerMinorMag*math.cos(majorAxisRadians)
        raOuterMinorY = outerMinorMag*math.sin(majorAxisRadians)

        laRadians = [0.0]*elementsCountAround
        laRadians[0] = laSeptumRadians
        deltaRadians = 2.0*math.pi/elementsCountAround
        laDeltaRadians = [ deltaRadians ]*elementsCountAround
        for i in range(1, elementsCountAround):
            laRadians[i] = laRadians[i - 1] + deltaRadians

        raRadians = [0.0]*elementsCountAround
        raRadians[0] = raSeptumRadians
        raDeltaRadians = [ deltaRadians ]*elementsCountAround
        for i in range(1, elementsCountAround):
            raRadians[i] = raRadians[i - 1] + deltaRadians

        baseRadiansUp = math.pi - totalArcUpRadians
        radiansPerElementUp = totalArcUpRadians/elementsCountUp
        scaleZ = lengthRatio/(1.0 - math.cos(totalArcUpRadians))

        laNodeId = [ [], [] ]
        raNodeId = [ [], [] ]

        for n3 in range(2):
            for n2 in range(elementsCountUp):

                radiansUp = baseRadiansUp + radiansPerElementUp*n2
                scalingUp = baseToEquatorRatio*math.sin(radiansUp)
                innerZ = -scaleZ*math.cos(radiansUp)

                laLayerNodeId = [-1]*elementsCountAround
                laNodeId[n3].append(laLayerNodeId)
                raLayerNodeId = [-1]*elementsCountAround
                raNodeId[n3].append(raLayerNodeId)

                # left atrium
                for n1 in range(elementsCountAround):
                    radiansAround = laRadians[n1]
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    inner = [
                        laCentreX + scalingUp*(cosRadiansAround*laInnerMajorX + sinRadiansAround*laInnerMinorX),
                        laCentreY + scalingUp*(cosRadiansAround*laInnerMajorY + sinRadiansAround*laInnerMinorY),
                        innerZ ]
                    outer = [
                        laCentreX + scalingUp*(cosRadiansAround*laOuterMajorX + sinRadiansAround*laOuterMinorX),
                        laCentreY + scalingUp*(cosRadiansAround*laOuterMajorY + sinRadiansAround*laOuterMinorY),
                        innerZ ]  # error
                    if (n3 == 1) and (n1 == 0):
                        continue  # right septum node created in next loop
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    laLayerNodeId[n1] = nodeIdentifier
                    cache.setNode(node)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                    if n3 == 0:
                        dx_ds1 = [
                            scalingUp*laDeltaRadians[n1]*(-sinRadiansAround*laInnerMajorX + cosRadiansAround*laInnerMinorX),
                            scalingUp*laDeltaRadians[n1]*(-sinRadiansAround*laInnerMajorY + cosRadiansAround*laInnerMinorY),
                            0.0 ]
                    else:
                        dx_ds1 = [
                            scalingUp*laDeltaRadians[n1]*(-sinRadiansAround*laOuterMajorX + cosRadiansAround*laOuterMinorX),
                            scalingUp*laDeltaRadians[n1]*(-sinRadiansAround*laOuterMajorY + cosRadiansAround*laOuterMinorY),
                            0.0 ]
                    dx_ds2 = [ 0.0, 0.0, 0.1 ]
                    dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                    if n1 == 0:
                        dx_ds3 = [ raSeptumX - laSeptumX, raSeptumY - laSeptumY, 0.0 ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier += 1

                # right atrium
                for n1 in range(elementsCountAround):
                    radiansAround = raRadians[n1]
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    inner = [
                        raCentreX + scalingUp*(cosRadiansAround*raInnerMajorX + sinRadiansAround*raInnerMinorX),
                        raCentreY + scalingUp*(cosRadiansAround*raInnerMajorY + sinRadiansAround*raInnerMinorY),
                        innerZ ]
                    outer = [
                        raCentreX + scalingUp*(cosRadiansAround*raOuterMajorX + sinRadiansAround*raOuterMinorX),
                        raCentreY + scalingUp*(cosRadiansAround*raOuterMajorY + sinRadiansAround*raOuterMinorY),
                        innerZ ]  # error
                    if (n3 == 1) and (n1 == 0):
                        continue  # right septum node created in previous loop
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    raLayerNodeId[n1] = nodeIdentifier
                    cache.setNode(node)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                    if n3 == 0:
                        dx_ds1 = [
                            scalingUp*raDeltaRadians[n1]*(-sinRadiansAround*raInnerMajorX + cosRadiansAround*raInnerMinorX),
                            scalingUp*raDeltaRadians[n1]*(-sinRadiansAround*raInnerMajorY + cosRadiansAround*raInnerMinorY),
                            0.0 ]
                    else:
                        dx_ds1 = [
                            scalingUp*raDeltaRadians[n1]*(-sinRadiansAround*raOuterMajorX + cosRadiansAround*raOuterMinorX),
                            scalingUp*raDeltaRadians[n1]*(-sinRadiansAround*raOuterMajorY + cosRadiansAround*raOuterMinorY),
                            0.0 ]
                    dx_ds2 = [ 0.0, 0.0, 0.1 ]
                    dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                    if n1 == 0:
                        dx_ds3 = [ raSeptumX - laSeptumX, raSeptumY - laSeptumY, 0.0 ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier += 1

        # fix up septum outer nodes:
        for n2 in range(elementsCountUp):
            laNodeId[1][n2][0] = raNodeId[0][n2][0]
            raNodeId[1][n2][0] = laNodeId[0][n2][0]

        if True:
            # show centre/axes of atria
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, -scaleZ*math.cos(totalArcUpRadians) ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laInnerMajorX, laInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laInnerMinorX, laInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, 0.1 ])
            nodeIdentifier += 1

            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, -scaleZ*math.cos(totalArcUpRadians) ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raInnerMajorX, raInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raInnerMinorX, raInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, 0.1 ])
            nodeIdentifier += 1

        #################
        # Create elements
        #################

        elementIdentifier = 1

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        eftRightOut = tricubichermite.createEftSplitXi1RightOut()
        elementtemplateRightOut = mesh.createElementtemplate()
        elementtemplateRightOut.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateRightOut.defineField(coordinates, -1, eftRightOut)

        eftRightStraight = tricubichermite.createEftSplitXi1RightStraight()
        elementtemplateRightStraight = mesh.createElementtemplate()
        elementtemplateRightStraight.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateRightStraight.defineField(coordinates, -1, eftRightStraight)

        for e2 in range(elementsCountUp - 1):

            # left atrium elements
            for e1 in range(elementsCountAround):
                if e1 == 0:
                    eft1 = eftRightOut
                    elementtemplate1 = elementtemplateRightOut
                elif e1 == (elementsCountAround - 1):
                    eft1 = eftRightStraight
                    elementtemplate1 = elementtemplateRightStraight
                else:
                    eft1 = eft
                    elementtemplate1 = elementtemplate

                en = (e1 + 1) % elementsCountAround
                nids = [ laNodeId[0][e2][e1], laNodeId[0][e2][en], laNodeId[0][e2 + 1][e1], laNodeId[0][e2 + 1][en], \
                         laNodeId[1][e2][e1], laNodeId[1][e2][en], laNodeId[1][e2 + 1][e1], laNodeId[1][e2 + 1][en] ]

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = 1
                print('create element la', elementIdentifier, result, result2, result3, nids)
                elementIdentifier += 1

        fm.endChange()

