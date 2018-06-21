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
            'Number of elements up septum' : 2,
            'Number of elements around' : 6,
            'Number of extra elements along septum' : 0,
            'Total arc up degrees' : 120.0,
            'Septum arc up ratio' : 0.5,
            'Length ratio' : 0.5,
            'Base septum thickness' : 0.1,
            'Free wall thickness' : 0.03,
            'Major axis rotation degrees' : 30.0,
            'Base inner major axis length' : 0.5,
            'Base inner minor axis length' : 0.4,
            'Vena cava inner diameter' : 0.12,
            'Vena cava wall thickness' : 0.012,
            'Pulmonary vein inner diameter' : 0.08,
            'Pulmonary vein wall thickness' : 0.008,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements up septum',
            'Number of elements around',
            'Number of extra elements along septum',
            'Total arc up degrees',
            'Septum arc up ratio',
            'Length ratio',
            'Base septum thickness',
            'Free wall thickness',
            'Major axis rotation degrees',
            'Base inner major axis length',
            'Base inner minor axis length',
            'Vena cava inner diameter',
            'Vena cava wall thickness',
            'Pulmonary vein inner diameter',
            'Pulmonary vein wall thickness'
            #,'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if options['Number of elements up'] < 3:
            options['Number of elements up'] = 3
        if options['Number of elements up septum'] > (options['Number of elements up'] - 2):
            options['Number of elements up septum'] = options['Number of elements up'] - 2
        if options['Number of elements up septum'] < 1:
            options['Number of elements up septum'] = 1
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        if options['Number of extra elements along septum'] > (options['Number of elements around'] - 4):
            options['Number of extra elements along septum'] = options['Number of elements around'] - 4
        if options['Number of extra elements along septum'] < 0:
            options['Number of extra elements along septum'] = 0
        if options['Total arc up degrees'] < 10.0:
            options['Total arc up degrees'] = 10.0
        elif options['Total arc up degrees'] > 170.0:
            options['Total arc up degrees'] = 170.0
        if options['Septum arc up ratio'] > 0.75:
            options['Septum arc up ratio'] = 0.75
        elif options['Septum arc up ratio'] < 0.25:
            options['Septum arc up ratio'] = 0.25
        for key in options:
            if key in [
                'Length ratio',
                'Base septum thickness',
                'Free wall thickness',
                'Base inner major axis length',
                'Base inner minor axis length',
                'Vena cava inner diameter',
                'Vena cava wall thickness',
                'Pulmonary vein inner diameter',
                'Pulmonary vein wall thickness']:
                if options[key] < 1.0E-6:
                    options[key] = 1.0E-6

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountUp = options['Number of elements up']
        elementsCountUpSeptum = options['Number of elements up septum']
        elementsCountAround = options['Number of elements around']
        #print('elementsCountAround', elementsCountAround)
        elementsCountAlongSeptum = options['Number of extra elements along septum']
        atriaAxisRadians = math.radians(options['Major axis rotation degrees'])
        totalArcUpRadians = math.radians(options['Total arc up degrees'])
        septumArcUpRadians = totalArcUpRadians*options['Septum arc up ratio']
        lengthRatio = options['Length ratio']
        baseSeptumThickness = options['Base septum thickness']
        freeWallThickness = options['Free wall thickness']
        majorAxisRadians = math.radians(options['Major axis rotation degrees'])
        baseToEquatorRatio = 1.0/math.sin(totalArcUpRadians)
        innerMajorMag = 0.5*options['Base inner major axis length']
        innerMinorMag = 0.5*options['Base inner minor axis length']
        vcInnerDiameter = options['Vena cava inner diameter']
        vcWallThickness = options['Vena cava wall thickness']
        pvInnerDiameter = options['Pulmonary vein inner diameter']
        pvWallThickness = options['Pulmonary vein wall thickness']
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
        laSeptumY = -0.5*baseSeptumThickness
        laCentreX = laSeptumX - laInnerMajorX*math.cos(laSeptumRadians) - laInnerMinorX*math.sin(laSeptumRadians)
        laCentreY = laSeptumY - laInnerMajorY*math.cos(laSeptumRadians) - laInnerMinorY*math.sin(laSeptumRadians)

        raSeptumX = 0.0
        raSeptumY = 0.5*baseSeptumThickness

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

        upRadians = [ [ 0.0 ]*(elementsCountUp +  1), [ 0.0 ]*(elementsCountUp +  1) ]
        deltaUpRadians = [ [ 0.0 ]*(elementsCountUp +  1), [ 0.0 ]*(elementsCountUp +  1) ]
        baseRadiansUp = math.pi - totalArcUpRadians
        deltaRadians1 = septumArcUpRadians/elementsCountUpSeptum
        deltaRadiansN = (totalArcUpRadians - septumArcUpRadians - 0.5*deltaRadians1) / (elementsCountUp - elementsCountUpSeptum - 0.5)
        for i in range(elementsCountUpSeptum + 1):
            upRadians[0][i] = baseRadiansUp + i*deltaRadians1
            deltaUpRadians[0][i] = deltaRadians1
        for i in range(elementsCountUp - elementsCountUpSeptum):
            upRadians[0][elementsCountUp - i] = math.pi - i*deltaRadiansN
            deltaUpRadians[0][elementsCountUp - i] = deltaRadiansN
        # GRC temp: should make separate range to get inlet incline at bottom:
        for i in range(elementsCountUp + 1):
            upRadians[1][i] = upRadians[0][i]
            deltaUpRadians[1][i] = deltaUpRadians[0][i]

        #print('upRadians',upRadians[0])
        #print('deltaUpRadians',deltaUpRadians[0])

        innerScaleZ = (lengthRatio - freeWallThickness)/(1.0 - math.cos(totalArcUpRadians))
        outerScaleZ = lengthRatio/(1.0 - math.cos(totalArcUpRadians))

        laNodeId = [ [], [] ]
        raNodeId = [ [], [] ]
        laApexNodeId = [ -1 ]*2
        raApexNodeId = [ -1 ]*2

        for n3 in range(2):
            for n2 in range(elementsCountUp):

                # GRC support separate inner and outer radians up
                radiansUp = upRadians[n3][n2]
                scalingUp = baseToEquatorRatio*math.sin(radiansUp)
                innerZ = -innerScaleZ*math.cos(radiansUp)
                outerZ = -outerScaleZ*math.cos(radiansUp)

                laLayerNodeId = [-1]*elementsCountAround
                laNodeId[n3].append(laLayerNodeId)
                raLayerNodeId = [-1]*elementsCountAround
                raNodeId[n3].append(raLayerNodeId)

                # regular nodes up atria
                for i in range(2):
                    if i == 0:
                        # left
                        centreX, centreY = laCentreX, laCentreY
                        aRadians = laRadians
                        aDeltaRadians = laDeltaRadians
                        layerNodeId = laLayerNodeId
                        innerMajorX, innerMajorY = laInnerMajorX, laInnerMajorY
                        innerMinorX, innerMinorY = laInnerMinorX, laInnerMinorY
                        outerMajorX, outerMajorY = laOuterMajorX, laOuterMajorY
                        outerMinorX, outerMinorY = laOuterMinorX, laOuterMinorY
                    else:
                        # right
                        centreX, centreY = raCentreX, raCentreY
                        aRadians = raRadians
                        aDeltaRadians = raDeltaRadians
                        layerNodeId = raLayerNodeId
                        innerMajorX, innerMajorY = raInnerMajorX, raInnerMajorY
                        innerMinorX, innerMinorY = raInnerMinorX, raInnerMinorY
                        outerMajorX, outerMajorY = raOuterMajorX, raOuterMajorY
                        outerMinorX, outerMinorY = raOuterMinorX, raOuterMinorY

                    for n1 in range(elementsCountAround):
                        radiansAround = aRadians[n1]
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        inner = [
                            centreX + scalingUp*(cosRadiansAround*innerMajorX + sinRadiansAround*innerMinorX),
                            centreY + scalingUp*(cosRadiansAround*innerMajorY + sinRadiansAround*innerMinorY),
                            innerZ ]
                        outer = [
                            centreX + scalingUp*(cosRadiansAround*outerMajorX + sinRadiansAround*outerMinorX),
                            centreY + scalingUp*(cosRadiansAround*outerMajorY + sinRadiansAround*outerMinorY),
                            outerZ ]
                        if (n3 == 1) and (n2 <= elementsCountUpSeptum) and (n1 == 0):
                            continue  # right septum node created in next loop
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        layerNodeId[n1] = nodeIdentifier
                        cache.setNode(node)
                        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                        if n3 == 0:
                            dx_ds1 = [
                                scalingUp*aDeltaRadians[n1]*(-sinRadiansAround*innerMajorX + cosRadiansAround*innerMinorX),
                                scalingUp*aDeltaRadians[n1]*(-sinRadiansAround*innerMajorY + cosRadiansAround*innerMinorY),
                                0.0 ]
                            dx_ds2 = [
                                deltaUpRadians[n3][n2]*math.cos(radiansUp)*(cosRadiansAround*innerMajorX + sinRadiansAround*innerMinorX),
                                deltaUpRadians[n3][n2]*math.cos(radiansUp)*(cosRadiansAround*innerMajorY + sinRadiansAround*innerMinorY),
                                deltaUpRadians[n3][n2]*math.sin(radiansUp)*innerScaleZ ]
                        else:
                            dx_ds1 = [
                                scalingUp*aDeltaRadians[n1]*(-sinRadiansAround*outerMajorX + cosRadiansAround*outerMinorX),
                                scalingUp*aDeltaRadians[n1]*(-sinRadiansAround*outerMajorY + cosRadiansAround*outerMinorY),
                                0.0 ]
                            dx_ds2 = [
                                deltaUpRadians[n3][n2]*math.cos(radiansUp)*(cosRadiansAround*outerMajorX + sinRadiansAround*outerMinorX),
                                deltaUpRadians[n3][n2]*math.cos(radiansUp)*(cosRadiansAround*outerMajorY + sinRadiansAround*outerMinorY),
                                deltaUpRadians[n3][n2]*math.sin(radiansUp)*outerScaleZ ]
                        dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                        nodeIdentifier += 1

            # apexes
            for i in range(2):
                if i == 0:
                    # left
                    centreX, centreY = laCentreX, laCentreY
                    apexNodeId = laApexNodeId
                    innerMajorX, innerMajorY = laInnerMajorX, laInnerMajorY
                    innerMinorX, innerMinorY = laInnerMinorX, laInnerMinorY
                    outerMajorX, outerMajorY = laOuterMajorX, laOuterMajorY
                    outerMinorX, outerMinorY = laOuterMinorX, laOuterMinorY
                else:
                    # right
                    centreX, centreY = raCentreX, raCentreY
                    aRadians = raRadians
                    apexNodeId = raApexNodeId
                    innerMajorX, innerMajorY = raInnerMajorX, raInnerMajorY
                    innerMinorX, innerMinorY = raInnerMinorX, raInnerMinorY
                    outerMajorX, outerMajorY = raOuterMajorX, raOuterMajorY
                    outerMinorX, outerMinorY = raOuterMinorX, raOuterMinorY

                x = [ centreX, centreY, innerScaleZ if (n3 == 0) else outerScaleZ ]
                if n3 == 0:
                    dx_ds1 = [
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*innerMajorX,
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*innerMajorY,
                        0.0 ]
                    dx_ds2 = [
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*innerMinorX,
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*innerMinorY,
                        0.0 ]
                else:
                    dx_ds1 = [
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*outerMajorX,
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*outerMajorY,
                        0.0 ]
                    dx_ds2 = [
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*outerMinorX,
                        baseToEquatorRatio*deltaUpRadians[n3][-1]*outerMinorY,
                        0.0 ]
                dx_ds3 = [ 0.0, 0.0, freeWallThickness ]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                apexNodeId[n3] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # transfer inner septum nodes to outer on opposite side, set derivative 3 to be node difference
        for n2 in range(elementsCountUpSeptum + 1):
            laNodeId[1][n2][0] = raNodeId[0][n2][0]
            raNodeId[1][n2][0] = laNodeId[0][n2][0]
            node2 = nodes.findNodeByIdentifier(laNodeId[1][n2][0])
            cache.setNode(node2)
            result, x_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            node1 = nodes.findNodeByIdentifier(laNodeId[0][n2][0])
            cache.setNode(node1)
            result, x_i = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            dx_ds3 = [ (x_o[i] - x_i[i]) for i in range(3) ]
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            dx_ds3 = [ -v for v in dx_ds3 ]
            cache.setNode(node2)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # create extra node(s) at top of septum
        n2 = elementsCountUpSeptum + 1
        node1 = nodes.findNodeByIdentifier(laNodeId[1][n2][0])
        cache.setNode(node1)
        result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        d1 = [ -d for d in d1 ]
        node2 = nodes.findNodeByIdentifier(raNodeId[1][n2][0])
        cache.setNode(node2)
        result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        xi = 0.5
        vc = interpolateCubicHermite(v1, d1, v2, d2, xi )
        dc = interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi )
        x = [ vc[0], vc[1], vc[2] ]
        # get magnitude of dx_ds1 from arc around apex to next node
        node = nodes.findNodeByIdentifier(apexNodeId[0])
        cache.setNode(node)
        result, ac = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        node = nodes.findNodeByIdentifier(laNodeId[1][n2 - 1][1])
        cache.setNode(node)
        result, vb = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        a = [ (vc[i] - ac[i]) for i in range(3) ]
        b = [ (vb[i] - ac[i]) for i in range(3) ]
        mag_a = math.sqrt(a[0]*a[0] + a[1]*a[1] + a[1]*a[1])
        mag_b = math.sqrt(b[0]*b[0] + b[1]*b[1] + b[1]*b[1])
        arcRadians = math.acos((a[0]*b[0] + a[1]*b[1] + a[1]*b[1]) / (mag_a*mag_b))
        mag = 0.5*(mag_a + mag_b)*arcRadians
        dx_ds1 = [ mag, 0.0, 0.0 ]
        dx_ds2 = [ 0.5*d for d in dc ]
        dx_ds3 = [ 0.0, 0.0, vc[2] + innerScaleZ*math.cos(math.pi - totalArcUpRadians + septumArcUpRadians) ]
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        septumNodeId = nodeIdentifier
        cache.setNode(node)
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        #print('septum', node.isValid(), result, ' nodes', laNodeId[1][n2][0], raNodeId[1][n2][0])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        nodeIdentifier += 1

        if False:
            # show centre/axes of atria
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, innerScaleZ*math.cos(totalArcUpRadians) ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laInnerMajorX, laInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laInnerMinorX, laInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, lengthRatio - freeWallThickness ])
            nodeIdentifier += 1

            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, innerScaleZ*math.cos(totalArcUpRadians) ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raInnerMajorX, raInnerMajorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raInnerMinorX, raInnerMinorY, 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, lengthRatio - freeWallThickness ])
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

        eftSeptumSplitRight = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftSeptumSplitRight, [1], [])
        remapEftNodeValueLabel(eftSeptumSplitRight, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftSeptumSplitRight, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
        #scaleEftNodeValueLabels(eftSeptumSplitRight, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
        elementtemplateSeptumSplitRight = mesh.createElementtemplate()
        elementtemplateSeptumSplitRight.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateSeptumSplitRight.defineField(coordinates, -1, eftSeptumSplitRight)

        eftSeptumSplitCentre = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftSeptumSplitCentre, [1], [])
        remapEftNodeValueLabel(eftSeptumSplitCentre, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eftSeptumSplitCentre, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eftSeptumSplitCentre, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
        elementtemplateSeptumSplitCentre = mesh.createElementtemplate()
        elementtemplateSeptumSplitCentre.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateSeptumSplitCentre.defineField(coordinates, -1, eftSeptumSplitCentre)

        # regular rows within septum
        for e2 in range(elementsCountUpSeptum):

            for i in range(2):
                if i == 0:  # left
                    nodeId = laNodeId
                    otherNodeId = raNodeId
                else:  # right
                    nodeId = raNodeId
                    otherNodeId = laNodeId

                for e1 in range(elementsCountAround):

                    en = (e1 + 1) % elementsCountAround
                    nids = [ nodeId[0][e2][e1], nodeId[0][e2][en], nodeId[0][e2 + 1][e1], nodeId[0][e2 + 1][en], \
                             nodeId[1][e2][e1], nodeId[1][e2][en], nodeId[1][e2 + 1][e1], nodeId[1][e2 + 1][en] ]

                    if e1 == 0:
                        eft1 = eftSeptumSplitRight
                        elementtemplate1 = elementtemplateSeptumSplitRight
                        nids[4] = otherNodeId[1][e2][-1]
                        nids[6] = otherNodeId[1][e2 + 1][-1]
                    elif e1 == (elementsCountAround - 1):
                        eft1 = eftSeptumSplitCentre
                        elementtemplate1 = elementtemplateSeptumSplitCentre
                    else:
                        eft1 = eft
                        elementtemplate1 = elementtemplate

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if eft1.getNumberOfLocalScaleFactors() == 1:
                        result3 = element.setScaleFactors(eft1, [ -1.0 ])
                    else:
                        result3 = 1
                    #print('create element', 'la' if i == 0 else 'ra', element.isValid(), elementIdentifier, result2, result3, nids)
                    elementIdentifier += 1

        # septum transition interior collapsed elements

        n2 = elementsCountUpSeptum
        nids = [
            [ laNodeId[0][n2][ 0], laNodeId[0][n2][1], raNodeId[1][n2][-1], laNodeId[1][n2][1], septumNodeId ],
            [ laNodeId[0][n2][-1], laNodeId[0][n2][0], laNodeId[1][n2][-1], laNodeId[1][n2][0], septumNodeId ],
            [ raNodeId[0][n2][ 0], raNodeId[0][n2][1], laNodeId[1][n2][-1], raNodeId[1][n2][1], septumNodeId ],
            [ raNodeId[0][n2][-1], raNodeId[0][n2][0], raNodeId[1][n2][-1], raNodeId[1][n2][0], septumNodeId ]
        ]

        for e in range(len(nids)):
            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            if e in [ 0, 2 ]:
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4, 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ] )
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                if e == 0:
                    #remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1]) ])
                    #remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                    #remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                ln_map = [ 1, 2, 1, 2, 3, 4, 5, 4 ]
                remapEftLocalNodes(eft1, 5, ln_map)
            else:
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4, 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ] )
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                if e == 1:
                    #remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 1, 2, 3, 4, 3, 5 ]
                remapEftLocalNodes(eft1, 5, ln_map)

            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)

            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids[e])
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 1
            #print('create element st', elementIdentifier, result, result2, result3, nids[e])
            elementIdentifier += 1

        # semi-regular rows above septum but below apex, including second septum transition elements

        for e2 in range(elementsCountUpSeptum, elementsCountUp - 1):

            for i in range(2):
                if i == 0:  # left
                    nodeId = laNodeId
                    otherNodeId = raNodeId
                else:  # right
                    nodeId = raNodeId
                    otherNodeId = laNodeId

                for e1 in range(elementsCountAround):

                    en = (e1 + 1) % elementsCountAround
                    nids = [ nodeId[0][e2][e1], nodeId[0][e2][en], nodeId[0][e2 + 1][e1], nodeId[0][e2 + 1][en], \
                             nodeId[1][e2][e1], nodeId[1][e2][en], nodeId[1][e2 + 1][e1], nodeId[1][e2 + 1][en] ]

                    if (e2 == elementsCountUpSeptum) and ((e1 == 0) or (e1 == elementsCountAround - 1)):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if e1 == 0:
                            nids[4] = septumNodeId
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1 if (i == 0) else 0]), (Node.VALUE_LABEL_D_DS2, [1 if (i == 0) else 0]) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1 if (i == 0) else 0]) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1 if (i == 0) else 0]), (Node.VALUE_LABEL_D_DS3, [0]) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                        else:
                            nids[5] = septumNodeId
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1 if (i == 0) else 0]), (Node.VALUE_LABEL_D_DS2, [0 if (i == 0) else 1]) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1 if (i == 0) else 0]) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1 if (i == 0) else 0]), (Node.VALUE_LABEL_D_DS3, [0]) ])
                        elementtemplate1 = mesh.createElementtemplate()
                        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        result = elementtemplate1.defineField(coordinates, -1, eft1)
                    else:
                        eft1 = eft
                        elementtemplate1 = elementtemplate

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if eft1.getNumberOfLocalScaleFactors() == 1:
                        result3 = element.setScaleFactors(eft1, [ -1.0 ])
                    else:
                        result3 = 1
                    #print('create element', 'la' if i == 0 else 'ra', element.isValid(), elementIdentifier, result2, result3, nids)
                    if e2 == elementsCountUp - 2:
                        if i == 0:
                            if e1 == 0:
                                rapvElementId = elementIdentifier
                            elif e1 == (elementsCountAround - 1):
                                rppvElementId = elementIdentifier
                            elif e1 == (elementsCountAround//2 - 1):
                                lapvElementId = elementIdentifier
                            elif e1 == ((elementsCountAround + 1)//2):
                                lppvElementId = elementIdentifier
                        else:  # i == 1:
                            if e1 == 1:
                                ivcElementId = elementIdentifier
                            elif e1 == (elementsCountAround - 2):
                                svcElementId = elementIdentifier

                    elementIdentifier += 1

        for i in range(2):
            nodeId = laNodeId if (i == 0) else raNodeId
            apexNodeId = laApexNodeId if (i == 0) else raApexNodeId
            aRadians = laRadians if (i == 0) else raRadians
            deltaRadians = laDeltaRadians if (i == 0) else raDeltaRadians

            # create top apex elements
            # scale factor identifiers follow convention of offsetting by 100 for each 'version'
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1)%elementsCountAround
                eft1 = tricubichermite.createEftShellApexTop(va*100, vb*100)
                elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                n2 = elementsCountUp - 1
                nodeIdentifiers = [ nodeId[0][n2][va], nodeId[0][n2][vb], apexNodeId[0], nodeId[1][n2][va], nodeId[1][n2][vb], apexNodeId[1] ]
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                # set general linear map coefficients
                scalefactors = [
                    -1.0,
                    -math.cos(aRadians[va]), -math.sin(aRadians[va]), deltaRadians[va],
                    -math.cos(aRadians[vb]), -math.sin(aRadians[vb]), deltaRadians[vb],
                    -math.cos(aRadians[va]), -math.sin(aRadians[va]), deltaRadians[va],
                    -math.cos(aRadians[vb]), -math.sin(aRadians[vb]), deltaRadians[vb]
                ]
                result = element.setScaleFactors(eft1, scalefactors)
                elementIdentifier = elementIdentifier + 1

        # add right atria inlets (venae cavae)

        for elementId in [ ivcElementId, svcElementId ]:
            element = mesh.findElementByIdentifier(elementId)
            vcLength = vcInnerDiameter*0.5
            tricubichermite.replaceElementWithInlet4(element, elementIdentifier, nodetemplate, nodeIdentifier, vcLength, vcInnerDiameter*0.5, vcWallThickness)
            elementIdentifier += 4
            nodeIdentifier += 8

        # add left atria inlets (pulmonary veins)

        for elementId in [ lapvElementId, lppvElementId, rapvElementId, rppvElementId ]:
            element = mesh.findElementByIdentifier(elementId)
            pvLength = pvInnerDiameter*0.5
            tricubichermite.replaceElementWithInlet4(element, elementIdentifier, nodetemplate, nodeIdentifier, pvLength, pvInnerDiameter*0.5, pvWallThickness)
            elementIdentifier += 4
            nodeIdentifier += 8

        fm.endChange()

