"""
Generates 3-D Left and Right ventricles mesh starting from modified sphere shell mesh.
"""

from __future__ import division
import math
import scaffoldmaker.utils.vector as vector
from scaffoldmaker.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_heartventricles2:
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Ventricles 2'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around LV free wall' : 5,
            'Number of elements around septum' : 5,
            'Number of elements up' : 4,
            'Number of elements up septum' : 3,
            'Total height' : 1.0,
            'LV outer radius' : 0.5,
            'LV free wall thickness' : 0.12,
            'LV apex thickness' : 0.06,
            'RV height' : 0.7,
            'RV arc around degrees' : 160.0,
            'RV free wall thickness' : 0.04,
            'RV width' : 0.4,
            'RV extra cross radius base' : 0.05,
            'Septum thickness' : 0.08,
            'Septum base radial displacement' : 0.1,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements surface' : 1,
            'Refine number of elements through LV wall' : 1,
            'Refine number of elements through RV wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around LV free wall',
            'Number of elements around septum',
            'Number of elements up',
            'Number of elements up septum',
            'Total height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV height',
            'RV arc around degrees',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Septum thickness',
            'Septum base radial displacement',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall',
            'Number of elements up septum']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around LV free wall',
            'Number of elements around septum',
            'Number of elements up']:
            if options[key] < 2:
                options[key] = 2
        if options['Number of elements up septum'] >= options['Number of elements up']:
            options['Number of elements up septum'] = options['Number of elements up'] - 1
        for key in [
            'Total height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV height',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Septum thickness',
            'Septum base radial displacement']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['RV height'] > 0.99*options['Total height']:
            options['RV height'] = 0.99*options['Total height']
        if options['RV arc around degrees'] < 45.0:
            options['RV arc around degrees'] = 45.0
        elif options['RV arc around degrees'] > 270.0:
            options['RV arc around degrees'] = 270.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundSeptum = options['Number of elements around septum']
        elementsCountUp = options['Number of elements up']
        elementsCountUpSeptum = options['Number of elements up septum']
        totalHeight = options['Total height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        lvApexThickness = options['LV apex thickness']
        rvHeight = options['RV height']
        rvArcAroundRadians = math.radians(options['RV arc around degrees'])
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        rvExtraCrossRadiusBase = options['RV extra cross radius base']
        septumThickness = options['Septum thickness']
        septumBaseRadialDisplacement = options['Septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        #################
        # Create nodes
        #################

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplateApex = nodetemplate

        nodeIdentifier = 1

        # LV nodes

        # get element spacing up
        # approximate correction to move up to next inner node, since bottom of RV is halfway up element
        rvHeightCorrection = (elementsCountUpSeptum - 1)/(elementsCountUpSeptum - 0.5)
        rvThetaUp = math.acos(rvHeight*rvHeightCorrection/totalHeight)
        rvOuterHeight = 0.5*(rvHeight + totalHeight)
        radialDisplacementStartRadiansUp = math.acos(rvOuterHeight/totalHeight)
        elementsCountUpApex = elementsCountUp - elementsCountUpSeptum
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundSeptum
        radiansPerElementAroundLVFreeWall = (2.0*math.pi - rvArcAroundRadians)/elementsCountAroundLVFreeWall

        norl = elementsCountAroundLV
        nowl = elementsCountUp*elementsCountAroundLV + 1
        for n3 in range(2):
            a = totalHeight
            b = lvOuterRadius
            bVS = lvOuterRadius - lvFreeWallThickness + septumThickness
            if n3 == 0:
                a -= lvApexThickness
                b -= lvFreeWallThickness
                bVS = b

            lengthUp = 0.25*getApproximateEllipsePerimeter(a, b)
            lengthUpApex = getEllipseArcLength(a, b, 0.0, rvThetaUp)
            lengthUpSeptum = lengthUp - lengthUpApex
            elementSizeUpSeptum = lengthUpSeptum/(elementsCountUpSeptum - 1)
            if n3 == 1:
                elementSizeUpApex = (lengthUpApex - (2.0/3.0)*elementSizeUpSeptum)/(elementsCountUpApex + (2.0/3.0))
                elementSizeUpSeptumTransition = (2.0/3.0)*(elementSizeUpApex + elementSizeUpSeptum)
            else:
                elementSizeUpApex = (lengthUpApex - 0.5*elementSizeUpSeptum)/(elementsCountUpApex + 0.5)
                elementSizeUpSeptumTransition = 0.5*(elementSizeUpApex + elementSizeUpSeptum)
            print('length total apex septum', lengthUp, lengthUpApex, lengthUpSeptum, ' size apex septum', elementSizeUpApex, elementSizeUpSeptum)

            # apex node, noting s1, s2 is x, -y to get out outward s3
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -a ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ elementSizeUpApex, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ 0.0, -elementSizeUpApex, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, -lvApexThickness ])
            nodeIdentifier = nodeIdentifier + 1

            radiansUp = 0.0
            for n2 in range(elementsCountUp):

                if n2 <= elementsCountUpApex:
                    if n2 == elementsCountUpApex:
                        arcLengthUp = elementSizeUpSeptumTransition
                    else:
                        arcLengthUp = elementSizeUpApex
                    dArcLengthUp = elementSizeUpApex
                else:
                    arcLengthUp = elementSizeUpSeptum
                    dArcLengthUp = elementSizeUpSeptum
                radiansUp = updateEllipseAngleByArcLength(a, b, radiansUp, arcLengthUp)
                if n2 == (elementsCountUp - 1):
                    radiansUp = 0.5*math.pi
                #print(n3, n2, ' -> ', math.degrees(radiansUp), 'dArcLengthUp', dArcLengthUp)
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                z = -a*cosRadiansUp
                radius = b*sinRadiansUp
                bSeptum = b
                if (n3 == 1) and (n2 >= elementsCountUpApex):
                    bSeptum -= (lvFreeWallThickness - septumThickness)
                radiusSeptum = bSeptum*sinRadiansUp

                # get radial displacement of centre of septum, a function of radiansUp
                xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
                radialDisplacement = interpolateCubicHermite([0.0], [0.0], [septumBaseRadialDisplacement], [0.0], xiUp)[0]

                sx, sd1 = getSeptumPoints(rvArcAroundRadians, radiusSeptum, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundSeptum, z, n3)

                # do finite difference in xi2
                dxi = 1.0E-3
                dzr_dRadiansUp = [ a*sinRadiansUp, bSeptum*cosRadiansUp ]
                ds_dRadiansUp = vector.magnitude(dzr_dRadiansUp)
                dzr_ds = vector.normalise(dzr_dRadiansUp)
                ds_dxi = dArcLengthUp
                dz = dxi*ds_dxi*dzr_ds[0]
                dr = dxi*ds_dxi*dzr_ds[1]
                dxiUp_dxi = ds_dxi/(ds_dRadiansUp*(0.5*math.pi - radialDisplacementStartRadiansUp))
                dRadialDisplacement = dxi*dxiUp_dxi*interpolateCubicHermiteDerivative([0.0], [0.0], [septumBaseRadialDisplacement], [0.0], xiUp)[0]
                tx, td1 = getSeptumPoints(rvArcAroundRadians, radiusSeptum + dr, radialDisplacement + dRadialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundSeptum, z + dz, n3)

                for n1 in range(elementsCountAroundLV):

                    if (n1 > 0) and (n1 < elementsCountAroundSeptum):
                        n1mid = (elementsCountAroundSeptum + 1)//2
                        if n1 < n1mid:
                            nr = elementsCountAroundSeptum//2 - n1
                            x = [ sx[nr][0], -sx[nr][1], sx[nr][2] ]
                            dx_ds1 = [ -sd1[nr][0], sd1[nr][1], sd1[nr][2] ]
                            xt = [ tx[nr][0], -tx[nr][1], tx[nr][2] ]
                        else:
                            nr = n1 - n1mid
                            x = sx[nr]
                            dx_ds1 = sd1[nr]
                            xt = tx[nr]
                        dx_ds2 = [ (xt[c] - x[c])/dxi for c in range(3) ]
                    else:
                        if n1 == 0:
                            radiansAround = -0.5*rvArcAroundRadians
                        else:
                            radiansAround = 0.5*rvArcAroundRadians + (n1 - elementsCountAroundSeptum)*radiansPerElementAroundLVFreeWall
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [ radius*cosRadiansAround, radius*sinRadiansAround, z ]
                        dx_ds1 = [ radius*radiansPerElementAroundLVFreeWall*-sinRadiansAround, radius*radiansPerElementAroundLVFreeWall*cosRadiansAround, 0.0 ]
                        dx_ds2 = [ cosRadiansAround*ds_dxi*dzr_ds[1], sinRadiansAround*ds_dxi*dzr_ds[1], ds_dxi*dzr_ds[0] ]

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    if n3 == 1:
                        # get derivative3 from difference in coordinates
                        innerNode = nodes.findNodeByIdentifier(nodeIdentifier - nowl)
                        cache.setNode(innerNode)
                        result, xInner = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                        dx_ds3 = [ (x[c] - xInner[c]) for c in range(3) ]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier = nodeIdentifier + 1

        # RV nodes

        nidr = nodeIdentifier
        elementsCountUpRV = elementsCountUpSeptum + 1
        elementsCountAroundRV = elementsCountAroundSeptum + 2
        rvOuterWidthBase = rvWidth - septumBaseRadialDisplacement + rvFreeWallThickness

        rxOuter = []
        rd1Outer = []
        rd2Outer = []
        rd3Outer = []
        elementSizeUpRv = (lengthUp - (elementsCountUpApex + 0.5)*elementSizeUpApex)/(elementsCountUpRV - 0.5)
        elementSizeUpRvTransition = 0.5*(elementSizeUpApex + elementSizeUpRv)
        distUp = elementSizeUpApex*elementsCountUpApex + elementSizeUpRvTransition
        radiansUp = updateEllipseAngleByArcLength(a, b, 0.0, distUp)

        print('RV element size apex', elementSizeUpApex, ', transition', elementSizeUpRvTransition, ', rv', elementSizeUpRv, \
            ', total', elementsCountUpApex*elementsCountUpApex + elementSizeUpRvTransition + (elementsCountUpRV - 1)*elementSizeUpRv, ' vs ', lengthUp)

        for n2 in range(elementsCountUpRV):
            if n2 == (elementsCountUpRV - 1):
                radiansUp = 0.5*math.pi
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)

            radiansAround = -0.5*rvArcAroundRadians
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)

            lvOuterRadius = b*math.sin(radiansUp)
            ax = [ lvOuterRadius*cosRadiansAround, lvOuterRadius*sinRadiansAround, -a*cosRadiansUp ]
            z = ax[2]
            dEndMag = lvOuterRadius*radiansPerElementAroundLVFreeWall
            ad2 = [ b*cosRadiansUp*cosRadiansAround, b*cosRadiansUp*sinRadiansAround, a*sinRadiansUp ]
            scale = elementSizeUpRv/vector.magnitude(ad2)
            ad2 = [ d*scale for d in ad2 ]

            xiUp = 1.0 + z/rvOuterHeight

            rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUp, rvOuterWidthBase, rvExtraCrossRadiusBase)
            #print('z',z,'xiUp', xiUp, 'xiUpFast', xiUpFast, 'xiUpSlow', xiUpSlow)
            #print('rvAddWidthRadius =', rvAddWidthRadius, ', rvAddCrossRadius=', rvAddCrossRadius)
            rx, rd1 = getRvOuterPoints(rvArcAroundRadians, lvOuterRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundRV, dEndMag, z)

            if n2 < (elementsCountUpRV - 1):
                # get points at node position plus increment of xi2, to calculate derivative 2
                dxi = 1.0E-3
                xPlus = [ ax[i] + ad2[i]*dxi for i in range(3) ]
                lvOuterRadiusPlus = math.sqrt(xPlus[0]*xPlus[0] + xPlus[1]*xPlus[1])
                dEndMagPlus = dEndMag*lvOuterRadiusPlus/lvOuterRadius
                xiUpPlus = 1.0 + xPlus[2]/rvOuterHeight
                rvAddWidthRadiusPlus, rvAddCrossRadiusPlus = getRVOuterSize(xiUpPlus, rvOuterWidthBase, rvExtraCrossRadiusBase)
                px, pd1 = getRvOuterPoints(rvArcAroundRadians, lvOuterRadiusPlus, rvAddWidthRadiusPlus, rvAddCrossRadiusPlus, elementsCountAroundRV, dEndMagPlus, xPlus[2])

                rd2 = [ [(px[n][c] - rx[n][c])/dxi for c in range(3)] for n in range(len(rx)) ]
            else:
                rd2 = [ad2]*len(rx)

            rd3 = []
            for n in range(len(rx)):
                normal = vector.crossproduct3(rd1[n], rd2[n])
                mag = rvFreeWallThickness/vector.magnitude(normal)
                rd3.append([ c*mag for c in normal ])

            rxOuter.append(rx)
            rd1Outer.append(rd1)
            rd2Outer.append(rd2)
            rd3Outer.append(rd3)
            radiansUp = updateEllipseAngleByArcLength(a, b, radiansUp, elementSizeUpRv)

        # get inner RV nodes from outer
        rxInner = []
        rd1Inner = []
        rd2Inner = []
        rd3Inner = rd3Outer
        for n2 in range(elementsCountUpRV):
            ix = []
            id1 = []
            id2 = []
            rx = rxOuter[n2]
            rd1 = rd1Outer[n2]
            rd2 = rd2Outer[n2]
            rd3 = rd3Outer[n2]
            for n in range(elementsCountAroundRV - 1):
                ix.append([ (rx[n][c] - rd3[n][c]) for c in range(3) ] )
                unitRadial = vector.normalise(rd3[n])

                # calculate inner d1 from curvature around
                curvature = 0.0
                count = 0
                if n > 0:
                    curvature -= getCubicHermiteCurvature(rx[n - 1], rd1[n - 1], rx[n], rd1[n], unitRadial, 1.0)
                    count += 1
                if n < (elementsCountAroundRV - 2):
                    curvature -= getCubicHermiteCurvature(rx[n], rd1[n], rx[n + 1], rd1[n + 1], unitRadial, 0.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*rvFreeWallThickness
                id1.append([ factor*c for c in rd1[n] ])

                # calculate inner d2 from curvature up
                curvature = 0.0
                count = 0.0
                if n2 > 0:
                    curvature -= getCubicHermiteCurvature(rxOuter[n2 - 1][n], rd2Outer[n2 - 1][n], rx[n], rd2[n], unitRadial, 1.0)
                    count += 1
                if n2 < (elementsCountUpRV - 1):
                    curvature -= getCubicHermiteCurvature(rx[n], rd2[n], rxOuter[n2 + 1][n], rd2Outer[n2 + 1][n], unitRadial, 0.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*rvFreeWallThickness
                id2.append([ factor*c for c in rd2Outer[n2][n] ])

            rxInner.append(ix)
            rd1Inner.append(id1)
            rd2Inner.append(id2)

        zero = [ 0.0, 0.0, 0.0 ]
        rv_nids = [[],[]]
        for n3 in range(2):
            for n2 in range(elementsCountUpRV):
                nids = []
                rv_nids[n3].append(nids)
                if n3 == 0:
                    rx = rxInner[n2]
                    rd1 = rd1Inner[n2]
                    rd2 = rd2Inner[n2]
                    rd3 = rd3Inner[n2]
                else:
                    rx = rxOuter[n2]
                    rd1 = rd1Outer[n2]
                    rd2 = rd2Outer[n2]
                    rd3 = rd3Outer[n2]
                for n1 in range(elementsCountAroundRV - 1):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rx[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3[n1])
                    nid = nodeIdentifier
                    nodeIdentifier += 1
                    nids.append(nid)


        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)
        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        norr = elementsCountAroundRV - 1
        nowr = norr*elementsCountUpRV

        # LV elements

        # calculate values used for apex scale factors
        if elementsCountAroundSeptum > 2:
            radiansPerElementAroundSeptum = (rvArcAroundRadians - radiansPerElementAroundLVFreeWall)/(elementsCountAroundSeptum - 1)
            radiansPerElementAroundSeptumTransition = 0.5*(radiansPerElementAroundLVFreeWall + radiansPerElementAroundSeptum)
        else:
            radiansPerElementAroundSeptum = radiansPerElementAroundSeptumTransition = radiansPerElementAroundSeptum/elementsCountAroundSeptum
        radiansAround = -0.5*rvArcAroundRadians

        for e2 in range(elementsCountUp):

            for e1 in range(elementsCountAroundLV):

                eft1 = eft
                scalefactors = None
                if e2 == 0:
                    # create bottom apex elements, varying eft scale factor identifiers around apex
                    # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                    va = e1
                    vb = (e1 + 1)%elementsCountAroundLV
                    nids = [ 1       , 2 + va       , 2 + vb,
                             1 + nowl, 2 + va + nowl, 2 + vb + nowl ]
                    eft1 = tricubichermite.createEftShellApexBottom(va*100, vb*100)
                    # calculate general linear map coefficients
                    if e1 < elementsCountAroundSeptum:
                        deltaRadians = dRadians1 = dRadians2 = radiansPerElementAroundSeptum
                        if e1 == 0:
                            deltaRadians = radiansPerElementAroundSeptumTransition
                            dRadians1 = radiansPerElementAroundLVFreeWall
                        elif e1 == (elementsCountAroundSeptum - 1):
                            deltaRadians = radiansPerElementAroundSeptumTransition
                            dRadians2 = radiansPerElementAroundLVFreeWall
                    else:
                        deltaRadians = dRadians1 = dradians2 = radiansPerElementAroundLVFreeWall
                    radiansAroundNext = radiansAround + deltaRadians
                    radians1 = -radiansAround
                    radians2 = -radiansAroundNext
                    scalefactors = [ -1.0,
                        math.cos(radians1), math.sin(radians1), dRadians1,
                        math.cos(radians2), math.sin(radians2), dRadians2,
                        math.cos(radians1), math.sin(radians1), dRadians1,
                        math.cos(radians2), math.sin(radians2), dRadians2
                    ]
                    radiansAround = radiansAroundNext
                else:
                    bnil = 2 + norl*(e2 - 1) + e1
                    bnjl = 2 + norl*(e2 - 1) + (e1 + 1)%elementsCountAroundLV
                    nids = [ bnil       , bnjl       , bnil + norl       , bnjl + norl, \
                             bnil + nowl, bnjl + nowl, bnil + norl + nowl, bnjl + norl + nowl ]
                    if (e2 >= elementsCountUpApex) and (e1 < elementsCountAroundSeptum):
                        if (e2 == elementsCountUpApex) or (e1 == 0):
                            nids[4] = nidr + norr*(e2 - elementsCountUpApex) + e1
                        if e1 == 0:
                            nids[6] = nids[4] + norr
                        if (e2 == elementsCountUpApex) or (e1 == (elementsCountAroundSeptum - 1)):
                            nids[5] = nidr + norr*(e2 - elementsCountUpApex) + e1 + 1
                        if e1 == (elementsCountAroundSeptum - 1):
                            nids[7] = nids[5] + norr
                        if e2 == elementsCountUpApex:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            elif e1 == (elementsCountAroundSeptum - 1):
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            else:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                        else:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]
                            elif e1 == (elementsCountAroundSeptum - 1):
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                scalefactors = [ -1.0 ]

                result = elementtemplate1.defineField(coordinates, -1, eft1)

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() > 0:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                print('create element lv', elementIdentifier, result, result2, result3, nids)
                elementIdentifier = elementIdentifier + 1

        # RV elements

        nidl = 2 + (elementsCountUpApex - 1)*norl

        for e2 in range(elementsCountUpRV):

            for e1 in range(elementsCountAroundRV):

                eft1 = eft
                scalefactors = None

                if e2 == 0:
                    if (e1 == 0) or (e1 == (elementsCountAroundRV - 1)):
                        # skip as two fewer elements on bottom row
                        continue
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    if e1 == 1:
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    elif e1 == (elementsCountAroundRV - 2):
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    else:
                        remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    bnil = nidl + norl*e2 + e1 - 1
                    bnjl = bnil + 1
                    bnir = nidr + norr*e2 + e1 - 1
                    bnjr = bnir + 1
                    nids = [ bnil       , bnjl       , bnir       , bnjr, \
                             bnil + nowl, bnjl + nowl, bnir + nowr, bnjr + nowr ]
                    scalefactors = [ -1.0 ]
                else:
                    if e1 == 0:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if e2 == 1:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        bnil = nidl + norl*(e2 - 1) + e1
                        bnir = nidr + norr*(e2 - 1) + e1
                        nids = [ bnil       , bnir       , bnil + norl       , bnir + norr, \
                                 bnil + nowl, bnir + nowr, bnil + norl + nowl, bnir + norr + nowr ]
                        scalefactors = [ -1.0 ]
                    elif e1 == (elementsCountAroundRV - 1):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if e2 == 1:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        bnil = nidl + norl*(e2 - 1) + elementsCountAroundSeptum
                        bnir = nidr + norr*(e2 - 1) + elementsCountAroundRV - 2
                        nids = [ bnir       , bnil       , bnir + norr       , bnil + norl, \
                                 bnir + nowr, bnil + nowl, bnir + norr + nowr, bnil + norl + nowl ]
                        scalefactors = [ -1.0 ]
                    else:
                        bnir = nidr + norr*(e2 - 1) + e1 - 1
                        bnjr = bnir + 1
                        nids = [ bnir       , bnjr       , bnir + norr       , bnjr + norr, \
                                 bnir + nowr, bnjr + nowr, bnir + norr + nowr, bnjr + norr + nowr ]

                result = elementtemplate1.defineField(coordinates, -1, eft1)

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() > 0:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                print('create element rv', elementIdentifier, result, result2, result3, nids)
                elementIdentifier = elementIdentifier + 1


        fm.endChange()


    @staticmethod
    def refineMesh(meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        Stops at end of ventricles, hence can be called from ventriclesbase.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAround = options['Number of elements around LV free wall']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountAroundSeptum = options['Number of elements around septum']
        elementsCountBelowSeptum = options['Number of elements below septum']

        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']

        startRvElementIdentifier = elementsCountAround*elementsCountUp*elementsCountThroughLVWall + 1
        elementsCountUpSeptum = elementsCountUp - elementsCountBelowSeptum + 1
        limitRvElementIdentifier = startRvElementIdentifier + (elementsCountAroundSeptum + 4)*elementsCountUpSeptum + elementsCountAroundSeptum + 2

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            elementId = element.getIdentifier()
            if elementId < startRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < limitRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitRvElementIdentifier - 1):
                return  # finish on last so can continue in ventriclesbase
            element = meshrefinement._sourceElementiterator.next()

    @staticmethod
    def generateMesh(region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            MeshType_3d_heartventricles2.generateBaseMesh(region, options)
            return
        baseRegion = region.createRegion()
        MeshType_3d_heartventricles2.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region)
        MeshType_3d_heartventricles2.refineMesh(meshrefinement, options)


def getSeptumPoints(septumArcRadians, lvRadius, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundSeptum, z, n3):
        '''
        Symmetric 2-cubic interpolation of n septum elements around arc.
        :return: x[], dx_ds1[]
        '''
        radiansPerElementAroundLVFreeWall = (2.0*math.pi - septumArcRadians)/elementsCountAroundLVFreeWall
        # get cubic curve with arc length scaling across to centre of septum
        radiansAround = 0.5*septumArcRadians
        circleArcLength = lvRadius*radiansAround
        v1 = [ lvRadius - radialDisplacement, 0.0, z ]
        d1 = [ 0.0, circleArcLength, 0.0 ]
        v2 = [ lvRadius*math.cos(radiansAround), lvRadius*math.sin(radiansAround), z ]
        d2 = [ -circleArcLength*math.sin(radiansAround), circleArcLength*math.cos(radiansAround), 0.0 ]
        cubicArcLength = computeCubicHermiteArcLength(v1, d1, v2, d2, False)
        scale = cubicArcLength/circleArcLength
        d1 = [ d*scale for d in d1 ]
        d2 = [ d*scale for d in d2 ]
        elementLengthAroundSeptum = (2.0*cubicArcLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundSeptum - 1)
        elementLengthAroundTransition = 0.5*(lvRadius*radiansPerElementAroundLVFreeWall + elementLengthAroundSeptum)
        lengthAroundStart = (elementsCountAroundSeptum % 2)*elementLengthAroundSeptum*0.5
        x = []
        dx_ds1 = []
        for n1 in range((elementsCountAroundSeptum + 1)//2):
                xi = (n1*elementLengthAroundSeptum + lengthAroundStart)/cubicArcLength
                pos1 = interpolateCubicHermite(v1, d1, v2, d2, xi)
                x.append([ v for v in pos1 ])
                deriv1 = interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
                scale = elementLengthAroundSeptum/vector.magnitude(deriv1)
                dx_ds1.append([ d*scale for d in deriv1 ])
        return x, dx_ds1


def getRvOuterPoints(rvArcAroundRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAround, dEndMag, z):
    '''
    Get array of points and derivatives around RV, centred at +x axis.
    LV is assumed to be circular, centred at x, y = (0,0)
    :param rvArcAroundRadians: Angle in radians around outside of RV to tangent nodes where it meets LV.
    :param lvRadius: Radius of the LV.
    :param rvAddWidthRadius: Radius to add to LV to get RV width at centre.
    :param rvAddCrossRadius: Additional radius to add only laterally around septum.
    :param elementsCountAround: Number of elements around RV.
    :param dEndMag: Magnitude of the external derivative at each end.
    :param z: Z coordinate to give to all values of x[]. dx_ds1[2] is all zero.
    :return: Arrays x[], dx_ds1[], starting at first point after end at lowest angle around.
    '''
    #print('\ngetRvOuterPoints', rvArcAroundRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAround, dEndMag, z)
    startRadians = -0.5*rvArcAroundRadians
    ellipseStartRadians = startRadians + 0.5*rvArcAroundRadians/elementsCountAround
    b = lvRadius + rvAddCrossRadius  # cross radius
    theta = math.asin(-lvRadius*math.sin(ellipseStartRadians)/b) - math.pi
    #print('phi',ellipseStartRadians)
    #print('theta',theta)
    a = (lvRadius*(math.cos(ellipseStartRadians) - 1.0) - rvAddWidthRadius)/(math.cos(theta) - 1.0)  # width radius
    cx = lvRadius + rvAddWidthRadius - a
    #print(ellipseStartRadians,'a',a,'b',b,'cx',cx)
    # get cubic curve joining start point with half ellipse, first with unit derivatives then compute arc length
    x1 = ( lvRadius*math.cos(startRadians), lvRadius*math.sin(startRadians) )
    d1 = ( -math.sin(startRadians), math.cos(startRadians) )
    x2 = ( cx, -b )
    d2 = ( 1.0, 0.0 )
    cubicArcLength = computeCubicHermiteArcLength(x1, d1, x2, d2, True)
    d1 = ( d1[0]*cubicArcLength, d1[1]*cubicArcLength )
    d2 = ( d2[0]*cubicArcLength, d2[1]*cubicArcLength )
    halfEllipsePerimeter = 0.5*getApproximateEllipsePerimeter(a, b)
    totalLength = 2.0*cubicArcLength + halfEllipsePerimeter
    #print('length cubic', cubicArcLength, 'halfEllipsePerimeter',halfEllipsePerimeter,'total',totalLength)
    midLength = (totalLength - dEndMag)/(elementsCountAround - 1)
    endLength = 0.5*(dEndMag + midLength)
    #print('midLength',midLength,'endLength',endLength)
    length = 0.0
    angle = 0.0
    if (elementsCountAround % 2) == 1:
        length = -0.5*midLength
        angle = updateEllipseAngleByArcLength(a, b, angle, length)
    x = []
    dx_ds1 = []
    quarterEllipsePerimeter = halfEllipsePerimeter*0.5
    for n in range(elementsCountAround//2):
        if length >= -quarterEllipsePerimeter:
            v = ( cx + a*math.cos(angle), b*math.sin(angle) )
            t = ( -a*math.sin(angle), b*math.cos(angle) )
            angle = updateEllipseAngleByArcLength(a, b, angle, -midLength)
        else:
            xi = 1.0 + (length + quarterEllipsePerimeter)/cubicArcLength
            #print('interpolate',x1,d1,x2,d2, xi)
            v = interpolateCubicHermite(x1, d1, x2, d2, xi)
            t = interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        x.insert(0, [ v[0], v[1], z ] )
        scale = midLength/math.sqrt(t[0]*t[0] + t[1]*t[1])
        dx_ds1.insert(0, [ t[0]*scale, t[1]*scale, 0.0 ])
        length -= midLength
    length += midLength - endLength
    #print('length',length,' remainder',totalLength*0.5+length)
    # mirror in y (about x axis) to get other side
    for n in range((elementsCountAround - 1)//2 - 1, -1, -1):
        x.append([ x[n][0], -x[n][1], z ])
        dx_ds1.append([ -dx_ds1[n][0], dx_ds1[n][1], 0.0 ])
    return x, dx_ds1

def getRVOuterSize(xiUp, rvWidthRadius, rvExtraCrossRadiusBase):
    xiUpFast = 1.0 - (1.0 - xiUp)*(1.0 - xiUp)
    xiUpSlow = xiUp*xiUp
    rvAddWidthRadius = interpolateCubicHermite([0.0], [0.0], [rvWidthRadius], [0.0], xiUpFast)[0]
    rvAddCrossRadius = interpolateCubicHermite([0.0], [0.0], [rvExtraCrossRadiusBase], [0.0], xiUpSlow)[0]
    return rvAddWidthRadius, rvAddCrossRadius
