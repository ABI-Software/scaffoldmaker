"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valves regions.
"""

from __future__ import division
import math
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles2 import MeshType_3d_heartventricles2
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventriclesbase2(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Ventricles with Base 2'

    @staticmethod
    def getDefaultOptions():
        options = MeshType_3d_heartventricles2.getDefaultOptions()
        # only works with particular numbers of elements around
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around septum'] = 7
        # works best with particular numbers of elements up
        options['Number of elements up apex'] = 1
        options['Number of elements up septum'] = 3
        # additional options
        options['Number of elements around atria'] = 8
        options['Atrial septum thickness'] = 0.075
        options['Atria major axis rotation degrees'] = 30.0
        options['Base height'] = 0.15
        options['Base thickness'] = 0.05
        options['LV outlet outer diameter'] = 0.3
        options['RV outlet outer diameter'] = 0.3
        options['Outlet element length'] = 0.1
        options['Outlet incline degrees'] = 10.0
        options['Outlet spacing'] = 0.03
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventricles2.getOrderedOptionNames()
        optionNames += [
            'Number of elements around atria',
            'Atrial septum thickness',
            'Atria major axis rotation degrees',
            'Base height',
            'Base thickness',
            'LV outlet outer diameter',
            'RV outlet outer diameter',
            'Outlet element length',
            'Outlet incline degrees',
            'Outlet spacing']
        # want refinement options last
        for optionName in [
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_heartventricles2.checkOptions(options)
        # only works with particular numbers of elements around
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around septum'] = 7
        # need even number of refine surface elements for elements with hanging nodes to conform
        #if (options['Refine number of elements surface'] % 2) == 1:
        #    options['Refine number of elements surface'] += 1
        if options['Number of elements around atria'] < 6:
            options['Number of elements around atria'] = 6
        for key in [
            'Atrial septum thickness',
            'Base height',
            'Base thickness',
            'LV outlet outer diameter',
            'RV outlet outer diameter',
            'Outlet element length',
            'Outlet spacing']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Atria major axis rotation degrees'] < -75.0:
            options['Atria major axis rotation degrees'] = -75.0
        elif options['Atria major axis rotation degrees'] > 75.0:
            options['Atria major axis rotation degrees'] = 75.0

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
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundSeptum
        elementsCountUpApex = options['Number of elements up apex']
        elementsCountUpSeptum = options['Number of elements up septum']
        elementsCountUpLV = elementsCountUpApex + elementsCountUpSeptum
        elementsCountUpRV = elementsCountUpSeptum + 1
        elementsCountAroundRV = elementsCountAroundSeptum + 2
        elementsCountAtrialSeptum = elementsCountAroundSeptum - 5
        elementsCountAroundAtria = options['Number of elements around atria']
        totalHeight = options['Total height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        lvApexThickness = options['LV apex thickness']
        rvHeight = options['RV height']
        rvArcAroundRadians = math.radians(options['RV arc around degrees'])
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        rvExtraCrossRadiusBase = options['RV extra cross radius base']
        vSeptumThickness = options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = options['Ventricular septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aSeptumThickness = options['Atrial septum thickness']
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        lvOutletRadius = options['LV outlet outer diameter']*0.5
        rvOutletRadius = options['RV outlet outer diameter']*0.5
        outletElementLength = options['Outlet element length']
        outletInclineRadians = math.radians(options['Outlet incline degrees'])
        outletSpacing = options['Outlet spacing']

        # generate heartventricles2 model to add base plane to
        MeshType_3d_heartventricles2.generateBaseMesh(region, options)

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

        nodeIdentifier = startNodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

        # node offsets for row, wall in LV, plus first LV node on inside top
        norl = elementsCountAroundLV
        nowl = 1 + elementsCountUpLV*norl
        nidl = nowl - norl + 1
        # node offsets for row, wall in RV, plus first RV node on inside top
        norr = (elementsCountAroundRV - 1)
        nowr = elementsCountUpRV*norr
        nidr = nowl*2 + 1 + nowr - norr
        #print('nidl',nidl,'nidr',nidr)

        # LV outlet
        elementsCountAroundOutlet = 6
        # GRC Set properly:
        defaultOutletScale3 = 0.5
        node = nodes.findNodeByIdentifier(nidl + elementsCountAtrialSeptum + 1)
        cache.setNode(node)
        result, cx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        node = nodes.findNodeByIdentifier(nidl)
        cache.setNode(node)
        result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        node = nodes.findNodeByIdentifier(nidr)
        cache.setNode(node)
        result, bx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        bx = [ 0.5*(ax[i] + bx[i]) for i in range(2) ]
        bx.append(ax[2])
        ax = [ (bx[c] - cx[c]) for c in range(3) ]
        axMag = vector.magnitude(ax)
        ax = vector.normalise(ax)
        baseRotationRadians = math.atan2(ax[1], ax[0])
        #print('baseRotationRadians', baseRotationRadians)

        cosOutletInclineRadians = math.cos(outletInclineRadians)
        sinOutletInclineRadians = math.sin(outletInclineRadians)
        lvOutletCentre = [
            cx[0] - ax[0]*lvOutletRadius,
            cx[1] - ax[1]*lvOutletRadius,
            baseHeight + baseThickness + sinOutletInclineRadians*lvOutletRadius ]
        loAxis1 = [ lvOutletRadius*ax[c] for c in range(3) ]
        loAxis2 = [ -loAxis1[1]*cosOutletInclineRadians, loAxis1[0]*cosOutletInclineRadians, -lvOutletRadius*sinOutletInclineRadians ]
        loAxis3 = vector.crossproduct3(loAxis1, loAxis2)
        scale = outletElementLength/vector.magnitude(loAxis3)
        loAxis3 = [ v*scale for v in loAxis3 ]

        radiansPerElementAroundOutlet = 2.0*math.pi/elementsCountAroundOutlet
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = loAxis3
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        lvOutletNodeId = []
        for n1 in range(elementsCountAroundOutlet):
            radiansAround = n1*radiansPerElementAroundOutlet
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            outletScale3 = outletSpacing/lvOutletRadius if (n1 == 3) else defaultOutletScale3
            for c in range(3):
                x[c] = lvOutletCentre[c] + loAxis1[c]*cosRadiansAround + loAxis2[c]*sinRadiansAround
                dx_ds1[c] = radiansPerElementAroundOutlet*(loAxis1[c]*-sinRadiansAround + loAxis2[c]*cosRadiansAround)
                dx_ds3[c] = outletScale3*(loAxis1[c]*cosRadiansAround + loAxis2[c]*sinRadiansAround)
            #if n1 == 4:
            #    d3inclineRadians = math.pi/6.0
            #    dx_ds3 = [ dx_ds1[0]*math.cos(d3inclineRadians),
            #               dx_ds1[1]*math.cos(d3inclineRadians),
            #               dx_ds1[2]*math.sin(d3inclineRadians) ]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            lvOutletNodeId.append(nodeIdentifier)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
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
        outletCentreSpacing = lvOutletRadius + outletSpacing + rvOutletRadius
        rvOutletCentre = [ (lvOutletCentre[c] - outletCentreSpacing*ax[c]) for c in range(3) ]
        roAxis1 = [ rvOutletRadius*ax[c] for c in range(3) ]
        roAxis2 = [ -roAxis1[1]*cosOutletInclineRadians, roAxis1[0]*cosOutletInclineRadians, rvOutletRadius*sinOutletInclineRadians ]
        roAxis3 = vector.crossproduct3(roAxis1, roAxis2)
        scale = outletElementLength/vector.magnitude(roAxis3)
        roAxis3 = [ v*scale for v in roAxis3 ]

        dx_ds2 = roAxis3
        rvOutletNodeId = []
        for n1 in range(elementsCountAroundOutlet):
            radiansAround = n1*radiansPerElementAroundOutlet
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            outletScale3 = outletSpacing/rvOutletRadius if (n1 == 0) else defaultOutletScale3
            for c in range(3):
                x[c] = rvOutletCentre[c] + roAxis1[c]*cosRadiansAround + roAxis2[c]*sinRadiansAround
                dx_ds1[c] = radiansPerElementAroundOutlet*(roAxis1[c]*-sinRadiansAround + roAxis2[c]*cosRadiansAround)
                dx_ds3[c] = outletScale3*(roAxis1[c]*cosRadiansAround + roAxis2[c]*sinRadiansAround)
            #if (n1 == 1) or (n1 == (elementsCountAroundOutlet - 1)):
            #    d3inclineRadians = math.pi/6.0
            #    dx_ds3 = [ dx_ds1[0]*math.cos(d3inclineRadians),
            #               dx_ds1[1]*math.cos(d3inclineRadians),
            #               dx_ds1[2]*math.sin(d3inclineRadians) ]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            rvOutletNodeId.append(nodeIdentifier)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier += 1

        if False:
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
                node = nodes.createNode(nodeIdentifier, nodetemplate)
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

        #print('\nelementsCountAroundAtria = ',elementsCountAroundAtria)
        # GRC parameter?
        atriumInletSlopeRadians = math.pi/6.0
        # GRC change rvFreeWallThickness to new parameter?
        atriumInletSlopeLength = rvFreeWallThickness
        atriumInletSlopeHeight = rvFreeWallThickness*math.tan(atriumInletSlopeRadians)

        # GRC revisit:
        aInnerMajorMag = 0.8*(lvOuterRadius - lvFreeWallThickness - 0.5*vSeptumBaseRadialDisplacement)
        aInnerMinorMag = 0.9*(lvOuterRadius - lvFreeWallThickness - lvOutletRadius)
        #print('inner mag major', aInnerMajorMag, 'minor', aInnerMinorMag)
        aOuterMajorMag = aInnerMajorMag + atriumInletSlopeLength
        aOuterMinorMag = aInnerMinorMag + atriumInletSlopeLength
        #print('outer mag major', aOuterMajorMag, 'minor', aOuterMinorMag)

        laMajorAxisRadians = baseRotationRadians + 0.5*math.pi - aMajorAxisRadians
        laInnerMajor = [  aInnerMajorMag*math.cos(laMajorAxisRadians), aInnerMajorMag*math.sin(laMajorAxisRadians), 0.0 ]
        laInnerMinor = [ -aInnerMinorMag*math.sin(laMajorAxisRadians), aInnerMinorMag*math.cos(laMajorAxisRadians), 0.0 ]
        laOuterMajor = [  aOuterMajorMag*math.cos(laMajorAxisRadians), aOuterMajorMag*math.sin(laMajorAxisRadians), 0.0 ]
        laOuterMinor = [ -aOuterMinorMag*math.sin(laMajorAxisRadians), aOuterMinorMag*math.cos(laMajorAxisRadians), 0.0 ]

        raMajorAxisRadians = baseRotationRadians - 0.5*math.pi + aMajorAxisRadians
        raInnerMajor = [  aInnerMajorMag*math.cos(raMajorAxisRadians), aInnerMajorMag*math.sin(raMajorAxisRadians), 0.0 ]
        raInnerMinor = [ -aInnerMinorMag*math.sin(raMajorAxisRadians), aInnerMinorMag*math.cos(raMajorAxisRadians), 0.0 ]
        raOuterMajor = [  aOuterMajorMag*math.cos(raMajorAxisRadians), aOuterMajorMag*math.sin(raMajorAxisRadians), 0.0 ]
        raOuterMinor = [ -aOuterMinorMag*math.sin(raMajorAxisRadians), aOuterMinorMag*math.cos(raMajorAxisRadians), 0.0 ]

        # get la angle intersecting with cruxLeft, thence laCentre

        rotRadians = baseRotationRadians + 0.5*math.pi
        cosRotRadians = math.cos(rotRadians)
        sinRotRadians = math.sin(rotRadians)
        cruxLeftModX = (cruxLeft[0] - cruxCentre[0])*cosRotRadians + (cruxLeft[1] - cruxCentre[1])*sinRotRadians
        cruxLeftModY = (cruxLeft[0] - cruxCentre[0])*-sinRotRadians + (cruxLeft[1] - cruxCentre[1])*cosRotRadians

        #cruxLeftModX = (cruxLeft[0] - cruxCentre[0])  # GRC temp
        #cruxLeftModY = (cruxLeft[1] - cruxCentre[1])  # GRC temp

        axInnerMod = aInnerMajorMag*math.cos(aMajorAxisRadians)
        bxInnerMod = aInnerMinorMag*math.sin(aMajorAxisRadians)
        #print('Inner axMod', axInnerMod, 'bxMod', bxInnerMod)
        laSeptumRadians = math.atan2(bxInnerMod, axInnerMod)
        #print('laSeptumRadians', laSeptumRadians)
        raSeptumRadians = -laSeptumRadians
        laCentreModX = -0.5*aSeptumThickness - axInnerMod*math.cos(laSeptumRadians) - bxInnerMod*math.sin(laSeptumRadians)

        axMod = aOuterMajorMag*math.cos(aMajorAxisRadians)
        ayMod = aOuterMajorMag*-math.sin(aMajorAxisRadians)
        bxMod = aOuterMinorMag*math.sin(aMajorAxisRadians)
        byMod = aOuterMinorMag*math.cos(aMajorAxisRadians)
        #print('Outer axMod', axMod, 'bxMod', bxMod)

        dX = cruxLeftModX - laCentreModX
        #print('laCentreModX', laCentreModX, 'cruxLeftModX', cruxLeftModX, 'dX', dX)
        # iterate with Newton's method to get laCruxLeftRadians
        laCruxLeftRadians = math.pi*0.5
        iters = 0
        fTol = aOuterMajorMag*1.0E-10
        while True:
            cosAngle = math.cos(laCruxLeftRadians)
            sinAngle = math.sin(laCruxLeftRadians)
            f = axMod*cosAngle + bxMod*sinAngle - dX
            if math.fabs(f) < fTol:
                break;
            df = -axMod*sinAngle + bxMod*cosAngle
            #print(iters, '. theta', laCruxLeftRadians, 'f', f,'df',df,'-->',laCruxLeftRadians - f/df)
            laCruxLeftRadians -= f/df
            iters += 1
            if iters == 100:
                print('No convergence!')
                break
        #print(iters,'iters : laCruxLeftRadians', laCruxLeftRadians)
        laCentreModY = cruxLeftModY - ayMod*math.cos(laCruxLeftRadians) - byMod*math.sin(laCruxLeftRadians)

        #print('laCentreMod', laCentreModX, laCentreModY)
        laCentreX = cruxCentre[0] + laCentreModX*cosRotRadians + laCentreModY*-sinRotRadians
        laCentreY = cruxCentre[1] + laCentreModX*sinRotRadians + laCentreModY*cosRotRadians

        raCruxLeftRadians = -laCruxLeftRadians
        raCentreX = cruxCentre[0] - laCentreModX*cosRotRadians + laCentreModY*-sinRotRadians
        raCentreY = cruxCentre[1] - laCentreModX*sinRotRadians + laCentreModY*cosRotRadians

        atrialPerimeterLength = getApproximateEllipsePerimeter(aOuterMajorMag, aOuterMinorMag)
        atrialSeptumCentreToCruxLeftLength = getEllipseArcLength(aOuterMajorMag, aOuterMinorMag, laSeptumRadians, laCruxLeftRadians)
        atrialSeptumElementLength = atrialSeptumCentreToCruxLeftLength/(1.0 + elementsCountAtrialSeptum*0.5)
        atrialFreeWallElementLength = (atrialPerimeterLength - atrialSeptumElementLength*(elementsCountAtrialSeptum + 3)) \
            / (elementsCountAroundAtria - elementsCountAtrialSeptum - 3)
        atrialTransitionElementLength = 0.5*(atrialSeptumElementLength + atrialFreeWallElementLength)

        #print('lengths septum', atrialSeptumElementLength, 'transition', atrialTransitionElementLength, 'freewall', atrialFreeWallElementLength)
        #print('total length', atrialSeptumElementLength*(elementsCountAtrialSeptum + 2) + 2*atrialTransitionElementLength \
        #    + (elementsCountAroundAtria - elementsCountAtrialSeptum - 4)*atrialFreeWallElementLength, 'vs.', atrialPerimeterLength)

        laRadians = []
        laOuterDerivatives = []
        radiansAround = laSeptumRadians
        if (elementsCountAtrialSeptum % 2) == 1:
            radiansAround = updateEllipseAngleByArcLength(aOuterMajorMag, aOuterMinorMag, radiansAround, 0.5*atrialSeptumElementLength)
        outerDerivative = atrialSeptumElementLength
        lan1CruxLimit = elementsCountAtrialSeptum//2 + 1
        lan1SeptumLimit = elementsCountAroundAtria - lan1CruxLimit - 1
        #print('lan1CruxLimit', lan1CruxLimit, 'lan1SeptumLimit', lan1SeptumLimit)
        for n1 in range(elementsCountAroundAtria):
            laRadians.append(radiansAround)
            laOuterDerivatives.append(outerDerivative)
            if (n1 < lan1CruxLimit) or (n1 > lan1SeptumLimit):
                elementLength = atrialSeptumElementLength
                outerDerivative = atrialSeptumElementLength
            elif n1 == lan1CruxLimit:
                elementLength = atrialTransitionElementLength
                outerDerivative = atrialFreeWallElementLength
            elif n1 == lan1SeptumLimit:
                elementLength = atrialTransitionElementLength
                outerDerivative = atrialSeptumElementLength
            else:
                elementLength = atrialFreeWallElementLength
                outerDerivative = atrialFreeWallElementLength
            #print(n1,': elementLength', elementLength, 'outerDerivative', outerDerivative)
            radiansAround = updateEllipseAngleByArcLength(aOuterMajorMag, aOuterMinorMag, radiansAround, elementLength)
        laInnerDerivatives = []
        finalArcLength = prevArcLength = getEllipseArcLength(aInnerMajorMag, aInnerMinorMag, laRadians[-1] - 2.0*math.pi, laRadians[0])
        for n1 in range(elementsCountAroundAtria):
            if n1 == (elementsCountAroundAtria - 1):
                nextArcLength = finalArcLength
            else:
                nextArcLength = getEllipseArcLength(aInnerMajorMag, aInnerMinorMag, laRadians[n1], laRadians[n1 + 1])
            if laOuterDerivatives[n1] is atrialSeptumElementLength:
                arcLength = min(prevArcLength, nextArcLength)
            else:
                arcLength = max(prevArcLength, nextArcLength)
            laInnerDerivatives.append(arcLength)
            prevArcLength = nextArcLength

        #print('raRadians', laRadians)
        #print('laOuterDerivatives', laOuterDerivatives)
        #print('laInnerDerivatives', laInnerDerivatives)

        raRadians = []
        raInnerDerivatives = []
        raOuterDerivatives = []
        for n1 in range(elementsCountAroundAtria):
            raRadians.append(2.0*math.pi - laRadians[-n1])
            raInnerDerivatives.append(laInnerDerivatives[-n1])
            raOuterDerivatives.append(laOuterDerivatives[-n1])
        # fix first one so not out by 2pi
        raRadians[0] = raSeptumRadians

        laNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        laNodeId[1][1] = lvOutletNodeId[0]
        laNodeId[1][2] = lvOutletNodeId[-1]

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = laRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    laCentreX + cosRadiansAround*laInnerMajor[0] + sinRadiansAround*laInnerMinor[0],
                    laCentreY + cosRadiansAround*laInnerMajor[1] + sinRadiansAround*laInnerMinor[1],
                    cruxCentre[2] - atriumInletSlopeHeight ]
                outer = [
                    laCentreX + cosRadiansAround*laOuterMajor[0] + sinRadiansAround*laOuterMinor[0],
                    laCentreY + cosRadiansAround*laOuterMajor[1] + sinRadiansAround*laOuterMinor[1],
                    cruxCentre[2] ]

                if (n3 == 1) and ((n1 <= lan1CruxLimit) or (n1 > (lan1SeptumLimit + 1))):
                    continue  # already have a node from crux or will get from right atrial septum
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                laNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        -sinRadiansAround*laInnerMajor[0] + cosRadiansAround*laInnerMinor[0],
                        -sinRadiansAround*laInnerMajor[1] + cosRadiansAround*laInnerMinor[1],
                        0.0 ]
                    scale1 = laInnerDerivatives[n1]
                else:
                    dx_ds1 = [
                        -sinRadiansAround*laOuterMajor[0] + cosRadiansAround*laOuterMinor[0],
                        -sinRadiansAround*laOuterMajor[1] + cosRadiansAround*laOuterMinor[1],
                        0.0 ]
                    scale1 = laOuterDerivatives[n1]
                scale1 /= vector.magnitude(dx_ds1)
                dx_ds1 = [ d*scale1 for d in dx_ds1 ]
                dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                dx_ds2 = [
                    dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                    dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                    dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                # GRC check scaling here:
                scale2 = inner[2]/vector.magnitude(dx_ds2)
                dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if False:
            # show axes of left atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, cruxCentre[2] - atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laInnerMajor[0], laInnerMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laInnerMinor[0], laInnerMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, cruxCentre[2] ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laOuterMajor[0], laOuterMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laOuterMinor[0], laOuterMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        ran1SeptumLimit = elementsCountAtrialSeptum//2
        ran1CruxLimit = elementsCountAroundAtria - ran1SeptumLimit - 1
        raNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        raNodeId[1][0] = laNodeId[0][0]
        raNodeId[1][1] = laNodeId[0][-1]
        raNodeId[1][-2] = lvOutletNodeId[1]
        raNodeId[1][-1] = lvOutletNodeId[0]

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = raRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    raCentreX + cosRadiansAround*raInnerMajor[0] + sinRadiansAround*raInnerMinor[0],
                    raCentreY + cosRadiansAround*raInnerMajor[1] + sinRadiansAround*raInnerMinor[1],
                    cruxCentre[2] - atriumInletSlopeHeight ]
                outer = [
                    raCentreX + cosRadiansAround*raOuterMajor[0] + sinRadiansAround*raOuterMinor[0],
                    raCentreY + cosRadiansAround*raOuterMajor[1] + sinRadiansAround*raOuterMinor[1],
                    cruxCentre[2] ]

                if (n3 == 1) and ((n1 <= ran1SeptumLimit) or (n1 >= ran1CruxLimit)):
                    continue  # already have a node from crux or will get from left atrial septum
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                raNodeId[n3][n1] = nodeIdentifier
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, inner if (n3 == 0) else outer)
                if n3 == 0:
                    dx_ds1 = [
                        -sinRadiansAround*raInnerMajor[0] + cosRadiansAround*raInnerMinor[0],
                        -sinRadiansAround*raInnerMajor[1] + cosRadiansAround*raInnerMinor[1],
                        0.0 ]
                    scale1 = raInnerDerivatives[n1]
                else:
                    dx_ds1 = [
                        -sinRadiansAround*raOuterMajor[0] + cosRadiansAround*raOuterMinor[0],
                        -sinRadiansAround*raOuterMajor[1] + cosRadiansAround*raOuterMinor[1],
                        0.0 ]
                    scale1 = raOuterDerivatives[n1]
                scale1 /= vector.magnitude(dx_ds1)
                dx_ds1 = [ d*scale1 for d in dx_ds1 ]
                dx_ds3 = [ outer[0] - inner[0], outer[1] - inner[1], outer[2] - inner[2] ]
                dx_ds2 = [
                    dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                    dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                    dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                # GRC check scaling here:
                scale2 = inner[2]/vector.magnitude(dx_ds2)
                dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if False:
            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, cruxCentre[2] - atriumInletSlopeHeight ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raInnerMajor[0], raInnerMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raInnerMinor[0], raInnerMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, cruxCentre[2] ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raOuterMajor[0], raOuterMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raOuterMinor[0], raOuterMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, cruxCentre[2] ])
            nodeIdentifier += 1

        laNodeId[1][-1] = raNodeId[0][1]
        laNodeId[1][0] = raNodeId[0][0]

        # fix up shared atrial septum nodes and derivatives
        dx_ds2 = [ 0.0, 0.0, cruxCentre[2] - atriumInletSlopeHeight ]
        for n1 in range(-1, 2):
            node1 = nodes.findNodeByIdentifier(laNodeId[0][n1])
            cache.setNode(node1)
            result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
            node2 = nodes.findNodeByIdentifier(laNodeId[1][n1])
            cache.setNode(node2)
            result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
            dx_ds3 = [ x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] ]
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ -d for d in dx_ds3 ])
            cache.setNode(node1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        fm.endChange()

    @staticmethod
    def refineMesh(meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)

        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountAcrossSeptum = options['Number of elements across septum']
        elementsCountBelowSeptum = options['Number of elements below septum']

        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']

        MeshType_3d_heartventricles2.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startBaseLvElementIdentifier = element.getIdentifier()
        startBaseRvElementIdentifier = startBaseLvElementIdentifier + 16
        startLvOutletElementIdentifier = startBaseRvElementIdentifier + 11
        limitLvOutletElementIdentifier = startLvOutletElementIdentifier + 6

        startHangingElementIdentifier = startBaseRvElementIdentifier + 4
        limitHangingElementIdentifier = startHangingElementIdentifier + 4

        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            elementId = element.getIdentifier()
            if elementId < startBaseRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < startLvOutletElementIdentifier:
                numberInXi3 = refineElementsCountThroughRVWall
                if (elementId >= startHangingElementIdentifier) and (elementId < limitHangingElementIdentifier):
                    numberInXi1 //= 2
            elif elementId < limitLvOutletElementIdentifier:
                numberInXi3 = 1
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitLvOutletElementIdentifier - 1):
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
            MeshType_3d_heartventriclesbase2.generateBaseMesh(region, options)
            return
        baseRegion = region.createRegion()
        MeshType_3d_heartventriclesbase2.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region)
        MeshType_3d_heartventriclesbase2.refineMesh(meshrefinement, options)
