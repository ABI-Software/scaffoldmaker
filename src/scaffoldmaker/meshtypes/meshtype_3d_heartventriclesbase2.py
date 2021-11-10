"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valve regions.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles2 import MeshType_3d_heartventricles2
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getEllipseArcLength, getEllipseRadiansToX, updateEllipseAngleByArcLength
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import interpolateNodesCubicHermite


class MeshType_3d_heartventriclesbase2(Scaffold_base):
    '''
    Generates a 3-D heart ventricles with base plane model, ready to attach the
    atria, mitral and tricuspid valves, with LV + RV outlets ready to attach
    aorta and pulmonary trunk and their valve regions.
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles with Base 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = MeshType_3d_heartventricles2.getDefaultOptions(parameterSetName)
        # only works with particular numbers of elements around
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around ventricular septum'] = 7
        options['Number of elements around atria'] = 8
        options['Number of elements around atrial septum'] = 2
        # works best with particular numbers of elements up
        options['Number of elements up LV apex'] = 1
        options['Number of elements up ventricular septum'] = 4
        # additional options
        options['Atria base inner major axis length'] = 0.55
        options['Atria base inner minor axis length'] = 0.42
        options['Atria major axis rotation degrees'] = 40.0
        options['Atrial septum thickness'] = 0.06
        options['Atrial base wall thickness'] = 0.05
        options['Atrial base slope degrees'] = 30.0
        options['Base height'] = 0.12
        options['Base thickness'] = 0.06
        options['Fibrous ring thickness'] = 0.01
        options['LV outlet inner diameter'] = 0.3
        options['LV outlet wall thickness'] = 0.025
        options['RV outlet inner diameter'] = 0.27
        options['RV outlet wall thickness'] = 0.025
        options['Ventricles outlet element length'] = 0.1
        options['Ventricles outlet incline degrees'] = 15.0
        options['Ventricles outlet spacing'] = 0.04
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventricles2.getOrderedOptionNames()
        optionNames.insert(4, 'Number of elements around atria')
        optionNames.insert(5, 'Number of elements around atrial septum')
        optionNames += [
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atrial septum thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Base height',
            'Base thickness',
            'Fibrous ring thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Ventricles outlet element length',
            'Ventricles outlet incline degrees',
            'Ventricles outlet spacing']
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
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = MeshType_3d_heartventricles2.checkOptions(options)
        # only works with particular numbers of elements around
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around ventricular septum'] = 7
        options['Number of elements around atria'] = 8
        options['Number of elements around atrial septum'] = 2
        for key in [
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atrial septum thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Base height',
            'Base thickness',
            'Fibrous ring thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Ventricles outlet element length',
            'Ventricles outlet spacing']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Atria major axis rotation degrees'] < -75.0:
            options['Atria major axis rotation degrees'] = -75.0
        elif options['Atria major axis rotation degrees'] > 75.0:
            options['Atria major axis rotation degrees'] = 75.0
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundVSeptum = options['Number of elements around ventricular septum']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundVSeptum
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpVSeptum = options['Number of elements up ventricular septum']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpVSeptum
        elementsCountUpRV = elementsCountUpVSeptum + 1
        elementsCountAroundRV = elementsCountAroundVSeptum + 2
        elementsCountAroundAtria = options['Number of elements around atria']
        elementsCountAtrialSeptum = options['Number of elements around atrial septum']
        lvOuterHeight = options['LV outer height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        lvApexThickness = options['LV apex thickness']
        rvHeight = options['RV inner height']
        rvArcAroundRadians = math.radians(options['RV arc around degrees'])
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        rvExtraCrossRadiusBase = options['RV extra cross radius base']
        vSeptumThickness = options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = options['Ventricular septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']
        aBaseInnerMajorMag = 0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = 0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aSeptumThickness = options['Atrial septum thickness']
        aBaseWallThickness = options['Atrial base wall thickness']
        aBaseSlopeRadians = math.radians(options['Atrial base slope degrees'])
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        fibrousRingThickness = options['Fibrous ring thickness']
        lvOutletInnerRadius = options['LV outlet inner diameter']*0.5
        lvOutletWallThickness = options['LV outlet wall thickness']
        lvOutletOuterRadius = lvOutletInnerRadius + lvOutletWallThickness
        rvOutletInnerRadius = options['RV outlet inner diameter']*0.5
        rvOutletWallThickness = options['RV outlet wall thickness']
        rvOutletOuterRadius = rvOutletInnerRadius + rvOutletWallThickness
        vOutletElementLength = options['Ventricles outlet element length']
        vOutletInclineRadians = math.radians(options['Ventricles outlet incline degrees'])
        vOutletSpacing = options['Ventricles outlet spacing']

        # generate heartventricles2 model to add base plane to
        annotationGroups = MeshType_3d_heartventricles2.generateBaseMesh(region, options)
        lvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left ventricle myocardium"))
        rvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right ventricle myocardium"))
        vSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interventricular septum"))
        conusArteriosusGroup = AnnotationGroup(region, get_heart_term("conus arteriosus"))
        annotationGroups += [ conusArteriosusGroup ]
        # av boundary nodes are put in left and right fibrous ring groups only so they can be found by heart1
        lFibrousRingGroup = AnnotationGroup(region, get_heart_term("left fibrous ring"))
        rFibrousRingGroup = AnnotationGroup(region, get_heart_term("right fibrous ring"))

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)
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
        # LV/RV outlet elements are linear through the wall, hence their nodes do not have D_DS3 parameters
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

        # node offsets for row, wall in LV, plus first LV node on inside top
        norl = elementsCountAroundLV
        nowl = 1 + elementsCountUpLV*norl
        nidl = nowl - norl + 1
        # make constants for end of septum and LV
        nsdl = nidl + elementsCountAroundVSeptum
        nedl = nidl + elementsCountAroundLV
        # node offsets for row, wall in RV, plus first RV node on inside top
        norr = elementsCountAroundRV - 1
        nowr = elementsCountUpRV*norr
        nidr = nowl*2 + 1 + nowr - norr
        #print('nidl',nidl,'nidr',nidr)

        # LV outlet
        elementsCountAroundOutlet = 6
        # GRC Set properly:
        defaultOutletScale3 = 0.5
        nidca = nidl + nowl + elementsCountAroundVSeptum - 1
        nidcb = nidr + elementsCountAroundVSeptum - 1
        #print('px nodes', nidca, nidcb)
        cache.setNode(nodes.findNodeByIdentifier(nidca))
        result, pxa = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        cache.setNode(nodes.findNodeByIdentifier(nidcb))
        result, pxb = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        px = [ 0.5*(pxa[c] + pxb[c]) for c in range(3) ]
        node = nodes.findNodeByIdentifier(nidl)
        cache.setNode(node)
        result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        node = nodes.findNodeByIdentifier(nidr)
        cache.setNode(node)
        result, bx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        bx = [ 0.5*(ax[i] + bx[i]) for i in range(2) ]
        bx.append(ax[2])
        ax = [ (bx[c] - px[c]) for c in range(3) ]
        ax = vector.normalise(ax)
        baseRotationRadians = math.atan2(ax[1], ax[0])
        # get crux location
        outletSpacingRadians = 0.25*math.pi  # GRC make option?
        outletSpacingHorizontal = vOutletSpacing*math.cos(outletSpacingRadians)
        outletSpacingVertical = vOutletSpacing*math.sin(outletSpacingRadians)
        cruxOffset = rvOutletOuterRadius + outletSpacingHorizontal + 2.0*lvOutletOuterRadius
        cx = [ (px[c] + ax[c]*cruxOffset) for c in range(3) ]

        #print('baseRotationRadians', baseRotationRadians)

        # calculate atria base inlet slope
        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)
        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)

        cosOutletInclineRadians = math.cos(vOutletInclineRadians)
        sinOutletInclineRadians = math.sin(vOutletInclineRadians)
        lvOutletCentre = [
            cx[0] - ax[0]*lvOutletOuterRadius,
            cx[1] - ax[1]*lvOutletOuterRadius,
            baseHeight + baseThickness - aBaseSlopeHeight + sinOutletInclineRadians*lvOutletOuterRadius ]

        radiansPerElementAroundOutlet = 2.0*math.pi/elementsCountAroundOutlet
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        lvOutletNodeId = []
        for n3 in range(2):
            radius = lvOutletInnerRadius if (n3 == 0) else lvOutletOuterRadius
            loAxis1 = [ radius*ax[c] for c in range(3) ]
            loAxis2 = [ -loAxis1[1]*cosOutletInclineRadians, loAxis1[0]*cosOutletInclineRadians, -radius*sinOutletInclineRadians ]
            loAxis3 = vector.crossproduct3(loAxis1, loAxis2)
            scale = vOutletElementLength/vector.magnitude(loAxis3)
            dx_ds2 = [ v*scale for v in loAxis3 ]
            outletNodeId = []
            for n1 in range(elementsCountAroundOutlet):
                radiansAround = n1*radiansPerElementAroundOutlet
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                outletScale3 = vOutletSpacing/radius if (n1 == 3) else defaultOutletScale3
                for c in range(3):
                    x[c] = lvOutletCentre[c] + loAxis1[c]*cosRadiansAround + loAxis2[c]*sinRadiansAround
                    dx_ds1[c] = radiansPerElementAroundOutlet*(loAxis1[c]*-sinRadiansAround + loAxis2[c]*cosRadiansAround)
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3 if (n3 == 0) else nodetemplate)
                outletNodeId.append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if n3 == 1:
                    dx_ds3 = [ outletScale3*(loAxis1[c]*cosRadiansAround + loAxis2[c]*sinRadiansAround) for c in range(3) ]
                    if n1 in [ 2, 4 ]:
                        if n1 == 2:
                            dx_ds3[2] = -dx_ds3[2]
                        else:
                            dx_ds3[2] = -2.0*dx_ds3[2]
                        scale = radiansPerElementAroundOutlet*rvOutletOuterRadius/vector.magnitude(dx_ds3)
                        dx_ds3 = [ d*scale for d in dx_ds3 ]
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
                    elif n1 == 3:
                        lvOutletOuterSpaceX = [ x[0], x[1], x[2] ]
                nodeIdentifier += 1
            lvOutletNodeId.append(outletNodeId)

        # RV outlet - for bicubic-linear tube connection
        outletCentreSpacing = lvOutletOuterRadius + outletSpacingHorizontal + rvOutletOuterRadius
        rvOutletCentre = [ (lvOutletCentre[c] - outletCentreSpacing*ax[c]) for c in range(3) ]
        # add outletSpacingVertical rotated by vOutletInclineRadians
        unitCrossX = vector.normalise([-ax[1], ax[0]])
        rvOutletCentre[0] -= outletSpacingVertical*sinOutletInclineRadians*unitCrossX[0]
        rvOutletCentre[1] -= outletSpacingVertical*sinOutletInclineRadians*unitCrossX[1]
        rvOutletCentre[2] += outletSpacingVertical*cosOutletInclineRadians

        rvOutletNodeId = []
        for n3 in range(2):
            radius = rvOutletInnerRadius if (n3 == 0) else rvOutletOuterRadius
            roAxis1 = [ radius*ax[c] for c in range(3) ]
            roAxis2 = [ -roAxis1[1]*cosOutletInclineRadians, roAxis1[0]*cosOutletInclineRadians, radius*sinOutletInclineRadians ]
            roAxis3 = vector.crossproduct3(roAxis1, roAxis2)
            scale = vOutletElementLength/vector.magnitude(roAxis3)
            dx_ds2 = [ v*scale for v in roAxis3 ]
            outletNodeId = []
            for n1 in range(elementsCountAroundOutlet):
                radiansAround = n1*radiansPerElementAroundOutlet
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                outletScale3 = vOutletSpacing/radius if (n1 == 0) else defaultOutletScale3
                for c in range(3):
                    x[c] = rvOutletCentre[c] + roAxis1[c]*cosRadiansAround + roAxis2[c]*sinRadiansAround
                    dx_ds1[c] = radiansPerElementAroundOutlet*(roAxis1[c]*-sinRadiansAround + roAxis2[c]*cosRadiansAround)
                hasDerivative3 = (n3 == 1) and (n1 in [ 0, 5 ])
                node = nodes.createNode(nodeIdentifier, nodetemplate if hasDerivative3 else nodetemplateLinearS3)
                outletNodeId.append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if hasDerivative3:
                    dx_ds3 = [ outletScale3*(roAxis1[c]*cosRadiansAround + roAxis2[c]*sinRadiansAround) for c in range(3) ]
                    if n1 in [ 1, 5 ]:
                        if n1 == 1:
                            dx_ds3[2] = -dx_ds3[2]
                        else:
                            dx_ds3[2] = 4.0*dx_ds3[2]
                            dx_ds3 = [ (dx_ds1[c] + dx_ds3[c]) for c in range(3) ]
                        mag3 = radiansPerElementAroundOutlet*rvOutletOuterRadius
                        scale = mag3/vector.magnitude(dx_ds3)
                        dx_ds3 = [ d*scale for d in dx_ds3 ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if n1 == 0:
                        rvOutletOuterSpaceX = [ x[0], x[1], x[2] ]
                nodeIdentifier += 1
            rvOutletNodeId.append(outletNodeId)

        # fix derivative 3 between lv, rv outlets
        cache.setNode(nodes.findNodeByIdentifier(lvOutletNodeId[1][3]))
        dx_ds3 = [ (rvOutletOuterSpaceX[c] - lvOutletOuterSpaceX[c]) for c in range(3)]
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        cache.setNode(nodes.findNodeByIdentifier(rvOutletNodeId[1][0]))
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ -d for d in dx_ds3])

        # Atria nodes

        # Previously computed these:
        #aInnerMajorMag = 0.9*(lvOuterRadius - lvFreeWallThickness - 0.5*vSeptumBaseRadialDisplacement)
        #aInnerMinorMag = 1.0*(lvOuterRadius - lvFreeWallThickness - lvOutletOuterRadius)
        aInnerMajorMag = aBaseInnerMajorMag
        aInnerMinorMag = aBaseInnerMinorMag
        #print('inner mag major', aInnerMajorMag, 'minor', aInnerMinorMag)
        aOuterMajorMag = aInnerMajorMag + aBaseSlopeLength
        aOuterMinorMag = aInnerMinorMag + aBaseSlopeLength
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

        laCruxLeftRadians = getEllipseRadiansToX(axMod, bxMod, cruxLeftModX - laCentreModX, math.pi*0.5)
        laCentreModY = cruxLeftModY - ayMod*math.cos(laCruxLeftRadians) - byMod*math.sin(laCruxLeftRadians)

        #print('laCentreMod', laCentreModX, laCentreModY)
        laCentreX = cruxCentre[0] + laCentreModX*cosRotRadians + laCentreModY*-sinRotRadians
        laCentreY = cruxCentre[1] + laCentreModX*sinRotRadians + laCentreModY*cosRotRadians

        raCruxLeftRadians = -laCruxLeftRadians
        raCentreX = cruxCentre[0] - laCentreModX*cosRotRadians + laCentreModY*-sinRotRadians
        raCentreY = cruxCentre[1] - laCentreModX*sinRotRadians + laCentreModY*cosRotRadians

        aCentreOuterZ = baseHeight + baseThickness
        aCentreInnerZ = aCentreOuterZ - aBaseSlopeHeight

        atrialPerimeterLength = getApproximateEllipsePerimeter(aOuterMajorMag, aOuterMinorMag)
        atrialSeptumCentreToCruxLeftLength = getEllipseArcLength(aOuterMajorMag, aOuterMinorMag, laSeptumRadians, laCruxLeftRadians)
        atrialSeptumElementLength = atrialSeptumCentreToCruxLeftLength/(1.0 + elementsCountAtrialSeptum*0.5)
        atrialFreeWallElementLength = (atrialPerimeterLength - atrialSeptumElementLength*(elementsCountAtrialSeptum + 2)) \
            / (elementsCountAroundAtria - elementsCountAtrialSeptum - 2)
        atrialTransitionElementLength = 0.5*(atrialSeptumElementLength + atrialFreeWallElementLength)

        #print('lengths septum', atrialSeptumElementLength, 'transition', atrialTransitionElementLength, 'freewall', atrialFreeWallElementLength)
        #print('total length', atrialSeptumElementLength*(elementsCountAtrialSeptum + 1) + 2*atrialTransitionElementLength \
        #    + (elementsCountAroundAtria - elementsCountAtrialSeptum - 3)*atrialFreeWallElementLength, 'vs.', atrialPerimeterLength)

        laRadians = []
        laOuterDerivatives = []
        radiansAround = laSeptumRadians
        if (elementsCountAtrialSeptum % 2) == 1:
            radiansAround = updateEllipseAngleByArcLength(aOuterMajorMag, aOuterMinorMag, radiansAround, 0.5*atrialSeptumElementLength)
        outerDerivative = atrialSeptumElementLength
        lan1CruxLimit = elementsCountAtrialSeptum//2 + 1
        lan1SeptumLimit = elementsCountAroundAtria - (elementsCountAtrialSeptum + 1)//2 - 1
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

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = laRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    laCentreX + cosRadiansAround*laInnerMajor[0] + sinRadiansAround*laInnerMinor[0],
                    laCentreY + cosRadiansAround*laInnerMajor[1] + sinRadiansAround*laInnerMinor[1],
                    aCentreInnerZ ]
                outer = [
                    laCentreX + cosRadiansAround*laOuterMajor[0] + sinRadiansAround*laOuterMinor[0],
                    laCentreY + cosRadiansAround*laOuterMajor[1] + sinRadiansAround*laOuterMinor[1],
                    aCentreOuterZ ]

                if (n3 == 1) and ((n1 <= lan1CruxLimit) or (n1 > (lan1SeptumLimit + 2))):
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
                if (n1 < lan1CruxLimit) or (n1 > lan1SeptumLimit):
                    dx_ds2 = [ 0.0, 0.0, aCentreInnerZ ]
                else:
                    dx_ds2 = [
                        dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                        dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                        dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                    if n1 == (lan1CruxLimit + 1):
                        # make derivative 2 on lv crest larger and less inclined
                        dx_ds2[2] *= 0.5 if (n3 == 0) else 0.25
                        mag2 = 1.5*(baseHeight + baseThickness)
                    else:
                        # GRC check scaling here:
                        mag2 = inner[2]
                    scale2 = mag2/vector.magnitude(dx_ds2)
                    dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if False:
            # show axes of left atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, aCentreInnerZ ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laInnerMajor[0], laInnerMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laInnerMinor[0], laInnerMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, aCentreInnerZ ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ laCentreX, laCentreY, aCentreOuterZ ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ laOuterMajor[0], laOuterMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ laOuterMinor[0], laOuterMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, aCentreOuterZ ])
            nodeIdentifier += 1

        ran1SeptumLimit = elementsCountAtrialSeptum//2
        ran1CruxLimit = elementsCountAroundAtria - ran1SeptumLimit - 1
        raNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        raNodeId[1][0] = laNodeId[0][0]
        raNodeId[1][-2] = lvOutletNodeId[1][1]
        raNodeId[1][-1] = lvOutletNodeId[1][0]

        for n3 in range(2):
            for n1 in range(elementsCountAroundAtria):
                radiansAround = raRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                inner = [
                    raCentreX + cosRadiansAround*raInnerMajor[0] + sinRadiansAround*raInnerMinor[0],
                    raCentreY + cosRadiansAround*raInnerMajor[1] + sinRadiansAround*raInnerMinor[1],
                    aCentreInnerZ ]
                outer = [
                    raCentreX + cosRadiansAround*raOuterMajor[0] + sinRadiansAround*raOuterMinor[0],
                    raCentreY + cosRadiansAround*raOuterMajor[1] + sinRadiansAround*raOuterMinor[1],
                    aCentreOuterZ ]

                if (n3 == 1) and ((n1 < ran1SeptumLimit) or (n1 >= ran1CruxLimit)):
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
                if (n1 <= ran1SeptumLimit) or (n1 >= ran1CruxLimit):
                    dx_ds2 = [ 0.0, 0.0, aCentreInnerZ ]
                else:
                    dx_ds2 = [
                        dx_ds3[1]*dx_ds1[2] - dx_ds3[2]*dx_ds1[1],
                        dx_ds3[2]*dx_ds1[0] - dx_ds3[0]*dx_ds1[2],
                        dx_ds3[0]*dx_ds1[1] - dx_ds3[1]*dx_ds1[0] ]
                    if n1 == (ran1CruxLimit - 1):
                        # make derivative 2 on sv crest larger and less inclined
                        dx_ds2[2] *= 0.5 if (n3 == 0) else 0.25
                        mag2 = 1.5*(baseHeight + baseThickness)
                    else:
                        # GRC check scaling here:
                        mag2 = inner[2]
                    scale2 = mag2/vector.magnitude(dx_ds2)
                    dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        if False:
            # show axes of right atrium
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, aCentreInnerZ ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raInnerMajor[0], raInnerMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raInnerMinor[0], raInnerMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, aCentreInnerZ ])
            nodeIdentifier += 1

            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ raCentreX, raCentreY, aCentreOuterZ ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ raOuterMajor[0], raOuterMajor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ raOuterMinor[0], raOuterMinor[1], 0.0 ])
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, aCentreOuterZ ])
            nodeIdentifier += 1

        laNodeId[1][0] = raNodeId[0][0]
        laNodeId[1][1] = lvOutletNodeId[1][ 0]
        laNodeId[1][2] = lvOutletNodeId[1][-1]

        #print('laNodeId[0]', laNodeId[0])
        #print('laNodeId[1]', laNodeId[1])
        #print('raNodeId[0]', raNodeId[0])
        #print('raNodeId[1]', raNodeId[1])

        # compute dx_ds3 around atria from differences of nodes
        for i in range(2):
            aNodeId = laNodeId if (i == 0) else raNodeId
            for n1 in range(elementsCountAroundAtria):
                nid2 = aNodeId[1][n1]
                node2 = nodes.findNodeByIdentifier(nid2)
                cache.setNode(node2)
                result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
                nid1 = aNodeId[0][n1]
                node1 = nodes.findNodeByIdentifier(nid1)
                cache.setNode(node1)
                result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
                dx_ds3 = [ x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] ]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                if (i == 1) and ((n1 == 0) or (nid2 == cruxCentreNodeId)):
                    continue
                if nid2 in [ cruxLeftNodeId, cruxRightNodeId ]:
                    dx_ds3 = [ -d for d in dx_ds3 ]
                cache.setNode(node2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # fix crux centre dx_ds3:
        cache.setNode(nodes.findNodeByIdentifier(laNodeId[0][1]))
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        cache.setNode(nodes.findNodeByIdentifier(raNodeId[0][-1]))
        result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        cache.setNode(nodes.findNodeByIdentifier(cruxCentreNodeId))
        result, xc = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        d1 = [ (x1[c] - xc[c]) for c in range(3) ]
        d2 = [ (x2[c] - xc[c]) for c in range(3) ]
        dx_ds3 = [ d1[0] + d2[0], d1[1] + d2[1], d1[2] + d2[2] ]
        scale = vector.magnitude(d1)/vector.magnitude(dx_ds3)
        dx_ds3 = [ d*scale for d in dx_ds3 ]
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3 )

        # fix la-lo bridge dx_ds2 on inside la
        cache.setNode(nodes.findNodeByIdentifier(laNodeId[0][2]))
        result, dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3 )
        result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3 )
        dx_ds2 = vector.crossproduct3(dx_ds3, dx_ds1)
        scale2 = 0.5*vector.magnitude(dx_ds3)/vector.magnitude(dx_ds2)
        dx_ds2 = [ scale2*d for d in dx_ds2 ]
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2 )

        # create node on top of lv 'crest'
        nida = nsdl + nowl
        nidb = lvOutletNodeId[1][-2]
        #print('lv crest interpolated from nodes', nida, nidb)
        node1 = nodes.findNodeByIdentifier(nida)
        node2 = nodes.findNodeByIdentifier(nidb)
        x, dx_ds2, dx_ds1, dx_ds3 = interpolateNodesCubicHermite(cache, coordinates, 0.5, baseThickness, \
            node1, Node.VALUE_LABEL_D_DS2,  1.0, Node.VALUE_LABEL_D_DS1, 1.0, \
            node2, Node.VALUE_LABEL_D_DS3, -1.0, Node.VALUE_LABEL_D_DS1, 1.0)

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        lv_crest_nid1 = nodeIdentifier
        nodeIdentifier += 1

        # create nodes on bottom and top of rv supraventricular crest
        nida = nidr + nowr + 4
        nidb = lvOutletNodeId[1][2]
        #print('rv crest interpolated from nodes', nida, nidb)
        node = nodes.findNodeByIdentifier(nida)
        cache.setNode(node)
        result, xa = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        result, d1a = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3 )
        result, d2a = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3 )
        node = nodes.findNodeByIdentifier(nidb)
        cache.setNode(node)
        result, xb = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        result, d1b = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3 )
        result, d2b = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3 )
        d2b = [ -2.0*d for d in d2b ]
        scale = 4.0*(baseHeight + baseThickness)/vector.magnitude(d2a)
        d2a = [ scale*d for d in d2a ]
        xi = 0.5
        xr = 1.0 - xi
        x = interp.interpolateCubicHermite(xa, d2a, xb, d2b, xi)
        dx_ds1 = [ (xr*d1a[c] + xi*d1b[c]) for c in range(3) ]
        dx_ds2 = interp.interpolateCubicHermiteDerivative(xa, d2a, xb, d2b, xi)
        dx_ds2 = [ xr*d for d in dx_ds2 ]
        radialVector = vector.normalise(vector.crossproduct3(dx_ds1, dx_ds2))
        dx_ds3 = [ baseThickness*d for d in radialVector ]

        x_inner = [ (x[c] - dx_ds3[c]) for c in range(3) ]
        curvatureScale = 1.0 - baseThickness*interp.getCubicHermiteCurvature(xa, d2a, x, dx_ds2, radialVector, 1.0)
        dx_ds2_inner = [ curvatureScale*d for d in dx_ds2 ]

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x_inner)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2_inner)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        rv_crest_nid1 = nodeIdentifier
        nodeIdentifier += 1

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        rv_crest_nid2 = nodeIdentifier
        nodeIdentifier += 1

        # create node on bottom of lv 'bridge' between la and lo
        nida = laNodeId[0][2]
        nidb = lvOutletNodeId[0][-1]
        #print('lv bridge interpolated from nodes', nida, nidb)
        node1 = nodes.findNodeByIdentifier(nida)
        node2 = nodes.findNodeByIdentifier(nidb)
        x, dx_ds2, dx_ds1, dx_ds3 = interpolateNodesCubicHermite(cache, coordinates, 0.4, lvOutletWallThickness, \
            node1, Node.VALUE_LABEL_D_DS2, -1.0, Node.VALUE_LABEL_D_DS1, -1.0, \
            node2, Node.VALUE_LABEL_D_DS2,  1.0, Node.VALUE_LABEL_D_DS1,  1.0)
        # dx_ds1 needs to be larger
        dx_ds1 = [ 2.0*d  for d in dx_ds1 ]

        node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
        lv_bridge_nid1 = nodeIdentifier
        nodeIdentifier += 1

        # create nodes on top of fibrous ring, between ventricles and atria
        tlaNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        traNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        for n3 in range(2):

            for i in range(2):
                if i == 0:
                    baNodeId = laNodeId
                    taNodeId = tlaNodeId
                else:
                    baNodeId = raNodeId
                    taNodeId = traNodeId

                for n1 in range(elementsCountAroundAtria):
                    if (n3 == 1) and \
                        (((i == 0) and ((n1 < (lan1CruxLimit - 1)) or (n1 > (lan1SeptumLimit + 2)))) or \
                         ((i == 1) and ((n1 < ran1SeptumLimit) or (n1 > ran1CruxLimit)))):
                        # get node from other side below
                        continue

                    node = nodes.findNodeByIdentifier(baNodeId[n3][n1])
                    cache.setNode(node)
                    result, x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
                    result, dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3 )
                    result, dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3 )
                    result, dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3 )
                    x[2] += fibrousRingThickness
                    if (n3 == 1) and (((i == 0) and ((n1 == 1) or (n1 == 2))) or ((i == 1) and (n1 == (elementsCountAroundAtria - 2)))):
                        dx_ds1 = [ -d for d in dx_ds1 ]
                        dx_ds3 = [ -d for d in dx_ds3 ]

                    taNodeId[n3][n1] = nodeIdentifier
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    nodeIdentifier += 1

        # tie up overlaps
        tlaNodeId[1][ 0] = traNodeId[0][0]
        traNodeId[1][-1] = tlaNodeId[1][1]
        traNodeId[1][ 0] = tlaNodeId[0][0]

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)
        conusArteriosusMeshGroup = conusArteriosusGroup.getMeshGroup(mesh)
        lFibrousRingMeshGroup = lFibrousRingGroup.getMeshGroup(mesh)
        rFibrousRingMeshGroup = rFibrousRingGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # LV base elements
        for e in range(19):
            eft1 = eft
            nids = None
            meshGroups = [ lvMeshGroup ]

            if e == 0:
                # 8-node atrial septum element 1
                nids = [ nidl +  0, nidl +  1, laNodeId[0][-1],   laNodeId[0][ 0], nidr        +  0, nidl + nowl +  1, raNodeId[0][ 1], laNodeId[1][ 0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ] )
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                scaleEftNodeValueLabels(eft1, [ 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 1:
                # 8-node atrial septum element 2
                nids = [ nidl +  1, nidl +  2, laNodeId[0][ 0],   laNodeId[0][ 1], nidl + nowl +  1, nidl + nowl +  2, laNodeId[1][ 0], raNodeId[0][-1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                scaleEftNodeValueLabels(eft1, [ 7 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 2:
                # 6-node crux element, multiple collapses
                nids = [ nidl +  2, lvOutletNodeId[0][0], laNodeId[0][ 1], lvOutletNodeId[1][0], nidl + nowl +  2, raNodeId[0][-1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                # remap parameters before collapsing nodes
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, 2, 4, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                # must set DS3 before DS1 to allow latter to equal former
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, [])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, [])
                #remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                #remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                # must set DS3 before DS1 to allow latter to equal former
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                ln_map = [ 1, 2, 3, 4, 5, 4, 6, 4 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 3:
                # 6 node collapsed vs-ra shim element
                nids = [ lvOutletNodeId[1][0], lvOutletNodeId[1][1], nidl + nowl + 2, nidl + nowl + 3, raNodeId[0][-1], raNodeId[0][-2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                #remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                #remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                ln_map = [ 1, 2, 1, 2, 3, 4, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e <= 6:
                # 8-node ventricular septum elements
                n = e - 4
                nids = [ nidl        + n + 2, nidl        + n + 3, lvOutletNodeId[0][n], lvOutletNodeId[0][n + 1], nidl + nowl + n + 2, nidl + nowl + n + 3, lvOutletNodeId[1][n], lvOutletNodeId[1][n + 1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                if e == 4:
                    # the following will not be overridden by following remappings on the same node x value
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 7:
                # 8-node ventricular septum element past lo-ro junction
                nids = [ nidl + 5, nidl + 6, lvOutletNodeId[0][3], lvOutletNodeId[0][4], nidl + nowl + 5, nidl + nowl + 6, lvOutletNodeId[1][3], lvOutletNodeId[1][4] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 8:
                # 7-node collapsed final ventricular septum element by rv junction past lo-ro junction, triangle on top & outside
                nids = [ nidl + 6, nidl + 7, lvOutletNodeId[0][4], lv_crest_nid1, nidl + nowl + 6, nidr + norr - 1, lvOutletNodeId[1][4] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 7, 4 ]
                remapEftLocalNodes(eft1, 7, ln_map)
                meshGroups += [ rvMeshGroup, vSeptumMeshGroup ]
            elif e == 9:
                # 5-node multiple collapse element junction of rv, septum at ro
                nids = [ nidr + norr - 1, nsdl, nidr + nowr + norr - 1, nsdl + nowl, lv_crest_nid1 ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                ln_map = [ 1, 2, 1, 2, 3, 4, 5, 5 ]
                remapEftLocalNodes(eft1, 5, ln_map)
                meshGroups += [ rvMeshGroup ]
            elif e == 10:
                # 7-node collapsed LV free wall element 1, by lv crest and la
                nids = [ nsdl, nsdl + 1, laNodeId[0][3], nsdl + nowl, nsdl + nowl + 1, lv_crest_nid1, laNodeId[1][3] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 1, 3, 4, 5, 6, 7]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 11:
                # regular LV free wall element 2, by la
                nids = [ nedl -  4, nedl -  3, laNodeId[0][-5],   laNodeId[0][-4], nedl + nowl -  4, nedl + nowl -  3, laNodeId[1][-5], laNodeId[1][-4] ]
            elif e == 12:
                # regular LV free wall element 3, by la
                nids = [ nedl -  3, nedl -  2, laNodeId[0][-4],   laNodeId[0][-3], nedl + nowl -  3, nedl + nowl -  2, laNodeId[1][-4], laNodeId[1][-3] ]
            elif e == 13:
                # regular LV free wall element 4, by la
                nids = [ nedl -  2, nedl -  1, laNodeId[0][-3],   laNodeId[0][-2], nedl + nowl -  2, nedl + nowl -  1, laNodeId[1][-3], laNodeId[1][-2] ]
            elif e == 14:
                # regular LV free wall element 5, by la
                nids = [ nedl -  1, nidl +  0, laNodeId[0][-2],   laNodeId[0][-1], nedl + nowl -  1, nidl + nowl +  0, laNodeId[1][-2], laNodeId[1][-1] ]
            elif e == 15:
                # 8-node on bottom of lv crest between la and lo
                nids = [ nidl + 6, nidl + 7, lv_bridge_nid1, laNodeId[0][3], lvOutletNodeId[0][4], lv_crest_nid1, lvOutletNodeId[0][5], laNodeId[1][3] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
            elif e == 16:
                # 6-node wedge element on top of lv crest beween la and lo
                nids = [ lv_crest_nid1, laNodeId[1][3], lvOutletNodeId[0][4], lvOutletNodeId[0][5], lvOutletNodeId[1][4], lvOutletNodeId[1][5] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [])
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 17:
                # 6-node wedge element on lv bridge beween la and lo
                nids = [ laNodeId[0][3], laNodeId[0][2], lv_bridge_nid1, laNodeId[1][3], lvOutletNodeId[1][5], lvOutletNodeId[0][5] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, 6, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1, 2, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 1, 3, 4, 5, 4, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 18:
                # 8-node element on lv bridge between la and lo, connecting to septum
                nids = [ lv_bridge_nid1, nidl + elementsCountAtrialSeptum, lvOutletNodeId[0][5], lvOutletNodeId[0][0], \
                         laNodeId[0][2], laNodeId[0][1], lvOutletNodeId[1][5], lvOutletNodeId[1][0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ]) # , ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary, to swap with D_DS2
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # swap from above
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 1
            #print('create element lv', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # RV base elements
        for e in range(15):
            eft1 = eft
            nids = None
            meshGroups = [ rvMeshGroup ]

            if e == 0:
                # lv-rv junction
                nids = [ nidl + 0, nidr + 0, laNodeId[0][-1], raNodeId[0][ 1], nidl + nowl + 0, nidr + nowr + 0, laNodeId[1][-1], raNodeId[1][ 1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                meshGroups += [ lvMeshGroup ]
            elif e == 1:
                # regular rv free wall element 1
                nids = [ nidr        + 0, nidr + 1, raNodeId[0][ 1], raNodeId[0][ 2], nidr + nowr + 0, nidr + nowr + 1, raNodeId[1][ 1], raNodeId[1][ 2] ]
            elif e == 2:
                # regular rv free wall element 2
                nids = [ nidr        + 1, nidr + 2, raNodeId[0][ 2], raNodeId[0][ 3], nidr + nowr + 1, nidr + nowr + 2, raNodeId[1][ 2], raNodeId[1][ 3] ]
            elif e == 3:
                # regular rv free wall element 3
                nids = [ nidr        + 2, nidr + 3, raNodeId[0][ 3], raNodeId[0][ 4], nidr + nowr + 2, nidr + nowr + 3, raNodeId[1][ 3], raNodeId[1][ 4] ]
            elif e == 4:
                # supraventricular crest outer 1
                nids = [ nidr        + 3, nidr + 4, raNodeId[0][ 4],      rv_crest_nid1, nidr + nowr + 3, nidr + nowr + 4, raNodeId[1][ 4],      rv_crest_nid2 ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
            elif e == 5:
                # supraventricular crest outer 2, outer infundibulum 1
                nids = [ nidr + 4, nidr + 5, rv_crest_nid1, rvOutletNodeId[0][2], nidr + nowr + 4, nidr + nowr + 5, rv_crest_nid2, rvOutletNodeId[1][2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 6:
                # outer infundibulum 2
                nids = [ nidr + 5, nidr + 6, rvOutletNodeId[0][2], rvOutletNodeId[0][3], nidr + nowr + 5, nidr + nowr + 6, rvOutletNodeId[1][2], rvOutletNodeId[1][3] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 7:
                # outer infundibulum 3
                nids = [ nidr        + 6, nidr + 7, rvOutletNodeId[0][3], rvOutletNodeId[0][4], nidr + nowr + 6, nidr + nowr + 7, rvOutletNodeId[1][3], rvOutletNodeId[1][4] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 8:
                # 7-node collapsed rv crest inner 1, by RA-LV outlet junction
                nids = [ raNodeId[0][-3], nidl + nowl +  4, raNodeId[0][-2], nidl + nowl + 3, raNodeId[1][-3], lvOutletNodeId[1][2], lvOutletNodeId[1][1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 7, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 9:
                # 8-node rv crest inner 2
                nids = [ raNodeId[0][-4], rv_crest_nid1, raNodeId[0][-3], nidl + nowl +  4, raNodeId[1][-4], rv_crest_nid2, raNodeId[1][-3], lvOutletNodeId[1][2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] )])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
            elif e == 10:
                # 8-node wedge rv crest inner 3
                nids = [ rv_crest_nid1, rvOutletNodeId[0][2], nidl + nowl +  4, rvOutletNodeId[0][1], rv_crest_nid2, rvOutletNodeId[1][2], lvOutletNodeId[1][2], rvOutletNodeId[1][1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 2, 3, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 11:
                # 8-node rv crest inner 4 by rv outlet
                nids = [ nidl + nowl +  4, rvOutletNodeId[0][1], nidl + nowl +  5, rvOutletNodeId[0][0], lvOutletNodeId[1][2], rvOutletNodeId[1][1], lvOutletNodeId[1][3], rvOutletNodeId[1][0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary, to swap with D_DS2
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # swap from above
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary, to swap with D_DS2
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 12:
                # 8-node ro-septum 1 past lo-ro junction
                nids = [ nidl + nowl + 6, nidl + nowl + 5, rvOutletNodeId[0][-1], rvOutletNodeId[0][0], \
                         lvOutletNodeId[1][4], lvOutletNodeId[1][3], rvOutletNodeId[1][-1], rvOutletNodeId[1][0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary, to swap with D_DS2
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # swap from above
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 13:
                # 8-node ro-septum 2 past lo-ro junction
                nids = [ nidr + elementsCountAroundVSeptum, nidl + nowl + 6, rvOutletNodeId[0][-2], rvOutletNodeId[0][-1],
                         lv_crest_nid1, lvOutletNodeId[1][4], rvOutletNodeId[1][-2], rvOutletNodeId[1][-1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 14:
                # 5-node collapsed 'open corner of book' ro closure
                nids = [ nidr + elementsCountAroundVSeptum, rvOutletNodeId[0][-2], nidr + nowr + elementsCountAroundVSeptum, lv_crest_nid1, rvOutletNodeId[1][-2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 1, 2, 2, 3, 4, 5, 5 ]
                remapEftLocalNodes(eft1, 5, ln_map)
                meshGroups += [ conusArteriosusMeshGroup ]

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 1
            #print('create element rv', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        #print('blaNodeId', laNodeId)
        #print('tlaNodeId', tlaNodeId)
        #print('braNodeId', raNodeId)
        #print('traNodeId', traNodeId)

        # fibrous ring
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives, linearAxis = 2, d_ds1 = Node.VALUE_LABEL_D_DS1, d_ds2 = Node.VALUE_LABEL_D_DS3)
        eftFibrousRing = bicubichermitelinear.createEftBasic()
        for i in range(2):
            if i == 0:
                baNodeId = laNodeId
                taNodeId = tlaNodeId
                meshGroupsSide = [ lFibrousRingMeshGroup ]
            else:
                baNodeId = raNodeId
                taNodeId = traNodeId
                meshGroupsSide = [ rFibrousRingMeshGroup ]
            for e1 in range(elementsCountAroundAtria + 2):

                if (i == 1) and (e1 < 4):
                    continue

                if e1 < 4:
                    meshGroups = [ lFibrousRingMeshGroup, rFibrousRingMeshGroup ]
                else:
                    meshGroups = meshGroupsSide

                eft1 = eftFibrousRing
                if (e1 == 0) or (e1 == 3):
                    eft1 = bicubichermitelinear.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    if e1 == 0:
                        nids = [ laNodeId[0][-1], raNodeId[0][ 1], tlaNodeId[0][-1], traNodeId[0][ 1], \
                                 laNodeId[1][-1], raNodeId[1][ 1], tlaNodeId[1][-1], traNodeId[1][ 1] ]
                    else:
                        nids = [ raNodeId[0][-1], laNodeId[0][ 1], traNodeId[0][-1], tlaNodeId[0][ 1], \
                                 laNodeId[1][ 1], tlaNodeId[1][ 1] ]
                        # 6 node fibrous ring crux, collapse xi1 on xi3 == 1
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                elif e1 == 1:
                    nids = [ laNodeId[0][-1], laNodeId[0][ 0], tlaNodeId[0][-1], tlaNodeId[0][ 0], \
                                raNodeId[0][ 1], raNodeId[0][ 0], traNodeId[0][ 1], traNodeId[0][ 0] ]
                    eft1 = bicubichermitelinear.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                elif e1 == 2:
                    nids = [ laNodeId[0][ 0], laNodeId[0][ 1], tlaNodeId[0][ 0], tlaNodeId[0][ 1], \
                             raNodeId[0][ 0], raNodeId[0][-1], traNodeId[0][ 0], traNodeId[0][-1] ]
                    eft1 = bicubichermitelinear.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                else:
                    ea = e1 - 3
                    eb = ea + 1 - elementsCountAroundAtria
                    nids = [ baNodeId[0][ea], baNodeId[0][eb], taNodeId[0][ea], taNodeId[0][eb], \
                             baNodeId[1][ea], baNodeId[1][eb], taNodeId[1][ea], taNodeId[1][eb] ]
                    if ((i == 0) and ((e1 == 4) or (e1 == 5))) or ((i == 1) and (e1 >= elementsCountAroundAtria)):
                        eft1 = bicubichermitelinear.createEftBasic()
                        setEftScaleFactorIds(eft1, [1], [])
                        if e1 == 4:
                            scaleEftNodeValueLabels(eft1, [ 6 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        elif e1 == 5:
                            scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            #remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                        elif e1 == elementsCountAroundAtria:
                            scaleEftNodeValueLabels(eft1, [ 6 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            #remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:  # e1 == (elementsCountAroundAtria + 1)
                            scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = 1
                #print(i, e1, 'create element fr', elementIdentifier, result, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']
        MeshType_3d_heartventricles2.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startBaseLvElementIdentifier = element.getIdentifier()
        startBaseRvElementIdentifier = startBaseLvElementIdentifier + 19
        limitBaseRvElementIdentifier = startBaseRvElementIdentifier + 15
        limitFibrousRingElementIdentifier = limitBaseRvElementIdentifier + 16
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            elementId = element.getIdentifier()
            if elementId < startBaseRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < limitBaseRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughRVWall
            else:  # elementId < limitFibrousRingElementIdentifier:
                numberInXi2 = 1
                numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitFibrousRingElementIdentifier - 1):
                return  # finish on last so can continue in full heart mesh
            element = meshrefinement._sourceElementiterator.next()
