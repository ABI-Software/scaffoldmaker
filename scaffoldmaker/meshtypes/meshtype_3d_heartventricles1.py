"""
Generates 3-D Left and Right ventricles mesh starting from modified sphere shell mesh.
"""

from __future__ import division
import math
from scaffoldmaker.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventricles1:
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Ventricles 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around' : 12,
            'Number of elements up' : 4,
            'Number of elements through LV wall' : 1,
            'Number of elements across septum' : 5,
            'Number of elements below septum' : 2,
            'LV wall thickness' : 0.15,
            'LV wall thickness ratio apex' : 0.5,
            'LV wall thickness ratio base': 0.8,
            'LV base flatten ratio': 0.75,
            'LV base flatten angle degrees': 0.0,
            'RV free wall thickness' : 0.05,
            'RV width' : 0.2,
            'Length ratio' : 2.0,
            'Element length ratio equator/apex' : 1.0,
            'Septum arc angle degrees' : 125.0,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements surface' : 1,
            'Refine number of elements through LV wall' : 1,
            'Refine number of elements through RV wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements up',
            'Number of elements through LV wall',
            'Number of elements across septum',
            'Number of elements below septum',
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'LV wall thickness ratio base',
            'LV base flatten ratio',
            'LV base flatten angle degrees',
            'RV free wall thickness',
            'RV width',
            'Length ratio',
            'Element length ratio equator/apex',
            'Septum arc angle degrees',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements through LV wall',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through RV wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 4:
            options['Number of elements up'] = 4
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        if options['Number of elements across septum'] < 2:
            options['Number of elements across septum'] = 2
        elif options['Number of elements across septum'] > (options['Number of elements around'] - 2):
            options['Number of elements across septum'] = options['Number of elements around'] - 2
        if options['Number of elements below septum'] < 2:
            options['Number of elements below septum'] = 2
        elif options['Number of elements below septum'] > (options['Number of elements up'] - 1):
            options['Number of elements below septum'] = options['Number of elements up'] - 1
        if options['LV wall thickness'] < 0.0:
            options['LV wall thickness'] = 0.0
        elif options['LV wall thickness'] > 0.5:
            options['LV wall thickness'] = 0.5
        if options['LV wall thickness ratio apex'] < 0.0:
            options['LV wall thickness ratio apex'] = 0.0
        if options['LV wall thickness ratio base'] < 0.0:
            options['LV wall thickness ratio base'] = 0.0
        if options['LV base flatten ratio'] < 0.0:
            options['LV base flatten ratio'] = 0.0
        if options['RV free wall thickness'] < 0.0:
            options['RV free wall thickness'] = 0.0
        if options['RV width'] < 0.0:
            options['RV width'] = 0.0
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6
        if options['Septum arc angle degrees'] < 45.0:
            options['Septum arc angle degrees'] = 45.0
        elif options['Septum arc angle degrees'] > 270.0:
            options['Septum arc angle degrees'] = 270.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountAcrossSeptum = options['Number of elements across septum']
        elementsCountBelowSeptum = options['Number of elements below septum']
        LVWallThickness = options['LV wall thickness']
        LVWallThicknessRatioApex = options['LV wall thickness ratio apex']
        LVWallThicknessRatioBase = options['LV wall thickness ratio base']
        LVBaseFlattenRatio = options['LV base flatten ratio']
        LVBaseFlattenAngleRadians = options['LV base flatten angle degrees']*math.pi/180.0

        lengthRatio = options['Length ratio']
        RVWidthTop = options['RV width']
        RVFreeWallThickness = options['RV free wall thickness']
        septumArcAngleRadians = options['Septum arc angle degrees']*math.pi/180.0
        useCrossDerivatives = options['Use cross derivatives']

        # generate a half sphere shell which will be edited
        sphereShellOptions = MeshType_3d_sphereshell1.getDefaultOptions()
        sphereShellOptions['Number of elements up'] = elementsCountUp*2
        sphereShellOptions['Number of elements around'] = elementsCountAround
        sphereShellOptions['Number of elements through wall'] = elementsCountThroughLVWall
        sphereShellOptions['Exclude top rows'] = elementsCountUp
        sphereShellOptions['Wall thickness'] = LVWallThickness
        sphereShellOptions['Wall thickness ratio apex'] = LVWallThicknessRatioApex
        sphereShellOptions['Length ratio'] = lengthRatio
        sphereShellOptions['Element length ratio equator/apex'] = options['Element length ratio equator/apex']
        MeshType_3d_sphereshell1.generateBaseMesh(region, sphereShellOptions)

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fm.findMeshByDimension(3)

        cache = fm.createFieldcache()

        nor = elementsCountAround
        now = 1 + elementsCountUp*nor

        # Resize elements around LV to get desired septum arc angle
        radiansPerElementOrig = 2.0*math.pi/elementsCountAround
        radiansPerElementSeptum = septumArcAngleRadians/elementsCountAcrossSeptum
        # want LV-RV 'transition' elements to be mean size of Septum and LV FreeWall elements
        radiansRemaining = 2.0*math.pi - septumArcAngleRadians
        elementsCountAroundLVFreeWall = elementsCountAround - elementsCountAcrossSeptum - 2
        radiansPerElementLVFreeWall = (radiansRemaining - radiansPerElementSeptum) / (elementsCountAroundLVFreeWall + 1)
        radiansPerElementTransition = 0.5*(radiansPerElementSeptum + radiansPerElementLVFreeWall)
        #print('Element size ratio LVFreeWall / Septum', radiansPerElementLVFreeWall/radiansPerElementSeptum)
        xyz_scale = fm.createFieldConstant([1.0, 1.0, lengthRatio/2.0 - LVWallThickness*LVWallThicknessRatioApex])
        coordinates_scale = fm.createFieldMultiply(coordinates, xyz_scale)
        sp = fm.createFieldCoordinateTransformation(coordinates_scale)
        sp.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_SPHERICAL_POLAR)
        r = fm.createFieldComponent(sp, 1)
        theta = fm.createFieldComponent(sp, 2)
        phi = fm.createFieldComponent(sp, 3)
        zero = fm.createFieldConstant([0.0])
        one = fm.createFieldConstant([1.0])  # also used as 'true'
        two_pi = fm.createFieldConstant([2.0*math.pi])
        radians_per_element_orig = fm.createFieldConstant([radiansPerElementOrig])
        half_radians_per_element_orig = fm.createFieldConstant([0.5*radiansPerElementOrig])
        radians_per_element_transition = fm.createFieldConstant([radiansPerElementTransition])

        theta_continuous = fm.createFieldIf(fm.createFieldLessThan(theta, half_radians_per_element_orig), fm.createFieldAdd(theta, two_pi), theta)
        theta_offset = fm.createFieldSubtract(theta_continuous, radians_per_element_orig)
        theta_offset_septum_start = fm.createFieldConstant([-0.5*radiansPerElementOrig])
        theta_offset_septum_end = fm.createFieldConstant([(elementsCountAcrossSeptum + 0.5)*radiansPerElementOrig])
        in_septum = fm.createFieldAnd(fm.createFieldGreaterThan(theta_offset, theta_offset_septum_start), fm.createFieldLessThan(theta_offset, theta_offset_septum_end))
        septum_scale = fm.createFieldConstant([radiansPerElementSeptum/radiansPerElementOrig])
        thetaNewSeptumStart = radiansPerElementOrig + (radiansPerElementOrig - radiansPerElementSeptum)*elementsCountAcrossSeptum/2.0
        theta_new_septum_start = fm.createFieldConstant([thetaNewSeptumStart])
        theta_new_septum = fm.createFieldAdd(fm.createFieldMultiply(theta_offset, septum_scale), theta_new_septum_start)
        lvfreewall_scale = fm.createFieldConstant([radiansPerElementLVFreeWall/radiansPerElementOrig])
        theta_offset_lvfreewall_start = fm.createFieldConstant([radiansPerElementOrig*(elementsCountAcrossSeptum + 1.0)])
        theta_new_lvfreewall_start = fm.createFieldConstant([thetaNewSeptumStart + septumArcAngleRadians + radiansPerElementTransition])
        theta_new_lvfreewall = fm.createFieldAdd(
            fm.createFieldMultiply(fm.createFieldSubtract(theta_offset, theta_offset_lvfreewall_start), lvfreewall_scale),
            theta_new_lvfreewall_start)
        theta_new = fm.createFieldIf(in_septum, theta_new_septum, theta_new_lvfreewall)

        sp_new = fm.createFieldConcatenate([r, theta_new, phi])
        sp_new.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_SPHERICAL_POLAR)
        coordinates_new_scale = fm.createFieldCoordinateTransformation(sp_new)
        coordinates_new_scale.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        coordinates_new = fm.createFieldDivide(coordinates_new_scale, xyz_scale)
        nonApexNodesetGroupField = fm.createFieldNodeGroup(nodes)
        nonApexNodesetGroup = nonApexNodesetGroupField.getNodesetGroup()
        nonApexNodesetGroup.addNodesConditional(one)  # add all
        # remove apex nodes:
        for i in range(elementsCountThroughLVWall + 1):
            node = nodes.findNodeByIdentifier(i*now + 1)
            nonApexNodesetGroup.removeNode(node)
        fieldassignment = coordinates.createFieldassignment(coordinates_new)
        fieldassignment.setNodeset(nonApexNodesetGroup)
        fieldassignment.assign()

        baseNodesetGroup = None
        baseNodeGroupField = fm.createFieldNodeGroup(nodes)
        z = fm.createFieldComponent(coordinates, 3)
        is_base = fm.createFieldGreaterThan(z, fm.createFieldConstant([-0.0001]))
        baseNodesetGroup = baseNodeGroupField.getNodesetGroup()
        baseNodesetGroup.addNodesConditional(is_base)
        #print('baseNodesetGroup.getSize()', baseNodesetGroup.getSize())

        if LVWallThicknessRatioBase != 1.0:
            # make LV walls thinner at base
            # get inside node at middle of RV
            now = 1 + elementsCountUp*elementsCountAround
            midRVnid = now - elementsCountAround + 2 + (elementsCountAcrossSeptum // 2)
            #print('midRVnid', midRVnid)
            midRVnode = nodes.findNodeByIdentifier(midRVnid)
            cp_coordinates = fm.createFieldCoordinateTransformation(coordinates)
            cp_coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
            radius = fm.createFieldComponent(cp_coordinates, 1)
            theta = fm.createFieldComponent(cp_coordinates, 2)
            z = fm.createFieldComponent(cp_coordinates, 3)
            cache.setNode(midRVnode)
            #result, cp = cp_coordinates.evaluateReal(cache, 3)
            #coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, innerRadius = radius.evaluateReal(cache, 1)
            #print('innerRadius', innerRadius)
            ir = fm.createFieldConstant([innerRadius*0.9999])
            radius_gt_ir = fm.createFieldGreaterThan(radius, ir)
            radius_minus_ir = fm.createFieldSubtract(radius, ir)
            thickness_scale = fm.createFieldConstant([LVWallThicknessRatioBase])
            delta_radius = fm.createFieldMultiply(radius_minus_ir, thickness_scale)
            new_radius = fm.createFieldAdd(ir, delta_radius)
            new_cp_coordinates = fm.createFieldConcatenate([new_radius, theta, z])
            new_cp_coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
            new_coordinates = fm.createFieldCoordinateTransformation(new_cp_coordinates)
            new_coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)

            fieldassignment = coordinates.createFieldassignment(new_coordinates)
            result = fieldassignment.setNodeset(baseNodesetGroup)
            fieldassignment.assign()

        if LVBaseFlattenRatio != 1.0:
            # flatten LV normal to mid-RV-aorta-mitral axis plus additional flatten angle
            septumCentreRadians = (1.0 + elementsCountAcrossSeptum/2.0)*radiansPerElementOrig
            flattenAxisRadians = septumCentreRadians + LVBaseFlattenAngleRadians - 0.5*math.pi

            half = fm.createFieldConstant([0.5])
            minus_one = fm.createFieldConstant([-1.0])
            minus_two = fm.createFieldConstant([-2.0])

            # flatten ratio = new inner radius / original inner radius
            beta_base = fm.createFieldConstant([1.0 - LVBaseFlattenRatio])
            z = fm.createFieldComponent(coordinates, 3)
            z2 = fm.createFieldMultiply(z, z)
            zero_dist_node = nodes.findNodeByIdentifier(now - nor)
            cache.setNode(zero_dist_node)
            result, z_zero_dist_value = z.evaluateReal(cache, 1)
            #print('zero dist', result, z_zero_dist_value)
            z_zerodist = fm.createFieldConstant([z_zero_dist_value])
            z_zerodist2 = fm.createFieldMultiply(z_zerodist, z_zerodist)
            z_a = fm.createFieldDivide(one, z_zerodist2)
            z_b = fm.createFieldDivide(minus_two, z_zerodist)
            zfact = fm.createFieldAdd(fm.createFieldAdd(fm.createFieldMultiply(z_a, z2), fm.createFieldMultiply(z_b, z)), one)
            beta = fm.createFieldMultiply(zfact, beta_base)  # 1 - squash factor
            alpha = fm.createFieldSubtract(one, beta)  # z-dependent squash factor

            psi = fm.createFieldConstant([flattenAxisRadians])
            ri = fm.createFieldConstant([0.5 - LVWallThickness])
            theta_minus_psi = fm.createFieldSubtract(theta, psi)
            cos_theta_minus_psi = fm.createFieldCos(theta_minus_psi)
            sin_theta_minus_psi = fm.createFieldSin(theta_minus_psi)
            r_minus_ri = fm.createFieldSubtract(r, ri)

            rf = fm.createFieldAdd(fm.createFieldMultiply(alpha, r), fm.createFieldMultiply(beta, r_minus_ri))
            x_new = fm.createFieldMultiply(rf, cos_theta_minus_psi)
            y_new = fm.createFieldMultiply(r, sin_theta_minus_psi)
            r2_new = fm.createFieldAdd(fm.createFieldMultiply(x_new, x_new), fm.createFieldMultiply(y_new, y_new))
            r_new = fm.createFieldSqrt(r2_new)
            theta_minus_psi_raw = fm.createFieldAtan2(y_new, x_new)
            theta_wrap = fm.createFieldAnd(
                fm.createFieldLessThan(theta_minus_psi, minus_one),
                fm.createFieldGreaterThan(theta_minus_psi_raw, one))
            theta_minus_psi_fix = fm.createFieldIf(theta_wrap, fm.createFieldAdd(theta_minus_psi, two_pi), theta_minus_psi)
            # above theta is too great; average with theta_minus_psi_raw
            theta_minus_psi_new = fm.createFieldMultiply(half, fm.createFieldAdd(theta_minus_psi_fix, theta_minus_psi_raw))
            theta_new = fm.createFieldAdd(theta_minus_psi_new, psi)
            sp_new = fm.createFieldConcatenate([r_new, theta_new, phi])
            sp_new.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_SPHERICAL_POLAR)
            coordinates_new_scale = fm.createFieldCoordinateTransformation(sp_new)
            coordinates_new_scale.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
            coordinates_new = fm.createFieldDivide(coordinates_new_scale, xyz_scale)

            fieldassignment = coordinates.createFieldassignment(coordinates_new)
            result = fieldassignment.setNodeset(baseNodesetGroup)
            fieldassignment.assign()

        baseNodesetGroup = None

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        crossAngle = math.pi/8
        sinCrossAngle = math.sin(crossAngle)
        cosCrossAngle = math.cos(crossAngle)
        edgeRotateCrossAngle = math.pi/24
        sinEdgeRotateCrossAngle = math.sin(edgeRotateCrossAngle)
        cosEdgeRotateCrossAngle = math.cos(edgeRotateCrossAngle)

        # create RV nodes and modify adjoining LV nodes
        elementsCountUpRV = elementsCountUp - elementsCountBelowSeptum
        lv_bni_base = 3 + elementsCountThroughLVWall*now + (elementsCountBelowSeptum - 1)*elementsCountAround
        nodeIdentifier = startNodeIdentifier = getMaximumNodeIdentifier(nodes) + 1
        baseExtraRVWidth = LVWallThickness*(1.0 - LVWallThicknessRatioBase)
        rv_nids = []
        for n3 in range(2):  # only 1 element through RV free wall

            onInside = n3 == 0

            for n2 in range(elementsCountUpRV + 1):

                #if n2 > (elementsCountUpRV - elementsCountUpBase):
                #    theta = math.pi*0.5
                #else:
                #    theta = (n2 / (elementsCountUpRV - elementsCountUpBase))*math.pi*0.5
                #theta = (n2 / elementsCountUpRV)*math.pi*0.5
                RVWidth = RVWidthTop  #*math.sin(theta)
                #RVWidth = RVWidthTop

                #if n2 == 0:
                #    edgeOffsetInner = RVWidth*0.25
                #else:
                edgeOffsetInner = RVWidth*0.65
                if n2 == elementsCountUpRV:
                    RVWidth += baseExtraRVWidth
                    edgeOffsetInner = edgeOffsetInner + baseExtraRVWidth*0.5

                onBottomEdge = (n2 == 0)

                for n1 in range(elementsCountAcrossSeptum + 1):

                    onSideEdge1 = (n1 == 0)
                    onSideEdge2 = (n1 == elementsCountAcrossSeptum)
                    onSideEdge = onSideEdge1 or onSideEdge2

                    existingNodeIdentifier = lv_bni_base + n2*nor + n1

                    baseNode = nodes.findNodeByIdentifier(existingNodeIdentifier)
                    cache.setNode(baseNode)
                    result, base_x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                    result, base_dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                    result, base_dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                    result, base_dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                    #print('node ', existingNodeIdentifier, 'dx_ds3', result, base_dx_ds3)
                    mag = math.sqrt(base_dx_ds3[0]*base_dx_ds3[0] + base_dx_ds3[1]*base_dx_ds3[1] + base_dx_ds3[2]*base_dx_ds3[2])
                    unitOutward = [ base_dx_ds3[0]/mag, base_dx_ds3[1]/mag, base_dx_ds3[2]/mag ]
                    baseRadiusAround = unitOutward[0]*base_x[0] + unitOutward[1]*base_x[1] + unitOutward[2]*base_x[2]
                    rotateEdge = False
                    if onInside:
                        if onBottomEdge and onSideEdge:
                            offset = 0.0
                        elif onBottomEdge or onSideEdge:
                            offset = edgeOffsetInner
                            rotateEdge = True
                        else:
                            offset = RVWidth
                    else:
                        if onBottomEdge and onSideEdge:
                            offset = RVFreeWallThickness
                            #rotateEdge = True
                        elif onBottomEdge or onSideEdge:
                            offset = edgeOffsetInner  # + RVFreeWallThickness
                            rotateEdge = True
                        else:
                            offset = RVWidth + RVFreeWallThickness
                    x = [
                        base_x[0] + unitOutward[0]*offset,
                        base_x[1] + unitOutward[1]*offset,
                        base_x[2] + unitOutward[2]*offset,
                    ]
                    scale1 = (baseRadiusAround + offset)/baseRadiusAround
                    dx_ds1 = [
                        base_dx_ds1[0]*scale1,
                        base_dx_ds1[1]*scale1,
                        base_dx_ds1[2]*scale1
                    ]
                    # GRC not sure if appropriate to use scale1 here:
                    dx_ds2 = [
                        base_dx_ds2[0]*scale1,
                        base_dx_ds2[1]*scale1,
                        base_dx_ds2[2]*scale1
                    ]
                    scale3 = RVFreeWallThickness
                    dx_ds3 = [
                        unitOutward[0]*scale3,
                        unitOutward[1]*scale3,
                        unitOutward[2]*scale3
                    ]

                    if rotateEdge:
                        if onSideEdge1:
                            rotatedOutward = [
                                cosEdgeRotateCrossAngle*base_dx_ds3[0] - sinEdgeRotateCrossAngle*base_dx_ds1[0],
                                cosEdgeRotateCrossAngle*base_dx_ds3[1] - sinEdgeRotateCrossAngle*base_dx_ds1[1],
                                cosEdgeRotateCrossAngle*base_dx_ds3[2] - sinEdgeRotateCrossAngle*base_dx_ds1[2]
                            ]
                        elif onSideEdge2:
                            rotatedOutward = [
                                cosEdgeRotateCrossAngle*base_dx_ds3[0] + sinEdgeRotateCrossAngle*base_dx_ds1[0],
                                cosEdgeRotateCrossAngle*base_dx_ds3[1] + sinEdgeRotateCrossAngle*base_dx_ds1[1],
                                cosEdgeRotateCrossAngle*base_dx_ds3[2] + sinEdgeRotateCrossAngle*base_dx_ds1[2]
                            ]
                        else:  # onBottomEdge:
                            rotatedOutward = [
                                cosEdgeRotateCrossAngle*base_dx_ds3[0] - sinEdgeRotateCrossAngle*base_dx_ds2[0],
                                cosEdgeRotateCrossAngle*base_dx_ds3[1] - sinEdgeRotateCrossAngle*base_dx_ds2[1],
                                cosEdgeRotateCrossAngle*base_dx_ds3[2] - sinEdgeRotateCrossAngle*base_dx_ds2[2]
                            ]
                        mag = math.sqrt(rotatedOutward[0]*rotatedOutward[0] + rotatedOutward[1]*rotatedOutward[1] + rotatedOutward[2]*rotatedOutward[2])
                        unitRotatedOutward = [ rotatedOutward[0]/mag, rotatedOutward[1]/mag, rotatedOutward[2]/mag ]
                        scale3r = RVFreeWallThickness  # RVFreeWallThickness + edgeOffsetInner
                        if not onInside:
                            x[0] += unitRotatedOutward[0]*scale3r
                            x[1] += unitRotatedOutward[1]*scale3r
                            x[2] += unitRotatedOutward[2]*scale3r
                        dx_ds3 = [
                            unitRotatedOutward[0]*scale3r,
                            unitRotatedOutward[1]*scale3r,
                            unitRotatedOutward[2]*scale3r
                        ]
                        if onBottomEdge:
                            rscale2 = math.sqrt(dx_ds2[0]*dx_ds2[0] + dx_ds2[1]*dx_ds2[1] + dx_ds2[2]*dx_ds2[2])
                            tang = [
                                unitRotatedOutward[1]*dx_ds1[2] - unitRotatedOutward[2]*dx_ds1[1],
                                unitRotatedOutward[2]*dx_ds1[0] - unitRotatedOutward[0]*dx_ds1[2],
                                unitRotatedOutward[0]*dx_ds1[1] - unitRotatedOutward[1]*dx_ds1[0]
                            ]
                            rscale2 /= math.sqrt(tang[0]*tang[0] + tang[1]*tang[1] + tang[2]*tang[2])
                            dx_ds2 = [
                                tang[0]*rscale2,
                                tang[1]*rscale2,
                                tang[2]*rscale2
                            ]
                        if onSideEdge:
                            rscale1 = math.sqrt(dx_ds1[0]*dx_ds1[0] + dx_ds1[1]*dx_ds1[1] + dx_ds1[2]*dx_ds1[2])
                            tang = [
                                dx_ds2[1]*unitRotatedOutward[2] - dx_ds2[2]*unitRotatedOutward[1],
                                dx_ds2[2]*unitRotatedOutward[0] - dx_ds2[0]*unitRotatedOutward[2],
                                dx_ds2[0]*unitRotatedOutward[1] - dx_ds2[1]*unitRotatedOutward[0]
                            ]
                            rscale1 /= math.sqrt(tang[0]*tang[0] + tang[1]*tang[1] + tang[2]*tang[2])
                            dx_ds1 = [
                                tang[0]*rscale1,
                                tang[1]*rscale1,
                                tang[2]*rscale1
                            ]

                    if onInside and onBottomEdge and onSideEdge:
                        node = baseNode
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        nodeIdentifier += 1
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    rv_nids.append(node.getIdentifier())

        # create RV elements and modify adjoining LV element fields
        elementIdentifier = elementsCountThroughLVWall*elementsCountUp*elementsCountAround + 1
        scalefactors5 = [ -1.0, sinCrossAngle, cosCrossAngle, sinCrossAngle, cosCrossAngle ]
        scalefactors9 = [ -1.0, 0.5, 0.25, 0.125, 0.75, sinCrossAngle, cosCrossAngle, sinCrossAngle, cosCrossAngle ]
        eow = elementsCountUp*elementsCountAround
        RVSeptumElementIdBase = eow*(elementsCountThroughLVWall - 1) + elementsCountBelowSeptum*elementsCountAround + 2

        # Add RV elements

        # RV LV joiner bottom elements: h-shape
        eftRVBottom = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVBottom, [1], [])
        # d/dxi2 is reversed on inside
        for ln in [1, 2]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVBottom, n*8 + 3, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
        # d/dxi3 = -d/dS2 on LV side
        for ln in [1, 2, 5, 6]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVBottom, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
        #print('eftRVBottom', eftRVBottom.validate())

        # RV LV joiner side 1 elements: h-shape
        eftRVSide1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVSide1, [1], [])
        # d/dxi1 is reversed on inside
        for ln in [1, 3]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVSide1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        # d/dxi3 = -d/dS1 on LV side
        for ln in [1, 3, 5, 7]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVSide1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        #print('eftRVSide1', eftRVSide1.validate())

        # RV LV joiner side 2 elements: h-shape
        eftRVSide2 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVSide2, [1], [])
        # d/dxi1 is reversed on inside
        for ln in [2, 4]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVSide2, n*8 + 2, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
        # d/dxi3 = d/dS1 on LV side
        for ln in [2, 4, 6, 8]:
            n = ln - 1
            mapEftFunction1Node1Term(eftRVSide2, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
        #print('eftRVSide2', eftRVSide2.validate())

        rv_nor = (elementsCountAcrossSeptum + 1)
        rv_now = (elementsCountUpRV + 1)*rv_nor

        for n2 in range(-1, elementsCountUpRV):

            for n1 in range(-1, elementsCountAcrossSeptum + 1):

                bni = lv_bni_base + n2*elementsCountAround + n1
                rv_bni = n2*rv_nor + n1
                eft1 = None

                if n2 == -1:
                    if n1 == -1:
                        # pinched to zero thickness on 3 sides
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        remapEftNodeValueLabel(eft1, [ 1, 2, 3, 5, 6, 7 ], Node.VALUE_LABEL_D_DS3, [])
                        ln_map = [ 1, 2, 3, 4, 1, 2, 3, 5 ]
                        remapEftLocalNodes(eft1, 5, ln_map)
                        nodeIdentifiers = [ bni, bni + 1, bni + nor, bni + nor + 1, rv_nids[rv_bni + rv_nor + rv_now + 1] ]
                        #print('eft1 corner a', eft1.validate())
                    elif n1 == elementsCountAcrossSeptum:
                        # pinched to zero thickness on 3 sides
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        remapEftNodeValueLabel(eft1, [ 1, 2, 4, 5, 6, 8 ], Node.VALUE_LABEL_D_DS3, [])
                        ln_map = [ 1, 2, 3, 4, 1, 2, 5, 4 ]
                        remapEftLocalNodes(eft1, 5, ln_map)
                        nodeIdentifiers = [ bni, bni + 1, bni + nor, bni + nor + 1, rv_nids[rv_bni + rv_nor + rv_now] ]
                        #print('eft1 corner b', eft1.validate())
                    else:
                        nodeIdentifiers = [
                            bni + nor, bni + nor + 1, rv_nids[rv_bni + rv_nor], rv_nids[rv_bni + rv_nor + 1],
                            bni, bni + 1, rv_nids[rv_bni + rv_nor + rv_now], rv_nids[rv_bni + rv_nor + 1 + rv_now]
                        ]
                        if n1 == 0:
                            eft1 = tricubichermite.createEftBasic()
                            # GRC to fix: uses repeated nodes instead of reducing to 7
                            #eft1.setNumberOfLocalNodes(7)
                            setEftScaleFactorIds(eft1, [1], [])
                            # d/dxi2 is zero on left side, reversed on inner right side
                            eft1.setFunctionNumberOfTerms(0*8 + 3, 0)
                            eft1.setFunctionNumberOfTerms(2*8 + 3, 0)
                            mapEftFunction1Node1Term(eft1, 1*8 + 3, 2, Node.VALUE_LABEL_D_DS2, 1, [1])
                            # d/dxi3 = -d/dS2 on LV side
                            for ln in [1, 2, 5, 6]:
                                n = ln - 1
                                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
                        elif n1 == (elementsCountAcrossSeptum - 1):
                            eft1 = tricubichermite.createEftBasic()
                            # GRC to fix: uses repeated nodes instead of reducing to 7
                            #eft1.setNumberOfLocalNodes(7)
                            setEftScaleFactorIds(eft1, [1], [])
                            # d/dxi2 is zero on right side, reversed on inner left side
                            eft1.setFunctionNumberOfTerms(1*8 + 3, 0)
                            eft1.setFunctionNumberOfTerms(3*8 + 3, 0)
                            mapEftFunction1Node1Term(eft1, 0*8 + 3, 1, Node.VALUE_LABEL_D_DS2, 1, [1])
                            # d/dxi3 = -d/dS2 on LV side
                            for ln in [1, 2, 5, 6]:
                                n = ln - 1
                                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [1])
                        else:
                            eft1 = eftRVBottom
                else:
                    if n1 == -1:
                        nodeIdentifiers = [
                            bni + 1, rv_nids[rv_bni + 1], bni + nor + 1, rv_nids[rv_bni + rv_nor + 1],
                            bni, rv_nids[rv_bni + 1 + rv_now], bni + nor, rv_nids[rv_bni + rv_nor + 1 + rv_now]
                        ]
                        if n2 == 0:
                            eft1 = tricubichermite.createEftBasic()
                            # GRC to fix: uses repeated nodes instead of reducing to 7
                            #eft1.setNumberOfLocalNodes(7)
                            setEftScaleFactorIds(eft1, [1], [])
                            # d/dxi1 is zero on bottom side, reversed on inner top side
                            eft1.setFunctionNumberOfTerms(0*8 + 2, 0)
                            eft1.setFunctionNumberOfTerms(1*8 + 2, 0)
                            mapEftFunction1Node1Term(eft1, 2*8 + 2, 3, Node.VALUE_LABEL_D_DS1, 1, [1])
                            # d/dxi3 = -d/dS1 on LV side
                            for ln in [1, 3, 5, 7]:
                                n = ln - 1
                                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [1])
                            #print('eft1 next', eft1.validate())
                        else:
                            eft1 = eftRVSide1
                    elif n1 == elementsCountAcrossSeptum:
                        nodeIdentifiers = [
                            rv_nids[rv_bni], bni, rv_nids[rv_bni + rv_nor], bni + nor,
                            rv_nids[rv_bni + rv_now], bni + 1, rv_nids[rv_bni + rv_nor + rv_now], bni + nor + 1
                        ]
                        if n2 == 0:
                            eft1 = tricubichermite.createEftBasic()
                            # GRC to fix: uses repeated nodes instead of reducing to 7
                            #eft1.setNumberOfLocalNodes(7)
                            setEftScaleFactorIds(eft1, [1], [])
                            # d/dxi1 is zero on bottom side, reversed on inner top side
                            eft1.setFunctionNumberOfTerms(0*8 + 2, 0)
                            eft1.setFunctionNumberOfTerms(1*8 + 2, 0)
                            mapEftFunction1Node1Term(eft1, 3*8 + 2, 4, Node.VALUE_LABEL_D_DS1, 1, [1])
                            # d/dxi3 = d/dS1 on LV side
                            for ln in [2, 4, 6, 8]:
                                n = ln - 1
                                mapEftFunction1Node1Term(eft1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [])
                            #print('eft1 next', eft1.validate())
                        else:
                            eft1 = eftRVSide2
                    else:
                        eft1 = eft
                        nodeIdentifiers = [
                            rv_nids[rv_bni], rv_nids[rv_bni + 1], rv_nids[rv_bni + rv_nor], rv_nids[rv_bni + rv_nor + 1],
                            rv_nids[rv_bni + rv_now], rv_nids[rv_bni + 1 + rv_now], rv_nids[rv_bni + rv_nor + rv_now], rv_nids[rv_bni + rv_nor + 1 + rv_now]
                        ]

                useElementtemplate = mesh.createElementtemplate()
                useElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result1 = useElementtemplate.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, useElementtemplate)
                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft1, [-1.0])
                #print('RV element create', elementIdentifier, result1, result2, nodeIdentifiers)
                elementIdentifier += 1

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
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountAcrossSeptum = options['Number of elements across septum']
        elementsCountBelowSeptum = options['Number of elements below septum']

        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']

        startRvElementIdentifier = elementsCountAround*elementsCountUp*elementsCountThroughLVWall + 1
        limitRvElementIdentifier = startRvElementIdentifier + (elementsCountAcrossSeptum + 2)*(elementsCountUp - elementsCountBelowSeptum + 1)

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
            MeshType_3d_heartventricles1.generateBaseMesh(region, options)
            return
        baseRegion = region.createRegion()
        MeshType_3d_heartventricles1.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region)
        MeshType_3d_heartventricles1.refineMesh(meshrefinement, options)
