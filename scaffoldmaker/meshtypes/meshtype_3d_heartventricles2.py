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
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


def getRvOuterPoints(rvArcRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAround, dEndMag, z):
    '''
    Get array of points and derivatives around RV, centred at +x axis.
    LV is assumed to be circular, centred at x, y = (0,0)
    :param rvArcRadians: Angle in radians around outside of RV to tangent nodes where it meets LV.
    :param lvRadius: Radius of the LV.
    :param rvAddWidthRadius: Radius to add to LV to get RV width at centre.
    :param rvAddCrossRadius: Additional radius to add only laterally around septum.
    :param elementsCountAround: Number of elements around RV.
    :param dEndMag: Magnitude of the external derivative at each end.
    :param z: Z coordinate to give to all values of x[]. dx_ds1[2] is all zero.
    :return: Arrays x[], dx_ds1[], starting at first point after end at lowest angle around.
    '''
    print('\ngetRvOuterPoints', rvArcRadians, lvRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAround, dEndMag, z)
    startRadians = -0.5*rvArcRadians
    ellipseStartRadians = startRadians + 0.5*rvArcRadians/elementsCountAround
    b = lvRadius + rvAddCrossRadius  # cross radius
    theta = math.asin(-lvRadius*math.sin(ellipseStartRadians)/b) - math.pi
    print('phi',ellipseStartRadians)
    print('theta',theta)
    a = (lvRadius*(math.cos(ellipseStartRadians) - 1.0) - rvAddWidthRadius)/(math.cos(theta) - 1.0)  # width radius
    cx = lvRadius + rvAddWidthRadius - a
    print(ellipseStartRadians,'a',a,'b',b,'cx',cx)
    # get cubic curve joining start point with half ellipse, first with unit derivatives then compute arc length
    x1 = ( lvRadius*math.cos(startRadians), lvRadius*math.sin(startRadians) )
    d1 = ( -math.sin(startRadians), math.cos(startRadians) )
    x2 = ( cx, -b )
    d2 = ( 1.0, 0.0 )
    cubicArcLength = computeCubicHermiteArcLength(x1, d1, x2, d2)
    d1 = ( d1[0]*cubicArcLength, d1[1]*cubicArcLength )
    d2 = ( d2[0]*cubicArcLength, d2[1]*cubicArcLength )
    halfEllipsePerimeter = 0.5*getApproximateEllipsePerimeter(a, b)
    totalLength = 2.0*cubicArcLength + halfEllipsePerimeter
    print('length cubic', cubicArcLength, 'halfEllipsePerimeter',halfEllipsePerimeter,'total',totalLength)
    midLength = (totalLength - dEndMag)/(elementsCountAround - 1)
    endLength = 0.5*(dEndMag + midLength)
    print('midLength',midLength,'endLength',endLength)
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
            print('interpolate',x1,d1,x2,d2, xi)
            v = interpolateCubicHermite(x1, d1, x2, d2, xi)
            t = interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        x.insert(0, [ v[0], v[1], z ] )
        scale = midLength/math.sqrt(t[0]*t[0] + t[1]*t[1])
        dx_ds1.insert(0, [ t[0]*scale, t[1]*scale, 0.0 ])
        length -= midLength
    length += midLength - endLength
    print('length',length,' remainder',totalLength*0.5+length)
    # mirror in y (about x axis) to get other side
    for n in range((elementsCountAround - 1)//2 - 1, -1, -1):
        x.append([ x[n][0], -x[n][1], z ])
        dx_ds1.append([ -dx_ds1[n][0], dx_ds1[n][1], 0.0 ])
    return x, dx_ds1

def getRVOuterSize(xiUp, rvWidthRadius, RVExtraCrossRadiusBase):
    xiUpFast = 1.0 - (1.0 - xiUp)*(1.0 - xiUp)
    xiUpSlow = xiUp*xiUp
    rvAddWidthRadius = interpolateCubicHermite([0.0], [0.0], [rvWidthRadius], [0.0], xiUpFast)[0]
    rvAddCrossRadius = interpolateCubicHermite([0.0], [0.0], [RVExtraCrossRadiusBase], [0.0], xiUpSlow)[0]
    return rvAddWidthRadius, rvAddCrossRadius

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
            'Number of elements around' : 12,
            'Number of elements up' : 4,
            'Number of elements through LV wall' : 1,
            'Number of elements around septum' : 5,
            'Number of elements below septum' : 2,
            'Length ratio' : 1.8,
            'Element length ratio equator/apex' : 1.0,
            'LV wall thickness' : 0.15,
            'LV wall thickness ratio apex' : 0.5,
            'LV wall thickness ratio base': 1.0,
            'LV base flatten ratio': 1.0,
            'LV base flatten angle degrees': 0.0,
            'RV free wall thickness' : 0.05,
            'RV width' : 0.3,
            'RV extra cross radius base' : 0.1,
            'Septum arc angle degrees' : 120.0,
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
            'Number of elements around septum',
            'Number of elements below septum',
            'Length ratio',
            'Element length ratio equator/apex',
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'LV wall thickness ratio base',
            'LV base flatten ratio',
            'LV base flatten angle degrees',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
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
        if options['Number of elements around septum'] < 2:
            options['Number of elements around septum'] = 2
        elif options['Number of elements around septum'] > (options['Number of elements around'] - 2):
            options['Number of elements around septum'] = options['Number of elements around'] - 2
        if options['Number of elements below septum'] < 2:
            options['Number of elements below septum'] = 2
        elif options['Number of elements below septum'] > (options['Number of elements up'] - 1):
            options['Number of elements below septum'] = options['Number of elements up'] - 1
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6
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
        if options['RV extra cross radius base'] < 0.0:
            options['RV extra cross radius base'] = 0.0
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
        elementsCountAroundSeptum = options['Number of elements around septum']
        elementsCountBelowSeptum = options['Number of elements below septum']
        LVWallThickness = options['LV wall thickness']
        LVWallThicknessRatioApex = options['LV wall thickness ratio apex']
        LVWallThicknessRatioBase = options['LV wall thickness ratio base']
        LVBaseFlattenRatio = options['LV base flatten ratio']
        LVBaseFlattenAngleRadians = options['LV base flatten angle degrees']*math.pi/180.0

        lengthRatio = options['Length ratio']
        RVFreeWallThickness = options['RV free wall thickness']
        RVWidthTop = options['RV width']
        RVExtraCrossRadiusBase = options['RV extra cross radius base']

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

        # Rotate model so RV & Septum are centred and symmetric about +x axis
        radiansPerElementOrig = 2.0*math.pi/elementsCountAround
        septumCentreRadians = (1.0 + elementsCountAroundSeptum/2.0)*radiansPerElementOrig
        cosSeptumCentre = math.cos(-septumCentreRadians)
        sinSeptumCentre = math.sin(-septumCentreRadians)
        rotation_matrix = fm.createFieldConstant([ cosSeptumCentre, -sinSeptumCentre, 0.0, sinSeptumCentre, cosSeptumCentre, 0.0, 0.0, 0.0, 1.0 ])
        coordinates_new = fm.createFieldMatrixMultiply(3, rotation_matrix, coordinates)
        fieldassignment = coordinates.createFieldassignment(coordinates_new)
        fieldassignment.assign()
        del coordinates_new
        del rotation_matrix
        del septumCentreRadians  # since now zero

        # Resize elements around LV to get desired septum arc angle
        radiansPerElementSeptum = septumArcAngleRadians/elementsCountAroundSeptum
        # want LV-RV 'transition' elements to be mean size of Septum and LV FreeWall elements
        radiansRemaining = 2.0*math.pi - septumArcAngleRadians
        elementsCountAroundLVFreeWall = elementsCountAround - elementsCountAroundSeptum - 2
        radiansPerElementLVFreeWall = (radiansRemaining - radiansPerElementSeptum) / (elementsCountAroundLVFreeWall + 1)
        radiansPerElementTransition = 0.5*(radiansPerElementSeptum + radiansPerElementLVFreeWall)
        #print('Element size ratio LVFreeWall / Septum', radiansPerElementLVFreeWall/radiansPerElementSeptum)

        cp = fm.createFieldCoordinateTransformation(coordinates)
        z = fm.createFieldComponent(coordinates, 3)
        cp.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
        r = fm.createFieldComponent(cp, 1)
        theta = fm.createFieldComponent(cp, 2)
        zero = fm.createFieldConstant([0.0])
        one = fm.createFieldConstant([1.0])  # also used as 'true'
        two_pi = fm.createFieldConstant([2.0*math.pi])
        pi = fm.createFieldConstant([math.pi])
        thetaSeptumLimit = 0.5*radiansPerElementOrig*(elementsCountAroundSeptum + 1)
        theta_septum_start = fm.createFieldConstant([-thetaSeptumLimit])
        theta_septum_end = fm.createFieldConstant([+thetaSeptumLimit])
        septum_scale = fm.createFieldConstant([radiansPerElementSeptum/radiansPerElementOrig])
        freewall_scale = fm.createFieldConstant([radiansPerElementLVFreeWall/radiansPerElementOrig])

        in_septum = fm.createFieldAnd(fm.createFieldGreaterThan(theta, theta_septum_start), fm.createFieldLessThan(theta, theta_septum_end))
        theta_septum = fm.createFieldMultiply(theta, septum_scale)
        is_negative = fm.createFieldLessThan(theta, zero)
        theta_freewall_negative = fm.createFieldSubtract(fm.createFieldMultiply(fm.createFieldAdd(theta, pi), freewall_scale), pi)
        theta_freewall_positive = fm.createFieldAdd(fm.createFieldMultiply(fm.createFieldSubtract(theta, pi), freewall_scale), pi)
        theta_freewall = fm.createFieldIf(is_negative, theta_freewall_negative, theta_freewall_positive)
        theta_new = fm.createFieldIf(in_septum, theta_septum, theta_freewall)

        cp_new = fm.createFieldConcatenate([r, theta_new, z])
        cp_new.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
        coordinates_new = fm.createFieldCoordinateTransformation(cp_new)
        coordinates_new.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)

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
            midRVnid = now - elementsCountAround + 2 + (elementsCountAroundSeptum // 2)
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
            flattenAxisRadians = LVBaseFlattenAngleRadians - 0.5*math.pi

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

        #################
        # Create nodes
        #################

        elementsCountUpRV = elementsCountUp - elementsCountBelowSeptum + 1
        nodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

        norl = elementsCountAround
        nowl = 1 + elementsCountUp*norl
        nidl = elementsCountThroughLVWall*nowl + 2 + (elementsCountBelowSeptum - 2)*norl
        node = nodes.findNodeByIdentifier(nidl)
        cache.setNode(node)
        
        result, rvBottomX = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, rvBottomD2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        rvOuterHeight = rvBottomX[2]
        rvOuterDz = rvBottomD2[2]
        print('node',nidl,'rvOuterHeight',rvOuterHeight,'rvOuterDz',rvOuterDz)

        rvArcRadians = septumArcAngleRadians + 2*radiansPerElementTransition

        elementsCountAroundRV = elementsCountAroundSeptum + 2
        rvOuterWidthBase = RVWidthTop + RVFreeWallThickness

        rxOuter = []
        rd1Outer = []
        rd2Outer = []
        rd3Outer = []
        for n2 in range(elementsCountUpRV):
            nida = nidl + norl*(n2 + 1)
            node = nodes.findNodeByIdentifier(nida)
            cache.setNode(node)
            result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)

            # get node at level, hence z
            z = ax[2]
            lvOuterRadius = math.sqrt(ax[0]*ax[0] + ax[1]*ax[1])
            dEndMag = math.sqrt(sum(d*d for d in ad1))

            xiUp = 1.0 - z/rvOuterHeight

            rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUp, rvOuterWidthBase, RVExtraCrossRadiusBase)
            #print('z',z,'xiUp', xiUp, 'xiUpFast', xiUpFast, 'xiUpSlow', xiUpSlow)
            #print('rvAddWidthRadius =', rvAddWidthRadius, ', rvAddCrossRadius=', rvAddCrossRadius)
            rx, rd1 = getRvOuterPoints(rvArcRadians, lvOuterRadius, rvAddWidthRadius, rvAddCrossRadius, elementsCountAroundRV, dEndMag, z)

            if n2 < (elementsCountUpRV - 1):
                # get points at node position plus increment of xi2, to calculate derivative 2
                dxi = 1.0E-3
                xPlus = [ ax[i] + ad2[i]*dxi for i in range(3) ]
                lvOuterRadiusPlus = math.sqrt(xPlus[0]*xPlus[0] + xPlus[1]*xPlus[1])
                dEndMagPlus = dEndMag*lvOuterRadiusPlus/lvOuterRadius
                xiUpPlus = 1.0 - xPlus[2]/rvOuterHeight
                rvAddWidthRadiusPlus, rvAddCrossRadiusPlus = getRVOuterSize(xiUpPlus, rvOuterWidthBase, RVExtraCrossRadiusBase)
                px, pd1 = getRvOuterPoints(rvArcRadians, lvOuterRadiusPlus, rvAddWidthRadiusPlus, rvAddCrossRadiusPlus, elementsCountAroundRV, dEndMagPlus, xPlus[2])

                rd2 = [ [(px[n][c] - rx[n][c])/dxi for c in range(3)] for n in range(len(rx)) ]
            else:
                rd2 = [ad2]*len(rx)

            rd3 = []
            for n in range(len(rx)):
                normal = vector.crossproduct3(rd1[n], rd2[n])
                mag = RVFreeWallThickness/vector.magnitude(normal)
                rd3.append([ c*mag for c in normal ])

            rxOuter.append(rx)
            rd1Outer.append(rd1)
            rd2Outer.append(rd2)
            rd3Outer.append(rd3)

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
                factor = 1.0 - curvature*RVFreeWallThickness
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
                factor = 1.0 - curvature*RVFreeWallThickness
                id2.append([ factor*c for c in rd2Outer[n2][n] ])

            rxInner.append(ix)
            rd1Inner.append(id1)
            rd2Inner.append(id2)

        zero = [ 0.0, 0.0, 0.0 ]
        rv_nids = [[],[]]
        for n3 in range(2):
            for n2 in range(-2, elementsCountUpRV):
                nids = []
                rv_nids[n3].append(nids)
                n2r = max(0, n2)
                if n3 == 0:
                    rx = rxInner[n2r]
                    rd1 = rd1Inner[n2r]
                    rd2 = rd2Inner[n2r]
                    rd3 = rd3Inner[n2r]
                else:
                    rx = rxOuter[n2r]
                    rd1 = rd1Outer[n2r]
                    rd2 = rd2Outer[n2r]
                    rd3 = rd3Outer[n2r]
                for n1 in range(-2, elementsCountAroundRV + 1):
                    n1r = min(max(0, n1), elementsCountAroundRV - 2)
                    nid = -1
                    if n2 < 0:
                        nid = nidl + norl + 1
                        if n1 == -2:
                            pass
                        elif n1 == elementsCountAroundRV:
                            nid += elementsCountAroundRV - 2
                        elif n2 == -2:
                            nid += n1
                        elif n3 == 0:  # and (n2 == -1)
                            if n1 != n1r:
                                # create but calculate later
                                x = dx_dxi1 = dx_dxi2 = dx_dxi3 = zero
                                nid = -1
                            else:
                                # interpolate node inside RV bottom edge
                                nid += n1
                                #print('bottom nid', nid)
                                node = nodes.findNodeByIdentifier(nid)
                                cache.setNode(node)
                                result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                                result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                                ad2 = [ -d for d in ad2 ]
                                nid -= norl
                                node = nodes.findNodeByIdentifier(nid)
                                cache.setNode(node)
                                result, cx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                #print('  bottom nid', nid)
                                bx = rx[n1r]
                                bd1 = rd1[n1r]
                                bd2 = rd2[n1r]
                                ad2 = [ 1.5*d for d in ad2 ]
                                bd2 = [ 1.5*d for d in bd2 ]
                                x = list(interpolateCubicHermite(ax, ad2, bx, bd2, 0.5))
                                dx_dxi2 = interpolateCubicHermiteDerivative(ax, ad2, bx, bd2, 0.5)
                                dx_dxi2 = [ (1.0/3.0)*d for d in dx_dxi2 ]
                                dx_dxi1 = [ 0.5*(ad1[c] + bd1[c]) for c in range(3) ]
                                dx_dxi3 = [ (cx[c] - x[c]) for c in range(3) ]
                                nid = -1
                        else: # n3 == 1
                            nid += n1 - norl
                    else:  # 0 <= n2
                        if n1 != n1r:
                            nid = nidl + norl*(n2 + 1) + (0 if (n1 < 0) else elementsCountAroundRV)
                            if n1 == -2:
                                nid += 1
                            elif n1 == elementsCountAroundRV:
                                nid -= 1
                            elif n3 == 0:
                                # interpolate node inside RV side edge
                                #print('nid', nid)
                                node = nodes.findNodeByIdentifier(nid)
                                cache.setNode(node)
                                result, cx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                if n1 < 0:
                                    nid += 1
                                else:
                                    nid -= 1
                                #print('  nid', nid)
                                node = nodes.findNodeByIdentifier(nid)
                                cache.setNode(node)
                                result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                                result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                                ad1 = [ -d for d in ad1 ]
                                bx = rx[n1r]
                                bd1 = rd1[n1r]
                                bd2 = rd2[n1r]
                                ad1 = [ 1.5*d for d in ad1 ]
                                bd1 = [ 1.5*d for d in bd1 ]
                                if n1 < 0:
                                    x = list(interpolateCubicHermite(ax, ad1, bx, bd1, 0.5))
                                    dx_dxi1 = interpolateCubicHermiteDerivative(ax, ad1, bx, bd1, 0.5)
                                else:
                                    x = list(interpolateCubicHermite(bx, bd1, ax, ad1, 0.5))
                                    dx_dxi1 = interpolateCubicHermiteDerivative(bx, bd1, ax, ad1, 0.5)
                                dx_dxi1 = [ (1.0/3.0)*d for d in dx_dxi1 ]
                                dx_dxi2 = [ 0.5*(ad2[c] + bd2[c]) for c in range(3) ]
                                dx_dxi3 = [ (cx[c] - x[c]) for c in range(3) ]
                                nid = -1
                        else:
                            x = rx[n1r]
                            dx_dxi1 = rd1[n1r]
                            dx_dxi2 = rd2[n1r]
                            dx_dxi3 = rd3[n1r]

                    if nid < 0:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_dxi1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_dxi2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_dxi3)
                        nid = nodeIdentifier
                        nodeIdentifier += 1
                    nids.append(nid)

        # set parameters for inside corners of RV
        fix_nids = [
            [ rv_nids[0][1][1], rv_nids[0][1][2], rv_nids[0][2][1], rv_nids[1][1][1] ],
            [ rv_nids[0][1][-2], rv_nids[0][1][-3], rv_nids[0][2][-2], rv_nids[1][1][-2] ]
        ]
        #print('fix_nids', fix_nids)
        for i in range(len(fix_nids)):
            nid = fix_nids[i][0]
            nida = fix_nids[i][1]
            nidb = fix_nids[i][2]
            nidc = fix_nids[i][3]
            node = nodes.findNodeByIdentifier(nida)
            cache.setNode(node)
            result, ax = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, ad1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, ad2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            node = nodes.findNodeByIdentifier(nidb)
            cache.setNode(node)
            result, bx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, bd1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, bd2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            node = nodes.findNodeByIdentifier(nidc)
            cache.setNode(node)
            result, cx = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            if i == 0:
                ad1 = [ -d for d in ad1 ]
            ad1 = [ 1.5*d for d in ad1 ]
            bd2 = [ 1.5*d for d in bd2 ]
            x = list(interpolateCubicHermite(ax, ad1, bx, bd2, 0.5))
            dx_dxi2 = interpolateCubicHermiteDerivative(ax, ad1, bx, bd2, 0.5)
            dx_dxi2 = [ (1.0/3.0)*d for d in dx_dxi2 ]
            dx_dxi1 = [ 0.5*(ad2[c] + bd1[c]) for c in range(3) ]
            dx_dxi3 = [ (cx[c] - x[c]) for c in range(3) ]
            node = nodes.findNodeByIdentifier(nid)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_dxi1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_dxi2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_dxi3)

        #################
        # Create elements
        #################

        # redefine elements on boundary of RV to be slightly concave

        eftSplitBottomInner = tricubichermite.createEftSplitXi2RightStraight()
        eftSplitSideInner = tricubichermite.createEftSplitXi1RightStraight()
        elementtemplateSplitSideInner = mesh.createElementtemplate()
        elementtemplateSplitSideInner.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateSplitSideInner.defineField(coordinates, -1, eftSplitSideInner)

        eor = elementsCountAround
        eow = eor*elementsCountUp
        eidl = eow*(elementsCountThroughLVWall - 1) + eor*(elementsCountBelowSeptum - 1) + 1
        print('eidl', eidl)

        for e1 in range(0, elementsCountAroundRV):
            eid = eidl + e1
            element = mesh.findElementByIdentifier(eid)
            eft1 = element.getElementfieldtemplate(coordinates, -1)
            nids  = getElementNodeIdentifiers(element, eft1)
            if e1 == 0:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
            elif e1 == (elementsCountAroundRV - 1):
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
            else:
                eft1 = eftSplitBottomInner
            elementtemplate1 = mesh.createElementtemplate()
            elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element.merge(elementtemplate1)
            element.setNodesByIdentifier(eft1, nids)
            element.setScaleFactors(eft1, [ -1.0 ])

        for e2 in range(1, elementsCountUpRV):
            for e1 in [ 0, elementsCountAroundRV - 1]:
                eid = eidl + e2*eor + e1
                element = mesh.findElementByIdentifier(eid)
                eft1 = element.getElementfieldtemplate(coordinates, -1)
                nids  = getElementNodeIdentifiers(element, eft1)
                element.merge(elementtemplateSplitSideInner)
                element.setNodesByIdentifier(eftSplitSideInner, nids)
                element.setScaleFactors(eftSplitSideInner, [ -1.0 ])

        # create RV elements

        elementIdentifier = eow*elementsCountThroughLVWall + 1

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        eftOuterExit1 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterExit1, [1], [])
        remapEftNodeValueLabel(eftOuterExit1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
        remapEftNodeValueLabel(eftOuterExit1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ])
        remapEftNodeValueLabel(eftOuterExit1, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eftOuterExit1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftOuterExit1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
        ln_map = [ 1, 2, 3, 4, 1, 5, 3, 6 ]
        remapEftLocalNodes(eftOuterExit1, 6, ln_map)

        eftOuterExit2 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterExit2, [1], [])
        remapEftNodeValueLabel(eftOuterExit2, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]) ])

        eftOuterEntry1 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterEntry1, [1], [])
        remapEftNodeValueLabel(eftOuterEntry1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
        remapEftNodeValueLabel(eftOuterEntry1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ])
        remapEftNodeValueLabel(eftOuterEntry1, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eftOuterEntry1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []) ])
        remapEftNodeValueLabel(eftOuterEntry1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
        ln_map = [ 1, 2, 3, 4, 5, 2, 6, 4 ]
        remapEftLocalNodes(eftOuterEntry1, 6, ln_map)

        eftOuterEntry2 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterEntry2, [1], [])
        remapEftNodeValueLabel(eftOuterEntry2, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []) ])

        eftOuterApexExit1 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterApexExit1, [1], [])
        remapEftNodeValueLabel(eftOuterApexExit1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
        remapEftNodeValueLabel(eftOuterApexExit1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ])
        remapEftNodeValueLabel(eftOuterApexExit1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eftOuterApexExit1, [ 7, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftOuterApexExit1, [ 7, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
        ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
        remapEftLocalNodes(eftOuterApexExit1, 6, ln_map)

        eftOuterApexExit2 = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftOuterApexExit2, [1], [])
        remapEftNodeValueLabel(eftOuterApexExit2, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]) ])

        for e2 in range(elementsCountUpRV + 1):
            for e1 in range(elementsCountAroundRV + 2):
                eft1 = eft

                nids = [ rv_nids[0][e2][e1], rv_nids[0][e2][e1 + 1], rv_nids[0][e2 + 1][e1], rv_nids[0][e2 + 1][e1 + 1], \
                         rv_nids[1][e2][e1], rv_nids[1][e2][e1 + 1], rv_nids[1][e2 + 1][e1], rv_nids[1][e2 + 1][e1 + 1] ]

                if e2 == 0:
                    if (e1 == 0) or (e1 >= elementsCountAroundRV):
                        continue
                    if e1 == 1:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
                        remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, [ ])
                        remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
                        remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
                        remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
                        remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                        ln_map = [ 1, 2, 3, 4, 5, 6, 3, 7 ]
                        remapEftLocalNodes(eft1, 7, ln_map)
                        nids = [ rv_nids[0][1][2], rv_nids[0][1][1], rv_nids[0][0][0], rv_nids[0][2][1], \
                                 rv_nids[1][1][2], rv_nids[1][1][1], rv_nids[1][2][1] ]
                    else:
                        eft1 = eftOuterApexExit1
                        nids.pop(4)
                        nids.pop(4)
                elif e2 == 1:
                    if (e1 < 2) or (e1 >= elementsCountAroundRV):
                        continue
                    eft1 = eftOuterApexExit2
                else:
                    if e1 == 0:
                        eft1 = eftOuterExit1
                        nids.pop(4)
                        nids.pop(5)
                    if e1 == 1:
                        eft1 = eftOuterExit2
                    elif e1 == elementsCountAroundRV:
                        eft1 = eftOuterEntry2
                    elif e1 == (elementsCountAroundRV + 1):
                        eft1 = eftOuterEntry1
                        nids.pop(5)
                        nids.pop(6)

                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                result = elementtemplate1.defineField(coordinates, -1, eft1)

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = 7
                print('create element rv', elementIdentifier, result, result2, result3, nids)
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
        elementsCountAroundSeptum = options['Number of elements around septum']
        elementsCountBelowSeptum = options['Number of elements below septum']

        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through RV wall']

        startRvElementIdentifier = elementsCountAround*elementsCountUp*elementsCountThroughLVWall + 1
        elementsCountUpRV = elementsCountUp - elementsCountBelowSeptum + 1
        limitRvElementIdentifier = startRvElementIdentifier + (elementsCountAroundSeptum + 4)*elementsCountUpRV + elementsCountAroundSeptum + 2

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
