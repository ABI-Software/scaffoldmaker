"""
Generates a 3-D lens mesh with variable numbers of elements
around, up the central axis, and radially.
"""

from __future__ import division
import math
from scaffoldmaker.meshtypes.meshtype_3d_solidsphere1 import MeshType_3d_solidsphere1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_lens1:
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Lens 1'

    @staticmethod
    def getDefaultOptions():
        options = MeshType_3d_solidsphere1.getDefaultOptions()
        options['Number of elements around'] = 8
        optionsLens = {
            'Axial thickness' : 4.0,
            'Anterior radius of curvature' : 10.0,
            'Posterior radius of curvature' : 6.0,
            'Spherical radius fraction' : 0.7
        }
        options.update(optionsLens)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_solidsphere1.getOrderedOptionNames()
        optionNames.remove('Diameter')
        for optionName in [
            'Spherical radius fraction',
            'Posterior radius of curvature',
            'Anterior radius of curvature',
            'Axial thickness']:
            optionNames.insert(3, optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_solidsphere1.checkOptions(options)
        for key in [
            'Spherical radius fraction',
            'Posterior radius of curvature',
            'Anterior radius of curvature',
            'Axial thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Spherical radius fraction'] > 1.0:
            options['Spherical radius fraction'] = 1.0

    @classmethod
    def sphereToLens(cls, region, options):
        """
        Map coordinates of the sphere within a width limit olim to arcs with radii 
        rAnt and rPos, respectively. Elements maintain constant size radially and match
        up at dlim. Outside of rlim, affine transformation applied to transform sphere 
        coordinates to be tangential to morphed spherical surfaces
        return: lensRC
        """
        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        sphereCoordinates = getOrCreateCoordinateField(fm)

        radiusSphere = options['Diameter']*0.5
        radiusAnt = options['Anterior radius of curvature']
        radiusPos = options['Posterior radius of curvature']
        lensThickness = options['Axial thickness']
        sphericalRadiusFraction = options['Spherical radius fraction']
        
        # Estimate dLim 
        if lensThickness*0.5 == radiusAnt and lensThickness*0.5 == radiusPos:
            fm.endChange()
            return sphereCoordinates
        A = radiusAnt + radiusPos - lensThickness
        B = radiusAnt*radiusAnt - (radiusPos*radiusPos - radiusAnt*radiusAnt - A*A)**2/(4*A*A)
        if B < 0 or A == 0:
            dLim = min(radiusAnt, radiusPos)
            # print('math domain error - use min radius of curvature')
        elif  math.sqrt(B) < 0.2*min(radiusAnt, radiusPos):
            dLim = min(radiusAnt, radiusPos)
            # print('small intersection - use min radius')
        else:
            dLim = dLim = math.sqrt(B)
        dLimitConstant = dLim*options['Spherical radius fraction']

        # Define constants
        half = fm.createFieldConstant([0.5])
        one = fm.createFieldConstant([1.0])
        rSphere = fm.createFieldConstant([radiusSphere])
        rAnt = fm.createFieldConstant([radiusAnt])
        rPos = fm.createFieldConstant([radiusPos])
        halfLensThickness = fm.createFieldConstant([lensThickness*0.5])

        rLimit = fm.createFieldConstant([0.4])
        psiLimit = fm.createFieldAsin(fm.createFieldDivide(rLimit, rSphere))
        dLimit = fm.createFieldConstant([dLimitConstant])

        # Convert to cylindrical coordinates
        sphereCP = fm.createFieldCoordinateTransformation(sphereCoordinates)
        sphereCP.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
        r = fm.createFieldComponent(sphereCP, 1)
        theta = fm.createFieldComponent(sphereCP, 2)
        z = fm.createFieldComponent(sphereCP, 3)

        # Get anterior weight wa ranging from 0.0 on posterior to 1.0 on anterior, posterior weight wp is reverse
        a = fm.createFieldSqrt(fm.createFieldSubtract(fm.createFieldMultiply(rSphere,rSphere), fm.createFieldMultiply(r,r)))
        psi = fm.createFieldAcos(fm.createFieldDivide(a, rSphere))
        xi = fm.createFieldDivide(z, a) # xi ranges from -1 to 1
        wa = fm.createFieldAdd(fm.createFieldMultiply(half, xi), half)
        wp = fm.createFieldSubtract(one, wa)

        # Spherical anterior mapping
        phiLimitAnt = fm.createFieldAsin(fm.createFieldDivide(dLimit, rAnt))
        dphiAnt_dpsi = fm.createFieldDivide(phiLimitAnt, psiLimit)
        phiAnt = fm.createFieldMultiply(dphiAnt_dpsi, psi)
        rSphereMapAnt = fm.createFieldMultiply(rAnt, fm.createFieldSin(phiAnt))
        zSphereMapAnt = fm.createFieldMultiply(rAnt, fm.createFieldCos(phiAnt))
        displAnt = fm.createFieldSubtract(halfLensThickness, rAnt)
        zSphereMapDisplAnt = fm.createFieldAdd(zSphereMapAnt, displAnt)

        # Spherical posterior mapping
        phiLimitPos = fm.createFieldAsin(fm.createFieldDivide(dLimit, rPos))
        dphiPos_dpsi = fm.createFieldDivide(phiLimitPos, psiLimit)
        phiPos = fm.createFieldMultiply(dphiPos_dpsi, psi)
        rSphereMapPos = fm.createFieldMultiply(rPos, fm.createFieldSin(phiPos))
        zSphereMapPos = fm.createFieldMultiply(rPos, fm.createFieldCos(phiPos))
        displPos = fm.createFieldSubtract(rPos, halfLensThickness)
        zSphereMapDisplPos = fm.createFieldSubtract(displPos, zSphereMapPos)

        # Interpolate between surfaces using weights wa and wp
        rSphereMapLens = fm.createFieldAdd(fm.createFieldMultiply(wa, rSphereMapAnt), fm.createFieldMultiply(wp, rSphereMapPos))
        zSphereMapLens = fm.createFieldAdd(fm.createFieldMultiply(wa, zSphereMapDisplAnt), fm.createFieldMultiply(wp, zSphereMapDisplPos))
        sphereMapLensCP = fm.createFieldConcatenate([rSphereMapLens, theta, zSphereMapLens])
        sphereMapLensCP.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
        sphereMapLensRC = fm.createFieldCoordinateTransformation(sphereMapLensCP)
        sphereMapLensRC.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)

        # Affine transformation beyond rLimit / dLimit
        # Form triangle past rLimit tangent to circle |>, rAffLimit out, zAffLimit half up
        cos_psiLimit = fm.createFieldCos(psiLimit)
        zAffLimit = fm.createFieldMultiply(rSphere, cos_psiLimit)
        tanHalfPiMinusPsi = fm.createFieldDivide(zAffLimit, rLimit)
        rAffLimit = fm.createFieldMultiply(zAffLimit, tanHalfPiMinusPsi)

        # Get unit coordinates across source affine triangle
        r_minus_rLimit = fm.createFieldSubtract(r, rLimit)
        s0Aff = fm.createFieldDivide(r_minus_rLimit, cos_psiLimit)
        psiAff = fm.createFieldDivide(s0Aff, rSphere)
        wrAff = fm.createFieldDivide(r_minus_rLimit, rAffLimit)
        one_minus_wrAff = fm.createFieldSubtract(one, wrAff)
        zAff = fm.createFieldMultiply(zAffLimit, one_minus_wrAff)
        wzAff = fm.createFieldDivide(z, zAff)
        waAff = fm.createFieldAdd(fm.createFieldMultiply(half, wzAff), half)
        wpAff = fm.createFieldSubtract(one, waAff)

        # Get anterior function for affine transformed r, z
        cos_phiLimitAnt = fm.createFieldCos(phiLimitAnt)
        zAffLimitAnt = fm.createFieldMultiply(rAnt, cos_phiLimitAnt)
        zAffDisplAnt = fm.createFieldAdd(displAnt, zAffLimitAnt)
        phiAffAnt = fm.createFieldMultiply(dphiAnt_dpsi, psiAff)
        sAffAnt = fm.createFieldMultiply(rAnt, phiAffAnt)
        oAffAnt = fm.createFieldMultiply(sAffAnt, cos_phiLimitAnt)
        aAffAnt = fm.createFieldMultiply(sAffAnt, fm.createFieldSin(phiLimitAnt))
        rAffAnt = fm.createFieldAdd(dLimit, oAffAnt)
        zAffAnt = fm.createFieldSubtract(zAffDisplAnt, aAffAnt)

        # Get posterior function for affine transformed r,z
        cos_phiLimitPos = fm.createFieldCos(phiLimitPos)
        zAffLimitPos = fm.createFieldMultiply(rPos, cos_phiLimitPos)
        zAffDisplPos = fm.createFieldSubtract(displPos, zAffLimitPos)
        phiAffPos = fm.createFieldMultiply(dphiPos_dpsi, psiAff)
        sAffPos = fm.createFieldMultiply(rPos, phiAffPos)
        oAffPos = fm.createFieldMultiply(sAffPos, cos_phiLimitPos)
        aAffPos = fm.createFieldMultiply(sAffPos, fm.createFieldSin(phiLimitPos))
        rAffPos = fm.createFieldAdd(dLimit, oAffPos)
        zAffPos = fm.createFieldAdd(zAffDisplPos, aAffPos)

        # Interpolate between anterior and posterior using weights waAff and wpAff
        rAffMapLens = fm.createFieldAdd(fm.createFieldMultiply(waAff, rAffAnt), fm.createFieldMultiply(wpAff, rAffPos))
        zAffMapLens = fm.createFieldAdd(fm.createFieldMultiply(waAff, zAffAnt), fm.createFieldMultiply(wpAff, zAffPos))
        affMapLensCP = fm.createFieldConcatenate([rAffMapLens, theta, zAffMapLens])
        affMapLensCP.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
        affMapLensRC = fm.createFieldCoordinateTransformation(affMapLensCP)
        affMapLensRC.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)

        # Combine spherical and affine mapping
        rGreaterThanRLimit = fm.createFieldGreaterThan(r, rLimit)
        lensRC = fm.createFieldIf(rGreaterThanRLimit, affMapLensRC, sphereMapLensRC)

        # # For debugging
        # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodeiter = nodes.createNodeiterator()
        # node = nodeiter.next()
        # print('node.isValid()',node.isValid())
        # while node.isValid():
            # cache.setNode(node)
            # resultnew, newx = lensRC.evaluateReal(cache, 3)
            # #print(node.getIdentifier(), ':', resultold, oldx, '-->', resultnew, newx)
            # print(newx)
            # node = nodeiter.next()

        fm.endChange()

        return lensRC

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        options['Diameter'] = 1.0

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()

        # generate solidsphere1 model
        MeshType_3d_solidsphere1.generateBaseMesh(region, options)
        sphereCoordinates = getOrCreateCoordinateField(fm)

        # Morph sphere surface to lens surface
        lensRC = cls.sphereToLens(region, options)

        # Assign Field
        fieldassignment = sphereCoordinates.createFieldassignment(lensRC)
        result = fieldassignment.assign()
        # print('fieldassignment', result)

        fm.endChange()

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountRadial = options['Refine number of elements radial']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountRadial)