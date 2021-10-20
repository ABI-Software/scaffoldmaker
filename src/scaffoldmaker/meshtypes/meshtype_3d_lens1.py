"""
Generates a 3-D lens mesh with variable numbers of elements
around, up the central axis, and radially.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from scaffoldmaker.meshtypes.meshtype_3d_solidsphere1 import MeshType_3d_solidsphere1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_lens1(Scaffold_base):
    '''
    Generates a 3-D spherical lens mesh with variable numbers
    of elements around, up the central axis, and radially.
    The spherical lens is created by a function morphing a
    sphere to a lens
    '''
    @staticmethod
    def getName():
        return '3D Lens 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = MeshType_3d_solidsphere1.getDefaultOptions(parameterSetName)
        options['Number of elements around'] = 8
        optionsLens = {
            'Axial thickness' : 4.0,
            'Anterior radius of curvature' : 10.0,
            'Posterior radius of curvature' : 6.0,
            'Sphere spherical radius fraction' : 0.8,
            'Lens spherical radius fraction' : 0.7
        }
        options.update(optionsLens)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_solidsphere1.getOrderedOptionNames()
        optionNames.remove('Diameter')
        for optionName in [
            'Lens spherical radius fraction',
            'Sphere spherical radius fraction',
            'Posterior radius of curvature',
            'Anterior radius of curvature',
            'Axial thickness']:
            optionNames.insert(3, optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_solidsphere1.checkOptions(options)
        for key in [
            'Axial thickness',
            'Sphere spherical radius fraction',
            'Lens spherical radius fraction']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'Anterior radius of curvature',
            'Posterior radius of curvature']:
            if options[key] < options['Axial thickness']*0.5:
                options[key] = options['Axial thickness']*0.5
        for key in [
            'Sphere spherical radius fraction',
            'Lens spherical radius fraction']:
            if options[key] > 1.0:
                options[key] = 1.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        options['Diameter'] = 1.0
        radiusSphere = options['Diameter']*0.5
        radiusAnt = options['Anterior radius of curvature']
        radiusPos = options['Posterior radius of curvature']
        lensThickness = options['Axial thickness']
        sphereSphericalRadiusFraction = options['Sphere spherical radius fraction']
        lensSphericalRadiusFraction = options['Lens spherical radius fraction']

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()

        # generate solidsphere with unit diameter
        MeshType_3d_solidsphere1.generateBaseMesh(region, options)
        sphereCoordinates = findOrCreateFieldCoordinates(fm)

        # Morph sphere surface to lens surface
        lensRC = getSphereToLensCoordinates(sphereCoordinates, radiusSphere, radiusAnt, radiusPos, lensThickness,
            sphereSphericalRadiusFraction, lensSphericalRadiusFraction)

        # Assign Field
        fieldassignment = sphereCoordinates.createFieldassignment(lensRC)
        result = fieldassignment.assign()

        fm.endChange()
        return []

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountRadial = options['Refine number of elements radial']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountRadial)


def getSphereToLensCoordinates(sphereCoordinates, radiusSphere, radiusAnt, radiusPos, lensThickness,
        sphereSphericalRadiusFraction = 0.8, lensSphericalRadiusFraction = 0.7):
    """
    Define a field giving positions on a spherical lens as functions of positions
    in the sphere. Map coordinates of the sphere within a width limit olim to arcs
    with radii rAnt and rPos, respectively. Elements maintain constant size radially
    and match up at dlim. Outside of rlim, affine transformation applied to transform
    sphere coordinates to be tangential to morphed spherical surfaces
    param sphereCoordinates: Coordinate field defined over sphere
    param radiusSphere: Radius of solid sphere
    param radiusAnt: Radius of curvature on anterior lens surface
    param radiusPos: Radius of curvature on posterior lens surface
    param lensThickness: Axial thickness of lens
    param sphereSphericalRadiusFraction: Proportion of sphere radius mapped to spherical
    part of lens
    param lensSphericalRadiusFraction: Proportion of lens radius with spherical shape
    return: Zinc Field giving lens coordinates
    """

    fm = sphereCoordinates.getFieldmodule()

    # Estimate dLim
    rMax = max(radiusAnt, radiusPos)
    rMin = min(radiusAnt, radiusPos)
    zMax = rMax - math.sqrt(rMax*rMax - rMin*rMin)
    if zMax > lensThickness - rMin:
        thicknessMin = lensThickness*(2*rMax - lensThickness)/(2*(rMax + rMin - lensThickness))
        dLim = math.sqrt(rMin*rMin - (rMin - thicknessMin)*(rMin - thicknessMin))
    else:
        dLim = rMin
    dLimitConstant = dLim*lensSphericalRadiusFraction

    # Define constants
    half = fm.createFieldConstant([0.5])
    one = fm.createFieldConstant([1.0])
    rSphere = fm.createFieldConstant([radiusSphere])
    rAnt = fm.createFieldConstant([radiusAnt])
    rPos = fm.createFieldConstant([radiusPos])
    halfLensThickness = fm.createFieldConstant([lensThickness*0.5])

    psiLimit = fm.createFieldConstant([math.asin(sphereSphericalRadiusFraction)])
    rLimit = fm.createFieldConstant([sphereSphericalRadiusFraction*radiusSphere])
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

    return lensRC
