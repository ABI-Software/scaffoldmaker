"""
Generates a 3-D heart arterial root scaffold with semilunar valve,
for attaching to a 6-element-around bicubic-linear orifice.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
import scaffoldmaker.utils.vector as vector
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartarterialroot1(object):
    '''
    Generates a 3-D heart arterial root scaffold with semilunar valve,
    for attaching to a 6-element-around bicubic-linear orifice.
    '''

    @staticmethod
    def getName():
        return '3D Heart Arterial Root 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Unit scale' : 1.0,
            'Outer height' : 0.5,
            'Inner depth' : 0.2,
            'Cusp height' : 0.5,
            'Inner diameter': 1.0,
            'Sinus radial displacement': 0.05,
            'Wall thickness': 0.1,
            'Cusp thickness' : 0.04,
            'Aortic not pulmonary' : True,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Unit scale',
            'Outer height',
            'Inner depth',
            'Cusp height',
            'Inner diameter',
            'Sinus radial displacement',
            'Wall thickness',
            'Cusp thickness',
            'Aortic not pulmonary',
            'Refine',
            'Refine number of elements surface'
        ]

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        for key in [
            'Unit scale',
            'Outer height',
            'Inner depth',
            'Cusp height',
            'Inner diameter',
            'Sinus radial displacement',
            'Wall thickness',
            'Cusp thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'Refine number of elements surface']:
            if options[key] < 1:
                options[key] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options, baseCentre=[ 0.0, 0.0, 0.0 ], axisSide1=[ 0.0, -1.0, 0.0 ], axisUp=[ 0.0, 0.0, 1.0]):
        """
        Generate the base bicubic-linear Hermite mesh. See also generateMesh().
        Optional extra parameters allow centre and axes to be set.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :param baseCentre: Centre of valve on ventriculo-arterial junction.
        :param axisSide: Unit vector in first side direction where angle around starts.
        :param axisUp: Unit vector in outflow direction of valve.
        :return: list of AnnotationGroup
         """
        unitScale = options['Unit scale']
        outerHeight = unitScale*options['Outer height']
        innerDepth = unitScale*options['Inner depth']
        cuspHeight = unitScale*options['Cusp height']
        innerRadius = unitScale*0.5*options['Inner diameter']
        sinusRadialDisplacement = unitScale*options['Sinus radial displacement']
        wallThickness = unitScale*options['Wall thickness']
        cuspThickness = unitScale*options['Cusp thickness']
        aorticNotPulmonary = options['Aortic not pulmonary']
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        if aorticNotPulmonary:
            arterialRootGroup = AnnotationGroup(region, 'root of aorta', FMANumber = 3740, lyphID = 'Lyph ID unknown')
            cuspGroups = [
                AnnotationGroup(region, 'right cusp of aortic valve',     FMANumber = 7252, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'left cusp of aortic valve',      FMANumber = 7251, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'posterior cusp of aortic valve', FMANumber = 7253, lyphID = 'Lyph ID unknown') ]
        else:
            arterialRootGroup = AnnotationGroup(region, 'root of pulmonary trunk', FMANumber = 8612, lyphID = 'Lyph ID unknown')
            cuspGroups = [
                AnnotationGroup(region, 'right cusp of pulmonary valve',    FMANumber = 7250, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'anterior cusp of pulmonary valve', FMANumber = 7249, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'left cusp of pulmonary valve',     FMANumber = 7247, lyphID = 'Lyph ID unknown') ]

        allGroups = [ arterialRootGroup ]  # groups that all elements in scaffold will go in
        annotationGroups = allGroups + cuspGroups

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
        # most nodes in this scaffold do not have a DS3 derivative
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        elementsCountAround = 6
        axisSide2 = vector.crossproduct3(axisUp, axisSide1)
        startRadians = (math.pi/3.0) if aorticNotPulmonary else 0.0
        outerRadius = innerRadius + wallThickness

        # lower points
        lowerx, lowerd1 = getArterialRootLowerPoints(baseCentre, axisSide1, axisSide2, axisUp, innerRadius, outerRadius, innerDepth, sinusRadialDisplacement, startRadians, elementsCountAround)

        # upper points
        topCentre = [ (baseCentre[c] + axisUp[c]*outerHeight) for c in range(3) ]
        ix, id1 = createCirclePoints(topCentre,
            [ axisSide1[c]*innerRadius for c in range(3) ], [ axisSide2[c]*innerRadius for c in range(3) ],
            elementsCountAround, startRadians)
        ox, od1 = createCirclePoints(topCentre,
            [ axisSide1[c]*outerRadius for c in range(3) ], [ axisSide2[c]*outerRadius for c in range(3) ],
            elementsCountAround, startRadians)
        upperx, upperd1 = [ ix, ox ], [ id1, od1 ]

        for n3 in range(2):
            for n1 in range(elementsCountAround):
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lowerx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lowerd1[n3][n1])
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lowerd2[n3][n1])
                nodeIdentifier += 1

        for n3 in range(2):
            for n1 in range(elementsCountAround):
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, upperx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, upperd1[n3][n1])
                #coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, upperd2[n3][n1])
                nodeIdentifier += 1


        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        allMeshGroups = [ allGroup.getMeshGroup(mesh) for allGroup in allGroups ]
        cuspMeshGroups = [ cuspGroup.getMeshGroup(mesh) for cuspGroup in cuspGroups ]

        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives()

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)

        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        scalefactors5hanging = [ -1.0, 0.5, 0.25, 0.125, 0.75 ]

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
        MeshType_3d_heartventricles1.refineMesh(meshrefinement, options)
        numberInXi1 = refineElementsCountSurface
        numberInXi2 = refineElementsCountSurface
        numberInXi3 = 1
        for i in range(0):  # eventually: (24):  # fixed number of elements
            element = meshrefinement._sourceElementiterator.next()
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()


def getArterialRootLowerPoints(baseCentre, axisSide1, axisSide2, axisUp, innerRadius, outerRadius, innerDepth, sinusRadialDisplacement,
        startRadians, elementsCountAround = 6):
    '''
    :param startRadians: Angle from axisSide1 toward axisSide2 where first cusp begins.
    :return: Coordinates, derivatives x[2][6], d1[2][6] for inner:outer and 6 around.
    '''
    assert elementsCountAround == 6, 'getArterialRootLowerPoints.  Only supports 6 elements around'
    ix, id1 = createCirclePoints([ (baseCentre[c] - axisUp[c]*innerDepth) for c in range(3) ],
        [ axisSide1[c]*innerRadius for c in range(3) ], [ axisSide2[c]*innerRadius for c in range(3) ],
        elementsCountAround, startRadians)

    # every second outer points is displaced by sinusRadialDisplacement to make sinus dilatation
    ox = []
    od1 = []
    outerRadiusPlus = outerRadius + sinusRadialDisplacement
    radiansPerElementAround = 2.0*math.pi/elementsCountAround
    radiansAround = startRadians
    for n in range(elementsCountAround):
        radius = outerRadiusPlus if (n % 2) else outerRadius
        rcosRadiansAround = radius*math.cos(radiansAround)
        rsinRadiansAround = radius*math.sin(radiansAround)
        ox.append([ (baseCentre[c] + rcosRadiansAround*axisSide1[c] + rsinRadiansAround*axisSide2[c]) for c in range(3) ])
        od1.append([ radiansPerElementAround*(-rsinRadiansAround*axisSide1[c] + rcosRadiansAround*axisSide2[c]) for c in range(3) ])
        radiansAround += radiansPerElementAround
    # smooth to get derivative in sinus
    pd1 = smoothCubicHermiteDerivativesLine(ox[0:2], od1[0:2], fixStartDerivative=True, fixEndDirection=True)
    od1[1] = pd1[1]
    magSinus = vector.magnitude(pd1[1])
    for n in range(3, elementsCountAround, 2):
        od1[n] = vector.setMagnitude(od1[n], magSinus)

    return [ ix, ox ], [ id1, od1 ]
