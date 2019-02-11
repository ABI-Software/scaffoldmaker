"""
Generates a 3-D vertebra model including body, spinous process, transverse processes,
articular processes, vertebral arch and pedicle.
"""

from __future__ import division

import math

from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.torusmesh import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.vector import *


class MeshType_3d_vertebra1(Scaffold_base):
    """
    Generates a 3-D vertebral mesh
    """

    @staticmethod
    def getName():
        return '3D Vertebra 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {'Number of elements around': 4, 'Number of elements up': 4, 'Number of elements through wall': 1,
            'Number of elements radial': 1, 'Major diameter': 2.0, 'Minor diameter': 1.0, 'Height': 1.0,
            'Use cross derivatives': False, 'Refine': False, 'Refine number of elements around': 1,
            'Refine number of elements up': 1, 'Refine number of elements radial': 1}

    @staticmethod
    def getOrderedOptionNames():
        return ['Number of elements around', 'Number of elements radial', 'Number of elements up',
            'Number of elements through wall', 'Major diameter', 'Minor diameter', 'Height', 'Use cross derivatives',
            'Refine', 'Refine number of elements around', 'Refine number of elements up',
            'Refine number of elements radial']

    @staticmethod
    def checkOptions(options):
        for key in ['Number of elements radial', 'Refine number of elements around', 'Refine number of elements up',
            'Refine number of elements radial']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 2:
            options['Number of elements around'] = 2
        if options['Number of elements through wall'] < 1:
            options['Number of elements around'] = 1
        if options['Major diameter'] < 0.0:
            options['Diameter'] = 0.0
        if options['Minor diameter'] < 0.0:
            options['Diameter'] = 0.0
        if options['Height'] < 0.0:
            options['Diameter'] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generates the base tricubic Hermite mesh for this MeshType. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountRadial = options['Number of elements radial']
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']
        axisA = options['Major diameter']
        axisB = options['Minor diameter']
        height = options['Height']
        useCubicHermiteThroughWall = True

        zero = [0.0, 0.0, 0.0]

        """ Get nodes along the centre line """
        cx = [zero, [0.0, 0.0, height]]
        cd1 = [[0.0, 0.0, height / elementsCountUp], [0.0, 0.0, height / elementsCountUp]]
        sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, elementsCountUp)
        # sd2 = interpolateSampleLinear(cd2, se, sxi)
        # sd3 = interpolateSampleLinear(cd3, se, sxi)
        # st2 = interpolateSampleLinear(t2, se, sxi)
        # st2d = interpolateSampleLinear(t2d, se, sxi)
        # st3 = interpolateSampleLinear(t3, se, sxi)
        # st3d = interpolateSampleLinear(t3d, se, sxi)

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)

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

        cache = fm.createFieldcache()

        mesh = fm.findMeshByDimension(3)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        """ Create nodes """
        nodeIdentifier = 1

        # Create bottom node
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, height / elementsCountUp])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, axisB / elementsCountThroughWall, 0.0])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                      [axisA / elementsCountThroughWall, 0.0, 0.0])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

        # Create node along the centre line between bottom and top nodes
        for n2 in range(1, elementsCountUp):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [0, 0, n2 * height / elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, height / elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, axisB / elementsCountThroughWall, 0.0])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                          [axisA / elementsCountThroughWall, 0.0, 0.0])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        # Create top node
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [0, 0, height])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, height/elementsCountUp])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, axisB/elementsCountThroughWall, 0.0])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [axisA/elementsCountThroughWall, 0.0, 0.0])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

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
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
