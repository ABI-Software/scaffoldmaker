"""
Generates a solid cylinder using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor and length directions.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, ConeBaseProgression, Tapered, \
    CylinderEnds


class MeshType_3d_solidcylinder1(Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor and length directions.
    """
    @staticmethod
    def getName():
        return '3D Solid Cylinder 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements along': 1,
            'Lower half': False,
            'Length': 1.0,
            'Major radius': 1.0,
            'Major radius geometric progression change': True,
            'Major radius end ratio': 0.92,
            'Minor radius': 1.0,
            'Minor radius geometric progression change': True,
            'Minor radius end ratio': 0.92,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements along',
            'Lower half',
            'Length',
            'Major radius',
            'Major radius geometric progression change',
            'Major radius end ratio',
            'Minor radius',
            'Minor radius geometric progression change',
            'Minor radius end ratio',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @staticmethod
    def checkOptions(options):
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1
        if options['Number of elements along'] < 1:
            options['Number of elements along'] = 1

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        full = not options['Lower half']
        length = options['Length']
        majorRadius = options['Major radius']
        majorGeometric = options['Major radius geometric progression change']
        minorGeometric = options['Minor radius geometric progression change']
        minorRadius = options['Minor radius']
        majorRadiusEndRatio = options['Major radius end ratio']
        minorRadiusEndRatio = options['Minor radius end ratio']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]

        if majorGeometric:
            majorRatio = math.pow(majorRadiusEndRatio, 1.0 / elementsCountAlong)
            majorProgression = ConeBaseProgression.GEOMETRIC_PROGRESSION
        else:
            majorRatio = (majorRadiusEndRatio*majorRadius - majorRadius)/elementsCountAlong
            majorProgression = ConeBaseProgression.ARITHMETIC_PROGRESSION

        if minorGeometric:
            minorRatio = math.pow(minorRadiusEndRatio, 1.0/elementsCountAlong)
            minorProgression = ConeBaseProgression.GEOMETRIC_PROGRESSION
        else:
            minorRatio = (minorRadiusEndRatio * minorRadius - minorRadius) / elementsCountAlong
            minorProgression = ConeBaseProgression.ARITHMETIC_PROGRESSION

        radiusChanges = Tapered(majorRatio, majorProgression, minorRatio, minorProgression)
        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, [0.0, 0.0, 0.0],
                            vector.setMagnitude(axis3, length), vector.setMagnitude(axis1, majorRadius), minorRadius)
        cylinder1 = CylinderMesh(fm, coordinates, base, elementsCountAlong,
                                 cylinderShape=cylinderShape, tapered=radiusChanges,
                                 useCrossDerivatives=False)

        annotationGroup = []
        return annotationGroup

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong,
                                                       refineElementsCountAcrossMajor)
