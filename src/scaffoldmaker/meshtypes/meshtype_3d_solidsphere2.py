"""
Generates a solid sphere (spheroid/ellipsoid in general) using a ShieldMesh of all cube elements,
 with variable numbers of elements across axes and shell directions.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape


class MeshType_3d_solidsphere2(Scaffold_base):
    """
Generates a solid sphere using a ShieldMesh of all cube elements,
with variable numbers of elements across axes and shell directions.
    """

    @staticmethod
    def getName():
        return '3D Solid Sphere 2'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            'Number of elements across axis 1': 4,
            'Number of elements across axis 2': 4,
            'Number of elements across axis 3': 4,
            'Number of elements across shell': 0,
            'Shell element thickness proportion': 1.0,
            'Number of elements across transition': 1,
            'Radius1': 1.0,
            'Radius2': 1.0,
            'Radius3': 1.0,
            'Crop number of elements in direction 1': [0, 0],
            'Crop number of elements in direction 2': [0, 0],
            'Crop number of elements in direction 3': [0, 0],
            'Box derivatives': [1, 2, 3],
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements across axis 1',
            'Number of elements across axis 2',
            'Number of elements across axis 3',
            'Number of elements across shell',
            'Shell element thickness proportion',
            # 'Number of elements across transition',
            'Radius1',
            'Radius2',
            'Radius3',
            'Crop number of elements in direction 1',
            'Crop number of elements in direction 2',
            'Crop number of elements in direction 3',
            'Box derivatives',
            'Refine',
            'Refine number of elements'
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False

        min1, min2, min3 = 4, 4, 4
        co1, co2, co3 = 1, 1, 1

        if options['Number of elements across axis 1'] < min1:
            options['Number of elements across axis 1'] = min1
        if options['Number of elements across axis 2'] < min2:
            options['Number of elements across axis 2'] = min2
        if options['Number of elements across axis 3'] < min3:
            options['Number of elements across axis 3'] = min3

        if options['Number of elements across axis 1'] % 2:
            options['Number of elements across axis 1'] += co1
        if options['Number of elements across axis 2'] % 2:
            options['Number of elements across axis 2'] += co2
        if options['Number of elements across axis 3'] % 2:
            options['Number of elements across axis 3'] += co3

        for radius in ['Radius1', 'Radius2', 'Radius3', ]:
            if options[radius] <= 0:
                options[radius] = 1.0

        if not all(abs(d) in [1, 2, 3] for d in options['Box derivatives']):
            options['Box derivatives'] = [1, 2, 3]
        if len(options['Box derivatives']) > len(set(options['Box derivatives'])):
            options['Box derivatives'] = [1, 2, 3]

        maxelems = [options['Number of elements across axis 1'],
                    options['Number of elements across axis 2'],
                    options['Number of elements across axis 3']]

        cropElements = [
            options['Crop number of elements in direction 1'],
            options['Crop number of elements in direction 2'],
            options['Crop number of elements in direction 3'],
        ]

        for i in range(3):
            for j in [0, 1]:
                if not (1 + options['Number of elements across shell'] < cropElements[i][j] < maxelems[i]):
                    options['Crop number of elements in direction {}'.format(i + 1)][j] = 0

        elementsCount_min = min(options['Number of elements across axis 1'],
                                options['Number of elements across axis 2'],
                                options['Number of elements across axis 3'],)

        max_shell_elements = elementsCount_min//2 - 2

        if options['Number of elements across shell'] > max_shell_elements:
            dependentChanges = True
            options['Number of elements across shell'] = max_shell_elements

        if options['Shell element thickness proportion'] < 0.15:
            options['Shell element thickness proportion'] = 1.0

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        elementsCountAcrossAxis1 = options['Number of elements across axis 1']
        elementsCountAcrossAxis2 = options['Number of elements across axis 2']
        elementsCountAcrossAxis3 = options['Number of elements across axis 3']

        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        shellProportion = options['Shell element thickness proportion']
        radius = [options['Radius1'], options['Radius2'], options['Radius3']]
        useCrossDerivatives = options['Use cross derivatives']

        cropElements = [
            options['Crop number of elements in direction 1'],
            options['Crop number of elements in direction 2'],
            options['Crop number of elements in direction 3'],
        ]
        rangeOfRequiredElements = [
            [cropElements[0][0], elementsCountAcrossAxis1 - cropElements[0][1]],
            [cropElements[1][0], elementsCountAcrossAxis2 - cropElements[1][1]],
            [cropElements[2][0], elementsCountAcrossAxis3 - cropElements[2][1]],
        ]
        sphereBoxDerivatives = [-options['Box derivatives'][0], options['Box derivatives'][1],
                                options['Box derivatives'][2]]  # To make the values more intuitive for the user but
        # consistent with [back, right, up]
        # sphereBoxDerivatives = [1, 3, 2]  # consistent with default derivatives of cylinder mesh.
        # This is the default value that is used for base sphere.
        sphere_shape = SphereShape.SPHERE_SHAPE_FULL

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        mesh = fm.findMeshByDimension(3)
        boxGroup = AnnotationGroup(region, ("box group", ""))
        boxMeshGroup = boxGroup.getMeshGroup(mesh)
        transitionGroup = AnnotationGroup(region, ("transition group", ""))
        transitionMeshGroup = transitionGroup.getMeshGroup(mesh)
        meshGroups = [boxMeshGroup, transitionMeshGroup]
        annotationGroups = [boxGroup, transitionGroup]

        centre = [0.0, 0.0, 0.0]
        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        axes = [vector.scaleVector(axis1, radius[0]),
                vector.scaleVector(axis2, radius[1]),
                vector.scaleVector(axis3, radius[2])]
        elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]

        sphere1 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                             elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False,  meshGroups=meshGroups)

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCount = options['Refine number of elements']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)
