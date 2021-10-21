"""
Generates a solid sphere (spheroid/ellipsoid in general) using a ShieldMesh of all cube elements,
 with variable numbers of elements across axes and shell directions.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node
from opencmiss.zinc.field import Field
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils import vector
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup


class MeshType_3d_solidsphere2(Scaffold_base):
    """
Generates a solid sphere using a ShieldMesh of all cube elements,
with variable numbers of elements across axes and shell directions.
    """
    # centralPathDefaultScaffoldPackages = {
    #     'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'D2 derivatives': True,
    #             'D3 derivatives': True,
    #             'Length': 3.0,
    #             'Number of elements': 3
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
    #              Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
    #                 [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
    #                 [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
    #                 [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
    #                 [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    #             ])
    #     })
    # }

    @staticmethod
    def getName():
        return '3D Solid Sphere 2'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        # centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            # 'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across axis 1': 4,
            'Number of elements across axis 2': 4,
            'Number of elements across axis 3': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Radius1': 1.0,
            'Radius2': 1.0,
            'Radius3': 1.0,
            'Shell element thickness proportion': 1.0,
            'Octant': False,
            'Hemisphere': False,
            'Full': True,
            'Box derivatives': [-1, 2, 3],
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            # 'Central path',
            'Number of elements across axis 1',
            'Number of elements across axis 2',
            'Number of elements across axis 3',
            # 'Number of elements across shell',
            # 'Number of elements across transition',
            'Radius1',
            'Radius2',
            'Radius3',
            # 'Shell element thickness proportion',
            'Octant',
            'Hemisphere',
            'Full',
            'Box derivatives',
            'Refine',
            'Refine number of elements'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        # if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
        #     options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Octant']:
            dependentChanges = True
            options['Hemisphere'] = False
            options['Full'] = False
        else:
            if options['Hemisphere']:
                dependentChanges = True
                options['Full'] = False
            else:
                options['Full'] = True

        if options['Octant']:
            min1, min2, min3 = 2, 2, 2
            co1, co2, co3 = 0, 0, 0
        elif options['Hemisphere']:
            dependentChanges = True
            min1, min2, min3 = 4, 4, 2
            co1, co2, co3 = 1, 1, 0
        else:
            dependentChanges = True
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
            options['Box derivatives'] = [1, 3, 2]
        if len(options['Box derivatives']) > len(set(options['Box derivatives'])):
            options['Box derivatives'] = [1, 3, 2]

        # if options['Number of elements across transition'] < 1:
        #     options['Number of elements across transition'] = 1
        # Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        # if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
        #     dependentChanges = True
        #     options['Number of elements across shell'] = Rcrit
        #     options['Number of elements across transition'] = 1
        #
        # if options['Shell element thickness proportion'] < 0.15:
        #     options['Shell element thickness proportion'] = 1.0

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        # centralPath = options['Central path']
        elementsCountAcrossAxis1 = options['Number of elements across axis 1']
        elementsCountAcrossAxis2 = options['Number of elements across axis 2']
        elementsCountAcrossAxis3 = options['Number of elements across axis 3']

        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        shellProportion = options['Shell element thickness proportion']
        radius = [options['Radius1'], options['Radius2'], options['Radius3']]
        useCrossDerivatives = options['Use cross derivatives']
        sphereBoxDerivatives = options['Box derivatives']

        if options['Octant']:
            sphere_shape = SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP
        elif options['Hemisphere']:
            sphere_shape = SphereShape.SPHERE_SHAPE_HALF_AAP
        else:
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

        # sphereBoxDerivatives = [1, 3, 2]  # consistent with default derivatives of cylinder mesh.
        centre = [0.0, 0.0, 0.0]
        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        axes = [vector.scaleVector(axis1, radius[0]), vector.scaleVector(axis2, radius[1]), vector.scaleVector(axis3, radius[2])]
        elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]

        sphere1 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                     elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                     sphereShape=sphere_shape, useCrossDerivatives=False, boxDerivatives=sphereBoxDerivatives, meshGroups=meshGroups)

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
