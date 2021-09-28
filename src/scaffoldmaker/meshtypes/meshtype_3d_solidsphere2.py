"""
Generates a solid sphere using a ShieldMesh of all cube elements,
 with variable numbers of elements across axes and shell directions.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
# from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node
from opencmiss.zinc.field import Field
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D, ShieldShape3D



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
            'Number of elements across axis 1': 2,
            'Number of elements across axis 2': 2,
            'Number of elements across axis 3': 2,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            # 'Number of elements along': 1,
            'Derivatives configuration1': True,
            'Derivatives configuration2': False,
            'Shell element thickness proportion': 1.0,
            'Octant': True,
            'Hemisphere': False,
            'Full': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            # 'Central path',
            'Number of elements across axis 1',
            'Number of elements across axis 2',
            'Number of elements across axis 3',
            'Number of elements across shell',
            'Number of elements across transition',
            # 'Number of elements along',
            'Derivatives configuration1',
            'Derivatives configuration2',
            'Shell element thickness proportion',
            'Octant',
            'Hemisphere',
            'Full',
            'Refine',
            'Refine number of elements across'
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

        # if options['Number of elements across major'] < 4:
        #     options['Number of elements across major'] = 4
        # if options['Number of elements across major'] % 2:
        #     options['Number of elements across major'] += 1
        #
        # if options['Number of elements across minor'] < 4:
        #     options['Number of elements across minor'] = 4
        # if options['Number of elements across minor'] % 2:
        #     options['Number of elements across minor'] += 1
        # if options['Number of elements along'] < 1:
        #     options['Number of elements along'] = 1
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
        # full = not options['Lower half']
        octant = options['Octant']
        hemisphere = options['Hemisphere']
        full = options['Full']
        elementsCountAcrossAxis1 = options['Number of elements across axis 1']
        elementsCountAcrossAxis2 = options['Number of elements across axis 2']
        elementsCountAcrossAxis3 = options['Number of elements across axis 3']
        # if not full:
        #     elementsCountAcrossMajor //= 2
        # elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        # elementsCountAlong = options['Number of elements along']
        shellProportion = options['Shell element thickness proportion']
        useCrossDerivatives = options['Use cross derivatives']
        first = options['Derivatives configuration1']
        second = options['Derivatives configuration2']
        if first:
            boxMapping = [1, 3, 2]
        if second and not first:
            boxMapping = [-1, 2, 3]

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        centre = [0.0, 0.0, 0.0]

        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        axes = [axis1, axis2, axis3]
        elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]

        sphere1 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                     elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                     sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

        if (hemisphere or full) and not octant:
            axis1 = [0.0, -1.0, 0.0]
            axis2 = [1.0, 0.0, 0.0]
            axis3 = [0.0, 0.0, 1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis2, elementsCountAcrossAxis1, elementsCountAcrossAxis3]
            boxMapping = [3, -1, 2]
            sphere2 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

            axis1 = [-1.0, 0.0, 0.0]
            axis2 = [0.0, -1.0, 0.0]
            axis3 = [0.0, 0.0, 1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]
            boxMapping = [-1, -3, 2]
            sphere3 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

            axis1 = [0.0, 1.0, 0.0]
            axis2 = [-1.0, 0.0, 0.0]
            axis3 = [0.0, 0.0, 1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis2, elementsCountAcrossAxis1, elementsCountAcrossAxis3]
            boxMapping = [-3, 1, 2]
            sphere4 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

        if full and not octant and not hemisphere:
            axis1 = [0.0, 1.0, 0.0]
            axis2 = [1.0, 0.0, 0.0]
            axis3 = [0.0, 0.0, -1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis2, elementsCountAcrossAxis1, elementsCountAcrossAxis3]
            boxMapping = [-3, -1, -2]
            sphere5 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

            axis2 = [0.0, -1.0, 0.0]
            axis1 = [1.0, 0.0, 0.0]
            axis3 = [0.0, 0.0, -1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]
            boxMapping = [1, -3, -2]
            sphere6 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

            axis2 = [-1.0, 0.0, 0.0]
            axis1 = [0.0, -1.0, 0.0]
            axis3 = [0.0, 0.0, -1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis2, elementsCountAcrossAxis1, elementsCountAcrossAxis3]
            boxMapping = [3, 1, -2]
            sphere7 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

            axis2 = [0.0, 1.0, 0.0]
            axis1 = [-1.0, 0.0, 0.0]
            axis3 = [0.0, 0.0, -1.0]
            axes = [axis1, axis2, axis3]
            elementsCountAcross = [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]
            boxMapping = [-1, 3, -2]
            sphere8 = SphereMesh(fm, coordinates, centre, axes, elementsCountAcross,
                         elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                         sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxMapping=boxMapping)

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
        refineElementsCountAcross = options['Refine number of elements across']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcross, refineElementsCountAcross, refineElementsCountAcross)
