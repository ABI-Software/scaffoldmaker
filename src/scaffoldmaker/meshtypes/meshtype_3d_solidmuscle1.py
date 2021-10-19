"""
Generates a one-headed muscle using a CylinderMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node


class MeshType_3d_solidmuscle1(Scaffold_base):
    """
Generates a one-headed muscle using a CylinderMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    centralPathDefaultScaffoldPackages = {
        'Brachioradialis 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 6
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[-0.27, 0.44, 0.87], [6.76, 7.80, 44.13], [4.58, -1.00, -0.40], [0.00, 0.00, 0.00], [0.13, 1.52, -0.27], [0.00, 0.00, 0.00]],
                    [[3.25, 5.21, 45.57], [0.39, 1.74, 45.29], [4.61, -0.73, -0.01], [0.00, 0.00, 0.00], [0.25, 2.60, -0.10], [0.00, 0.00, 0.00]],
                    [[0.46, 3.96, 90.94], [-0.65, -2.04, 59.65], [7.90, -0.75, 0.05], [0.00, 0.00, 0.00], [0.37, 3.98, 0.14], [0.00, 0.00, 0.00]],
                    [[3.28, 1.02, 164.51], [5.43, 7.26, 72.36], [16.13, -2.18, -0.95], [0.00, 0.00, 0.00], [0.96, 7.43, -0.81], [0.00, 0.00, 0.00]],
                    [[11.03, 18.03, 234.00], [11.35, 10.96, 50.55], [8.23, -1.76, -1.46], [0.00, 0.00, 0.00], [1.09, 6.46, -1.65], [0.00, 0.00, 0.00]],
                    [[21.86, 23.91, 265.44], [9.66, 3.35, 33.17], [3.99, -1.33, -1.01], [0.00, 0.00, 0.00], [1.16, 6.60, -0.99], [0.00, 0.00, 0.00]],
                    [[30.07, 24.60, 299.90], [6.77, -1.97, 35.76], [1.86, -0.92, -0.40], [0.00, 0.00, 0.00], [0.53, 2.86, 0.05], [0.00, 0.00, 0.00]]
                ])
        }),
        'Long head of right biceps femoris 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[1.11, 0.78, 0.89], [-16.00, 0.00, 44.00], [-2.00, 5.00, -1.00], [0.00, 0.00, 0.00], [-4.00, -1.00, -1.00], [0.00, 0.00, 0.00]],
                    [[-10.89, 4.78, 45.89], [-6.00, 7.00, 45.00], [-2.00, 8.00, -2.00], [0.00, 0.00, 0.00], [-5.00, -2.00, -1.00], [0.00, 0.00, 0.00]],
                    [[-11.89, 14.78, 89.89], [-8.00, 6.00, 48.00], [-1.00, 15.00, -2.00], [0.00, 0.00, 0.00], [-14.00, -1.00, -2.00], [0.00, 0.00, 0.00]],
                    [[-26.89, 16.78, 139.89], [-12.00, 2.00, 49.00], [0.00, 25.00, -1.00], [0.00, 0.00, 0.00], [-16.00, 0.00, -4.00], [0.00, 0.00, 0.00]],
                    [[-35.89, 18.78, 187.89], [-8.00, 3.00, 49.00], [-1.00, 26.00, -1.00], [0.00, 0.00, 0.00], [-20.00, -1.00, -3.00], [0.00, 0.00, 0.00]],
                    [[-42.89, 21.78, 237.89], [-7.00, 5.00, 48.00], [-1.00, 24.00, -2.00], [0.00, 0.00, 0.00], [-19.00, -1.00, -3.00], [0.00, 0.00, 0.00]],
                    [[-48.89, 28.78, 284.89], [-5.00, 10.00, 47.00], [-2.00, 13.00, -3.00], [0.00, 0.00, 0.00], [-14.00, -2.00, -1.00], [0.00, 0.00, 0.00]],
                    [[-52.89, 42.78, 330.89], [-3.00, 17.00, 46.00], [-2.00, 6.00, -2.00], [0.00, 0.00, 0.00], [-8.00, -3.00, 1.00], [0.00, 0.00, 0.00]]
                ])
        }),

    }

    @staticmethod
    def getName():
        return '3D Solid Muscle 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Brachioradialis 1',
            'Long head of right biceps femoris 1',
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Brachioradialis 1'

        if 'Brachioradialis 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Brachioradialis 1']
        elif 'Long head of right biceps femoris 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Long head of right biceps femoris 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Number of elements along': 6,
            'Shell element thickness proportion': 1.0,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements across transition',
            'Number of elements along',
            'Shell element thickness proportion',
            'Lower half',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
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
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
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
        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1
        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1

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

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        elementsCountAlong = options['Number of elements along']
        shellProportion = options['Shell element thickness proportion']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            elementsCountAcrossTransition,
                            shellProportion,
                            [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

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
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)
