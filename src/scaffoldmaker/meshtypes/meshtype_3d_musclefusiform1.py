"""
Generates a fusiform muscle using a CylinderMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division

import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class MeshType_3d_musclefusiform1(Scaffold_base):
    """
Generates a fusiform muscle using a CylinderMesh of all cube elements,
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
                    [[-0.27, 0.44, 0.87], [6.61, 7.71, 43.82], [4.56, -1.02, -0.51], [-1.56, 0.44, 0.72], [0.29, 1.48, -0.30], [0.18, 0.93, 0.18]],
                    [[3.25, 5.21, 45.57], [0.37, 1.78, 45.29], [4.61, -0.73, -0.01], [1.66, 0.13, 0.29], [0.40, 2.58, -0.10], [0.04, 1.24, 0.22]],
                    [[0.46, 3.96, 90.94], [-0.69, -2.01, 59.61], [7.90, -0.75, 0.06], [5.17, -0.55, -0.35], [0.37, 3.97, 0.14], [0.20, 2.18, -0.21]],
                    [[3.28, 1.02, 164.51], [5.37, 7.23, 72.38], [16.12, -2.18, -0.99], [0.07, -0.49, -0.76], [0.95, 7.43, -0.81], [0.35, 1.21, -0.89]],
                    [[11.03, 18.03, 234.00], [11.41, 10.94, 50.55], [8.22, -1.76, -1.46], [-5.41, 0.42, 0.14], [1.09, 6.46, -1.64], [0.55, -0.34, 0.05]],
                    [[21.86, 23.91, 265.44], [9.61, 3.36, 33.13], [3.98, -1.33, -1.02], [-3.21, 0.42, 0.52], [1.83, 6.41, -1.18], [0.10, -1.87, 0.76]],
                    [[30.07, 24.60, 299.90], [6.77, -1.97, 35.60], [1.85, -0.92, -0.40], [-1.04, 0.40, 0.71], [1.27, 2.61, -0.09], [-1.23, -5.71, 1.40]]
                ])
        }),
        'Biceps femoris 1': ScaffoldPackage(MeshType_1d_path1, {
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
                    [[1.11, 0.78, 0.89], [-17.32, 0.85, 44.44], [-2.04, 5.00, -0.89], [-0.58, 1.13, -0.73], [-3.62, -1.72, -1.38], [1.97, 0.30, 2.07]],
                    [[-10.89, 4.78, 45.89], [-6.48, 7.14, 45.04], [-2.06, 8.08, -1.58], [0.54, 5.02, -0.64], [-5.26, -1.44, -0.53], [-5.24, 0.26, -0.37]],
                    [[-11.89, 14.78, 89.89], [-7.65, 6.42, 47.75], [-0.97, 14.98, -2.17], [1.04, 8.35, 0.22], [-13.97, -1.21, -2.08], [-5.63, 0.61, -1.67]],
                    [[-26.89, 16.78, 139.89], [-11.94, 2.01, 49.14], [0.00, 25.00, -1.02], [-0.01, 5.35, 0.32], [-16.03, -0.16, -3.89], [-3.03, 0.13, -0.53]],
                    [[-35.89, 18.78, 187.89], [-8.02, 2.49, 49.02], [-0.92, 25.98, -1.47], [-0.47, -0.50, -0.79], [-19.97, -0.89, -3.22], [-1.57, -0.42, 0.72]],
                    [[-42.89, 21.78, 237.89], [-6.50, 5.06, 48.56], [-0.91, 23.94, -2.62], [-0.55, -6.63, -0.82], [-19.08, -0.99, -2.45], [3.07, -0.76, 1.14]],
                    [[-48.89, 28.78, 284.89], [-5.02, 10.53, 46.68], [-1.99, 12.97, -3.14], [-0.54, -9.05, 0.11], [-13.94, -2.37, -0.97], [5.44, -0.79, 1.45]],
                    [[-52.89, 42.78, 330.89], [-2.97, 17.41, 45.15], [-1.98, 5.86, -2.39], [0.55, -5.16, 1.38], [-8.19, -2.58, 0.46], [6.06, 0.38, 1.39]]
                ])
        }),

    }

    @staticmethod
    def getName():
        return '3D Muscle Fusiform 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Brachioradialis 1',
            'Biceps femoris 1',
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Brachioradialis 1'

        if 'Brachioradialis 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Brachioradialis 1']
        elif 'Biceps femoris 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Biceps femoris 1']

        options = {
            'Base parameter set': parameterSetName,
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Number of elements along': 10,
            'Shell element thickness proportion': 1.0,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across': 4,
            'Refine number of elements along': 4
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
            'Refine',
            'Refine number of elements across',
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
        parameterSetName = options['Base parameter set']

        centralPath = options['Central path']
        full = True
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

        mesh = fm.findMeshByDimension(3)
        annotationGroups = []

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

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcross = options['Refine number of elements across']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcross, refineElementsCountAlong, refineElementsCountAcross)
