"""
Generates a spinal cord scaffold using the solid cylinder elements,
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


class MeshType_3d_spinalcord1(Scaffold_base):
    """
Generates a spinal cord mesh using a solid cylinder
    """
    centralPathDefaultScaffoldPackages = {
        'Spinal cord 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 5.0,
                'Number of elements': 5
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[-0.4620, 0.2860, 2.8640], [-0.4351, -0.2354, -0.4089], [-0.0124, 0.0557, -0.0236],
                     [-0.0053, -0.0051, 0.0104], [0.0542, -0.0086, -0.0487], [0.0075, 0.0420, 0.0744]],
                    [[-0.6560, 0.0958, 2.3777], [0.1213, -0.1049, -0.4937], [-0.0137, 0.0477, -0.0142],
                     [0.0026, -0.0109, 0.0084], [0.0459, 0.0161, 0.0098], [-0.0241, 0.0074, 0.0425]],
                    [[-0.3348, 0.1003, 2.0765], [0.5537, 0.0257, -0.3639], [-0.0079, 0.0344, -0.0065],
                     [0.0075, -0.0019, 0.0046], [0.0090, 0.0095, 0.0393], [-0.0058, -0.0082, 0.0097]],
                    [[0.4236, 0.1551, 1.7489], [0.3797, -0.1409, -0.9449], [0.0028, 0.0540, -0.0077],
                     [0.0033, -0.0010, 0.0012], [0.0623, -0.0016, 0.0114], [0.0170, -0.0040, -0.0239]],
                    [[0.1487, -0.1261, 0.7841], [-0.2629, -0.2135, -1.2602], [-0.0037, 0.0258, -0.0034],
                     [-0.0026, -0.0235, 0.0039], [0.0317, 0.0036, -0.0073], [-0.0255, 0.0023, -0.0061]],
                    [[-0.0608, -0.2126, -0.7292], [0.1978, 0.0575, -1.5941], [-0.0014, 0.0082, -0.0000],
                     [0.0072, -0.0118, 0.0028], [0.0125, 0.0022, 0.0022], [-0.0128, -0.0051, 0.0251]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Spinal Cord 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Spinal cord 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements along': 12,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        if parameterSetName == 'Default':
            parameterSetName = 'Rat 1'
        options['Base parameter set'] = parameterSetName
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements along',
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
        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit

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
        isRat = 'Rat' in parameterSetName

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell, [0.0, 0.0, 0.0],
                            cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
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
