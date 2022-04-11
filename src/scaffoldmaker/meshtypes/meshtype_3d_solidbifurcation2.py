"""
Generates a solid cylinder using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division
import math
import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_1d_stickman1 import MeshType_1d_stickman1
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.bifurcation3d2 import BifurcationMesh
from scaffoldmaker.meshtypes.meshtype_1d_stickman1 import extractPathParametersFromRegion
from scaffoldmaker.utils import vector


class MeshType_3d_solidbifurcation2(Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    centralPathDefaultScaffoldPackages = {
        # 'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
        #     'scaffoldSettings': {
        #         'Coordinate dimensions': 3,
        #         'D2 derivatives': True,
        #         'D3 derivatives': True,
        #         'Length': 3.0,
        #         'Number of elements': 3
        #     }#,
        #     # 'meshEdits': exnodeStringFromNodeValues(
        #     #     [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
        #     #      Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
        #     #         [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        #     #         [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        #     #         [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        #     #         [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
        #     #     ])
        # }),
        'control curves': ScaffoldPackage(MeshType_1d_stickman1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
            }
        }),
        'control curves1': ScaffoldPackage(MeshType_1d_stickman1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
            }
        })
    }

    @staticmethod
    def getName():
        return '3D Solid Bifurcation 2'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['control curves']
        centralPathOption1 = cls.centralPathDefaultScaffoldPackages['control curves1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Central path1': copy.deepcopy(centralPathOption1),
            'Armpit': [1.2, 0.0, 1.0],
            'Torso radius': 1.0,
            'Left arm radius': 1.0,
            'Right arm radius': 1.0,
            'Neck radius': 0.8,
            'Neck radius 2': 0.8,
            'Neck length': 0.4,
            'Neck number of elements': 2,
            'Neck shoulder point': [0.7, 0.0, 2.8],
            'Head length': 1.0,
            'Head number of elements': 5,
            'Head radius': 1.0,
            'Shoulder height': 2.2,
            'Shoulder joint': [1.1, 0.0, 2.4],
            'Shoulder point': [1.0, 0.0, 2.8],
            'Shoulder start': [0.5, 0.0, 2.2],
            'Neck height': 3.6,
            'Right arm angle': 0.0,
            'Left arm angle': 0.0,
            'Right shoulder length': 0.95,
            'Right arm length': 2.0,
            'Right arm number of elements': 5,
            'Right wrist radius': 0.8,
            'Lower torso length': 5.5,
            'Lower torso number of elements': 4,
            'Lower torso radii': [1.3, 1.0],
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Number of elements along': 1,
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
            'Central path1',
            'Armpit',
            'Torso radius',
            'Left arm radius',
            'Right arm radius',
            'Neck radius',
            'Neck radius 2',
            'Neck length',
            'Neck number of elements',
            'Neck shoulder point',
            'Head length',
            'Head number of elements',
            'Head radius',
            'Shoulder joint',
            'Shoulder height',
            'Shoulder point',
            'Shoulder start',
            'Neck height',
            'Right arm angle',
            'Left arm angle',
            'Right shoulder length',
            'Right arm length',
            'Right arm number of elements',
            'Right wrist radius',
            'Lower torso length',
            'Lower torso number of elements',
            'Lower torso radii',
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

        # cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)
        #
        # cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF
        #
        # base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
        #                     elementsCountAcrossTransition,
        #                     shellProportion,
        #                     [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
        #                     cylinderCentralPath.minorRadii[0])
        # cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
        #                          cylinderShape=cylinderShape,
        #                          cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)
        torso_radius = options['Torso radius']
        left_arm_radius = options['Left arm radius']
        right_arm_radius = options['Right arm radius']
        neck_radius = options['Neck radius']
        neck_radius2 = options['Neck radius 2']
        neck_length = options['Neck length']
        neck_number_of_elements = options['Neck number of elements']
        shoulder_height = options['Shoulder height']
        neck_height = options['Neck height']
        right_arm_angle = options['Right arm angle']
        left_arm_angle = options['Left arm angle']
        right_shoulder_length = options['Right shoulder length']
        right_arm_length = options['Right arm length']
        rightArmNumberOfElements = options['Right arm number of elements']
        righ_wrist_radius = options['Right wrist radius']
        shoulder_joint = options['Shoulder joint']
        armpit = options['Armpit']
        neck_shoulder = options['Neck shoulder point']
        shoulder_point = options['Shoulder point']
        shoulder_start = options['Shoulder start']
        head_length = options['Head length']
        head_number_of_elements = options['Head number of elements']
        head_radius = options['Head radius']
        lower_torso_length = options['Lower torso length']
        lower_torso_number_of_elements = options['Lower torso number of elements']
        lower_torso_radii = options['Lower torso radii']



        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE])[0]

        deltacx = vector.vectorRejection(vector.addVectors([cx[3], cx[1]], [1, -1]), [0.0, 1.0, 0.0])
        right_arm_angle = vector.angleBetweenVectors([-1.0, 0.0, 0.0], deltacx)
        if vector.crossproduct3(deltacx, [1.0, 0.0, 0.0])[1] > 0:
            right_arm_angle = -right_arm_angle

        deltacx = vector.vectorRejection(vector.addVectors([cx[7], cx[1]], [1, -1]), [0.0, 1.0, 0.0])
        left_arm_angle = vector.angleBetweenVectors([1.0, 0.0, 0.0], deltacx)
        if vector.crossproduct3(deltacx, [1.0, 0.0, 0.0])[1] > 0:
            left_arm_angle = -left_arm_angle

        # bifurcation1 = BifurcationMesh(fm, coordinates, region, torso_radius, left_arm_radius, right_arm_radius,
        #                                neck_radius, shoulder_height, neck_height, right_arm_angle, right_arm_length,
        #                                shoulder_joint, armpit, neck_shoulder, shoulder_point, shoulder_start)
        bifurcation1 = BifurcationMesh(fm, coordinates, region, torso_radius, left_arm_radius, right_arm_radius,
                                       neck_radius, shoulder_height, neck_height, right_arm_angle,left_arm_angle,
                                       right_shoulder_length, right_arm_length, rightArmNumberOfElements,
                                       righ_wrist_radius, neck_radius2, neck_length, neck_number_of_elements,
                                       head_length, head_number_of_elements, head_radius, armpit, lower_torso_length,
                                       lower_torso_number_of_elements, lower_torso_radii)

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
