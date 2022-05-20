"""
Generates a whole body scaffold using a mesh of all cube elements,
 with features to control arms, torso, neck and head.
"""

from __future__ import division
import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_1d_stickman1 import MeshType_1d_stickman1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation3d2 import TrifurcationMesh, BranchType, PathNodes, BifurcationMesh
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_solidbifurcation2(Scaffold_base):
    """
Generates a whole body scaffold using a mesh of all cube elements,
 with features to control arms, torso, neck and head.
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
            'Head length': 1.0,
            'Head number of elements': 5,
            'Head radius': 1.0,
            'Left arm radius': 1.0,
            'Lower torso length': 5.5,
            'Lower torso number of elements': 4,
            'Lower torso radii': [1.3, 1.0],
            'Neck height': 3.6,
<<<<<<< HEAD
            'Right arm angle': 0.0,
<<<<<<< HEAD
            'Right arm length': 1.7,
=======
            'Left arm angle': 0.0,
            'Right shoulder length': 0.95,
=======
            'Neck length': 0.4,
            'Neck number of elements': 2,
            'Neck radius': 0.8,
            'Neck radius 2': 0.8,
>>>>>>> 5bc5787 (change the order of options)
            'Right arm length': 2.0,
            'Right arm number of elements': 5,
            'Right arm radius': 1.0,
            'Right shoulder length': 0.95,
            'Right wrist radius': 0.8,
<<<<<<< HEAD
            'Lower torso length': 5.5,
            'Lower torso number of elements': 4,
            'Lower torso radii': [1.3, 1.0],
>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
=======
            'Shoulder height': 2.2,
            'Torso radius': 1.0,
>>>>>>> 5bc5787 (change the order of options)
            'Number of elements across major': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Shell element thickness proportion': 1.0,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Central path1',
            'Armpit',
            'Head length',
            'Head number of elements',
            'Head radius',
            'Left arm radius',
            'Lower torso length',
            'Lower torso number of elements',
            'Lower torso radii',
            'Neck height',
            'Neck length',
            'Neck number of elements',
            'Neck radius',
            'Neck radius 2',
            'Right arm length',
            'Right arm number of elements',
            'Right arm radius',
            'Right shoulder length',
            'Right wrist radius',
            'Shoulder height',
            'Torso radius',
            'Number of elements across major',
            'Number of elements across shell',
            'Number of elements across transition',
            'Shell element thickness proportion',
            'Refine',
            'Refine number of elements across major',
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

        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1
        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across major']-4)//2
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
        elementsCountAcrossMajor = options['Number of elements across major']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        shellProportion = options['Shell element thickness proportion']
        useCrossDerivatives = options['Use cross derivatives']

        torso_radius = options['Torso radius']
        left_arm_radius = options['Left arm radius']
        right_arm_radius = options['Right arm radius']
        neck_radius = options['Neck radius']
        neck_radius2 = options['Neck radius 2']
        neck_length = options['Neck length']
        neck_number_of_elements = options['Neck number of elements']
        shoulder_height = options['Shoulder height']
        neck_height = options['Neck height']
        right_shoulder_length = options['Right shoulder length']
        right_arm_length = options['Right arm length']
        rightArmNumberOfElements = options['Right arm number of elements']
        righ_wrist_radius = options['Right wrist radius']
        armpit = options['Armpit']
        head_length = options['Head length']
        head_number_of_elements = options['Head number of elements']
        head_radius = options['Head radius']
        lower_torso_length = options['Lower torso length']
        lower_torso_number_of_elements = options['Lower torso number of elements']
        lower_torso_radii = options['Lower torso radii']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

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

        trifurcation1 = TrifurcationMesh(fm, coordinates, region, torso_radius, left_arm_radius, right_arm_radius,
                                       neck_radius, shoulder_height, neck_height, right_arm_angle,left_arm_angle,
                                       right_shoulder_length, armpit, [elementsCountAcrossMajor, elementsCountAcrossMajor, 2])

        trifurcation1.create_branch_cylinder([[right_arm_radius]*2, [righ_wrist_radius]*2],
                                            right_arm_length, [elementsCountAcrossMajor, elementsCountAcrossMajor, rightArmNumberOfElements],
                                            branch_type=BranchType.LEFT_ARM)
        trifurcation1.create_branch_cylinder([[right_arm_radius]*2, [righ_wrist_radius]*2],
                                            right_arm_length, [elementsCountAcrossMajor, elementsCountAcrossMajor, rightArmNumberOfElements],
                                            branch_type=BranchType.RIGHT_ARM)
        neck_cylinder = trifurcation1.create_branch_cylinder([[neck_radius2]*2, [neck_radius2]*2], neck_length,
                                                            [elementsCountAcrossMajor,elementsCountAcrossMajor, neck_number_of_elements], branch_type=BranchType.NECK)

        neck_cyliner_shield = neck_cylinder._shield
        pn = PathNodes(neck_cyliner_shield, [[neck_radius2]*2, [head_radius, neck_radius2]],
                       head_length/head_number_of_elements, [elementsCountAcrossMajor,elementsCountAcrossMajor, 1])
        path_list = pn.get_path_list()
        path_list[1][0] = vector.addVectors(
            [path_list[1][0], vector.setMagnitude([0.0, -1.0, 0.0],
                                                  head_length/head_number_of_elements)], [1, 1])

        cw, d1w, d2w = path_list[1][:3]
        d3w = path_list[1][4]
        for ni in range(2, head_number_of_elements + 1):
            cw = vector.addVectors([cw, vector.setMagnitude(d1w, head_length/head_number_of_elements)], [1, 1])
            if ni == 3:
                path_list.append([cw, d1w, vector.scaleVector(d2w, 1.1), [0.0, 0.0, 0.0], d3w, [0.0, 0.0, 0.0]])
            else:
                path_list.append([cw, d1w, d2w, [0.0, 0.0, 0.0], d3w, [0.0, 0.0, 0.0]])

        head_cylinder = trifurcation1.create_branch_cylinder([[neck_radius2] * 2, [head_radius, neck_radius2]],
                                                            head_length/head_number_of_elements,
                                                            [elementsCountAcrossMajor,elementsCountAcrossMajor, head_number_of_elements], path_list=path_list,
                                                            part1=neck_cyliner_shield, branch_type=4)

        cap = trifurcation1.create_branch_cap(head_cylinder, head_radius)

        lower_torso_cylinder = trifurcation1.create_branch_cylinder([[torso_radius]*2, lower_torso_radii],
                                                                   lower_torso_length,
                                                                   [elementsCountAcrossMajor,elementsCountAcrossMajor, lower_torso_number_of_elements],
                                                                   part1=trifurcation1._torso_upper_part, branch_type=4,
                                                                   attach_bottom=False)

        bifurcation1 = BifurcationMesh(fm, coordinates, region, [0, 0, -lower_torso_length], lower_torso_radii)

        # trifurcation1.smooth_all_derivatives()

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
        refineElementsCountAlong = options['Refine number of elements across major']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)
