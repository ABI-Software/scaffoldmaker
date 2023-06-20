"""
Generates a whole body scaffold using a mesh of all cube elements,
 with features to control arms, torso, neck and head.
"""

from __future__ import division
import copy

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_1d_stickman1 import MeshType_1d_stickman1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation3d import TrifurcationMesh, BranchType, PathNodes, BifurcationMesh
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_wholebody2(Scaffold_base):
    """
Generates a whole body scaffold using a mesh of all cube elements,
 with features to control arms, torso, neck and head.
    """
    centralPathDefaultScaffoldPackages = {
        'Stickman': ScaffoldPackage(MeshType_1d_stickman1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Arm length': 3.5,
                'Leg length': 6.0,
                'Left arm angle': 0.0,
                'Left leg angle': 1.56,
                'Right arm angle': 0.0,
                'Right leg angle': 1.56
            }
        }),
        'Left arm': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Number of elements': 3
            }  # ,
            # 'meshEdits': exnodeStringFromNodeValues(
            #     [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            #      Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
            #         [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
            #     ])
        }),
        'Left leg': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Number of elements': 3
            }
        }),
        'Right arm': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Number of elements': 3
            }
        }),
        'Right leg': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Number of elements': 3
            }
        }),
    }

    @staticmethod
    def getName():
        return '3D Whole Body 2'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        stickmanOption = cls.centralPathDefaultScaffoldPackages['Stickman']
        leftArmOption = cls.centralPathDefaultScaffoldPackages['Left arm']
        leftLegOption = cls.centralPathDefaultScaffoldPackages['Left leg']
        rightArmOption = cls.centralPathDefaultScaffoldPackages['Right arm']
        rightLegOption = cls.centralPathDefaultScaffoldPackages['Right leg']
        options = {
            'Stickman': copy.deepcopy(stickmanOption),
            'Left arm central path': copy.deepcopy(leftArmOption),
            'Left leg central path': copy.deepcopy(leftLegOption),
            'Right arm central path': copy.deepcopy(rightArmOption),
            'Right leg central path': copy.deepcopy(rightLegOption),
            'Armpit': [1.6, 0.0, 1.2],
            'Head length': 1.5,
            'Head number of elements': 5,
            'Head radius': 1.0,
            'Left arm radius': 0.6,
            'Lower torso length': 2.6,
            'Lower torso number of elements': 4,
            'Lower torso radii': [1.3, 1.6],
            'Neck height': 3.3,
            'Neck length': 0.2,
            'Neck number of elements': 1,
            'Neck radius': 0.8,
            'Neck radius 2': 0.8,
            'Right arm length': 3.5,
            'Right arm number of elements': 5,
            'Right arm radius': 0.6,
            'Right shoulder length': 1.4,
            'Right wrist radius': 0.4,
            'Shoulder height': 2.1,
            'Torso radius': 1.3,
            'Number of elements across major': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Shell element thickness proportion': 1.0,
            'Pre-fit configuration': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Stickman',
            'Left arm central path',
            'Left leg central path',
            'Right arm central path',
            'Right leg central path',
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
            'Refine',
            'Refine number of elements across major',
                ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Stickman':
            return [MeshType_1d_stickman1]
        elif optionName in ['Left arm central path', 'Left leg central path',
                            'Right arm central path', 'Right leg central path']:
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Stickman':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        elif optionName in ['Left arm central path', 'Left leg central path',
                            'Right arm central path', 'Right leg central path']:
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
        if optionName == 'Stickman':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        elif optionName in ['Left arm central path', 'Left leg central path',
                            'Right arm central path', 'Right leg central path']:
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Stickman'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Stickman'):
            options['Stickman'] = cls.getOptionScaffoldPackage('Stickman', MeshType_1d_stickman1)
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

        stickman = options['Stickman']
        left_arm_central_path = options['Left arm central path']
        left_leg_central_path = options['Left leg central path']
        right_arm_central_path = options['Right arm central path']
        right_leg_central_path = options['Right leg central path']
        elementsCountAcrossMajor = options['Number of elements across major']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        shellProportion = options['Shell element thickness proportion']
        useCrossDerivatives = options['Use cross derivatives']

        armpit = options['Armpit']
        head_length = options['Head length']
        head_number_of_elements = options['Head number of elements']
        head_radius = options['Head radius']
        left_arm_radius = options['Left arm radius']
        lower_torso_length = options['Lower torso length']
        lower_torso_number_of_elements = options['Lower torso number of elements']
        lower_torso_radii = options['Lower torso radii']
        neck_height = options['Neck height']
        neck_length = options['Neck length']
        neck_number_of_elements = options['Neck number of elements']
        neck_radius = options['Neck radius']
        neck_radius2 = options['Neck radius 2']
        righ_wrist_radius = options['Right wrist radius']
        rightArmNumberOfElements = options['Right arm number of elements']
        right_arm_length = options['Right arm length']
        right_arm_radius = options['Right arm radius']
        right_shoulder_length = options['Right shoulder length']
        shoulder_height = options['Shoulder height']
        torso_radius = options['Torso radius']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        tmpRegion = region.createRegion()
        scaffoldSettings = stickman.getScaffoldSettings()
        scaffoldSettings['Arm length'] = right_arm_length
        stickman.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE])[0]

        def limb_angle(p1, p2):
            p1p2 = vector.vectorRejection(vector.addVectors([p2, p1], [1, -1]), [0.0, 1.0, 0.0])
            angle = vector.angleBetweenVectors([p2[0], 0.0, 0.0], p1p2)
            if vector.crossproduct3(p1p2, [1.0, 0.0, 0.0])[1] > 0:
                angle = -angle
            return angle

        # get the angles from the stickman
        right_arm_angle = limb_angle(cx[1], cx[3])
        left_arm_angle = limb_angle(cx[1], cx[5])
        right_leg_angle = limb_angle(cx[2], cx[4])
        left_leg_angle = limb_angle(cx[2], cx[6])

        # set the angles and regenerate the stickman with the new settings.
        if stickman.getMeshEdits():
            scaffoldSettings['Left arm angle'] = left_arm_angle
            scaffoldSettings['Left leg angle'] = left_leg_angle
            scaffoldSettings['Right arm angle'] = right_arm_angle
            scaffoldSettings['Right leg angle'] = right_leg_angle
            stickman.setMeshEdits(None)
            stickman.generate(tmpRegion)

        # create torso, head and arms
        trifurcation1 = TrifurcationMesh(fm, coordinates, region, torso_radius, left_arm_radius, right_arm_radius,
                                         neck_radius, shoulder_height, neck_height, right_arm_angle, left_arm_angle,
                                         right_shoulder_length, armpit,
                                         [elementsCountAcrossMajor, elementsCountAcrossMajor, 2])

        trifurcation1.create_branch_cylinder([[right_arm_radius]*2, [righ_wrist_radius]*2], right_arm_length,
                                             [elementsCountAcrossMajor, elementsCountAcrossMajor,
                                              rightArmNumberOfElements],
                                             branch_type=BranchType.LEFT_ARM)
        trifurcation1.create_branch_cylinder([[right_arm_radius]*2, [righ_wrist_radius]*2], right_arm_length,
                                             [elementsCountAcrossMajor, elementsCountAcrossMajor,
                                              rightArmNumberOfElements],
                                             branch_type=BranchType.RIGHT_ARM)
        neck_cylinder = trifurcation1.create_branch_cylinder([[neck_radius2]*2, [neck_radius2]*2], neck_length,
                                                             [elementsCountAcrossMajor, elementsCountAcrossMajor,
                                                              neck_number_of_elements],
                                                             branch_type=BranchType.NECK)

        # neck
        neck_cylinder_shield = neck_cylinder.getShield()
        pn = PathNodes(neck_cylinder_shield, [[neck_radius2]*2, [head_radius, neck_radius2]],
                       head_length/head_number_of_elements, [elementsCountAcrossMajor, elementsCountAcrossMajor, 1])
        path_list = pn.get_path_list()
        path_list[1][0] = vector.addVectors(
            [path_list[1][0], vector.setMagnitude([0.0, -1.0, 0.0],
                                                  head_length/head_number_of_elements)], [1, 1])

        # extend the neck and create the head using a cylinder with its central path
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
                                                             [elementsCountAcrossMajor, elementsCountAcrossMajor,
                                                              head_number_of_elements], path_list=path_list,
                                                             part1=neck_cylinder_shield, branch_type=4)

        cap = trifurcation1.create_branch_cap(head_cylinder, head_radius)

        lower_torso_cylinder = trifurcation1.create_branch_cylinder([[torso_radius]*2, lower_torso_radii],
                                                                    lower_torso_length,
                                                                    [elementsCountAcrossMajor, elementsCountAcrossMajor,
                                                                     lower_torso_number_of_elements],
                                                                    part1=trifurcation1._torso_upper_part,
                                                                    branch_type=4, attach_bottom=False)

        # create the legs
        bifurcation1 = BifurcationMesh(fm, coordinates, region, [0, 0, -lower_torso_length], lower_torso_radii,
                                       right_leg_angle, left_leg_angle,
                                       left_leg_central_path, right_leg_central_path, part1=lower_torso_cylinder)

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
