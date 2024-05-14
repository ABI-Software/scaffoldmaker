"""
Generates a solid bone using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division

import copy

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldStoredString, \
    findOrCreateFieldStoredMeshLocation
from cmlibs.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.bone_terms import get_bone_term

class MeshType_3d_bone1 (Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    centralPathDefaultScaffoldPackages = {
        'Bone 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (2, [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (3, [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (4, [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
                ])
        }),
        'Ulna 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[ 5.356588e-01, 1.344625e+00,-1.057989e+00], [-2.179484e-01,-6.019707e-01, 8.118885e-01], [ 1.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 1.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (2, [[ 5.342651e-01, 2.042545e-01, 7.154451e-01], [ 0.000000e+00, 0.000000e+00, 1.000000e+00], [ 1.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 1.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (3, [[ 4.256270e-01,-7.595472e-01, 2.423991e+00], [ 0.000000e+00, 0.000000e+00, 1.000000e+00], [ 1.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 1.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (4, [[ 8.636938e-01,-2.271568e+00, 3.615462e+00], [ 3.065797e-01,-1.318752e+00, 8.066436e-01], [ 1.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 0.000000e+00, 1.000000e+00, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]])
                ])
        }),
        'Radius 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[-1.046980e+00,-6.643628e-01,-1.123687e+00], [ 0.000000e+00, 0.000000e+00, 1.000000e+00], [ 2.606645e-01, 3.753671e-02, 1.551701e-02], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [-6.214456e-02, 2.659024e-01, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (2, [[-9.877670e-01,-7.028367e-01,-4.943481e-01], [ 8.919182e-02,-5.293884e-02, 7.619579e-01], [ 2.444035e-01, 9.618276e-03, 1.657146e-02], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [-4.259481e-02, 2.119403e-01, 3.984784e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (3, [[-9.341713e-01,-6.580989e-01, 1.537421e-01], [ 1.589471e-01,-5.608755e-02, 8.047271e-01], [ 2.196727e-01, 2.198117e-02, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [-2.156681e-02, 2.155317e-01, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (4, [[-8.515940e-01,-6.576398e-01, 6.944961e-01], [ 6.039198e-03, 5.546106e-02, 2.692810e-01], [ 2.174544e-01, 6.164221e-02, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [-6.445803e-02, 2.273877e-01, 0.000000e+00], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]])
                ])
        }),
        'Tibia 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[ 5.478874e-02,-2.031468e-02,-9.043539e+00], [ 1.707559e-01,-5.112576e-02, 5.609759e+00], [ 7.943801e-01,-7.964818e-02,-7.322865e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 9.087504e-02, 1.136004e+00, 7.587069e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (2, [[ 6.058305e-02,-9.069170e-02,-3.895807e+00], [-1.109146e-01, 5.163967e-02, 5.480916e+00], [ 5.653736e-01,-1.230512e-01, 1.260054e-02], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 1.389491e-01, 6.380904e-01,-3.200067e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (3, [[-1.877598e-01, 1.107110e-02, 1.514034e+00], [-4.509590e-02, 6.698397e-02, 5.530504e+00], [ 7.153731e-01,-1.353970e-01, 7.473066e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 1.388781e-01, 7.333377e-01,-7.749571e-03], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]]),
                    (4, [[-2.093406e-02, 4.163051e-02, 7.155221e+00], [ 3.783990e-01,-5.859759e-03, 5.746579e+00], [ 1.346095e+00,-1.457858e-01,-8.878590e-02], [ 0.000000e+00, 0.000000e+00, 0.000000e+00], [ 2.337453e-01, 2.166286e+00,-1.318264e-02], [ 0.000000e+00, 0.000000e+00, 0.000000e+00]])
                ])
        }),

    }

    @staticmethod
    def getName():
        return '3D Bone 1'
    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Bone 1',
            'Ulna 1',
            'Radius 1',
            'Tibia 1',
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):


        if parameterSetName == 'Default':
            parameterSetName = 'Bone 1'

        if 'Bone 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Bone 1']
        elif 'Ulna 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Ulna 1']
        elif 'Radius 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Radius 1']
        elif 'Tibia 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Tibia 1']





        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 8,
            'Number of elements across minor': 8,
            'Number of elements across shell': 2,
            'Number of elements across transition': 1,
            'Number of elements along': 4,
            'Shell element thickness proportion': 1.0,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1,
            'Upper scale': 1.0,
            'Upper scale_Z': 1.0,
            'Lower scale': 1.0,
            'Lower scale_Z': 1.0
        }

        options['Base parameter set'] = parameterSetName
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
            # 'Lower half',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along',
            'Upper scale',
            'Upper scale_Z',
            'Lower scale',
            'Lower scale_Z'
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

        if options['Number of elements across major'] < 6:
            options['Number of elements across major'] = 6
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        #if options['Number of elements across minor'] < 6:
        #    options['Number of elements across minor'] = 6
        #if options['Number of elements across minor'] % 2:
        #    options['Number of elements across minor'] += 1
        options['Number of elements across minor'] = options['Number of elements across major' ]
        #if options['Number of elements across minor'] != options['Number of elements across major']:
        #    options['number of elements across minor'] = options['number of elements across major']
        if options['Number of elements along'] < 1:
            options['Number of elements along'] = 1
        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1

        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1
            if options['Number of elements across major'] == 8:
                options['Number of elements across shell'] = 2

        if options['Shell element thickness proportion'] < 0.15:
            options['Shell element thickness proportion'] = 1.0

        parameterSetName = options['Base parameter set']
        isRadius = 'Radius 1' in parameterSetName
        isUlna = 'Ulna 1' in parameterSetName
        isTibia = 'Tibia 1' in parameterSetName

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
        isRadius = 'Radius 1' in parameterSetName
        isUlna = 'Ulna 1' in parameterSetName
        isTibia = 'Tibia 1' in parameterSetName
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

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        sphere_shape = SphereShape.SPHERE_SHAPE_FULL
        sphere_base = cylinder1._ellipses[0]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = options['Lower scale_Z']
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                vector.setMagnitude(vector.crossproduct3(sphere_base.majorAxis, sphere_base.minorAxis),
                                    sphere_radius_3)]
        elementsCountAcross = [cylinder1._elementsCountAcrossMajor, cylinder1._elementsCountAcrossMinor,
                               cylinder1._elementsCountAcrossMajor]
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, -1]]
        sphereBoxDerivatives = [1, 3, 2]

        sphere1 = SphereMesh(fm, coordinates, [0.0,0.0,0.0], axes, elementsCountAcross,
                             elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False)

        hemisphere = ShieldMesh3D(elementsCountAcross,
                                  elementsCountAcrossShell + elementsCountAcrossTransition - 1)  # 0

        # get hemisphere nodes from both cylinder end and top of the sphere and mix them
        hemisphere._boxDerivatives = sphere1._shield3D._boxDerivatives

        hemisphere._boxMapping = sphere1._shield3D._boxMapping
        hemisphere._box_deriv_mapping = sphere1._shield3D._box_deriv_mapping
        hemisphere._element_needs_scale_factor = sphere1._shield3D._element_needs_scale_factor
        hemisphere._xi_mapping = sphere1._shield3D._xi_mapping
        hemisphere._xi_signs = sphere1._shield3D._xi_signs

        hemisphere.px = sphere1._shield3D.px
        hemisphere.pd1 = sphere1._shield3D.pd1
        hemisphere.pd2 = sphere1._shield3D.pd2
        hemisphere.pd3 = sphere1._shield3D.pd3
        nodesIdCylinderProximalEnd = []
        for n3 in range(elementsCountAcross[2] + 1):
            for n2 in range(elementsCountAcross[0] + 1):
                for n1 in range(elementsCountAcross[1] + 1):
                    if n3 < elementsCountAcross[2] // 2:
                        if sphere1._shield3D.px[n3][n2][n1]:
                            hemisphere.px[n3][n2][n1] = vector.addVectors([sphere1._shield3D.px[n3][n2][n1],
                                                                            sphere_centre], [1, 1])
                    #cylinder end
                    elif n3 == elementsCountAcross[2] // 2:
                        # find nodes on the triple line. Note that cylinder and sphere have a little bit different
                        # numbering for nodes on the triple line
                        n2c, n1c = n2, n1
                        if n2 < 1 and n1 == n2:
                            n1c = 1
                        elif n2 < 1 and n1 == elementsCountAcross[1] - n2:
                            n1c = elementsCountAcross[1] - 1
                        elif n2 > elementsCountAcross[1] - 1:
                            if n1 == elementsCountAcross[1] - n2:
                                n1c = 1
                            elif n1 == n2:
                                n1c = elementsCountAcross[1] - 1
                        if cylinder1._shield.nodeId[0][n2c][n1c]:
                            nodesIdCylinderProximalEnd.append(cylinder1._shield.nodeId[0][n2c][n1c])

        # generate hemisphere extra nodes.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]],
                                   [0, int(0.5 * elementsCountAcross[2]) - 1]]
        radius_scale = options['Lower scale']
        for n2 in range(elementsCountAcross[0] + 1):
            for n1 in range(elementsCountAcross[1] + 1):
                if hemisphere.px[0][n2][n1]:
                    x = hemisphere.px[0][n2][n1]
                    hemisphere.px[0][n2][n1] = [radius_scale * x[0], radius_scale * x[1], x[2]]
        #
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = hemisphere.generateNodes(fm, coordinates, nodeIdentifier,
                                                  rangeOfRequiredElements)

        # replace the hemispherical nodes at interface with nodes on cylinder end
        count = 0
        hemisphereNodesToDelete = []
        for n2 in range(elementsCountAcross[0] + 1):
            for n1 in range(elementsCountAcross[1] + 1):
                if hemisphere.px[elementsCountAcross[2] // 2][n2][n1]:
                    hemisphereNodesToDelete.append(hemisphere.nodeId[elementsCountAcross[2] // 2][n2][n1])
                    hemisphere.nodeId[elementsCountAcross[2] // 2][n2][n1] = nodesIdCylinderProximalEnd[count]
                    count += 1


        # generate hemisphere elements.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]],
                                   [0, int(0.5 * elementsCountAcross[2])]]
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = hemisphere.generateElements(fm, coordinates, elementIdentifier,
                                                        rangeOfRequiredElements)

        # ***********************************************************************************************************************
        # top half of the hemisphere
        sphere_base = cylinder1._ellipses[-1]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = options['Upper scale_Z']
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                vector.setMagnitude(vector.crossproduct3(sphere_base.majorAxis, sphere_base.minorAxis),
                                    sphere_radius_3)]
        elementsCountAcross = [cylinder1._elementsCountAcrossMajor, cylinder1._elementsCountAcrossMinor,
                               cylinder1._elementsCountAcrossMajor]
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]],[0,-1]]

        sphereBoxDerivatives = [1, 3, 2]

        sphere2 = SphereMesh(fm, coordinates, [0.0, 0.0, 0.0], axes, elementsCountAcross,
                             elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False)

        hemisphere2 = ShieldMesh3D(elementsCountAcross,
                                   elementsCountAcrossShell + elementsCountAcrossTransition - 1)

        # get hemisphere nodes from both cylinder end and top of the sphere and mix them
        hemisphere2._boxDerivatives = sphere2._shield3D._boxDerivatives
        hemisphere2._boxMapping = sphere2._shield3D._boxMapping
        hemisphere2._box_deriv_mapping = sphere2._shield3D._box_deriv_mapping
        hemisphere2._element_needs_scale_factor = sphere2._shield3D._element_needs_scale_factor
        hemisphere2._xi_mapping = sphere2._shield3D._xi_mapping
        hemisphere2._xi_signs = sphere2._shield3D._xi_signs

        hemisphere2.px = sphere2._shield3D.px
        hemisphere2.pd1 = sphere2._shield3D.pd1
        hemisphere2.pd2 = sphere2._shield3D.pd2
        hemisphere2.pd3 = sphere2._shield3D.pd3
        nodesIdCylinderDistalEnd = []
        for n3 in range(elementsCountAcross[2] + 1):
            for n2 in range(elementsCountAcross[0] + 1):
                for n1 in range(elementsCountAcross[1] + 1):
                    if n3 > elementsCountAcross[2] // 2:
                        if sphere2._shield3D.px[n3][n2][n1]:
                            hemisphere2.px[n3][n2][n1] = vector.addVectors([sphere2._shield3D.px[n3][n2][n1],
                                                                            sphere_centre], [1, 1])
                    #cylinder end
                    elif n3 == elementsCountAcross[2] // 2:
                        # find nodes on the triple line. Note that cylinder and sphere have a little bit different
                        # numbering for nodes on the triple line
                        n2c, n1c = n2, n1
                        if n2 < 1 and n1 == n2:
                            n1c = 1
                        elif n2 < 1 and n1 == elementsCountAcross[1] - n2:
                            n1c = elementsCountAcross[1] - 1
                        elif n2 > elementsCountAcross[1] - 1:
                            if n1 == elementsCountAcross[1] - n2:
                                n1c = 1
                            elif n1 == n2:
                                n1c = elementsCountAcross[1] - 1
                        if cylinder1._shield.nodeId[-1][n2c][n1c]:
                            nodesIdCylinderDistalEnd.append(cylinder1._shield.nodeId[-1][n2c][n1c])

        # ******************************************************************************************************************************
        # generate hemisphere extra nodes.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]],
                                   [int(0.5 * elementsCountAcross[2])+1,elementsCountAcross[2]]]
        radius_scale = options['Upper scale']
        for n2 in range(elementsCountAcross[0] + 1):
            for n1 in range(elementsCountAcross[1] + 1):
                if hemisphere2.px[elementsCountAcross[2]][n2][n1]:
                    x = hemisphere2.px[elementsCountAcross[2]][n2][n1]
                    hemisphere2.px[elementsCountAcross[2]][n2][n1] = [radius_scale * x[0], radius_scale * x[1],
                                                                      x[2]]
        #
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = hemisphere2.generateNodes(fm, coordinates, nodeIdentifier,
                                                   rangeOfRequiredElements)

        # replace the hemispherical nodes at interface with nodes on cylinder end
        count = 0
        hemisphereNodesToDelete2 = []
        for n2 in range(elementsCountAcross[0] + 1):
            for n1 in range(elementsCountAcross[1] + 1):
                if hemisphere2.px[elementsCountAcross[2] // 2][n2][n1]:
                    hemisphereNodesToDelete2.append(hemisphere2.nodeId[elementsCountAcross[2] // 2][n2][n1])
                    hemisphere2.nodeId[elementsCountAcross[2] // 2][n2][n1] = nodesIdCylinderDistalEnd[count]
                    count += 1

        # generate hemisphere elements.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]],
                                   [int(0.5 * elementsCountAcross[2]),elementsCountAcross[2]]]
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = hemisphere2.generateElements(fm, coordinates, elementIdentifier,
                                                         rangeOfRequiredElements)





        cancellousGroup = AnnotationGroup(region, get_bone_term("Cancellous bone"))

        cancellousMeshGroup = cancellousGroup.getMeshGroup(mesh)


        corticalGroup = AnnotationGroup(region, get_bone_term("Cortical bone"))

        corticalMeshGroup = corticalGroup.getMeshGroup(mesh)
        annotationGroup = [cancellousGroup,corticalGroup]


        array2 = [1,2,7,14,15,22,27,28]
        array1 = [1,2,5,10,11,16,19,20]

        More_elements = array1 if elementsCountAcrossMajor == 6 else array2 if elementsCountAcrossMajor == 8 else None
        addition = elementsCountAcrossMajor-2 if elementsCountAcrossMajor == 6 else elementsCountAcrossMajor if elementsCountAcrossMajor == 8 else None



        number_of_elements_per_layer = max(More_elements) #(elementsCountAcrossMajor-1)*4

        outer_element_trunk = [x+number_of_elements_per_layer*i for i in range(0,elementsCountAlong) for x in More_elements]
        print(outer_element_trunk)
        max_value = max(outer_element_trunk)
        dome1=[max_value + i for i in range(1,5)]
        print(dome1)
        dome1_final = dome1 + [x+max(dome1)+addition for x in More_elements]
        print(dome1_final)
        max_value = max(dome1_final)

        dome2 = [x + max_value for x in More_elements]
        max_value = max(dome2)
        print(dome2)
        dome2_final = dome2 + [max_value + i+addition for i in range(1,5)]
        print(dome2_final)
        final_outer = outer_element_trunk+dome1_final+dome2_final


        final_inner = []

        for number in range (1,max(final_outer)-1):
            if number not in final_outer:
                final_inner.append(number)



        for key in final_inner:
            element = mesh.findElementByIdentifier(int(key))
            cancellousMeshGroup.addElement(element)

        for key in final_outer:
            element = mesh.findElementByIdentifier(int(key))
            corticalMeshGroup.addElement(element)

        top_elem = max(final_outer)
        bottom_elem = number_of_elements_per_layer*elementsCountAlong+1
        lateral_elem = int(number_of_elements_per_layer*(elementsCountAlong/2-1)+More_elements[3])
        medial_elem = int(number_of_elements_per_layer*(elementsCountAlong/2-1)+More_elements[4])


        # markers with element number and xi position
        allMarkersRadius = {"Greater tuberosity of radius": {"elementID": top_elem, "xi": [1.0, 0.0, 0.5]},
                      "Lateral epicondyle of radius": {"elementID": bottom_elem, "xi": [1.0, 0.0, 1.0]},
                      "Medial epicondyle of radius": {"elementID": medial_elem, "xi": [0.0, 0.0, 0.0]},
                      "Epiphysial plate of radius": {"elementID": lateral_elem, "xi": [0.0, 0.0, 0.0]},
        }

        allMarkersUlna = {"Greater tuberosity of radius": {"elementID": top_elem, "xi": [1.0, 0.0, 0.5]},
                      "Lateral epicondyle of radius": {"elementID": bottom_elem, "xi": [1.0, 0.0, 1.0]},
        }

        allMarkersTibia = {"Intercondylar eminence of tibia": {"elementID": top_elem, "xi": [0.0,1.0,0.0]},
                           "Inferior articular surface of tibia": {"elementID": bottom_elem, "xi": [1.0,0.0,1.0]},
                           "Medial surface of tibia": {"elementID": medial_elem, "xi": [1.0,1.0,1.0]},
                           "Lateral surface of tibia": {"elementID": lateral_elem, "xi":[1.0,1.0,1.0]},

        }

        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        markerPoints = markerGroup.getOrCreateNodesetGroup(nodes)
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)+10

        #annotation fiducial point
        if isRadius:
            for key in allMarkersRadius:
                xi = allMarkersRadius[key]["xi"]
                addMarker = {"name": key, "xi": allMarkersRadius[key]["xi"]}
                markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
                nodeIdentifier += 1
                cache.setNode(markerPoint)
                markerName.assignString(cache, addMarker["name"])
                elementID = allMarkersRadius[key]["elementID"]
                element = mesh.findElementByIdentifier(elementID)
                markerLocation.assignMeshLocation(cache, element, addMarker["xi"])

        elif isUlna:
            for key in allMarkersUlna:
                xi = allMarkersUlna[key]["xi"]
                addMarker = {"name": key, "xi": allMarkersUlna[key]["xi"]}
                markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
                nodeIdentifier += 1
                cache.setNode(markerPoint)
                markerName.assignString(cache, addMarker["name"])
                elementID = allMarkersRadius[key]["elementID"]
                element = mesh.findElementByIdentifier(elementID)
                markerLocation.assignMeshLocation(cache, element, addMarker["xi"])

        elif isTibia:
            for key in allMarkersTibia:
                xi = allMarkersTibia[key]["xi"]
                addMarker = {"name": key, "xi": allMarkersTibia[key]["xi"]}
                markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
                nodeIdentifier += 1
                cache.setNode(markerPoint)
                markerName.assignString(cache, addMarker["name"])
                elementID = allMarkersTibia[key]["elementID"]
                element = mesh.findElementByIdentifier(elementID)
                markerLocation.assignMeshLocation(cache, element, addMarker["xi"])


#         # AnnotationGroup.createMarkerNode(startNodeIdentifier=1, materialCoordinatesField: FieldFiniteElement = None, materialCoordinates = None, element = None, xi = [0.0, 0.0, 0.0])
#
#
#         # markerTermNameBoneCoordinatesMap = {
#         #     'tibial tuberosity': [-5.076472492200136e+01, -4.592226612078402e+01, -1.261953033704384e+03],
#         # }
#         # annotationGroups = []
#         # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
#         # nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
#         # print(nodeIdentifier)
#         # for termName, boneCoordinatesValues in markerTermNameBoneCoordinatesMap.items():
#         #     annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, ('tibia point', 'FM:63'))
#         #     annotationGroup.createMarkerNode(nodeIdentifier, coordinates, boneCoordinatesValues)
#         #     nodeIdentifier += 1

        # smoothing = DerivativeSmoothing(region, coordinates)
        # smoothing.smooth(True)

        #annotationGroup = []
        return annotationGroup, None

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

        return
