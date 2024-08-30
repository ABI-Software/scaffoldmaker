"""
Generates a solid bone using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division

import copy

from cmlibs.maths.vectorops import add_vectors, set_magnitude, cross
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters


class MeshType_3d_bone1 (Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    centralPathDefaultScaffoldPackages = {
        'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (2, [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (3, [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
                    (4, [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
                ]])
        })
    }

    @staticmethod
    def getName():
        return '3D Bone 1'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 6,
            'Number of elements across minor': 6,
            'Number of elements across shell': 1,
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

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fm.findMeshByDimension(3)

        sphere_shape = SphereShape.SPHERE_SHAPE_FULL
        sphere_base = cylinder1._ellipses[0]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = options['Lower scale_Z']
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                set_magnitude(cross(sphere_base.majorAxis, sphere_base.minorAxis),
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
                            hemisphere.px[n3][n2][n1] = add_vectors([sphere1._shield3D.px[n3][n2][n1], sphere_centre], [1, 1])
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

        # ******************************************************************************************************************************
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
                set_magnitude(cross(sphere_base.majorAxis, sphere_base.minorAxis),
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
                            hemisphere2.px[n3][n2][n1] = add_vectors([sphere2._shield3D.px[n3][n2][n1], sphere_centre], [1, 1])
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


        annotationGroup = []
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
