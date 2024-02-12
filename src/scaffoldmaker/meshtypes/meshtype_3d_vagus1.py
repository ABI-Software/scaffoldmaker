"""
Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
"""

import copy

from cmlibs.maths.vectorops import add, sub
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, get_group_list
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldfitter.fitter import Fitter


class MeshType_3d_vagus1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @staticmethod
    def getName():
        return "3D Vagus 1"

    @staticmethod
    def getParameterSetNames():
        return [
            'Human Left Trunk 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk': 20,
            'Number of elements along the branch': 3
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along the trunk',
            'Number of elements along the branch',
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options['Number of elements along the trunk'] < 10:
            options['Number of elements along the trunk'] = 10

        if options['Number of elements along the branch'] < 3:
            options['Number of elements along the branch'] = 3
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        marker_TermNameVagusLengthList = {
            "centroid of level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
            "centroid of level of superior border of the jugular foramen on the vagal trunk": 8.6342,
            "centroid of level of inferior border of the jugular foramen on the vagal trunk": 16.7227,
            "centroid of level of C1 transverse process on the vagal trunk": 32.1129,
            "centroid of level of angle of mandible on the vagal trunk": 42.2450,
            "centroid of level of tubercles of the greater horn of hyoid bone on the vagal trunk": 45.6122,
            "centroid of level of carotid bifurcation on the vagal trunk": 48.3581,
            "centroid of level of laryngeal prominence on the vagal trunk": 68.8431,
            "centroid of level of superior border of clavicle on the vagal trunk": 117.5627,
            "centroid of level of jugular notch on the vagal trunk": 124.6407,
            "centroid of level of sternal angle on the vagal trunk": 151.2352,
            "centroid of 1 cm superior to esophageal plexus on the vagal trunk": 165.5876,
            "centroid of level of esophageal hiatus on the vagal trunk": 254.32879,
            "centroid of level of aortic hiatus on the vagal trunk": 291.3695,
            "centroid of level of end of trunk on Japanese dataset": 312.5 # note this term is also not on the list of annotations
        }

        lengthToDiameterRatio = 312.5 # calculated from total length of nerve/average diameter of nerve
        rescaledmarker_TermNameVagusLengthList = {}
        for term in marker_TermNameVagusLengthList:
            rescaledmarker_TermNameVagusLengthList[term] = marker_TermNameVagusLengthList[term] / 312.5 * lengthToDiameterRatio

        # extracting the data from datafile
        data_region = region.getParent().findChildByName('data')
        assert data_region.isValid()
        fieldmodule = data_region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)

        markers = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        marker_coordinates = fieldmodule.findFieldByName("marker_data_coordinates").castFiniteElement()
        marker_names = fieldmodule.findFieldByName("marker_data_name")
        assert marker_coordinates.isValid() and (marker_coordinates.getNumberOfComponents() == 3)

        group_list = get_group_list(fieldmodule)
        group_map = {}
        branch_groups = []
        for group in group_list:
            group_name = group.getName()
            group_map[group_name] = group
            if 'trunk' in group_name:
                trunk_group = group_name
            if 'branch' and 'left A' in group_name:
                branch_groups.append(group_name)

        marker_data = {}
        marker_group = group_map.get("marker")
        if marker_group:
            marker_nodes = marker_group.getNodesetGroup(markers)
            marker_node_iter = marker_nodes.createNodeiterator()
            marker_node = marker_node_iter.next()
            while marker_node.isValid():
                fieldcache.setNode(marker_node)
                result, x = marker_coordinates.evaluateReal(fieldcache, 3)
                marker_name = marker_names.evaluateString(fieldcache)
                marker_data[marker_name] = x
                marker_node = marker_node_iter.next()
            #print(marker_data.keys())

        left_vagus_trunk_group = group_map.get(trunk_group)
        if left_vagus_trunk_group:
            trunk_data_coordinates_list = []
            trunk_nodes = left_vagus_trunk_group.getNodesetGroup(nodes)
            node_iter = trunk_nodes.createNodeiterator()
            node = node_iter.next()
            while node.isValid():
                fieldcache.setNode(node)
                result, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                trunk_data_coordinates_list.append(x)
                node = node_iter.next()
            print(trunk_group, len(trunk_data_coordinates_list))

        branch_data = {}
        for branch_name in branch_groups:
            branch_group = group_map.get(branch_name)
            if branch_group:
                branch_data_coordinates_list = []
                branch_nodes = branch_group.getNodesetGroup(nodes)
                node_iter = branch_nodes.createNodeiterator()
                node = node_iter.next()
                while node.isValid():
                    fieldcache.setNode(node)
                    result, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                    branch_data_coordinates_list.append(x)
                    node = node_iter.next()
                branch_data[branch_name] = branch_data_coordinates_list
                print(branch_name, len(branch_data_coordinates_list))


        use_marker_names = [
            "centroid of level of angle of mandible on the vagal trunk",
            "centroid of level of aortic hiatus on the vagal trunk"
        ]
        assert [name in rescaledmarker_TermNameVagusLengthList for name in use_marker_names]
        x1, y1, z1 = marker_data[use_marker_names[0]]
        x2, y2, z2 = marker_data[use_marker_names[1]]
        t1 = rescaledmarker_TermNameVagusLengthList[use_marker_names[0]]
        t2 = rescaledmarker_TermNameVagusLengthList[use_marker_names[1]]
        dx, dy, dz = [(x2 - x1) / (t2 - t1), (y2 - y1) / (t2 - t1), (z2 - z1) / (t2 - t1)]

        elementsAlongTrunk = options['Number of elements along the trunk']
        step = lengthToDiameterRatio/(elementsAlongTrunk - 1)
        trunk_nodes = []
        vagus_trunk_nodes = []
        for i in range(elementsAlongTrunk):
            trunk_nodes.append([
                x1 + dx * (i * step - t1),
                y1 + dy * (i * step - t1),
                z1 + dz * (i * step - t1)])
            vagus_trunk_nodes.append([0, 0, i*step])

        elementsAlongBranch = options['Number of elements along the branch']
        branch_nodes = {}
        vagus_branch_nodes = {}
        for branch_name in branch_data.keys():
            # find closest points in the branch and in trunk nodes to determine branch approximate start
            min_distance_squared = float('inf')
            closest_pair = (None, None)
            for ind, trunk_node in enumerate(trunk_nodes):
                for branch_node in branch_data[branch_name]:
                    distance_squared = (trunk_node[0] - branch_node[0])**2 + \
                                       (trunk_node[1] - branch_node[1])**2 + \
                                       (trunk_node[2] - branch_node[2])**2
                    if distance_squared < min_distance_squared:
                        min_distance_squared = distance_squared
                        closest_pair = (trunk_node, branch_node)
                        indexTrunkNode = ind

            branch_start_name = branch_name
            branch_start = trunk_nodes[indexTrunkNode]
            vagus_branch_start = vagus_trunk_nodes[indexTrunkNode]

            # for now while branch direction is unknown
            #branch_end = [closest_pair[1][0] + 15, closest_pair[1][1] - 15, closest_pair[1][2] - 15]
            branch_data_starts_from_end = (branch_data[branch_name][-1][0] - branch_start[0]) ** 2 + \
                                          (branch_data[branch_name][-1][1] - branch_start[1]) ** 2 + \
                                          (branch_data[branch_name][-1][2] - branch_start[2]) ** 2 \
                                          < \
                                          (branch_data[branch_name][0][0] - branch_start[0]) ** 2 + \
                                          (branch_data[branch_name][0][1] - branch_start[1]) ** 2 + \
                                          (branch_data[branch_name][0][2] - branch_start[2]) ** 2
            branch_end = branch_data[branch_name][0] if branch_data_starts_from_end else branch_data[branch_name][-1]

            branch_coordinates = []
            vagus_branch_coordinates = []
            b_dx, b_dy, b_dz = [branch_end[0] - branch_start[0], branch_end[1] - branch_start[1], branch_end[2] - branch_start[2]]
            for i in range(elementsAlongBranch):
                branch_coordinates.append([
                    branch_start[0] + b_dx * i / (elementsAlongBranch - 1),
                    branch_start[1] + b_dy * i / (elementsAlongBranch - 1),
                    branch_start[2] + b_dz * i / (elementsAlongBranch - 1),
                ])
                vagus_branch_coordinates.append([
                    vagus_branch_start[0] + b_dx * i / (elementsAlongBranch - 1),
                    vagus_branch_start[1] + b_dy * i / (elementsAlongBranch - 1),
                    vagus_branch_start[2] + b_dz * i / (elementsAlongBranch - 1),
                ])
            branch_nodes[branch_start_name] = branch_coordinates
            vagus_branch_nodes[branch_start_name] = vagus_branch_coordinates



        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        annotationGroups = []
        vagusTrunkGroup = AnnotationGroup(region, (trunk_group, 'None'))
        annotationGroups.append(vagusTrunkGroup)
        branchGroups = {}
        for branch_name in branch_data.keys():
            branchGroup = AnnotationGroup(region, (branch_name, 'None'))
            branchGroups[branch_name] = branchGroup
            annotationGroups.append(branchGroup)

        # Geometric coordinates
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)

        mesh = fieldmodule.findMeshByDimension(1)
        elementbasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(elementbasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        elementtemplate.defineField(coordinates, -1, eft)

        nodeIdentifier = 1
        elementIdentifier = 1

        ld1 = [1.0, 0.0, 0.0]
        ld2 = [0.0, 1.0, 0.0]
        ld3 = [0.0, 0.0, 1.0]

        # build trunk
        for n in range(elementsAlongTrunk):
            lx = trunk_nodes[n]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, ld2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, ld3)

            if n > 0:
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [nodeIdentifier - 1, nodeIdentifier])
                meshGroup = vagusTrunkGroup.getMeshGroup(mesh)
                meshGroup.addElement(element)
                elementIdentifier += 1

            nodeIdentifier += 1

        # build branches
        for branch_name in branch_nodes.keys():
            for n in range(elementsAlongBranch):
                lx = branch_nodes[branch_name][n]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, ld2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, ld3)

                if n > 0:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, [nodeIdentifier - 1, nodeIdentifier])
                    meshGroup = branchGroups[branch_name].getMeshGroup(mesh)
                    meshGroup.addElement(element)
                    elementIdentifier += 1

                nodeIdentifier += 1

        # set marker nodes
        for marker_name, marker_coordinate in marker_data.items():
            if 'centroid of' in marker_name:
                annotationGroup = findOrCreateAnnotationGroupForTerm(
                    annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
                annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
                nodeIdentifier += 1




        # Material coordinates
        vagusCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="vagus coordinates")
        vagusNodetemplate = nodes.createNodetemplate()
        vagusNodetemplate.defineField(vagusCoordinates)
        vagusNodetemplate.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        vagusNodetemplate.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        vagusNodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        vagusNodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        vagusNodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        vagusNodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)

        eftVagus = mesh.createElementfieldtemplate(elementbasis)
        vagusElementtemplate = mesh.createElementtemplate()
        vagusElementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        vagusElementtemplate.defineField(vagusCoordinates, -1, eftVagus)

        vagusNodeIdentifier = 1
        vagusElementIdentifier = 1

        # build trunk in material coordinates
        for n in range(elementsAlongTrunk):
            lx = vagus_trunk_nodes[n]
            node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
            node.merge(vagusNodetemplate)
            fieldcache.setNode(node)
            vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
            vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
            vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, ld2)
            vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, ld3)

            if n > 0:
                element = mesh.findElementByIdentifier(vagusElementIdentifier)
                element.merge(vagusElementtemplate)
                element.setNodesByIdentifier(eftVagus, [vagusNodeIdentifier - 1, vagusNodeIdentifier])
                vagusElementIdentifier += 1

            vagusNodeIdentifier += 1

        # build branches in material coordinates
        for branch_name in vagus_branch_nodes.keys():
            for n in range(elementsAlongBranch):
                lx = vagus_branch_nodes[branch_name][n]
                node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
                node.merge(vagusNodetemplate)
                fieldcache.setNode(node)
                vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
                vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, ld2)
                vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, ld3)

                if n > 0:
                    element = mesh.findElementByIdentifier(vagusElementIdentifier)
                    element.merge(vagusElementtemplate)
                    element.setNodesByIdentifier(eft, [vagusNodeIdentifier - 1, vagusNodeIdentifier])
                    vagusElementIdentifier += 1

                vagusNodeIdentifier += 1








        return annotationGroups, None
