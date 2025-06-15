import os
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.field import find_or_create_field_group
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_nerve1 import MeshType_3d_nerve1
from scaffoldmaker.utils.read_vagus_data import VagusInputData

here = os.path.abspath(os.path.dirname(__file__))


class VagusScaffoldTestCase(unittest.TestCase):

    def test_input_vagus_data(self):
        """
        Reading vagus input data
        """

        data_file = os.path.join(here, "resources", "vagus_data1.exf")

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        assert(region.isValid())
        data_region = region.getParent().createChild('data')
        assert(data_region.isValid())
        result = data_region.readFile(data_file)
        assert result == RESULT_OK

        vagus_data = VagusInputData(data_region)
        self.assertEqual(vagus_data.get_side_label(), 'right')

        marker_data = vagus_data.get_level_markers()
        self.assertEqual(len(marker_data), 2)
        self.assertTrue('right level of superior border of the clavicle on the vagus nerve' in marker_data)

        orientation_data = vagus_data.get_orientation_data()
        self.assertEqual(len(orientation_data), 1)
        self.assertEqual(len(orientation_data['orientation anterior']), 3)

        trunk_group_name = vagus_data.get_trunk_group_name()
        self.assertEqual(trunk_group_name, "right vagus X nerve trunk")
        trunk_coordinates = vagus_data.get_trunk_coordinates()
        self.assertEqual(len(trunk_coordinates), 32)
        annotation_term_map = vagus_data.get_annotation_term_map()
        self.assertTrue(trunk_group_name in annotation_term_map)
        self.assertEqual(annotation_term_map[trunk_group_name], 'http://purl.obolibrary.org/obo/UBERON_0035020')

        branch_data = vagus_data.get_branch_data()
        self.assertEqual(len(branch_data), 3)
        self.assertTrue('right B branch' in branch_data)
        self.assertEqual(len(branch_data['right A branch']), 5)
        self.assertTrue('right A branch of branch B' in branch_data)
        self.assertEqual(len(branch_data['right A branch of branch B']), 2)

        branch_parents = vagus_data.get_branch_parent_map()
        self.assertEqual(branch_parents['right B branch'], trunk_group_name)
        self.assertEqual(branch_parents['right A branch of branch B'], 'right B branch')

        branch_common_groups = vagus_data.get_branch_common_group_map()
        self.assertEqual(len(branch_common_groups), 1)
        self.assertTrue('right branches' in branch_common_groups)
        self.assertTrue('right B branch' in branch_common_groups['right branches'])

    def test_no_input_file(self):
        """
        No input file.
        """

        scaffold = MeshType_3d_nerve1
        options = scaffold.getDefaultOptions("Human Right Vagus 1")

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')

        # check that it doesn't crash with  no input file
        annotation_groups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(annotation_groups, [])

    def test_vagus_box(self):
        """
        Test creation of vagus scaffold.
        """
        scaffold = MeshType_3d_nerve1
        scaffoldname = MeshType_3d_nerve1.getName()
        self.assertEqual(scaffoldname, '3D Nerve 1')
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human Left Vagus 1', 'Human Right Vagus 1'])
        options = scaffold.getDefaultOptions("Human Right Vagus 1")
        self.assertEqual(len(options), 2)
        self.assertEqual(options.get('Number of elements along the trunk'), 30)
        options['Number of elements along the trunk'] = 8

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        data_region = region.getParent().createChild('data')
        data_file = os.path.join(here, "resources", "vagus_data1.exf")
        result = data_region.readFile(data_file)

        # check annotation groups
        annotation_groups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(len(annotation_groups), 13)

        # check all meshes are created
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(14, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(130, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(270, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(19, nodes.getSize())  # including 4 marker points
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check all markers are added
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        marker_group = fieldmodule.findFieldByName("marker").castGroup()
        marker_nodes = marker_group.getNodesetGroup(nodes)
        self.assertEqual(4, marker_nodes.getSize())
        marker_name = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(marker_name.isValid())
        marker_location = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(marker_location.isValid())

        # check vagus material coordinates
        vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()
        trunk_group_name = "right vagus X nerve trunk"
        trunk_vagus_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        minimums, maximums = evaluateFieldNodesetRange(trunk_vagus_group, nodes)
        self.assertEqual(minimums, 0.0)
        self.assertEqual(maximums, 1.0)

        # check annotation group: added and connected to the right node segment
        trunk_annotation_group = findAnnotationGroupByName(annotation_groups, trunk_group_name)
        trunk_nodeset_group = trunk_annotation_group.getNodesetGroup(nodes)
        self.assertEqual(trunk_nodeset_group.getSize(), 9)

        branchB_annotation_group = findAnnotationGroupByName(annotation_groups, "right B branch")
        self.assertNotEqual(branchB_annotation_group, 'None')
        branchB_nodeset_group = branchB_annotation_group.getNodesetGroup(nodes)
        self.assertEqual(branchB_nodeset_group.getSize(), 4)
        branchB_mesh_group = branchB_annotation_group.getMeshGroup(mesh3d)
        self.assertEqual(branchB_mesh_group.getSize(), 2)

        node_iter = branchB_nodeset_group.createNodeiterator()
        node = node_iter.next()
        self.assertEqual(node.getIdentifier(), 6)  # first trunk node


if __name__ == "__main__":
    unittest.main()
