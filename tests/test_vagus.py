import copy
import os
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.field import find_or_create_field_group
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field

from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.read_vagus_data import load_vagus_data, VagusInputData
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    get_nodeset_path_field_parameters

from scaffoldmaker.meshtypes.meshtype_3d_vagus_box1 import MeshType_3d_vagus_box1

from testutils import assertAlmostEqualList

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
        assert 'right level of superior border of the clavicle on the vagus nerve' in marker_data

        orientation_data = vagus_data.get_orientation_data()
        self.assertEqual(len(orientation_data), 1)
        self.assertEqual(len(orientation_data['orientation anterior']), 3)

        trunk_group_name = vagus_data.get_trunk_group_name()
        self.assertEqual(trunk_group_name, "right vagus X nerve trunk")
        trunk_coordinates = vagus_data.get_trunk_coordinates()
        self.assertEqual(len(trunk_coordinates), 32)
        annotation_term_map = vagus_data.get_annotation_term_map()
        assert trunk_group_name in annotation_term_map
        self.assertEqual(annotation_term_map[trunk_group_name], 'http://purl.obolibrary.org/obo/UBERON_0035020')

        branch_data = vagus_data.get_branch_data()
        self.assertEqual(len(branch_data), 3)
        assert 'right B branch' in branch_data
        self.assertEqual(len(branch_data['right A branch']), 5)
        assert 'right A branch of branch B' in branch_data
        self.assertEqual(len(branch_data['right A branch of branch B']), 2)

        branch_parents = vagus_data.get_branch_parent_map()
        self.assertEqual(branch_parents['right B branch'], trunk_group_name)
        self.assertEqual(branch_parents['right A branch of branch B'], 'right B branch')

        branch_common_groups = vagus_data.get_branch_common_group_map()
        self.assertEqual(len(branch_common_groups), 1)
        assert 'right branches' in branch_common_groups
        assert 'right B branch' in branch_common_groups['right branches']

    def test_vagus_box(self):
        """
        Test creation of vagus scaffold.
        """

        scaffold = MeshType_3d_vagus_box1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human Left Vagus 1', 'Human Right Vagus 1'])
        options = scaffold.getDefaultOptions("Human Right Vagus 1")
        self.assertEqual(len(options), 2)
        self.assertEqual(options.get('Number of points along the trunk'), 30)

        data_file = os.path.join(here, "resources", "vagus_data1.exf")
        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        data_region = region.getParent().createChild('data')
        result = data_region.readFile(data_file)

        annotation_groups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(len(annotation_groups), 13)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(35, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(319, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(627, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(40, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        marker_group = fieldmodule.findFieldByName("marker").castGroup()
        marker_nodes = marker_group.getNodesetGroup(nodes)
        self.assertEqual(4, marker_nodes.getSize())
        marker_name = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(marker_name.isValid())
        marker_location = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(marker_location.isValid())

        vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()
        trunk_group_name = "right vagus X nerve trunk"
        trunk_vagus_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        minimums, maximums = evaluateFieldNodesetRange(trunk_vagus_group, nodes)
        self.assertEqual(minimums, 0.0)
        self.assertEqual(maximums, 1.0)


if __name__ == "__main__":
    unittest.main()
