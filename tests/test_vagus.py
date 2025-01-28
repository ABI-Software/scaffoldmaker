import copy
import os
import unittest

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

        marker_data = vagus_data.get_level_markers()
        self.assertEqual(len(marker_data), 2)
        orientation_data = vagus_data.get_orientation_data()
        self.assertEqual(len(orientation_data), 1)
        self.assertEqual(len(orientation_data['orientation anterior']), 3)

        trunk_group_name = vagus_data.get_trunk_group_name()
        self.assertEqual(trunk_group_name, "left vagus X nerve trunk")
        branch_data = vagus_data.get_branch_data()
        self.assertEqual(len(branch_data), 3)
        annotation_term_map = vagus_data.get_annotation_term_map()
        assert trunk_group_name in annotation_term_map

        branch_common_groups = vagus_data.get_branch_common_group_map()
        print(branch_common_groups)


    def test_vagus_box(self):
        """
        Test creation of vagus scaffold.
        """

        scaffold = MeshType_3d_vagus_box1
        parameterSetNames = scaffold.getParameterSetNames()
        options = scaffold.getDefaultOptions()
        self.assertEqual(2, len(options))

        data_file = os.path.join(here, "resources", "vagus_data1.exf")
        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        assert (region.isValid())
        data_region = region.getParent().createChild('data')
        assert (data_region.isValid())
        result = data_region.readFile(data_file)
        assert result == RESULT_OK

        annotationGroups = scaffold.generateBaseMesh(region, options)[0]


if __name__ == "__main__":
    unittest.main()
