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

    def test_vagus_box1(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_vagus_box1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Human Trunk 1'])
        options = scaffold.getDefaultOptions()
        self.assertEqual(2, len(options))

        data_file = os.path.join(here, "resources", "vagus_data1.exf")
        print(data_file)

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        assert(region.isValid())
        data_region = region.getParent().createChild('data')
        assert(data_region.isValid())
        result = data_region.readFile(data_file)
        assert result == RESULT_OK

        vagus_data = VagusInputData(data_region)
        trunk_group_name = vagus_data.get_trunk_group_name()
        assert trunk_group_name != None


        #annotationGroups = scaffold.generateBaseMesh(region, options)[0]




if __name__ == "__main__":
    unittest.main()