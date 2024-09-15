import copy
import os
import unittest

from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.meshtype_3d_vagus1 import MeshType_3d_vagus1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    get_nodeset_path_field_parameters

from testutils import assertAlmostEqualList

here = os.path.abspath(os.path.dirname(__file__))

class VagusScaffoldTestCase(unittest.TestCase):

    def test_vagus1(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_vagus1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Human Left Trunk 1'])
        options = scaffold.getDefaultOptions()
        self.assertEqual(1, len(options))

        data_file = os.path.join(here, "resources", "vagus_contoursWithMarkers.exf")
        print(data_file)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        data_region = region.getParent().createChild("data")
        print(region.getParent().getNextSibling().getName())


        assert data_region.isValid()


        result = data_region.readFile(data_file)
        #assert result == RESULT_OK
        print(result)




if __name__ == "__main__":
    unittest.main()