import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine

from testutils import assertAlmostEqualList


class NetworkScaffoldTestCase(unittest.TestCase):

    def test_network_layout(self):
        """
        Test creation of network layout scaffold.
        """
        scaffold = MeshType_1d_network_layout1
        options = scaffold.getDefaultOptions()
        self.assertEqual(1, len(options))
        self.assertEqual("1-2", options.get("Structure"))
        options["Structure"] = "1-2-3,3-4,3.2-5"

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups, networkMesh = scaffold.generateBaseMesh(region, options)
        self.assertEqual(0, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(4, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(5, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        interactiveFunctions = scaffold.getInteractiveFunctions()
        functionOptions = None
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Smooth derivatives...":
                functionOptions = interactiveFunction[1]
                break
        functionOptions["Update directions"] = True
        scaffold.smoothDerivatives(region, options, None, functionOptions, "meshEdits")

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.25, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [3.0, 0.25, 0.0], 1.0E-6)

        networkSegments = networkMesh.getNetworkSegments()
        self.assertEqual(3, len(networkSegments))
        self.assertEqual([1, 2, 3], networkSegments[0].getNodeIdentifiers())
        self.assertEqual([1, 1, 1], networkSegments[0].getNodeVersions())
        self.assertEqual([3, 4], networkSegments[1].getNodeIdentifiers())
        self.assertEqual([1, 1], networkSegments[1].getNodeVersions())
        self.assertEqual([3, 5], networkSegments[2].getNodeIdentifiers())
        self.assertEqual([2, 1], networkSegments[2].getNodeVersions())

        # get path parameters with versions
        nx, nd1 = get_nodeset_path_ordered_field_parameters(
            nodes, coordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1],
            networkSegments[2].getNodeIdentifiers(), networkSegments[2].getNodeVersions())
        self.assertEqual(2, len(nx))
        assertAlmostEqualList(self, nx[0], [2.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, nx[1], [3.0, 0.25, 0.0], 1.0E-6)
        expected_nd = [nx[1][c] - nx[0][c] for c in range(3)]
        assertAlmostEqualList(self, nd1[0], expected_nd, 1.0E-6)
        assertAlmostEqualList(self, nd1[1], expected_nd, 1.0E-6)

if __name__ == "__main__":
    unittest.main()
