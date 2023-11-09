import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_2d_tubenetwork1 import MeshType_2d_tubenetwork1
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters

from testutils import assertAlmostEqualList


class NetworkScaffoldTestCase(unittest.TestCase):

    def test_network_layout(self):
        """
        Test creation of network layout scaffold.
        """
        scaffold = MeshType_1d_network_layout1
        options = scaffold.getDefaultOptions()
        self.assertEqual(2, len(options))
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

    def test_2d_tube_network_sphere_cube(self):
        """
        Test sphere cube is generated correctly.
        """
        scaffold = MeshType_2d_tubenetwork1
        options = scaffold.getDefaultOptions("Sphere cube")
        self.assertEqual(4, len(options))
        options["Elements count around"] = 8
        options["Target element aspect ratio"] = 2.0
        options["Serendipity"] = True

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            annotationGroups, _ = scaffold.generateBaseMesh(region, options)
        self.assertEqual(0, len(annotationGroups))

        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(32 * 12, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(8 * 7 * 12 + 4 * 3 * 8, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(8 * 3 * 12 + (2 + 3 * 3) * 8, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5658330029836134, -0.5787347924219615, -0.5903131690984363], 1.0E-8)
        assertAlmostEqualList(self, maximums, [0.5658331365887045, 0.578734792421953, 0.5903125661076012], 1.0E-8)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, 3.879343151083442, delta=1.0E-8)


if __name__ == "__main__":
    unittest.main()
