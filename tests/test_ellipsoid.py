from cmlibs.maths.vectorops import magnitude
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_ellipsoid1 import MeshType_3d_ellipsoid1
from testutils import assertAlmostEqualList
import unittest


class EllipsoidScaffoldTestCase(unittest.TestCase):

    def test_ellipsoid_2D(self):
        """
        Test creation of 2-D ellipsoid surface.
        """
        Scaffold = MeshType_3d_ellipsoid1
        parameterSetNames = Scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default"])
        options = Scaffold.getDefaultOptions("Default")
        self.assertEqual(13, len(options))
        self.assertEqual(4, options["Number of elements across axis 1"])
        self.assertEqual(6, options["Number of elements across axis 2"])
        self.assertEqual(8, options["Number of elements across axis 3"])
        self.assertFalse(options["2D surface only"])
        self.assertEqual(1, options["Number of transition elements"])
        self.assertEqual(1.0, options["Axis length x"])
        self.assertEqual(1.5, options["Axis length y"])
        self.assertEqual(2.0, options["Axis length z"])
        self.assertEqual(0.0, options["Axis 2 x-rotation degrees"])
        self.assertEqual(90.0, options["Axis 3 x-rotation degrees"])
        self.assertEqual(0.6, options["Advanced n-way derivative factor"])
        self.assertFalse(options["Refine"])
        self.assertEqual(4, options["Refine number of elements"])
        # set test options
        options["2D surface only"] = True

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = Scaffold.generateMesh(region, options)[0]
        self.assertEqual(0, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(0, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(88, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(176, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(90, nodes.getSize())

        # check coordinates range, sphere volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        TOL = 1.0E-6
        assertAlmostEqualList(self, minimums, [-1.0, -1.5, -2.0], TOL)
        assertAlmostEqualList(self, maximums, [1.0, 1.5, 2.0], TOL)
        # test symmetry of 3-way points
        fieldcache = fieldmodule.createFieldcache()
        node_3way1 = nodes.findNodeByIdentifier(90)
        fieldcache.setNode(node_3way1)
        result, x_3way1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, d1_3way1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, d2_3way1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        node_3way2 = nodes.findNodeByIdentifier(78)
        fieldcache.setNode(node_3way2)
        result, x_3way2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, d1_3way2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, d2_3way2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        assertAlmostEqualList(self, x_3way2, [x_3way1[0], -x_3way1[1], x_3way1[2]], TOL)
        assertAlmostEqualList(self, d1_3way2, [d1_3way1[0], -d1_3way1[1], d1_3way1[2]], TOL)
        assertAlmostEqualList(self, d2_3way2, [-d2_3way1[0], d2_3way1[1], -d2_3way1[2]], TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 27.86848567909992, delta=TOL)

if __name__ == "__main__":
    unittest.main()
