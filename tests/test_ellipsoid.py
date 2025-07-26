from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
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
        self.assertEqual(12, len(options))
        self.assertEqual(4, options.get("Number of elements across axis 1"))
        self.assertEqual(6, options.get("Number of elements across axis 2"))
        self.assertEqual(8, options.get("Number of elements across axis 3"))
        self.assertFalse(options.get("2D surface only"))
        self.assertEqual(1, options.get("Number of transition elements"))
        self.assertEqual(1.0, options.get("Axis length x"))
        self.assertEqual(1.5, options.get("Axis length y"))
        self.assertEqual(2.0, options.get("Axis length z"))
        self.assertEqual(0.0, options.get("Axis 2 x-rotation degrees"))
        self.assertEqual(90.0, options.get("Axis 3 x-rotation degrees"))
        self.assertFalse(options.get("Refine"))
        self.assertEqual(4, options.get("Refine number of elements"))
        # set test options
        options["2D surface only"] = True
        options["Axis 2 x-rotation degrees"] = -45.0
        options["Axis 3 x-rotation degrees"] = 45.0

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
        assertAlmostEqualList(self, minimums, [-1.0, -1.4882638616539217, -1.9598574478664705], 1.0E-6)
        assertAlmostEqualList(self, maximums, [1.0, 1.4882638616539217, 1.9598574478664705], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 27.843226402805207, delta=1.0E-4)


if __name__ == "__main__":
    unittest.main()
