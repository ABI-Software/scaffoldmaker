import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_solidcylinder1 import MeshType_3d_solidcylinder1

from testutils import assertAlmostEqualList


class CylinderScaffoldTestCase(unittest.TestCase):

    def test_cylinder1(self):
        """
        Test creation of cylinder scaffold.
        """
        scaffold = MeshType_3d_solidcylinder1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default"])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(14, len(options))
        self.assertEqual(4, options.get("Number of elements across major"))
        self.assertEqual(4, options.get("Number of elements across minor"))
        self.assertEqual(0, options.get("Number of elements across shell"))
        self.assertEqual(1, options.get("Number of elements across transition"))
        self.assertEqual(1, options.get("Number of elements along"))
        self.assertEqual(1.0, options.get("Shell element thickness proportion"))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)
        self.assertEqual(0, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(12, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(52, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(73, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(34, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, cylinder volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-1.0, -1.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 1.0, 1.0, 3.0 ], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 9.414866630615249, delta=1.0E-3)


if __name__ == "__main__":
    unittest.main()
