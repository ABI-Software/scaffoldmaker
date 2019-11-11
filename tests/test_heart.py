import unittest
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_3d_heart1 import MeshType_3d_heart1
from scaffoldmaker.utils.zinc_utils import evaluateFieldRange

def assertAlmostEqualList(testcase, actualList, expectedList, delta):
    assert len(actualList) == len(expectedList)
    for actual, expected in zip(actualList, expectedList):
        testcase.assertAlmostEqual(actual, expected, delta=delta)

class HeartScaffoldTestCase(unittest.TestCase):

    def test_heart1(self):
        """
        Test creation of heart scaffold.
        """
        parameterSetNames = MeshType_3d_heart1.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Pig 1", "Rat 1",
            "Unit Human 1", "Unit Mouse 1", "Unit Pig 1", "Unit Rat 1" ]);
        options = MeshType_3d_heart1.getDefaultOptions("Human 1")
        self.assertEqual(116, len(options))
        self.assertEqual(0.9, options.get("LV outer height"))
        self.assertEqual(80.0, options.get("Unit scale"))
        self.assertEqual(7, options.get("Number of elements around LV free wall"))
        self.assertEqual(7, options.get("Number of elements around RV free wall"))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_heart1.generateBaseMesh(region, options)

        self.assertEqual(18, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        if annotationGroups is not None:
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(289, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1117, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1356, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(523, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(7, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ -50.7876375290527, -57.76590573823474, -91.6 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 43.81084359764995, 39.03925080604259, 40.71693637558552 ], 1.0E-6)

if __name__ == "__main__":
    unittest.main()
