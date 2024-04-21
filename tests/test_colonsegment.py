import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class ColonSegmentScaffoldTestCase(unittest.TestCase):

    def test_humancolonsegment1(self):
        """
        Test creation of human colon segment scaffold.
        """
        parameterSetNames = MeshType_3d_colonsegment1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cattle 1", "Human 1", "Human 2", "Mouse 1", "Pig 1"])
        options = MeshType_3d_colonsegment1.getDefaultOptions("Human 1")
        self.assertEqual(32, len(options))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(2, options.get("Number of elements around tenia coli"))
        self.assertEqual(4, options.get("Number of elements along segment"))
        self.assertEqual(4, options.get("Number of elements through wall"))
        self.assertEqual(None, options.get("Start inner radius"))
        self.assertEqual(None, options.get("Start inner radius derivative"))
        self.assertEqual(None, options.get("End inner radius"))
        self.assertEqual(None, options.get("End inner radius derivative"))
        self.assertEqual(3.0, options.get("Segment length mid derivative factor"))
        self.assertEqual(50.0, options.get("Segment length"))
        self.assertEqual(3, options.get("Number of tenia coli"))
        self.assertEqual(11.1, options.get("Start tenia coli width"))
        self.assertEqual(0.4, options.get("End tenia coli width derivative"))
        self.assertEqual(1.6, options.get("Wall thickness"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colonsegment1.generateBaseMesh(region, options)[0]
        self.assertEqual(8, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(504, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1746, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(2007, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(765, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-2.172286248499807e-15, -59.07417787753859, -55.60850335002243], 1.0E-6)
        assertAlmostEqualList(self, maximums, [50.0, 50.60638395601446, 55.60850335002243], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [392.25861666443416, 50.0, 2.2], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 20870.63159749863, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 39988.68541401406, delta=1.0E-6)

    def test_mousecolonsegment1(self):
        """
        Test creation of mouse colon segment scaffold.
        """
        options = MeshType_3d_colonsegment1.getDefaultOptions("Mouse 1")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colonsegment1.generateBaseMesh(region, options)[0]
        self.assertEqual(7, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(160, mesh3d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)

        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 467.98576432245915, delta=1.0E-6)
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 467.98905126160287, delta=1.0E-3)


if __name__ == "__main__":
    unittest.main()
