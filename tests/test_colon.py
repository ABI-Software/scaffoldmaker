import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colon1 import MeshType_3d_colon1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, exnodeStringFromNodeValues

from testutils import assertAlmostEqualList


class ColonScaffoldTestCase(unittest.TestCase):

    def test_colon1(self):
        """
        Test creation of colon scaffold.
        """
        scaffold = MeshType_3d_colon1
        parameterSetNames = MeshType_3d_colon1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cattle 1", "Human 1", "Human 2", "Human 3", "Mouse 1",
                                             "Mouse 2", "Pig 1"])

        options = scaffold.getDefaultOptions("Human 1")
        self.assertEqual(12, len(options))
        centralPath = options['Central path']
        segmentProfile = options.get("Segment profile")
        segmentSettings = segmentProfile.getScaffoldSettings()
        self.assertEqual(8, segmentSettings.get("Number of elements around haustrum"))
        self.assertEqual(0.5, segmentSettings.get("Corner inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Haustrum inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Segment length end derivative factor"))
        self.assertEqual(3, segmentSettings.get("Number of tenia coli"))
        self.assertEqual(0.6, segmentSettings.get("Tenia coli thickness"))
        self.assertEqual(30, options.get("Number of segments"))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(None, options.get("Transverse length"))
        self.assertEqual(None, options.get("Proximal inner radius"))
        self.assertEqual(10.0, options.get("Proximal-transverse tenia coli width"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE])[0]
        self.assertEqual(9, len(cx))
        assertAlmostEqualList(self, cx[0], [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, cx[1], [-47.4, 188.6, 0.0], 1.0E-6)
        del tmpRegion

        annotationGroups = scaffold.generateBaseMesh(region, options)
        self.assertEqual(11, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        if annotationGroups is not None:
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(15120, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(48726, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(52119, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(18513, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-102.92476320962298, -80.82196978241228, -65.70156604792459], 1.0E-6)
        assertAlmostEqualList(self, maximums, [519.9276475845722, 451.090297174708, 53.94971987636461], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [397.25007484305684, 1490.368038490398, 2.2], 1.0E-6)

        colonCoordinates = fieldmodule.findFieldByName("colon coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(colonCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.6, 0.0, -0.6], 1.0E-4)
        assertAlmostEqualList(self, maximums, [0.6, 24.0, 0.625], 1.0E-4)

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
        self.assertAlmostEqual(surfaceArea, 503394.01077174983, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 944619.4065131614, delta=1.0E-6)

    def test_mousecolon1(self):
        """
        Test creation of mouse colon scaffold.
        """
        options = MeshType_3d_colon1.getDefaultOptions("Mouse 2")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colon1.generateBaseMesh(region, options)
        self.assertEqual(10, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(1600, mesh3d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        colonCoordinates = fieldmodule.findFieldByName("colon coordinates").castFiniteElement()
        self.assertTrue(colonCoordinates.isValid())

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)
            colonSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, colonCoordinates, faceMeshGroup)
            colonSurfaceAreaField.setNumbersOfPoints(4)
            colonVolumeField = fieldmodule.createFieldMeshIntegral(one, colonCoordinates, mesh3d)
            colonVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 651.2624843023726, delta=1.0E-6)
        result, colonSurfaceArea = colonSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonSurfaceArea, 90.45789313556386, delta=1.0E-6)
        result, colonVolume = colonVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonVolume, 8.290058800222006, delta=1.0E-6)


if __name__ == "__main__":
    unittest.main()
