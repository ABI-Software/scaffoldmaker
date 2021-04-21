import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_stomach1 import MeshType_3d_stomach1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace
from testutils import assertAlmostEqualList

class StomachScaffoldTestCase(unittest.TestCase):

    def test_stomach1(self):
        """
        Test creation of stomach scaffold.
        """
        parameterSetNames = MeshType_3d_stomach1.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Rat 1" ])
        options = MeshType_3d_stomach1.getDefaultOptions("Rat 1")
        self.assertEqual(17, len(options))
        self.assertEqual(12, options.get("Number of elements around esophagus"))
        self.assertEqual(14, options.get("Number of elements around duodenum"))
        self.assertEqual(2, options.get("Number of elements between annulus and duodenum"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(1, options.get("Number of radial elements in annulus"))
        self.assertEqual(0.5, options.get("Wall thickness"))
        self.assertEqual(True, options.get("Limiting ridge"))
        ostiumOptions = options['Gastro-esophagal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        self.assertEqual(1, ostiumSettings.get("Number of vessels"))
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))
        self.assertEqual(4.0, ostiumSettings.get("Ostium diameter"))
        self.assertEqual(3.5, ostiumSettings.get("Ostium length"))
        self.assertEqual(0.5, ostiumSettings.get("Ostium wall thickness"))
        self.assertEqual(1.25, ostiumSettings.get("Vessel inner diameter"))
        self.assertEqual(0.5, ostiumSettings.get("Vessel wall thickness"))
        self.assertEqual(0.0, ostiumSettings.get("Vessel angle 1 degrees"))
        self.assertEqual(0.55, options.get("Gastro-esophagal junction position along factor"))
        self.assertEqual(0.2, options.get("Annulus derivative factor"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_stomach1.generateBaseMesh(region, options)
        self.assertEqual(8, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(158, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(643, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(823, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(341, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-17.977095466754566, -15.437578891036395, -8.706694306455596], 1.0E-6)
        assertAlmostEqualList(self, maximums, [17.943222678234513, 15.191836205539767, 8.725549026319936], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 2436.4955183926895, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 1168.9418891341106, delta=1.0E-6)

if __name__ == "__main__":
    unittest.main()
