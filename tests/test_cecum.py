import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.cecum_terms import get_cecum_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_3d_cecum1 import MeshType_3d_cecum1
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class CecumScaffoldTestCase(unittest.TestCase):

    def test_cecum1(self):
        """
        Test creation of cecum scaffold.
        """
        parameterSetNames = MeshType_3d_cecum1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Human 2", "Human 3", "Pig 1"])
        options = MeshType_3d_cecum1.getDefaultOptions("Human 3")

        centralPath = options.get("Central path")
        centralPathSettings = centralPath.getScaffoldSettings()
        self.assertEqual("1-2-3.2, 4-3-5", centralPathSettings["Structure"])

        self.assertEqual(28, len(options))
        self.assertEqual(1, options.get("Number of segments"))
        self.assertEqual(2, options.get("Number of elements around tenia coli"))
        self.assertEqual(8, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(0.5, options.get("Corner inner radius factor"))
        self.assertEqual(0.4, options.get("Haustrum inner radius factor"))
        self.assertEqual(3.0, options.get("Segment length mid derivative factor"))
        self.assertEqual(3, options.get("Number of tenia coli"))
        self.assertEqual(10.0, options.get("Start tenia coli width"))
        self.assertEqual(0.0, options.get("End tenia coli width derivative"))
        self.assertEqual(1.6, options.get("Wall thickness"))
        ostiumOptions = options['Ileocecal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_cecum1.generateBaseMesh(region, options)[0]
        self.assertEqual(5, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(308, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1164, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1412, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(558, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-110.98815807992678, -144.1649444946355, 854.4533092097239], 1.0E-6)
        assertAlmostEqualList(self, maximums, [-55.3885994876074, -77.17764626881537, 900.1398770921413], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 8014.802826468518, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 12562.715396659782, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "caecum": 308,
            "ileum": 16,
            "ileocecal junction": 8
        }

        for name in expectedSizes3d:
            if name == "caecum":
                term = get_cecum_term(name)
            else:
                term = get_smallintestine_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)


if __name__ == "__main__":
    unittest.main()
