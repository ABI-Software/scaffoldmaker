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
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Human 2", "Pig 1"])
        options = MeshType_3d_cecum1.getDefaultOptions("Human 2")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3.2, 4-3-5", networkLayoutSettings["Structure"])

        self.assertEqual(29, len(options))
        self.assertEqual(1, options.get("Number of segments"))
        self.assertEqual(2, options.get("Number of elements around tenia coli"))
        self.assertEqual(12, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(0.536, options.get("Corner outer radius factor"))
        self.assertEqual(0.464, options.get("Haustrum outer radius factor"))
        self.assertEqual(3.0, options.get("Segment length mid derivative factor"))
        self.assertEqual(3, options.get("Number of tenia coli"))
        self.assertEqual(10.0, options.get("Start tenia coli width"))
        self.assertEqual(0.0, options.get("End tenia coli width derivative"))
        self.assertEqual(1.6, options.get("Wall thickness"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_cecum1.generateBaseMesh(region, options)[0]
        self.assertEqual(7, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(436, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1642, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1982, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(776, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-112.41968159644387, -146.34422797225153, 852.6082677676069], 1.0E-6)
        assertAlmostEqualList(self, maximums, [-54.22370705290374, -77.56, 899.9973429272325], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 8546.983090282285, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 13790.25181377472, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "caecum": 436,
            "ileum": 12,
            "ileocecal junction": 12
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
