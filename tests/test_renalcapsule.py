import math
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager

from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.meshtypes.meshtype_3d_renal_capsule1 import MeshType_3d_renal_capsule1


from testutils import assertAlmostEqualList


class RenalCapsulecaffoldTestCase(unittest.TestCase):

    def test_renalcapsule(self):
        """
        Test creation of renal capsule scaffold.
        """
        scaffold = MeshType_3d_renal_capsule1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1"])
        options = scaffold.getDefaultOptions("Human 1")

        self.assertEqual(9, len(options))
        self.assertEqual(12, options["Elements count around"])
        self.assertEqual(1, options["Elements count through shell"])
        self.assertEqual([0], options["Annotation elements counts around"])
        self.assertEqual(4.0, options["Target element density along longest segment"])
        self.assertEqual(2, options["Number of elements across core box minor"])
        self.assertEqual(1, options["Number of elements across core transition"])
        self.assertEqual([0], options["Annotation numbers of elements across core box minor"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(3, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(288, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(920, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(994, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(363, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-4
        assertAlmostEqualList(self, minimums, [-1.583346623141804, -0.9520066012170885, -0.75], tol)
        assertAlmostEqualList(self, maximums, [1.583349401938375, 0.7499999999986053, 0.75], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 4.844335733470136, delta=tol)
            self.assertAlmostEqual(surfaceArea, 15.289023987470623, delta=tol)

        # check some annotation groups:

        expectedSizes3d = {
            "core": (176, 2.880428953529323),
            "shell": (112, 1.9640545714444255),
            "kidney capsule": (288, 4.844483524973759)
            }
        for name in expectedSizes3d:
            term = get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0], size, name)
            volumeMeshGroup = annotationGroup.getMeshGroup(mesh3d)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, volumeMeshGroup)
            volumeField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(volume, expectedSizes3d[name][1], delta=tol)

        expectedSizes2d = {
            "shell": (448, 37.8673195697525),
            "kidney capsule": (920, 66.45790167672409)
            }
        for name in expectedSizes2d:
            term = get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name][0], size, name)
            surfaceMeshGroup = annotationGroup.getMeshGroup(mesh2d)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, expectedSizes2d[name][1], delta=tol)

if __name__ == "__main__":
    unittest.main()
