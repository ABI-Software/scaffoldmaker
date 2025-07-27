import math
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager

from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.annotation.ureter_terms import get_ureter_term
from scaffoldmaker.meshtypes.meshtype_3d_renal_pelvis1 import MeshType_3d_renal_pelvis1


from testutils import assertAlmostEqualList


class RenalPelviscaffoldTestCase(unittest.TestCase):

    def test_renalpelvis(self):
        """
        Test creation of renal pelvis scaffold.
        """
        scaffold = MeshType_3d_renal_pelvis1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1"])
        options = scaffold.getDefaultOptions("Human 1")

        self.assertEqual(9, len(options))
        self.assertEqual(8, options["Elements count around"])
        self.assertEqual(1, options["Elements count through shell"])
        self.assertEqual([0], options["Annotation elements counts around"])
        self.assertEqual(4.0, options["Target element density along longest segment"])
        self.assertEqual(False, options["Use linear through shell"])
        self.assertEqual(True, options["Use outer trim surfaces"])
        self.assertEqual(False, options["Show trim surfaces"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(7, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(672, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(2532, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(3063, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1204, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-4
        assertAlmostEqualList(self, minimums, [0.7789485582121838, -2.250268470689302, -0.1318007396451522], tol)
        assertAlmostEqualList(self, maximums, [4.188222198878539, 2.250268470689302, 0.1318007396451522], tol)

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

            self.assertAlmostEqual(volume, 0.41239028202973105, delta=tol)
            self.assertAlmostEqual(surfaceArea, 15.14363706045464, delta=tol)

        # check some annotation groups:

        expectedSizes3d = {
            "core": (240, 0.19386743947037463),
            "major calyx": (64, 0.025523445233479505),
            "minor calyx": (160, 0.039595947108606436),
            "renal pelvis": (272, 0.10994753925440572),
            "renal pyramid": (480, 0.3199459642092786),
            "shell": (432, 0.21852600996918725),
            "ureter": (64, 0.049524462117630445)
            }
        for name in expectedSizes3d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
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
            "major calyx": (272, 3.0306454318904343),
            "minor calyx": (692, 4.858039077961195),
            "renal pelvis": (1132, 11.555677753247785),
            "renal pyramid": (1780, 18.695951082914355),
            "ureter": (264, 4.2985089237162715)
            }
        for name in expectedSizes2d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
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
