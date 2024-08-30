import math
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager

from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm, findAnnotationGroupByName
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.meshtypes.meshtype_3d_wholebody2 import MeshType_3d_wholebody2

from testutils import assertAlmostEqualList


class WholeBody2ScaffoldTestCase(unittest.TestCase):

    def test_wholebody2_core(self):
        """
        Test creation of Whole-body scaffold with solid core.
        """
        scaffold = MeshType_3d_wholebody2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Human 2 Coarse", "Human 2 Medium",
                                             "Human 2 Fine"])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(15, len(options))
        self.assertEqual(12, options["Number of elements around head"])
        self.assertEqual(12, options["Number of elements around torso"])
        self.assertEqual(8, options["Number of elements around arm"])
        self.assertEqual(8, options["Number of elements around leg"])
        self.assertEqual(1, options["Number of elements through wall"])
        self.assertEqual(5.0, options["Target element density along longest segment"])
        self.assertEqual(False, options["Show trim surfaces"])
        self.assertEqual(True, options["Use Core"])
        self.assertEqual(6, options["Number of elements across head core major"])
        self.assertEqual(6, options["Number of elements across torso core major"])
        self.assertEqual(4, options["Number of elements across arm core major"])
        self.assertEqual(4, options["Number of elements across leg core major"])
        self.assertEqual(1, options["Number of elements across core transition"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(5, len(annotationGroups))  # Needs updating as we add more annotation groups

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(784, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(2562, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(2809, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1032, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-6
        assertAlmostEqualList(self, minimums, [0.0, -1.976, -0.6113507244855241], tol)
        assertAlmostEqualList(self, maximums, [8.748998777022969, 1.968, 0.7035640657448223], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            isExteriorXi3_1 = fieldmodule.createFieldAnd(
                isExterior, fieldmodule.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExteriorXi3_1, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 6.419528010968665, delta=tol)
            self.assertAlmostEqual(surfaceArea, 35.880031911102506, delta=tol)

        # check some annotationGroups:
        expectedSizes2d = {
            "skin epidermis": (308, 35.880031911102506)
            }
        for name in expectedSizes2d:
            term = get_body_term(name)
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


    def test_wholebody2_tube(self):
        """
        Test creation of Whole-body scaffold without solid core.
        """
        scaffold = MeshType_3d_wholebody2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Human 2 Coarse", "Human 2 Medium",
                                             "Human 2 Fine"])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(15, len(options))
        self.assertEqual(12, options["Number of elements around head"])
        self.assertEqual(12, options["Number of elements around torso"])
        self.assertEqual(8, options["Number of elements around arm"])
        self.assertEqual(8, options["Number of elements around leg"])
        self.assertEqual(1, options["Number of elements through wall"])
        self.assertEqual(5.0, options["Target element density along longest segment"])
        self.assertEqual(False, options["Show trim surfaces"])
        self.assertEqual(True, options["Use Core"])
        self.assertEqual(6, options["Number of elements across head core major"])
        self.assertEqual(6, options["Number of elements across torso core major"])
        self.assertEqual(4, options["Number of elements across arm core major"])
        self.assertEqual(4, options["Number of elements across leg core major"])
        self.assertEqual(1, options["Number of elements across core transition"])

        options["Use Core"] = False
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(5, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(308, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1254, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1603, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(654, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-6
        assertAlmostEqualList(self, minimums, [0.0, -1.976, -0.6113507244855241], tol)
        assertAlmostEqualList(self, maximums, [8.748998777022969, 1.968, 0.7035640657448223], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            isExteriorXi3_0 = fieldmodule.createFieldAnd(
                isExterior, fieldmodule.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
            isExteriorXi3_1 = fieldmodule.createFieldAnd(
                isExterior, fieldmodule.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            outerSurfaceAreaField = fieldmodule.createFieldMeshIntegral(isExteriorXi3_1, coordinates, mesh2d)
            outerSurfaceAreaField.setNumbersOfPoints(4)
            result, outerSurfaceArea = outerSurfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            innerSurfaceAreaField = fieldmodule.createFieldMeshIntegral(isExteriorXi3_0, coordinates, mesh2d)
            innerSurfaceAreaField.setNumbersOfPoints(4)
            result, innerSurfaceArea = innerSurfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 3.1692513375640607, delta=tol)
            self.assertAlmostEqual(outerSurfaceArea, 35.880031911102506, delta=tol)
            self.assertAlmostEqual(innerSurfaceArea, 25.851093527837623, delta=tol)

        # check some annotationGroups:
        expectedSizes2d = {
            "skin epidermis": (308, 35.880031911102506)
            }
        for name in expectedSizes2d:
            term = get_body_term(name)
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
