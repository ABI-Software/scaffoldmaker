import unittest
import copy
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_3d_solidsphere2 import MeshType_3d_solidsphere2
from testutils import assertAlmostEqualList
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class SphereScaffoldTestCase(unittest.TestCase):

    def test_sphere1(self):
        """
        Test creation of Sphere scaffold.
        """
        scaffold = MeshType_3d_solidsphere2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default"])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(16, len(options))
        self.assertEqual(4, options.get("Number of elements across axis 1"))
        self.assertEqual(4, options.get("Number of elements across axis 2"))
        self.assertEqual(4, options.get("Number of elements across axis 3"))
        self.assertEqual(0, options.get("Number of elements across shell"))
        self.assertEqual(1, options.get("Number of elements across transition"))
        self.assertEqual(1.0, options.get("Radius1"))
        self.assertEqual(1.0, options.get("Radius2"))
        self.assertEqual(1.0, options.get("Radius3"))
        self.assertEqual(1.0, options.get("Shell element thickness proportion"))
        self.assertEqual([-1, 2, 3], options.get("Box derivatives"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)
        self.assertEqual(2, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(32, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(108, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(128, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(53, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, sphere volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-1.0, -1.0, -1.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [1.0, 1.0, 1.0], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceGroup = AnnotationGroup(region, ("sphere surface", ""))
            is_exterior = fieldmodule.createFieldIsExterior()
            surfaceMeshGroup = surfaceGroup.getMeshGroup(mesh2d)
            surfaceMeshGroup.addElementsConditional(is_exterior)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 12.460033954564986, delta=2.0E-1)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 4.132033912594377, delta=1.0E-3)

        # check some annotationGroups:
        expectedSizes3d = {
            "box group": 8,
            "transition group": 24,
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, (name, ''))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 8x8x8 and check result
        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, [])
        scaffold.refineMesh(meshrefinement, options)

        refineCoordinates = refineFieldmodule.findFieldByName("coordinates").castFiniteElement()

        refineFieldmodule.defineAllFaces()
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(256, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(816, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(880, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(321, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # obtain surface area and volume again. If they are the same as previous, mesh is conform.
        with ChangeManager(refineFieldmodule):
            one = refineFieldmodule.createFieldConstant(1.0)
            surfaceGroup = AnnotationGroup(refineRegion, ("sphere surface", ""))
            is_exterior = refineFieldmodule.createFieldIsExterior()
            surfaceMeshGroup = surfaceGroup.getMeshGroup(mesh2d)
            surfaceMeshGroup.addElementsConditional(is_exterior)

            surfaceAreaField = refineFieldmodule.createFieldMeshIntegral(one, refineCoordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = refineFieldmodule.createFieldMeshIntegral(one, refineCoordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        refineFieldcache = refineFieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(refineFieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 12.460033954564986, delta=5.0E-1)
        result, volume = volumeField.evaluateReal(refineFieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 4.132033912594377, delta=3.0E-1)


        # check ellipsoid.
        scaffold1 = MeshType_3d_solidsphere2
        parameterSetNames = scaffold1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default"])
        options = scaffold1.getDefaultOptions("Default")
        options['Number of elements across axis 1'] = 4
        options['Number of elements across axis 2'] = 6
        options['Number of elements across axis 3'] = 8

        options['Radius1'] = 0.5
        options['Radius2'] = 0.8
        options['Radius3'] = 1.0

        context2 = Context("Test2")
        region = context2.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold1.generateMesh(region, options)
        self.assertEqual(2, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(136, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(452, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(510, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(195, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, sphere volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5, -0.8, -1.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.5, 0.8, 1.0], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceGroup = AnnotationGroup(region, ("sphere surface", ""))
            is_exterior = fieldmodule.createFieldIsExterior()
            surfaceMeshGroup = surfaceGroup.getMeshGroup(mesh2d)
            surfaceMeshGroup.addElementsConditional(is_exterior)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 7.3015173697377245, delta=2.0E-1)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 1.6741674010981926, delta=1.0E-3)

if __name__ == "__main__":
    unittest.main()
