import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.group import identifier_ranges_to_string, mesh_group_add_identifier_ranges, \
    mesh_group_to_identifier_ranges
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_tubenetwork1 import MeshType_3d_tubenetwork1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage

from testutils import assertAlmostEqualList

class CapScaffoldTestCase(unittest.TestCase):

    def test_3d_cap_tube_network_default(self):
        """
        Test default 3-D tube network with cap at both ends without core is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Default")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change the network layout to have cap at both ends of the tube
        networkLayoutSettings["Structure"] = "(1-2)"
        self.assertEqual("(1-2)", networkLayoutSettings["Structure"])

        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertTrue(settings["Use outer trim surfaces"])
        self.assertFalse(settings["Show trim surfaces"])
        self.assertFalse(settings["Core"])
        self.assertEqual(2, settings["Number of elements across core box minor"])
        self.assertEqual(1, settings["Number of elements across core transition"])
        self.assertEqual([0], settings["Annotation numbers of elements across core box minor"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(8 * 4 + 8 * 3 * 2, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 * 5 + 8 * 2 * 2 + 2) * 2, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.150, -0.100, -0.100], X_TOL)
        assertAlmostEqualList(self, maximums, [1.150, 0.100, 0.100], X_TOL)

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
            volumeField.setNumbersOfPoints(4)
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

            self.assertAlmostEqual(volume, 0.014526098773694766, delta=X_TOL)
            self.assertAlmostEqual(outerSurfaceArea, 0.8247195017487451, delta=X_TOL)
            self.assertAlmostEqual(innerSurfaceArea, 0.6421069263444478, delta=X_TOL)

    def test_3d_cap_tube_network_default_core(self):
        """
        Test default 3-D tube network with cap at both ends and with a solid core is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Default")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change the network layout to have cap at both ends of the tube
        networkLayoutSettings["Structure"] = "(1-2)"
        self.assertEqual("(1-2)",networkLayoutSettings["Structure"])

        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertTrue(settings["Use outer trim surfaces"])
        self.assertFalse(settings["Show trim surfaces"])
        self.assertFalse(settings["Core"])
        self.assertEqual(2, settings["Number of elements across core box minor"])
        self.assertEqual(1, settings["Number of elements across core transition"])
        self.assertEqual([0], settings["Annotation numbers of elements across core box minor"])
        settings["Core"] = True

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual((8 + 8 + 4) * 4 + ((8 + 8 + 4) * 2 + 4 + 4) * 2, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 + 8 + 9) * 5 + ((8 + 8 + 9) + 9 * 3) * 2, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.150, -0.100, -0.100], X_TOL)
        assertAlmostEqualList(self, maximums, [1.150, 0.100, 0.100], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(4)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 0.03862588577889281, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 0.8148830981004554, delta=X_TOL)

    def test_3d_cap_tube_network_bifurcation(self):
        """
        Test bifurcation 3-D tube network with cap at both ends is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change the network layout to have cap at both ends of the tube
        networkLayoutSettings["Structure"] = "(1-2.1,2.2-3),2.3-4)"
        self.assertEqual("(1-2.1,2.2-3),2.3-4)",networkLayoutSettings["Structure"])

        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertTrue(settings["Use outer trim surfaces"])
        self.assertFalse(settings["Show trim surfaces"])
        self.assertFalse(settings["Core"])
        self.assertEqual(2, settings["Number of elements across core box minor"])
        self.assertEqual(1, settings["Number of elements across core transition"])
        self.assertEqual([0], settings["Annotation numbers of elements across core box minor"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual((8 * 4 + 8 * 3) * 3 , mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 * 4 * 3 + 3 * 3 + 2) * 2 + (8 * 2 * 2 + 2) * 3, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.1500000000000000, -0.6246941344953365, -0.1000000000000000], X_TOL)
        assertAlmostEqualList(self, maximums, [2.1433222518126906, 0.6246941344953366, 0.1000000000000000], X_TOL)

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
            volumeField.setNumbersOfPoints(4)
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

            self.assertAlmostEqual(volume, 0.039770715282355665, delta=X_TOL)
            self.assertAlmostEqual(outerSurfaceArea, 2.2235629407464303, delta=X_TOL)
            self.assertAlmostEqual(innerSurfaceArea, 1.7691117020959595, delta=X_TOL)

    def test_3d_tube_network_bifurcation_core(self):
        """
        Test bifurcation 3-D tube network with solid core and cap at both ends is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change the network layout to have cap at both ends of the tube
        networkLayoutSettings["Structure"] = "(1-2.1,2.2-3),2.3-4)"
        self.assertEqual("(1-2.1,2.2-3),2.3-4)",networkLayoutSettings["Structure"])

        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertFalse(settings["Show trim surfaces"])
        self.assertFalse(settings["Core"])
        self.assertEqual(2, settings["Number of elements across core box minor"])
        self.assertEqual(1, settings["Number of elements across core transition"])
        self.assertEqual([0], settings["Annotation numbers of elements across core box minor"])
        settings["Core"] = True

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual((8 * 4 * 3) * 2 + (4 * 4 * 3) +((8 + 8 + 4) * 2 + 4 + 4) * 3, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 * 4 * 3 + 3 * 3 + 2) * 2 + (9 * 4 * 3 + 3 * 4) + ((8 + 8 + 9) + 9 * 3) * 3, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.1500000000000000, -0.6172290095800493, -0.1000000000000000], X_TOL)
        assertAlmostEqualList(self, maximums, [2.1395896893550477, 0.6172290095800493, 0.1000000000000000], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(4)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 0.1096893209534004, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 2.2086099452749344, delta=X_TOL)

    def test_3d_tube_network_converging_bifurcation_core(self):
        """
        Test converging bifurcation 3-D tube network with solid core and 12, 12, 8 elements around.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change the network layout to have cap at both ends of the tube
        networkLayoutSettings["Structure"] = "(1-3.1,(2-3.2,3.3-4)"
        self.assertEqual("(1-3.1,(2-3.2,3.3-4)",networkLayoutSettings["Structure"])

        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertFalse(settings["Show trim surfaces"])
        self.assertFalse(settings["Core"])
        self.assertEqual(2, settings["Number of elements across core box minor"])
        self.assertEqual(1, settings["Number of elements across core transition"])
        self.assertEqual([0], settings["Annotation numbers of elements across core box minor"])
        settings["Core"] = True
        settings["Number of elements around"] = 12
        settings["Annotation numbers of elements around"] = [8]

        context = Context("Test")
        region = context.getDefaultRegion()

        # add a user-defined annotation group to network layout to vary elements count around. Must generate first
        tmpRegion = region.createRegion()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        networkLayoutScaffoldPackage.generate(tmpRegion)

        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("segment 3", "SEGMENT:3"))
        group = annotationGroup1.getGroup()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[3, 3]])
        self.assertEqual(1, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("3", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(3, len(annotationGroups))
        self.assertTrue(findAnnotationGroupByName(annotationGroups, "core") is not None)
        self.assertTrue(findAnnotationGroupByName(annotationGroups, "shell") is not None)

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(544, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(680, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        # check annotation group transferred to 3D tube
        annotationGroup = findAnnotationGroupByName(annotationGroups, "segment 3")
        self.assertTrue(annotationGroup is not None)
        self.assertEqual("SEGMENT:3", annotationGroup.getId())
        self.assertEqual(128, annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(3)).getSize())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.14453769545049167, -0.6221797157655231, -0.1000000000000000], X_TOL)
        assertAlmostEqualList(self, maximums, [2.15, 0.6221795584199485, 0.1000000000000000], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(4)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 0.10923734449979247, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 2.205681406661578, delta=X_TOL)

if __name__ == "__main__":
    unittest.main()
