import math
import unittest

from cmlibs.maths.vectorops import magnitude
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.group import identifier_ranges_to_string, mesh_group_add_identifier_ranges, \
    mesh_group_to_identifier_ranges
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_2d_tubenetwork1 import MeshType_2d_tubenetwork1
from scaffoldmaker.meshtypes.meshtype_3d_boxnetwork1 import MeshType_3d_boxnetwork1
from scaffoldmaker.meshtypes.meshtype_3d_tubenetwork1 import MeshType_3d_tubenetwork1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters

from testutils import assertAlmostEqualList


class NetworkScaffoldTestCase(unittest.TestCase):

    def test_network_layout(self):
        """
        Test creation of network layout scaffold.
        """
        scaffold = MeshType_1d_network_layout1
        options = scaffold.getDefaultOptions()
        self.assertEqual(3, len(options))
        self.assertEqual("1-2", options.get("Structure"))
        options["Structure"] = "1-2-3,3-4,3.2-5"

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups, networkMesh = scaffold.generateBaseMesh(region, options)
        self.assertEqual(0, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(4, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(5, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        interactiveFunctions = scaffold.getInteractiveFunctions()
        functionOptions = None
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Smooth derivatives...":
                functionOptions = interactiveFunction[1]
                break
        functionOptions["Update directions"] = True
        scaffold.smoothDerivatives(region, options, None, functionOptions, "meshEdits")

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [3.0, 0.5, 0.0], 1.0E-6)

        networkSegments = networkMesh.getNetworkSegments()
        self.assertEqual(3, len(networkSegments))
        self.assertEqual([1, 2, 3], networkSegments[0].getNodeIdentifiers())
        self.assertEqual([1, 1, 1], networkSegments[0].getNodeVersions())
        self.assertEqual([3, 4], networkSegments[1].getNodeIdentifiers())
        self.assertEqual([1, 1], networkSegments[1].getNodeVersions())
        self.assertEqual([3, 5], networkSegments[2].getNodeIdentifiers())
        self.assertEqual([2, 1], networkSegments[2].getNodeVersions())

        # get path parameters with versions
        nx, nd1 = get_nodeset_path_ordered_field_parameters(
            nodes, coordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1],
            networkSegments[2].getNodeIdentifiers(), networkSegments[2].getNodeVersions())
        self.assertEqual(2, len(nx))
        assertAlmostEqualList(self, nx[0], [2.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, nx[1], [3.0, 0.5, 0.0], 1.0E-6)
        expected_nd = [nx[1][c] - nx[0][c] for c in range(3)]
        assertAlmostEqualList(self, nd1[0], expected_nd, 1.0E-6)
        assertAlmostEqualList(self, nd1[1], expected_nd, 1.0E-6)

    def test_2d_tube_network_bifurcation(self):
        """
        Test 2D tube bifurcation is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_2d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        self.assertFalse(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(6, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        settings["Target element density along longest segment"] = 3.3
        MeshType_2d_tubenetwork1.checkOptions(settings)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(88, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(99, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5894427190999916, -0.10000000000000002], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, 1.9298310795249618, delta=X_TOL)

    def test_2d_tube_network_snake(self):
        """
        Test 2D tube snake has radial elements.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_2d_tubenetwork1, defaultParameterSetName="Snake")
        settings = scaffoldPackage.getScaffoldSettings()
        self.assertEqual(12.0, settings["Target element density along longest segment"])
        MeshType_2d_tubenetwork1.checkOptions(settings)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(96, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(104, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.1, -0.5196340284402325, -0.1], X_TOL)
        assertAlmostEqualList(self, maximums, [4.1, 0.5196340284402319, 0.1], X_TOL)

        with ChangeManager(fieldmodule):
            # check range of d2 shows element sizes vary from inside to outside of curves
            d2 = fieldmodule.createFieldNodeValue(coordinates, Node.VALUE_LABEL_D_DS2, 1)
            mag_d2 = fieldmodule.createFieldMagnitude(d2)
            min_mag_d2, max_mag_d2 = evaluateFieldNodesetRange(mag_d2, nodes)

            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(min_mag_d2, 0.41678801141467386, delta=X_TOL)
            self.assertAlmostEqual(max_mag_d2, 0.6251820171220115, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 3.883499820061501, delta=X_TOL)

    def test_2d_tube_network_sphere_cube(self):
        """
        Test 2D sphere cube is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_2d_tubenetwork1, defaultParameterSetName="Sphere cube")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        self.assertFalse(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(6, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])

        context = Context("Test")
        region = context.getDefaultRegion()

        # add a user-defined annotation group to network layout. Must generate first
        tmpRegion = region.createRegion()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        networkLayoutScaffoldPackage.generate(tmpRegion)

        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("bob", "BOB:1"))
        group = annotationGroup1.getGroup()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[1, 1], [5, 5]])
        self.assertEqual(2, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("1,5", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(1, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(32 * 12, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(8 * 7 * 12 + 4 * 3 * 8, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(8 * 3 * 12 + (2 + 3 * 3) * 8, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        # check annotation group transferred to 2D tube
        annotationGroup = annotationGroups[0]
        self.assertEqual("bob", annotationGroup.getName())
        self.assertEqual("BOB:1", annotationGroup.getId())
        self.assertEqual(64, annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(2)).getSize())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5665335420558559, -0.5965021612011158, -0.5986833971069179], X_TOL)
        assertAlmostEqualList(self, maximums, [0.5665335420558559, 0.5965021612011159, 0.5986833971069178], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, 4.040053575361621, delta=X_TOL)

    def test_2d_tube_network_trifurcation(self):
        """
        Test 2D tube triifurcation is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_2d_tubenetwork1, defaultParameterSetName="Trifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        self.assertFalse(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(6, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        MeshType_2d_tubenetwork1.checkOptions(settings)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(112, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(126, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -1.0707106781186548, -0.10000000000000002], X_TOL)
        assertAlmostEqualList(self, maximums, [2.0707106781186546, 1.0707106781186548, 0.1], X_TOL)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, 2.791976505648992, delta=X_TOL)

    def test_2d_tube_network_vase(self):
        """
        Test 2D tube vase has near constant length elements despite radius changes.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_2d_tubenetwork1, defaultParameterSetName="Vase")
        settings = scaffoldPackage.getScaffoldSettings()
        self.assertEqual(12.0, settings["Target element density along longest segment"])
        MeshType_2d_tubenetwork1.checkOptions(settings)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(96, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(104, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-1.5, -1.5, 0.0], X_TOL)
        assertAlmostEqualList(self, maximums, [1.5, 1.5, 4.0], X_TOL)

        with ChangeManager(fieldmodule):
            # check range of d2 shows near constant element sizes
            d2 = fieldmodule.createFieldNodeValue(coordinates, Node.VALUE_LABEL_D_DS2, 1)
            mag_d2 = fieldmodule.createFieldMagnitude(d2)
            min_mag_d2, max_mag_d2 = evaluateFieldNodesetRange(mag_d2, nodes)

            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(min_mag_d2, 0.38259775266954776, delta=X_TOL)
            self.assertAlmostEqual(max_mag_d2, 0.3825977526695479, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 28.820994366312384, delta=X_TOL)

    def test_3d_tube_network_bifurcation(self):
        """
        Test bifurcation 3-D tube network is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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
        self.assertEqual(8 * 4 * 3, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 * 4 * 3 + 3 * 3 + 2) * 2, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5894427190999916, -0.10000000000000002], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

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

            self.assertAlmostEqual(volume, 0.0349284962620723, delta=X_TOL)
            self.assertAlmostEqual(outerSurfaceArea, 1.9284739468709864, delta=X_TOL)
            self.assertAlmostEqual(innerSurfaceArea, 1.5595435883893172, delta=X_TOL)

    def test_3d_tube_network_bifurcation_core(self):
        """
        Test bifurcation 3-D tube network with solid core is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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

        self.assertEqual((8 * 4 * 3) * 2 + (4 * 4 * 3), mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual((8 * 4 * 3 + 3 * 3 + 2) * 2 + (9 * 4 * 3 + 3 * 4), nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5894427190999916, -0.10000000000000002], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

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

            self.assertAlmostEqual(volume, 0.09883609668362349, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 2.0226236083210507, delta=X_TOL)

    def test_3d_tube_network_converging_bifurcation_core(self):
        """
        Test converging bifurcation 3-D tube network with solid core and 12, 12, 8 elements around.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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
        self.assertEqual(336, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(460, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        # check annotation group transferred to 3D tube
        annotationGroup = findAnnotationGroupByName(annotationGroups, "segment 3")
        self.assertTrue(annotationGroup is not None)
        self.assertEqual("SEGMENT:3", annotationGroup.getId())
        self.assertEqual(80, annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(3)).getSize())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5894427190999916, -0.10000000000000002], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

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

            self.assertAlmostEqual(volume, 0.09854306375590477, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 2.0220707879327655, delta=X_TOL)

    def test_3d_tube_network_line_core_transition2(self):
        """
        Test line 3-D tube network with solid core and 2 transition elements.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Default")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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
        settings["Number of elements across core transition"] = 2

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)

        self.assertEqual(112, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(165, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.1, -0.1], X_TOL)
        assertAlmostEqualList(self, maximums, [1.0, 0.1, 0.1], X_TOL)

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

            self.assertAlmostEqual(volume, 0.0313832204833548, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 0.6907602069977625, delta=X_TOL)

    def test_3d_tube_network_sphere_cube(self):
        """
        Test sphere cube 3-D tube network is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Sphere cube")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertFalse(settings["Show trim surfaces"])
        settings["Number of elements through shell"] = 2
        settings["Use linear through shell"] = True

        context = Context("Test")
        region = context.getDefaultRegion()

        # set custom inner coordinates
        tmpRegion = region.createRegion()
        networkLayoutScaffoldPackage.generate(tmpRegion)
        networkMesh = networkLayoutScaffoldPackage.getConstructionObject()
        functionOptions = {
            "To field": {"coordinates": False, "inner coordinates": True},
            "From field": {"coordinates": True, "inner coordinates": False},
            "Mode": {"Scale": True, "Offset": False},
            "D2 value": 0.8,
            "D3 value": 0.8}
        editGroupName = "meshEdits"
        MeshType_1d_network_layout1.assignCoordinates(tmpRegion, networkLayoutSettings, networkMesh,
                                                      functionOptions, editGroupName=editGroupName)
        # put edited coordinates into scaffold package
        sir = tmpRegion.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        sir.setResourceGroupName(srm, editGroupName)
        sir.setResourceFieldNames(srm, ["coordinates", "inner coordinates"])
        tmpRegion.write(sir)
        result, meshEditsString = srm.getBuffer()
        self.assertEqual(RESULT_OK, result)
        networkLayoutScaffoldPackage.setMeshEdits(meshEditsString)

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(32 * 12 * 2, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(8 * 3 * 12 * 3 + (2 + 3 * 3) * 8 * 3, nodes.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(32 * 12 * 5 + 24 * 12 * 2 + 12 * 8 * 2, mesh2d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5665130262270113, -0.5965021612011158, -0.5986773876363235], X_TOL)
        assertAlmostEqualList(self, maximums, [0.5665130262270113, 0.5965021612011159, 0.5986773876363234], X_TOL)

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

            self.assertAlmostEqual(volume, 0.07331492814968889, delta=X_TOL)
            self.assertAlmostEqual(outerSurfaceArea, 4.04004649646839, delta=X_TOL)
            self.assertAlmostEqual(innerSurfaceArea, 3.325814823258647, delta=X_TOL)

    def test_3d_tube_network_sphere_cube_core(self):
        """
        Test sphere cube 3-D tube network with solid core is generated correctly.
        Use different number of elements around on some segments to mix it up.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Sphere cube")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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
        settings["Number of elements through shell"] = 2
        settings["Annotation numbers of elements around"] = [12]
        settings["Core"] = True

        context = Context("Test")
        region = context.getDefaultRegion()

        # set custom inner coordinates
        tmpRegion = region.createRegion()
        networkLayoutScaffoldPackage.generate(tmpRegion)
        networkMesh = networkLayoutScaffoldPackage.getConstructionObject()
        functionOptions = {
            "To field": {"coordinates": False, "inner coordinates": True},
            "From field": {"coordinates": True, "inner coordinates": False},
            "Mode": {"Scale": True, "Offset": False},
            "D2 value": 0.75,
            "D3 value": 0.75}
        editGroupName = "meshEdits"
        MeshType_1d_network_layout1.assignCoordinates(tmpRegion, networkLayoutSettings, networkMesh,
                                                      functionOptions, editGroupName=editGroupName)
        # put edited coordinates into scaffold package
        sir = tmpRegion.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        sir.setResourceGroupName(srm, editGroupName)
        sir.setResourceFieldNames(srm, ["coordinates", "inner coordinates"])
        tmpRegion.write(sir)
        result, meshEditsString = srm.getBuffer()
        self.assertEqual(RESULT_OK, result)
        networkLayoutScaffoldPackage.setMeshEdits(meshEditsString)

        # add a user-defined annotation group to network layout to vary elements count around. Must generate first
        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("group1", ""))
        group = annotationGroup1.getGroup()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[2, 2], [5, 5], [8, 8], [10, 10]])
        self.assertEqual(4, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("2,5,8,10", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(3, len(annotationGroups))
        self.assertTrue(findAnnotationGroupByName(annotationGroups, "core") is not None)
        self.assertTrue(findAnnotationGroupByName(annotationGroups, "shell") is not None)
        self.assertTrue(findAnnotationGroupByName(annotationGroups, "group1") is not None)

        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(1600, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1836, nodes.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(5024, mesh2d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5762017364104554, -0.5965021612011158, -0.5909194631968028], X_TOL)
        assertAlmostEqualList(self, maximums, [0.5762017364104554, 0.5965021612011158, 0.5909194631968027], X_TOL)

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

            self.assertAlmostEqual(volume, 0.20985014947157313, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 4.033518137808064, delta=X_TOL)

        expectedSizes3d = {
            "core": (704, 0.12129037249755871),
            "shell": (896, 0.08855977697401322)
            }
        for name in expectedSizes3d:
            annotationGroup = findAnnotationGroupByName(annotationGroups, name)
            size = annotationGroup.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0], size, name)
            volumeMeshGroup = annotationGroup.getMeshGroup(mesh3d)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, volumeMeshGroup)
            volumeField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(volume, expectedSizes3d[name][1], delta=X_TOL)

    def test_3d_tube_network_trifurcation_cross(self):
        """
        Test trifurcation cross 3-D tube network is generated correctly with variable elements count around.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Trifurcation cross")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        self.assertTrue(networkLayoutSettings["Define inner coordinates"])
        self.assertEqual(13, len(settings))
        self.assertEqual(8, settings["Number of elements around"])
        self.assertEqual(1, settings["Number of elements through shell"])
        self.assertEqual([0], settings["Annotation numbers of elements around"])
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])
        self.assertFalse(settings["Use linear through shell"])
        self.assertFalse(settings["Show trim surfaces"])
        settings["Annotation numbers of elements around"] = [10]  # requires annotation group below

        context = Context("Test")
        region = context.getDefaultRegion()

        # add a user-defined annotation group to network layout to vary elements count around. Must generate first
        tmpRegion = region.createRegion()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        networkLayoutScaffoldPackage.generate(tmpRegion)

        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("straight", "STRAIGHT:1"))
        group = annotationGroup1.getGroup()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[1, 1], [4, 4]])
        self.assertEqual(2, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("1,4", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(1, len(annotationGroups))

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(144, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(320, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        # check annotation group transferred to 3D tube
        annotationGroup = annotationGroups[0]
        self.assertEqual("straight", annotationGroup.getName())
        self.assertEqual("STRAIGHT:1", annotationGroup.getId())
        self.assertEqual(80, annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(3)).getSize())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.0447213595499958, -0.5894427190999916, -0.1], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

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

            self.assertAlmostEqual(volume, 0.047176020014987795, delta=X_TOL)
            self.assertAlmostEqual(outerSurfaceArea, 2.5987990744712053, delta=X_TOL)
            self.assertAlmostEqual(innerSurfaceArea, 2.113990124755675, delta=X_TOL)

    def test_3d_tube_network_trifurcation_cross_core(self):
        """
        Test trifurcation cross 3-D tube network with solid core is generated correctly with
        variable elements count around.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Trifurcation cross")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
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
        settings["Annotation numbers of elements around"] = [12]  # requires annotation group below
        settings["Annotation numbers of elements across core box minor"] = [2]

        context = Context("Test")
        region = context.getDefaultRegion()

        # add a user-defined annotation group to network layout to vary elements count around. Must generate first
        tmpRegion = region.createRegion()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        networkLayoutScaffoldPackage.generate(tmpRegion)

        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("straight", "STRAIGHT:1"))
        group = annotationGroup1.getGroup()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[1, 1], [4, 4]])
        self.assertEqual(2, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("1,4", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(3, len(annotationGroups))

        fieldmodule = region.getFieldmodule()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(416, mesh3d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(569, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        # check annotation group transferred to 3D tube
        annotationGroup = findAnnotationGroupByName(annotationGroups, "straight")
        self.assertTrue(annotationGroup is not None)
        self.assertEqual("STRAIGHT:1", annotationGroup.getId())
        self.assertEqual(256, annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(3)).getSize())

        X_TOL = 1.0E-6

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.0447213595499958, -0.5894427190999916, -0.1], X_TOL)
        assertAlmostEqualList(self, maximums, [2.044721359549996, 0.5894427190999916, 0.10000000000000002], X_TOL)

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

            self.assertAlmostEqual(volume, 0.1346381132294423, delta=X_TOL)
            self.assertAlmostEqual(surfaceArea, 2.7234652881096166, delta=X_TOL)

    def test_3d_box_network_bifurcation(self):
        """
        Test 3-D box network bifurcation is generated correctly.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_boxnetwork1, defaultParameterSetName="Bifurcation")
        settings = scaffoldPackage.getScaffoldSettings()
        self.assertEqual(3, len(settings))
        self.assertEqual(4.0, settings["Target element density along longest segment"])
        self.assertEqual([0], settings["Annotation numbers of elements along"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(12, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(63, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(108, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(13, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        X_TOL = 1.0E-6
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, -0.5, 0.0], X_TOL)
        assertAlmostEqualList(self, maximums, [2.0, 0.5, 0.0], X_TOL)

        L2 = math.sqrt(1.25)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(1)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            expectedVolume = 0.2 * 0.2 * (1.0 + 2 * L2)
            self.assertAlmostEqual(volume, expectedVolume, delta=X_TOL)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(1)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            expectedSurfaceArea = 6 * 0.2 * 0.2 + 4 * 0.2 * (1.0 + 2 * L2)
            self.assertAlmostEqual(surfaceArea, expectedSurfaceArea, delta=X_TOL)

    def test_3d_box_network_smooth(self):
        """
        Test 3-D box network derivative smoothing is working between segments sharing a version.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_boxnetwork1)
        settings = scaffoldPackage.getScaffoldSettings()
        settings["Target element density along longest segment"] = 1.0
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        networkLayoutSettings["Structure"] = "1-2-3,3-4"  # 2 unequal-sized segments

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(3, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        node2 = nodes.findNodeByIdentifier(2)
        self.assertTrue(node2.isValid())

        # test magnitude of d1 between segments is harmonic mean of element sizes
        fieldcache = fieldmodule.createFieldcache()
        fieldcache.setNode(node2)
        result, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        self.assertEqual(result, RESULT_OK)
        d1Mag = magnitude(d1)
        self.assertAlmostEqual(4.0 / 3.0, d1Mag, delta=1.0E-12)

    def test_3d_tube_network_loop(self):
        """
        Test loop 3-D tube network is generated correctly.
        This has one segment which loops back on itself so nodes are common at start and end.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Loop")
        settings = scaffoldPackage.getScaffoldSettings()
        settings["Target element density along longest segment"] = 8.0

        context = Context("Test")
        region = context.getDefaultRegion()
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(8 * 8, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(8 * 8 * 4, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(8 * 8 * 2, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.6, -0.6, -0.1], 1.0E-8)
        assertAlmostEqualList(self, maximums, [0.6, 0.6, 0.1], 1.0E-8)

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

            self.assertAlmostEqual(volume, 0.03534439013604324, delta=1.0E-6)
            self.assertAlmostEqual(outerSurfaceArea, 1.9683574196198823, delta=1.0E-6)
            self.assertAlmostEqual(innerSurfaceArea, 1.5748510621127434, delta=1.0E-6)

    def test_3d_tube_network_loop_core(self):
        """
        Test loop 3-D tube network with solid core is generated correctly.
        This has one segment which loops back on itself so nodes are common at start and end.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Loop")
        settings = scaffoldPackage.getScaffoldSettings()
        settings["Target element density along longest segment"] = 8.0
        settings["Core"] = True

        context = Context("Test")
        region = context.getDefaultRegion()
        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(160, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(512, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(200, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.6, -0.6, -0.1], 1.0E-8)
        assertAlmostEqualList(self, maximums, [0.6, 0.6, 0.1], 1.0E-8)

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

            self.assertAlmostEqual(volume, 0.0982033864405135, delta=1.0E-6)
            self.assertAlmostEqual(surfaceArea, 1.9683574196198823, delta=1.0E-6)

    def test_3d_tube_network_loop_two_segments(self):
        """
        Test loop 3-D tube network is generated with 2 segments with fixed element boundary between them.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_tubenetwork1, defaultParameterSetName="Loop")
        settings = scaffoldPackage.getScaffoldSettings()
        networkLayoutScaffoldPackage = settings["Network layout"]
        networkLayoutSettings = networkLayoutScaffoldPackage.getScaffoldSettings()
        # change structure to make two segments but use regular loop parameters:
        networkLayoutSettings["Structure"] = "1-2-3-4-5-6-7,7-8-1"
        settings["Target element density along longest segment"] = 7.0

        context = Context("Test")
        region = context.getDefaultRegion()

        # add a user-defined annotation group to network layout. Must generate first
        tmpRegion = region.createRegion()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        networkLayoutScaffoldPackage.generate(tmpRegion)

        annotationGroup1 = networkLayoutScaffoldPackage.createUserAnnotationGroup(("bob", "BOB:1"))
        group = annotationGroup1.getGroup()
        mesh1d = tmpFieldmodule.findMeshByDimension(1)
        meshGroup = group.createMeshGroup(mesh1d)
        mesh_group_add_identifier_ranges(meshGroup, [[7, 8]])
        self.assertEqual(2, meshGroup.getSize())
        self.assertEqual(1, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual("7-8", identifier_ranges_string)
        networkLayoutScaffoldPackage.updateUserAnnotationGroups()

        self.assertTrue(region.isValid())
        scaffoldPackage.generate(region)
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(1, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(10 * 8, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(10 * 8 * 4, mesh2d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(10 * 8 * 2, nodes.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())

        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5846409928643533, -0.6, -0.1], 1.0E-8)
        assertAlmostEqualList(self, maximums, [0.6, 0.5846409928643533, 0.1], 1.0E-8)

        bob = fieldmodule.findFieldByName("bob").castGroup()
        self.assertTrue(bob.isValid())
        bobNodes = bob.getNodesetGroup(nodes)
        self.assertTrue(bobNodes.isValid())
        self.assertEqual(4 * 8 * 2, bobNodes.getSize())
        bobMinimums, bobMaximums = evaluateFieldNodesetRange(coordinates, bobNodes)
        assertAlmostEqualList(self, bobMinimums, [0.0, -0.6, -0.1], 1.0E-8)
        assertAlmostEqualList(self, bobMaximums, [0.6, 0.0, 0.1], 1.0E-8)

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

            self.assertAlmostEqual(volume, 0.0353515741325893, delta=1.0E-6)
            self.assertAlmostEqual(outerSurfaceArea, 1.9681077595642782, delta=1.0E-6)
            self.assertAlmostEqual(innerSurfaceArea, 1.5745958498454014, delta=1.0E-6)


if __name__ == "__main__":
    unittest.main()
