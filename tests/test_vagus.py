import os
import unittest

from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.meshtypes.meshtype_3d_nerve1 import MeshType_3d_nerve1, get_left_vagus_marker_locations_list
from scaffoldmaker.utils.read_vagus_data import VagusInputData

from testutils import assertAlmostEqualList


here = os.path.abspath(os.path.dirname(__file__))


class VagusScaffoldTestCase(unittest.TestCase):

    def test_input_vagus_data(self):
        """
        Reading vagus input data
        """

        data_file = os.path.join(here, "resources", "vagus_test_data1.exf")

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')
        assert(region.isValid())
        data_region = region.getParent().createChild('data')
        assert(data_region.isValid())
        result = data_region.readFile(data_file)
        assert result == RESULT_OK

        vagus_data = VagusInputData(data_region)
        self.assertEqual(vagus_data.get_side_label(), 'left')

        marker_data = vagus_data.get_level_markers()
        self.assertEqual(len(marker_data), 4)
        self.assertTrue('left level of superior border of the clavicle on the vagus nerve' in marker_data)

        orientation_data = vagus_data.get_orientation_data()
        self.assertEqual(len(orientation_data), 8)
        expected_orientation_info = [
            ("orientation left", 1),
            ("orientation left anterior", 2),
            ("orientation anterior", 12),
            ("orientation right anterior", 1),
            ("orientation right", 1),
            ("orientation right posterior", 1),
            ("orientation posterior", 1),
            ("orientation left posterior", 1),
        ]
        for orientation_direction_name, expected_count in expected_orientation_info:
            self.assertEqual(len(orientation_data[orientation_direction_name]), expected_count)

        trunk_group_name = vagus_data.get_trunk_group_name()
        self.assertEqual(trunk_group_name, "left vagus nerve")
        trunk_coordinates = vagus_data.get_trunk_coordinates()
        self.assertEqual(len(trunk_coordinates), 201)
        annotation_term_map = vagus_data.get_annotation_term_map()
        self.assertTrue(trunk_group_name in annotation_term_map)
        # self.assertEqual(annotation_term_map[trunk_group_name], 'http://purl.obolibrary.org/obo/UBERON_0035020')

        branch_data = vagus_data.get_branch_data()
        self.assertEqual(len(branch_data), 4)
        self.assertTrue("left superior laryngeal nerve A" in branch_data)
        self.assertEqual(len(branch_data["left superior laryngeal nerve A"]), 42)
        self.assertTrue("branch of left superior laryngeal nerve A" in branch_data)
        self.assertEqual(len(branch_data["branch of left superior laryngeal nerve A"]), 22)
        left_thoracic_cardiopulmonary_branches = (
            "left thoracic cardiopulmonary branch A", "left thoracic cardiopulmonary branch B")
        for branch_name in left_thoracic_cardiopulmonary_branches:
            self.assertTrue(branch_name in branch_data)

        branch_parents = vagus_data.get_branch_parent_map()
        self.assertEqual(branch_parents["left superior laryngeal nerve A"], trunk_group_name)
        self.assertEqual(branch_parents["branch of left superior laryngeal nerve A"], "left superior laryngeal nerve A")
        for branch_name in left_thoracic_cardiopulmonary_branches:
            self.assertEqual(branch_parents[branch_name], trunk_group_name)

        branch_common_groups = vagus_data.get_branch_common_group_map()
        self.assertEqual(len(branch_common_groups), 1)
        self.assertTrue("left thoracic cardiopulmonary branches" in branch_common_groups)
        for branch_name in left_thoracic_cardiopulmonary_branches:
            self.assertTrue(branch_name in branch_common_groups["left thoracic cardiopulmonary branches"])

    def test_no_input_file(self):
        """
        No input file.
        """

        scaffold = MeshType_3d_nerve1
        options = scaffold.getDefaultOptions("Human Right Vagus 1")

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')

        # check that it doesn't crash with  no input file
        annotation_groups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(annotation_groups, [])

    def test_vagus_nerve_1(self):
        """
        Test creation of vagus nerve scaffold with simple, synthetic data similar to REVA data.
        """
        scaffold = MeshType_3d_nerve1
        scaffoldname = MeshType_3d_nerve1.getName()
        self.assertEqual(scaffoldname, '3D Nerve 1')
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human Left Vagus 1', 'Human Right Vagus 1'])
        options = scaffold.getDefaultOptions("Human Left Vagus 1")
        self.assertEqual(len(options), 6)
        self.assertEqual(options.get('Number of elements along the trunk pre-fit'), 20)
        self.assertEqual(options.get('Number of elements along the trunk'), 50)
        self.assertEqual(options.get('Trunk proportion'), 1.0)
        self.assertEqual(options.get('Trunk fit number of iterations'), 5)
        self.assertEqual(options.get('Default trunk diameter mm'), 3.0)
        self.assertEqual(options.get('Branch diameter trunk proportion'), 0.5)
        # change options to make test fast and consistent, with minor effect on result:
        options['Number of elements along the trunk pre-fit'] = 10
        options['Number of elements along the trunk'] = 25
        options['Trunk fit number of iterations'] = 2

        context = Context("Test")
        root_region = context.getDefaultRegion()
        region = root_region.createChild('vagus')
        data_region = root_region.createChild('data')
        data_file = os.path.join(here, "resources", "vagus_test_data1.exf")
        self.assertEqual(data_region.readFile(data_file), RESULT_OK)

        # check annotation groups
        annotation_groups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(len(annotation_groups), 16)

        # (parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3,
        #  expected_surface_area, expected_volume)
        expected_group_info = {
            "left vagus nerve": (
                None, 25,
                [-1269.8293183474027, -6359.8794918062185, -69.79596358103436],
                [2163.5065559034797, -1112.4694863693492, 121.49477098496129],
                [49.69809961146075, 258.19137709267125, 1479.1407307248198],
                248424548.42757902,
                32973169952.33672),
            "left superior laryngeal nerve A": (
                "left vagus nerve", 3,
                [5888.361608956813, -4417.766403817663, -201.35981775061967],
                [-1488.8488795282458, 791.7429838047342, 83.0590111566718],
                [40.03987180084414, 14.864881979714482, 576.0260276199538],
                9712700.895177323,
                554689189.2182714),
            "branch of left superior laryngeal nerve A": (
                "left superior laryngeal nerve A", 2,
                [5095.698554751851, -1435.5245960200803, -4.95378510864149],
                [-1310.3674569985124, 264.340042537044, 53.75707581244285],
                [9.731322981637277, -12.682910300344247, 299.5737724490257],
                4680192.321708083,
                236055059.57983524),
            "left thoracic cardiopulmonary branch A": (
                "left vagus nerve", 2,
                [20651.94765027247, -2955.4531267489256, -609.718731363599],
                [-11.845579664618569, -1713.7381754113042, -52.47409739912885],
                [-8.83290713867791, 10.741342747816361, -348.80482226706044],
                6206690.003928819,
                329878058.1182652),
            "left thoracic cardiopulmonary branch B": (
                "left vagus nerve", 2,
                [22198.762145694476, -3181.553031186403, -628.400239786329],
                [886.9465489726451, 756.0670892895126, -89.42734555851769],
                [1.0856990421307273, 39.326727353878596, 343.25743550415484],
                4386798.475452815,
                230441652.0883617)
        }
        groups_count = len(expected_group_info)

        # check all meshes and groups are created and of appropriate size
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        expected_elements_count = 34
        self.assertEqual(expected_elements_count, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(expected_elements_count * 9 + groups_count, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(expected_elements_count * 17 + groups_count * 8, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(expected_elements_count + 1 + 8, nodes.getSize())  # including 6 marker points
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())
        meshDerivative1 = mesh3d.getChartDifferentialoperator(order=1, term=1)
        meshDerivative3 = mesh3d.getChartDifferentialoperator(order=1, term=3)
        TOL = 0.01  # coordinates and derivatives
        STOL = 0.1  # surface area
        VTOL = 1.0  # volume
        MTOL = 1.0E-7  # material coordinate
        one = fieldmodule.createFieldConstant(1.0)
        for group_name in expected_group_info.keys():
            parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3, \
                expected_surface_area, expected_volume = expected_group_info[group_name]
            group = fieldmodule.findFieldByName(group_name).castGroup()
            mesh_group3d = group.getMeshGroup(mesh3d)
            self.assertEqual(expected_elements_count, mesh_group3d.getSize())
            mesh_group2d = group.getMeshGroup(mesh2d)
            expected_face_count = expected_elements_count * 9 + 1
            self.assertEqual(expected_face_count, mesh_group2d.getSize())
            mesh_group1d = group.getMeshGroup(mesh1d)
            expected_line_count = expected_elements_count * 17 + 8
            self.assertEqual(expected_line_count, mesh_group1d.getSize())
            nodeset_group = group.getNodesetGroup(nodes)
            expected_node_count = expected_elements_count + (2 if parent_group_name else 1)
            self.assertEqual(expected_node_count, nodeset_group.getSize())
            if parent_group_name:
                # check first 2 nodes are in parent nodeset group
                parent_group = fieldmodule.findFieldByName(parent_group_name).castGroup()
                parent_nodeset_group = parent_group.getNodesetGroup(nodes)
                nodeiterator = nodeset_group.createNodeiterator()
                for n in range(2):
                    node = nodeiterator.next()
                    self.assertTrue(parent_nodeset_group.containsNode(node))
            element = mesh_group3d.createElementiterator().next()
            self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5]))
            result, start_x = coordinates.evaluateReal(fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
            result, start_d1 = coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
            result, start_d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_d3, expected_start_d3, delta=TOL)
            # note surface area is merely sum of all surfaces including epineurium, box and interior surfaces
            surface_area_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group2d)
            surface_area_field.setNumbersOfPoints(4)
            volume_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group3d)
            volume_field.setNumbersOfPoints(3)
            fieldcache.clearLocation()  # otherwise integrates over element only
            result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            result, volume = volume_field.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(expected_surface_area, surface_area, delta=STOL)
            self.assertAlmostEqual(expected_volume, volume, delta=VTOL)

        # check sampled trunk d3 for orientation and radius fit, all at element centre
        xi_centre = [0.5, 0.5, 0.5]
        # (element_identifier, expected_d3)
        expected_d3_info = [
            (2, [-58.6320893646309, 252.10632092292698, 1160.1294488820856]),
            (4, [-466.62766145492685, 684.2084398271177, 792.1096579340568]),
            (6, [-39.95762551866224, 626.8288519565403, 194.81718412021738]),
            (8, [-3.574742930372935, 202.7684810518603, 665.5762493112433]),
            (10, [-35.304634943338584, -271.8616602272858, 641.3892257020541]),
            (12, [-185.8222316950974, -564.9761535332963, 246.0793858722639]),
            (14, [1.0315516085406102, -474.6578157024839, 117.94936286585241]),
            (16, [-3.5248364119734674, -465.9424087280286, 105.13518670585816])]
        for element_identifier, expected_d3 in expected_d3_info:
            element = mesh3d.findElementByIdentifier(element_identifier)
            self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, xi_centre))
            result, d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, d3, expected_d3, delta=TOL)

        # check surface area of epineurium, length of centroids
        expected_elements_count = 34
        group = fieldmodule.findFieldByName("vagus epineurium").castGroup()
        mesh_group2d = group.getMeshGroup(mesh2d)
        self.assertEqual(expected_elements_count * 4, mesh_group2d.getSize())
        surface_area_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group2d)
        surface_area_field.setNumbersOfPoints(4)
        fieldcache.clearLocation()
        result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(72152126.79268955, surface_area, delta=STOL)
        group = fieldmodule.findFieldByName("vagus centroid").castGroup()
        mesh_group1d = group.getMeshGroup(mesh1d)
        self.assertEqual(expected_elements_count, mesh_group1d.getSize())
        length_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group1d)
        length_field.setNumbersOfPoints(4)
        result, length = length_field.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(75773.76066920695, length, delta=STOL)

        # check all markers are added
        marker_group = fieldmodule.findFieldByName("marker").castGroup()
        marker_nodes = marker_group.getNodesetGroup(nodes)
        self.assertEqual(8, marker_nodes.getSize())
        marker_name_field = fieldmodule.findFieldByName("marker_name").castStoredString()
        marker_location_field = fieldmodule.findFieldByName("marker_location").castStoredMeshLocation()
        expected_marker_info = get_left_vagus_marker_locations_list()
        expected_elements_count = 25
        node_iter = marker_nodes.createNodeiterator()
        node = node_iter.next()
        for expected_marker_name, expected_material_coordinate1 in expected_marker_info.items():
            fieldcache.setNode(node)
            marker_name = marker_name_field.evaluateString(fieldcache)
            self.assertEqual(expected_marker_name, marker_name)
            element, xi = marker_location_field.evaluateMeshLocation(fieldcache, 3)
            material_coordinate1 = (element.getIdentifier() - 1 + xi[0]) / expected_elements_count
            self.assertAlmostEqual(expected_material_coordinate1, material_coordinate1, delta=MTOL)
            node = node_iter.next()

        # check vagus material coordinates
        vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()
        self.assertTrue(vagus_coordinates.isValid())
        # (expected_start_x, expected_start_d1, expected_start_d3, expected_surface_area, expected_volume)
        expected_group_material_info = {
            'left vagus nerve': (
                [0.0, 0.0, 0.0],
                [0.04, 0.0, 0.0],
                [0.0, 0.0, 0.012],
                0.07044881379783888,
                0.00014399999999999916),
            'left superior laryngeal nerve A': (
                [0.13157231647692647, 0.000881946577575204, 0.00046342826229338254],
                [-0.007050424464461594, 0.012375232591217526, 0.011945671613974565],
                [-0.00014529821484897398, -0.004122294839398416, 0.004267844200042748],
                0.002002117608655125,
                2.002583526823736e-06),
            'branch of left superior laryngeal nerve A': (
                [0.11477792076676165, 0.02817015018376029, 0.025803094067409528],
                [-0.014572453104873985, -0.010775661380116038, -0.013815688836866404],
                [-0.0007112995737230954, -0.004342817670527603, 0.004084144597586936],
                0.0015430135555127994,
                1.4763762746144708e-06),
            'left thoracic cardiopulmonary branch A': (
                [0.38136799563066826, -0.00033316105639673993, 1.7775357276363264e-05],
                [0.0045715391201361375, -0.026876114055372165, 0.011301086871669301],
                [-0.00012166893939763446, -0.002302094674267968, -0.005536330847622523],
                0.002076223576695471,
                2.1241213756367916e-06),
            'left thoracic cardiopulmonary branch B': (
                [0.4067480269851232, 0.0010023426903899028, -0.0015549456000039472],
                [0.010640904036129741, 0.011289784094864912, -0.01261655629280709],
                [6.442917430682371e-07, 0.00446818676605633, 0.00400443416074313],
                0.0014410574883490532,
                1.4381927734194566e-06)}
        TOL = 1.0E-7  # coordinates and derivatives
        STOL = 1.0E-9  # surface area
        VTOL = 1.0E-11  # volume
        for group_name in expected_group_info.keys():
            expected_start_x, expected_start_d1, expected_start_d3, expected_surface_area, expected_volume = \
                expected_group_material_info[group_name]
            group = fieldmodule.findFieldByName(group_name).castGroup()
            mesh_group3d = group.getMeshGroup(mesh3d)
            mesh_group2d = group.getMeshGroup(mesh2d)
            element = mesh_group3d.createElementiterator().next()
            self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5]))
            result, start_x = vagus_coordinates.evaluateReal(fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
            result, start_d1 = vagus_coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
            result, start_d3 = vagus_coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, start_d3, expected_start_d3, delta=TOL)
            # note surface area is merely sum of all surfaces including epineurium, box and interior surfaces
            surface_area_field = fieldmodule.createFieldMeshIntegral(one, vagus_coordinates, mesh_group2d)
            surface_area_field.setNumbersOfPoints(4)
            volume_field = fieldmodule.createFieldMeshIntegral(one, vagus_coordinates, mesh_group3d)
            volume_field.setNumbersOfPoints(3)
            fieldcache.clearLocation()  # otherwise integrates over element only
            result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            result, volume = volume_field.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(expected_surface_area, surface_area, delta=STOL)
            self.assertAlmostEqual(expected_volume, volume, delta=VTOL)


if __name__ == "__main__":
    unittest.main()
