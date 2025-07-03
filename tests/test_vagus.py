import os
import unittest

from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.group import mesh_group_to_identifier_ranges, nodeset_group_to_identifier_ranges
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import findAnnotationGroupByName
from scaffoldmaker.annotation.annotation_utils import annotation_term_id_to_url
from scaffoldmaker.annotation.vagus_terms import vagus_branch_terms, vagus_marker_terms
from scaffoldmaker.meshtypes.meshtype_3d_nerve1 import MeshType_3d_nerve1, get_left_vagus_marker_locations_list
from scaffoldmaker.utils.interpolation import get_curve_from_points, getCubicHermiteCurvesLength
from scaffoldmaker.utils.read_vagus_data import VagusInputData

from testutils import assertAlmostEqualList


here = os.path.abspath(os.path.dirname(__file__))


def reorder_vagus_test_data1(testcase, region):
    """
    Break up and reorder vagus test data to test trunk ordering code.
    """
    fieldmodule = region.getFieldmodule()
    mesh1d = fieldmodule.findMeshByDimension(1)
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    trunk_group = fieldmodule.findFieldByName("left vagus X nerve trunk").castGroup()
    coordinates = fieldmodule.findFieldByName("coordinates")
    IDENTIFIER_OFFSET = 1000
    UNUSED_IDENTIFIER = 100000
    with ChangeManager(fieldmodule):
        for element_identifier in range(1, 101):
            element = mesh1d.findElementByIdentifier(element_identifier)
            eft = element.getElementfieldtemplate(coordinates, -1)
            local_node_count = eft.getNumberOfLocalNodes()
            node_identifiers = [element.getNode(eft, ln).getIdentifier() for ln in range(1, local_node_count + 1)]
            node_identifiers.reverse()
            testcase.assertEqual(RESULT_OK, element.setNodesByIdentifier(eft, node_identifiers))
        for node_identifier in range(1, 51):
            other_node_identifier = 101 - node_identifier
            node = nodes.findNodeByIdentifier(node_identifier)
            other_node = nodes.findNodeByIdentifier(other_node_identifier)
            testcase.assertEqual(RESULT_OK, other_node.setIdentifier(UNUSED_IDENTIFIER))
            testcase.assertEqual(RESULT_OK, node.setIdentifier(other_node_identifier))
            testcase.assertEqual(RESULT_OK, other_node.setIdentifier(node_identifier))
        for element_identifier in range(76, 101):
            other_element_identifier = 151 - element_identifier
            element = mesh1d.findElementByIdentifier(element_identifier)
            other_element = mesh1d.findElementByIdentifier(other_element_identifier)
            testcase.assertEqual(RESULT_OK, other_element.setIdentifier(UNUSED_IDENTIFIER))
            testcase.assertEqual(RESULT_OK, element.setIdentifier(other_element_identifier))
            testcase.assertEqual(RESULT_OK, other_element.setIdentifier(element_identifier))
        # for node_identifier in range(126, 151):
        #     node = nodes.findNodeByIdentifier(node_identifier)
        #     testcase.assertEqual(RESULT_OK, node.setIdentifier(node_identifier + IDENTIFIER_OFFSET))
        for element_identifier in range(101, 104):
            mesh1d.destroyElement(mesh1d.findElementByIdentifier(element_identifier))


class VagusScaffoldTestCase(unittest.TestCase):


    def test_vagus_terms(self):
        """
        Test that all vagus terms are UBERON or ILX, and urls.
        """
        for term in (vagus_branch_terms + vagus_marker_terms):
            term_url = annotation_term_id_to_url(term)[1]
            self.assertTrue((not term_url) or
                            term_url.startswith("http://purl.obolibrary.org/obo/UBERON_") or
                            term_url.startswith("http://uri.interlex.org/base/ilx_"), "Invalid vagus term" + str(term))


    def test_input_vagus_data(self):
        """
        Reading vagus input data
        """

        data_file = os.path.join(here, "resources", "vagus_test_data1.exf")

        # test with original ordered, and reordered vagus data file
        for i in range(2):
            context = Context("Test")
            base_region = context.getDefaultRegion()
            region = base_region.createChild('vagus')
            assert(region.isValid())
            data_region = region.getParent().createChild('data')
            assert(data_region.isValid())
            result = data_region.readFile(data_file)
            assert result == RESULT_OK
            if i == 1:
                reorder_vagus_test_data1(self, data_region)

            data_fieldmodule = data_region.getFieldmodule()
            data_mesh1d = data_fieldmodule.findMeshByDimension(1)
            data_nodes = data_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            trunk_group = data_fieldmodule.findFieldByName("left vagus X nerve trunk").castGroup()
            data_coordinates = data_fieldmodule.findFieldByName("coordinates")
            data_trunk_mesh_group = trunk_group.getMeshGroup(data_mesh1d)
            data_trunk_nodeset_group = trunk_group.getNodesetGroup(data_nodes)
            mesh_ranges = mesh_group_to_identifier_ranges(data_trunk_mesh_group)
            nodeset_ranges = nodeset_group_to_identifier_ranges(data_trunk_nodeset_group)
            if i == 0:
                self.assertEqual([[1, 200]], mesh_ranges)
                self.assertEqual([[1, 201]], nodeset_ranges)
                self.assertEqual(200, data_trunk_mesh_group.getSize())
                expected_element_info = {
                    1: [1, 2],
                    2: [2, 3],
                    51: [51, 52],
                    101: [101, 102],
                    102: [102, 103],
                    103: [103, 104]
                }
            else:
                self.assertEqual([[1, 100], [104, 200]], mesh_ranges)
                # self.assertEqual([[1, 125], [151, 201], [1126, 1150]], nodeset_ranges)
                self.assertEqual(197, data_trunk_mesh_group.getSize())
                expected_element_info = {
                    1: [99, 100],
                    2: [98, 99],
                    51: [101, 1],
                    101: None,
                    102: None,
                    103: None
                }
            for element_id, expected_node_ids in expected_element_info.items():
                element = data_mesh1d.findElementByIdentifier(element_id)
                eft = element.getElementfieldtemplate(data_coordinates, -1)
                if expected_node_ids:
                    node_ids = [element.getNode(eft, n + 1).getIdentifier() for n in range(eft.getNumberOfLocalNodes())]
                    self.assertEqual(expected_node_ids, node_ids)
                else:
                    self.assertFalse(eft.isValid())

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

            # do a simple fit to the trunk data coordinates to check trunk ordering is working
            trunk_data_coordinates = vagus_data.get_trunk_coordinates()
            px = [e[0] for e in trunk_data_coordinates]
            self.assertEqual(201, len(px))
            bx, bd1 = get_curve_from_points(px, number_of_elements=10)
            length = getCubicHermiteCurvesLength(bx, bd1)
            self.assertAlmostEqual(31726.825262197974, length, delta=1.0E-3)

            branch_data = vagus_data.get_branch_data()
            self.assertEqual(len(branch_data), 4)
            self.assertTrue("left superior laryngeal nerve" in branch_data)
            self.assertEqual(len(branch_data["left superior laryngeal nerve"]), 42)
            self.assertTrue("left A branch of superior laryngeal nerve" in branch_data)
            self.assertEqual(len(branch_data["left A branch of superior laryngeal nerve"]), 22)
            left_thoracic_cardiopulmonary_branches = (
                "left A thoracic cardiopulmonary branch of vagus nerve",
                "left B thoracic cardiopulmonary branch of vagus nerve")
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertTrue(branch_name in branch_data)

            branch_parents = vagus_data.get_branch_parent_map()
            self.assertEqual(branch_parents["left superior laryngeal nerve"], trunk_group_name)
            self.assertEqual(branch_parents["left A branch of superior laryngeal nerve"], "left superior laryngeal nerve")
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertEqual(branch_parents[branch_name], trunk_group_name)

            branch_common_groups = vagus_data.get_branch_common_group_map()
            self.assertEqual(len(branch_common_groups), 2)
            self.assertTrue("left branch of superior laryngeal nerve" in branch_common_groups)
            self.assertTrue("left A branch of superior laryngeal nerve" in \
                            branch_common_groups["left branch of superior laryngeal nerve"])
            self.assertTrue("left thoracic cardiopulmonary branch of vagus nerve" in branch_common_groups)
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertTrue(branch_name in branch_common_groups["left thoracic cardiopulmonary branch of vagus nerve"])

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
        self.assertEqual(len(options), 7)
        self.assertEqual(options.get('Base parameter set'), 'Human Left Vagus 1')
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

        # test with original ordered, and reordered vagus data file
        for i in range(2):

            context = Context("Test")
            root_region = context.getDefaultRegion()
            region = root_region.createChild('vagus')
            data_region = root_region.createChild('data')
            data_file = os.path.join(here, "resources", "vagus_test_data1.exf")
            self.assertEqual(data_region.readFile(data_file), RESULT_OK)
            if i == 1:
                reorder_vagus_test_data1(self, data_region)

            # check annotation groups
            annotation_groups = scaffold.generateMesh(region, options)[0]
            self.assertEqual(len(annotation_groups), 20)

            # (term_id, parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3,
            #  expected_surface_area, expected_volume)
            expected_group_info = {
                "left vagus nerve": (
                    "http://uri.interlex.org/base/ilx_0785628", None, 25,
                    [-1269.8048516184547, -6359.977051431916, -69.78642824721726],
                    [2163.657939271601, -1111.9771974322234, 121.45057496461462],
                    [49.68213484225328, 258.2220400479382, 1479.1356481323735],
                    248430668.23162192,
                    32973643021.906002),
                "left superior laryngeal nerve": (
                    "http://uri.interlex.org/base/ilx_0788780", "left vagus nerve", 3,
                    [5923.104657597037, -4450.247919770724, -196.91175665569304],
                    [-1473.665051675918, 858.0807042974036, 37.61890734384076],
                    [29.08457359713975, 24.68196779670643, 576.3537772211523],
                    9788808.200600598,
                    558709236.3511268),
                "left A branch of superior laryngeal nerve": (
                    "http://uri.interlex.org/base/ilx_0795823", "left superior laryngeal nerve", 2,
                    [5105.456364262517, -1456.268405569011, 0.18793093373064806],
                    [-1289.581295107282, 381.4601337342457, 17.4939305617649],
                    [2.990626464088564, -3.646186476940329, 299.9629335060045],
                    4696615.9900110625,
                    236649716.13221297),
                "left A thoracic cardiopulmonary branch of vagus nerve": (
                    "http://uri.interlex.org/base/ilx_0794192", "left vagus nerve", 2,
                    [20637.123232118385, -2947.094130818923, -608.014306886659],
                    [99.3811573594032, -1713.8817535655435, -61.05879554434691],
                    [-8.860965678419234, 11.91104819082625, -348.7579634464091],
                    6201808.834555441,
                    328598403.66705346),
                "left B thoracic cardiopulmonary branch of vagus nerve": (
                    "http://uri.interlex.org/base/ilx_0794193", "left vagus nerve", 1,
                    [22164.37254634065, -3219.4137858081867, -620.4335804280934],
                    [1775.1656728388923, 1620.6261382213854, -217.236772843207],
                    [2.3594134224331356, 43.375834887613564, 342.871791578891],
                    4658643.2842846215,
                    267729615.1145657)
            }
            groups_count = len(expected_group_info)

            # check all meshes and groups are created and of appropriate size
            fieldmodule = region.getFieldmodule()
            fieldcache = fieldmodule.createFieldcache()
            coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
            self.assertTrue(coordinates.isValid())
            self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
            mesh3d = fieldmodule.findMeshByDimension(3)
            expected_elements_count = 33
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
            XTOL = 0.01  # coordinates and derivatives
            LTOL = 0.0001  # length
            STOL = 0.1  # surface area
            VTOL = 1.0  # volume
            MTOL = 1.0E-7  # material coordinate
            trunk_group_name = "left vagus trunk"
            one = fieldmodule.createFieldConstant(1.0)
            for group_name in expected_group_info.keys():
                term_id, parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3, \
                    expected_surface_area, expected_volume = expected_group_info[group_name]
                group = fieldmodule.findFieldByName(group_name).castGroup()
                annotation_group = findAnnotationGroupByName(annotation_groups, group_name)
                self.assertEqual(term_id, annotation_group.getId())
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
                branch_of_branch = False
                if parent_group_name:
                    # check first 2 nodes are in parent nodeset group
                    parent_group = fieldmodule.findFieldByName(parent_group_name).castGroup()
                    parent_nodeset_group = parent_group.getNodesetGroup(nodes)
                    nodeiterator = nodeset_group.createNodeiterator()
                    for n in range(2):
                        node = nodeiterator.next()
                        self.assertTrue(parent_nodeset_group.containsNode(node))
                    branch_of_branch = parent_group_name != trunk_group_name
                element = mesh_group3d.createElementiterator().next()
                self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5]))
                result, start_x = coordinates.evaluateReal(fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d1 = coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                TOL = 10.0 * XTOL if branch_of_branch else XTOL
                assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
                assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
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
                self.assertAlmostEqual(expected_surface_area, surface_area, delta=350.0 if branch_of_branch else STOL)
                self.assertAlmostEqual(expected_volume, volume, delta=20000.0 if branch_of_branch else VTOL)

            # check sampled trunk d3 for orientation and radius fit, all at element centre
            xi_centre = [0.5, 0.5, 0.5]
            # (element_identifier, expected_d3)
            expected_d3_info = [
                (2, [-58.61783307229305, 252.12576421778604, 1160.1311747920795]),
                (4, [-466.6104180485502, 684.1844267896083, 792.1344245659976]),
                (6, [-39.96816960556015, 626.8518820786068, 194.7956076398493]),
                (8, [-3.5724520665194177, 202.78915581284485, 665.5701348977416]),
                (10, [-35.32843904524995, -271.90830737446174, 641.3643386269923]),
                (12, [-185.7857471471036, -564.9499544682276, 246.03766692690448]),
                (14, [1.024675995239022, -474.64623017953363, 117.92958351592786]),
                (16, [-3.5041316672463836, -465.94920208179235, 105.11830084069737])]

            for element_identifier, expected_d3 in expected_d3_info:
                element = mesh3d.findElementByIdentifier(element_identifier)
                self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, xi_centre))
                result, d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                assertAlmostEqualList(self, d3, expected_d3, delta=XTOL)

            # check volume of trunk, surface area of epineurium, length of centroids, coordinates and straight coordinates
            straight_coordinates = fieldmodule.findFieldByName("straight coordinates").castFiniteElement()
            self.assertTrue(straight_coordinates.isValid())
            STOL = 100.0
            for coordinate_field in (coordinates, straight_coordinates):
                group = fieldmodule.findFieldByName("left vagus nerve").castGroup()
                mesh_group3d = group.getMeshGroup(mesh3d)
                volume_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group3d)
                volume_field.setNumbersOfPoints(3)
                fieldcache.clearLocation()
                result, volume = volume_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                expected_volume = 32973643021.906002 if (coordinate_field is coordinates) else 33282940849.74868
                self.assertAlmostEqual(expected_volume, volume, delta=STOL)
                expected_elements_count = 33
                group = fieldmodule.findFieldByName("vagus epineurium").castGroup()
                mesh_group2d = group.getMeshGroup(mesh2d)
                self.assertEqual(expected_elements_count * 4, mesh_group2d.getSize())
                surface_area_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group2d)
                surface_area_field.setNumbersOfPoints(4)
                fieldcache.clearLocation()
                result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                expected_surface_area = 72273975.5732966 if (coordinate_field is coordinates) else 72581935.90689336
                self.assertAlmostEqual(expected_surface_area, surface_area, delta=STOL)
                group = fieldmodule.findFieldByName("vagus centroid").castGroup()
                mesh_group1d = group.getMeshGroup(mesh1d)
                self.assertEqual(expected_elements_count, mesh_group1d.getSize())
                length_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group1d)
                length_field.setNumbersOfPoints(4)
                result, length = length_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                self.assertAlmostEqual(75894.09718530288, length, delta=LTOL)

            # check all markers are added
            marker_group = fieldmodule.findFieldByName("marker").castGroup()
            marker_nodes = marker_group.getNodesetGroup(nodes)
            self.assertEqual(8, marker_nodes.getSize())
            marker_name_field = fieldmodule.findFieldByName("marker_name").castStoredString()
            marker_location_field = fieldmodule.findFieldByName("marker_location").castStoredMeshLocation()
            expected_marker_info = get_left_vagus_marker_locations_list()
            expected_elements_count = 25
            expected_marker_node_identifier = 10001
            node_iter = marker_nodes.createNodeiterator()
            node = node_iter.next()
            for expected_marker_name, expected_material_coordinate3 in expected_marker_info.items():
                self.assertEqual(expected_marker_node_identifier, node.getIdentifier())
                fieldcache.setNode(node)
                marker_name = marker_name_field.evaluateString(fieldcache)
                self.assertEqual(expected_marker_name, marker_name)
                element, xi = marker_location_field.evaluateMeshLocation(fieldcache, 3)
                material_coordinate3 = (element.getIdentifier() - 1 + xi[0]) / expected_elements_count
                self.assertAlmostEqual(expected_material_coordinate3, material_coordinate3, delta=MTOL)
                expected_marker_node_identifier += 1
                node = node_iter.next()

            # check vagus material coordinates
            vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()
            self.assertTrue(vagus_coordinates.isValid())
            # (expected_start_x, expected_start_d1, expected_start_d3, expected_surface_area, expected_volume)
            expected_group_material_info = {
                'left vagus nerve': (
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.04],
                    [0.0, 0.012, 0.0],
                    0.07044881379783888,
                    0.00014399999999999916),
                'left superior laryngeal nerve': (
                    [0.0004810143452341199, 0.00016215537359726758, 0.1315566022306105],
                    [0.012918159918245761, 0.011858038889900713, -0.005993045531008481],
                    [-0.004009185375038346, 0.004374444695059173, -0.00014749081697862376],
                    0.0019961825738478577,
                    1.9965460129518142e-06),
                'left A branch of superior laryngeal nerve': (
                    [0.029038466228280567, 0.02563897682357003, 0.11742199136150967],
                    [-0.009398043228514652, -0.013415035881397646, -0.017342847114134856],
                    [-0.004251880840935028, 0.004133184586737899, -0.0009537255926423682],
                    0.0016020806733303852,
                    1.5300531912537515e-06),
                'left A thoracic cardiopulmonary branch of vagus nerve': (
                    [-0.00023305107492456, -5.219544526299368e-06, 0.381097023646834],
                    [-0.02667084814216207, 0.011012987852591622, 0.006327322474065759],
                    [-0.002266194165367729, -0.005551192787385297, -0.0001201586042652858],
                    0.002078138926255945,
                    2.1261681901588285e-06),
                'left B thoracic cardiopulmonary branch of vagus nerve': (
                    [0.0005132914610188526, -0.0009107494333354117, 0.4063642319952978],
                    [0.02359392799517729, -0.02679536894825348, 0.020864764064142793],
                    [0.004500416564857745, 0.0039681076697987636, 1.0469921152278516e-08],
                    0.0014496559346632183,
                    1.4855508165393375e-06)}
            XTOL = 1.0E-7  # coordinates and derivatives
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
                result, start_d1 = vagus_coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d3 = vagus_coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                branch_of_branch = (group_name == 'left A branch of superior laryngeal nerve')
                TOL = (10.0 * XTOL) if branch_of_branch else XTOL
                assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
                assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
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

            # test combined groups
            branch_common_map = {
                'left branch of superior laryngeal nerve': [
                    'left A branch of superior laryngeal nerve'],
                'left thoracic cardiopulmonary branch of vagus nerve': [
                    'left A thoracic cardiopulmonary branch of vagus nerve',
                    'left B thoracic cardiopulmonary branch of vagus nerve']}
            for common_branch_name, variant_branch_name_list in branch_common_map.items():
                common_mesh_group = fieldmodule.findFieldByName(common_branch_name).castGroup().getMeshGroup(mesh3d)
                sum_variant_group_sizes = 0
                for variant_branch_name in variant_branch_name_list:
                    variant_mesh_group = fieldmodule.findFieldByName(variant_branch_name).castGroup().getMeshGroup(mesh3d)
                    sum_variant_group_sizes += variant_mesh_group.getSize()
                    elem_iter = variant_mesh_group.createElementiterator()
                    element = elem_iter.next()
                    while element.isValid():
                        self.assertTrue(common_mesh_group.containsElement(element))
                        element = elem_iter.next()
                self.assertEqual(common_mesh_group.getSize(), sum_variant_group_sizes)

            # test sizes of cervical and thoracic parts
            expected_section_info = {
                "left cervical vagus nerve": ("http://uri.interlex.org/base/ilx_0794142", 8),
                "left thoracic vagus nerve": ("http://uri.interlex.org/base/ilx_0787543", 17)}
            for section_name, info in expected_section_info.items():
                expected_id, expected_mesh_size = info
                annotation_group = findAnnotationGroupByName(annotation_groups, section_name)
                self.assertEqual(expected_id, annotation_group.getId())
                self.assertEqual(expected_mesh_size, annotation_group.getMeshGroup(mesh3d).getSize())


if __name__ == "__main__":
    unittest.main()
