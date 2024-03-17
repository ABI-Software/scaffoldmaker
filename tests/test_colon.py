import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_3d_colon1 import MeshType_3d_colon1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, \
    exnode_string_from_nodeset_field_parameters, get_nodeset_path_field_parameters

from testutils import assertAlmostEqualList


class ColonScaffoldTestCase(unittest.TestCase):

    def test_colon1(self):
        """
        Test creation of colon scaffold.
        """
        parameterSetNames = MeshType_3d_colon1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cattle 1", "Human 1", "Human 2", "Human 3", "Mouse 1",
                                             "Mouse 2", "Pig 1"])

        testNetworkLayout = ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[0.00, 0.00, 0.00], [-50.70, 178.20, 0.00], [-37.97, -9.49, -18.98], [-6.86, -11.39, -2.36], [-18.61, -3.98, 39.12], [-14.00, -1.00, -12.00]]),
                    (2, [[-47.40, 188.60, 0.00], [-19.30, 177.10, 0.00], [-35.79, -6.51, -13.01], [11.23, 17.36, 14.31], [-12.66, -3.99, 36.28], [-4.00, 19.00, 22.00]]),
                    (3, [[-4.40, 396.50, 0.00], [206.00, 40.10, 0.00], [-13.89, 27.78, 11.11], [13.54, -1.87, 21.51], [-6.05, -12.50, 29.93], [-6.00, 0.00, 51.00]]),
                    (4, [[130.00, 384.10, 0.00], [130.80, -40.50, 0.00], [-5.35, 4.28, 31.06], [5.83, -8.41, 8.86], [-15.28, -27.78, 2.51], [0.00, 1.00, 24.00]])]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_colon_term('colon')[0],
                    'ontId': get_colon_term('colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_colon_term('ascending colon')[0],
                    'ontId': get_colon_term('ascending colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3',
                    'name': get_colon_term('descending colon')[0],
                    'ontId': get_colon_term('descending colon')[1]
                }]
        })

        segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Human 1')
        options = {
            'Base parameter set': 'Human 1',
            'Network layout': testNetworkLayout,
            'Segment profile': segmentProfileOption,
            'Number of segments': 3,
            'Start phase': 0.0,
            'Proximal tenia coli width': 8.0,
            'Proximal-transverse tenia coli width': 6.0,
            'Transverse-distal tenia coli width': 5.0,
            'Distal tenia coli width': 5.0,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        self.assertEqual(14, len(options))
        networkLayout = options['Network layout']
        segmentProfile = options.get("Segment profile")
        segmentSettings = segmentProfile.getScaffoldSettings()
        self.assertEqual(8, segmentSettings.get("Number of elements around haustrum"))
        self.assertEqual(0.536, segmentSettings.get("Corner outer radius factor"))
        self.assertEqual(0.464, segmentSettings.get("Haustrum outer radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Segment length end derivative factor"))
        self.assertEqual(3, segmentSettings.get("Number of tenia coli"))
        self.assertEqual(0.6, segmentSettings.get("Tenia coli thickness"))
        self.assertEqual(3, options.get("Number of segments"))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(None, options.get("Transverse length"))
        self.assertEqual(None, options.get("Proximal inner radius"))
        self.assertEqual(6.0, options.get("Proximal-transverse tenia coli width"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        tmpRegion = region.createRegion()
        networkLayout.generate(tmpRegion)
        tmpFieldmodule = tmpRegion.getFieldmodule()
        cx = get_nodeset_path_field_parameters(
            tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES),
            tmpFieldmodule.findFieldByName('coordinates'),
            [Node.VALUE_LABEL_VALUE])[0]
        self.assertEqual(4, len(cx))
        assertAlmostEqualList(self, cx[0], [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, cx[1], [ -47.40, 188.60, 0.00 ], 1.0E-6)
        del tmpFieldmodule
        del tmpRegion

        annotationGroups = MeshType_3d_colon1.generateBaseMesh(region, options)[0]
        self.assertEqual(11, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        if annotationGroups is not None:
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(1512, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(4986, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(5463, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1989, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-100.08668914502144, -11.342531039430702, -60.09737764484168], 1.0E-6)
        assertAlmostEqualList(self, maximums, [138.56897619597683, 447.93686761958236, 50.13111447364834], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [377.06179507784395, 561.3383898388232, 2.2], 1.0E-6)

        colonCoordinates = fieldmodule.findFieldByName("colon coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(colonCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.599920500334689, 0.0, -0.6], 1.0E-4)
        assertAlmostEqualList(self, maximums, [0.599920500334689, 24.0, 0.625], 1.0E-4)

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
        self.assertAlmostEqual(surfaceArea, 164285.41543554227, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 294925.4567043401, delta=1.0E-6)

    def test_mousecolon1(self):
        """
        Test creation of mouse colon scaffold.
        """
        options = MeshType_3d_colon1.getDefaultOptions("Mouse 2")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = MeshType_3d_colon1.generateBaseMesh(region, options)[0]
        self.assertEqual(10, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(1600, mesh3d.getSize())
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        colonCoordinates = fieldmodule.findFieldByName("colon coordinates").castFiniteElement()
        self.assertTrue(colonCoordinates.isValid())

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)
            colonSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, colonCoordinates, faceMeshGroup)
            colonSurfaceAreaField.setNumbersOfPoints(4)
            colonVolumeField = fieldmodule.createFieldMeshIntegral(one, colonCoordinates, mesh3d)
            colonVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 650.8919830877059, delta=1.0E-6)
        result, colonSurfaceArea = colonSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonSurfaceArea, 90.45785287687968, delta=1.0E-6)
        result, colonVolume = colonVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonVolume, 8.290050554292705, delta=1.0E-6)


if __name__ == "__main__":
    unittest.main()
