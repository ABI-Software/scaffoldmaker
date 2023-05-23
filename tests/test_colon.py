import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_3d_colon1 import MeshType_3d_colon1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, \
    exnode_string_from_nodeset_field_parameters, extractPathParametersFromRegion

from testutils import assertAlmostEqualList


class ColonScaffoldTestCase(unittest.TestCase):

    def test_colon1(self):
        """
        Test creation of colon scaffold.
        """
        parameterSetNames = MeshType_3d_colon1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cattle 1", "Human 1", "Human 2", "Mouse 1", "Mouse 2", "Pig 1",
                                             "Pig 2"])
        centralPathDefaultScaffoldPackages = {
            'Test line': ScaffoldPackage(MeshType_1d_path1, {
                'scaffoldSettings': {
                    'Coordinate dimensions': 3,
                    'D2 derivatives': True,
                    'Length': 1.0,
                    'Number of elements': 1
                },
                'meshEdits': exnode_string_from_nodeset_field_parameters(
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                     Node.VALUE_LABEL_D2_DS1DS2], [
                        (1, [[163.7, -25.2, 12.2], [-21.7, 50.1, -18.1], [0.0, 0.0, 5.0], [0.0, 0.0, 0.5]]),
                        (2, [[117.2, 32.8, -2.6], [-64.3, 34.4, -3.9], [0.0, 0.0, 5.0], [0.0, 0.0, 0.5]])
                    ])
            })
        }
        centralPathOption = centralPathDefaultScaffoldPackages['Test line']
        segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Human 1')
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Segment profile': segmentProfileOption,
            'Number of segments': 3,
            'Start phase': 0.0,
            'Proximal length': 25.0,
            'Transverse length': 25.0,
            'Distal length': 25.0,
            'Proximal inner radius': 20.0,
            'Proximal tenia coli width': 8.0,
            'Proximal-transverse inner radius': 18.0,
            'Proximal-transverse tenia coli width': 6.0,
            'Transverse-distal inner radius': 16.0,
            'Transverse-distal tenia coli width': 5.0,
            'Distal inner radius': 15.0,
            'Distal tenia coli width': 5.0,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        self.assertEqual(19, len(options))
        centralPath = options['Central path']
        segmentProfile = options.get("Segment profile")
        segmentSettings = segmentProfile.getScaffoldSettings()
        self.assertEqual(8, segmentSettings.get("Number of elements around haustrum"))
        self.assertEqual(0.5, segmentSettings.get("Corner inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Haustrum inner radius factor"))
        self.assertEqual(0.5, segmentSettings.get("Segment length end derivative factor"))
        self.assertEqual(3, segmentSettings.get("Number of tenia coli"))
        self.assertEqual(0.6, segmentSettings.get("Tenia coli thickness"))
        self.assertEqual(3, options.get("Number of segments"))
        self.assertEqual(0.0, options.get("Start phase"))
        self.assertEqual(25.0, options.get("Transverse length"))
        self.assertEqual(20.0, options.get("Proximal inner radius"))
        self.assertEqual(6.0, options.get("Proximal-transverse tenia coli width"))
        self.assertEqual(16.0, options.get("Transverse-distal inner radius"))
        self.assertEqual(5.0, options.get("Distal tenia coli width"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx = extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE])[0]
        self.assertEqual(2, len(cx))
        assertAlmostEqualList(self, cx[0], [163.7, -25.2, 12.2], 1.0E-6)
        assertAlmostEqualList(self, cx[1], [117.2, 32.8, -2.6], 1.0E-6)
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
        assertAlmostEqualList(self, minimums, [108.02506479907721, -36.405037279268456, -25.89741158484918], 1.0E-6)
        assertAlmostEqualList(self, maximums, [185.46457506076914, 48.1011574894518, 34.05259862880112], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [186.72988844629867, 77.41781871321301, 2.2], 1.0E-6)

        colonCoordinates = fieldmodule.findFieldByName("colon coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(colonCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.6, 0.0, -0.6], 1.0E-4)
        assertAlmostEqualList(self, maximums, [0.6, 24.0, 0.625], 1.0E-4)

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
        self.assertAlmostEqual(surfaceArea, 14342.540002125375, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 25983.483155342656, delta=1.0E-6)

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
        self.assertAlmostEqual(flatSurfaceArea, 629.4883774904393, delta=1.0E-6)
        result, colonSurfaceArea = colonSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonSurfaceArea, 90.4578820802557, delta=1.0E-6)
        result, colonVolume = colonVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(colonVolume, 8.290058800222006, delta=1.0E-6)


if __name__ == "__main__":
    unittest.main()
