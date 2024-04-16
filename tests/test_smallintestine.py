import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_3d_smallintestine1 import MeshType_3d_smallintestine1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace, \
    exnode_string_from_nodeset_field_parameters, get_nodeset_path_field_parameters
from testutils import assertAlmostEqualList


class SmallIntestineScaffoldTestCase(unittest.TestCase):

    def test_smallintestine1(self):
        """
        Test creation of small intestine scaffold.
        """
        parameterSetNames = MeshType_3d_smallintestine1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cattle 1", "Human 1", "Mouse 1"])
        testNetworkLayout = ScaffoldPackage(MeshType_1d_network_layout1, {
                'scaffoldSettings': {
                    "Structure": "1-2-3-4"
                },
                'meshEdits': exnode_string_from_nodeset_field_parameters(
                    ['coordinates'],
                    [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3  ], [[
                    (1,  [ [  -2.30, 18.50,  -4.40 ], [ -4.20, -0.80,   3.70 ], [  0.00,  0.60,  0.00 ], [  0.00,  0.11,  0.00 ], [ -0.33,  0.01, -0.50 ], [ 0.00, 0.00, 0.50 ] ] ),
                    (2,  [ [  -8.60, 16.30,  -0.40 ], [ -7.10, -2.70,   1.60 ], [  0.00,  0.73,  0.00 ], [  0.00,  0.14,  0.00 ], [  0.08,  0.09, -0.72 ], [ 0.00, 0.00, 0.50 ] ] ),
                    (3,  [ [ -18.30, 12.60,  -1.50 ], [ -6.40, -1.70,  -3.80 ], [  0.00,  0.90,  0.00 ], [  0.00,  0.13,  0.00 ], [  0.61,  0.04, -0.65 ], [ 0.00, 0.00, 0.50 ] ] ),
                    (4,  [ [ -15.60, 13.70,  -6.10 ], [  7.00,  2.10,  -1.80 ], [  0.00,  1.00,  0.00 ], [  0.00,  0.05,  0.00 ], [  0.50,  0.08,  0.86 ], [ 0.00, 0.00, 0.50 ] ] )]]),
                    
                'userAnnotationGroups': [
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1-3',
                        'name': get_smallintestine_term('small intestine')[0],
                        'ontId': get_smallintestine_term('small intestine')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1',
                        'name': get_smallintestine_term('duodenum')[0],
                        'ontId': get_smallintestine_term('duodenum')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '2',
                        'name': get_smallintestine_term('jejunum')[0],
                        'ontId': get_smallintestine_term('jejunum')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '3',
                        'name': get_smallintestine_term('ileum')[0],
                        'ontId': get_smallintestine_term('ileum')[1]
                    }]
            })

        options = MeshType_3d_smallintestine1.getDefaultOptions("Mouse 1")
        options['Network layout'] = testNetworkLayout
        options['Number of segments'] = 3
        self.assertEqual(16, len(options))
        networkLayout = options['Network layout']
        self.assertEqual(3, options.get("Number of segments"))
        self.assertEqual(8, options.get("Number of elements around"))
        self.assertEqual(3, options.get("Number of elements along segment"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(None, options.get("Duodenum length"))
        self.assertEqual(None, options.get("Jejunum-ileum inner radius"))
        self.assertEqual(0.1, options.get("Wall thickness"))

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
        assertAlmostEqualList(self, cx[0], [-2.3, 18.5, -4.4], 1.0E-6)
        assertAlmostEqualList(self, cx[2], [-18.3, 12.6, -1.5], 1.0E-6)
        del tmpFieldmodule
        del tmpRegion

        annotationGroups = MeshType_3d_smallintestine1.generateBaseMesh(region, options)[0]
        self.assertEqual(4, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(72, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(296, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(384, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(160, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-19.956658197690174, 11.509153871441058, -7.068492932596849], 1.0E-6)
        assertAlmostEqualList(self, maximums, [-1.9033816741586018, 19.093964008729433, 0.6249529551386831], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-1.2566370614359177, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [4.790928796724435, 25.314180157278805, 0.1], 1.0E-6)

        smallintestineCoordinates = fieldmodule.findFieldByName("small intestine coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(smallintestineCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.65, 0.0, -0.65], 1.0E-4)
        assertAlmostEqualList(self, maximums, [0.65, 200.0, 0.65], 1.0E-4)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            faceMeshGroup = createFaceMeshGroupExteriorOnFace(fieldmodule, Element.FACE_TYPE_XI3_1)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, faceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            flatSurfaceAreaField = fieldmodule.createFieldMeshIntegral(one, flatCoordinates, faceMeshGroup)
            flatSurfaceAreaField.setNumbersOfPoints(4)

        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 127.90324495617577, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 12.011251644596864, delta=1.0E-6)
        result, flatSurfaceArea = flatSurfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(flatSurfaceArea, 144.42091907104523, delta=1.0E-3)


if __name__ == "__main__":
    unittest.main()
