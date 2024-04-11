import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_3d_uterus2 import MeshType_3d_uterus2
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class StomachScaffoldTestCase(unittest.TestCase):

    def test_uterus1(self):
        """
        Test creation of uterus scaffold.
        """
        scaffold = MeshType_3d_uterus2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Mouse 1"])
        options = scaffold.getDefaultOptions("Human 1")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3, 4-5-6, 3-7-8-11.1, 6-9-10-11.2, 11.3-12-13-14,14-15-16,16-17-18",
                         networkLayoutSettings["Structure"])

        self.assertEqual(14, len(options))
        self.assertEqual(10, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around horns"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(5.0, options.get("Target element density along longest segment"))
        self.assertEqual(True, options.get("Serendipity"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(248, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1005, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1278, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(522, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-9.360152113397383, -0.05, -8.943209696521798], 1.0E-6)
        assertAlmostEqualList(self, maximums, [9.360152113397215, 12.920227960030479, 1.278732071803069], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 264.7433300356289, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 188.23971703799126, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "body of uterus": 130,
            "left uterine tube": 24,
            "uterine cervix": 20,
            "vagina": 50,
            "uterus": 248
            }

        for name in expectedSizes3d:
            term = get_uterus_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 4x4x4 and check result
        # first remove any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if (not annotationGroup.hasMeshGroup(mesh3d)) and \
                    (annotationGroup.hasMeshGroup(mesh2d) or annotationGroup.hasMeshGroup(mesh1d)):
                removeAnnotationGroups.append(annotationGroup)

        for annotationGroup in removeAnnotationGroups:
            annotationGroups.remove(annotationGroup)
        self.assertEqual(17, len(annotationGroups))

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements along'] = 4
        options['Refine number of elements around'] = 4
        options['Refine number of elements through wall'] = 4
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        oldAnnotationGroups = copy.copy(annotationGroups)
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(refineRegion, options, annotationGroups)
        for annotation in annotationGroups:
            if annotation not in oldAnnotationGroups:
                annotationGroup.addSubelements()
        self.assertEqual(29, len(annotationGroups))
#
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(15872, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(51792, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(56016, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(20097, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*64, size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(2, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "junction of left round ligament with uterus")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(3712, element.getIdentifier())
        assertAlmostEqualList(self, xi, [1.0, 1.0, 1.0], 1.0E-04)


if __name__ == "__main__":
    unittest.main()
