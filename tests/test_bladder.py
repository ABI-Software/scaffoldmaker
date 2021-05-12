import copy
import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_3d_bladderurethra1 import MeshType_3d_bladderurethra1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace
from testutils import assertAlmostEqualList

class BladderScaffoldTestCase(unittest.TestCase):

    def test_bladderurethra1(self):
        """
        Test creation of bladder scaffold.
        """
        scaffold = MeshType_3d_bladderurethra1
        parameterSetNames = MeshType_3d_bladderurethra1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cat 1", "Rat 1", "Human 1"])
        options = scaffold.getDefaultOptions("Cat 1")
        self.assertEqual(26, len(options))
        self.assertEqual(12, options.get("Number of elements along bladder"))
        self.assertEqual(12, options.get("Number of elements around"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(0.5, options.get("Wall thickness"))
        self.assertEqual(45, options.get("Neck angle degrees"))
        ostiumOptions = options['Ureter']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        self.assertEqual(1, ostiumSettings.get("Number of vessels"))
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))
        self.assertEqual(2.2, ostiumSettings.get("Ostium diameter"))
        self.assertEqual(0.5, ostiumSettings.get("Ostium length"))
        self.assertEqual(0.5, ostiumSettings.get("Ostium wall thickness"))
        self.assertEqual(0.8, ostiumSettings.get("Vessel inner diameter"))
        self.assertEqual(0.25, ostiumSettings.get("Vessel wall thickness"))
        self.assertEqual(0.0, ostiumSettings.get("Vessel angle 1 degrees"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)
        self.assertEqual(13, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(240, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(960, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1201, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(487, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-15.48570141314588, -12.992184072505665, -0.5], 1.0E-6)
        assertAlmostEqualList(self, maximums, [15.485696373577879, 13.837600783918127, 127.68599990411343], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 5335.535739164492, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 2555.81694845211, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "dome of the bladder": 108,
            "neck of urinary bladder": 36,
            "urethra": 96,
            "urinary bladder": 144
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_bladder_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getFieldNodeGroup(nodes).getNodesetGroup()
        self.assertEqual(5, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of urinary bladder")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(1, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 0.0, 0.0 ], 1.0E-10)
        apexGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("apex of urinary bladder"))
        self.assertTrue(apexGroup.getNodesetGroup(nodes).containsNode(node))

        # refine 4x4x4 and check result
        # first remove any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if (not annotationGroup.hasMeshGroup(mesh3d)) and (annotationGroup.hasMeshGroup(mesh2d) or annotationGroup.hasMeshGroup(mesh1d)):
                removeAnnotationGroups.append(annotationGroup)

        for annotationGroup in removeAnnotationGroups:
            annotationGroups.remove(annotationGroup)
        self.assertEqual(13, len(annotationGroups))

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
        self.assertEqual(37, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(15360, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(49920, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(53764, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(19210, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_bladder_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*64, size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getFieldNodeGroup(refinedNodes).getNodesetGroup()
        self.assertEqual(5, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of urinary bladder")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(1, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 0.0, 0.0 ], 1.0E-10)

if __name__ == "__main__":
    unittest.main()
