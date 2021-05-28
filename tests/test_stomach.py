import copy
import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from scaffoldmaker.meshtypes.meshtype_3d_stomach1 import MeshType_3d_stomach1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace
from testutils import assertAlmostEqualList

class StomachScaffoldTestCase(unittest.TestCase):

    def test_stomach1(self):
        """
        Test creation of stomach scaffold.
        """
        scaffold = MeshType_3d_stomach1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Rat 1" ])
        options = scaffold.getDefaultOptions("Rat 1")
        self.assertEqual(15, len(options))
        self.assertEqual(12, options.get("Number of elements around esophagus"))
        self.assertEqual(14, options.get("Number of elements around duodenum"))
        self.assertEqual(2, options.get("Number of elements between cardia and duodenum"))
        self.assertEqual(1, options.get("Number of elements across cardia"))
        self.assertEqual(0.5, options.get("Wall thickness"))
        self.assertEqual(True, options.get("Limiting ridge"))
        ostiumOptions = options['Gastro-esophagal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        self.assertEqual(1, ostiumSettings.get("Number of vessels"))
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))
        self.assertEqual(4.0, ostiumSettings.get("Ostium diameter"))
        self.assertEqual(3.5, ostiumSettings.get("Ostium length"))
        self.assertEqual(0.5, ostiumSettings.get("Ostium wall thickness"))
        self.assertEqual(1.25, ostiumSettings.get("Vessel inner diameter"))
        self.assertEqual(0.5, ostiumSettings.get("Vessel wall thickness"))
        self.assertEqual(0.0, ostiumSettings.get("Vessel angle 1 degrees"))
        self.assertEqual(0.55, options.get("Gastro-esophagal junction position along factor"))
        self.assertEqual(0.2, options.get("Cardia derivative factor"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)
        self.assertEqual(13, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(158, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(643, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(823, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(342, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-18.535034704449842, -14.989550140413265, -8.726949180639382], 1.0E-6)
        assertAlmostEqualList(self, maximums, [18.39909561216101, 15.187335864125645, 8.727618869278126], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 2452.212186914634, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 1175.367658574881, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "body of stomach" : 28,
            "cardia of stomach" : 12,
            "duodenum" : 14,
            "esophagus" : 24,
            "fundus of stomach" : 38,
            "pyloric antrum" : 28,
            "pylorus": 14,
            "stomach": 158
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_stomach_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

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
        options['Refine number of elements surface'] = 4
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
        self.assertEqual(16, len(annotationGroups))
#
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(10112, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(33232, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(36028, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(12936, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_stomach_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*64, size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getFieldNodeGroup(refinedNodes).getNodesetGroup()
        self.assertEqual(4, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "gastro-esophagal junction on lesser curvature")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(1213, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

if __name__ == "__main__":
    unittest.main()
