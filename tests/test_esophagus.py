import copy
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.esophagus_terms import get_esophagus_term
from scaffoldmaker.meshtypes.meshtype_3d_esophagus1 import MeshType_3d_esophagus1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList


class EsophagusScaffoldTestCase(unittest.TestCase):

    def test_esophagus1(self):
        """
        Test creation of esophagus scaffold.
        """
        scaffold = MeshType_3d_esophagus1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1"])
        options = scaffold.getDefaultOptions("Human 1")
        self.assertEqual(15, len(options))
        self.assertEqual(8, options.get("Number of elements around"))
        self.assertEqual(20, options.get("Number of elements along"))
        self.assertEqual(4, options.get("Number of elements through wall"))
        self.assertEqual(3.2, options.get("Wall thickness"))
        self.assertEqual(0.35, options.get("Mucosa relative thickness"))
        self.assertEqual(0.15, options.get("Submucosa relative thickness"))
        self.assertEqual(0.25, options.get("Circular muscle layer relative thickness"))
        self.assertEqual(0.25, options.get("Longitudinal muscle layer relative thickness"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)
        self.assertEqual(12, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(640, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(2112, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(2312, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(844, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-11.461359765715862, -117.02328007059768, 1124.746516933097], 1.0E-6)
        assertAlmostEqualList(self, maximums, [14.41809738372959, -62.10186570804781, 1403.745932596971], 1.0E-6)

        flatCoordinates = fieldmodule.findFieldByName("flat coordinates").castFiniteElement()
        self.assertTrue(flatCoordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(flatCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [0.0, 0.0, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [31.765221935414147, 287.8507353356311, 3.2], 1.0E-6)

        esophagusCoordinates = fieldmodule.findFieldByName("esophagus coordinates").castFiniteElement()
        minimums, maximums = evaluateFieldNodesetRange(esophagusCoordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.65, 0.0, -0.65], 1.0E-4)
        assertAlmostEqualList(self, maximums, [0.65, 15.0, 0.65], 1.0E-4)

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
        self.assertAlmostEqual(surfaceArea, 13994.064785386161, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 35763.04473815047, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "abdominal part of esophagus": 96,
            "cervical part of esophagus": 160,
            "esophagus": 640,
            "esophagus mucosa": 160,
            "esophagus smooth muscle circular layer": 160,
            "esophagus smooth muscle longitudinal layer": 160,
            "submucosa of esophagus": 160,
            "thoracic part of esophagus": 384
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_esophagus_term(name))
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
        self.assertEqual(12, len(annotationGroups))

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements around'] = 4
        options['Refine number of elements along'] = 4
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
        self.assertEqual(14, len(annotationGroups))
#
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(40960, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(125952, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(129056, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(44068, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_esophagus_term(name))
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
        node = findNodeWithName(markerNodes, markerName,
                                "distal point of lower esophageal sphincter serosa on the greater curvature of stomach")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(40509, element.getIdentifier())
        assertAlmostEqualList(self, xi, [0.0, 1.0, 1.0], 1.0E-10)

if __name__ == "__main__":
    unittest.main()
