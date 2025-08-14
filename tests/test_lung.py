import copy
import math
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term, lung_terms
from scaffoldmaker.meshtypes.meshtype_3d_lung1 import MeshType_3d_lung1
from scaffoldmaker.meshtypes.meshtype_3d_lung2 import MeshType_3d_lung2
from scaffoldmaker.meshtypes.meshtype_3d_lung3 import MeshType_3d_lung3
from scaffoldmaker.utils.meshrefinement import MeshRefinement

from testutils import assertAlmostEqualList, check_annotation_term_ids


class LungScaffoldTestCase(unittest.TestCase):

    def test_lung_annotations(self):
        """
        Test nomenclature of the lung terms. 
        """
        for term_ids in lung_terms:
            self.assertTrue(check_annotation_term_ids(term_ids), "Invalid primary term id or order not UBERON < ILX < FMA for lung annotation term ids " + str(term_ids)) 
            
    def test_lung1(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_lung1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Pig 1", "Rat 1" ])
        options = scaffold.getDefaultOptions(["Human 1"])
        self.assertEqual(3, len(options))
        self.assertFalse(scaffold.checkOptions(options))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(20, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(88, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(292, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(334, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(141, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [ 14.153, 62.444, -353.40564 ], 1.0E-6)
        assertAlmostEqualList(self, maximums, [ 297.502, 245.3003, -23.486 ], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            upperRightFissureGroup = fieldmodule.findFieldByName('upper lobe of right lung').castGroup()
            self.assertTrue(upperRightFissureGroup.isValid())
            upperRightFissureMeshGroup = upperRightFissureGroup.getMeshGroup(mesh2d)
            self.assertTrue(upperRightFissureMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, upperRightFissureMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 138894.27893716586, delta=1.0E-2)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 7225489.982550352, delta=1.0E-2)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 88,
            "left lung": 44,
            "right lung": 44,
            "upper lobe of left lung": 24,
            "lower lobe of left lung": 20,
            "upper lobe of right lung": 16,
            "middle lobe of right lung": 8,
            "lower lobe of right lung": 20,
            # "right lung accessory lobe": 0
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "horizontal fissure of right lung": 4,
            "oblique fissure of left lung": 8,
            "oblique fissure of right lung": 8
            }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getNodesetGroup(nodes)
        self.assertEqual(9, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(37, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 0.0, 1.0 ], 1.0E-10)

        # refine 4 and check result
        # first remove any face (but not point) annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if (not annotationGroup.hasMeshGroup(mesh3d)) and (annotationGroup.hasMeshGroup(mesh2d) or annotationGroup.hasMeshGroup(mesh1d)):
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            annotationGroups.remove(annotationGroup)
        self.assertEqual(17, len(annotationGroups))

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine'] = True
        options['Refine number of elements'] = 4
        refineNumberOfElements = options['Refine number of elements']
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
        self.assertEqual(20, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(5632, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(17344, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(17848, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(6147, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*(refineNumberOfElements**3), size, name)
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name]*(refineNumberOfElements**2), size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(9, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(2353, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 0.0, 1.0 ], 1.0E-10)

    def test_lung2_human(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_lung2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Rat 1", "Pig 1", "Material" ])
        options = scaffold.getDefaultOptions(["Human 1"])
        self.assertEqual(37, len(options))
        self.assertFalse(scaffold.checkOptions(options))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        # Need to do the following manually to save originalAnnotationGroups which has some temporary groups
        fieldmodule = region.getFieldmodule()
        # annotationGroups = scaffold.generateMesh(region, options)
        with ChangeManager(fieldmodule):
            annotationGroups = scaffold.generateBaseMesh(region, options)[0]
            fieldmodule.defineAllFaces()
            originalAnnotationGroups = copy.copy(annotationGroups)
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
            scaffold.defineFaceAnnotations(region, options, annotationGroups)
            for annotationGroup in annotationGroups:
                if annotationGroup not in originalAnnotationGroups:
                    annotationGroup.addSubelements()

        self.assertEqual(54, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(88, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(292, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(334, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(139, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        material_coordinates = fieldmodule.findFieldByName("lung coordinates").castFiniteElement()
        self.assertTrue(material_coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5604046150578359, -0.46588462707467704, -0.20848015499159522], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.5604046150578359, 0.3440517958126045, 0.9604261163849366], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            upperRightFissureGroup = fieldmodule.findFieldByName('upper lobe of right lung').castGroup()
            self.assertTrue(upperRightFissureGroup.isValid())
            upperRightFissureMeshGroup = upperRightFissureGroup.getMeshGroup(mesh2d)
            self.assertTrue(upperRightFissureMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, upperRightFissureMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            left_lung_group = getAnnotationGroupForTerm(annotationGroups, get_lung_term("left lung"))
            left_lung_mesh_group = left_lung_group.getMeshGroup(mesh3d)
            materialVolumeField = fieldmodule.createFieldMeshIntegral(one, material_coordinates, left_lung_mesh_group)
            materialVolumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 1.838244693520495, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume,  0.3859394886952645, delta=1.0E-6)
        result, material_volume = materialVolumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        halfEllipsoidVolume = 0.5 * (4.0 / 3.0) * math.pi * 0.25 * 0.5 * 1.0  # 3 axis lengths
        self.assertAlmostEqual(halfEllipsoidVolume, 0.2617993877991494, delta=1.0E-12)
        # not quite the same as the above due to sharpening on the ventral edge
        self.assertAlmostEqual(material_volume, 0.2597062752904966, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 88,
            "left lung": 44,
            "right lung": 44,
            "upper lobe of left lung": 24,
            "lower lobe of left lung": 20,
            "upper lobe of right lung": 16,
            "middle lobe of right lung": 8,
            "lower lobe of right lung": 20,
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "base of left lung surface": 12,
            "base of lower lobe of left lung surface": 8,
            "base of upper lobe of left lung surface": 4,
            "base of right lung surface": 12,
            "base of lower lobe of right lung surface": 8,
            "base of middle lobe of right lung surface": 4,
            "horizontal fissure of right lung": 4,
            "horizontal fissure of middle lobe of right lung": 4,
            "horizontal fissure of upper lobe of right lung": 4,
            'lateral surface of lower lobe of left lung': 10,
            'lateral surface of upper lobe of left lung': 12,
            'lateral surface of lower lobe of right lung': 10,
            'lateral surface of middle lobe of right lung': 4,
            'lateral surface of upper lobe of right lung': 8,
            'medial surface of lower lobe of left lung': 10,
            'medial surface of upper lobe of left lung': 12,
            'medial surface of lower lobe of right lung': 10,
            'medial surface of middle lobe of right lung': 4,
            'medial surface of upper lobe of right lung': 8,
            "oblique fissure of left lung": 8,
            "oblique fissure of right lung": 8,
            "oblique fissure of lower lobe of left lung": 8,
            "oblique fissure of upper lobe of left lung": 8,
            "oblique fissure of lower lobe of right lung": 8,
            # "oblique fissure of middle lobe of right lung": 4,
            # "oblique fissure of upper lobe of right lung": 4,
            "left lung surface": 56,
            "lower lobe of left lung surface": 36,
            "upper lobe of left lung surface": 36,
            "right lung surface": 56,
            "lower lobe of right lung surface": 36,
            "middle lobe of right lung surface": 20,
            "upper lobe of right lung surface": 24,
        }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)
        expectedSizes1d = {
            "anterior border of left lung": 5,
            "anterior border of right lung": 5
            }
        for name in expectedSizes1d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh1d).getSize()
            self.assertEqual(expectedSizes1d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getNodesetGroup(nodes)
        self.assertEqual(7, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(40, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

        # refine 2x2x2 and check result
        # need to use original annotation groups to get temporaries
        annotationGroups = originalAnnotationGroups
        self.assertEqual(19, len(annotationGroups))  # including 4 temporary groups

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine'] = True
        options['Refine number of elements'] = 2
        refineNumberOfElements = options['Refine number of elements']
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
        self.assertEqual(54, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(704, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(2224, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(2364, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(853, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*(refineNumberOfElements**3), size, name)
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name]*(refineNumberOfElements**2), size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(7, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(319, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

    def test_lung2_mouse(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_lung2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Rat 1", "Pig 1", "Material" ])
        options = scaffold.getDefaultOptions(["Mouse 1"])
        self.assertEqual(37, len(options))
        self.assertFalse(scaffold.checkOptions(options))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        # Need to do the following manually to save originalAnnotationGroups which has some temporary groups
        fieldmodule = region.getFieldmodule()
        # annotationGroups = scaffold.generateMesh(region, options)
        with ChangeManager(fieldmodule):
            annotationGroups = scaffold.generateBaseMesh(region, options)[0]
            fieldmodule.defineAllFaces()
            originalAnnotationGroups = copy.copy(annotationGroups)
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
            scaffold.defineFaceAnnotations(region, options, annotationGroups)
            for annotationGroup in annotationGroups:
                if annotationGroup not in originalAnnotationGroups:
                    annotationGroup.addSubelements()

        self.assertEqual(54, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(108, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(366, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(429, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(187, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5965062723149795, -0.9035466860119679, -0.6098596059713717], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.6642156880908214, 0.009214327892980281, 0.669464028759718], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            upperRightFissureGroup = fieldmodule.findFieldByName('upper lobe of right lung').castGroup()
            self.assertTrue(upperRightFissureGroup.isValid())
            upperRightFissureMeshGroup = upperRightFissureGroup.getMeshGroup(mesh2d)
            self.assertTrue(upperRightFissureMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, upperRightFissureMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 2.077889388115966, delta=1.0E-2)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume,  0.3799935491917908, delta=1.0E-2)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 108,
            "left lung": 44,
            "right lung": 44,
            "upper lobe of right lung": 16,
            "middle lobe of right lung": 8,
            "lower lobe of right lung": 20,
            "right lung accessory lobe": 20
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "base of left lung surface": 12,
            "base of right lung surface": 12,
            "base of lower lobe of right lung surface": 8,
            "base of middle lobe of right lung surface": 4,
            "horizontal fissure of right lung": 4,
            "horizontal fissure of middle lobe of right lung": 4,
            "horizontal fissure of upper lobe of right lung": 4,
            'lateral surface of lower lobe of right lung': 10,
            'lateral surface of middle lobe of right lung': 4,
            'lateral surface of upper lobe of right lung': 8,
            'medial surface of lower lobe of right lung': 10,
            'medial surface of middle lobe of right lung': 4,
            'medial surface of upper lobe of right lung': 8,
            "oblique fissure of right lung": 8,
            "oblique fissure of lower lobe of right lung": 8,
            # "oblique fissure of middle lobe of right lung": 4,
            # "oblique fissure of upper lobe of right lung": 4,
            "left lung surface": 56,
            "right lung surface": 56,
            "lower lobe of right lung surface": 36,
            "middle lobe of right lung surface": 20,
            "upper lobe of right lung surface": 24,
        }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)
        expectedSizes1d = {
            "anterior border of left lung": 5,
            "anterior border of right lung": 5
            }
        for name in expectedSizes1d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh1d).getSize()
            self.assertEqual(expectedSizes1d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getNodesetGroup(nodes)
        self.assertEqual(13, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(40, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

        # refine 2x2x2 and check result
        # need to use original annotation groups to get temporaries
        annotationGroups = originalAnnotationGroups
        self.assertEqual(24, len(annotationGroups))  # including 4 temporary groups

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine'] = True
        options['Refine number of elements'] = 2
        refineNumberOfElements = options['Refine number of elements']
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
        self.assertEqual(54, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(864, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(2760, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(2970, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1090, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name]*(refineNumberOfElements**3), size, name)
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name]*(refineNumberOfElements**2), size, name)
        # for name in expectedSizes1d:
        #     group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
        #     size = group.getMeshGroup(mesh1d).getSize()
        #     self.assertEqual(expectedSizes1d[name]*refineNumberOfElements, size, name)


        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(13, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(319, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

    def test_lung2_human_openFissures(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_lung2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Rat 1", "Pig 1", "Material" ])
        options = scaffold.getDefaultOptions(["Human 1"])
        options['Open fissures'] = True
        self.assertEqual(37, len(options))
        self.assertFalse(scaffold.checkOptions(options))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(54, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(88, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(312, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(384, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(172, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5604046150578359, -0.46588462707467704, -0.20848015499159522], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.5604046150578359, 0.3440517958126045, 0.9604261163849366], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            upperRightFissureGroup = fieldmodule.findFieldByName('upper lobe of right lung').castGroup()
            self.assertTrue(upperRightFissureGroup.isValid())
            upperRightFissureMeshGroup = upperRightFissureGroup.getMeshGroup(mesh2d)
            self.assertTrue(upperRightFissureMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, upperRightFissureMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 1.8382446935204972, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume,  0.3859394886952645, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 88,
            "left lung": 44,
            "right lung": 44,
            "upper lobe of left lung": 24,
            "lower lobe of left lung": 20,
            "upper lobe of right lung": 16,
            "middle lobe of right lung": 8,
            "lower lobe of right lung": 20,
            # "right lung accessory lobe": 0
        }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "base of left lung surface": 12,
            "base of lower lobe of left lung surface": 8,
            "base of upper lobe of left lung surface": 4,
            "base of right lung surface": 12,
            "base of lower lobe of right lung surface": 8,
            "base of middle lobe of right lung surface": 4,
            "horizontal fissure of right lung": 8,
            "horizontal fissure of middle lobe of right lung": 4,
            "horizontal fissure of upper lobe of right lung": 4,
            'lateral surface of lower lobe of left lung': 10,
            'lateral surface of upper lobe of left lung': 12,
            'lateral surface of lower lobe of right lung': 10,
            'lateral surface of middle lobe of right lung': 4,
            'lateral surface of upper lobe of right lung': 8,
            'medial surface of lower lobe of left lung': 10,
            'medial surface of upper lobe of left lung': 12,
            'medial surface of lower lobe of right lung': 10,
            'medial surface of middle lobe of right lung': 4,
            'medial surface of upper lobe of right lung': 8,
            "oblique fissure of left lung": 16,
            "oblique fissure of right lung": 16,
            "oblique fissure of lower lobe of left lung": 8,
            "oblique fissure of upper lobe of left lung": 8,
            "oblique fissure of lower lobe of right lung": 8,
            # "oblique fissure of middle lobe of right lung": 4,
            # "oblique fissure of upper lobe of right lung": 4,
            "left lung surface": 56,
            "lower lobe of left lung surface": 36,
            "upper lobe of left lung surface": 36,
            "right lung surface": 56,
            "lower lobe of right lung surface": 36,
            "middle lobe of right lung surface": 20,
            "upper lobe of right lung surface": 24,
        }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)
        expectedSizes1d = {
            "anterior border of left lung": 5,
            "anterior border of right lung": 5
        }
        for name in expectedSizes1d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh1d).getSize()
            self.assertEqual(expectedSizes1d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getNodesetGroup(nodes)
        self.assertEqual(7, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(40, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

    def test_lung2_mouse_openFissures(self):
        """
        Test creation of heart scaffold.
        """
        scaffold = MeshType_3d_lung2
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, [ "Default", "Human 1", "Mouse 1", "Rat 1", "Pig 1", "Material" ])
        options = scaffold.getDefaultOptions(["Mouse 1"])
        options['Open fissures'] = True
        self.assertEqual(37, len(options))
        self.assertFalse(scaffold.checkOptions(options))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(54, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(108, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(378, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(459, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(207, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5965062723149795, -0.9035466860119679, -0.6098596059713717], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.6642156880908214, 0.009214327892980281, 0.669464028759718], 1.0E-6)
        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            upperRightFissureGroup = fieldmodule.findFieldByName('upper lobe of right lung').castGroup()
            self.assertTrue(upperRightFissureGroup.isValid())
            upperRightFissureMeshGroup = upperRightFissureGroup.getMeshGroup(mesh2d)
            self.assertTrue(upperRightFissureMeshGroup.isValid())
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, upperRightFissureMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 2.077889388115967, delta=1.0E-2)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume,  0.3799935491917908, delta=1.0E-2)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 108,
            "left lung": 44,
            "right lung": 44,
            "upper lobe of right lung": 16,
            "middle lobe of right lung": 8,
            "lower lobe of right lung": 20,
            "right lung accessory lobe": 20
            }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)
        expectedSizes2d = {
            "base of left lung surface": 12,
            "base of right lung surface": 12,
            "base of lower lobe of right lung surface": 8,
            "base of middle lobe of right lung surface": 4,
            "horizontal fissure of right lung": 8,
            "horizontal fissure of middle lobe of right lung": 4,
            "horizontal fissure of upper lobe of right lung": 4,
            'lateral surface of lower lobe of right lung': 10,
            'lateral surface of middle lobe of right lung': 4,
            'lateral surface of upper lobe of right lung': 8,
            'medial surface of lower lobe of right lung': 10,
            'medial surface of middle lobe of right lung': 4,
            'medial surface of upper lobe of right lung': 8,
            "oblique fissure of right lung": 16,
            "oblique fissure of lower lobe of right lung": 8,
            # "oblique fissure of middle lobe of right lung": 4,
            # "oblique fissure of upper lobe of right lung": 4,
            "left lung surface": 56,
            "right lung surface": 56,
            "lower lobe of right lung surface": 36,
            "middle lobe of right lung surface": 20,
            "upper lobe of right lung surface": 24,
        }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)
        expectedSizes1d = {
            "anterior border of left lung": 5,
            "anterior border of right lung": 5
            }
        for name in expectedSizes1d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh1d).getSize()
            self.assertEqual(expectedSizes1d[name], size, name)

        # test finding a marker in scaffold
        markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        markerNodes = markerGroup.getNodesetGroup(nodes)
        self.assertEqual(13, markerNodes.getSize())
        markerName = fieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = fieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        # test apex marker point
        cache = fieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertEqual(40, element.getIdentifier())
        assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

    def test_lung3_human(self):
        """
        Test creation of human lung scaffold using lung3 mesh.
        """
        scaffold = MeshType_3d_lung3
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames,
                         ["Default", "Human 1 Coarse", "Human 1 Medium", "Human 1 Fine",
                          "Ellipsoid Coarse", "Ellipsoid Medium", "Ellipsoid Fine"])
        options = scaffold.getDefaultOptions("Human 1 Coarse")
        self.assertEqual(20, len(options))
        self.assertFalse(scaffold.checkOptions(options))
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        # Need to do the following manually to save originalAnnotationGroups which has some temporary groups
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            annotationGroups = scaffold.generateBaseMesh(region, options)[0]
            fieldmodule.defineAllFaces()
            originalAnnotationGroups = copy.copy(annotationGroups)
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()
            scaffold.defineFaceAnnotations(region, options, annotationGroups)
            for annotationGroup in annotationGroups:
                if annotationGroup not in originalAnnotationGroups:
                    annotationGroup.addSubelements()
        self.assertEqual(63, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(192, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(640, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(728, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(282, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        tol = 1.0E-6
        # check coordinates range, outside surface area and volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        # material_coordinates = fieldmodule.findFieldByName("lung coordinates").castFiniteElement()
        # self.assertTrue(material_coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5150077744090888, -0.39797018819161917, 0.035835394723902014], tol)
        assertAlmostEqualList(self, maximums, [0.5150077744090888, 0.3290391566888997, 0.9617023416663923], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            leftLungGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("left lung"))
            leftLungSurfaceMeshGroup = leftLungGroup.getMeshGroup(mesh2d)
            self.assertTrue(leftLungSurfaceMeshGroup.isValid())

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, leftLungSurfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)

        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 4.887248605193398, delta=tol)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 0.301840881252553, delta=tol)

        # check some annotationGroups:
        expectedSizes3d = {
            "lung": 192,
            "left lung": 96,
            "right lung": 96,
            "upper lobe of left lung": 48,
            "lower lobe of left lung": 48,
            "upper lobe of right lung": 24,
            "middle lobe of right lung": 24,
            "lower lobe of right lung": 48,
        }
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        expectedSizes2d = {
            "base of left lung surface": 16,
            "base of lower lobe of left lung surface": 8,
            "base of upper lobe of left lung surface": 8,
            "base of right lung surface": 16,
            "base of lower lobe of right lung surface": 8,
            "base of middle lobe of right lung surface": 8,
            "horizontal fissure of right lung": 10,
            "horizontal fissure of middle lobe of right lung": 10,
            "horizontal fissure of upper lobe of right lung": 10,
            'lateral surface of lower lobe of left lung': 16,
            'lateral surface of upper lobe of left lung': 16,
            'lateral surface of lower lobe of right lung': 16,
            'lateral surface of middle lobe of right lung': 8,
            'lateral surface of upper lobe of right lung': 8,
            'medial surface of lower lobe of left lung': 16,
            'medial surface of upper lobe of left lung': 16,
            'medial surface of lower lobe of right lung': 16,
            'medial surface of middle lobe of right lung': 8,
            'medial surface of upper lobe of right lung': 8,
            "oblique fissure of left lung": 20,
            "oblique fissure of right lung": 20,
            "oblique fissure of lower lobe of left lung": 20,
            "oblique fissure of upper lobe of left lung": 20,
            "oblique fissure of lower lobe of right lung": 20,
            "oblique fissure of middle lobe of right lung": 10,
            "oblique fissure of upper lobe of right lung": 10,
            "left lung surface": 64,
            "lower lobe of left lung surface": 52,
            "upper lobe of left lung surface": 52,
            "right lung surface": 64,
            "lower lobe of right lung surface": 52,
            "middle lobe of right lung surface": 36,
            "upper lobe of right lung surface": 36,
        }
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)

        expectedSizes1d = {
            "antero-posterior edge of upper lobe of left lung": 8,
            "antero-posterior edge of upper lobe of right lung": 4
        }
        for name in expectedSizes1d:
            group = getAnnotationGroupForTerm(annotationGroups, [name, ""])
            size = group.getMeshGroup(mesh1d).getSize()
            self.assertEqual(expectedSizes1d[name], size, name)

        # # test finding a marker in scaffold
        # markerGroup = fieldmodule.findFieldByName("marker").castGroup()
        # markerNodes = markerGroup.getNodesetGroup(nodes)
        # self.assertEqual(7, markerNodes.getSize())
        # markerName = fieldmodule.findFieldByName("marker_name")
        # self.assertTrue(markerName.isValid())
        # markerLocation = fieldmodule.findFieldByName("marker_location")
        # self.assertTrue(markerLocation.isValid())
        # # test apex marker point
        # cache = fieldmodule.createFieldcache()
        # node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        # self.assertTrue(node.isValid())
        # cache.setNode(node)
        # element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        # self.assertEqual(40, element.getIdentifier())
        # assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

        # refine 2x2x2 and check result
        # need to use original annotation groups to get temporaries
        annotationGroups = originalAnnotationGroups
        self.assertEqual(16, len(annotationGroups))  # including temporary groups

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine'] = True
        options['Refine number of elements'] = 2
        refineNumberOfElements = options['Refine number of elements']
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
        self.assertEqual(63, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(1536, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(4864, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(5168, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1842, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name] * (refineNumberOfElements ** 3), size, name)
        for name in expectedSizes2d:
            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(name))
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name] * (refineNumberOfElements ** 2), size, name)

        # # test finding a marker in refined scaffold
        # markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        # refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        # self.assertEqual(7, markerNodes.getSize())
        # markerName = refineFieldmodule.findFieldByName("marker_name")
        # self.assertTrue(markerName.isValid())
        # markerLocation = refineFieldmodule.findFieldByName("marker_location")
        # self.assertTrue(markerLocation.isValid())
        # # test apex marker point
        # cache = refineFieldmodule.createFieldcache()
        # node = findNodeWithName(markerNodes, markerName, "apex of left lung")
        # self.assertTrue(node.isValid())
        # cache.setNode(node)
        # element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        # self.assertEqual(319, element.getIdentifier())
        # assertAlmostEqualList(self, xi, [ 0.0, 1.0, 1.0 ], 1.0E-10)

if __name__ == "__main__":
    unittest.main()
