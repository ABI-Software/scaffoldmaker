import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term, uterus_terms
from scaffoldmaker.meshtypes.meshtype_3d_uterus1 import MeshType_3d_uterus1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import createFaceMeshGroupExteriorOnFace

from testutils import assertAlmostEqualList, check_annotation_term_ids


class UterusScaffoldTestCase(unittest.TestCase):

    def test_uterus_annotations(self):
        """
        Test nomenclature of the uterus terms. 
        """
        for term_ids in uterus_terms:
            self.assertTrue(check_annotation_term_ids(term_ids), "Invalid primary term id or order not UBERON < ILX < FMA for uterus annotation term ids " + str(term_ids)) 

    def test_uterus1_human_1shell(self):
        """
        Test creation of human uterus scaffold with 1 shell count.
        """
        scaffold = MeshType_3d_uterus1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human 1', 'Human 2', 'Human Pregnant 1', 'Human Pregnant 2',
                                             'Mouse 1', 'Rat 1'])
        options = scaffold.getDefaultOptions("Human 1")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3-4-5-6-7-8-23.1,9-10-11-12-13-14-15-16-23.2,#17-18-19-20-21-22-23.3,"
                         "23.4-24-25-26-27-28-29,29-30-31,31-32-33-34-35-36-37-38",
                         networkLayoutSettings["Structure"])

        self.assertEqual(15, len(options))
        self.assertEqual(20, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around oviduct/uterine horn"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(True, options.get("Use linear through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(320, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1298, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1653, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(678, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-2.9999999999999964, -14.0, -8.201149866258765], 1.0E-6)
        assertAlmostEqualList(self, maximums, [13.994479106259377, 14.0, 2.9919276924165974], 1.0E-6)

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
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 336.64153452445015, delta=5.0E-2)
        self.assertAlmostEqual(volume, 255.69163498084598, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(45, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 32,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            'left broad ligament of uterus': 12,
            'right broad ligament of uterus': 12,
            "uterus": 240
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 2x2x2 and check result
        # first remove faces/lines and any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getDimension() in [1, 2]:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            if "cervix" in annotationGroup.getName():
                continue
            annotationGroups.remove(annotationGroup)
        self.assertEqual(22, len(annotationGroups))
        # must keep all faces and lines as used for refinement

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(22, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2560, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(9032, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(10418, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(3949, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        expectedSizes3d = {
            "fundus of uterus": 32,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            "uterus": 240
        }
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(4, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "junction of left round ligament with uterus")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertTrue(element.isValid())

    def test_uterus1_human_0shell(self):
        """
        Test creation of human uterus scaffold with 0 shell count.
        """
        scaffold = MeshType_3d_uterus1
        options = scaffold.getDefaultOptions("Human 1")
        options["Number of elements through wall"] = 0

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(0, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(320, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(658, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(341, nodes.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-2.9999999999999964, -14.0, -8.201149866258765], 1.0E-6)
        assertAlmostEqualList(self, maximums, [13.994479106259377, 14.0, 2.9919276924165974], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 335.4880848987993, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(29, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 32,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            'left broad ligament of uterus': 12,
            'right broad ligament of uterus': 12,
            "uterus": 240
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

    def test_uterus2_human_1shell(self):
        """
        Test creation of human 2 uterus scaffold with 1 shell count.
        """
        scaffold = MeshType_3d_uterus1
        parameterSetNames = scaffold.getParameterSetNames()
        options = scaffold.getDefaultOptions("Human 2")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3-4-10.1,5-6-7-8-10.2,(9-10.3,10.4-11-12-13,13-14,14-15-16-17-18",
                         networkLayoutSettings["Structure"])

        self.assertEqual(15, len(options))
        self.assertEqual(20, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around oviduct/uterine horn"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(True, options.get("Use linear through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(340, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1378, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1753, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(718, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-3.15, -14.0, -8.201149866258763], 1.0E-6)
        assertAlmostEqualList(self, maximums, [13.994479106259377, 14.0, 3.000008384729209], 1.0E-6)

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
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 340.0350366489269, delta=5.0E-2)
        self.assertAlmostEqual(volume, 251.3644411904487, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(45, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 52,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            'left broad ligament of uterus': 12,
            'right broad ligament of uterus': 12,
            "uterus": 260
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 2x2x2 and check result
        # first remove faces/lines and any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getDimension() in [1, 2]:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            if "cervix" in annotationGroup.getName():
                continue
            annotationGroups.remove(annotationGroup)
        self.assertEqual(22, len(annotationGroups))
        # must keep all faces and lines as used for refinement

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(22, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2720, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(9592, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(11058, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(4189, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        expectedSizes3d = {
            "fundus of uterus": 52,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            "uterus": 260
        }
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(4, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "junction of left round ligament with uterus")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertTrue(element.isValid())

    def test_uterus1_human_pregnant_1shell(self):
        """
        Test creation of human pregnant 1 uterus scaffold with 1 shell count.
        """
        scaffold = MeshType_3d_uterus1
        options = scaffold.getDefaultOptions("Human Pregnant 1")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3-4-5-6-7-8-9-31.1,10-11-12-13-14-15-16-17-18-31.2,"
                         "#19-20-21-22-23-24-25-26-27-28-29-30-31.3,"
                         "31.4-32-33-34-35-36-37-38-39-40-41-42-43,43-44,44-45-46-47-48",
                         networkLayoutSettings["Structure"])

        self.assertEqual(15, len(options))
        self.assertEqual(20, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around oviduct/uterine horn"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(True, options.get("Use linear through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(320, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1298, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1653, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(678, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-5.999999999999985, -16.0, -8.201149866258765], 1.0E-6)
        assertAlmostEqualList(self, maximums, [20.994479106259377, 16.0, 5.983806305763669], 1.0E-6)

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
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 809.956584291775, delta=5.0E-2)
        self.assertAlmostEqual(volume, 267.60969130550086, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(45, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 32,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            'left broad ligament of uterus': 12,
            'right broad ligament of uterus': 12,
            "uterus": 240
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 2x2x2 and check result
        # first remove faces/lines and any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getDimension() in [1, 2]:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            if "cervix" in annotationGroup.getName():
                continue
            annotationGroups.remove(annotationGroup)
        self.assertEqual(22, len(annotationGroups))
        # must keep all faces and lines as used for refinement

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(22, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2560, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(9032, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(10418, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(3949, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        expectedSizes3d = {
            "fundus of uterus": 32,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            "uterus": 240
        }
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(4, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "junction of left round ligament with uterus")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertTrue(element.isValid())

    def test_uterus2_human_pregnant_1shell(self):
        """
        Test creation of human pregnant 2 uterus scaffold with 1 shell count.
        """
        scaffold = MeshType_3d_uterus1
        options = scaffold.getDefaultOptions("Human Pregnant 2")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3-4-5-6-7-8-9-20.1,10-11-12-13-14-15-16-17-18-20.2,"
                         "(19-20.3,20.4-21-22-23,23-24,24-25-26-27-28",
                         networkLayoutSettings["Structure"])

        self.assertEqual(15, len(options))
        self.assertEqual(20, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around oviduct/uterine horn"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(True, options.get("Use linear through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(17, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(340, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1378, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1753, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(718, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-6.3, -16.0, -8.201149866258765], 1.0E-6)
        assertAlmostEqualList(self, maximums, [20.994479106259377, 16.0, 5.979109874196457], 1.0E-6)

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
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 813.7381551716634, delta=5.0E-2)
        self.assertAlmostEqual(volume, 280.18220494062393, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(45, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 52,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            'left broad ligament of uterus': 12,
            'right broad ligament of uterus': 12,
            "uterus": 260
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 2x2x2 and check result
        # first remove faces/lines and any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getDimension() in [1, 2]:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            if "cervix" in annotationGroup.getName():
                continue
            annotationGroups.remove(annotationGroup)
        self.assertEqual(22, len(annotationGroups))
        # must keep all faces and lines as used for refinement

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(22, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2720, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(9592, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(11058, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(4189, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        expectedSizes3d = {
            "fundus of uterus": 52,
            "body of uterus": 128,
            "left oviduct": 40,
            "vagina": 80,
            "uterus": 260
        }
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(4, markerNodes.getSize())
        markerName = refineFieldmodule.findFieldByName("marker_name")
        self.assertTrue(markerName.isValid())
        markerLocation = refineFieldmodule.findFieldByName("marker_location")
        self.assertTrue(markerLocation.isValid())
        cache = refineFieldmodule.createFieldcache()
        node = findNodeWithName(markerNodes, markerName, "junction of left round ligament with uterus")
        self.assertTrue(node.isValid())
        cache.setNode(node)
        element, xi = markerLocation.evaluateMeshLocation(cache, 3)
        self.assertTrue(element.isValid())

    def test_uterus_rat_1shell(self):
        """
        Test creation of rat uterus scaffold with 1 shell count.
        """
        scaffold = MeshType_3d_uterus1
        options = scaffold.getDefaultOptions("Rat 1")

        networkLayout = options.get("Network layout")
        networkLayoutSettings = networkLayout.getScaffoldSettings()
        self.assertEqual("1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-36.1,"
                         "17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-36.2,#33-34-35-36.3,"
                         "36.4-37-38-39,39-40,40-41-42",
                         networkLayoutSettings["Structure"])

        self.assertEqual(15, len(options))
        self.assertEqual(12, options.get("Number of elements around"))
        self.assertEqual(8, options.get("Number of elements around oviduct/uterine horn"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(False, options.get("Use linear through wall"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(16, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(296, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(1174, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1470, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(591, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-3.186063418775099, -2.67365658085329, -0.25], 1.0E-6)
        assertAlmostEqualList(self, maximums, [1.4, 2.67365658085329, 0.25], 1.0E-6)

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
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 14.809112350090114, delta=5.0E-2)
        self.assertAlmostEqual(volume, 1.4588529735165636, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(33, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 8,
            "body of uterus": 64,
            "left uterine horn": 96,
            "vagina": 24,
            "uterus": 264
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        # refine 2x2x2 and check result
        # first remove faces/lines and any surface annotation groups as they are re-added by defineFaceAnnotations
        removeAnnotationGroups = []
        for annotationGroup in annotationGroups:
            if annotationGroup.getDimension() in [1, 2]:
                removeAnnotationGroups.append(annotationGroup)
        for annotationGroup in removeAnnotationGroups:
            if "cervix" in annotationGroup.getName():
                continue
            annotationGroups.remove(annotationGroup)
        self.assertEqual(19, len(annotationGroups))
        # must keep all faces and lines as used for refinement

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(19, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(2368, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(8248, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(9412, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(3531, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        expectedSizes3d = {
            "fundus of uterus": 8,
            "body of uterus": 64,
            "left uterine horn": 96,
            "vagina": 24,
            "uterus": 264
        }
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(0, markerNodes.getSize())

    def test_uterus_rat_0shell(self):
        """
        Test creation of rat uterus scaffold with 0 shell count.
        """
        scaffold = MeshType_3d_uterus1
        options = scaffold.getDefaultOptions("Rat 1")
        options["Number of elements through wall"] = 0

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(13, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(0, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(260, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(534, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(273, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-3.186063418775099, -2.67365658085329, -0.25], 1.0E-6)
        assertAlmostEqualList(self, maximums, [1.4, 2.67365658085329, 0.25], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 14.808339326907051, delta=5.0E-2)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(20, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "fundus of uterus": 8,
            "body of uterus": 36,
            "left uterine horn": 96,
            "vagina": 24,
            "uterus": 236
            }

        meshes = [mesh1d, mesh2d, mesh3d]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name], size, name)

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        self.assertEqual(20, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(0, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(1040, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(2108, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1067, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        meshes = [mesh1d, mesh2d, mesh3d]
        sizeScales = [2, 4, 8]
        for name in expectedSizes3d:
            term = get_uterus_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(meshes[annotationGroup.getDimension() - 1]).getSize()
            self.assertEqual(expectedSizes3d[name] * sizeScales[annotationGroup.getDimension() - 1], size, name)

        # test finding a marker in refined scaffold
        markerGroup = refineFieldmodule.findFieldByName("marker").castGroup()
        refinedNodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNodes = markerGroup.getNodesetGroup(refinedNodes)
        self.assertEqual(0, markerNodes.getSize())


if __name__ == "__main__":
    unittest.main()
