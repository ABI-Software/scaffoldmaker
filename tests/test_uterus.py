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


class UterusScaffoldTestCase(unittest.TestCase):

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

        self.assertEqual(12, len(options))
        self.assertEqual(12, options.get("Number of elements around"))
        self.assertEqual(12, options.get("Number of elements around horns"))
        self.assertEqual(1, options.get("Number of elements through wall"))
        self.assertEqual(5.0, options.get("Target element density along longest segment"))
        self.assertEqual(True, options.get("Use linear through wall"))
        # test with fewer elements
        options["Number of elements around"] = 8
        options["Number of elements around horns"] = 8
        options["Target element density along longest segment"] = 3.0

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(14, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(136, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(556, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(715, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(296, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-9.361977045958657, -0.048, -8.920886838294868], 1.0E-6)
        assertAlmostEqualList(self, maximums, [9.36197704595863, 12.846826512420462, 1.09], 1.0E-6)

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
        self.assertAlmostEqual(surfaceArea, 262.06878921216844, delta=1.0E-6)
        self.assertAlmostEqual(volume, 182.50188979389358, delta=1.0E-6)

        fieldmodule.defineAllFaces()
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(region, options, annotationGroups)
        self.assertEqual(34, len(annotationGroups))

        # check some annotation groups
        expectedSizes3d = {
            "body of uterus": 72,
            "left uterine tube": 16,
            "uterine cervix": 8,
            "vagina": 24,
            'left broad ligament of uterus': 8,
            'right broad ligament of uterus': 8,
            "uterus": 136
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
            annotationGroups.remove(annotationGroup)
        self.assertEqual(14, len(annotationGroups))
        # also remove all faces and lines as not needed for refinement
        mesh2d.destroyAllElements()
        mesh1d.destroyAllElements()

        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements along'] = 2
        options['Refine number of elements around'] = 2
        options['Refine number of elements through wall'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, annotationGroups)
        scaffold.refineMesh(meshrefinement, options)
        annotationGroups = meshrefinement.getAnnotationGroups()

        refineFieldmodule.defineAllFaces()
        oldAnnotationGroups = copy.copy(annotationGroups)
        for annotationGroup in annotationGroups:
            annotationGroup.addSubelements()
        scaffold.defineFaceAnnotations(refineRegion, options, annotationGroups)
        for annotationGroup in annotationGroups:
            if annotationGroup not in oldAnnotationGroups:
                annotationGroup.addSubelements()
        self.assertEqual(34, len(annotationGroups))
#
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(1088, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(3856, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(4470, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1703, nodes.getSize())
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
        self.assertTrue(element.isValid())


if __name__ == "__main__":
    unittest.main()
