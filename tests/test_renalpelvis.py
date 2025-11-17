import copy
import math
import unittest

from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager

from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.annotation.ureter_terms import get_ureter_term
from scaffoldmaker.meshtypes.meshtype_3d_renal_pelvis1 import MeshType_3d_renal_pelvis1
from scaffoldmaker.utils.meshrefinement import MeshRefinement

from testutils import assertAlmostEqualList


class RenalPelviScaffoldTestCase(unittest.TestCase):

    def test_renal_pelvis_human(self):
        """
        Test creation of human renal pelvis scaffold.
        """
        scaffold = MeshType_3d_renal_pelvis1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Rat 1"])
        options = scaffold.getDefaultOptions("Human 1")

        self.assertEqual(12, len(options))
        self.assertEqual(8, options["Elements count around"])
        self.assertEqual(1, options["Elements count through shell"])
        self.assertEqual([0], options["Annotation elements counts around"])
        self.assertEqual(4.0, options["Target element density along longest segment"])
        self.assertEqual(False, options["Use linear through shell"])
        self.assertEqual(True, options["Use outer trim surfaces"])
        self.assertEqual(False, options["Show trim surfaces"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

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

        self.assertEqual(8, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(936, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(3428, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(4073, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1582, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-4
        assertAlmostEqualList(self, minimums, [-1.8502096232770495, -2.9281582229743637, -0.19871303721933536], tol)
        assertAlmostEqualList(self, maximums, [2.8357557021117126, 1.1993129591292884, 0.19871303721933536], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 0.7814475699705258, delta=tol)
            self.assertAlmostEqual(surfaceArea, 20.83194524004843, delta=tol)

        # check some annotation groups:

        expectedSizes3d = {
            "core": (360, 0.3846094864179342),
            "major calyx": (48, 0.029128480517289796),
            "minor calyx": (160, 0.036662643737561694),
            "renal pelvis": (256, 0.16709749228397455),
            "renal pyramid": (680, 0.6143476950036447),
            "shell": (576, 0.39683570086969006),
            "ureter": (64, 0.10778687283916291)
            }
        for name in expectedSizes3d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0], size, name)
            volumeMeshGroup = annotationGroup.getMeshGroup(mesh3d)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, volumeMeshGroup)
            volumeField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(volume, expectedSizes3d[name][1], delta=tol)

        expectedSizes2d = {
            "major calyx": (208, 2.978894400047078),
            "minor calyx": (692, 4.4296179632944614),
            "renal pelvis": (1068, 13.807110887897247),
            "renal pyramid": (2440, 25.945226755465498),
            "ureter": (264, 7.203840960105056)
            }
        for name in expectedSizes2d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name][0], size, name)
            surfaceMeshGroup = annotationGroup.getMeshGroup(mesh2d)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, expectedSizes2d[name][1], delta=tol)

        # refine 2x2x2 and check result
        annotationGroups = originalAnnotationGroups

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

        self.assertEqual(8, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(7488, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(24944, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(27474, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(10019, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0] * (refineNumberOfElements ** 3), size, name)


    def test_renal_pelvis_rat(self):
        """
        Test creation of rat renal pelvis scaffold.
        """
        scaffold = MeshType_3d_renal_pelvis1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Human 1", "Rat 1"])
        options = scaffold.getDefaultOptions("Rat 1")

        self.assertEqual(12, len(options))
        self.assertEqual(8, options["Elements count around"])
        self.assertEqual(1, options["Elements count through shell"])
        self.assertEqual([0], options["Annotation elements counts around"])
        self.assertEqual(4.0, options["Target element density along longest segment"])
        self.assertEqual(False, options["Use linear through shell"])
        self.assertEqual(True, options["Use outer trim surfaces"])
        self.assertEqual(False, options["Show trim surfaces"])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

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
        self.assertEqual(8, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(148, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(564, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(691, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(276, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # Check coordinates range, volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        tol = 1.0E-4
        assertAlmostEqualList(self, minimums, [-0.29836657719423826, -2.221051441787816, -0.19868582233021914], tol)
        assertAlmostEqualList(self, maximums, [2.128648920925165, 1.2000000000000002, 0.19868582233021914], tol)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            isExterior = fieldmodule.createFieldIsExterior()
            mesh2d = fieldmodule.findMeshByDimension(2)
            fieldcache = fieldmodule.createFieldcache()

            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(isExterior, coordinates, mesh2d)
            surfaceAreaField.setNumbersOfPoints(4)
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)

            self.assertAlmostEqual(volume, 0.12005376065418445, delta=tol)
            self.assertAlmostEqual(surfaceArea, 5.096473140665664, delta=tol)

        # check some annotation groups:

        expectedSizes3d = {
            "core": (36, 0.03842050612475523),
            "major calyx": (8, 0.005543006096630809),
            "minor calyx": (16, 0.007802597971432378),
            "renal pelvis": (80, 0.05861664077426631),
            "renal pyramid": (68, 0.06143679397579967),
            "shell": (112, 0.08163292862531038),
            "ureter": (64, 0.050814042802833914)
            }
        for name in expectedSizes3d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0], size, name)
            volumeMeshGroup = annotationGroup.getMeshGroup(mesh3d)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, volumeMeshGroup)
            volumeField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, volume = volumeField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(volume, expectedSizes3d[name][1], delta=tol)

        expectedSizes2d = {
            "major calyx": (40, 0.5472137937926038),
            "minor calyx": (72, 0.816589508937787),
            "renal pelvis": (328, 5.171627911264767),
            "renal pyramid": (244, 2.595915075303359),
            "ureter": (264, 4.371043844773501)
            }
        for name in expectedSizes2d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            annotationGroup = getAnnotationGroupForTerm(annotationGroups, term)
            size = annotationGroup.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name][0], size, name)
            surfaceMeshGroup = annotationGroup.getMeshGroup(mesh2d)
            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            fieldcache = fieldmodule.createFieldcache()
            result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
            self.assertEqual(result, RESULT_OK)
            self.assertAlmostEqual(surfaceArea, expectedSizes2d[name][1], delta=tol)

        # refine 2x2x2 and check result
        annotationGroups = originalAnnotationGroups

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

        self.assertEqual(8, len(annotationGroups))

        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(1184, mesh3d.getSize())
        mesh2d = refineFieldmodule.findMeshByDimension(2)
        self.assertEqual(4032, mesh2d.getSize())
        mesh1d = refineFieldmodule.findMeshByDimension(1)
        self.assertEqual(4526, mesh1d.getSize())
        nodes = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(1679, nodes.getSize())
        datapoints = refineFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check some refined annotationGroups:
        for name in expectedSizes3d:
            term = get_ureter_term(name) if name == "ureter" else get_kidney_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh3d).getSize()
            self.assertEqual(expectedSizes3d[name][0] * (refineNumberOfElements ** 3), size, name)


if __name__ == "__main__":
    unittest.main()
