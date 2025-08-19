import unittest
from testutils import assertAlmostEqualList

from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm
from scaffoldmaker.meshtypes.meshtype_3d_wholebody1 import MeshType_3d_wholebody1
from scaffoldmaker.annotation.body_terms import get_body_term, body_terms
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class WholeBodyScaffoldTestCase(unittest.TestCase):

    def test_wholebody1(self):
        """
        Test creation of Whole-body scaffold.
        """
        scaffold = MeshType_3d_wholebody1
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human Coarse', 'Human Fine', 'Rat Coarse', 'Rat Fine'])
        options = scaffold.getDefaultOptions("Default")
        self.assertEqual(17, len(options))
        self.assertEqual(6, options['Number of elements across major'])
        self.assertEqual(6, options['Number of elements across minor'])
        self.assertEqual(1, options['Number of elements across shell'])
        self.assertEqual(1, options['Number of elements across transition'])
        self.assertEqual(5, options['Number of elements in abdomen'])
        self.assertEqual(3, options['Number of elements in thorax'])
        self.assertEqual(1, options['Number of elements in neck'])
        self.assertEqual(2, options['Number of elements in head'])
        self.assertEqual(0.2, options['Shell thickness proportion'])
        self.assertEqual(True, options['Discontinuity on the core boundary'])
        self.assertEqual(False, options['Lower half'])
        self.assertEqual(False, options['Use cross derivatives'])
        self.assertEqual(False, options['Refine'])
        self.assertEqual(1, options['Refine number of elements across major'])
        self.assertEqual(1, options['Refine number of elements along'])

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        annotationGroups = scaffold.generateMesh(region, options)[0]
        self.assertEqual(373, len(annotationGroups))

        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(220, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(724, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(803, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(662, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        # check coordinates range, sphere volume
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-0.5, -0.5, 0.0], 1.0E-6)
        assertAlmostEqualList(self, maximums, [0.5, 0.5, 3.5], 1.0E-6)

        with ChangeManager(fieldmodule):
            one = fieldmodule.createFieldConstant(1.0)
            surfaceGroup = AnnotationGroup(region, ("sphere surface", ""))
            is_exterior = fieldmodule.createFieldIsExterior()
            surfaceMeshGroup = surfaceGroup.getMeshGroup(mesh2d)
            surfaceMeshGroup.addElementsConditional(is_exterior)

            surfaceAreaField = fieldmodule.createFieldMeshIntegral(one, coordinates, surfaceMeshGroup)
            surfaceAreaField.setNumbersOfPoints(4)
            volumeField = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh3d)
            volumeField.setNumbersOfPoints(3)
        fieldcache = fieldmodule.createFieldcache()
        result, surfaceArea = surfaceAreaField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(surfaceArea, 12.55907051130264, delta=1.0E-6)
        result, volume = volumeField.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(volume, 2.7460390972760136, delta=1.0E-6)

        # check some annotationGroups:
        expectedSizes2d = {
            "skin epidermis outer surface": 88
            }
        for name in expectedSizes2d:
            term = get_body_term(name)
            group = getAnnotationGroupForTerm(annotationGroups, term)
            size = group.getMeshGroup(mesh2d).getSize()
            self.assertEqual(expectedSizes2d[name], size, name)

        # refine 8x8x8 and check result
        refineRegion = region.createRegion()
        refineFieldmodule = refineRegion.getFieldmodule()
        options['Refine number of elements across major'] = 2
        options['Refine number of elements along'] = 2
        meshrefinement = MeshRefinement(region, refineRegion, [])
        scaffold.refineMesh(meshrefinement, options)

        refineFieldmodule.defineAllFaces()
        mesh3d = refineFieldmodule.findMeshByDimension(3)
        self.assertEqual(1760, mesh3d.getSize())


if __name__ == "__main__":
    unittest.main()
