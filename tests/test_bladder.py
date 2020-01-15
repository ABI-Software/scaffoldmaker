import copy
import unittest
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1
from scaffoldmaker.meshtypes.meshtype_3d_bladder1 import MeshType_3d_bladder1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from testutils import assertAlmostEqualList

class BladderScaffoldTestCase(unittest.TestCase):

    def test_bladder1(self):
        """
        Test creation of bladder scaffold.
        """
        parameterSetNames = MeshType_3d_bladder1.getParameterSetNames()
        self.assertEqual(parameterSetNames, ["Default", "Cat 1", "Rat 1"])
        ostiumDefaultScaffoldPackages = {
            'Test ostium': ScaffoldPackage(MeshType_3d_ostium1, {
                'scaffoldSettings': {
                    'Number of vessels': 1,
                    'Number of elements across common': 2,
                    'Number of elements around ostium': 8,
                    'Number of elements along': 2,
                    'Number of elements through wall': 1,  # not implemented for > 1
                    'Unit scale': 1.0,
                    'Outlet': False,
                    'Ostium diameter': 0.3,
                    'Ostium length': 0.2,
                    'Ostium wall thickness': 0.05,
                    'Ostium inter-vessel distance': 0.8,
                    'Ostium inter-vessel height': 0.0,
                    'Use linear through ostium wall': True,
                    'Vessel end length factor': 1.0,
                    'Vessel inner diameter': 0.15,
                    'Vessel wall thickness': 0.04,
                    'Vessel angle 1 degrees': 0.0,
                    'Vessel angle 1 spread degrees': 0.0,
                    'Vessel angle 2 degrees': 0.0,
                    'Use linear through vessel wall': True,
                    # 'Use cross derivatives' : False,
                    'Refine': False,
                    'Refine number of elements around': 4,
                    'Refine number of elements along': 4,
                    'Refine number of elements through wall': 1
                },
            })
        }
        ostiumOption = ostiumDefaultScaffoldPackages['Test ostium']
        options = {
            'Number of elements up neck': 8,
            'Number of elements up body': 16,
            'Number of elements around': 8,  # should be even
            'Number of elements through wall': 1,
            'Number of elements around ostium': 8,  # implemented for 8
            'Number of elements radially on annulus': 1,
            'Height': 5.0,
            'Major diameter': 6.0,
            'Minor diameter': 6.0,
            'Bladder wall thickness': 0.05,
            'Urethra diameter': 1.0,
            'Ureter': copy.deepcopy(ostiumOption),
            'Ostium position around': 0.15,
            'Ostium position up': 0.25,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements up': 4,
            'Refine number of elements through wall': 1
        }
        self.assertEqual(19, len(options))
        ostiumSettings = ostiumOption.getScaffoldSettings()
        self.assertEqual(1, ostiumSettings.get("Number of vessels"))
        self.assertEqual(8, ostiumSettings.get("Number of elements around ostium"))
        self.assertEqual(1, ostiumSettings.get("Number of elements through wall"))
        self.assertEqual(0.3, ostiumSettings.get("Ostium diameter"))
        self.assertEqual(8, options.get("Number of elements up neck"))
        self.assertEqual(16, options.get("Number of elements up body"))
        self.assertEqual(8, options.get("Number of elements around"))
        self.assertEqual(1.0, options.get("Urethra diameter"))

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        annotationGroups = MeshType_3d_bladder1.generateBaseMesh(region, options)
        self.assertEqual(2, len(annotationGroups))
        fieldmodule = region.getFieldmodule()
        self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
        if annotationGroups is not None:
            for annotationGroup in annotationGroups:
                annotationGroup.addSubelements()

        mesh3d = fieldmodule.findMeshByDimension(3)
        self.assertEqual(232, mesh3d.getSize())
        mesh2d = fieldmodule.findMeshByDimension(2)
        self.assertEqual(936, mesh2d.getSize())
        mesh1d = fieldmodule.findMeshByDimension(1)
        self.assertEqual(1183, mesh1d.getSize())
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(478, nodes.getSize())
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        self.assertEqual(0, datapoints.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [-2.996386368615517, -2.996386368615517, -6.464466094067262], 1.0E-6)
        assertAlmostEqualList(self, maximums, [2.996386368615517, 2.996386368615517, 5.0], 1.0E-6)

if __name__ == "__main__":
    unittest.main()
