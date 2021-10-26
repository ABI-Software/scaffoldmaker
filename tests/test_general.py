import unittest
from opencmiss.maths.vectorops import magnitude
from opencmiss.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_3d_box1 import MeshType_3d_box1
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.scaffolds import Scaffolds
from testutils import assertAlmostEqualList

from scaffoldmaker.utils.zinc_utils import identifier_ranges_from_string, identifier_ranges_to_string, \
    mesh_group_add_identifier_ranges, mesh_group_to_identifier_ranges, \
    nodeset_group_add_identifier_ranges, nodeset_group_to_identifier_ranges


class GeneralScaffoldTestCase(unittest.TestCase):

    def test_transformation(self):
        """
        Test transformation of a box scaffold with scaffold package.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_box1)

        tmpScale = scaffoldPackage.getScale()
        TOL = 1.0E-7
        FINETOL = 1.0E-12
        assertAlmostEqualList(self, tmpScale, [ 1.0, 1.0, 1.0 ], delta=FINETOL)
        newScale = [ 2.0, 1.5, 0.5 ]
        scaffoldPackage.setScale(newScale)
        tmpScale = scaffoldPackage.getScale()
        assertAlmostEqualList(self, tmpScale, newScale, delta=FINETOL)

        tmpRotation = scaffoldPackage.getRotation()
        assertAlmostEqualList(self, tmpRotation, [ 0.0, 0.0, 0.0 ], delta=FINETOL)
        newRotation = [ 30.0, -10.0, 90.0 ]
        scaffoldPackage.setRotation(newRotation)
        tmpRotation = scaffoldPackage.getRotation()
        assertAlmostEqualList(self, tmpRotation, newRotation, delta=FINETOL)

        tmpTranslation = scaffoldPackage.getTranslation()
        assertAlmostEqualList(self, tmpTranslation, [ 0.0, 0.0, 0.0 ], delta=FINETOL)
        newTranslation = [ 0.5, 1.2, -0.1 ]
        scaffoldPackage.setTranslation(newTranslation)
        tmpTranslation = scaffoldPackage.getTranslation()
        assertAlmostEqualList(self, tmpTranslation, newTranslation, delta=FINETOL)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())

        scaffoldPackage.generate(region)

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.assertEqual(8, nodes.getSize())

        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        self.assertTrue(coordinates.isValid())
        minimums, maximums = evaluateFieldNodesetRange(coordinates, nodes)
        assertAlmostEqualList(self, minimums, [  2.744244002293470e-01,  6.367511648575830e-01, -1.000000000000000e-01 ], delta=TOL)
        assertAlmostEqualList(self, maximums, [  2.455737063904887e+00,  2.184807753012208e+00,  1.724507984852172e+00 ], delta=TOL)

        node = nodes.findNodeByIdentifier(8)
        self.assertTrue(node.isValid())
        fieldcache = fieldmodule.createFieldcache()
        fieldcache.setNode(node)
        result, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(self, x , [  2.230161464134234e+00,  1.621558917869791e+00,  1.724507984852172e+00 ], delta=TOL)
        # derivative magnitudes must also equal scale
        result, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(self, d1, [  1.705737064039425e+00,  9.848077530127952e-01,  3.472963553408093e-01 ], delta=TOL)
        self.assertAlmostEqual(newScale[0], magnitude(d1), delta=TOL)
        result, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(self, d2, [ -2.255755995328457e-01, -1.302361332111701e-01,  1.477211629352659e+00 ], delta=TOL)
        self.assertAlmostEqual(newScale[1], magnitude(d2), delta=TOL)
        result, d3 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(self, d3, [  2.499999998128999e-01, -4.330127019169794e-01,  0.000000000000000e+00 ], delta=TOL)
        self.assertAlmostEqual(newScale[2], magnitude(d3), delta=TOL)

    def test_user_annotation_groups(self):
        """
        Test user annotation group on heartatria1 scaffold with scaffold package.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_heartatria1)
        # can't add user annotation groups until generate is called()
        try:
            annotationGroup = scaffoldPackage.createUserAnnotationGroup()
            self.assertTrue(False)  # should never get here as above raises expection
        except:
            pass

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        fieldmodule = region.getFieldmodule()

        scaffoldPackage.generate(region)

        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(22, len(annotationGroups))

        endocardium_of_la = scaffoldPackage.findAnnotationGroupByName('endocardium of left atrium')
        self.assertTrue(isinstance(endocardium_of_la, AnnotationGroup))
        self.assertFalse(scaffoldPackage.isUserAnnotationGroup(endocardium_of_la))
        self.assertFalse(scaffoldPackage.deleteAnnotationGroup(endocardium_of_la))  # can't delete auto annotation groups

        annotationGroup1 = scaffoldPackage.createUserAnnotationGroup()
        self.assertEqual('group1', annotationGroup1.getName())  # default name
        self.assertIsNone(annotationGroup1.getId())
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup1))
        self.assertEqual(-1, annotationGroup1.getDimension())  # -1 = empty
        group = annotationGroup1.getGroup()
        self.assertTrue(group.isValid())
        mesh2d = fieldmodule.findMeshByDimension(2)
        meshGroup = group.createFieldElementGroup(mesh2d).getMeshGroup()
        mesh_group_add_identifier_ranges(meshGroup, [[1,2],[4,4]])
        self.assertEqual(3, meshGroup.getSize())
        self.assertEqual(2, annotationGroup1.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup))
        self.assertEqual('1-2,4', identifier_ranges_string)

        annotationGroup2 = scaffoldPackage.createUserAnnotationGroup(('bob', 'BOB:1'))
        self.assertEqual('bob', annotationGroup2.getName())
        self.assertEqual('BOB:1', annotationGroup2.getId())
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup2))
        group = annotationGroup2.getGroup()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodesetGroup = group.createFieldNodeGroup(nodes).getNodesetGroup()
        nodeset_group_add_identifier_ranges(nodesetGroup, identifier_ranges_from_string('1,3-5,7'))
        self.assertEqual(5, nodesetGroup.getSize())
        self.assertEqual(0, annotationGroup2.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(nodeset_group_to_identifier_ranges(nodesetGroup))
        self.assertEqual('1,3-5,7', identifier_ranges_string)

        annotationGroup3 = scaffoldPackage.createUserAnnotationGroup()
        self.assertEqual('group2', annotationGroup3.getName())  # default name
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup3))
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(25, len(annotationGroups))

        # rename group1 to fred
        self.assertTrue(annotationGroup1.setName('fred'))
        self.assertTrue(annotationGroup1.setId('FRED:1'))
        self.assertEqual('fred', annotationGroup1.getName())
        self.assertEqual('FRED:1', annotationGroup1.getId())

        self.assertTrue(scaffoldPackage.deleteAnnotationGroup(annotationGroup3))
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(24, len(annotationGroups))

        # test serialisation
        dct = scaffoldPackage.toDict()
        self.assertEqual('3D Heart Atria 1', dct['scaffoldTypeName'])
        scaffoldType = Scaffolds().findScaffoldTypeByName(dct['scaffoldTypeName'])

        scaffoldPackage2 = ScaffoldPackage(scaffoldType, dct)
        region2 = context.createRegion()
        fieldmodule2 = region2.getFieldmodule()

        scaffoldPackage2.generate(region2)

        annotationGroups2 = scaffoldPackage2.getAnnotationGroups()
        self.assertEqual(24, len(annotationGroups2))

        annotationGroup1 = scaffoldPackage2.findAnnotationGroupByName('fred')
        self.assertEqual('fred', annotationGroup1.getName())
        self.assertEqual('FRED:1', annotationGroup1.getId())
        self.assertTrue(scaffoldPackage2.isUserAnnotationGroup(annotationGroup1))
        self.assertEqual(2, annotationGroup1.getDimension())
        mesh2d2 = fieldmodule2.findMeshByDimension(2)
        meshGroup2 = annotationGroup1.getMeshGroup(mesh2d2)
        self.assertEqual(3, meshGroup2.getSize())
        identifier_ranges_string = identifier_ranges_to_string(mesh_group_to_identifier_ranges(meshGroup2))
        self.assertEqual('1-2,4', identifier_ranges_string)

        annotationGroup2 = scaffoldPackage2.findAnnotationGroupByName('bob')
        self.assertEqual('bob', annotationGroup2.getName())
        self.assertEqual('BOB:1', annotationGroup2.getId())
        self.assertTrue(scaffoldPackage2.isUserAnnotationGroup(annotationGroup2))
        self.assertEqual(0, annotationGroup2.getDimension())
        nodes2 = fieldmodule2.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodesetGroup2 = annotationGroup2.getNodesetGroup(nodes2)
        self.assertEqual(5, nodesetGroup.getSize())
        identifier_ranges_string = identifier_ranges_to_string(nodeset_group_to_identifier_ranges(nodesetGroup2))
        self.assertEqual('1,3-5,7', identifier_ranges_string)


if __name__ == "__main__":
    unittest.main()
