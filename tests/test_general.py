import math
import unittest

from cmlibs.maths.vectorops import dot, magnitude, normalize, sub
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange, findNodeWithName
from cmlibs.utils.zinc.group import identifier_ranges_from_string, identifier_ranges_to_string, \
    mesh_group_add_identifier_ranges, mesh_group_to_identifier_ranges, \
    nodeset_group_add_identifier_ranges, nodeset_group_to_identifier_ranges
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationMarkerNameField
from scaffoldmaker.meshtypes.meshtype_3d_box1 import MeshType_3d_box1
from scaffoldmaker.meshtypes.meshtype_3d_brainstem import MeshType_3d_brainstem1
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_stomach1 import MeshType_3d_stomach1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.scaffolds import Scaffolds
from scaffoldmaker.utils.geometry import getEllipsoidPlaneA, getEllipsoidPolarCoordinatesFromPosition, \
    getEllipsoidPolarCoordinatesTangents
from scaffoldmaker.utils.tracksurface import TrackSurface

from testutils import assertAlmostEqualList


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
        self.assertEqual(33, len(annotationGroups))

        endocardium_of_la = scaffoldPackage.findAnnotationGroupByName('left atrium endocardium')
        self.assertTrue(isinstance(endocardium_of_la, AnnotationGroup))
        self.assertFalse(scaffoldPackage.isUserAnnotationGroup(endocardium_of_la))
        self.assertFalse(scaffoldPackage.deleteAnnotationGroup(endocardium_of_la))  # can't delete auto annotation groups

        annotationGroup1 = scaffoldPackage.createUserAnnotationGroup()
        self.assertEqual('group1', annotationGroup1.getName())  # default name
        self.assertEqual('None', annotationGroup1.getId())
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup1))
        self.assertEqual(-1, annotationGroup1.getDimension())  # -1 = empty
        group = annotationGroup1.getGroup()
        self.assertTrue(group.isValid())
        mesh2d = fieldmodule.findMeshByDimension(2)
        meshGroup = group.createMeshGroup(mesh2d)
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
        nodesetGroup = group.createNodesetGroup(nodes)
        nodeset_group_add_identifier_ranges(nodesetGroup, identifier_ranges_from_string('1,3-5,7'))
        self.assertEqual(5, nodesetGroup.getSize())
        self.assertEqual(0, annotationGroup2.getDimension())
        identifier_ranges_string = identifier_ranges_to_string(nodeset_group_to_identifier_ranges(nodesetGroup))
        self.assertEqual('1,3-5,7', identifier_ranges_string)

        annotationGroup3 = scaffoldPackage.createUserAnnotationGroup()
        self.assertEqual('group2', annotationGroup3.getName())  # default name
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup3))
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(36, len(annotationGroups))

        # rename group1 to fred
        self.assertTrue(annotationGroup1.setName('fred'))
        self.assertTrue(annotationGroup1.setId('FRED:1'))
        self.assertEqual('fred', annotationGroup1.getName())
        self.assertEqual('FRED:1', annotationGroup1.getId())

        self.assertTrue(scaffoldPackage.deleteAnnotationGroup(annotationGroup3))
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(35, len(annotationGroups))

        # test serialisation
        dct = scaffoldPackage.toDict()
        self.assertEqual('3D Heart Atria 1', dct['scaffoldTypeName'])
        scaffoldType = Scaffolds().findScaffoldTypeByName(dct['scaffoldTypeName'])

        scaffoldPackage2 = ScaffoldPackage(scaffoldType, dct)
        region2 = context.createRegion()
        fieldmodule2 = region2.getFieldmodule()

        scaffoldPackage2.generate(region2)

        annotationGroups2 = scaffoldPackage2.getAnnotationGroups()
        self.assertEqual(35, len(annotationGroups2))

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

    def test_user_marker_points(self):
        """
        Test user marker point on brainstem1 scaffold which defined "brainstem coordinates".
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_brainstem1)

        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(3)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        scaffoldPackage.generate(region)
        # 1 higher than last node in scaffold. Make marker nodes from this number
        nextNodeIdentifier = scaffoldPackage.getNextNodeIdentifier()

        brainstemCoordinatesField = fieldmodule.findFieldByName("brainstem coordinates")
        self.assertTrue(brainstemCoordinatesField.isValid())

        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(18, len(annotationGroups))
        TOL = 1.0E-6  # coarse to handle find xi tolerances

        # check a built-in non-marker annotation group
        ponsGroup = scaffoldPackage.findAnnotationGroupByName('pons')
        self.assertFalse(ponsGroup.isMarker())
        self.assertRaises(AssertionError, lambda: ponsGroup.getMarkerMaterialCoordinates())
        self.assertRaises(AssertionError, lambda: ponsGroup.getMarkerLocation())
        self.assertRaises(AssertionError, lambda: ponsGroup.createMarkerNode(nextNodeIdentifier))

        # check a built-in marker annotation group
        brainstemVentralCranialPointGroup = \
            scaffoldPackage.findAnnotationGroupByName('brainstem ventral midline cranial point')
        self.assertIsNotNone(brainstemVentralCranialPointGroup)
        self.assertTrue(brainstemVentralCranialPointGroup.isMarker())
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = \
            brainstemVentralCranialPointGroup.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField)
        assertAlmostEqualList(self, [0.0, -1.0, 8.0], brainstemCoordinatesValueOut, delta=TOL)
        elementOut, xiOut = brainstemVentralCranialPointGroup.getMarkerLocation()
        self.assertEqual(235, elementOut.getIdentifier())
        assertAlmostEqualList(self, [1.0, 1.0, 0.0], xiOut, delta=TOL)
        self.assertRaises(AssertionError, lambda: brainstemVentralCranialPointGroup.createMarkerNode(nextNodeIdentifier))

        # check a non-existant annotation group
        bobGroup = scaffoldPackage.findAnnotationGroupByName("bob")
        self.assertIsNone(bobGroup)

        # now make a marker annotation named "bob" at the default location
        bobGroup = scaffoldPackage.createUserAnnotationGroup(('bob', 'BOB:1'), isMarker=True)
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(bobGroup))
        self.assertFalse(bobGroup.isMarker())
        node = bobGroup.createMarkerNode(nextNodeIdentifier)
        bobNodeIdentifier = node.getIdentifier()
        self.assertEqual(nextNodeIdentifier, bobNodeIdentifier)
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertIsNone(brainstemCoordinatesFieldOut)
        self.assertIsNone(brainstemCoordinatesValueOut)
        elementOut, xiOut = bobGroup.getMarkerLocation()
        self.assertEqual(1, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.0, 0.0, 0.0], xiOut, delta=TOL)
        # now request brainstem coordinates and let the annotation group determine its values from element:xi
        bobGroup.setMarkerMaterialCoordinates(brainstemCoordinatesField)
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField)
        # these should be precisely cos(45) but are not due to ellipse approximations
        cos45 = math.cos(math.pi / 4.0)
        assertAlmostEqualList(self, [cos45, -cos45, 0], brainstemCoordinatesValueOut, delta=TOL)
        # set element:xi location and check brainstem coordinates change
        bobGroup.setMarkerLocation(mesh.findElementByIdentifier(33), [0.0, 0.5, 0.0])
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField)
        assertAlmostEqualList(self, [cos45, -cos45, 1.5], brainstemCoordinatesValueOut, delta=TOL)
        # assign brainstem coordinates and check element:xi has moved
        bobGroup.setMarkerMaterialCoordinates(brainstemCoordinatesField, [-0.1, -0.5, 2.2])
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField)
        assertAlmostEqualList(self, [-0.1, -0.5, 2.2], brainstemCoordinatesValueOut, delta=TOL)
        elementOut, xiOut = bobGroup.getMarkerLocation()
        self.assertEqual(82, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.30536897479419056, 0.2, 0.486013009421796], xiOut, delta=TOL)

        # now make a marker annotation named "fred" with brainstem coordinates from the start
        fredGroup = scaffoldPackage.createUserAnnotationGroup(('fred', 'FRED:1'))
        # AnnotationGroup.createMarkerNode increments nextNodeIdentifier to one not used by existing node
        node = fredGroup.createMarkerNode(nextNodeIdentifier, brainstemCoordinatesField, [0.5, 0.5, 4])
        fredNodeIdentifier = node.getIdentifier()
        self.assertEqual(nextNodeIdentifier + 1, fredNodeIdentifier)
        del node
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = fredGroup.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField)
        assertAlmostEqualList(self, [0.5, 0.5, 4], brainstemCoordinatesValueOut, delta=TOL)
        elementOut, xiOut = fredGroup.getMarkerLocation()
        self.assertEqual(105, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.3452673123795837, 1.0, 0.6634646029995092], xiOut, delta=TOL)

        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(20, len(annotationGroups))

        # test renaming group assigns the marker_name field:
        fredGroup.setName("freddy")
        self.assertEqual("freddy", fredGroup.getName())
        node = fredGroup.getMarkerNode()
        markerName = getAnnotationMarkerNameField(fieldmodule)
        fieldcache = fieldmodule.createFieldcache()
        fieldcache.setNode(node)
        markerNodeName = markerName.evaluateString(fieldcache)
        self.assertEqual("freddy", markerNodeName)
        del node

        # test deleting a marker annotation group
        scaffoldPackage.deleteAnnotationGroup(fredGroup)

        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(19, len(annotationGroups))

        # check node fred has gone
        node = nodes.findNodeByIdentifier(fredNodeIdentifier)
        self.assertFalse(node.isValid())

        # re-recreate fred with just element:xi location
        fredGroup = scaffoldPackage.createUserAnnotationGroup(('fred', 'FRED:1'), isMarker=True)
        element = mesh.findElementByIdentifier(105)
        node = fredGroup.createMarkerNode(nextNodeIdentifier, element=element,
                                          xi=[0.3452673123795837, 1.0, 0.6634646029995092])
        self.assertEqual(fredNodeIdentifier, node.getIdentifier())
        del node
        elementOut, xiOut = fredGroup.getMarkerLocation()
        self.assertEqual(105, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.3452673123795837, 1.0, 0.6634646029995092], xiOut, delta=TOL)

        # check the total number of groups
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(20, len(annotationGroups))

        # test serialisation
        dct = scaffoldPackage.toDict()
        scaffoldType = Scaffolds().findScaffoldTypeByName(dct['scaffoldTypeName'])

        scaffoldPackage2 = ScaffoldPackage(scaffoldType, dct)
        region2 = context.createRegion()
        fieldmodule2 = region2.getFieldmodule()

        scaffoldPackage2.generate(region2)

        brainstemCoordinatesField2 = fieldmodule2.findFieldByName("brainstem coordinates")
        self.assertTrue(brainstemCoordinatesField2.isValid())

        annotationGroups2 = scaffoldPackage2.getAnnotationGroups()
        self.assertEqual(20, len(annotationGroups2))

        # check user markers have been defined correctly for scaffoldPackage2

        bobGroup2 = scaffoldPackage2.findAnnotationGroupByName('bob')
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = bobGroup2.getMarkerMaterialCoordinates()
        self.assertEqual(brainstemCoordinatesFieldOut, brainstemCoordinatesField2)
        assertAlmostEqualList(self, [-0.1, -0.5, 2.2], brainstemCoordinatesValueOut, delta=TOL)
        elementOut, xiOut = bobGroup2.getMarkerLocation()
        self.assertEqual(82, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.30536897479419256, 0.20000000000000018, 0.48601300942181097], xiOut, delta=TOL)

        fredGroup2 = scaffoldPackage2.findAnnotationGroupByName('fred')
        brainstemCoordinatesFieldOut, brainstemCoordinatesValueOut = fredGroup2.getMarkerMaterialCoordinates()
        self.assertIsNone(brainstemCoordinatesFieldOut)
        self.assertIsNone(brainstemCoordinatesValueOut)
        elementOut, xiOut = fredGroup2.getMarkerLocation()
        self.assertEqual(105, elementOut.getIdentifier())
        assertAlmostEqualList(self, [0.3452673123795837, 1.0, 0.6634646029995092], xiOut, delta=TOL)

    def test_deletion(self):
        """
        Test deletion of element ranges on a stomach scaffold with scaffold package.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_stomach1)
        TOL = 1.0E-7
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertTrue(region.isValid())
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        mesh3d = fieldmodule.findMeshByDimension(3)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        scaffoldPackage.generate(region)
        coordinates = fieldmodule.findFieldByName("coordinates")
        stomachCoordinatesField = fieldmodule.findFieldByName("stomach coordinates")

        # Check existing markers
        bodyAntrumGCLuminal = \
            scaffoldPackage.findAnnotationGroupByName('body-antrum junction along the greater curvature on luminal '
                                                      'surface')
        self.assertTrue(isinstance(bodyAntrumGCLuminal, AnnotationGroup))
        self.assertTrue(bodyAntrumGCLuminal.isMarker())
        bodyAntrumGCLuminalNode = bodyAntrumGCLuminal.getMarkerNode()
        bodyAntrumGCLuminalNodeIdentifier = bodyAntrumGCLuminalNode.getIdentifier()
        bodyAntrumGCLuminalMaterialCoordinates = bodyAntrumGCLuminal.getMarkerMaterialCoordinates()[1]
        element, xi = bodyAntrumGCLuminal.getMarkerLocation()
        fieldcache.setMeshLocation(element, xi)
        bodyAntrumGCLuminalCoordinates = coordinates.evaluateReal(fieldcache, coordinates.getNumberOfComponents())[1]

        # check a non-existent annotation group
        bobGroup = scaffoldPackage.findAnnotationGroupByName("bob")
        self.assertIsNone(bobGroup)

        # now make a marker annotation named "bob" at a location to be deleted
        bobGroup = scaffoldPackage.createUserAnnotationGroup(('bob', 'BOB:1'), isMarker=True)
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(bobGroup))
        self.assertFalse(bobGroup.isMarker())
        nextNodeIdentifier = scaffoldPackage.getNextNodeIdentifier()
        node = bobGroup.createMarkerNode(nextNodeIdentifier, stomachCoordinatesField, [0.900712, 0.291771, 0.391829])
        bobNodeIdentifier = node.getIdentifier()
        self.assertEqual(nextNodeIdentifier, bobNodeIdentifier)
        stomachCoordinatesFieldOut, stomachCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertEqual([0.900712, 0.291771, 0.391829], stomachCoordinatesValueOut)

        # check bob is made before deletion
        bob = scaffoldPackage.findAnnotationGroupByName('bob')
        self.assertTrue(isinstance(bob, AnnotationGroup))
        self.assertTrue(bob.isMarker())

        # make a user marker named "fred" at a location not being deleted, but without material coordinates
        # a bug was causing such markers to be deleted everywhere if any elements were deleted
        fredGroup = scaffoldPackage.createUserAnnotationGroup(('fred', 'FRED:1'))  # not setting isMarker here
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(fredGroup))
        self.assertFalse(fredGroup.isMarker())
        element300 = mesh3d.findElementByIdentifier(300)
        node = fredGroup.createMarkerNode(nextNodeIdentifier, element=element300, xi=[0.5, 0.5, 1.0])
        fredNodeIdentifier = node.getIdentifier()
        self.assertEqual(nextNodeIdentifier + 1, fredNodeIdentifier)

        # check fred is made before deletion
        fred = scaffoldPackage.findAnnotationGroupByName('fred')
        self.assertTrue(isinstance(fred, AnnotationGroup))
        self.assertTrue(fred.isMarker())

        # delete element ranges for body
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(73, len(annotationGroups))
        scaffoldPackage.deleteElementsInRanges(region, [[313, 496]])
        self.assertEqual(600, mesh3d.getSize())
        element = mesh3d.findElementByIdentifier(400)
        self.assertFalse(element.isValid())

        # check that body-antrum junction marker is still present after deletion
        self.assertIn(bodyAntrumGCLuminal, annotationGroups)
        self.assertTrue(bodyAntrumGCLuminal.isMarker())
        node = nodes.findNodeByIdentifier(bodyAntrumGCLuminalNodeIdentifier)
        self.assertTrue(node.isValid())
        fieldcache.setNode(bodyAntrumGCLuminal.getMarkerNode())
        markerMaterialCoordinates = bodyAntrumGCLuminal.getMarkerMaterialCoordinates()[1]
        element, xi = bodyAntrumGCLuminal.getMarkerLocation()
        fieldcache.setMeshLocation(element, xi)
        result, markerCoordinates = coordinates.evaluateReal(fieldcache, coordinates.getNumberOfComponents())
        assertAlmostEqualList(self, markerCoordinates, bodyAntrumGCLuminalCoordinates, delta=TOL)
        assertAlmostEqualList(self, markerMaterialCoordinates, bodyAntrumGCLuminalMaterialCoordinates, delta=TOL)

        # check that bob is deleted
        annotationGroups = scaffoldPackage.getAnnotationGroups()
        self.assertEqual(72, len(annotationGroups))
        bob = scaffoldPackage.findAnnotationGroupByName('bob')
        self.assertNotIn(bob, annotationGroups)
        node = nodes.findNodeByIdentifier(bobNodeIdentifier)
        self.assertFalse(node.isValid())

        # check that fred is still around
        self.assertIn(fred, annotationGroups)
        node = nodes.findNodeByIdentifier(fredNodeIdentifier)
        self.assertTrue(node.isValid())

    def test_utils_ellipsoid(self):
        """
        Test ellipsoid functions converting between coordinates.
        """
        a = 1.0
        b = 0.75
        c = 2.0
        aa = a * a
        bb = b * b
        cc = c * c
        TOL = 1.0E-6
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, [1.0, 0.0, 0.0])
        self.assertAlmostEqual(u, 0.0, delta=TOL)
        self.assertAlmostEqual(v, 0.5 * math.pi, delta=TOL)
        x, dx_du, dx_dv = getEllipsoidPolarCoordinatesTangents(a, b, c, u, v)
        assertAlmostEqualList(self, x, [1.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_du, [0.0, 0.75, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_dv, [0.0, 0.0, 2.0], delta=TOL)

        # test a point not on the surface
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, [1.0, 0.75, 1.0])
        self.assertAlmostEqual(u, 0.6795893730435009, delta=TOL)
        self.assertAlmostEqual(v, 2.038067043010934, delta=TOL)
        x = getEllipsoidPolarCoordinatesTangents(a, b, c, u, v)[0]
        assertAlmostEqualList(self, x, [0.6944481794937621, 0.42082645661076656, 0.9009025175142785], delta=TOL)
        mag = (x[0] * x[0]) / aa + (x[1] * x[1]) / bb + (x[2] * x[2]) / cc
        self.assertAlmostEqual(mag, 1.0, delta=TOL)
        # test the nearest point found on the surface has same polar coordinates
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, x)
        self.assertAlmostEqual(u, 0.6795893730435009, delta=TOL)
        self.assertAlmostEqual(v, 2.038067043010935, delta=TOL)

        # test a point not on the surface
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, [-0.9, -0.25, -1.2])
        self.assertAlmostEqual(u, -2.8190484300147065, delta=TOL)
        self.assertAlmostEqual(v, 0.9561038921906282, delta=TOL)
        x = getEllipsoidPolarCoordinatesTangents(a, b, c, u, v)[0]
        assertAlmostEqualList(self, x, [-0.7748223543206793, -0.19421818481623945, -1.153414567503636], delta=TOL)
        mag = (x[0] * x[0]) / aa + (x[1] * x[1]) / bb + (x[2] * x[2]) / cc
        self.assertAlmostEqual(mag, 1.0, delta=TOL)
        # test the nearest point found on the surface has same polar coordinates
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, x)
        self.assertAlmostEqual(u, -2.8190484300147065, delta=TOL)
        self.assertAlmostEqual(v, 0.9561038921906281, delta=TOL)

        u_in = math.pi / 3.0
        v_in = math.pi / 2.0
        x, dx_du, dx_dv = getEllipsoidPolarCoordinatesTangents(a, b, c, u_in, v_in)
        assertAlmostEqualList(self, x, [0.5, 0.649519052838329, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_du, [-0.8660254037844386, 0.375, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_dv, [0.0, 0.0, 2.0], delta=TOL)
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, x)
        self.assertAlmostEqual(u, u_in, delta=TOL)
        self.assertAlmostEqual(v, v_in, delta=TOL)

        u_in = -0.7 * math.pi
        v_in = 0.3 * math.pi
        x, dx_du, dx_dv = getEllipsoidPolarCoordinatesTangents(a, b, c, u_in, v_in)
        assertAlmostEqualList(self, x, [-0.4755282581475767, -0.4908813728906053, -1.1755705045849463], delta=TOL)
        assertAlmostEqualList(self, dx_du, [0.6545084971874737, -0.35664619361068256, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_dv, [-0.3454915028125262, -0.3566461936106826, 1.618033988749895], delta=TOL)
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, x)
        self.assertAlmostEqual(u, u_in, delta=TOL)
        self.assertAlmostEqual(v, v_in, delta=TOL)

        u_in = 0.35 * math.pi
        v_in = 0.65 * math.pi
        x, dx_du, dx_dv = getEllipsoidPolarCoordinatesTangents(a, b, c, u_in, v_in)
        assertAlmostEqualList(self, x, [0.4045084971874737, 0.5954194696096774, 0.9079809994790935], delta=TOL)
        assertAlmostEqualList(self, dx_du, [-0.7938926261462366, 0.3033813728906053, 0.0], delta=TOL)
        assertAlmostEqualList(self, dx_dv, [-0.20610737385376343, -0.30338137289060524, 1.7820130483767358], delta=TOL)
        u, v = getEllipsoidPolarCoordinatesFromPosition(a, b, c, x)
        self.assertAlmostEqual(u, u_in, delta=TOL)
        self.assertAlmostEqual(v, v_in, delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, 0.0], [0.0, -b, 0.0])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, -b, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, b / 3.0, 0.0], [0.0, -b, 0.0])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, -b, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, 0.0], [0.0, 0.0, c])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, 0.0, c], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, c / 4.0], [0.0, 0.0, c])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, 0.0, c], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        z = 0.25 * c
        y = -math.sqrt((b * b) * (1.0 - (z * z) / (c * c)))
        x = math.sqrt((a * a) * (1.0 - (z * z) / (c * c)))
        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, z], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, 0.0, z], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, y, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis2, [x, 0.0, 0.0], delta=TOL)
        mag = (axis2[0] * axis2[0]) / aa + (z * z) / cc
        self.assertAlmostEqual(mag, 1.0, delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, 0.0], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, y, z], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.1, 0.1], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, 0.009782267266570388, 0.1436792247347149], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, -0.735966644680461, 0.3563207752652851], delta=TOL)
        assertAlmostEqualList(self, axis2, [0.9973309128094611, 0.0, 0.0], delta=TOL)
        mag = (axis2[0] * axis2[0]) / aa + (centre[1] * centre[1]) / bb + (centre[2] * centre[2]) / cc
        self.assertAlmostEqual(mag, 1.0, delta=TOL)
        # check original midx is on axis1
        dir1 = normalize(sub(centre, [0.0, 0.1, 0.1]))
        dir2 = normalize(axis1)
        dir3 = normalize(sub([0.0, y, z], centre))
        assertAlmostEqualList(self, dir1, dir2, delta=TOL)
        assertAlmostEqualList(self, dir1, dir3, delta=TOL)

        z = 0.8 * c
        y = math.sqrt((b * b) * (1.0 - (z * z) / (c * c)))
        x = math.sqrt((a * a) * (1.0 - (y * y) / (b * b)))
        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, y, 0.0], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, y, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, 0.0, z], delta=TOL)
        assertAlmostEqualList(self, axis2, [x, 0.0, 0.0], delta=TOL)
        mag = (axis2[0] * axis2[0]) / aa + (y * y) / bb
        self.assertAlmostEqual(mag, 1.0, delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, 0.0, 0.0], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, 0.0, 0.0], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, y, z], delta=TOL)
        assertAlmostEqualList(self, axis2, [a, 0.0, 0.0], delta=TOL)

        centre, axis1, axis2 = getEllipsoidPlaneA(a, b, c, [0.0, -0.3, -0.1], [0.0, y, z])
        assertAlmostEqualList(self, centre, [0.0, -0.10732946298984034, 0.3367198838896952], delta=TOL)
        assertAlmostEqualList(self, axis1, [0.0, 0.5573294629898402, 1.2632801161103049], delta=TOL)
        assertAlmostEqualList(self, axis2, [0.9752823267321079, 0.0, 0.0], delta=TOL)
        mag = (axis2[0] * axis2[0]) / aa + (centre[1] * centre[1]) / bb + (centre[2] * centre[2]) / cc
        self.assertAlmostEqual(mag, 1.0, delta=TOL)
        # check original midx is on axis1
        dir1 = normalize(sub(centre, [0.0, -0.3, -0.1]))
        dir2 = normalize(axis1)
        dir3 = normalize(sub([0.0, y, z], centre))
        assertAlmostEqualList(self, dir1, dir2, delta=TOL)
        assertAlmostEqualList(self, dir1, dir3, delta=TOL)

    def test_track_surface_intersection(self):
        """
        Test finding points on intersection between 2 track surfaces.
        """
        surf1_x = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        surf1_d1 = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        surf1_d2 = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        surf1 = TrackSurface(1, 1, surf1_x, surf1_d1, surf1_d2)

        surf2_x = [[0.1, 0.0, -0.2], [0.9, 0.0, 0.6], [0.1, 1.0, 0.6], [0.9, 1.0, -0.2]]
        surf2_d1 = [[0.9, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        surf2_d2 = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        surf2 = TrackSurface(1, 1, surf2_x, surf2_d1, surf2_d2)

        xi_tol = 1.0E-5
        p1, p1x = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.25, 0.1), surf2.createPositionProportion(0.5, 0.5))
        self.assertEqual(p1.e1, 0)
        self.assertEqual(p1.e2, 0)
        self.assertEqual(p1.xi1, 0.34801277492351584)
        self.assertEqual(p1.xi2, 0.12877821314270047)

        p2, p2x = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.5, 0.5), surf2.createPositionProportion(0.5, 0.5))
        self.assertEqual(p2.e1, 0)
        self.assertEqual(p2.e2, 0)
        self.assertEqual(p2.xi1, 0.71452743629564)
        self.assertEqual(p2.xi2, 0.7304540495036337)

        p3, p3x = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.85, 0.75), surf2.createPositionProportion(0.5, 0.5))
        self.assertEqual(p3.e1, 0)
        self.assertEqual(p3.e2, 0)
        self.assertEqual(p3.xi1, 0.837522778571395)
        self.assertEqual(p3.xi2, 0.6783464093929652)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # surf1.generateMesh(region)
        # surf2.generateMesh(region)
        # fieldmodule = region.getFieldmodule()
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        # nodetemplate.defineField(coordinates)
        # fieldcache = fieldmodule.createFieldcache()
        # pp = [p1, p2, p3]
        # px = [p1x, p2x, p3x]
        # for i in range(len(pp)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     coordinates.assignReal(fieldcache, px[i])
        # fieldmodule.defineAllFaces()
        # region.writeFile("C:\\Users\\gchr006\\tmp\\tracksurface_intersection.exf")


if __name__ == "__main__":
    unittest.main()
