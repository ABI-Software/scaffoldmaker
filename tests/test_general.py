import math
import unittest

from cmlibs.maths.vectorops import dot, magnitude, normalize, sub
# from cmlibs.utils.zinc.field import find_or_create_field_coordinates, find_or_create_field_group
from cmlibs.utils.zinc.finiteelement import evaluateFieldNodesetRange
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
from scaffoldmaker.utils.interpolation import evaluateCoordinatesOnCurve, getCubicHermiteCurvesLength, \
    getNearestLocationBetweenCurves, getNearestLocationOnCurve
from scaffoldmaker.utils.networkmesh import getPathRawTubeCoordinates, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition
# from scaffoldmaker.utils.zinc_utils import generateCurveMesh

from testutils import assertAlmostEqualList


class GeneralScaffoldTestCase(unittest.TestCase):

    def test_transformation(self):
        """
        Test transformation of a box scaffold with scaffold package.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_box1)

        tmpScale = scaffoldPackage.getScale()
        TOL = 1.0E-7
        FINE_TOL = 1.0E-12
        assertAlmostEqualList(self, tmpScale, [1.0, 1.0, 1.0], delta=FINE_TOL)
        newScale = [2.0, 1.5, 0.5]
        scaffoldPackage.setScale(newScale)
        tmpScale = scaffoldPackage.getScale()
        assertAlmostEqualList(self, tmpScale, newScale, delta=FINE_TOL)

        tmpRotation = scaffoldPackage.getRotation()
        assertAlmostEqualList(self, tmpRotation, [0.0, 0.0, 0.0], delta=FINE_TOL)
        newRotation = [30.0, -10.0, 90.0]
        scaffoldPackage.setRotation(newRotation)
        tmpRotation = scaffoldPackage.getRotation()
        assertAlmostEqualList(self, tmpRotation, newRotation, delta=FINE_TOL)

        tmpTranslation = scaffoldPackage.getTranslation()
        assertAlmostEqualList(self, tmpTranslation, [0.0, 0.0, 0.0], delta=FINE_TOL)
        newTranslation = [0.5, 1.2, -0.1]
        scaffoldPackage.setTranslation(newTranslation)
        tmpTranslation = scaffoldPackage.getTranslation()
        assertAlmostEqualList(self, tmpTranslation, newTranslation, delta=FINE_TOL)

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
        assertAlmostEqualList(
            self, minimums, [2.744244002293470e-01, 6.367511648575830e-01, -1.000000000000000e-01], delta=TOL)
        assertAlmostEqualList(
            self, maximums, [2.455737063904887e+00, 2.184807753012208e+00,  1.724507984852172e+00], delta=TOL)

        node = nodes.findNodeByIdentifier(8)
        self.assertTrue(node.isValid())
        fieldcache = fieldmodule.createFieldcache()
        fieldcache.setNode(node)
        result, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(
            self, x, [2.230161464134234e+00,  1.621558917869791e+00,  1.724507984852172e+00], delta=TOL)
        # derivative magnitudes must also equal scale
        result, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(
            self, d1, [1.705737064039425e+00,  9.848077530127952e-01,  3.472963553408093e-01], delta=TOL)
        self.assertAlmostEqual(newScale[0], magnitude(d1), delta=TOL)
        result, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(
            self, d2, [-2.255755995328457e-01, -1.302361332111701e-01,  1.477211629352659e+00], delta=TOL)
        self.assertAlmostEqual(newScale[1], magnitude(d2), delta=TOL)
        result, d3 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
        self.assertEqual(RESULT_OK, result)
        assertAlmostEqualList(
            self, d3, [2.499999998128999e-01, -4.330127019169794e-01,  0.000000000000000e+00], delta=TOL)
        self.assertAlmostEqual(newScale[2], magnitude(d3), delta=TOL)

    def test_user_annotation_groups(self):
        """
        Test user annotation group on heartatria1 scaffold with scaffold package.
        """
        scaffoldPackage = ScaffoldPackage(MeshType_3d_heartatria1)
        # can't add user annotation groups until generate is called()
        self.assertRaises(AssertionError, lambda: scaffoldPackage.createUserAnnotationGroup())

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
        self.assertFalse(scaffoldPackage.deleteAnnotationGroup(endocardium_of_la))  # can't delete auto annotation group

        annotationGroup1 = scaffoldPackage.createUserAnnotationGroup()
        self.assertEqual('group1', annotationGroup1.getName())  # default name
        self.assertEqual('None', annotationGroup1.getId())
        self.assertTrue(scaffoldPackage.isUserAnnotationGroup(annotationGroup1))
        self.assertEqual(-1, annotationGroup1.getDimension())  # -1 = empty
        group = annotationGroup1.getGroup()
        self.assertTrue(group.isValid())
        mesh2d = fieldmodule.findMeshByDimension(2)
        meshGroup = group.createMeshGroup(mesh2d)
        mesh_group_add_identifier_ranges(meshGroup, [[1, 2], [4, 4]])
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
        self.assertRaises(
            AssertionError, lambda: brainstemVentralCranialPointGroup.createMarkerNode(nextNodeIdentifier))

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
        node = bobGroup.createMarkerNode(nextNodeIdentifier, stomachCoordinatesField, [-0.165034, 0.12313, 0.484541])
        bobNodeIdentifier = node.getIdentifier()
        self.assertEqual(nextNodeIdentifier, bobNodeIdentifier)
        stomachCoordinatesFieldOut, stomachCoordinatesValueOut = bobGroup.getMarkerMaterialCoordinates()
        self.assertEqual([-0.165034, 0.12313, 0.484541], stomachCoordinatesValueOut)

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
        self.assertEqual(824, mesh3d.getSize())
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
        self.assertEqual(70, len(annotationGroups))
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

    def test_curve_nearest_intersections(self):
        """
        Test finding nearest points and intersections for curves.
        """
        curve1_x = [[0.0, 0.0, 0.0], [1.0, 0.0, 1.0]]
        curve1_d1 = [[1.0, 0.0, 1.0], [1.0, 0.0, 1.0]]
        curve2_x = [[0.0, -0.2, 0.7], [1.0, -0.3, 0.3]]
        curve2_d1 = [[1.0, -0.1, -0.6], [1.0, -0.1, 0.0]]
        curve3_x = [[0.0, 0.3, 0.5], [1.0, 0.7, 0.5]]
        curve3_d1 = [[1.0, -0.2, 0.0], [1.0, -0.2, 0.0]]
        curve4_x = [[0.4, 0.4, 0.0], [0.3, 0.3, 0.4]]
        curve4_d1 = [[0.0, 0.0, 0.5], [0.0, 0.0, 0.5]]
        loop1_x = [[0.1, 0.5, 0.5], [0.9, 0.5, 0.5]]
        loop1_d1 = [[0.0, -1.5, 0.0], [0.0, 1.5, 0.0]]

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        t1x = [0.1, 0.2, 0.3]
        p1, p1x = getNearestLocationOnCurve(curve1_x, curve1_d1, t1x)
        self.assertEqual(p1[0], 0)
        self.assertAlmostEqual(p1[1], 0.2, delta=XI_TOL)
        assertAlmostEqualList(self, [0.2, 0.0, 0.2], p1x, delta=X_TOL)
        p1t = evaluateCoordinatesOnCurve(curve1_x, curve1_d1, p1, derivative=True)[1]
        # nearest projection is normal to curve
        self.assertAlmostEqual(dot(sub(t1x, p1x), p1t), 0.0, delta=X_TOL)

        t2x = [0.1, 0.2, 0.7]
        p2, p2x = getNearestLocationOnCurve(loop1_x, loop1_d1, t2x, loop=True)
        self.assertEqual(p2[0], 0)
        self.assertAlmostEqual(p2[1], 0.19457868521558017, delta=XI_TOL)
        assertAlmostEqualList(self, [0.1790790431027201, 0.26492322619542114, 0.5], p2x, delta=X_TOL)
        p2t = evaluateCoordinatesOnCurve(loop1_x, loop1_d1, p2, loop=True, derivative=True)[1]
        # nearest projection is normal to curve
        self.assertAlmostEqual(dot(sub(t2x, p2x), p2t), 0.0, delta=X_TOL)

        t3x = [-0.1, 0.2, 0.4]
        p3, p3x = getNearestLocationOnCurve(curve3_x, curve3_d1, t3x, startLocation=(0, 0.7))
        self.assertEqual(p3[0], 0)
        self.assertAlmostEqual(p3[1], 0.0, delta=XI_TOL)
        assertAlmostEqualList(self, curve3_x[0], p3x, delta=X_TOL)

        p4, op4, p4intersects = getNearestLocationBetweenCurves(curve1_x, curve1_d1, curve2_x, curve2_d1)
        p4x, p4t = evaluateCoordinatesOnCurve(curve1_x, curve1_d1, p4, derivative=True)
        op4x = evaluateCoordinatesOnCurve(curve2_x, curve2_d1, op4)
        self.assertFalse(p4intersects)
        self.assertEqual(p4[0], 0)
        self.assertAlmostEqual(p4[1], 0.44315296477065264, delta=XI_TOL)
        # nearest projection is normal to curve
        self.assertAlmostEqual(dot(sub(op4x, p4x), p4t), 0.0, delta=X_TOL)

        p5, op5, p5intersects = getNearestLocationBetweenCurves(
            loop1_x, loop1_d1, curve1_x, curve1_d1, startLocation=(1, 0.9), nLoop=True)
        p5x, p5t = evaluateCoordinatesOnCurve(loop1_x, loop1_d1, p5, loop=True, derivative=True)
        op5x = evaluateCoordinatesOnCurve(curve1_x, curve1_d1, op5)
        self.assertFalse(p5intersects)
        self.assertEqual(p5[0], 0)
        self.assertAlmostEqual(p5[1], 0.5, delta=XI_TOL)
        # nearest projection is normal to curve
        self.assertAlmostEqual(dot(sub(op5x, p5x), p5t), 0.0, delta=X_TOL)

        p6, op6, p6intersects = getNearestLocationBetweenCurves(
            curve3_x, curve3_d1, loop1_x, loop1_d1, startLocation=(0, 0.0), oLoop=True)
        p6x = evaluateCoordinatesOnCurve(curve3_x, curve3_d1, p6)
        op6x = evaluateCoordinatesOnCurve(loop1_x, loop1_d1, op6, loop=True)
        self.assertTrue(p6intersects)
        self.assertEqual(p6[0], 0)
        self.assertAlmostEqual(p6[1], 0.14995988970099794, delta=XI_TOL)
        assertAlmostEqualList(self, p6x, op6x, delta=X_TOL)

        p7, op7, p7intersects = getNearestLocationBetweenCurves(
            loop1_x, loop1_d1, curve3_x, curve3_d1, startLocation=(0, 0.9), nLoop=True)
        p7x = evaluateCoordinatesOnCurve(loop1_x, loop1_d1, p7, loop=True)
        op7x = evaluateCoordinatesOnCurve(curve3_x, curve3_d1, op7)
        self.assertTrue(p7intersects)
        self.assertEqual(p7[0], 1)
        self.assertAlmostEqual(p7[1], 0.152207831857098, delta=XI_TOL)
        assertAlmostEqualList(self, p7x, op7x, delta=X_TOL)

        p8, op8, p8intersects = getNearestLocationBetweenCurves(
            curve4_x, curve4_d1, loop1_x, loop1_d1, startLocation=(0, 0.5), oLoop=True)
        p8x = evaluateCoordinatesOnCurve(curve4_x, curve4_d1, p8, loop=True)
        op8x, op8t = evaluateCoordinatesOnCurve(loop1_x, loop1_d1, op8, derivative=True)
        self.assertFalse(p8intersects)
        self.assertEqual(p8[0], 0)
        self.assertAlmostEqual(p8[1], 1.0, delta=XI_TOL)
        # nearest projection is normal to loop1
        self.assertAlmostEqual(dot(sub(p8x, op8x), op8t), 0.0, delta=X_TOL)

        # context = Context("Curve nearest and intersection")
        # region = context.getDefaultRegion()
        # generateCurveMesh(region, curve1_x, curve1_d1)
        # generateCurveMesh(region, curve2_x, curve2_d1)
        # generateCurveMesh(region, curve3_x, curve3_d1)
        # generateCurveMesh(region, curve4_x, curve4_d1)
        # generateCurveMesh(region, loop1_x, loop1_d1, loop=True)
        # fieldmodule = region.getFieldmodule()
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # pointCoordinates = create_field_coordinates(fieldmodule, "point_coordinates", managed=True)
        # nodetemplate.defineField(pointCoordinates)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [t1x, p1x, t2x, p2x, t3x, p3x, p4x, op4x, p5x, op5x, p6x, p7x, p8x, op8x]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     pointCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])

    def test_curve_track_surface_nearest_intersection(self):
        """
        Test finding nearest/intersection points on a curve and a track surface.
        """
        curve3_x = [[0.0, 0.3, 0.5], [1.0, 0.7, 0.5]]
        curve3_d1 = [[1.0, -0.2, 0.0], [1.0, -0.2, 0.0]]
        curve4_x = [[0.4, 0.4, 0.0], [0.3, 0.3, 0.4]]
        curve4_d1 = [[0.0, 0.0, 0.5], [0.0, 0.0, 0.5]]
        surf2_x = [[0.1, 0.0, -0.2], [0.9, 0.0, 0.5], [0.1, 1.0, 0.6], [0.9, 1.0, -0.3]]
        surf2_d1 = [[0.9, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        surf2_d2 = [[0.0, 1.0, 0.0], [0.0, 1.0, -0.8], [0.0, 1.0, 0.0], [0.0, 1.0, -0.8]]
        surf2 = TrackSurface(1, 1, surf2_x, surf2_d1, surf2_d2)

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        p1, cp1, p1intersects = surf2.findNearestPositionOnCurve(curve3_x, curve3_d1)
        self.assertFalse(p1intersects)
        self.assertEqual(cp1[0], 0)
        self.assertAlmostEqual(cp1[1], 0.4432272731754388, delta=XI_TOL)
        p1x, p1d1, p1d2 = surf2.evaluateCoordinates(p1, derivatives=True)
        cp1x = evaluateCoordinatesOnCurve(curve3_x, curve3_d1, cp1)
        # nearest projection is normal to surf2
        delta = sub(cp1x, p1x)
        self.assertAlmostEqual(dot(delta, p1d1), 0.0, delta=X_TOL)
        self.assertAlmostEqual(dot(delta, p1d2), 0.0, delta=X_TOL)

        p2, cp2, p2intersects = surf2.findNearestPositionOnCurve(curve4_x, curve4_d1)
        self.assertTrue(p2intersects)
        self.assertEqual(cp2[0], 0)
        self.assertAlmostEqual(cp2[1], 0.23080643876585413, delta=XI_TOL)
        p2x = surf2.evaluateCoordinates(p2)
        cp2x = evaluateCoordinatesOnCurve(curve4_x, curve4_d1, cp2)
        assertAlmostEqualList(self, p2x, cp2x, delta=X_TOL)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # generateCurveMesh(region, curve3_x, curve3_d1)
        # generateCurveMesh(region, curve4_x, curve4_d1)
        # surf2.generateMesh(region)
        # fieldmodule = region.getFieldmodule()
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # pointCoordinates = create_field_coordinates(fieldmodule, "point_coordinates", managed=True)
        # nodetemplate.defineField(pointCoordinates)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [p1x, cp1x, cp2x]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     pointCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])

    def test_track_surface_intersection(self):
        """
        Test finding points on intersection between 2 track surfaces.
        """
        surf1_x = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        surf1_d1 = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        surf1_d2 = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        surf1 = TrackSurface(1, 1, surf1_x, surf1_d1, surf1_d2)

        surf2_x = [[0.1, 0.0, -0.2], [0.9, 0.0, 0.5], [0.1, 1.0, 0.6], [0.9, 1.0, -0.3]]
        surf2_d1 = [[0.9, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        surf2_d2 = [[0.0, 1.0, 0.0], [0.0, 1.0, -0.8], [0.0, 1.0, 0.0], [0.0, 1.0, -0.8]]
        surf2 = TrackSurface(1, 1, surf2_x, surf2_d1, surf2_d2)

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        p1, op1, p1x, p1t, p1bdy = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.25, 0.1), surf2.createPositionProportion(0.5, 0.5))
        self.assertEqual(p1.e1, 0)
        self.assertEqual(p1.e2, 0)
        self.assertAlmostEqual(p1.xi1, 0.38287186861026046, delta=XI_TOL)
        self.assertAlmostEqual(p1.xi2, 0.14232301624230445, delta=XI_TOL)
        self.assertEqual(op1.e1, 0)
        self.assertEqual(op1.e2, 0)
        self.assertAlmostEqual(op1.xi1, 0.35438003143365127, delta=XI_TOL)
        self.assertAlmostEqual(op1.xi2, 0.14232301622713658, delta=XI_TOL)
        assertAlmostEqualList(self, [0.1700422367869369, -0.9854367750944224, 0.0], p1t, delta=X_TOL)

        p2, op2, p2x, p2t, p2bdy = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.7, 0.6), surf2.createPositionProportion(0.5, 0.5))
        self.assertEqual(p2.e1, 0)
        self.assertEqual(p2.e2, 0)
        self.assertAlmostEqual(p2.xi1, 0.7455699485954836, delta=XI_TOL)
        self.assertAlmostEqual(p2.xi2, 0.6652325809349918, delta=XI_TOL)
        self.assertEqual(op2.e1, 0)
        self.assertEqual(op2.e2, 0)
        self.assertAlmostEqual(op2.xi1, 0.830980843977797, delta=XI_TOL)
        self.assertAlmostEqual(op2.xi2, 0.6652325809312423, delta=XI_TOL)
        assertAlmostEqualList(self, [-0.8391083674789601, 0.5439642889260237, 0.0], p2t, delta=X_TOL)

        p3, op3, p3x, p3t, p3bdy = surf1.findIntersectionPoint(
            surf2, surf1.createPositionProportion(0.0, 0.326), surf2.createPositionProportion(0.0, 0.25))
        self.assertEqual(p3.e1, 0)
        self.assertEqual(p3.e2, 0)
        self.assertAlmostEqual(p3.xi1, 0.1, delta=XI_TOL)
        self.assertAlmostEqual(p3.xi2, 0.3263518272375643, delta=XI_TOL)
        self.assertEqual(op3.e1, 0)
        self.assertEqual(op3.e2, 0)
        self.assertAlmostEqual(op3.xi1, 0.0, delta=XI_TOL)
        self.assertAlmostEqual(op3.xi2, 0.3263518246534994, delta=XI_TOL)
        assertAlmostEqualList(self, [1.0, 0.0, 0.0], p3t, delta=X_TOL)

        p4, op4, p4x, p4t, p4bdy = surf1.findIntersectionPoint(
            surf2, TrackSurfacePosition(0, 0, 0.35032948194631897, 0.0),
            TrackSurfacePosition(0, 0, 0.35431492084144983, 0.006191352063889618))
        self.assertEqual(p4.e1, 0)
        self.assertEqual(p4.e2, 0)
        self.assertAlmostEqual(p4.xi1, 0.38097836477185065, delta=XI_TOL)
        self.assertAlmostEqual(p4.xi2, 0.0, delta=XI_TOL)
        self.assertEqual(op4.e1, 0)
        self.assertEqual(op4.e2, 0)
        self.assertAlmostEqual(op4.xi1, 0.3528986153613523, delta=XI_TOL)
        self.assertAlmostEqual(op4.xi2, 0.0, delta=XI_TOL)
        assertAlmostEqualList(self, [-0.1713415101967921, -0.9852116964812602, 0.0], p4t, delta=X_TOL)

        cx, cd1, cprops, loop = surf1.findIntersectionCurve(
            surf2, surf1.createPositionProportion(0.25, 0.1), MAX_MAG_DXI=0.2)
        self.assertEqual(len(cx), 9)
        self.assertFalse(loop)
        clength = getCubicHermiteCurvesLength(cx, cd1)
        self.assertAlmostEqual(clength, 0.5166105006376641, delta=X_TOL)
        assertAlmostEqualList(self, [0.1000000000000001, 0.3263518562674236, 0.0], cx[0], delta=X_TOL)
        assertAlmostEqualList(self, [0.3373116664675114, 0.2471498514765583, 0.0], cx[4], delta=X_TOL)
        assertAlmostEqualList(self, [0.3809786164663428, 0.0, 0.0], cx[8], delta=X_TOL)
        assertAlmostEqualList(self, [0.1, 0.3263518562674236], cprops[0], delta=XI_TOL)
        assertAlmostEqualList(self, [0.3373116664675114, 0.2471498514765583], cprops[4], delta=XI_TOL)
        assertAlmostEqualList(self, [0.3809786164663428, 0.0], cprops[8], delta=XI_TOL)

        dx, dd1, dprops, loop = surf1.findIntersectionCurve(
            surf2, surf1.createPositionProportion(0.85, 0.75), MAX_MAG_DXI=0.2)
        self.assertEqual(len(dx), 9)
        self.assertFalse(loop)
        dlength = getCubicHermiteCurvesLength(dx, dd1)
        self.assertAlmostEqual(dlength, 0.5433209057545304, delta=X_TOL)
        assertAlmostEqualList(self, [0.9, 0.625, 0.0], dx[0], delta=X_TOL)
        assertAlmostEqualList(self, [0.6714699015724952, 0.7450512483881269, 0.0], dx[4], delta=3.0*X_TOL)
        assertAlmostEqualList(self, [0.5797035235769858, 1.0, 0.0], dx[8], delta=3.0*X_TOL)
        assertAlmostEqualList(self, [0.9, 0.625], dprops[0], delta=XI_TOL)
        assertAlmostEqualList(self, [0.6714699015724951, 0.7450512483881269], dprops[4], delta=3.0*XI_TOL)
        assertAlmostEqualList(self, [0.5797035235769857, 1.0], dprops[8], delta=3.0*XI_TOL)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # surfaceGroupName = "surface"
        # surf1.generateMesh(region, group_name=surfaceGroupName)
        # surf2.generateMesh(region, group_name=surfaceGroupName)
        # coordinateFieldName = "curve_coordinates"
        # curveGroupName = "curve"
        # fieldmodule = region.getFieldmodule()
        # curveCoordinates = find_or_create_field_coordinates(fieldmodule, coordinateFieldName, managed=True)
        # curveGroup = find_or_create_field_group(fieldmodule, curveGroupName)
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(curveCoordinates)
        # nodetemplate.setValueNumberOfVersions(curveCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # curveNodesetGroup = curveGroup.getOrCreateNodesetGroup(nodes)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [p1x, p2x, p3x, p4x]
        # pd1 = [p1t, p2t, p3t, p4t]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n])
        #     curveNodesetGroup.addNode(node)
        # generateCurveMesh(region, cx, cd1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)
        # generateCurveMesh(region, dx, dd1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)

    def test_tube_intersections1(self):
        """
        Test tube intersections in a diverging bifurcation with one pair of tubes equal sized and continuous,
        but with one small outward tube which only intersects with the first tube.
        """
        elementsCountAround = 8
        elementsCountAlong = 6
        path1Params = [
            [[0.000, 0.000, 0.000], [1.000, 0.000, 0.000]],
            [[0.979, 0.123, 0.000], [1.010, -0.124, 0.000]],
            [[-0.031, 0.248, 0.000], [0.031, 0.348, 0.000]],
            [[0.063, 0.007, 0.000], [0.060, 0.307, 0.000]],
            [[0.000, 0.000, 0.250], [0.000, 0.000, 0.413]],
            [[0.000, 0.000, 0.126], [0.000, 0.000, 0.601]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path1Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube1Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        path2Params = [
            [[1.000, 0.000, 0.000], [2.000, -0.250, 0.000]],
            [[1.010, -0.124, 0.000], [0.980, -0.372, 0.000]],
            [[0.031, 0.348, 0.000], [0.089, 0.234, 0.000]],
            [[0.060, 0.307, 0.000], [-0.013, -0.101, -0.031]],
            [[0.000, 0.000, 0.413], [-0.000, 0.000, 0.250]],
            [[0.000, 0.000, 0.601], [-0.023, -0.019, -0.111]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path2Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube2Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        path3Params = [
            [[1.000, 0.000, 0.000], [0.488, 0.731, 0.107]],
            [[-0.512, 0.731, 0.107], [-0.512, 0.731, 0.107]],
            [[-0.034, -0.015, -0.055], [-0.052, -0.024, -0.088]],
            [[0.000, 0.000, 0.000], [-0.053, -0.053, -0.046]],
            [[-0.068, -0.056, 0.057], [-0.065, -0.054, 0.053]],
            [[0.000, 0.000, 0.000], [-0.086, -0.042, -0.109]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path3Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube3Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        startPosition = TrackSurfacePosition(12, 5, 0.8242639579553614, 0.999980904098904)
        p1x = tube1Surface.evaluateCoordinates(startPosition)
        otherStartPosition = TrackSurfacePosition(4, 0, 0.1580888007016199, 4.3205057553252083e-16)
        p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(tube2Surface, startPosition, otherStartPosition)
        self.assertEqual(p2.e1, 12)
        self.assertEqual(p2.e2, 5)
        self.assertAlmostEqual(p2.xi1, 0.8242639579553614, delta=XI_TOL)
        self.assertAlmostEqual(p2.xi2, 1.0, delta=XI_TOL)
        self.assertEqual(op2.e1, 4)
        self.assertEqual(op2.e2, 0)
        self.assertAlmostEqual(op2.xi1, 0.8242639579553614, delta=XI_TOL)
        self.assertAlmostEqual(op2.xi2, 0.0, delta=XI_TOL)
        assertAlmostEqualList(self, [0.045070823041903824, 0.5059563360833004, -0.8613807560466439], p2t, delta=X_TOL)

        # get loop intersection of fully connected tube1 and tube2
        ax, ad1, aprops, aloop = tube1Surface.findIntersectionCurve(tube2Surface, curveElementsCount=12)
        self.assertEqual(len(ax), 12)
        self.assertTrue(aloop)
        aCircumference = getCubicHermiteCurvesLength(ax, ad1, loop=True)
        self.assertAlmostEqual(aCircumference, 2.3973806470916053, delta=X_TOL)
        assertAlmostEqualList(self, [0.9710855559899819, -0.3245879521124608, 0.14882333760790006], ax[0], delta=X_TOL)
        assertAlmostEqualList(self, [1.0042792217936645, 0.04803771561920248, -0.4088273215322407], ax[4], delta=X_TOL)
        assertAlmostEqualList(self, [1.0253800137427358, 0.28491112201522784, 0.23710643763851688], ax[8], delta=X_TOL)
        assertAlmostEqualList(self, [1.437028796684921, 1.0], aprops[0], delta=XI_TOL)
        assertAlmostEqualList(self, [1.7702955198293207, 1.0], aprops[4], delta=XI_TOL)
        assertAlmostEqualList(self, [0.10347469710521648, 1.0], aprops[8], delta=XI_TOL)

        # get loop intersection of unconnected tube2 and tube3
        bx, bd1, bprops, bloop = tube2Surface.findIntersectionCurve(tube3Surface, curveElementsCount=12)
        self.assertIsNone(bx)
        self.assertIsNone(bd1)
        self.assertIsNone(bprops)
        self.assertFalse(bloop)

        # get loop intersection of large tube1 and small tube3
        cx, cd1, cprops, cloop = tube1Surface.findIntersectionCurve(
             tube3Surface, tube1Surface.createPositionProportion(0.1, 0.7))
        self.assertEqual(len(cx), 8)
        self.assertTrue(cloop)
        cCircumference = getCubicHermiteCurvesLength(cx, cd1, loop=True)
        self.assertAlmostEqual(cCircumference, 0.587857727905694, delta=X_TOL)
        assertAlmostEqualList(self, [0.7318758089726128, 0.28115396634786666, 0.13238674558014327], cx[0], delta=X_TOL)
        assertAlmostEqualList(self, [0.8408765202354725, 0.3182280500580874, -0.04549438718340043], cx[4], delta=X_TOL)
        assertAlmostEqualList(self, [1.0762363614295223, 0.7158314515175085], cprops[0], delta=XI_TOL)
        assertAlmostEqualList(self, [0.9766132119015503, 0.8198673158401282], cprops[4], delta=XI_TOL)

        # make trimmed tube3 starting at intersection with tube1
        tx, td1, td2, td12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong,
                                                     startSurface=tube1Surface)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(tx)):
            nx += tx[i]
            nd1 += td1[i]
            nd2 += td2[i]
            nd12 += td12[i]
        tube3TrimmedSurface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)
        tCircumference = getCubicHermiteCurvesLength(tx[0], td1[0], loop=True)
        self.assertAlmostEqual(tCircumference, 0.589158910546622, delta=X_TOL)
        tLength = getCubicHermiteCurvesLength([tx[n][0] for n in range(elementsCountAlong + 1)],
                                              [td2[n][0] for n in range(elementsCountAlong + 1)])
        self.assertAlmostEqual(tLength, 0.5004154200664181, delta=X_TOL)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # surfaceGroupName = "surface"
        # tube1Surface.generateMesh(region, group_name=surfaceGroupName)
        # tube2Surface.generateMesh(region, group_name=surfaceGroupName)
        # # tube3Surface.generateMesh(region, group_name=surfaceGroupName)
        # tube3TrimmedSurface.generateMesh(region, group_name=surfaceGroupName)
        # coordinateFieldName = "curve_coordinates"
        # curveGroupName = "curve"
        # fieldmodule = region.getFieldmodule()
        # curveCoordinates = find_or_create_field_coordinates(fieldmodule, coordinateFieldName, managed=True)
        # curveGroup = find_or_create_field_group(fieldmodule, curveGroupName)
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(curveCoordinates)
        # nodetemplate.setValueNumberOfVersions(curveCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # curveNodesetGroup = curveGroup.getOrCreateNodesetGroup(nodes)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [p1x, p2x]
        # pd1 = [[0.0, 0.0, 0.0], p2t]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n])
        #     curveNodesetGroup.addNode(node)
        # generateCurveMesh(region, ax, ad1, aloop, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)
        # generateCurveMesh(region, cx, cd1, cloop, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)

    def test_tube_intersections1_coarse(self):
        """
        Same as above but fewer along.
        """
        elementsCountAround = 8
        elementsCountAlong = 4
        path1Params = [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.9899453611864429, 0.1306534639262169, 0.1208354837258566], [0.986629229688665, -0.1356318297362336, -0.140819493622476]],
            [[-0.03223527987199958, 0.2478818360588229, -0.003934727904278583], [0.0333737119177766, 0.2477165835464935, -0.004763358991550144]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[-0.03029087281311803, 0.0, 0.2481581411604695], [0.03532398595224329, 8.623473827193372e-19, 0.2474918504041006]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path1Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube1Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        path2Params = [
            [[1.0, 0.0, 0.0], [1.938402743610923, -0.5417378094686381, -0.3033747813098204]],
            [[0.7703805269668589, -0.789032983143896, -0.3738697910125267], [1.133060151191206, -0.1265084623478793, 0.2366603646370917]],
            [[0.1645877843296779, 0.1824773036586876, -0.04596623651010114], [0.0291111109020875, 0.2482090649360781, -0.006693527142258089]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0897382556035165, -0.02243456390087916, 0.2322579079898364], [-0.04972023063292857, 0.01243005765823214, 0.2446904009813655]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path2Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube2Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        XI_TOL = 1.0E-6

        startPosition = TrackSurfacePosition(6, 3, 0.5717869946250334, 1.0)
        otherStartPosition = TrackSurfacePosition(6, 0, 0.11046416298376727, 0.07944602603404333)
        p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(tube2Surface, startPosition, otherStartPosition)
        self.assertEqual(p2.e1, 6)
        self.assertEqual(p2.e2, 3)
        self.assertAlmostEqual(p2.xi1, 0.4639613077782414, delta=XI_TOL)
        self.assertAlmostEqual(p2.xi2, 1.0, delta=XI_TOL)
        self.assertEqual(op2.e1, 6)
        self.assertEqual(op2.e2, 0)
        self.assertAlmostEqual(op2.xi1, 0.4873273111526366, delta=XI_TOL)
        self.assertAlmostEqual(op2.xi2, 0.0032771820253160657, delta=XI_TOL)
        self.assertIsNone(p2t)

    def test_tube_intersections2(self):
        """
        Test tube intersections in a diverging bifurcation case with similar tube sizes.
        """
        elementsCountAround = 8
        elementsCountAlong = 6
        path1Params = [
            [[0.000, 0.000, 0.000], [1.000, 0.000, 0.000]],
            [[0.990, 0.131, 0.121], [0.987, -0.136, -0.141]],
            [[-0.032, 0.248, -0.004], [0.033, 0.248, -0.005]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]],
            [[-0.030, 0.000, 0.248], [0.035, 0.000, 0.247]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path1Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube1Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        path2Params = [
            [[1.000, 0.000, 0.000], [1.938, -0.542, -0.303]],
            [[0.770, -0.789, -0.374], [1.133, -0.127, 0.237]],
            [[0.165, 0.182, -0.046], [0.029, 0.248, -0.007]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]],
            [[0.090, -0.022, 0.232], [-0.050, 0.012, 0.245]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path2Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube2Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        path3Params = [
            [[1.000, 0.000, 0.000], [1.658, 0.634, 0.017]],
            [[0.795, 0.473, -0.005], [0.441, 0.812, -0.048]],
            [[-0.128, 0.215, 0.000], [-0.216, 0.120, 0.037]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]],
            [[0.001, 0.000, 0.250], [0.038, -0.006, 0.247]],
            [[0.000, 0.000, 0.000], [0.000, 0.000, 0.000]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path3Params, elementsCountAround)
        sx, sd1, sd2, sd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(sx)):
            nx += sx[i]
            nd1 += sd1[i]
            nd2 += sd2[i]
            nd12 += sd12[i]
        tube3Surface = TrackSurface(elementsCountAround, elementsCountAlong, nx, nd1, nd2, nd12, loop1=True)

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        p1, op1, p1x, p1t, p1bdy = tube1Surface.findIntersectionPoint(
            tube2Surface, tube1Surface.createPositionProportion(0.25, 1.0),
            tube2Surface.createPositionProportion(0.25, 0.0))
        self.assertIsNotNone(p1)
        if elementsCountAround == 8:
            self.assertEqual(p1.e1, 2)
            self.assertEqual(p1.e2, 5)
            self.assertEqual(p1bdy, 2)
            self.assertAlmostEqual(p1.xi1, 0.551275863396822, delta=XI_TOL)
            self.assertAlmostEqual(p1.xi2, 1.0, delta=XI_TOL)
            self.assertEqual(op1.e1, 2)
            self.assertEqual(op1.e2, 0)
            self.assertAlmostEqual(op1.xi1, 0.5566551590679389, delta=XI_TOL)
            self.assertAlmostEqual(op1.xi2, 0.052773843347578696, delta=XI_TOL)
            assertAlmostEqualList(self, [-0.6607544357449996, -0.7021257304258945, -0.26537338660312426], p1t, delta=X_TOL)

        if elementsCountAround == 8:
            startPosition = TrackSurfacePosition(6, 5, 0.0, 1.0)
            p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(
                tube3Surface, startPosition, TrackSurfacePosition(6, 0, 0.0, 0.0))
            self.assertEqual(p2.e1, 6)
            self.assertEqual(p2.e2, 5)
            self.assertAlmostEqual(p2.xi1, 0.33908360205244037, delta=XI_TOL)
            self.assertAlmostEqual(p2.xi2, 0.899262977194283, delta=XI_TOL)
            self.assertEqual(op2.e1, 6)
            self.assertEqual(op2.e2, 0)
            self.assertAlmostEqual(op2.xi1, 0.4059622427354084, delta=XI_TOL)
            self.assertAlmostEqual(op2.xi2, 0.0, delta=XI_TOL)
            assertAlmostEqualList(self, [0.1412180735660879, 0.9590457759564066, 0.24553748251225724], p2t, delta=X_TOL)

            startPosition = TrackSurfacePosition(6, 5, 0.053597696166878706, 1.0)
            p3, op3, p3x, p3t, p3bdy = tube1Surface.findIntersectionPoint(
                tube3Surface, startPosition, TrackSurfacePosition(6, 0, 0.5777545444894798, 0.1230948585475792))
            self.assertEqual(p3.e1, 6)
            self.assertEqual(p3.e2, 5)
            self.assertAlmostEqual(p3.xi1, 0.3390958322811555, delta=XI_TOL)
            self.assertAlmostEqual(p3.xi2, 0.899250982796703, delta=XI_TOL)

            startPosition = TrackSurfacePosition(5, 5, 0.9639270978611663, 0.7487688426600227)
            p4, op4, p4x, p4t, p4bdy = tube1Surface.findIntersectionPoint(
                tube2Surface, startPosition, TrackSurfacePosition(5, 0, 0.49953899323872, 0.27331994814227073))
            self.assertEqual(p4.e1, 5)
            self.assertEqual(p4.e2, 5)
            self.assertAlmostEqual(p4.xi1, 0.9585925958894155, delta=XI_TOL)
            self.assertAlmostEqual(p4.xi2, 0.7720200484294564, delta=XI_TOL)

            startPosition = TrackSurfacePosition(6, 0, 0.27152105049897823,0.0)
            p5, op5, p5x, p5t, p5bdy = tube2Surface.findIntersectionPoint(
                tube3Surface, startPosition, TrackSurfacePosition(6, 0, 0.1648182942983345,0.4437057086482662))
            self.assertEqual(p5.e1, 6)
            self.assertEqual(p5.e2, 0)
            self.assertAlmostEqual(p5.xi1, 0.3265731742895994, delta=XI_TOL)
            self.assertAlmostEqual(p5.xi2, 0.06905957836045906, delta=XI_TOL)

            startPosition = TrackSurfacePosition(5, 5, 0.8414377968474449, 0.9143600731749943)
            otherStartPosition = TrackSurfacePosition(6, 0, 0.4076509361724341, 0.0017373594850997659)
            p6, op6, p6x, p6t, p6bdy = tube1Surface.findIntersectionPoint(tube3Surface, startPosition, otherStartPosition)
            self.assertEqual(p6.e1, 6)
            self.assertEqual(p6.e2, 5)
            self.assertAlmostEqual(p6.xi1, 0.3390904086633162, delta=XI_TOL)
            self.assertAlmostEqual(p6.xi2, 0.89925630182541039, delta=XI_TOL)

            startPosition = TrackSurfacePosition(2, 0, 0.6590628632288293, 0.0)
            otherStartPosition = TrackSurfacePosition(2, 0, 0.5297149091022884, 0.09186286116765788)
            p7, op7, p7x, p7t, p7bdy = tube2Surface.findIntersectionPoint(tube3Surface, startPosition, otherStartPosition)
            self.assertEqual(p7.e1, 2)
            self.assertEqual(p7.e2, 0)
            self.assertAlmostEqual(p7.xi1, 0.3488485362030178, delta=XI_TOL)
            self.assertAlmostEqual(p7.xi2, 0.10171143282929856, delta=XI_TOL)

        if elementsCountAround == 4:
            # startPosition = TrackSurfacePosition(3, 5, 0.36148219241126434, 1.0)
            # p1x = tube1Surface.evaluateCoordinates(startPosition)
            # p1t = [0.0, 0.0, 0.0]
            # p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(
            #     tube2Surface, startPosition, TrackSurfacePosition(2, 0, 0.9794494698544911, 0.1562967852059411))
            # self.assertIsNotNone(p2)
            # self.assertEqual(p2.e1, 6)
            # self.assertEqual(p2.e2, 5)
            # self.assertAlmostEqual(p2.xi1, 0.33908360205244037, delta=XI_TOL)
            # self.assertAlmostEqual(p2.xi2, 0.899262977194283, delta=XI_TOL)
            # self.assertEqual(op2.e1, 6)
            # self.assertEqual(op2.e2, 0)
            # self.assertAlmostEqual(op2.xi1, 0.4059622427354084, delta=XI_TOL)
            # self.assertAlmostEqual(op2.xi2, 0.0, delta=XI_TOL)
            # assertAlmostEqualList(self, [0.1412180735660879, 0.9590457759564066, 0.24553748251225724], p2t, delta=X_TOL)

            # startPosition = TrackSurfacePosition(1, 5, 0.3631033157903245, 1.0)
            # p1x = tube1Surface.evaluateCoordinates(startPosition)
            # p1t = [0.0, 0.0, 0.0]
            # p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(
            #     tube2Surface, startPosition, TrackSurfacePosition(1, 0, 0.33228274550722303,0.0))
            # self.assertIsNotNone(p2)

            # startPosition = TrackSurfacePosition(1, 5, 0.3631033157903245, 1.0)
            # p1x = tube1Surface.evaluateCoordinates(startPosition)
            # p1t = [0.0, 0.0, 0.0]
            # p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(
            #     tube3Surface, startPosition, TrackSurfacePosition(1, 0, 0.33228274550722303, 0.0))
            # self.assertIsNotNone(p2)

            startPosition = TrackSurfacePosition(3, 0, 0.0865879388733597, 0.0)
            p1x = tube1Surface.evaluateCoordinates(startPosition)
            p1t = [0.0, 0.0, 0.0]
            p2, op2, p2x, p2t, p2bdy = tube2Surface.findIntersectionPoint(
                tube3Surface, startPosition, TrackSurfacePosition(3, 0, 0.19438867305821184, 0.0))
            self.assertIsNotNone(p2)

        # get non-loop intersection of tube1 and tube2
        ax, ad1, aprops, aloop = tube1Surface.findIntersectionCurve(tube2Surface)
        self.assertEqual(len(ax), 9)
        self.assertFalse(aloop)
        aLength = getCubicHermiteCurvesLength(ax, ad1, loop=True)
        if elementsCountAround == 8:
            self.assertAlmostEqual(aLength, 1.291896441648134, delta=X_TOL)
            assertAlmostEqualList(self, [1.0179538858287902, -0.10372009294261042, 0.22619317098715946], ax[0], delta=X_TOL)
            assertAlmostEqualList(self, [0.8745386899334123, -0.21894738318303653, -0.07777014486867816], ax[4], delta=X_TOL)
            assertAlmostEqualList(self, [0.9779282561711389, 0.0820751155699032, -0.23449126177600493], ax[8], delta=X_TOL)
            assertAlmostEqualList(self, [0.31891164977233744, 1.0], aprops[0], delta=XI_TOL)
            assertAlmostEqualList(self, [0.5615258901138819, 0.9099206438856747], aprops[4], delta=XI_TOL)
            assertAlmostEqualList(self, [0.8039061628699322, 1.0], aprops[8], delta=XI_TOL)
        elif elementsCountAround == 6:
            self.assertAlmostEqual(aLength, 1.2935036652991265, delta=X_TOL)
        elif elementsCountAround == 4:
            self.assertAlmostEqual(aLength, 1.2777849511330288, delta=X_TOL)

        # # get non-loop intersection of tube1 and tube3
        bx, bd1, bprops, bloop = tube1Surface.findIntersectionCurve(tube3Surface)
        self.assertEqual(len(bx), 9)
        self.assertFalse(bloop)
        bLength = getCubicHermiteCurvesLength(bx, bd1, loop=True)
        if elementsCountAround == 8:
            self.assertAlmostEqual(bLength, 1.2644609503534696, delta=X_TOL)
            assertAlmostEqualList(self, [0.959094526660969, 0.06711431057004068, -0.2372654694905775], bx[0], delta=X_TOL)
            assertAlmostEqualList(self, [0.9486154683235082, 0.2545005003987398, 0.05281998774919236], bx[4], delta=X_TOL)
            assertAlmostEqualList(self, [1.0147890884040247, -0.023171620406177137, 0.24848738902991338], bx[8], delta=X_TOL)
            assertAlmostEqualList(self, [0.7923873358701699, 0.9832080394698671], bprops[0], delta=XI_TOL)
            assertAlmostEqualList(self, [1.0300414330440506, 0.9150675749144951], bprops[4], delta=XI_TOL)
            assertAlmostEqualList(self, [1.2664762695247593, 0.9841126335511698], bprops[8], delta=XI_TOL)
        elif elementsCountAround == 6:
            self.assertAlmostEqual(bLength, 1.2594725877612785, delta=X_TOL)
        elif elementsCountAround == 4:
            self.assertAlmostEqual(bLength, 1.2658594359610489, delta=X_TOL)

        # get non-loop intersection of tube2 and tube3
        cx, cd1, cprops, cloop = tube2Surface.findIntersectionCurve(tube3Surface)
        self.assertEqual(len(cx), 9)
        self.assertFalse(cloop)
        cLength = getCubicHermiteCurvesLength(cx, cd1, loop=True)
        if elementsCountAround == 8:
            self.assertAlmostEqual(cLength, 1.5064021938074896, delta=X_TOL)
            assertAlmostEqualList(self, [1.0515320702829902, -0.08500706675336937, 0.22935825287467218], cx[0], delta=X_TOL)
            assertAlmostEqualList(self, [1.371362500857284, -0.024407127528850117, -0.08285047871256904], cx[4], delta=X_TOL)
            assertAlmostEqualList(self, [0.9646355114600913, 0.05778940350574679, -0.240589269211589], cx[8], delta=X_TOL)
            assertAlmostEqualList(self, [0.29558389660701956, 0.016942064305594088], cprops[0], delta=XI_TOL)
            assertAlmostEqualList(self, [0.05134644399148322, 0.29537544509373714], cprops[4], delta=XI_TOL)
            assertAlmostEqualList(self, [0.7907974706501059, 0.011505022662710312], cprops[8], delta=XI_TOL)
        elif elementsCountAround == 6:
            self.assertAlmostEqual(cLength, 1.500631828350329, delta=X_TOL)
        elif elementsCountAround == 4:
            self.assertAlmostEqual(cLength, 1.5116074716918804, delta=X_TOL)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # surfaceGroupName = "surface"
        # tube1Surface.generateMesh(region, group_name=surfaceGroupName)
        # tube2Surface.generateMesh(region, group_name=surfaceGroupName)
        # tube3Surface.generateMesh(region, group_name=surfaceGroupName)
        # coordinateFieldName = "curve_coordinates"
        # curveGroupName = "curve"
        # fieldmodule = region.getFieldmodule()
        # curveCoordinates = find_or_create_field_coordinates(fieldmodule, coordinateFieldName, managed=True)
        # curveGroup = find_or_create_field_group(fieldmodule, curveGroupName)
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(curveCoordinates)
        # nodetemplate.setValueNumberOfVersions(curveCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # curveNodesetGroup = curveGroup.getOrCreateNodesetGroup(nodes)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [p1x, p2x, p3x, p4x, p5x, p6x, p7x]
        # pd1 = [p1t, p2t, p3t, p4t, p5t, p6t, p7t]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n])
        #     curveNodesetGroup.addNode(node)
        # generateCurveMesh(region, ax, ad1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)
        # generateCurveMesh(region, bx, bd1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)
        # generateCurveMesh(region, cx, cd1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)

    def test_tube_intersections3(self):
        """
        Regular converging bifurcation failed to find boundary intersection point.
        """
        elementsCountAround = 10
        # elementsCountAlong = 4
        path1Params = [
            [[0.0, -0.25, 0.0], [1.0, 0.0, 0.0]],
            [[0.9701425001453319, 0.24253562503633297, 0.0], [0.9701425001453319, 0.24253562503633297, 0.0]],
            [[-0.06063390625908324, 0.24253562503633297, 0.0], [-0.06063390625908324, 0.24253562503633297, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, -0.0, 0.25], [0.0, -0.0, 0.25]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path1Params, elementsCountAround)
        # px, pd1, pd2, pd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(px)):
            nx += px[i]
            nd1 += pd1[i]
            nd2 += pd2[i]
            nd12 += pd12[i]
        tube1Surface = TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True)

        path2Params = [
            [[0.0, 0.25, 0.0], [1.0, 0.0, 0.0]],
            [[0.9701425001453319, -0.24253562503633297, 0.0], [0.9701425001453319, -0.24253562503633297, 0.0]],
            [[0.06063390625908324, 0.24253562503633297, 0.0], [0.06063390625908324, 0.24253562503633297, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[-0.0, 0.0, 0.25], [-0.0, 0.0, 0.25]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(path2Params, elementsCountAround)
        # px, pd1, pd2, pd12 = resampleTubeCoordinates((px, pd1, pd2, pd12), elementsCountAlong)
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(px)):
            nx += px[i]
            nd1 += pd1[i]
            nd2 += pd2[i]
            nd12 += pd12[i]
        tube2Surface = TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True)

        # path3Params = [
        #     [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        #     [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        #     [[-0.0, 0.25, 0.0], [-0.0, 0.25, 0.0]],
        #     [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        #     [[0.0, -0.0, 0.25], [0.0, -0.0, 0.25]],
        #     [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]

        XI_TOL = 1.0E-6
        X_TOL = 1.0E-6

        startPosition = TrackSurfacePosition(12, 0, 0.7886706262895231, 1.0)
        p1x = tube1Surface.evaluateCoordinates(startPosition)
        otherStartPosition = TrackSurfacePosition(2, 0, 0.636882069585452, 0.9147548823734029)
        p2, op2, p2x, p2t, p2bdy = tube1Surface.findIntersectionPoint(
            tube2Surface, startPosition, otherStartPosition)
        self.assertEqual(p2.e1, 12)
        self.assertEqual(p2.e2, 0)
        self.assertAlmostEqual(p2.xi1, 0.4973473706811795, delta=XI_TOL)
        self.assertAlmostEqual(p2.xi2, 1.0, delta=XI_TOL)
        self.assertEqual(op2.e1, 2)
        self.assertEqual(op2.e2, 0)
        self.assertAlmostEqual(op2.xi1, 0.49765944425751174, delta=XI_TOL)
        self.assertAlmostEqual(op2.xi2, 0.9998039549126785, delta=XI_TOL)
        self.assertIsNone(p2t)

        startPosition = TrackSurfacePosition(0, 0, 0.08468387605632358, 0.0)
        p3x = tube1Surface.evaluateCoordinates(startPosition)
        otherStartPosition = TrackSurfacePosition(4, 0, 0.42675444143171326, 0.09458268381308738)
        p4, op4, p4x, p4t, p4bdy = tube1Surface.findIntersectionPoint(
            tube2Surface, startPosition, otherStartPosition)
        self.assertEqual(p4.e1, 0)
        self.assertEqual(p4.e2, 0)
        self.assertAlmostEqual(p4.xi1, 0.001112953121053552, delta=XI_TOL)
        self.assertAlmostEqual(p4.xi2, 0.030691469925266602, delta=XI_TOL)
        assertAlmostEqualList(self, [0.002978246340503309, 1.6957770277963287e-12, 0.999995565014533], p4t, delta=X_TOL)

        # get non-loop intersection of tube1 and tube2
        ax, ad1, aprops, aloop = tube1Surface.findIntersectionCurve(tube2Surface)
        self.assertEqual(len(ax), 9)
        self.assertFalse(aloop)
        aLength = getCubicHermiteCurvesLength(ax, ad1, loop=True)
        self.assertAlmostEqual(aLength, 2.7550717460511067, delta=X_TOL)
        assertAlmostEqualList(self, [0.9999027578036812, 0.00038896878527516776, -0.24987413008984302], ax[0], delta=X_TOL)
        assertAlmostEqualList(self, [-0.03077634623170676, 1.5966942528521302e-09, 7.945917568834134e-05], ax[4], delta=X_TOL)
        assertAlmostEqualList(self, [0.9997557619375088, 0.0009769522499643996, 0.24987245334458275], ax[8], delta=X_TOL)
        assertAlmostEqualList(self, [0.7502552895299374, 1.0], aprops[0], delta=XI_TOL)
        assertAlmostEqualList(self, [1.0000505870312695, 0.030691274736817955], aprops[4], delta=XI_TOL)
        assertAlmostEqualList(self, [1.249358801390615, 1.0], aprops[8], delta=XI_TOL)

        # context = Context("TrackSurface")
        # region = context.getDefaultRegion()
        # surfaceGroupName = "surface"
        # tube1Surface.generateMesh(region, group_name=surfaceGroupName)
        # tube2Surface.generateMesh(region, group_name=surfaceGroupName)
        # coordinateFieldName = "curve_coordinates"
        # curveGroupName = "curve"
        # fieldmodule = region.getFieldmodule()
        # curveCoordinates = find_or_create_field_coordinates(fieldmodule, coordinateFieldName, managed=True)
        # curveGroup = find_or_create_field_group(fieldmodule, curveGroupName)
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(curveCoordinates)
        # nodetemplate.setValueNumberOfVersions(curveCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # curveNodesetGroup = curveGroup.getOrCreateNodesetGroup(nodes)
        # fieldcache = fieldmodule.createFieldcache()
        # px = [p1x, p2x, p3x, p4x]
        # pd1 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], p4t]
        # for n in range(len(px)):
        #     node = nodes.createNode(-1, nodetemplate)
        #     fieldcache.setNode(node)
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, px[n])
        #     curveCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n])
        #     curveNodesetGroup.addNode(node)
        # generateCurveMesh(region, ax, ad1, coordinate_field_name=coordinateFieldName, group_name=curveGroupName)


if __name__ == "__main__":
    unittest.main()
