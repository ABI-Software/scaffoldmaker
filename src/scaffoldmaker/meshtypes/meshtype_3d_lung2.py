'''
Generates 3D lung surface mesh.
'''
import copy
import math
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, remapEftNodeValueLabelsVersion, setEftScaleFactorIds
from scaffoldmaker.utils.cylindermesh import createEllipsePerimeter
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.utils.zinc.general import ChangeManager
from scaffoldmaker.utils.geometry import createEllipsoidPoints, getEllipseRadiansToX, getEllipseArcLength, getApproximateEllipsePerimeter, \
    updateEllipseAngleByArcLength
from scaffoldmaker.utils.interpolation import DerivativeScalingMode
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.vector import magnitude, setMagnitude, crossproduct3, normalise
from opencmiss.utils.zinc.field import Field, findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from opencmiss.zinc.element import Element
from opencmiss.zinc.node import Node

class MeshType_3d_lung2(Scaffold_base):
    '''
    3D lung scaffold.
    '''

    """
        Generates an ellipsoid with a tear-shaped base for the lung mathematically,
         with x, y, z length and the position, angle of the oblique fissure. 
         Regions and markers of the lung are annotated.
    """

    @staticmethod
    def getName():
        return '3D Lung 2'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Human 1'

        options = {
            'Base parameter set': parameterSetName,
            'Length - Left/Right Lung': [0.5, 0.5],
            'Width - Left/Right Lung':  [0.25, 0.25],
            'Height - Left/Right Lung': [1.0, 1.0],
            'Distance from origin - Left/Right Lung': [0.5, 0.0, 0.0],
            'Fissure angle - Left/Right Lung': [45.0, 45.0],
            'Oblique proportion - Left/Right Lung': [0.8, 0.8],
            'Tilt apex along x-axis - Left/Right Lung': [0.0, 0.0],
            'Tilt apex along y-axis - Left/Right Lung': [0.0, 0.0],
            'Rotate around z-axis - Left/Right Lung': [0.0, 0.0],
            'Tilt diaphragm surface along x-axis - Left/Right Lung': [0.0, 0.0],
            'Tilt diaphragm surface along y-axis - Left/Right Lung': [0.0, 0.0],
            'Diaphragmatic curve radius - Left/Right Lung': [0.0, 0.0],
            'Bulge radius around y-axis - Left/Right Lung': [0.0, 0.0],
            'Bulge radius around z-axis - Left/Right Lung': [0.0, 0.0],
            'Medial curve radius - Left/Right Lung': [0.0, 0.0],
            'Sharpening edge - Left/Right Lung': [0.0, 0.0],
            'Tapering along z-axis - Left/Right Lung': [0.0, 0.0],
            'Open fissures - Left/Right Lung': False,
            'Length - Accessory lobe': 0.1,
            'Width - Accessory lobe': 0.25,
            'Height - Accessory lobe': 0.25,
            'Distance from origin - Accessory lobe': [0.0, -0.2, 0.0],
            'Refine': False,
            'Refine number of elements': 4
        }

        if 'Human 1' in parameterSetName:
            options['Distance from origin - Left/Right Lung'] = [0.35, 0.0, 0.0]
            options['Tilt apex along x-axis - Left/Right Lung'] = [5.0, 5.0]
            options['Tilt apex along y-axis - Left/Right Lung'] = [5.0, 5.0]
            options['Rotate around z-axis - Left/Right Lung'] = [20.0, -10.0]
            options['Tilt diaphragm surface along x-axis - Left/Right Lung'] = [15.0, 15.0]
            options['Diaphragmatic curve radius - Left/Right Lung'] = [1.0, 1.0]
            options['Bulge radius around z-axis - Left/Right Lung'] = [0.8, 0.8]
            options['Medial curve radius - Left/Right Lung'] = [2.0, 1.0]
            options['Sharpening edge - Left/Right Lung'] = [0.5, 0.0]
            options['Tapering along z-axis - Left/Right Lung'] = [0.2, 0.2]

        elif 'Mouse 1' in parameterSetName:
            options['Height - Left/Right Lung'] = [0.8, 0.8]
            options['Distance from origin - Left/Right Lung'] = [0.4, 0.0, 0.0]
            options['Tilt apex along x-axis - Left/Right Lung'] = [10, 10.0]
            options['Tilt apex along y-axis - Left/Right Lung'] = [5.0, 10.0]
            options['Rotate around z-axis - Left/Right Lung'] = [10.0, 20.0]
            options['Tilt diaphragm surface along x-axis - Left/Right Lung'] = [20.0, 20.0]
            options['Tilt diaphragm surface along y-axis - Left/Right Lung'] = [20.0, 10.0]
            options['Diaphragmatic curve radius - Left/Right Lung'] = [1.0, 1.0]
            options['Bulge radius around z-axis - Left/Right Lung'] = [0.6, 1.0]
            options['Medial curve radius - Left/Right Lung'] = [0.0, 2.0]
            options['Sharpening edge - Left/Right Lung'] = [1.8, 0.8]
            options['Tapering along z-axis - Left/Right Lung'] = [0.1, -0.2]
            options['Distance from origin - Accessory lobe'] = [0.0, -0.2, 0.5]

        elif 'Rat 1' in parameterSetName:
            options['Height - Left/Right Lung'] = [0.8, 0.8]
            options['Distance from origin - Left/Right Lung'] = [0.33, 0.0, 0.0]
            options['Tilt apex along x-axis - Left/Right Lung'] = [10, 10.0]
            options['Tilt apex along y-axis - Left/Right Lung'] = [15.0, 15.0]
            options['Rotate around z-axis - Left/Right Lung'] = [15.0, 15.0]
            options['Tilt diaphragm surface along x-axis - Left/Right Lung'] = [10.0, 10.0]
            options['Diaphragmatic curve radius - Left/Right Lung'] = [1.0, 1.0]
            options['Bulge radius around z-axis - Left/Right Lung'] = [0.6, 1.0]
            options['Medial curve radius - Left/Right Lung'] = [0.0, 2.0]
            options['Sharpening edge - Left/Right Lung'] = [1.6, 0.8]
            options['Tapering along z-axis - Left/Right Lung'] = [0.5, 0.5]
            options['Distance from origin - Accessory lobe'] = [0.0, -0.2, 0.5]

        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Length - Left/Right Lung',
            'Width - Left/Right Lung',
            'Height - Left/Right Lung',
            'Distance from origin - Left/Right Lung',
            'Fissure angle - Left/Right Lung',
            'Oblique proportion - Left/Right Lung',
            'Tilt apex along x-axis - Left/Right Lung',
            'Tilt apex along y-axis - Left/Right Lung',
            'Rotate around z-axis - Left/Right Lung',
            'Tilt diaphragm surface along x-axis - Left/Right Lung',
            'Tilt diaphragm surface along y-axis - Left/Right Lung',
            'Diaphragmatic curve radius - Left/Right Lung',
            'Bulge radius around y-axis - Left/Right Lung',
            'Bulge radius around z-axis - Left/Right Lung',
            'Medial curve radius - Left/Right Lung',
            'Sharpening edge - Left/Right Lung',
            'Tapering along z-axis - Left/Right Lung',
            'Open fissures - Left/Right Lung',
            'Length - Accessory lobe',
            'Width - Accessory lobe',
            'Height - Accessory lobe',
            'Distance from origin - Accessory lobe',
            'Refine',
            'Refine number of elements'
        ]
        return optionNames

    @classmethod
    def checkOptions(cls, options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False

        for i in range(2):
            if (options['Oblique proportion - Left/Right Lung'][i] < 0.7) or \
                    (options['Oblique proportion - Left/Right Lung'][i] > 0.99):
                options['Oblique proportion - Left/Right Lung'][i] = 0.8

        if options['Refine number of elements'] < 1:
            options['Refine number of elements'] = 1

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''

        parameterSetName = options['Base parameter set']
        isMouse = 'Mouse 1' in parameterSetName
        isHuman = 'Human 1' in parameterSetName
        isRat = 'Rat 1' in parameterSetName

        length = options['Length - Left/Right Lung']
        width = options['Width - Left/Right Lung']
        height = options['Height - Left/Right Lung']
        spaceFromCentre = options['Distance from origin - Left/Right Lung']
        length_1 = options['Length - Accessory lobe']
        width_1 = options['Width - Accessory lobe']
        height_1 = options['Height - Accessory lobe']
        spaceBelowCentre = options['Distance from origin - Accessory lobe']
        fissureAngle = options['Fissure angle - Left/Right Lung']
        obliqueProportion = options['Oblique proportion - Left/Right Lung']
        discontinuity = options['Open fissures - Left/Right Lung']
        tiltApex_xAxis = options['Tilt apex along x-axis - Left/Right Lung']
        tiltApex_yAxis = options['Tilt apex along y-axis - Left/Right Lung']
        rotate_ZAxis = options['Rotate around z-axis - Left/Right Lung']
        tiltDiap_xAxis = options['Tilt diaphragm surface along x-axis - Left/Right Lung']
        tiltDiap_yAxis = options['Tilt diaphragm surface along y-axis - Left/Right Lung']
        DiaphramaticCurveRadius = options['Diaphragmatic curve radius - Left/Right Lung']
        bulgeRadiusY = options['Bulge radius around y-axis - Left/Right Lung']
        bulgeRadiusZ = options['Bulge radius around z-axis - Left/Right Lung']
        medialCurveRadius = options['Medial curve radius - Left/Right Lung']
        sharpeningFactor = options['Sharpening edge - Left/Right Lung']
        taperingFactor = options['Tapering along z-axis - Left/Right Lung']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fm.findMeshByDimension(3)

        eftfactory = eftfactory_tricubichermite(mesh, None)
        eftRegular = eftfactory.createEftBasic()

        elementtemplateRegular = mesh.createElementtemplate()
        elementtemplateRegular.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateRegular.defineField(coordinates, -1, eftRegular)

        elementtemplateCustom = mesh.createElementtemplate()
        elementtemplateCustom.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Annotation groups
        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        annotationGroups = [leftLungGroup, lungGroup]
        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLungGroup)
        lowerRightLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of right lung"))
        lowerRightLungMeshGroup = lowerRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(lowerRightLungGroup)
        upperRightLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of right lung"))
        upperRightLungMeshGroup = upperRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(upperRightLungGroup)
        middleRightLungGroup = AnnotationGroup(region, get_lung_term("middle lobe of right lung"))
        middleRightLungMeshGroup = middleRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(middleRightLungGroup)
        mediastinumLeftGroup = AnnotationGroup(region, get_lung_term("anterior mediastinum of left lung"))
        mediastinumLeftGroupMeshGroup = mediastinumLeftGroup.getMeshGroup(mesh)
        annotationGroups.append(mediastinumLeftGroup)
        mediastinumRightGroup = AnnotationGroup(region, get_lung_term("anterior mediastinum of right lung"))
        mediastinumRightGroupMeshGroup = mediastinumRightGroup.getMeshGroup(mesh)
        annotationGroups.append(mediastinumRightGroup)

        # Nodeset group
        leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
        rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)

        # Marker points/groups
        leftApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_lung_term("apex of left lung"))
        rightApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_lung_term("apex of right lung"))
        leftDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_lung_term("dorsal base of left lung"))
        rightDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_lung_term("dorsal base of right lung"))
        leftVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_lung_term("ventral base of left lung"))
        rightVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_lung_term("ventral base of right lung"))
        rightLateralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_lung_term("laterodorsal tip of middle lobe of right lung"))
        leftMedialGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_lung_term("medial base of left lung"))
        rightMedialGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_lung_term("medial base of right lung"))

        if isHuman:
            # Annotation groups
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(upperLeftLungGroup)

        elif isMouse or isRat:
            # Annotation groups
            diaphragmaticLungGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
            diaphragmaticLungMeshGroup = diaphragmaticLungGroup.getMeshGroup(mesh)
            annotationGroups.append(diaphragmaticLungGroup)

            # Marker points
            accessoryApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                    get_lung_term("apex of right lung accessory lobe"))
            accessoryVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term("ventral base of right lung accessory lobe"))
            accessoryDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                      get_lung_term("dorsal base of right lung accessory lobe"))

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        cache = fm.createFieldcache()

        # common parameters in species
        leftLung = 0
        rightLung = 1

        # The number of the elements in the generic lungs
        # These counts are only values that work for nodeFieldParameters (KEEP THEM FIXED)
        lElementsCount1 = 2
        lElementsCount2 = 4
        lElementsCount3 = 3

        uElementsCount1 = 2
        uElementsCount2 = 4
        uElementsCount3 = 4

        if isHuman:
            # Create nodes
            nodeIndex = 0
            nodeIdentifier = 1
            lowerLeftNodeIds = []
            upperLeftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []

            # Left lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spaceFromCentre,
                length[0], width[0], height[0], fissureAngle[0], obliqueProportion[0],
                leftLung, cache, coordinates, nodes, nodetemplate, leftLungNodesetGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerLeftNodeIds, upperLeftNodeIds, nodeIndex, nodeIdentifier)

            # Right lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spaceFromCentre,
                length[1], width[1], height[1], fissureAngle[1], obliqueProportion[1],
                rightLung, cache, coordinates, nodes, nodetemplate, rightLungNodesetGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, nodeIndex, nodeIdentifier)

            # Create elements
            elementIdentifier = 1

            # Left lung elements
            elementIdentifier, leftUpperLobeElementID, leftLowerLobeElementID = \
                getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                leftLungMeshGroup, lowerLeftLungMeshGroup, None,
                upperLeftLungMeshGroup, mediastinumLeftGroupMeshGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerLeftNodeIds, upperLeftNodeIds, elementIdentifier)

            # Right lung elements
            elementIdentifier, rightUpperLobeElementID, rightLowerLobeElementID = \
                getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                rightLungMeshGroup, lowerRightLungMeshGroup, middleRightLungMeshGroup,
                upperRightLungMeshGroup, mediastinumRightGroupMeshGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, elementIdentifier)

            # Marker points
            lungNodesetGroup = lungGroup.getNodesetGroup(nodes)
            markerList = []

            lowerLeftElementCount = (lElementsCount1 * (lElementsCount2-1) * lElementsCount3 + lElementsCount1)

            idx = lowerLeftElementCount + (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount2)
            markerList.append({ "group" : leftApexGroup, "elementId" : idx, "xi" : [0.0, 1.0, 1.0]})

            idx = 1
            markerList.append({"group": leftDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 // 2)
            markerList.append({"group": leftMedialGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 - 1) + 1
            markerList.append({"group": leftVentralGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            upperLeftElementCount = (uElementsCount1 * uElementsCount2 * (uElementsCount3-1))
            leftLungElementCount = lowerLeftElementCount + upperLeftElementCount
            lowerRightElementCount = lowerLeftElementCount

            idx = leftLungElementCount + 1
            markerList.append({"group": rightDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = leftLungElementCount + lowerRightElementCount + \
                  (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount2)
            markerList.append({"group": rightApexGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})

            idx = leftLungElementCount + lElementsCount1 + 1
            markerList.append({"group": rightMedialGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            idx = leftLungElementCount + lElementsCount1 * lElementsCount2
            markerList.append({"group": rightVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * lElementsCount3) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 0.0, 1.0]})

        elif isMouse or isRat or test:
            # The number of the elements in the diaphragmatic animal lung
            diaphragmaticElementsCount1 = 2
            diaphragmaticElementsCount2 = 5
            diaphragmaticElementsCount3 = 2

            # Create nodes
            nodeIndex = 0
            nodeIdentifier = 1
            lowerLeftNodeIds = []
            upperLeftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []
            diaphragmaticNodeIds = []

            # Left lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spaceFromCentre,
                                                     length[0], width[0], height[0], fissureAngle[0], obliqueProportion[0],
                                                     leftLung, cache, coordinates,
                                                     nodes, nodetemplate, leftLungNodesetGroup,
                                                     lElementsCount1, lElementsCount2, lElementsCount3,
                                                     uElementsCount1, uElementsCount2, uElementsCount3,
                                                     lowerLeftNodeIds, upperLeftNodeIds, nodeIndex, nodeIdentifier)

            # Right lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spaceFromCentre,
                                                     length[1], width[1], height[1], fissureAngle[1], obliqueProportion[1],
                                                     rightLung, cache, coordinates,
                                                     nodes, nodetemplate, rightLungNodesetGroup,
                                                     lElementsCount1, lElementsCount2, lElementsCount3,
                                                     uElementsCount1, uElementsCount2, uElementsCount3,
                                                     lowerRightNodeIds, upperRightNodeIds, nodeIndex, nodeIdentifier)

            # Diaphragm lung nodes
            nodeIndex, nodeIdentifier = getDiaphragmaticLungNodes(spaceBelowCentre, cache, coordinates, nodes, nodetemplate,
                 diaphragmaticElementsCount1, diaphragmaticElementsCount2, diaphragmaticElementsCount3, width_1, length_1, height_1,
                 diaphragmaticNodeIds, nodeIndex, nodeIdentifier)

            # Create elements
            elementIdentifier = 1

            # Left lung elements
            elementIdentifier, leftUpperLobeElementID, leftLowerLobeElementID = \
                getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                                elementtemplateCustom, mesh, lungMeshGroup,
                                leftLungMeshGroup, None, None, None, mediastinumLeftGroupMeshGroup,
                                lElementsCount1, lElementsCount2, lElementsCount3,
                                uElementsCount1, uElementsCount2, uElementsCount3,
                                lowerLeftNodeIds, upperLeftNodeIds, elementIdentifier)

            # Right lung elements
            elementIdentifier, rightUpperLobeElementID, rightLowerLobeElementID = \
                getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                                elementtemplateCustom, mesh, lungMeshGroup,
                                rightLungMeshGroup, lowerRightLungMeshGroup, middleRightLungMeshGroup,
                                upperRightLungMeshGroup, mediastinumRightGroupMeshGroup,
                                lElementsCount1, lElementsCount2, lElementsCount3,
                                uElementsCount1, uElementsCount2, uElementsCount3,
                                lowerRightNodeIds, upperRightNodeIds, elementIdentifier)

            # Diaphragm lung elements
            getDiaphragmaticLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                elementtemplateCustom, mesh, lungMeshGroup,
                rightLungMeshGroup, diaphragmaticLungMeshGroup,
                diaphragmaticElementsCount1, diaphragmaticElementsCount2, diaphragmaticElementsCount3,
                diaphragmaticNodeIds, elementIdentifier)

            # Marker points
            lungNodesetGroup = lungGroup.getNodesetGroup(nodes)
            markerList = []

            lowerLeftElementCount = (lElementsCount1 * (lElementsCount2-1) * lElementsCount3 + lElementsCount1)

            idx = lowerLeftElementCount + (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)
            markerList.append({ "group" : leftApexGroup, "elementId" : idx, "xi" : [0.0, 0.0, 1.0]})

            idx = 1
            markerList.append({"group": leftDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 // 2)
            markerList.append({"group": leftMedialGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 - 1) + 1
            markerList.append({"group": leftVentralGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            upperLeftElementCount = (uElementsCount1 * uElementsCount2 * (uElementsCount3-1))
            leftLungElementCount = lowerLeftElementCount + upperLeftElementCount
            lowerRightElementCount = lowerLeftElementCount

            idx = leftLungElementCount + 1
            markerList.append({"group": rightDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = leftLungElementCount + lowerRightElementCount + \
                  (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount1//2)
            markerList.append({"group": rightApexGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx = leftLungElementCount + lElementsCount1 + 1
            markerList.append({"group": rightMedialGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            idx = leftLungElementCount + lElementsCount1 * lElementsCount2
            markerList.append({"group": rightVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * lElementsCount3) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 0.0, 1.0]})

            rightLungElementCount = leftLungElementCount

            idx_temp = diaphragmaticElementsCount1 * diaphragmaticElementsCount2 * (diaphragmaticElementsCount3 - 1) + 1
            idx = rightLungElementCount + leftLungElementCount + idx_temp
            markerList.append({"group": accessoryApexGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx_temp = diaphragmaticElementsCount1 * diaphragmaticElementsCount2 * (diaphragmaticElementsCount3 - 1) - 1
            idx = rightLungElementCount + leftLungElementCount + idx_temp
            markerList.append({"group": accessoryVentralGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            idx = rightLungElementCount + leftLungElementCount + 1
            markerList.append({"group": accessoryDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

        for marker in markerList:
            annotationGroup = marker["group"]
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            cache.setNode(markerPoint)
            markerLocation.assignMeshLocation(cache, mesh.findElementByIdentifier(marker["elementId"]),
                                              marker["xi"])
            markerName.assignString(cache, annotationGroup.getName())
            annotationGroup.getNodesetGroup(nodes).addNode(markerPoint)
            lungNodesetGroup.addNode(markerPoint)
            nodeIdentifier += 1

        # Transformation to the left and right lungs
        for i in range(2):
            if sharpeningFactor[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                sharpeningRidge(sharpeningFactor[i], fm, coordinates, LungNodeset, spaceFromCentre_temp)

            if taperingFactor[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                taperingZAxis(taperingFactor[i], fm, coordinates, LungNodeset, spaceFromCentre_temp)

            if DiaphramaticCurveRadius[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                concavingDiaphragmaticSurface(DiaphramaticCurveRadius[i], fm, coordinates, LungNodeset,
                                              spaceFromCentre_temp, length[i], width[i], height[i])

            if medialCurveRadius[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                width_temp = width[i] if i == 0 else -width[i]
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                concavingMedialSurface(medialCurveRadius[i], fm, coordinates, LungNodeset,
                                              spaceFromCentre_temp, length[i], width_temp, height[i])

            if bulgeRadiusY[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                bulgeRadiusY_temp = bulgeRadiusY[i] if i == 0 else -bulgeRadiusY[i]
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                bendingAroundYAxis(bulgeRadiusY_temp, fm, coordinates, LungNodeset, spaceFromCentre_temp)

            if bulgeRadiusZ[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                bulgeRadiusZ_temp = bulgeRadiusZ[i] if i == 0 else -bulgeRadiusZ[i]
                spaceFromCentre_temp = spaceFromCentre if i == 0 else [-spaceFromCentre[j] for j in range(3)]
                width_temp = width[i] if i == 0 else -width[i]
                bendingAroundZAxis(bulgeRadiusZ_temp, fm, coordinates, LungNodeset, spaceFromCentre_temp,
                                   length[i], width_temp, height[i])

            if tiltApex_xAxis[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                tiltApex_xAxis_temp = tiltApex_xAxis[i] if i == 0 else -tiltApex_xAxis[i]
                tiltLungs(tiltApex_xAxis_temp, 0, 0, 0, fm, coordinates, LungNodeset)

            if tiltApex_yAxis[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                tiltLungs(0, tiltApex_yAxis[i], 0, 0, fm, coordinates, LungNodeset)

            if rotate_ZAxis[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                rotate_ZAxis_temp = rotate_ZAxis[i] if i == 0 else -rotate_ZAxis[i]
                rotateLungs(rotate_ZAxis_temp, rotate_ZAxis_temp, rotate_ZAxis_temp, rotate_ZAxis_temp, fm, coordinates, LungNodeset)

            if tiltDiap_xAxis[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                tiltDiap_xAxis_temp = tiltDiap_xAxis[i] if i == 0 else -tiltDiap_xAxis[i]
                tiltLungs(0, 0, 0, tiltDiap_xAxis_temp, fm, coordinates, LungNodeset)

            if tiltDiap_yAxis[i] != 0:
                LungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                tiltLungs(0, 0, tiltDiap_yAxis[i], 0, fm, coordinates, LungNodeset)

        if discontinuity:
            # create discontinuity in d3 on the core boundary
            nodeIdentifier = createDiscontinuity(coordinates, nodes, mesh, cache, nodetemplate, leftUpperLobeElementID,
                            rightUpperLobeElementID, nodeIdentifier)

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCount = options['Refine number of elements']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        parameterSetName = options['Base parameter set']
        isMouse = 'Mouse 1' in parameterSetName
        isHuman = 'Human 1' in parameterSetName
        isRat = 'Rat 1' in parameterSetName

        # create fissure groups
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        # 1D Annotation
        is_exterior = fm.createFieldIsExterior()
        is_xi1_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)
        is_xi1_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)
        is_xi3_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)

        mediastanumLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("anterior mediastinum of left lung"))
        is_anteriorBorderLeftGroup = mediastanumLeftGroup.getFieldElementGroup(mesh1d)

        is_mediastanumLeftGroup_exterior = fm.createFieldAnd(is_anteriorBorderLeftGroup, is_exterior)
        is_mediastanumLeftGroup_exterior_xi1_0 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior, is_xi1_0)
        is_mediastanumLeftGroup_exterior_xi1_01 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi1_0, is_xi1_1)
        is_mediastanumLeftGroup_exterior_xi3_0 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior, is_xi3_0)

        is_upperLeftGroup_Diaphragm = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi3_0, is_mediastanumLeftGroup_exterior_xi1_01)

        is_ridge = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi1_01, fm.createFieldNot(is_upperLeftGroup_Diaphragm))
        anteriorBorderLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_lung_term("anterior border of left lung"))

        anteriorBorderLeftGroup.getMeshGroup(mesh1d).addElementsConditional(is_ridge)

        mediastanumRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("anterior mediastinum of right lung"))
        is_anteriorBorderRightGroup = mediastanumRightGroup.getFieldElementGroup(mesh1d)

        is_anteriorBorderRightGroup_exterior = fm.createFieldAnd(is_anteriorBorderRightGroup, is_exterior)
        is_anteriorBorderRightGroup_exterior_xi1_0 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior, is_xi1_0)
        is_anteriorBorderRightGroup_exterior_xi1_01 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi1_0, is_xi1_1)
        is_anteriorBorderRightGroup_exterior_xi3_0 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior, is_xi3_0)

        is_upperRightGroup_Diaphragm = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi3_0, is_anteriorBorderRightGroup_exterior_xi1_01)

        is_ridge = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi1_01, fm.createFieldNot(is_upperRightGroup_Diaphragm))
        anteriorBorderRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_lung_term("anterior border of right lung"))

        anteriorBorderRightGroup.getMeshGroup(mesh1d).addElementsConditional(is_ridge)

        if isHuman:
            upperLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of left lung"))
            lowerLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of left lung"))

            is_upperLeftGroup = upperLeftGroup.getFieldElementGroup(mesh2d)
            is_lowerLeftGroup = lowerLeftGroup.getFieldElementGroup(mesh2d)

            is_obliqueLeftGroup = fm.createFieldAnd(is_upperLeftGroup, is_lowerLeftGroup)
            obliqueLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of left lung"))
            obliqueLeftGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueLeftGroup)

        if isHuman or isMouse or isRat:
            upperRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of right lung"))
            middleRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("middle lobe of right lung"))
            lowerRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of right lung"))

            is_upperRightGroup = upperRightGroup.getFieldElementGroup(mesh2d)
            is_middleRightGroup = middleRightGroup.getFieldElementGroup(mesh2d)
            is_lowerRightGroup = lowerRightGroup.getFieldElementGroup(mesh2d)

            is_obliqueRightGroup = fm.createFieldAnd(fm.createFieldOr(is_middleRightGroup, is_upperRightGroup),
                                                     is_lowerRightGroup)
            is_horizontalRightGroup = fm.createFieldAnd(is_upperRightGroup, is_middleRightGroup)

            obliqueRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of right lung"))
            obliqueRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueRightGroup)
            horizontalRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("horizontal fissure of right lung"))
            horizontalRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_horizontalRightGroup)


def getLungNodes(spaceFromCentre, lengthUnit, widthUnit, heightUnit, fissureAngle, obliqueProportion,
                 lungSide, cache, coordinates, nodes, nodetemplate, lungSideNodesetGroup,
                 lElementsCount1, lElementsCount2, lElementsCount3,
                 uElementsCount1, uElementsCount2, uElementsCount3,
                 lowerNodeIds, upperNodeIds, nodeIndex, nodeIdentifier):
    """
    :param lowerNodeIds: nodeIdentifier array in the lower lobe filled by this function
        including indexing by [lElementsCount3 + 1][lElementsCount2 + 1][lElementsCount1 + 1]
    :param upperNodeIds: nodeIdentifier array in the upper lobe filled by this function
        including indexing by [uElementsCount3 + 1][uElementsCount2 + 1][uElementsCount1 + 1]
    :return: nodeIndex, nodeIdentifier
    """

    # Initialise parameters
    leftLung = 0
    d1 = [1.0, 0.0, 0.0]
    d2 = [0.0, 1.0, 0.0]
    d3 = [0.0, 0.0, 1.0]

    centre = [0.0, 0.0, 0.0]
    xAxis = [1.0, 0.0, 0.0]
    yAxis = [0.0, 1.0, 0.0]
    zAxis = [0.0, 0.0, 1.0]

    majorXAxis = [xAxis[i] * widthUnit for i in range(3)]
    minorYAxis = [yAxis[i] * -lengthUnit for i in range(3)]
    minorZAxis = [zAxis[i] * heightUnit for i in range(3)]

    fissureRadian = fissureAngle/180 * math.pi
    tanFissureRadian = math.tan(fissureRadian)

    # Create 5 points at the centre of 2D Ellipse base
    p1 = [yAxis[i] * -lengthUnit for i in range(3)]
    p2 = [yAxis[i] * lengthUnit for i in range(3)]
    p3 = [[p1[i] - (p1[i] - p2[i]) * obliqueProportion for i in range(3)]] # point of oblique fissure on diaphragm

    # Simplified simultaneous (linear and ellipsoid) equations
    C = tanFissureRadian * p3[0][1] # (y-offset) - y = mx + c

    # Quadratic equation - ax^2 + bx + c = 0
    a = (heightUnit ** 2) + ((lengthUnit ** 2) * (tanFissureRadian ** 2))
    b = (lengthUnit ** 2) * 2 * tanFissureRadian * C
    c = (lengthUnit ** 2) * (C ** 2) - ((lengthUnit ** 2) * (heightUnit ** 2))
    det = (b ** 2) - (4 * a * c) # determinant b^2 - 4ac
    assert det > 0, "No real solution"
    y = (-b + math.sqrt(det)) / (2 * a)
    z = tanFissureRadian*y + C
    p4 = [0.0, -y, z] # point at the oblique-ellispe intersection

    # points on 2D Ellipse base (two on the curve and one centre)
    x = p3[0][1]
    y = math.sqrt(1 - (x/lengthUnit) ** 2) * widthUnit
    p3.append([y, p3[0][1], p3[0][2]])
    p3.insert(0, [-y, p3[0][1], p3[0][2]])

    # points around on 2D Ellipse base
    # -------------------------------------------------- Lower lobe ------------------------------------------------
    lower_row1 = []
    lower_row1_d2 = []
    lower_row2 = []
    lower_row2_d2 = []
    obl = []
    obl_d2 = []
    lower_edge = []
    lower_edge_d2 = []
    lowerObl = []
    lowerObl_d2 = []

    for n2 in range(lElementsCount2 + 1):
        if n2 == 0:
            # row 1
            dx = p3[2][1]
            vec1 = minorYAxis
            vec2 = majorXAxis
            planeCentre = centre
            elementsCountAround = 4
            nodesCountAround = elementsCountAround + 1
            x = lower_row1
            nd2 = lower_row1_d2
        elif n2 == 1:
            # oblique
            vec1_1 = ([(p4[i] - p3[1][i]) / obliqueProportion for i in range(3)])
            vec1 = [vec1_1[i]/2 for i in range(3)]
            vec2 = majorXAxis
            planeCentre = [p4[i] - vec1[i] for i in range(3)]
            dx = magnitude([p3[1][i] - planeCentre[i] for i in range(3)])
            elementsCountAround = 4
            nodesCountAround = elementsCountAround + 1
            x = obl
            nd2 = obl_d2
        elif n2 == 2:
            # lower lobe
            dx = obl[0][1]
            elementsCountAround = 3
            vec1 = minorYAxis
            vec2 = [minorZAxis[i] * -1 for i in range(3)]
            planeCentre = centre
            nodesCountAround = elementsCountAround + 1
            x = lower_edge
            nd2 = lower_edge_d2
        elif n2 == 3:
            # 2nd oblique
            obliqueProportion_1 = (lengthUnit + lower_row1[3][1])/(lengthUnit*2)
            row1_temp = [0.0, lower_row1[3][1], lower_row1[3][2]]
            vec1_1 = ([(lower_edge[2][i] - row1_temp[i]) / obliqueProportion_1 for i in range(3)])
            vec1 = ([vec1_1[i]/2 for i in range(3)])
            vec2 = majorXAxis
            planeCentre = [lower_edge[2][i] - vec1[i] for i in range(3)]
            dx = magnitude([planeCentre[i] - row1_temp[i] for i in range(3)])
            elementsCountAround = 4
            nodesCountAround = elementsCountAround
            x = lowerObl
            nd2 = lowerObl_d2
        elif n2 == 4:
            # row 2
            row1_temp = [0.0, 0.0, lowerObl[3][2]]
            y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (lowerObl[3][2] / heightUnit) ** 2) * widthUnit
            row1_ellipse = [-y, 0.0, lowerObl[3][2]]
            vec1 = ([lower_edge[1][i] - row1_temp[i] for i in range(3)])
            vec2 = [-row1_ellipse[i] + row1_temp[i] for i in range(3)]
            planeCentre = row1_temp
            dx = lowerObl[3][1]
            elementsCountAround = 3
            nodesCountAround = elementsCountAround + 1
            x = lower_row2
            nd2 = lower_row2_d2

        magMajorAxis = magnitude(vec1)
        magMinorAxis = magnitude(vec2)
        unitMajorAxis = normalise(vec1)
        unitMinorAxis = normalise(vec2)

        radians = 0.0
        totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, -dx, 0.5 * math.pi)
        arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
        elementArcLength = arclengthAround / elementsCountAround
        radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

        x_temp = []
        nd2_temp = []
        for n1 in range(nodesCountAround):
            cosRadiansAround = math.cos(radians)
            sinRadiansAround = math.sin(radians)

            temp = ([(planeCentre[c] + cosRadiansAround * vec1[c] + sinRadiansAround * vec2[c]) for c in range(3)])
            x_temp.append(temp)

            ndab = setMagnitude([-sinRadiansAround * magMajorAxis, cosRadiansAround * magMinorAxis], -elementArcLength)
            nd2_temp.append(
                [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])
            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, -elementArcLength)

        # Gradual changes in node spacing
        if n2 == 0:
            # row 1
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 1, arcLengthDerivatives = True)
        elif n2 == 1:
            # oblique
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 1, arcLengthDerivatives = True)
        elif n2 == 2:
            # lower lobe
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2, lengthFractionStart = 0.3,
                                                             lengthFractionEnd = 1, arcLengthDerivatives = True)
        elif n2 == 3:
            # 2nd oblique
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2,  arcLengthDerivatives = True)
        elif n2 == 4:
            # rows 2
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp)-1, elementLengthStartEndRatio = 2, arcLengthDerivatives = True)

        [x.append(px[i]) for i in range(len(px))]
        [nd2.append(pd1[i]) for i in range(len(pd1))]

    # complete lower_row2 and lowerObl list
    lower_row2.append(obl[3])
    lower_row2_d2.append(obl_d2[3])
    lowerObl.append(lower_row1[3])
    lowerObl_d2.append(lower_row1[3])

    # smooth derivatives
    tx_d2 = []
    td2 = []
    tx_d3 = []
    td3 = []
    md2 = []
    md3 = []

    for n3 in range(lElementsCount3 + 1):
        tx_d2.append([])
        td2.append([])
        if n3 == 0:
            x = lower_row1
            d2 = lower_row1
        elif n3 == 1:
            x = lower_row2
            d2 = lower_row2
        elif n3 == 2:
            x = lowerObl
            d2 = lowerObl
        elif n3 == 3:
            x = obl
            d2 = obl
        for n2 in range(lElementsCount2 + 1):
            tx_d2[n3].append(x[n2])
            td2[n3].append(d2[n2])
            if n3 == 0:
                if n2 == 4:
                    reverse_obl = [obl[-i] for i in range(1, len(obl)+1)]
                    tx_d3.append(reverse_obl)
                    td3.append(reverse_obl)
                else:
                    tx_d3.append([lower_row1[n2], lower_row2[n2], lowerObl[n2], obl[n2]])
                    td3.append([lower_row1[n2], lower_row2[n2], lowerObl[n2], obl[n2]])
                md3.append(smoothCubicHermiteDerivativesLine(tx_d3[n2], td3[n2]))
        md2.append(smoothCubicHermiteDerivativesLine(tx_d2[n3], td2[n3]))

    # create nodes
    for n3 in range(lElementsCount3 + 1):
        lowerNodeIds.append([])
        for n2 in range(lElementsCount2 + 1):
            lowerNodeIds[n3].append([])
            for n1 in range(lElementsCount1 + 1):
                lowerNodeIds[n3][n2].append(None)
                d2 = md2[n3][n2]
                d3 = md3[n2][n3]

                # each i row
                if n3 == 0:
                    x = copy.deepcopy(lower_row1[n2])
                    if n2 < 4:
                        next_xd2 = copy.deepcopy(lower_row1[n2 + 1])
                elif n3 == 1:
                    x = copy.deepcopy(lower_row2[n2])
                    if n2 < 4:
                        next_xd2 = copy.deepcopy(lower_row2[n2 + 1])
                elif (n3 == 2) and (n2 < 3):
                    x = copy.deepcopy(lowerObl[n2])
                    if n2 < 4:
                        next_xd2 = copy.deepcopy(lowerObl[n2 + 1])
                elif (n3 == 3) and (n2 < 3):
                    x = copy.deepcopy(obl[n2])
                    if n2 < 2:
                        next_xd2 = copy.deepcopy(obl[n2 + 1])
                    else:
                        # 3-way point
                        d2 = [md3[n2][n3][0], md3[n2][n3][1], md3[n2][n3][2]]
                        d3 = [-md2[n3][n2][0], -md2[n3][n2][1], -md2[n3][n2][2]]
                else:
                    continue

                # skipping the first and last column and row - edge.
                if (n2 == 0) and (n1 != 1):
                    continue

                # symmetry
                if n1 == 0:
                    next_xd1 = [0.0, x[1], x[2]]
                    d1 = [next_xd1[i] - x[i] for i in range(3)]
                elif n1 == 1:
                    x = [0.0, x[1], x[2]]
                    next_x = [0.0, next_xd2[1], next_xd2[2]]
                    if n2 == 0:
                        d1 = [widthUnit * math.sin(0.5 * math.pi), 0.0, 0.0]
                        d2 = [next_x[i] - x[i] for i in range(3)]
                    # 3-way point
                    if (n3 == 3) and (n2 == 2):
                        d3 = [0.0, -md2[n3][n2][1], -md2[n3][n2][2]]
                    elif (n2 < 4):
                        next_xd2[0] = 0.0
                        d2 = [next_xd2[i] - x[i] for i in range(3)]
                        previous_d2 = d2
                    else:
                        d2 = previous_d2
                else:
                    x = [-x[0], x[1], x[2]]
                    d2 = [-d2[0], d2[1], d2[2]]
                    d3 = [-d3[0], d3[1], d3[2]]

                # translate right lung to the defined centre of the right lung
                if lungSide != leftLung:
                    x = [x[i] + spaceFromCentre[i] for i in range(3)]
                else:
                    x = [x[i] - spaceFromCentre[i] for i in range(3)]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                lowerNodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

                # Annotation
                lungSideNodesetGroup.addNode(node)
    # ---------------------------------------------------- Upper lobe --------------------------------------------
    upper_row1 = []
    upper_row2 = []
    upper_row3 = []
    upper_row4 = []
    upper_edge = []

    for n2 in range(uElementsCount1):
        if n2 == 0:
            # Upper lobe along the edge
            dx = 0.0
            vec1 = [minorYAxis[i] * -1 for i in range(3)]
            vec2 = minorZAxis
            planeCentre = centre
            elementsCountAround = 5
            nodesCountAround = elementsCountAround + 1
            x = upper_edge
        else:
            # Upper lobe along the edge on the lower lobe
            dx = obl[0][2]
            vec1 = minorZAxis
            vec2 = minorYAxis
            planeCentre = centre
            elementsCountAround = 3
            nodesCountAround = elementsCountAround + 1
            x = upper_edge

        radians = 0.0
        totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, dx, 0.5 * math.pi)
        arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
        elementArcLength = arclengthAround / elementsCountAround
        radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

        for n1 in range(nodesCountAround):
            cosRadians = math.cos(radians)
            sinRadians = math.sin(radians)
            x_temp = ([(planeCentre[c] + cosRadians * vec1[c] + sinRadians * vec2[c]) for c in range(3)])
            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, elementArcLength)
            if (n2 == 1) and (n1 == 0):
                continue
            x.insert(0, x_temp)

    # Interpolate the nodes between the upper lobe edge and fissures
    for n2 in range(uElementsCount3):
        if n2 == 0:
            # row 1
            dx = p3[2][1]
            vec1 = [minorYAxis[i] * -1 for i in range(3)]
            vec2 = majorXAxis
            planeCentre = centre
            elementsCountAround = 2
            nodesCountAround = elementsCountAround + 1
            x = upper_row1
        elif n2 == 1:
            # row 2
            obl_temp = [0.0, 0.0, obl[-2][2]]
            y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (obl[-2][2] / heightUnit) ** 2) * widthUnit
            obl_ellipse = [-y, 0.0, obl[-2][2]]
            vec1 = [upper_edge[-2][i] - obl_temp[i] for i in range(3)]
            vec2 = [-obl_ellipse[i] + obl_temp[i] for i in range(3)]
            planeCentre = obl_temp
            dx = obl[-2][1]
            elementsCountAround = 2
            nodesCountAround = elementsCountAround + 1
            x = upper_row2
        elif n2 == 2:
            # row 3
            obl_temp = [0.0, 0.0, obl[-3][2]]
            y = math.sqrt(1 - (0.0 / lengthUnit) ** 2 - (obl[-3][2] / heightUnit) ** 2) * widthUnit
            obl_ellipse = [-y, 0.0, obl[-3][2]]
            vec1 = [upper_edge[-3][i] - obl_temp[i] for i in range(3)]
            vec2 = [-obl_ellipse[i] + obl_temp[i] for i in range(3)]
            planeCentre = obl_temp
            dx = obl[-3][1]
            elementsCountAround = 2
            nodesCountAround = elementsCountAround + 1
            x = upper_row3

        radians = 0.0
        totalRadiansAround = getEllipseRadiansToX(magnitude(vec1), 0.0, dx, 0.5 * math.pi)
        arclengthAround = getEllipseArcLength(magnitude(vec1), magnitude(vec2), radians, totalRadiansAround)
        elementArcLength = arclengthAround / elementsCountAround
        radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), totalRadiansAround, -arclengthAround)

        x_temp = []
        nd2_temp = []
        for n1 in range(nodesCountAround):
            cosRadians = math.cos(radians)
            sinRadians = math.sin(radians)

            temp = ([(planeCentre[c] + cosRadians * vec1[c] + sinRadians * vec2[c]) for c in range(3)])
            x_temp.insert(0, temp)

            ndab = setMagnitude([-sinRadiansAround * magMajorAxis, cosRadiansAround * magMinorAxis], -elementArcLength)
            nd2_temp.append(
                [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])

            radians = updateEllipseAngleByArcLength(magnitude(vec1), magnitude(vec2), radians, -elementArcLength)

        # Gradual changes in node spacing
        if n2 == 0:
            # row 1
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                             elementLengthStartEndRatio=1.0, arcLengthDerivatives=True)
        elif n2 == 1:
            # row 2
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                             elementLengthStartEndRatio=1.0, arcLengthDerivatives=True)
        elif n2 == 2:
            # row 3
            px, pd1, pe, pxi, psf = sampleCubicHermiteCurves(x_temp, nd2_temp, len(x_temp) - 1,
                                                             elementLengthStartEndRatio=1.0, arcLengthDerivatives=True)
        elif n2 == 3:
            # row 4 - curve
            x = upper_row4
            px = []
            for j in range(5):
                x_temp_start = upper_edge[j + 1]
                if (j > 0) and (j < 4):
                    x_temp_end = obl[j] if j < 3 else upper_row3[1]
                    px_1, pd1, pe, pxi, psf = sampleCubicHermiteCurves([x_temp_start, x_temp_end], [x_temp_start, x_temp_end], 2,
                                                                     elementLengthStartEndRatio=1, arcLengthDerivatives=True)
                    x_ellipse = math.sqrt(1 - (px_1[1][1] / lengthUnit) ** 2 - (px_1[1][2] / heightUnit) ** 2) * widthUnit
                    x_temp = [-x_ellipse, px_1[1][1], px_1[1][2]]
                    px.append(x_temp)
                else:
                    px.append(x_temp_start)

        [x.append(px[i]) for i in range(len(px))]
        [nd2.append(pd1[i]) for i in range(len(pd1))]

    # smooth derivatives
    tx_d2 = []
    tx_d3 = []
    md2 = []
    md3 = []
    for n3 in range(uElementsCount3 + 1):
        if n3 == 0:
            tx_d2 = upper_row1
            tx_d3 = [upper_edge[-i] for i in range(1, 6)]
        elif n3 == 1:
            tx_d2 = upper_row2
            tx_d3 = [upper_row1[-n3-1], upper_row2[-n3-1], upper_row3[-n3-1], upper_row4[-n3-1], upper_edge[-5]]
        elif n3 == 2:
            tx_d2 = upper_row3
            tx_d3 = [obl[n3], upper_row4[n3], upper_edge[n3 + 1]]
        elif n3 == 3:
            tx_d2 = upper_row4
            tx_d3 = [obl[n3-2], upper_row4[n3-2], upper_edge[n3 - 1]]
        elif n3 == 4:
            # Apex row
            tx_d2 = [upper_edge[i] for i in range(1, 6)]
            tx_d3 = [upper_edge[i] for i in range(3)]

        md2.append(smoothCubicHermiteDerivativesLine(tx_d2, tx_d2))
        md3.append(smoothCubicHermiteDerivativesLine(tx_d3, tx_d3))

    # create nodes
    for i in range(uElementsCount3 + 1):
        upperNodeIds.append([])
        for j in range(uElementsCount2 + 1):
            upperNodeIds[i].append([])
            for k in range(uElementsCount1 + 1):
                upperNodeIds[i][j].append(None)

                # Oblique fissure nodes
                if (i < 2) and (j == 2):
                    upperNodeIds[i][j][k] = lowerNodeIds[i][-1][k]
                elif (i == 2) and (j < 4):
                    upperNodeIds[i][j][k] = lowerNodeIds[-1][j][k]

                # each i row
                if (i == 0) and (j > 2):
                    x = copy.deepcopy(upper_row1[j - 2])
                    d2 = md2[i][j-2]
                    idx = j-(j//2)*2
                    d3 = md3[idx][i]
                    if j < 4:
                        next_xd2 = copy.deepcopy(upper_row1[j - 1])
                elif (i == 1) and (j > 2):
                    x = copy.deepcopy(upper_row2[j - 2])
                    d2 = md2[i][j-2]
                    idx = j-(j//2)*2
                    d3 = md3[idx][i]
                    if j < 4:
                        next_xd2 = copy.deepcopy(upper_row2[j - 1])
                elif (i == 2) and (j > 2):
                    x = copy.deepcopy(upper_row3[j - 2])
                    d2 = md2[i][j-2]
                    idx = j-(j//2)*2
                    d3 = md3[idx][i]
                    if j < 4:
                        next_xd2 = copy.deepcopy(upper_row3[j - 1])
                elif (i == 3):
                    x = copy.deepcopy(upper_row4[j])
                    d2 = md2[i][j]
                    d3 = md3[-j-1][i] if j > 2 else md3[-j-1][i-2]
                    if j < 4:
                        next_xd2 = copy.deepcopy(upper_row4[j+1])
                elif (i == 4) and (j > 0) and (j < 4):
                    x = copy.deepcopy(upper_edge[j+1])
                    d2 = md2[i][j]
                    d3 = md3[-j-1][i] if j == 3 else md3[-j-1][i-3]
                    if j < 4:
                        next_xd2 = copy.deepcopy(upper_edge[j+2])
                else:
                    continue

                # skipping the first and last, k.
                if (((j == 0) or (j == 4)) and (k != 1)) or ((i == 4) and (k != 1)):
                    continue

                # symmetry
                if k == 0:
                    # d2 = md2[i][j]
                    pass
                elif k == 1:
                    if j == 4:
                        # Ridges
                        d2 = previous_d2
                    elif i == 4:
                        # Apex row
                        d3 = [0.0, d3[1], d3[2]]
                    else:
                        x = [0.0, x[1], x[2]]
                        next_xd2 = [0.0, next_xd2[1], next_xd2[2]]
                        d2 = [next_xd2[i] - x[i] for i in range(3)]
                        previous_d2 = d2
                        d3 = [0.0, d3[1], d3[2]]
                else:
                    x = [-x[0], x[1], x[2]]
                    d2 = [-d2[0], d2[1], d2[2]]
                    d3 = [-d3[0], d3[1], d3[2]]

                # translate right lung to the defined centre of the right lung
                if lungSide != leftLung:
                    x = [x[i] + spaceFromCentre[i] for i in range(3)]
                else:
                    x = [x[i] - spaceFromCentre[i] for i in range(3)]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                upperNodeIds[i][j][k] = nodeIdentifier
                nodeIdentifier += 1

                # Annotation
                lungSideNodesetGroup.addNode(node)

    return nodeIndex, nodeIdentifier

def getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular, elementtemplateCustom, mesh,
                    lungMeshGroup, lungSideMeshGroup, lowerLobeMeshGroup, middleLobeMeshGroup,
                    upperLobeMeshGroup, mediastinumGroupMeshGroup,
                    lElementsCount1, lElementsCount2, lElementsCount3,
                    uElementsCount1, uElementsCount2, uElementsCount3,
                    lowerNodeIds, upperNodeIds, elementIdentifier):
    """
    :param lowerNodeIds: Indexing by [lElementsCount3 + 1][lElementsCount2 + 1][lElementsCount1 + 1]
    :param upperNodeIds: Indexing by [uElementsCount3 + 1][uElementsCount2 + 1][uElementsCount1 + 1]
    :return: elementIdentifier
    """

    eftWedgeCollapseXi1_15 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
    eftWedgeCollapseXi1_26 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
    eftWedgeCollapseXi1_57 = eftfactory.createEftWedgeCollapseXi1Quadrant([5, 7])
    eftWedgeCollapseXi1_68 = eftfactory.createEftWedgeCollapseXi1Quadrant([6, 8])
    eftWedgeCollapseXi2_78 = eftfactory.createEftWedgeCollapseXi2Quadrant([7, 8])
    eftTetCollapseXi1Xi2_71 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(7, 1)
    eftTetCollapseXi1Xi2_82 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(8, 2)

    lowerLobeElementID = []
    upperLobeElementID = []
    # Lower lobe elements
    for e3 in range(lElementsCount3):
        lowerLobeElementID.append([])
        for e2 in range(lElementsCount2):
            lowerLobeElementID[e3].append([])
            for e1 in range(lElementsCount1):
                lowerLobeElementID[e3][e2].append(None)

                eft = eftRegular
                nodeIdentifiers = [
                    lowerNodeIds[e3][e2][e1], lowerNodeIds[e3][e2][e1 + 1], lowerNodeIds[e3][e2 + 1][e1],
                    lowerNodeIds[e3][e2 + 1][e1 + 1],
                    lowerNodeIds[e3 + 1][e2][e1], lowerNodeIds[e3 + 1][e2][e1 + 1],
                    lowerNodeIds[e3 + 1][e2 + 1][e1], lowerNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e2 == 0) and (e1 == 0):
                    # Back wedge elements
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftWedgeCollapseXi1_15
                elif (e2 == 0) and (e1 == (lElementsCount1 - 1)):
                    # Back wedge elements
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(1)
                    eft = eftWedgeCollapseXi1_26
                elif (e3 == 1) and (e2 == (lElementsCount2 - 2)):
                    # Middle wedge
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    eft = eftWedgeCollapseXi2_78
                elif (e3 == (lElementsCount3 - 1)) and (e2 == (lElementsCount2 - 3)):
                    # Remapped cube element 1
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                elif (e3 == (lElementsCount3 - 1)) and (e2 == (lElementsCount2 - 2)):
                    # Remapped cube element 2
                    nodeIdentifiers[2] = lowerNodeIds[e3 - 1][e2 + 1][e1]
                    nodeIdentifiers[3] = lowerNodeIds[e3 - 1][e2 + 1][e1 + 1]
                    nodeIdentifiers[6] = lowerNodeIds[e3 - 1][e2 + 2][e1]
                    nodeIdentifiers[7] = lowerNodeIds[e3 - 1][e2 + 2][e1 + 1]
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, [])])
                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                lowerLobeElementID[e3][e2][e1] = elementIdentifier
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if lowerLobeMeshGroup:
                    lowerLobeMeshGroup.addElement(element)

    # Upper lobe elements
    for e3 in range(uElementsCount3):
        upperLobeElementID.append([])
        for e2 in range(uElementsCount2):
            upperLobeElementID[e3].append([])
            for e1 in range(uElementsCount1):
                upperLobeElementID[e3][e2].append(None)
                is_mediastanum = False

                eft = eftRegular
                nodeIdentifiers = [
                    upperNodeIds[e3][e2][e1], upperNodeIds[e3][e2][e1 + 1], upperNodeIds[e3][e2 + 1][e1],
                    upperNodeIds[e3][e2 + 1][e1 + 1],
                    upperNodeIds[e3 + 1][e2][e1], upperNodeIds[e3 + 1][e2][e1 + 1],
                    upperNodeIds[e3 + 1][e2 + 1][e1], upperNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                    is_mediastanum = True
                    # Distal-front wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(2)
                    eft = eftfactory.createEftBasic()
                    nodes = [3, 4, 7, 8]
                    collapseNodes = [3, 7]
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e3 < (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    is_mediastanum = True
                    # Distal-back wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(3)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [3, 4, 7, 8]
                    collapseNodes = [4, 8]
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    ln_map = [1, 2, 3, 3, 4, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e3 == (uElementsCount3 - 2)) and (e2 == 0) and (e1 == 0):
                    # Medial-front wedge elements
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftWedgeCollapseXi1_15
                elif (e3 == (uElementsCount3 - 2)) and (e2 == 0) and (e1 == (uElementsCount1 - 1)):
                    # Medial-back wedge elements
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(1)
                    eft = eftWedgeCollapseXi1_26
                elif (e3 == (uElementsCount3 - 1)) and (0 < e2 < (uElementsCount2 - 1)) and (e1 == 0):
                    # Top-front wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    eft = eftWedgeCollapseXi1_57
                elif (e3 == (uElementsCount3 - 1)) and (0 < e2 < (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    # Top-back wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    eft = eftWedgeCollapseXi1_68
                elif (e3 == (uElementsCount3 - 1)) and (e2 == 0) and (e1 == 0):
                    # Top-front-medial tetrahedron wedge elements
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(0)
                    eft = eftTetCollapseXi1Xi2_82
                elif (e3 == (uElementsCount3 - 1)) and (e2 == 0) and (e1 == (uElementsCount1 - 1)):
                    # Top-back-medial tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(1)
                    eft = eftTetCollapseXi1Xi2_71
                elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (e1 == 0):
                    is_mediastanum = True
                    # Top-front-distal tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    nodeIdentifiers.pop(2)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    # remap parameters on xi3 = 1 before collapsing nodes
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, [3], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 4, 4, 4]
                    remapEftLocalNodes(eft, 4, ln_map)

                elif (e3 == (uElementsCount3 - 1)) and (e2 == (uElementsCount2 - 1)) and (
                        e1 == (uElementsCount1 - 1)):
                    is_mediastanum = True
                    # Top-front-distal tetrahedron wedge elements
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(5)
                    nodeIdentifiers.pop(3)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    # remap parameters on xi3 = 1 before collapsing nodes
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS2, [])
                    remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, [4], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                    ln_map = [1, 2, 3, 3, 4, 4, 4, 4]
                    remapEftLocalNodes(eft, 4, ln_map)

                elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 3)):
                    # Remapped cube element 1
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                elif (e3 == (uElementsCount3 - 2)) and (e2 == (uElementsCount2 - 2)):
                    # Remapped cube element 2
                    eft = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft, [1, 2], Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                upperLobeElementID[e3][e2][e1] = elementIdentifier
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)
                if middleLobeMeshGroup and (e3 < (uElementsCount3 - 2)):
                    middleLobeMeshGroup.addElement(element)
                elif upperLobeMeshGroup:
                    upperLobeMeshGroup.addElement(element)
                if is_mediastanum is True:
                    mediastinumGroupMeshGroup.addElement(element)

    return elementIdentifier, upperLobeElementID, lowerLobeElementID

def getDiaphragmaticLungNodes(spaceBelowCentre, cache, coordinates, nodes,
                              nodetemplate, elementsCount1, elementsCount2, elementsCount3,
                              elementsUnit1, elementsUnit2, elementsUnit3,
                              nodeIds, nodeIndex, nodeIdentifier):
    """
    Create a 3D triangular mesh from getDiaphragmaticLungNodes
    :parameter: elementsCount1 - x, elementsCount2 - y, elementsCount3 - z
    :return: nodeIndex, nodeIdentifier
    """

    unitRatio = 0.5 # the ratio of the unit length in the second row

    # Initialise parameters
    d1 = [elementsUnit1, 0.0, 0.0]
    d2 = [0.0, elementsUnit2, 0.0]
    d3 = [0.0, 0.0, elementsUnit3]

    # Diaphragmatic lobe nodes
    for n3 in range(elementsCount3 + 1):
        nodeIds.append([])
        for n2 in range(elementsCount2 + 1):
            nodeIds[n3].append([])
            for n1 in range(elementsCount1 + 1):
                nodeIds[n3][n2].append(None)

                if ((n1 == elementsCount1) or (n1 != 1)) and (n3 == elementsCount3):
                    continue

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)

                if n3 == (elementsCount3 - 1):
                    # second row
                    x = [elementsUnit1 * unitRatio * (n1 - 1), elementsUnit2 * (n2 - 2), elementsUnit3 * n3]
                    next_row_x = [0, elementsUnit2 * (n2 - 2), elementsUnit3 * (n3 + 1)]
                    d1 = [elementsUnit1 * unitRatio, 0.0, 0.0]
                else:
                    x = [elementsUnit1 * (n1 - 1), elementsUnit2 * (n2 - 2), elementsUnit3 * n3]
                    next_row_x = [elementsUnit1 * unitRatio * (n1 - 1), elementsUnit2 * (n2 - 2), elementsUnit3 * (n3 + 1)]

                d3 = [next_row_x[i] - x[i] for i in range(3)]

                # Translate the mesh
                x = [x[i] - spaceBelowCentre[i] for i in range(3)]

                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                nodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

    return nodeIndex, nodeIdentifier

def getDiaphragmaticLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular, elementtemplateCustom,
                    mesh, lungMeshGroup, lungSideMeshGroup, diaphragmaticLobeMeshGroup,
                    elementsCount1, elementsCount2, elementsCount3,
                    NodeIds, elementIdentifier):
    """
    Create a 3D triangular mesh from getDiaphragmaticLungNodes
    :parameter: elementsCount1 - x, elementsCount2 - y, elementsCount3 - z
    :return: elementIdentifier
    """

    # Diaphragmatic lobe elements
    for e3 in range(elementsCount3):
        for e2 in range(elementsCount2):
            for e1 in range(elementsCount1):
                eft = eftRegular
                nodeIdentifiers = [
                    NodeIds[e3][e2][e1], NodeIds[e3][e2][e1 + 1], NodeIds[e3][e2 + 1][e1],
                    NodeIds[e3][e2 + 1][e1 + 1],
                    NodeIds[e3 + 1][e2][e1], NodeIds[e3 + 1][e2][e1 + 1],
                    NodeIds[e3 + 1][e2 + 1][e1], NodeIds[e3 + 1][e2 + 1][e1 + 1]]

                if (e1 == 0) and (e3 == (elementsCount3 - 1)):
                    # wedge elements along crest
                    nodeIdentifiers.pop(6)
                    nodeIdentifiers.pop(4)
                    eft = eftfactory.createEftBasic()
                    nodes = [5, 6, 7, 8]
                    collapseNodes = [5, 7]
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])

                    ln_map = [1, 2, 3, 4, 5, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif (e1 == elementsCount1 - 1) and (e3 == (elementsCount3 - 1)):
                    nodeIdentifiers.pop(7)
                    nodeIdentifiers.pop(5)
                    eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft, [1], [])
                    nodes = [5, 6, 7, 8]
                    collapseNodes = [6, 8]
                    remapEftNodeValueLabel(eft, collapseNodes, Node.VALUE_LABEL_D_DS3,
                                           [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                    remapEftNodeValueLabel(eft, nodes, Node.VALUE_LABEL_D_DS1, [])
                    ln_map = [1, 2, 3, 4, 5, 5, 6, 6]
                    remapEftLocalNodes(eft, 6, ln_map)

                elif None in nodeIdentifiers:
                    continue

                if eft is eftRegular:
                    element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                else:
                    elementtemplateCustom.defineField(coordinates, -1, eft)
                    element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if eft.getNumberOfLocalScaleFactors() == 1:
                    element.setScaleFactors(eft, [-1.0])
                elementIdentifier += 1

                # Annotation
                lungMeshGroup.addElement(element)
                diaphragmaticLobeMeshGroup.addElement(element)
                lungSideMeshGroup.addElement(element)

    return elementIdentifier

def createDiscontinuity(coordinates, nodes, mesh, fieldcache, nodetemplate, leftUpperLobeElementID, rightUpperLobeElementID, nodeIdentifier):
    """
    Refine source mesh into separate region, with change of basis.
    :param meshrefinement: MeshRefinement, which knows source and target region.
    :param options: Dict containing options. See getDefaultOptions().
    """
    # Replicating nodes at fissures
    for e4 in range(2):
        upperNodeIds = []
        store_coordinate = []
        UpperLobeElementID = leftUpperLobeElementID if e4 == 0 else rightUpperLobeElementID
        if UpperLobeElementID == None:
            continue

        for e3 in range(len(UpperLobeElementID)):
            upperNodeIds.append([])
            for e2 in range(len(UpperLobeElementID[e3])):
                for e1 in range(len(UpperLobeElementID[e3][e2])):
                    # Fissures along the upper lobe
                    if ((e2 == 2) and (e3 < 2)) or ((e2 < 4) and (e3 == 2)):
                        elementIdentifier = UpperLobeElementID[e3][e2][e1]
                        element = mesh.findElementByIdentifier(elementIdentifier)
                        eft = element.getElementfieldtemplate(coordinates, -1)

                        # Exclude horizontal fissure in the left lobe
                        if (e4 == 0) and (e3 == 2) and (e2 > 1):
                            continue

                        if e3 < 2:
                            # Middle lobe
                            localNodeIndexes = [1, 2] if e1 == 0 else [2]
                            if e3 == 1:
                                if e1 == 0:
                                    localNodeIndexes.append(5)
                                localNodeIndexes.append(6)
                        else:
                            # Upper lobe
                            if e2 == 0:
                                # Dorsal wedge
                                localNodeIndexes = [1] if (e1 == 0) else []
                            elif e2 == 3:
                                # Ventral wedge
                                localNodeIndexes = [1, 2] if (e1 == 0) else [2, 3]
                            else:
                                # Regular element
                                localNodeIndexes = [1, 2] if (e1 == 0) else [2]

                        for localNodeIndex in localNodeIndexes:
                            node = element.getNode(eft, localNodeIndex)
                            nodetemplate.defineFieldFromNode(coordinates, node)
                            versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1,
                                                                                  Node.VALUE_LABEL_VALUE)

                            if versionsCount == 1:
                                fieldcache.setNode(node)
                                result0, x = coordinates.getNodeParameters(fieldcache, -1,
                                                                           Node.VALUE_LABEL_VALUE,
                                                                           1, 3)
                                result0, d1 = coordinates.getNodeParameters(fieldcache, -1,
                                                                            Node.VALUE_LABEL_D_DS1,
                                                                            1, 3)
                                result0, d2 = coordinates.getNodeParameters(fieldcache, -1,
                                                                            Node.VALUE_LABEL_D_DS2,
                                                                            1, 3)
                                result0, d3 = coordinates.getNodeParameters(fieldcache, -1,
                                                                            Node.VALUE_LABEL_D_DS3,
                                                                            1, 3)

                                if (localNodeIndex > 2) and (e3 == 1):
                                    # Store 3-way points in the middle lobe
                                    store_coordinate.append([x, d1, d2, d3])
                                    continue

                                node = nodes.createNode(nodeIdentifier, nodetemplate)
                                fieldcache.setNode(node)
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1,
                                                              x)
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                                              d1)
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                              d2)
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                              d3)
                                upperNodeIds[e3].append(nodeIdentifier)
                                nodeIdentifier += 1

                        # Create 3-way points in the middle lobe
                        if (e3 == 1) and (e2 == 2) and (e1 == 1):
                            for i in range(len(store_coordinate)):
                                node = nodes.createNode(nodeIdentifier, nodetemplate)
                                fieldcache.setNode(node)
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1,
                                                              store_coordinate[i][0])
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                                              store_coordinate[i][1])
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                              store_coordinate[i][2])
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                              store_coordinate[i][3])
                                upperNodeIds[e3].append(nodeIdentifier)
                                nodeIdentifier += 1

        # Change node index in the element
        for e3 in range(len(UpperLobeElementID)):
            # middle lobe and upper lobe idx
            temp_idx = upperNodeIds[0] + upperNodeIds[1] if e3 < 2 else upperNodeIds[2]
            # print("temp_idx", temp_idx)
            for e2 in range(len(UpperLobeElementID[e3])):
                for e1 in range(len(UpperLobeElementID[e3][e2])):

                    # Fissures along the upper lobe
                    if ((e2 == 2) and (e3 < 2)) or ((e2 < 4) and (e3 == 2)):
                        elementIdentifier = UpperLobeElementID[e3][e2][e1]
                        element = mesh.findElementByIdentifier(elementIdentifier)
                        eft = element.getElementfieldtemplate(coordinates, -1)
                        nodeIds = get_element_node_identifiers(element, eft)

                        # Exclude horizontal fissure in the left lobe
                        if (e4 == 0) and (e3 == 2) and (e2 > 2):
                            continue

                        if (e3 == 2) and (e2 == 0):
                            # Dorsal wedge
                            nd2 = e1 + 1
                            nodeIds[0] = temp_idx[0]
                            nodeIds[1:3] = [temp_idx[nd2], temp_idx[nd2 + 1]]
                        elif (e3 == 2) and (e2 < 3):
                            # upper and lower lobes
                            if e4 == 1:
                                # Right lung
                                nd1 = (e2 * e2) + e1
                                nd2 = (e2 * 3) + e1 + 1
                                nodeIds[:2] = [temp_idx[nd1], temp_idx[nd1 + 1]]
                                nodeIds[2:4] = [temp_idx[nd2], temp_idx[nd2 + 1]]
                            else:
                                # LeftLung
                                nd1 = (e2 * e2) + e1
                                nd2 = (e2 * 2) + e1 + 1 if (e2 < 2) else e2 + e1 + 1
                                temp_idx_1 = upperNodeIds[1]
                                nodeIds[:2] = [temp_idx[nd1], temp_idx[nd1 + 1]] if (e2 < 2) else \
                                    [temp_idx_1[nd2], temp_idx_1[nd2 + 1]]
                                if (e2 < 2):
                                    nodeIds[2:4] = [temp_idx_1[nd2], temp_idx_1[nd2 + 1]]
                        elif (e3 == 2) and (e2 == 3):
                            # Ventral wedge
                            nd2 = (e3 * 3) + e1 + 1
                            nodeIds[2] = temp_idx[-1]
                            nodeIds[:2] = [temp_idx[nd2], temp_idx[nd2 + 1]]
                        else:
                            # Middle lobe
                            nd1 = (e3 * 3) + e1
                            nd2 = ((e3 + 1) * 3) + e1
                            nodeIds[:2] = [temp_idx[nd1], temp_idx[nd1 + 1]]
                            nodeIds[4:6] = [temp_idx[nd2], temp_idx[nd2 + 1]]
                        element.setNodesByIdentifier(eft, nodeIds)

    return nodeIdentifier

def tiltLungs(tiltApex_xAxis, tiltApex_yAxis, tiltDiap_yAxis, tiltDiap_xAxis, fm, coordinates, lungNodesetGroup):
    """
    :param tiltDegree: [tilted degree for apex, for diaphragm]
    :param fm:
    :param coordinates:
    :param nodes:
    :return: transformed lungs
    """
    # FieldConstant - Matrix = [   x1,    x4, sh_zx,
    #                              x2,    x5, sh_zy,
    #                           sh_xz, sh_yz,    x9]
    sh_xz = tiltApex_xAxis / 180 * math.pi
    sh_yz = tiltApex_yAxis / 180 * math.pi
    sh_zy = tiltDiap_yAxis / 180 * math.pi
    sh_zx = tiltDiap_xAxis / 180 * math.pi
    shearMatrix = fm.createFieldConstant([1.0, 0.0, sh_xz, 0.0, 1.0, sh_yz, sh_zx, sh_zy, 1.0])
    newCoordinates = fm.createFieldMatrixMultiply(3, shearMatrix, coordinates)
    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def rotateLungs(tiltApex_xAxis, tiltApex_yAxis, tiltDiap_yAxis, tiltDiap_xAxis, fm, coordinates, lungNodesetGroup):
    """
    :param tiltDegree: [tilted degree for apex, for diaphragm]
    :param fm:
    :param coordinates:
    :param nodes:
    :return: transformed lungs
    """
    # FieldConstant - Matrix = [   x1,    x4, sh_zx,
    #                              x2,    x5, sh_zy,
    #                           sh_xz, sh_yz,    x9]
    sh_xx = tiltApex_xAxis / 180 * math.pi
    sh_xy = tiltApex_yAxis / 180 * math.pi
    sh_yx = tiltDiap_yAxis / 180 * math.pi
    sh_yy = tiltDiap_xAxis / 180 * math.pi
    shearMatrix = fm.createFieldConstant([math.cos(sh_xx), math.sin(sh_xy), 0.0,
                                          -math.sin(sh_yx), math.cos(sh_yy), 0.0,
                                          0.0, 0.0, 1.0])
    newCoordinates = fm.createFieldMatrixMultiply(3, shearMatrix, coordinates)
    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def concavingDiaphragmaticSurface(bulgeRadius, fm, coordinates, lungNodesetGroup, spaceFromCentre, length, width, height):
    """
    Quadratic transformation
    :param sharpRadius:
    :param fm:
    :param coordinates:
    :param lungNodesetGroup:
    :param spaceFromCentre:
    :return:
    """
    # Transformation matrix = [       1,                   |  [x,
    #                                 1,                   |   y,
    #                         ,k1x^2 + k2y^2 + (0)k3z^2]   |   z]

    offset = fm.createFieldConstant([spaceFromCentre[0], spaceFromCentre[1], spaceFromCentre[2] - height])
    origin = fm.createFieldAdd(coordinates, offset)
    absCoordinate = fm.createFieldAbs(origin)
    sqauredCoordinate = fm.createFieldMultiply(absCoordinate, absCoordinate)

    k1 = length/width * bulgeRadius  # quadratic coefficient for x
    k2 = width/length * bulgeRadius  # quadratic coefficient for y
    k3 = 0.0 * bulgeRadius  # quadratic coefficient for z

    scale = fm.createFieldConstant([k1, 0.0, 0.0])
    scaleFunction = fm.createFieldMultiply(sqauredCoordinate, scale)

    scale_1 = fm.createFieldConstant([0.0, k2, 0.0])
    scaleFunction_1 = fm.createFieldMultiply(sqauredCoordinate, scale_1)
    transformation_matrix_1 = fm.createFieldComponent(scaleFunction_1, [2, 1, 1])

    scale_2 = fm.createFieldConstant([0.0, 0.0, k3])
    scaleFunction_2 = fm.createFieldMultiply(absCoordinate, scale_2)
    transformation_matrix_2 = fm.createFieldComponent(scaleFunction_2, [3, 1, 1])

    twoVariablesFunction = fm.createFieldAdd(scaleFunction, transformation_matrix_1)
    threeVariablesFunction = fm.createFieldAdd(transformation_matrix_2, twoVariablesFunction)

    constant = fm.createFieldConstant([1.0, 1.0, 0.0])
    constantFunction = fm.createFieldAdd(threeVariablesFunction, constant)

    transformation_matrix = fm.createFieldComponent(constantFunction, [2, 2, 1])
    curve_coordinates = fm.createFieldMultiply(transformation_matrix, origin)

    ratio = (width ** 2) * bulgeRadius * (length/width) + 1

    reduceSize = fm.createFieldConstant([1.0, 1.0, 1.0/ratio])
    reduce_coordinates = fm.createFieldMultiply(curve_coordinates, reduceSize)

    translate_coordinates = fm.createFieldSubtract(reduce_coordinates, offset)
    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def concavingMedialSurface(bulgeRadius, fm, coordinates, lungNodesetGroup, spaceFromCentre, length, width, height):
    """
    Quadratic transformation
    :param sharpRadius:
    :param fm:
    :param coordinates:
    :param lungNodesetGroup:
    :param spaceFromCentre:
    :return:
    """
    # Transformation matrix = [       1,                   |  [x,
    #                                 1,                   |   y,
    #                         ,k1x^2 + k2y^2 + (0)k3z^2]   |   z]
    # the centre of the curve depends on the offset
    offset = fm.createFieldConstant([spaceFromCentre[0] + width, spaceFromCentre[1], spaceFromCentre[2] - height/2])
    origin = fm.createFieldAdd(coordinates, offset)
    absCoordinate = fm.createFieldAbs(origin)
    sqauredCoordinate = fm.createFieldMultiply(absCoordinate, absCoordinate)

    k1 = -0.7 * bulgeRadius  # quadratic coefficient for x
    k2 = 0.7 * bulgeRadius  # quadratic coefficient for y
    k3 = 0.7 * bulgeRadius # quadratic coefficient for z

    scale = fm.createFieldConstant([k1, 0.0, 0.0])
    scaleFunction = fm.createFieldMultiply(sqauredCoordinate, scale)

    scale_1 = fm.createFieldConstant([0.0, k2, 0.0])
    scaleFunction_1 = fm.createFieldMultiply(sqauredCoordinate, scale_1)
    transformation_matrix_1 = fm.createFieldComponent(scaleFunction_1, [2, 1, 1])

    scale_2 = fm.createFieldConstant([0.0, 0.0, k3])
    scaleFunction_2 = fm.createFieldMultiply(sqauredCoordinate, scale_2)
    transformation_matrix_2 = fm.createFieldComponent(scaleFunction_2, [3, 1, 1])

    twoVariablesFunction = fm.createFieldAdd(scaleFunction, transformation_matrix_1)
    threeVariablesFunction = fm.createFieldAdd(transformation_matrix_2, twoVariablesFunction)

    constant = fm.createFieldConstant([1.0, 1.0, 0.0])
    constantFunction = fm.createFieldAdd(threeVariablesFunction, constant)

    transformation_matrix = fm.createFieldComponent(constantFunction, [1, 2, 2])
    curve_coordinates = fm.createFieldMultiply(transformation_matrix, origin)
    translate_coordinates = fm.createFieldSubtract(curve_coordinates, offset)
    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def bendingAroundYAxis(bulgeRadius, fm, coordinates, lungNodesetGroup, spaceFromCentre):
    """
    :param bulgeRadius: the radius and the centre of curvature to transfrom the scaffold
    :param fm:
    :param coordinates:
    :param nodes:
    :return:
    """
    # cylindrical polar coordinates (x = r*cos(theta), y = r*sin(theta), z = z):
    # r = x - bulgeRadius
    # theta = -z / bulgeRadius
    # z = y
    xzyCoordinates = fm.createFieldComponent(coordinates, [1, 3, 2])
    scale = fm.createFieldConstant([1.0, -1.0 / bulgeRadius, 1.0])
    scaleCoordinates = fm.createFieldMultiply(xzyCoordinates, scale)
    offset_temp = [-bulgeRadius + spaceFromCentre[0], spaceFromCentre[1], spaceFromCentre[2]]
    offset = fm.createFieldConstant(offset_temp)
    polarCoordinates = fm.createFieldAdd(scaleCoordinates, offset)
    polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
    rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
    rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    newxzyCoordinates = fm.createFieldSubtract(rcCoordinates, offset)
    newCoordinates = fm.createFieldComponent(newxzyCoordinates, [1, 3, 2])
    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def bendingAroundZAxis(bulgeRadius, fm, coordinates, lungNodesetGroup, spaceFromCentre, length, width, height):
    """
    :param bulgeRadius: the radius and the centre of curvature to transfrom the scaffold
    :param fm:
    :param coordinates:
    :param nodes:
    :return:
    """
    # cylindrical polar coordinates (x = r*cos(theta), y = r*sin(theta), z = z):
    # r = x - bulgeRadius
    # theta = -y / bulgeRadius
    # z = z
    scale = fm.createFieldConstant([1.0, -1.0 / bulgeRadius, 1.0])
    scaleCoordinates = fm.createFieldMultiply(coordinates, scale)
    offset = fm.createFieldConstant([-bulgeRadius + spaceFromCentre[0], spaceFromCentre[1], spaceFromCentre[2]])
    polarCoordinates = fm.createFieldAdd(scaleCoordinates, offset)
    polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
    rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
    rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    newxyzCoordinates = fm.createFieldSubtract(rcCoordinates, offset)
    fieldassignment = coordinates.createFieldassignment(newxyzCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def sharpeningRidge(sharpeningFactor, fm, coordinates, lungNodesetGroup, spaceFromCentre):
    """
    Linear transformation
    :param sharpRadius:
    :param fm:
    :param coordinates:
    :param lungNodesetGroup:
    :param spaceFromCentre:
    :return:
    """
    # Transformation matrix = [     1, | [x,
    #                         k1y + 1, |  y,
    #                               1] |  z]
    offset = fm.createFieldConstant(spaceFromCentre)
    origin = fm.createFieldAdd(coordinates, offset)
    k1 = -sharpeningFactor # tempering factor < abs(0.18)
    scale = fm.createFieldConstant([0.0, k1, 0.0])
    scaleFunction = fm.createFieldMultiply(coordinates, scale)
    constant = fm.createFieldConstant([0.0, 1.0, 1.0])
    constantFunction = fm.createFieldAdd(scaleFunction, constant)
    transformation_matrix = fm.createFieldComponent(constantFunction, [2, 3, 3])
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translate_coordinates = fm.createFieldSubtract(taper_coordinates, offset)
    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def taperingZAxis(sharpeningFactor, fm, coordinates, lungNodesetGroup, spaceFromCentre):
    """
    Linear transformation
    :param sharpRadius:
    :param fm:
    :param coordinates:
    :param lungNodesetGroup:
    :param spaceFromCentre:
    :return:
    """
    # Transformation matrix = [     1, | [x,
    #                         k1y + 1, |  y,
    #                               1] |  z]
    offset = fm.createFieldConstant(spaceFromCentre)
    origin = fm.createFieldAdd(coordinates, offset)
    k1 = -sharpeningFactor  # tempering factor < abs(0.18)
    scale = fm.createFieldConstant([0.0, 0.0, k1])
    scaleFunction = fm.createFieldMultiply(coordinates, scale)
    constant = fm.createFieldConstant([0.0, 1.0, 1.0])
    constantFunction = fm.createFieldAdd(scaleFunction, constant)
    transformation_matrix = fm.createFieldComponent(constantFunction, [3, 3, 2])
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translate_coordinates = fm.createFieldSubtract(taper_coordinates, offset)
    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

