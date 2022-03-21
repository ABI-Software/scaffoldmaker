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
            'Rat 1',
            'Pig 1',
            'Material'
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Human 1'

        options = {
            'Base parameter set': parameterSetName,
            'Number of left lung lobes': 2,
            'Left-right lung spacing': 1.0,
            'Left lung width': 0.25,
            'Left lung depth': 0.5,
            'Left lung height': 1.0,
            'Left lung ventral edge sharpness factor': 0.0,
            'Left lung medial curvature': 0.0,
            'Right lung width': 0.25,
            'Right lung depth': 0.5,
            'Right lung height': 1.0,
            'Right lung ventral edge sharpness factor': 0.0,
            'Right lung medial curvature': 0.0,
            'Open fissures': False,
            'Accessory lobe': True,
            'Accessory lobe medial curvature': 0.0,
            'Accessory lobe length': 0.5,
            'Accessory lobe dorsal centre [x,y]': [0.0, 0.0],
            'Accessory lobe dorsal height': 0.5,
            'Accessory lobe dorsal width': 0.5,
            'Accessory lobe ventral height': 0.5,
            'Accessory lobe ventral width': 0.5,
            'Diaphragm centre x': 0.0,
            'Diaphragm centre y': 0.0,
            'Diaphragm curvature x': 0.0,
            'Diaphragm curvature y': 0.0,
            'Refine': False,
            'Refine number of elements': 4
        }

        if 'Human 1' in parameterSetName:
            options['Left-right lung spacing'] = 0.7
            options['Left lung ventral edge sharpness factor'] = 1.0
            options['Left lung medial curvature'] = 1.6
            options['Right lung ventral edge sharpness factor'] = 0.6
            options['Right lung medial curvature'] = 2.5
            options['Accessory lobe'] = False
            options['Diaphragm curvature x'] = 1.0
            options['Diaphragm curvature y'] = 1.0

        elif 'Mouse 1' in parameterSetName:
            options['Number of left lung lobes'] = 1
            options['Left lung height'] = 0.8
            options['Left lung ventral edge sharpness factor'] = 0.6
            options['Left lung medial curvature'] = 2.0
            options['Right lung height'] = 0.8
            options['Right lung ventral edge sharpness factor'] = 0.6
            options['Right lung medial curvature'] = 2.0
            options['Accessory lobe ventral width'] = 0.4
            options['Accessory lobe ventral height'] = 0.4
            options['Accessory lobe medial curvature'] = 0.5
            options['Diaphragm curvature y'] = 1.0

        elif 'Rat 1' in parameterSetName:
            options['Number of left lung lobes'] = 1
            options['Left-right lung spacing'] = 0.9
            options['Left lung ventral edge sharpness factor'] = 0.6
            options['Left lung medial curvature'] = 2.0
            options['Right lung ventral edge sharpness factor'] = 0.6
            options['Right lung medial curvature'] = 0.6
            options['Accessory lobe'] = True
            options['Accessory lobe dorsal centre [x,y]'] = [0.0, 0.2]
            options['Accessory lobe dorsal height'] = 0.4
            options['Accessory lobe dorsal width'] = 0.4
            options['Accessory lobe ventral width'] = 0.2
            options['Accessory lobe ventral height'] = 0.2
            options['Accessory lobe medial curvature'] = 1.11
            options['Diaphragm curvature x'] = -0.4

        elif 'Pig 1' in parameterSetName:
            options['Number of left lung lobes'] = 1
            options['Left-right lung spacing'] = 0.9
            options['Left lung ventral edge sharpness factor'] = 0.6
            options['Left lung medial curvature'] = 2.0
            options['Right lung ventral edge sharpness factor'] = 0.6
            options['Right lung medial curvature'] = 2
            options['Accessory lobe'] = True
            options['Accessory lobe dorsal centre [x,y]'] = [0.0, 0.2]
            options['Accessory lobe dorsal height'] = 0.5
            options['Accessory lobe dorsal width'] = 0.3
            options['Accessory lobe ventral width'] = 0.5
            options['Accessory lobe ventral height'] = 0.2
            options['Accessory lobe medial curvature'] = 0.9
            options['Accessory lobe length'] = 0.5
            options['Diaphragm curvature y'] = 1.0

        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Number of left lung lobes',
            'Left-right lung spacing',
            'Left lung width',
            'Left lung depth',
            'Left lung height',
            'Left lung ventral edge sharpness factor',
            'Left lung medial curvature',
            'Right lung width',
            'Right lung depth',
            'Right lung height',
            'Right lung ventral edge sharpness factor',
            'Right lung medial curvature',
            'Open fissures',
            'Accessory lobe',
            'Accessory lobe dorsal centre [x,y]',
            'Accessory lobe length',
            'Accessory lobe dorsal height',
            'Accessory lobe dorsal width',
            'Accessory lobe ventral height',
            'Accessory lobe ventral width',
            'Accessory lobe medial curvature',
            'Diaphragm centre x',
            'Diaphragm centre y',
            'Diaphragm curvature x',
            'Diaphragm curvature y',
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

        if options['Refine number of elements'] < 1:
            options['Refine number of elements'] = 1

        if options['Number of left lung lobes'] > 2:
            options['Number of left lung lobes'] = 2

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''
        # Generate two meshes: geometric[0] and lung[1] coordinates
        for coordinate in range(2):
            numberOfLeftLung = options['Number of left lung lobes']
            spacingBetweenLeftRight = options['Left-right lung spacing']/2.0 if coordinate == 0 else 0.5
            leftWidth = options['Left lung width'] if coordinate == 0 else 0.25
            leftDepth = options['Left lung depth'] if coordinate == 0 else 0.5
            leftHeight = options['Left lung height'] if coordinate == 0 else 1.0
            leftEdgeSharpFactor = options['Left lung ventral edge sharpness factor'] if coordinate == 0 else 0.0
            leftLungMedialCurvature = options['Left lung medial curvature'] if coordinate == 0 else 0.0
            rightWidth = options['Right lung width'] if coordinate == 0 else 0.25
            rightDepth = options['Right lung depth'] if coordinate == 0 else 0.5
            rightHeight = options['Right lung height'] if coordinate == 0 else 1.0
            rightEdgeSharpFactor = options['Right lung ventral edge sharpness factor'] if coordinate == 0 else 0.0
            rightLungMedialCurvature = options['Right lung medial curvature'] if coordinate == 0 else 0.0
            isOpenfissure = options['Open fissures']
            hasAccessoryLobe = options['Accessory lobe']
            accessoryLobeDorsalCentre = options['Accessory lobe dorsal centre [x,y]'] if coordinate == 0 else [0.0, 0.0]
            accessoryLobeLength = options['Accessory lobe length'] if coordinate == 0 else 0.5
            accessoryLobeDorsalHeight = options['Accessory lobe dorsal height'] if coordinate == 0 else 0.5
            accessoryLobeDorsalWidth = options['Accessory lobe dorsal width'] if coordinate == 0 else 0.5
            accessoryLobeVentralHeight = options['Accessory lobe ventral height'] if coordinate == 0 else 0.5
            accessoryLobeVentralWidth = options['Accessory lobe ventral width'] if coordinate == 0 else 0.5
            accessoryLobeMedialCurve = options['Accessory lobe medial curvature'] if coordinate == 0 else 0.0
            diaphragmCentreX = options['Diaphragm centre x'] if coordinate == 0 else 0.0
            diaphragmCentreY = options['Diaphragm centre y'] if coordinate == 0 else 0.0
            diaphragmCurvatureX = options['Diaphragm curvature x'] if coordinate == 0 else 0.0
            diaphragmCurvatureY = options['Diaphragm curvature y'] if coordinate == 0 else 0.0

            if coordinate == 0:
                fm_n0 = region.getFieldmodule()
                coordinates_n0 = findOrCreateFieldCoordinates(fm_n0)
                # point reference
                fm = fm_n0
                coordinates = coordinates_n0
            else:
                region_n1 = region.createRegion()
                fm_n1 = region_n1.getFieldmodule()
                coordinates_n1 = findOrCreateFieldCoordinates(fm_n1, name="lung coordinates")
                # point reference
                fm = fm_n1
                coordinates = coordinates_n1

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

            ####################
            # Annotation groups
            ####################
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

            if numberOfLeftLung == 2:
                lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
                lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
                annotationGroups.append(lowerLeftLungGroup)
                upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
                upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
                annotationGroups.append(upperLeftLungGroup)

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

            if hasAccessoryLobe:
                # Annotation groups
                rightLungAccessoryLobeGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
                rightLungAccessoryLobeMeshGroup = rightLungAccessoryLobeGroup.getMeshGroup(mesh)
                annotationGroups.append(rightLungAccessoryLobeGroup)
                rightLungAccessoryLobeNodesetGroup = rightLungAccessoryLobeGroup.getNodesetGroup(nodes)

                # Marker points
                accessoryApexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                        get_lung_term("apex of right lung accessory lobe"))
                accessoryVentralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                           get_lung_term("ventral base of right lung accessory lobe"))
                accessoryDorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                          get_lung_term("dorsal base of right lung accessory lobe"))

            # Nodeset group
            leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
            rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)
            lungNodesetGroup = lungGroup.getNodesetGroup(nodes)
            mediastinumLeftNodesetGroup = mediastinumLeftGroup.getNodesetGroup(nodes)
            mediastinumRightNodesetGroup = mediastinumRightGroup.getNodesetGroup(nodes)

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

            # Set internal parameters
            # The number of the elements in the generic lungs
            # These counts are only values that work for nodeFieldParameters (KEEP THEM FIXED)
            lElementsCount1 = 2
            lElementsCount2 = 4
            lElementsCount3 = 3

            uElementsCount1 = 2
            uElementsCount2 = 4
            uElementsCount3 = 4

            leftFissureAngle = math.atan(leftHeight/(leftDepth * 2.0)) * 180 / math.pi
            leftObliqueProportion = 0.8
            rightFissureAngle = math.atan(rightHeight/(rightDepth * 2.0)) * 180 / math.pi
            rightObliqueProportion = 0.8

            ###############
            # Create nodes
            ###############
            nodeIndex = 0
            nodeIdentifier = 1
            lowerLeftNodeIds = []
            upperLeftNodeIds = []
            lowerRightNodeIds = []
            upperRightNodeIds = []

            # Left lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spacingBetweenLeftRight,
                leftDepth, leftWidth, leftHeight, leftFissureAngle, leftObliqueProportion,
                leftLung, cache, coordinates, nodes, nodetemplate,
                mediastinumLeftNodesetGroup, leftLungNodesetGroup, lungNodesetGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerLeftNodeIds, upperLeftNodeIds, nodeIndex, nodeIdentifier)

            # Right lung nodes
            nodeIndex, nodeIdentifier = getLungNodes(spacingBetweenLeftRight,
                rightDepth, rightWidth, rightHeight, rightFissureAngle, rightObliqueProportion,
                rightLung, cache, coordinates, nodes, nodetemplate,
                mediastinumLeftNodesetGroup, rightLungNodesetGroup, lungNodesetGroup,
                lElementsCount1, lElementsCount2, lElementsCount3,
                uElementsCount1, uElementsCount2, uElementsCount3,
                lowerRightNodeIds, upperRightNodeIds, nodeIndex, nodeIdentifier)

            if hasAccessoryLobe:
                # The number of the elements in the diaphragmatic animal lung
                accesssoryElementsCount1 = 2
                accesssoryElementsCount2 = 5
                accesssoryElementsCount3 = 2
                diaphragmaticNodeIds = []

                # Diaphragm lung nodes
                nodeIndex, nodeIdentifier = getDiaphragmaticLungNodes(accessoryLobeDorsalCentre, cache, coordinates, nodes, nodetemplate,
                    rightLungAccessoryLobeNodesetGroup, lungNodesetGroup,
                    accesssoryElementsCount1, accesssoryElementsCount2, accesssoryElementsCount3,
                    accessoryLobeLength, accessoryLobeDorsalWidth, accessoryLobeDorsalHeight, accessoryLobeVentralWidth,
                    accessoryLobeVentralHeight, diaphragmaticNodeIds, nodeIndex, nodeIdentifier)

            ##################
            # Create elements
            ##################
            elementIdentifier = 1

            if numberOfLeftLung == 2:
                # Left lung elements
                elementIdentifier, leftUpperLobeElementID, leftLowerLobeElementID = \
                    getLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                    elementtemplateCustom, mesh, lungMeshGroup,
                    leftLungMeshGroup, lowerLeftLungMeshGroup, None,
                    upperLeftLungMeshGroup, mediastinumLeftGroupMeshGroup,
                    lElementsCount1, lElementsCount2, lElementsCount3,
                    uElementsCount1, uElementsCount2, uElementsCount3,
                    lowerLeftNodeIds, upperLeftNodeIds, elementIdentifier)
            else:
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

            if isOpenfissure:
                # create discontinuity in d3 on the core boundary
                nodeIdentifier = createDiscontinuity(coordinates, nodes, mesh, cache, nodetemplate, leftUpperLobeElementID,
                                rightUpperLobeElementID, nodeIdentifier)

            if hasAccessoryLobe:
                # Diaphragm lung elements
                getDiaphragmaticLungElements(coordinates, eftfactory, eftRegular, elementtemplateRegular,
                    elementtemplateCustom, mesh, lungMeshGroup,
                    rightLungMeshGroup, rightLungAccessoryLobeMeshGroup,
                    accesssoryElementsCount1, accesssoryElementsCount2, accesssoryElementsCount3,
                    diaphragmaticNodeIds, elementIdentifier)

            ########################
            # Write two coordinates
            ########################
            if coordinate == 1:
                sir = region_n1.createStreaminformationRegion()
                srm = sir.createStreamresourceMemory()
                region_n1.write(sir)
                result, buffer = srm.getBuffer()

                sir = region.createStreaminformationRegion()
                srm = sir.createStreamresourceMemoryBuffer(buffer)
                region.read(sir)

                del srm
                del sir
                del fm_n1
                del coordinates_n1
                del region_n1

            ####################################################################
            # Transformation and translation to the left (0) and right lungs (1)
            ####################################################################
            for i in range(2):
                lungNodeset = leftLungNodesetGroup if i == 0 else rightLungNodesetGroup
                edgeSharpFactor = leftEdgeSharpFactor if i == 0 else rightEdgeSharpFactor
                length = leftDepth if i == 0 else rightDepth
                spacing = spacingBetweenLeftRight if i == 0 else -spacingBetweenLeftRight
                lungMedialcurvature = -leftLungMedialCurvature if i == 0 else rightLungMedialCurvature

                # Transformation of the left and right lungs
                if edgeSharpFactor != 0.0:
                    sharpeningRidge(edgeSharpFactor, fm, coordinates, lungNodeset, spacing, length)

                if lungMedialcurvature != 0.0:
                    bendingAroundZAxis(lungMedialcurvature, fm, coordinates, lungNodeset, spacing)

            if accessoryLobeMedialCurve != 0.0:
                spacing = accessoryLobeDorsalCentre
                bendingAroundZAxis(-accessoryLobeMedialCurve, fm, coordinates, rightLungAccessoryLobeNodesetGroup, spacing)

            if (diaphragmCurvatureX != 0.0) or (diaphragmCurvatureY != 0.0):
                concavingDiaphragmaticSurface(diaphragmCurvatureX, diaphragmCurvatureY, fm, coordinates, diaphragmCentreX,
                                              diaphragmCentreY, lungNodesetGroup)

            ###############
            # Marker points
            ###############
            markerList = []

            lowerLeftElementCount = (lElementsCount1 * (lElementsCount2-1) * lElementsCount3 + lElementsCount1)

            idx = lowerLeftElementCount + (uElementsCount1 * uElementsCount2 * (uElementsCount3//2) + uElementsCount2)
            markerList.append({ "group" : leftApexGroup, "elementId" : idx, "xi" : [0.0, 1.0, 1.0]})

            idx = 1
            markerList.append({"group": leftDorsalGroup, "elementId": idx, "xi": [0.0, 0.0, 0.0]})

            idx = lElementsCount1 * (lElementsCount2 // 2)
            markerList.append({"group": leftMedialGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

            idx = lowerLeftElementCount + (uElementsCount1 * (uElementsCount2 - 2))
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

            idx = leftLungElementCount + lowerRightElementCount + (uElementsCount1 * (uElementsCount2 - 2))
            markerList.append({"group": rightVentralGroup, "elementId": idx, "xi": [0.0, 1.0, 0.0]})

            idx = leftLungElementCount + (lElementsCount1 * lElementsCount2 * lElementsCount3) + lElementsCount1
            markerList.append({"group": rightLateralGroup, "elementId": idx, "xi": [1.0, 0.0, 1.0]})

            if hasAccessoryLobe:
                rightLungElementCount = leftLungElementCount

                idx_temp = accesssoryElementsCount1 * accesssoryElementsCount2 * (accesssoryElementsCount3 - 1) + 1
                idx = rightLungElementCount + leftLungElementCount + idx_temp
                markerList.append({"group": accessoryApexGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

                idx_temp = accesssoryElementsCount1 * accesssoryElementsCount2 * (accesssoryElementsCount3 - 1)
                idx = rightLungElementCount + leftLungElementCount + idx_temp
                markerList.append({"group": accessoryVentralGroup, "elementId": idx, "xi": [1.0, 1.0, 0.0]})

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
        isPig = 'Pig 1' in parameterSetName
        isMaterial = 'Material' in parameterSetName

        # create fissure groups
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        # 1D Annotation
        is_exterior = fm.createFieldIsExterior()
        is_xi1_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)
        is_xi1_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)
        is_xi3_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)

        # edge markers
        mediastanumLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("anterior mediastinum of left lung"))
        is_anteriorBorderLeftGroup = mediastanumLeftGroup.getFieldElementGroup(mesh1d)

        is_mediastanumLeftGroup_exterior = fm.createFieldAnd(is_anteriorBorderLeftGroup, is_exterior)
        is_mediastanumLeftGroup_exterior_xi1_0 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior, is_xi1_0)
        is_mediastanumLeftGroup_exterior_xi1_01 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi1_0, is_xi1_1)
        is_mediastanumLeftGroup_exterior_xi3_0 = fm.createFieldAnd(is_mediastanumLeftGroup_exterior, is_xi3_0)

        is_upperLeftGroup_Diaphragm = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi3_0, is_mediastanumLeftGroup_exterior_xi1_01)

        is_ridge = fm.createFieldAnd(is_mediastanumLeftGroup_exterior_xi1_01, fm.createFieldNot(is_upperLeftGroup_Diaphragm))
        anteriorBorderLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("anterior border of left lung"))

        anteriorBorderLeftGroup.getMeshGroup(mesh1d).addElementsConditional(is_ridge)

        mediastanumRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("anterior mediastinum of right lung"))
        is_anteriorBorderRightGroup = mediastanumRightGroup.getFieldElementGroup(mesh1d)

        is_anteriorBorderRightGroup_exterior = fm.createFieldAnd(is_anteriorBorderRightGroup, is_exterior)
        is_anteriorBorderRightGroup_exterior_xi1_0 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior, is_xi1_0)
        is_anteriorBorderRightGroup_exterior_xi1_01 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi1_0, is_xi1_1)
        is_anteriorBorderRightGroup_exterior_xi3_0 = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior, is_xi3_0)

        is_upperRightGroup_Diaphragm = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi3_0, is_anteriorBorderRightGroup_exterior_xi1_01)

        is_ridge = fm.createFieldAnd(is_anteriorBorderRightGroup_exterior_xi1_01, fm.createFieldNot(is_upperRightGroup_Diaphragm))
        anteriorBorderRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("anterior border of right lung"))

        anteriorBorderRightGroup.getMeshGroup(mesh1d).addElementsConditional(is_ridge)

        # exterior surfaces
        leftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("left lung"))
        is_leftGroup = leftGroup.getFieldElementGroup(mesh2d)
        is_LeftGroup_exterior = fm.createFieldAnd(is_leftGroup, is_exterior)
        leftSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("left lung surface"))
        leftSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_LeftGroup_exterior)

        is_LeftGroup_exterior_xi3_0 = fm.createFieldAnd(is_LeftGroup_exterior, is_xi3_0)
        baseOfLeftSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("base of left lung surface"))
        baseOfLeftSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_LeftGroup_exterior_xi3_0)

        rightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("right lung"))
        is_rightGroup = rightGroup.getFieldElementGroup(mesh2d)
        is_rightGroup_exterior = fm.createFieldAnd(is_rightGroup, is_exterior)
        rightSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("right lung surface"))
        rightSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_rightGroup_exterior)

        upperRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of right lung"))
        is_upperRightGroup = upperRightGroup.getFieldElementGroup(mesh2d)
        is_upperRightGroup_exterior = fm.createFieldAnd(is_upperRightGroup, is_exterior)
        upperRightSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("upper lobe of right lung surface"))
        upperRightSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_upperRightGroup_exterior)

        middleRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("middle lobe of right lung"))
        is_middleRightGroup = middleRightGroup.getFieldElementGroup(mesh2d)
        is_middleRightGroup_exterior = fm.createFieldAnd(is_middleRightGroup, is_exterior)
        middleRightSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("middle lobe of right lung surface"))
        middleRightSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_middleRightGroup_exterior)

        lowerRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of right lung"))
        is_lowerRightGroup = lowerRightGroup.getFieldElementGroup(mesh2d)
        is_lowerRightGroup_exterior = fm.createFieldAnd(is_lowerRightGroup, is_exterior)
        lowerRightSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("lower lobe of right lung surface"))
        lowerRightSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_lowerRightGroup_exterior)

        is_rightGroup_exterior_xi3_0 = fm.createFieldAnd(is_rightGroup_exterior, is_xi3_0)
        baseOfRightSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("base of right lung surface"))
        baseOfRightSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_rightGroup_exterior_xi3_0)

        upperRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of right lung"))
        middleRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("middle lobe of right lung"))
        lowerRightGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of right lung"))

        is_upperRightGroup = upperRightGroup.getFieldElementGroup(mesh2d)
        is_middleRightGroup = middleRightGroup.getFieldElementGroup(mesh2d)
        is_lowerRightGroup = lowerRightGroup.getFieldElementGroup(mesh2d)

        is_obliqueRightGroup = fm.createFieldAnd(fm.createFieldOr(is_middleRightGroup, is_upperRightGroup), is_lowerRightGroup)
        is_horizontalRightGroup = fm.createFieldAnd(is_upperRightGroup, is_middleRightGroup)

        obliqueRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of right lung"))
        obliqueRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueRightGroup)
        horizontalRightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("horizontal fissure of right lung"))
        horizontalRightGroup.getMeshGroup(mesh2d).addElementsConditional(is_horizontalRightGroup)

        if isHuman or isMaterial:
            lowerLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("lower lobe of left lung"))
            is_lowerLeftGroup = lowerLeftGroup.getFieldElementGroup(mesh2d)
            is_lowerLeftGroup_exterior = fm.createFieldAnd(is_lowerLeftGroup, is_exterior)
            lowerLeftSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term("lower lobe of left lung surface"))
            lowerLeftSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_lowerLeftGroup_exterior)

            upperLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("upper lobe of left lung"))
            is_upperLeftGroup = upperLeftGroup.getFieldElementGroup(mesh2d)
            is_upperLeftGroup_exterior = fm.createFieldAnd(is_upperLeftGroup, is_exterior)
            upperLeftSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term("upper lobe of left lung surface"))
            upperLeftSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_upperLeftGroup_exterior)

            is_upperLeftGroup = upperLeftGroup.getFieldElementGroup(mesh2d)
            is_lowerLeftGroup = lowerLeftGroup.getFieldElementGroup(mesh2d)

            is_obliqueLeftGroup = fm.createFieldAnd(is_upperLeftGroup, is_lowerLeftGroup)
            obliqueLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("oblique fissure of left lung"))
            obliqueLeftGroup.getMeshGroup(mesh2d).addElementsConditional(is_obliqueLeftGroup)

        if isMouse or isRat or isMaterial or isPig:
            rightAccessoryGroup = getAnnotationGroupForTerm(annotationGroups, get_lung_term("right lung accessory lobe"))
            is_rightAccessoryGroup = rightAccessoryGroup.getFieldElementGroup(mesh2d)
            is_rightAccessoryGroup_exterior = fm.createFieldAnd(is_rightAccessoryGroup, is_exterior)
            rightAccessorySurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                            get_lung_term("right lung accessory lobe surface"))
            rightAccessorySurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_rightAccessoryGroup_exterior)

            is_rightAccessoryGroup_exterior_xi3_0 = fm.createFieldAnd(is_rightAccessoryGroup_exterior, is_xi3_0)
            baseOfRightAccessorySurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                                  get_lung_term("base of right lung accessory lobe surface"))
            baseOfRightAccessorySurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_rightAccessoryGroup_exterior_xi3_0)

def getLungNodes(spaceFromCentre, lengthUnit, widthUnit, heightUnit, fissureAngle, obliqueProportion,
                 lungSide, cache, coordinates, nodes, nodetemplate,
                 mediastinumNodesetGroup, lungSideNodesetGroup, lungNodesetGroup,
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
                    x = [x[0] + spaceFromCentre, x[1], x[2]]
                else:
                    x = [x[0] - spaceFromCentre, x[1], x[2]]

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
                lungNodesetGroup.addNode(node)

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
                    x = [x[0] + spaceFromCentre, x[1], x[2]]
                else:
                    x = [x[0] - spaceFromCentre, x[1], x[2]]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                upperNodeIds[i][j][k] = nodeIdentifier
                nodeIdentifier += 1

                # Annotation
                if j > 2:
                    mediastinumNodesetGroup.addNode(node)
                lungSideNodesetGroup.addNode(node)
                lungNodesetGroup.addNode(node)

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

def getDiaphragmaticLungNodes(centre, cache, coordinates, nodes, nodetemplate, lungSideNodesetGroup, lungNodesetGroup,
                              elementsCount1, elementsCount2, elementsCount3,
                              length, dorsalWidth, dorsalHeight, ventralWidth, ventralHeight,
                              nodeIds, nodeIndex, nodeIdentifier):
    """
    Create a 3D triangular mesh from getDiaphragmaticLungNodes
    :parameter: elementsCount1 - x, elementsCount2 - y, elementsCount3 - z
    :return: nodeIndex, nodeIdentifier
    """

    # Initialise parameters
    lengthUnit = length/elementsCount2

    # Determine a triangular face at the ventral and dorsal ends
    dorsalHeightUnit = dorsalHeight/elementsCount3
    dorsalWidthUnit = dorsalWidth/elementsCount1
    ventralHeightUnit = ventralHeight/elementsCount3
    ventralWidthUnit = ventralWidth/elementsCount1

    x_dorsal = []
    for n3 in range(elementsCount3 + 1):
        x_dorsal.append([])
        for n2 in range(3):
            x_dorsal[n3].append([])
            width = dorsalWidthUnit if n2 == 0 else ventralWidthUnit
            height = dorsalHeightUnit if n2 == 0 else ventralHeightUnit
            for n1 in range(elementsCount1 + 1):
                x_dorsal[n3][n2].append(None)
                if (n1 != 1) and (n3 == elementsCount3):
                    continue
                unitRatio = 0.5 if n3 == 1 else 1.0
                x_dorsal[n3][n2][n1] = [width * unitRatio * (n1 - 1), lengthUnit * (n2*elementsCount2), height * n3]

    # Interpolate triangular faces between two ends
    px = []
    pd2 = []
    for n3 in range(elementsCount3 + 1):
        px.append([])
        pd2.append([])
        for n1 in range(elementsCount1 + 1):
            px[n3].append(None)
            pd2[n3].append(None)
            x_temp = [x_dorsal[n3][0][n1], x_dorsal[n3][1][n1]]
            if None in x_temp:
                continue
            d2_temp = [(x_dorsal[n3][1][n1][i] - x_dorsal[n3][0][n1][i])/elementsCount2 for i in range(3)]
            px[n3][n1], pd2[n3][n1], pe, pxi, psf = sampleCubicHermiteCurves(x_temp, [d2_temp, d2_temp], elementsCount2)

    # Diaphragmatic lobe nodes
    for n3 in range(elementsCount3 + 1):
        nodeIds.append([])
        for n2 in range(elementsCount2 + 1):
            nodeIds[n3].append([])
            for n1 in range(elementsCount1 + 1):
                nodeIds[n3][n2].append(None)

                if (n1 != 1) and (n3 == elementsCount3):
                    continue

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)

                x = px[n3][n1][n2]
                d2 = pd2[n3][n1][n2]

                if n3 < elementsCount3:
                    d1 = [px[n3][n1][n2][i] - px[n3][n1-1][n2][i] for i in range(3)] if n1 == elementsCount1 else \
                    [px[n3][n1+1][n2][i] - px[n3][n1][n2][i] for i in range(3)]
                    d3 = [px[n3+1][1][n2][i] - px[n3][n1][n2][i] for i in range(3)] if n3 == 1 else \
                        [px[n3 + 1][n1][n2][i] - px[n3][n1][n2][i] for i in range(3)]
                else:
                    d1 = [px[n3-1][n1][n2][i] - px[n3-1][n1-1][n2][i] for i in range(3)]
                    d3 = [px[n3][n1][n2][i] - px[n3-1][n1][n2][i] for i in range(3)]

                # Translate the mesh
                x = [x[0] + centre[0], x[1]+centre[1], x[2]]

                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                nodeIds[n3][n2][n1] = nodeIdentifier
                nodeIdentifier += 1

                # Annotation
                lungSideNodesetGroup.addNode(node)
                lungNodesetGroup.addNode(node)

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

def concavingDiaphragmaticSurface(diaphragmCurvatureX, diaphragmCurvatureY, fm, coordinates, diaphragmCentreX,
                                  diaphragmCentreY, lungNodesetGroup):
    """
    # s_x = x - x_centre                        s_y = y - y_centre
    # kappa_x = diaphragmCurvatureX             kappa_y = diaphragmCurvatureY
    # theta_x = kappa_x * s_x                   theta_y = kappa_y * s_y
    # r_x = 1/kappa_x                           r_y = 1/kappa_y
    # s_zx = z + r_x                            s_zy = z + r_y
    # delta_zx = s_zx * (cos(theta_x) - 1.0)    delta_zy = s_zy * (cos(theta_y) - 1.0)
    # ----------------------------------------------------------------------------------------
    # x_new = s_zx * sin(theta_x)
    # y_new = s_zy * sin(theta_y)
    # z_new = z + deta_zx + deta_zy

    :param diaphragmCurvatureX:
    :param diaphragmCurvatureY:
    :param fm:
    :param coordinates:
    :param diaphragmCentre:
    :param lungNodesetGroup:
    :return:
    """

    # Initialise parameters
    diaphragmCentre = fm.createFieldConstant([diaphragmCentreX, diaphragmCentreY, 0.0])
    offset_coordinates = fm.createFieldSubtract(coordinates, diaphragmCentre)
    only_x = fm.createFieldConstant([1.0, 0.0, 0.0])
    only_y = fm.createFieldConstant([0.0, 1.0, 0.0])
    only_z = fm.createFieldConstant([0.0, 0.0, 1.0])

    x_new    = fm.createFieldMultiply(coordinates, only_x)
    y_new    = fm.createFieldMultiply(coordinates, only_y)
    delta_zx = fm.createFieldConstant([0.0, 0.0, 0.0])
    delta_zy = fm.createFieldConstant([0.0, 0.0, 0.0])

    # x-coordinates
    if diaphragmCurvatureX != 0.0:
        kappa_x = fm.createFieldConstant([diaphragmCurvatureX, 0.0, 0.0])
        r_x = fm.createFieldConstant([1/diaphragmCurvatureX, 0.0, 0.0])
        s_x = fm.createFieldMultiply(offset_coordinates, only_x)
        z_x = fm.createFieldMultiply(coordinates, only_z)
        z_x = fm.createFieldComponent(z_x, [3, 1, 1])
        s_zx = r_x # if no bulge s_zx = r_x
        theta_x = fm.createFieldMultiply(kappa_x, s_x)
        x_new = fm.createFieldMultiply(s_zx, fm.createFieldSin(theta_x))
        delta_zx = fm.createFieldMultiply(s_zx, fm.createFieldSubtract(fm.createFieldCos(theta_x),
                                                                        fm.createFieldConstant([1.0, 0.0, 0.0])))
        delta_zx = fm.createFieldComponent(delta_zx, [3, 3, 1])

    # y-coordinates
    if diaphragmCurvatureY != 0.0:
        kappa_y = fm.createFieldConstant([0.0, diaphragmCurvatureY, 0.0])
        r_y = fm.createFieldConstant([0.0, 1/diaphragmCurvatureY, 0.0])
        s_y = fm.createFieldMultiply(offset_coordinates, only_y)
        z_y = fm.createFieldMultiply(coordinates, only_z)
        z_y = fm.createFieldComponent(z_y, [1, 3, 1])
        s_zy = r_y # if no bulge s_zx = r_y
        theta_y = fm.createFieldMultiply(kappa_y, s_y)
        y_new = fm.createFieldMultiply(s_zy, fm.createFieldSin(theta_y))
        delta_zy = fm.createFieldMultiply(s_zy, fm.createFieldSubtract(fm.createFieldCos(theta_y),
                                                                       fm.createFieldConstant([0.0, 1.0, 0.0])))
        delta_zy = fm.createFieldComponent(delta_zy, [3, 3, 2])

    # z-coordinates
    z = fm.createFieldMultiply(coordinates, only_z)
    z_new = fm.createFieldAdd(delta_zx, delta_zy)
    z_new = fm.createFieldAdd(z_new, z)

    new_coordinates = fm.createFieldAdd(x_new, y_new)
    new_coordinates = fm.createFieldAdd(new_coordinates, z_new)

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def bendingAroundZAxis(curvature, fm, coordinates, lungNodesetGroup, spaceFromCentre):
    """
    :param bulgeRadius: the radius and the centre of curvature to transfrom the scaffold
    :param fm:
    :param coordinates:
    :param nodes:
    :return:
    """
    # cylindrical polar coordinates (x = r*cos(theta), y = r*sin(theta), z = z):
    # r = x - bulgeRadius
    # theta = y / bulgeRadius
    # z = z

    radius = 1/curvature
    scale = fm.createFieldConstant([1.0, curvature, 1.0])
    scaleCoordinates = fm.createFieldMultiply(coordinates, scale)
    if isinstance(spaceFromCentre, list):
        offset = fm.createFieldConstant([radius + spaceFromCentre[0], spaceFromCentre[1], 0.0])
    else:
        offset = fm.createFieldConstant([radius + spaceFromCentre, 0.0, 0.0])
    polarCoordinates = fm.createFieldAdd(scaleCoordinates, offset)
    polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
    rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
    rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    newxyzCoordinates = fm.createFieldSubtract(rcCoordinates, offset)
    fieldassignment = coordinates.createFieldassignment(newxyzCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

def sharpeningRidge(sharpeningFactor, fm, coordinates, lungNodesetGroup, spaceFromCentre, length):
    """
    Linear transformation
    :param sharpRadius:
    :param fm:
    :param coordinates:
    :param lungNodesetGroup:
    :param spaceFromCentre:
    :return:
    """
    # Transformation matrix = [ k1y + 1, | [x,
    #                                 1, |  y,
    #                                 1] |  z]
    offset = fm.createFieldConstant([spaceFromCentre, 0.75 * length, 0.0])
    origin = fm.createFieldAdd(coordinates, offset)
    k1 = -sharpeningFactor * 1.14 # 1.14 is the maximum value but the input range 0-1
    scale = fm.createFieldConstant([0.0, k1, 0.0])
    scaleFunction = fm.createFieldMultiply(origin, scale)
    constant = fm.createFieldConstant([0.0, 1.0, 1.0])
    constantFunction = fm.createFieldAdd(scaleFunction, constant)
    transformation_matrix = fm.createFieldComponent(constantFunction, [2, 3, 3])
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translate_coordinates = fm.createFieldSubtract(taper_coordinates, offset)
    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()
