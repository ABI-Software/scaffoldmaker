"""
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
"""

from __future__ import division

import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.vector import setMagnitude
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class MeshType_3d_wholebody1(Scaffold_base):
    """
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
    """
    cylinder1Settings = {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.5,
                'Number of elements': 1
            }

    axis1 = [0, 0, 1]
    axis2 = [1, 0, 0]
    axis3 = [0, 1, 0]
    centralPathDefaultScaffoldPackages = {
        'Default': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': cylinder1Settings,
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, 0.0, 0.0], setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]],
                    [setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Whole Body 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Coarse',
            'Fine']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Default']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 6,
            'Number of elements across minor': 6,
            'Number of elements across shell': 1,
            'Number of elements across transition': 1,
            'Number of elements in abdomen': 5,
            'Number of elements in thorax': 3,
            'Number of elements in neck': 1,
            'Number of elements in head': 2,
            'Shell thickness proportion': 0.2,
            'Discontinuity on the core boundary': True,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        if 'Coarse' in parameterSetName:
            pass
        if 'Fine' in parameterSetName:
            options['Number of elements across major'] = 10
            options['Number of elements across minor'] = 10
            options['Number of elements across shell'] = 1
            options['Number of elements across transition'] = 1
            options['Number of elements in abdomen'] = 10
            options['Number of elements in thorax'] = 6
            options['Number of elements in neck'] = 2
            options['Number of elements in head'] = 4
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements across transition',
            'Number of elements in abdomen',
            'Number of elements in thorax',
            'Number of elements in neck',
            'Number of elements in head',
            'Shell thickness proportion',
            'Discontinuity on the core boundary',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] != options['Number of elements across major']:
            options['Number of elements across minor'] = options['Number of elements across major']
            dependentChanges = True

        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1

        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1

        if options['Shell thickness proportion'] < 0.07 or options['Shell thickness proportion'] > 0.7:
            options['Shell thickness proportion'] = 2*options['Number of elements across shell']/options['Number of elements across major']

        if options['Number of elements in abdomen'] < 1:
            options['Number of elements in abdomen'] = 1
        if options['Number of elements in head'] < 1:
            options['Number of elements in head'] = 1
        if options['Number of elements in neck'] < 1:
            options['Number of elements in neck'] = 1
        if options['Number of elements in thorax'] < 1:
            options['Number of elements in thorax'] = 1

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: List of AnnotationGroup
        """

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        elementsCountAlongAbdomen = options['Number of elements in abdomen']
        elementsCountAlongHead = options['Number of elements in head']
        elementsCountAlongNeck = options['Number of elements in neck']
        elementsCountAlongThorax = options['Number of elements in thorax']
        shellRadiusProportion = options['Shell thickness proportion']
        shellProportion = 1/(1/shellRadiusProportion-1)*(elementsCountAcrossMajor/2/elementsCountAcrossShell - 1)
        discontinuity = options['Discontinuity on the core boundary']
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAlong = elementsCountAlongAbdomen + elementsCountAlongThorax + elementsCountAlongNeck + elementsCountAlongHead

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        mesh = fieldmodule.findMeshByDimension(3)

        bodyGroup = AnnotationGroup(region, get_body_term("body"))
        coreGroup = AnnotationGroup(region, get_body_term("core"))
        non_coreGroup = AnnotationGroup(region, get_body_term("non core"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        neckGroup = AnnotationGroup(region, get_body_term("neck core"))
        headGroup = AnnotationGroup(region, get_body_term("head core"))
        annotationGroups = [bodyGroup, coreGroup, non_coreGroup, abdomenGroup, thoraxGroup, neckGroup, headGroup]

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            elementsCountAcrossTransition, shellProportion,
                            [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fieldmodule, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # Groups of different parts of the body
        is_body = fieldmodule.createFieldConstant(1)
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        bodyMeshGroup.addElementsConditional(is_body)

        coreMeshGroup = coreGroup.getMeshGroup(mesh)

        # core group
        e1a = elementsCountAcrossShell
        e1z = elementsCountAcrossMinor - elementsCountAcrossShell - 1
        e2a = elementsCountAcrossShell
        e2b = e2a + elementsCountAcrossTransition
        e2z = elementsCountAcrossMajor - elementsCountAcrossShell - 1
        e2y = e2z - elementsCountAcrossTransition
        e1oc = elementsCountAcrossMinor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oc = elementsCountAcrossMajor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oCore = e2oc * e1oc + 2 * elementsCountAcrossTransition * (e2oc + e1oc)
        elementsCountAround = cylinder1.getElementsCountAround()
        e2oShell = elementsCountAround * elementsCountAcrossShell
        e2o = e2oCore + e2oShell
        elementId = cylinder1.getElementIdentifiers()
        for e3 in range(elementsCountAlong):
            for e2 in range(elementsCountAcrossMajor):
                for e1 in range(elementsCountAcrossMinor):
                    coreElement = ((e2 >= e2a) and (e2 <= e2z)) and ((e1 >= e1a) and (e1 <= e1z))
                    if coreElement:
                        elementIdentifier = elementId[e3][e2][e1]
                        if elementIdentifier:
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            coreMeshGroup.addElement(element)

        is_non_core = fieldmodule.createFieldNot(coreGroup.getGroup())
        non_coreMeshGroup = non_coreGroup.getMeshGroup(mesh)
        non_coreMeshGroup.addElementsConditional(is_non_core)

        abdomenMeshGroup = abdomenGroup.getMeshGroup(mesh)
        thoraxMeshGroup = thoraxGroup.getMeshGroup(mesh)
        neckMeshGroup = neckGroup.getMeshGroup(mesh)
        headMeshGroup = headGroup.getMeshGroup(mesh)
        meshGroups = [abdomenMeshGroup, thoraxMeshGroup, neckMeshGroup, headMeshGroup]

        abdomenRange = [1, elementsCountAlongAbdomen*e2o]
        thoraxRange = [abdomenRange[1]+1, abdomenRange[1]+elementsCountAlongThorax*e2o]
        neckRange = [thoraxRange[1]+1, thoraxRange[1] + elementsCountAlongNeck*e2o]
        headRange = [neckRange[1]+1, elementsCountAlong*e2o]
        groupsRanges = [abdomenRange, thoraxRange, neckRange, headRange]

        totalElements = e2o*elementsCountAlong
        for elementIdentifier in range(1, totalElements+1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            if coreMeshGroup.containsElement(element):
                ri = 0
                for groupRange in groupsRanges:
                    if (elementIdentifier >= groupRange[0]) and (elementIdentifier <= groupRange[1]):
                        meshGroups[ri].addElement(element)
                        break
                    ri += 1

        if discontinuity:
            # create discontinuity in d3 on the core boundary
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            elementtemplate = mesh.createElementtemplate()
            undefineNodetemplate = nodes.createNodetemplate()
            undefineNodetemplate.undefineField(coordinates)
            nodetemplate = nodes.createNodetemplate()
            fieldcache = fieldmodule.createFieldcache()
            with ChangeManager(fieldmodule):
                localNodeIndexes = [1, 2, 3, 4]
                valueLabel = Node.VALUE_LABEL_D_DS3
                for e3 in range(elementsCountAlong):
                    for e2 in range(elementsCountAcrossMajor):
                        for e1 in range(elementsCountAcrossMinor):
                            regularRowElement = (((e2 >= e2b) and (e2 <= e2y)) and ((e1 == e1a - 1) or (e1 == e1z + 1)))
                            non_coreFirstLayerElement = (e2 == e2a - 1) or regularRowElement or (e2 == e2z + 1)
                            elementIdentifier = elementId[e3][e2][e1]
                            if elementIdentifier and non_coreFirstLayerElement:
                                element = mesh.findElementByIdentifier(elementIdentifier)
                                eft = element.getElementfieldtemplate(coordinates, -1)
                                nodeIds = get_element_node_identifiers(element, eft)
                                for localNodeIndex in localNodeIndexes:
                                    node = element.getNode(eft, localNodeIndex)
                                    nodetemplate.defineFieldFromNode(coordinates, node)
                                    versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1, valueLabel)
                                    if versionsCount == 1:
                                        fieldcache.setNode(node)
                                        result0, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                        result0, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                                        result0, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                                        result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1, 3)
                                        result1 = node.merge(undefineNodetemplate)
                                        result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 2)
                                        result3 = node.merge(nodetemplate)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                        result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)
                                result1 = elementtemplate.defineField(coordinates, -1, eft)
                                result2 = element.merge(elementtemplate)
                                result3 = element.setNodesByIdentifier(eft, nodeIds)
        else:
            fieldcache = fieldmodule.createFieldcache()

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fieldmodule, "marker")
        markerName = findOrCreateFieldStoredString(fieldmodule, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fieldmodule, mesh, name="marker_location")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        #
        middleLeft = elementsCountAcrossMinor//2
        topElem = elementsCountAcrossMajor - 1
        middleRight = middleLeft - 1
        neckFirstElem = elementsCountAlongAbdomen+elementsCountAlongThorax
        thoraxFirstElem = elementsCountAlongAbdomen
        middleDown = elementsCountAcrossMajor//2 - 1

        bodyMarkerPoints = [
            {"name": "left hip joint", "elementId": elementId[1][0][middleLeft], "xi": [0.8, 0.5, 0.5]},
            {"name": "right hip joint", "elementId": elementId[1][topElem][middleLeft], "xi": [0.2, 0.5, 0.5]},
            {"name": "left shoulder joint", "elementId": elementId[neckFirstElem][0][middleRight], "xi": [0.8, 0.5, 0.5]},
            {"name": "right shoulder joint", "elementId": elementId[neckFirstElem][topElem][middleRight], "xi": [0.2, 0.5, 0.5]},
            {"name": "along left femur", "elementId": elementId[1][0][middleLeft], "xi": [0.2, 0.99, 0.5]},
            {"name": "along right femur", "elementId": elementId[1][topElem][middleLeft], "xi": [0.8, 0.99, 0.5]},
            {"name": "along left humerus", "elementId": elementId[neckFirstElem][0][middleRight], "xi": [0.5, 0.0, 0.5]},
            {"name": "along right humerus", "elementId": elementId[neckFirstElem][topElem][middleRight], "xi": [0.5, 0.0, 0.5]},
            {"name": "heart apex", "elementId": elementId[thoraxFirstElem][middleDown][middleRight], "xi": [0.5, 0.5, 0.0]},
            {"name": "atrial base", "elementId": elementId[thoraxFirstElem + 1][middleDown][middleRight], "xi": [0.7, 0.5, 0.0]},
            {"name": "aorta", "elementId": elementId[thoraxFirstElem + 1][middleDown][middleRight], "xi": [0.7, 0.5, 0.2]}
        ]

        nodeIdentifier = cylinder1._endNodeIdentifier + 1000
        for bodyMarkerPoint in bodyMarkerPoints:
            element = mesh.findElementByIdentifier(bodyMarkerPoint["elementId"])
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            fieldcache.setNode(markerPoint)
            markerName.assignString(fieldcache, bodyMarkerPoint["name"])
            markerLocation.assignMeshLocation(fieldcache, element, bodyMarkerPoint["xi"])
            nodeIdentifier += 1

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

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

        # create 2d surface mesh groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("body"))
        coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))
        non_coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("non core"))
        abdomenGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdomen"))
        thoraxGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thorax"))
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("neck core"))

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))
        coreBoundaryGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("core boundary"))
        diaphragmGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("diaphragm"))
        spinalCordGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("spinal cord"))

        is_exterior = fm.createFieldIsExterior()
        is_on_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_skin = fm.createFieldAnd(is_exterior, is_on_face_xi3_1)

        skinMeshGroup = skinGroup.getMeshGroup(mesh2d)
        skinMeshGroup.addElementsConditional(is_skin)

        is_core_boundary = fm.createFieldAnd(coreGroup.getGroup(), non_coreGroup.getGroup())
        coreBoundaryMeshGroup = coreBoundaryGroup.getMeshGroup(mesh2d)
        coreBoundaryMeshGroup.addElementsConditional(is_core_boundary)

        is_diaphragm = fm.createFieldAnd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        diaphragmMeshGroup = diaphragmGroup.getMeshGroup(mesh2d)
        diaphragmMeshGroup.addElementsConditional(is_diaphragm)

        # spinal cord
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()
        zero = fm.createFieldConstant(0)
        zero_m = fm.createFieldConstant(-0.01)
        zero_p = fm.createFieldConstant(0.01)
        comp2 = cls.axis2.index(max(cls.axis2)) + 1
        ax2_comp = fm.createFieldComponent(coordinates, comp2)
        ax2_gt_zero_m = fm.createFieldGreaterThan(ax2_comp, zero_m)
        ax2_lt_zero_p = fm.createFieldLessThan(ax2_comp, zero_p)
        ax2_gt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_gt_zero_m)
        ax2_lt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_lt_zero_p)
        is_ax2_zero = fm.createFieldAnd(ax2_lt_zero_xi10, ax2_gt_zero_xi10)
        comp3 = cls.axis3.index(max(cls.axis3)) + 1
        ax3_comp = fm.createFieldComponent(coordinates, comp3)
        ax3_positive = fm.createFieldGreaterThan(ax3_comp, zero)
        is_ax2_zero_ax3_positive = fm.createFieldAnd(is_ax2_zero, ax3_positive)
        is_abdomen_thorax = fm.createFieldAdd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        is_abdomen_thorax_neck = fm.createFieldAdd(is_abdomen_thorax, neckGroup.getGroup())
        is_abdomen_thorax_neck_boundary = fm.createFieldAnd(is_core_boundary, is_abdomen_thorax_neck)
        is_spinal_cord = fm.createFieldAnd(is_ax2_zero_ax3_positive, is_abdomen_thorax_neck_boundary)

        mesh1d = fm.findMeshByDimension(1)
        spinalCordMeshGroup = spinalCordGroup.getMeshGroup(mesh1d)
        spinalCordMeshGroup.addElementsConditional(is_spinal_cord)
