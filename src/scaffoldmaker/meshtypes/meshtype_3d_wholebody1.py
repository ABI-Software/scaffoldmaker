"""
Generates a solid cylinder using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node
from opencmiss.zinc.element import Element
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm, mergeAnnotationGroups
from scaffoldmaker.annotation.body_terms import get_body_term
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelsVersion


class MeshType_3d_wholebody1(Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    cylinder1Settings = {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.5,
                'Number of elements': 1
            }

    centralPathDefaultScaffoldPackages = {
        'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': cylinder1Settings,
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, 0.0, 0.0], [cylinder1Settings['Length'], 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0]],
                    [[cylinder1Settings['Length'], 0.0, 0.0], [cylinder1Settings['Length'], 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Whole Body 1'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 6,
            'Number of elements across minor': 6,
            'Number of elements across shell': 1,
            'Number of elements in abdomen': 5,
            'Number of elements in head': 2,
            'Number of elements in neck': 1,
            'Number of elements in thorax': 3,
            'Shell thickness': 0.17,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements in abdomen',
            'Number of elements in head',
            'Number of elements in neck',
            'Number of elements in thorax',
            'Shell thickness',
            'Lower half',
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

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1

        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit

        if options['Number of elements in abdomen'] < 1:
            options['Number of elements in abdomen'] = 1
        if options['Number of elements in head'] < 1:
            options['Number of elements in head'] = 1
        if options['Number of elements in neck'] < 1:
            options['Number of elements in neck'] = 1
        if options['Number of elements in thorax'] < 1:
            options['Number of elements in thorax'] = 1

        if options['Shell thickness'] < 0:
            options['Shell thickness'] = -options['Shell thickness']
        elif options['Shell thickness'] < 0.0001:
            if options['Number of elements across shell'] >= 1:
                options['Shell thickness'] = 0.2
                dependentChanges = True

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
        elementsCountAlongAbdomen = options['Number of elements in abdomen']
        elementsCountAlongHead = options['Number of elements in head']
        elementsCountAlongNeck = options['Number of elements in neck']
        elementsCountAlongThorax = options['Number of elements in thorax']
        shellThickness = options['Shell thickness']
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
        neckGroup = AnnotationGroup(region, get_body_term("neck"))
        headGroup = AnnotationGroup(region, get_body_term("head"))
        topLeftGroup = AnnotationGroup(region, get_body_term("top left"))
        topRightGroup = AnnotationGroup(region, get_body_term("top right"))
        non_coreFirstLayerGroup = AnnotationGroup(region, get_body_term("non core first layer"))
        annotationGroups = [bodyGroup, coreGroup, non_coreGroup, abdomenGroup, thoraxGroup, neckGroup, headGroup,
                            topLeftGroup, topRightGroup]

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            shellThickness,
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
        e1oa = elementsCountAcrossMinor - 2*elementsCountAcrossShell - 2
        e2oa = elementsCountAcrossMajor - 2*elementsCountAcrossShell
        e1ob = e1oa*elementsCountAcrossShell
        e1oc = e1ob + e1oa + elementsCountAcrossShell
        e2oCore = (e1oa+2)*(e2oa)-4
        e2oShell = cylinder1._elementsCountAround*elementsCountAcrossShell
        e2o = e2oCore + e2oShell
        e1oy = e2o - e1oa*(elementsCountAcrossShell+1)
        for e3 in range(elementsCountAlong):
            for e2 in range(e2oa):
                for e1 in range(1, e1oa+3):
                    if e2 == 0:
                        if e1 <= e1oa:
                            elementIdentifier = e3 * e2o + e1ob + e1
                        else:
                            continue
                    elif e2 == e2oa - 1:
                        if e1 <= e1oa:
                            elementIdentifier = e3 * e2o + e1oy + e1
                        else:
                            continue
                    else:
                        elementIdentifier = e3 * e2o + elementsCountAcrossMinor * (e2 - 1) + e1 + e1oc

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

        topLeftMeshGroup = topLeftGroup.getMeshGroup(mesh)
        topRightMeshGroup = topRightGroup.getMeshGroup(mesh)
        for e3 in range(elementsCountAlong - elementsCountAlongHead):
            elementIdentifier = e3 * e2o + e2o//2 - elementsCountAcrossShell
            element = mesh.findElementByIdentifier(elementIdentifier)
            topLeftMeshGroup.addElement(element)
            element = mesh.findElementByIdentifier(elementIdentifier + elementsCountAcrossMinor)
            topRightMeshGroup.addElement(element)

        # create discontinuity in d3 on the core boundary
        e1a = elementsCountAcrossShell
        e2a = elementsCountAcrossShell
        e1z = elementsCountAcrossMinor - elementsCountAcrossShell - 1
        e2z = elementsCountAcrossMajor - elementsCountAcrossShell - 1
        non_coreFirstLayerMeshGroup = non_coreFirstLayerGroup.getMeshGroup(mesh)
        for e3 in range(elementsCountAlong):
            for e2 in range(e2a-1, e2z + 2):
                for e1 in range(elementsCountAcrossMinor):
                    if e2 == e2a - 1:
                        if (e1 > e1a) and (e1 < e1z):
                            elementIdentifier = e3 * e2o + e1oa * (elementsCountAcrossShell - 1) + e1 - elementsCountAcrossShell
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            non_coreFirstLayerMeshGroup.addElement(element)
                    elif e2 == e2z + 1:
                        if (e1 > e1a) and (e1 < e1z):
                            elementIdentifier = e3 * e2o + e2o - e1oa - e1oa * (elementsCountAcrossShell - 1) + e1 - elementsCountAcrossShell
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            non_coreFirstLayerMeshGroup.addElement(element)
                    elif (e2 > e2a) and (e2 < e2z):
                        if e1 == e1a - 1:
                            elementIdentifier = e3 * e2o + e1oa * (elementsCountAcrossShell+1) + (e2-(e2a+1)) * elementsCountAcrossMinor + elementsCountAcrossShell
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            non_coreFirstLayerMeshGroup.addElement(element)
                        elif e1 == e1z + 1:
                            elementIdentifier = e3 * e2o + e1oa * (elementsCountAcrossShell+1) + (e2-(e2a+1)) * elementsCountAcrossMinor + elementsCountAcrossMinor - elementsCountAcrossShell + 1
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            non_coreFirstLayerMeshGroup.addElement(element)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        elementtemplate = mesh.createElementtemplate()
        undefineNodetemplate = nodes.createNodetemplate()
        undefineNodetemplate.undefineField(coordinates)
        nodetemplate = nodes.createNodetemplate()
        fieldcache = fieldmodule.createFieldcache()
        with ChangeManager(fieldmodule):
            elementIter = non_coreFirstLayerMeshGroup.createElementiterator()
            element = elementIter.next()
            localNodeIndexes = [1, 2, 3, 4]
            valueLabel = Node.VALUE_LABEL_D_DS3
            while element.isValid():
                print('Element', element.getIdentifier())
                eft = element.getElementfieldtemplate(coordinates, -1)
                nodeIds = get_element_node_identifiers(element, eft)
                for localNodeIndex in localNodeIndexes:
                    node = element.getNode(eft, localNodeIndex)
                    print('    Node', node.getIdentifier())
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
                        # print('Results:', result0, result1, result2, result3, result4, result5)
                        # if result3 == 0:
                        #     loggerMessageCount = logger.getNumberOfMessages()
                        #     if loggerMessageCount > 0:
                        #         for i in range(1, loggerMessageCount + 1):
                        #             print(logger.getMessageTypeAtIndex(i), logger.getMessageTextAtIndex(i))
                        #         logger.removeAllMessages()
                        #         exit()
                remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)
                result1 = elementtemplate.defineField(coordinates, -1, eft)
                result2 = element.merge(elementtemplate)
                result3 = element.setNodesByIdentifier(eft, nodeIds)
                # print('Element merge result',result1, result2, result3)
                # if (result1 != RESULT_OK) or (result2 != RESULT_OK):
                #     loggerMessageCount = logger.getNumberOfMessages()
                #     if loggerMessageCount > 0:
                #         for i in range(1, loggerMessageCount + 1):
                #             print(logger.getMessageTypeAtIndex(i), logger.getMessageTextAtIndex(i))
                #         logger.removeAllMessages()
                #         exit()
                element = elementIter.next()
        # region.writeFile('body_groups_discontinuity11.exf')

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
        topLeftGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("top left"))
        topRightGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("top right"))

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

        mesh1d = fm.findMeshByDimension(1)
        spinalCordMeshGroup = spinalCordGroup.getMeshGroup(mesh1d)
        is_top_middle_faces = fm.createFieldAnd(topLeftGroup.getGroup(), topRightGroup.getGroup())
        is_spinal_cord = fm.createFieldAnd(is_top_middle_faces, is_core_boundary)
        spinalCordMeshGroup.addElementsConditional(is_spinal_cord)