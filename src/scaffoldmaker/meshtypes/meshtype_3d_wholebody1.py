"""
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
"""

from __future__ import division

import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, get_element_node_identifiers
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation import heart_terms, bladder_terms, lung_terms, stomach_terms, brainstem_terms
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.annotationgroup import getAnnotationGroupForTerm
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.annotation.nerve_terms import get_nerve_term
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
                    [[0.0, 0.0, 0.0], setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5),
                     [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]],
                    [setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis1, cylinder1Settings['Length']),
                     setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]]
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
            'Human Coarse',
            'Human Fine',
            'Rat Coarse',
            'Rat Fine'
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Human Coarse'
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Default']
        options = {}
        options['Base parameter set'] = parameterSetName
        options['Central path'] = copy.deepcopy(centralPathOption)
        options['Number of elements across major'] = 6
        options['Number of elements across minor'] = 6
        options['Number of elements across shell'] = 1
        options['Number of elements across transition'] = 1
        options['Number of elements in abdomen'] = 5
        options['Number of elements in thorax'] = 3
        options['Number of elements in neck'] = 1
        options['Number of elements in head'] = 2
        options['Shell thickness proportion'] = 0.2
        options['Discontinuity on the core boundary'] = True
        options['Lower half'] = False
        options['Use cross derivatives'] = False
        options['Refine'] = False
        options['Refine number of elements across major'] = 1
        options['Refine number of elements along'] = 1

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

        baseParameterSetName = options['Base parameter set']
        isHuman = 'Human' in baseParameterSetName
        isRat = 'Rat' in baseParameterSetName

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

        elementsCountAlong = elementsCountAlongAbdomen + elementsCountAlongThorax + elementsCountAlongNeck +\
                             elementsCountAlongHead

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

        # body coordinates
        bodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="body coordinates")
        tmp_region = region.createRegion()
        tmp_fieldmodule = tmp_region.getFieldmodule()
        tmp_body_coordinates = findOrCreateFieldCoordinates(tmp_fieldmodule, name="body coordinates")
        tmp_cylinder = CylinderMesh(tmp_fieldmodule, tmp_body_coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)
        sir = tmp_region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        tmp_region.write(sir)
        result, buffer = srm.getBuffer()
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemoryBuffer(buffer)
        region.read(sir)

        del srm
        del sir
        del tmp_body_coordinates
        del tmp_fieldmodule
        del tmp_region

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
                                        result0, x = coordinates.getNodeParameters(fieldcache, -1,
                                                                                   Node.VALUE_LABEL_VALUE, 1, 3)
                                        result0, d1 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                    Node.VALUE_LABEL_D_DS1, 1, 3)
                                        result0, d2 = coordinates.getNodeParameters(fieldcache, -1,
                                                                                    Node.VALUE_LABEL_D_DS2, 1, 3)
                                        result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1, 3)
                                        result1 = node.merge(undefineNodetemplate)
                                        result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 2)
                                        result3 = node.merge(nodetemplate)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                Node.VALUE_LABEL_VALUE, 1, x)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                Node.VALUE_LABEL_D_DS1, 1, d1)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1,
                                                                                Node.VALUE_LABEL_D_DS2, 1, d2)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                        result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)
                                result1 = elementtemplate.defineField(coordinates, -1, eft)
                                result2 = element.merge(elementtemplate)
                                result3 = element.setNodesByIdentifier(eft, nodeIds)
        else:
            fieldcache = fieldmodule.createFieldcache()

        # organ landmarks groups
        apexOfHeart = heart_terms.get_heart_term('apex of heart')
        leftAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('left atrium epicardium venous midpoint')
        rightAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('right atrium epicardium venous midpoint')
        apexOfUrinaryBladder = bladder_terms.get_bladder_term('apex of urinary bladder')
        leftUreterJunctionWithBladder = bladder_terms.get_bladder_term('left ureter junction with bladder')
        rightUreterJunctionWithBladder = bladder_terms.get_bladder_term('right ureter junction with bladder')
        urethraJunctionWithBladderDorsal = bladder_terms.get_bladder_term('urethra junction of dorsal bladder neck')
        urethraJunctionWithBladderVentral = bladder_terms.get_bladder_term('urethra junction of ventral bladder neck')
        gastroesophagalJunctionOnLesserCurvature = stomach_terms.get_stomach_term(
            'esophagogastric junction along the lesser curvature on serosa')
        limitingRidgeOnGreaterCurvature = stomach_terms.get_stomach_term(
            'limiting ridge at the greater curvature on serosa')
        pylorusOnGreaterCurvature = stomach_terms.get_stomach_term(
            'gastroduodenal junction along the greater curvature on serosa')
        junctionBetweenFundusAndBodyOnGreaterCurvature = stomach_terms.get_stomach_term(
            "fundus-body junction along the greater curvature on serosa")
        apexOfLeftLung = lung_terms.get_lung_term('apex of left lung')
        ventralBaseOfLeftLung = lung_terms.get_lung_term('ventral base of left lung')
        dorsalBaseOfLeftLung = lung_terms.get_lung_term('dorsal base of left lung')
        apexOfRightLung = lung_terms.get_lung_term('apex of right lung')
        ventralBaseOfRightLung = lung_terms.get_lung_term('ventral base of right lung')
        dorsalBaseOfRightLung = lung_terms.get_lung_term('dorsal base of right lung')
        laterodorsalTipOfMiddleLobeOfRightLung = lung_terms.get_lung_term(
            'laterodorsal tip of middle lobe of right lung')
        apexOfRightLungAccessoryLobe = lung_terms.get_lung_term('apex of right lung accessory lobe')
        ventralBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('ventral base of right lung accessory lobe')
        dorsalBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('dorsal base of right lung accessory lobe')
        medialBaseOfLeftLung = lung_terms.get_lung_term("medial base of left lung")
        medialBaseOfRightLung = lung_terms.get_lung_term("medial base of right lung")
        brainstemDorsalMidlineCaudalPoint = brainstem_terms.get_brainstem_term(
            'brainstem dorsal midline caudal point')
        brainstemDorsalMidlineCranialPoint = brainstem_terms.get_brainstem_term(
            'brainstem dorsal midline cranial point')
        brainstemVentralMidlineCaudalPoint = brainstem_terms.get_brainstem_term(
            'brainstem ventral midline caudal point')
        brainstemVentralMidlineCranialPoint = brainstem_terms.get_brainstem_term(
            'brainstem ventral midline cranial point')

        # marker coordinates. In future we want to have only one table for all species.
        if isRat:
            bodyMarkerPoints = [
                {"group": ("left hip joint", ''), "x": [0.367, 0.266, 0.477]},
                {"group": ("right hip joint", ''), "x": [-0.367, 0.266, 0.477]},
                {"group": ("left shoulder joint", ''), "x": [0.456, -0.071, 2.705]},
                {"group": ("right shoulder joint", ''), "x": [-0.456, -0.071, 2.705]},
                {"group": ("along left femur", ''), "x": [0.456, 0.07, 0.633]},
                {"group": ("along right femur", ''), "x": [-0.456, 0.07, 0.633]},
                {"group": ("along left humerus", ''), "x": [0.423, -0.173, 2.545]},
                {"group": ("along right humerus", ''), "x": [-0.423, -0.173, 2.545]},
                {"group": apexOfUrinaryBladder, "x": [-0.124, -0.383, 0.434]},
                {"group": leftUreterJunctionWithBladder, "x": [-0.111, -0.172, 0.354]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.03, -0.196, 0.363]},
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.03, -0.26, 0.209]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.037, -0.304, 0.203]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [-0.032, 0.418, 2.713]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [-0.017, 0.203, 2.941]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [-0.028, 0.388, 2.72]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [-0.019, 0.167, 2.95]},
                {"group": apexOfHeart, "x": [0.096, -0.128, 1.601]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.127, -0.083, 2.079]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [0.039, -0.082, 2.075]},
                {"group": apexOfLeftLung, "x": [0.172, -0.175, 2.337]},
                {"group": ventralBaseOfLeftLung, "x": [0.274, -0.285, 1.602]},
                {"group": dorsalBaseOfLeftLung, "x": [0.037, 0.31, 1.649]},
                {"group": apexOfRightLung, "x": [-0.086, -0.096, 2.311]},
                {"group": ventralBaseOfRightLung, "x": [0.14, -0.357, 1.662]},
                {"group": dorsalBaseOfRightLung, "x": [-0.054, 0.304, 1.667]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.258, -0.173, 2.013]},
                {"group": apexOfRightLungAccessoryLobe, "x": [0.041, -0.063, 1.965]},
                {"group": ventralBaseOfRightLungAccessoryLobe, "x": [0.143, -0.356, 1.66]},
                {"group": dorsalBaseOfRightLungAccessoryLobe, "x": [0.121, -0.067, 1.627]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [0.12, 0.009, 1.446]},
                {"group": limitingRidgeOnGreaterCurvature, "x": [0.318, 0.097, 1.406]},
                {"group": pylorusOnGreaterCurvature, "x": [0.08, -0.111, 1.443]},
            ]
        elif isHuman:
            bodyMarkerPoints = [
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.0071, -0.2439, 0.1798]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.007, -0.2528, 0.1732]},
                {"group": leftUreterJunctionWithBladder, "x": [0.1074, 0.045, 0.1728]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.1058, 0.0533, 0.1701]},
                {"group": apexOfUrinaryBladder, "x": [0.005, 0.1286, 0.1264]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [0.0068, 0.427, 2.7389]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [0.008, -0.0231, 3.0778]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [0.0054, 0.3041, 2.7374]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [0.0025, -0.2308, 3.091]},
                {"group": apexOfHeart, "x": [0.1373, -0.1855, 1.421]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.0024, 0.1452, 1.8022]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [-0.0464, 0.0373, 1.7491]},
                {"group": apexOfLeftLung, "x": [0.0655, -0.0873, 2.3564]},
                {"group": apexOfRightLung, "x": [-0.088, -0.0363, 2.3518]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.2838, -0.0933, 1.9962]},
                {"group": ventralBaseOfLeftLung, "x": [0.219, -0.2866, 1.4602]},
                {"group": medialBaseOfLeftLung, "x": [0.0426, -0.0201, 1.4109]},
                {"group": ventralBaseOfRightLung, "x": [-0.2302, -0.2356, 1.3926]},
                {"group": medialBaseOfRightLung, "x": [-0.0363, 0.0589, 1.3984]},
                {"group": dorsalBaseOfLeftLung, "x": [0.1544, 0.2603, 1.3691]},
                {"group": dorsalBaseOfRightLung, "x": [0.0369, -0.2524, 0.912]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [-0.0062, -0.3259, 0.8586]},
                {"group": pylorusOnGreaterCurvature, "x": [-0.0761, -0.3189, 0.8663]},
                {"group": junctionBetweenFundusAndBodyOnGreaterCurvature, "x": [0.1884, -0.1839, 0.9639]},
            ]

        nodeIdentifier = cylinder1._endNodeIdentifier
        for bodyMarkerPoint in bodyMarkerPoints:
            annotationGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, bodyMarkerPoint["group"], isMarker=True)
            annotationGroup.createMarkerNode(nodeIdentifier, bodyCoordinates, bodyMarkerPoint["x"])
            nodeIdentifier += 1

        # Add annotations for the nerves fiducials
        markerTermNameNerveCoordinatesMap = {
            "C3_C5_spinal_n-phrenic_n": [-0.010138285213950535, 0.046746867202612036, 2.5973140511178574],
            "C3_dorsal_root_end": [0.04009154377867926, 0.42235778666136586, 2.679594014258995],
            "C3_spinal": [0.008354959186164506, 0.4020152910718021, 2.6617877641656733],
            "C3_ventral_root_end": [0.016423978370946554, 0.11322229502301429, 2.680340093002615],
            "C4_dorsal_root_end": [0.030566562004732283, 0.42557416070641385, 2.6345915939077065],
            "C4_spinal": [0.0086952058843628, 0.40299655708679416, 2.62113807182829],
            "C4_ventral_root_end": [0.01524840230383169, 0.11588302230464842, 2.635730500816299],
            "C5_dorsal_root_end": [0.03946112400871407, 0.4237861972345228, 2.5858033315964035],
            "C5_spinal": [0.014030332821629183, 0.4033186198288035, 2.570719946575644],
            "C5_ventral_root_end": [0.02843172889917061, 0.15556766343560036, 2.5870555102000035],
            "C6_dorsal_root_end": [0.0639892090705618, 0.4250139733240407, 2.5372114716214944],
            "C6_spinal": [-0.0020270150401892123, 0.40674309718172785, 2.525907820987904],
            "C7_dorsal_root_end": [0.0506051479493571, 0.42228546556253177, 2.4877814872260813],
            "C7_spinal": [0.004165487609905456, 0.4056608573753266, 2.481348272753478],
            "C8_dorsal_root_end": [0.02486514336892628, 0.4221816835541224, 2.4434641062411164],
            "C8_spinal": [0.0018577837929680627, 0.4050022805158228, 2.4395049908395947],
            "L1_dorsal_root_end": [0.025939970500889854, 0.4346763676202703, 1.0936604094280287],
            "L1_spinal": [-0.004384849913309982, 0.0903448478093626, 0.8341965328754388],
            "L1_ventral_root_end": [0.015875086446471454, 0.027504564692553815, 0.8353901771234414],
            "L2_dorsal_root_end": [0.028209865626316166, 0.43850427209263887, 0.924354735061185],
            "L2_spinal": [0.0007566375156007152, 0.4026926089912893, 0.8897290010308939],
            "L2_ventral_root_end": [0.018349629626238997, 0.023740947944245732, 0.7509913923857144],
            "L5_spinal": [-0.00024194527785240572, 0.40124171699896166, 0.5849167892990343],
            "L5_ventral_root_end": [0.027316100288407058, 0.2132181925709541, 0.522666125296055],
            "S1_dorsal_root_end": [0.03342074232118696, 0.42813337983082744, 0.510887411811316],
            "S1_spinal": [0.004476333170439515, 0.39790416286368574, 0.4924862211888632],
            "S1_ventral_root_end": [0.02998463713863766, 0.2721271778385812, 0.45788220894590553],
            "S43_C3": [0.03352929800850721, 0.14355138136811432, 2.6742758346881796],
            "S43_C4": [0.011594668400702182, 0.15471012667668957, 2.629138587552765],
            "S43_C5": [0.03258446072428602, 0.2207225149176348, 2.579755605097222],
            "S43_L1": [0.015223031785910151, 0.013366154052073422, 0.8289118942738474],
            "S43_L2": [0.01905455023120646, 0.044393943797411364, 0.7453945182875407],
            "S43_L5": [0.026368325254780067, 0.2384234029262272, 0.5239452458260677],
            "S43_S1": [0.02910187138647435, 0.30032992348851956, 0.45398653624920077],
            "S43_T1": [0.016525490747565513, 0.3116377417799617, 2.37385649603387],
            "S43_T10": [0.018939336580186632, 0.3513047023873593, 1.394070710443572],
            "S43_T11": [0.0218886915905585, 0.36381223761961173, 1.2979849083213453],
            "S43_T12": [0.011227647627429544, 0.14774776016617366, 0.972491706596798],
            "S43_T2": [0.01978885496937896, 0.33102468432269544, 2.2885178222413898],
            "S43_T3": [0.017703748271177946, 0.3589986564649193, 2.19544647500101],
            "S43_T4": [0.013342441492638248, 0.36681659323088667, 2.0833107117099066],
            "S43_T5": [0.01112620365204089, 0.36490791984142473, 1.9529163980372526],
            "S43_T6": [0.01416227521186865, 0.34151734048835847, 1.8379090940298832],
            "S43_T7": [0.01874273808681078, 0.3413596036485539, 1.7237832332696168],
            "S43_T8": [0.017778981148898813, 0.33957715738203503, 1.613070626300565],
            "S43_T9": [0.018291914360204317, 0.3420299037123524, 1.516234146158932],
            "S44_C3": [0.04585957304677293, 0.42153297056122085, 2.6764455912110874],
            "S44_C4": [0.031051106702403762, 0.4225486307381243, 2.630928990129954],
            "S44_C5": [0.03621132036706466, 0.4223711686793854, 2.5823923786978287],
            "S44_C6": [0.06631875543140325, 0.4214287040198945, 2.535353490518112],
            "S44_C7": [0.055150860229967354, 0.4185800600742488, 2.4859398542155984],
            "S44_C8": [0.02565027649697677, 0.4195642200102317, 2.4414866926400403],
            "S44_L1": [0.026661616545576448, 0.4328047508536241, 1.0850500347731722],
            "S44_L2": [0.030063978186969532, 0.4357279993952721, 0.9166354891216831],
            "S44_S1": [0.030135518094151586, 0.42256640115373906, 0.5038452436452542],
            "S44_T1": [0.022810405713903647, 0.4165206115175227, 2.38410684318936],
            "S44_T2": [0.01923695475768973, 0.4136717951646068, 2.2991646415606546],
            "S44_T3": [0.016944156961161673, 0.4146587729285394, 2.197661524197905],
            "S44_T4": [0.011491484961029166, 0.41918983463066073, 2.0918033784147827],
            "S44_T5": [0.00982132010124722, 0.423818991934699, 1.966284982394008],
            "S44_T6": [0.012483990836193622, 0.4301663171202181, 1.8579200975349006],
            "S45_C3": [0.044008339626270467, 0.4203782073436899, 2.673417168832901],
            "S45_C4": [0.0351021854659379, 0.41932993185300177, 2.627869222810391],
            "S45_C5": [0.03967329386724587, 0.41889454959459776, 2.5805212379306366],
            "S45_C6": [0.0670176797225276, 0.41835530509874536, 2.533130255037592],
            "S45_C7": [0.056433739971485056, 0.41529988763787407, 2.4850656030028566],
            "S45_C8": [0.027426776536023083, 0.4165741670797335, 2.440926707838184],
            "S45_L1": [0.0195959788760001, 0.20243330625518594, 0.8910682997995971],
            "S45_L2": [0.029243057379758332, 0.4326098116615097, 0.9085126370047192],
            "S45_S1": [0.0339963642334792, 0.4198531912166759, 0.503001130519617],
            "S45_T1": [0.019533363273623473, 0.4130270353258266, 2.379827295749163],
            "S45_T2": [0.02015832139065615, 0.4101052992604936, 2.2974592386270745],
            "S45_T3": [0.02215865237454608, 0.41118122484317304, 2.193568901235038],
            "S45_T4": [0.011681231290381483, 0.41578813697456546, 2.088311570862359],
            "S45_T5": [0.009251588701565915, 0.4202611008185268, 1.961031337445952],
            "S45_T6": [0.013092077515424303, 0.4251255573471926, 1.849390184197928],
            "S46_L1": [0.01815334803904297, 0.07573521840446609, 0.8244313022530534],
            "S46_L2": [0.022278003364078675, 0.10935694267580652, 0.7473338763090714],
            "S46_L3": [0.026586711691535796, 0.15265132037526932, 0.6894014078021441],
            "S46_L4": [0.030222640291741938, 0.2630649544375915, 0.6100367654335477],
            "S46_L5": [0.02672167368126046, 0.3092528798487019, 0.5112953579811701],
            "S46_S1": [0.032379260761676056, 0.3611231962473034, 0.4584490808540732],
            "S46_S2": [0.031053648110984974, 0.3660894429396693, 0.4033995494050888],
            "S46_T10": [0.020962255771327057, 0.39298210896855346, 1.3840025019881141],
            "S46_T11": [0.02102757935496368, 0.39369784349084386, 1.2848862676444577],
            "S46_T12": [0.011604851320986686, 0.1801313550000583, 0.9633231721280667],
            "S46_T2": [-0.003984587649489822, 0.40008608057457534, 2.2846747099033955],
            "S46_T3": [0.01395615786205882, 0.3962229870859504, 2.18333075854326],
            "S46_T4": [0.012029480870225505, 0.39833102361052314, 2.0743106056625473],
            "S46_T5": [0.011258222786843316, 0.3958160344741454, 1.9446986116252072],
            "S46_T6": [0.015453540921386804, 0.3953961485915715, 1.8297086461487904],
            "S46_T7": [0.016050633298270164, 0.39594111662889586, 1.709351022218289],
            "S46_T8": [0.022157310708381624, 0.3894505080747455, 1.6007911046159984],
            "S46_T9": [0.02173137595756408, 0.39183836215695067, 1.5003818350360076],
            "S48_L1": [0.018297801075471062, 0.06862228823503885, 0.8283256877690376],
            "S48_L2": [0.020201264131216666, 0.09264856006220012, 0.7470710759054864],
            "S48_T10": [0.022854569151858482, 0.3878148236648599, 1.3881003827164913],
            "S48_T11": [0.02137062406568845, 0.388631827322188, 1.2892502028612567],
            "S48_T12": [0.011475976089809908, 0.17818972506922556, 0.9685262026861576],
            "S48_T2": [-0.006151363002076865, 0.3910969290160922, 2.2892905995494717],
            "S48_T3": [0.015729238757681226, 0.3918482839913054, 2.1875678560303413],
            "S48_T4": [0.012722340347631993, 0.3952597530095745, 2.0783533908800598],
            "S48_T5": [0.010710630988212979, 0.3913272689841617, 1.9494137154183988],
            "S48_T6": [0.014661728181654085, 0.3898482169951014, 1.8319137149346336],
            "S48_T7": [0.017066251093088947, 0.3828533589422067, 1.711036218793946],
            "S48_T8": [0.020855867579673237, 0.3806356794143022, 1.6058805901382305],
            "S48_T9": [0.02207739912819795, 0.38783611890159564, 1.5049939605563534],
            "S49_L1": [0.019316494040743662, 0.08072199564346509, 0.8219737740847872],
            "S49_L2": [0.02087110088942092, 0.11823167223905497, 0.7444637157327187],
            "S49_T10": [0.02152608500033666, 0.39740987038914233, 1.3760893295943988],
            "S49_T11": [0.01989732740676114, 0.39707297364047234, 1.277167951026701],
            "S49_T12": [0.012338437853365856, 0.17219419804263555, 0.9524497292323838],
            "S49_T2": [-0.003884302304418874, 0.4004370914321539, 2.2778593053180693],
            "S49_T3": [0.012884393694254792, 0.3987738692340437, 2.1788581812815466],
            "S49_T4": [0.012352033383656183, 0.4005515407613809, 2.067045917709552],
            "S49_T5": [0.01144042680543557, 0.3987204786473204, 1.9386949200446075],
            "S49_T6": [0.017219730869396172, 0.39885582988532853, 1.8221971027934027],
            "S49_T7": [0.017678788295255266, 0.3983746129320296, 1.6989962812064117],
            "S49_T8": [0.02120578846059063, 0.3947404147252563, 1.5945322717344004],
            "S49_T9": [0.02134478137690721, 0.39367352862071825, 1.4947378197880996],
            "S50_C3_B": [0.04379813033725092, 0.41076373332281335, 2.668611105303296],
            "S50_C3_T": [0.03383507708453722, 0.27271773528595644, 2.6670860396516787],
            "S50_C4_B": [0.03264170274508996, 0.4145596485839877, 2.6249267131800242],
            "S50_C4_T": [0.012905492582429895, 0.25279094434325095, 2.6249856240620106],
            "S50_C5_B": [0.042746746054760976, 0.4148253017456548, 2.575581391124371],
            "S50_C5_T": [0.03949335600258939, 0.3219190611364926, 2.577766388707105],
            "S50_C6_B": [0.0600995274891082, 0.41543847234582765, 2.531483374270639],
            "S50_C7_B": [0.051478101084571924, 0.4121950188168679, 2.4842817170398273],
            "S50_C8_B": [0.026777609528986356, 0.4128426058308994, 2.4393064807005196],
            "S50_L1": [0.016594745403087398, 0.05236789435985481, 0.827375283729268],
            "S50_L1_B": [0.0197676410961049, 0.1799744883009029, 0.8759369734675617],
            "S50_L2": [0.02077468142050592, 0.07531672838139492, 0.7455383147113345],
            "S50_L2_B": [0.031865781651167845, 0.4286225770980082, 0.900762066184893],
            "S50_L5_B": [0.02639456283847429, 0.42819919771006476, 0.588085380035991],
            "S50_L5_T": [0.033143046498861815, 0.31172189479245965, 0.5467054843360633],
            "S50_S1_B": [0.03260681839900097, 0.41152528051542203, 0.4979911198628674],
            "S50_S1_T": [0.0329411164567827, 0.3447554781847321, 0.47043059487079786],
            "S50_T1": [0.020549131152726798, 0.369408069016317, 2.3767306863905753],
            "S50_T10": [0.02118062939089793, 0.3774923397394473, 1.3916775752106179],
            "S50_T11": [0.02187439989265712, 0.3801974399890965, 1.292302395954607],
            "S50_T12": [0.010596347491229808, 0.1742322734377156, 0.9725553088183632],
            "S50_T1_B": [0.019523686745962398, 0.4096040874068864, 2.3771822967364282],
            "S50_T2": [0.019771354028778238, 0.37694271134775154, 2.292475604481717],
            "S50_T2_B": [0.020303273482013266, 0.4072399515199125, 2.2962497509933972],
            "S50_T3": [0.0168501747789926, 0.3831257489648913, 2.1921779675631736],
            "S50_T3_B": [0.016744607175711662, 0.4094606736270623, 2.1936273092794987],
            "S50_T4": [0.013033634342149563, 0.3887360803675961, 2.0843733695189477],
            "S50_T4_B": [0.012979447082041915, 0.41226197200460846, 2.084493453493506],
            "S50_T5": [0.010712930707445107, 0.38380755308139747, 1.9539689066843213],
            "S50_T5_B": [0.009445485944845212, 0.41641270148333503, 1.9616755993515165],
            "S50_T6": [0.01537279654530296, 0.37990737821472914, 1.832718097905838],
            "S50_T6_B": [0.013440830560440944, 0.41929918866927146, 1.8432507545520709],
            "S50_T7": [0.01829009830287092, 0.3713300850499354, 1.721714848970642],
            "S50_T8": [0.020070384488674674, 0.3682775134117164, 1.6093274416738839],
            "S50_T9": [0.02144938557745213, 0.37157561341761214, 1.5122322794621406],
            "T10_spinal": [0.0004180048518173882, 0.4049514526144774, 1.390390044135262],
            "T10_ventral_root_end": [0.018706719553620195, 0.3331390017584663, 1.4018110449159857],
            "T11_spinal": [-0.004763546749091652, 0.40297890392061186, 1.2904750889067718],
            "T11_ventral_root_end": [0.021133991332534817, 0.3493495384803953, 1.3029401825425748],
            "T12_spinal": [-0.0025410701345697445, 0.4027367444389069, 1.1930700307209394],
            "T12_ventral_root_end": [0.015518380191592551, 0.3273063629701729, 1.1696194174428083],
            "T1_dorsal_root_end": [0.023822778951030205, 0.4192337130575748, 2.3871626726521145],
            "T1_spinal": [-0.00013920391178511604, 0.39657531703162663, 2.372882207920326],
            "T1_ventral_root_end": [0.018072779684073517, 0.2760413861979546, 2.381735512877754],
            "T2_dorsal_root_end": [0.020446328038535848, 0.4165673052554153, 2.301307102940742],
            "T2_spinal": [-0.0023400318832315545, 0.40194288267514494, 2.2967894447067168],
            "T2_ventral_root_end": [0.019092955663044343, 0.3046227084789411, 2.2971118557076626],
            "T3_dorsal_root_end": [0.016750099993353928, 0.4170030101652892, 2.2005458154830992],
            "T3_spinal": [0.0010470984266154225, 0.40260309419544765, 2.1911971061933624],
            "T3_ventral_root_end": [0.016749130396468023, 0.3384633819232629, 2.202318443810073],
            "T4_dorsal_root_end": [0.011741946056976255, 0.42294892926796307, 2.093604610204217],
            "T4_spinal": [-0.0008359681186187092, 0.4042247121462797, 2.0866228994859317],
            "T4_ventral_root_end": [0.012314011701061035, 0.35294201491359983, 2.0875456882790178],
            "T5_dorsal_root_end": [0.009369167664632699, 0.4260826070725257, 1.973095795034117],
            "T5_spinal": [-0.004580020028552358, 0.40490694160862456, 1.9589799761929745],
            "T5_ventral_root_end": [0.010652710344158848, 0.3464740499787812, 1.9600388949763532],
            "T6_dorsal_root_end": [0.011617869904320515, 0.4338612664350838, 1.8651631346679873],
            "T6_spinal": [-0.0014007673443994602, 0.4048554158357282, 1.8402695744326387],
            "T6_ventral_root_end": [0.013014639196911941, 0.32730338338297066, 1.8452243466884617],
            "T7_spinal": [0.000720759783012804, 0.4040180474647821, 1.7272997162170596],
            "T7_ventral_root_end": [0.017626968275536333, 0.32127725881591407, 1.7287897365114748],
            "T8_spinal": [-0.00030260611202322427, 0.4035331395779552, 1.6115449231958736],
            "T8_ventral_root_end": [0.01713935278215176, 0.3317252039935669, 1.620478039936626],
            "T9_spinal": [0.0020663997878669514, 0.40335813057109193, 1.5131447787035754],
            "T9_ventral_root_end": [0.01783882435252728, 0.32215925989901345, 1.5217821546779977],
            "ardell_1_branching_point": [0.1369243808878739, -0.014345676071479081, 1.7916987591437523],
            "ardell_1_start-1": [-2.3667508977702825e-06, 0.10449866786907987, 2.244578349073853],
            "ardell_1_start-2": [-0.004792910383749542, 0.16640878699549128, 2.140309803267339],
            "ardell_1_start-3": [0.022213571277388076, 0.19982430468925608, 1.9832306487963618],
            "ardell_1_start-4": [0.05263893224095571, 0.21096438092969103, 1.8817341813633324],
            "ardell_1_start-5": [0.010064702064377396, 0.22605626845353732, 1.7696418333325126],
            "ardell_1_start-6": [-0.01532522558245051, 0.2115534604157658, 1.6860891041195112],
            "ardell_1_start-7": [-0.021107608099796362, 0.19398632193300658, 1.6054710499221072],
            "ardell_1_start-8": [-0.00720084401991996, 0.18284889058676512, 1.5471281517594921],
            "ardell_1_start-9": [0.023623459547606696, 0.16669646291368984, 1.469306008161244],
            "ardell_2_branching_point": [0.06631618371791677, 0.07972506176270236, 1.8458736999218497],
            "ardell_2_start": [-0.049381269493019755, 0.04298081680877442, 2.113908504134355],
            "ardell_3_branching_point": [0.09036267655667948, 0.06572121909001478, 1.7726746289072972],
            "ardell_4_branching_point": [0.03588790684051415, 0.00010530996209381446, 1.8268258493096292],
            "ardell_5_branching_point": [0.15500608700246332, -0.07508956251317368, 1.728057761937825],
            "ardell_5_start-1": [0.03318075428152776, 0.15266961140506874, 1.7155592880066306],
            "ardell_5_start-2": [0.02817268349999068, 0.08322905614486505, 1.5966281469999084],
            "ardell_5_start-3": [0.04623157769218973, 0.0711630874603279, 1.5753716810125757],
            "ardell_5_start-4": [0.07930815179273962, 0.037327520208753955, 1.5072468193682795],
            "ardell_5_start-5": [0.09833866857039474, 0.15143719818494145, 1.6017204570500967],
            "ardell_5_start-6": [0.16620455994167818, 0.118757624664246, 1.5480602453915788],
            "ardell_6_branching_point": [0.1224768684378978, 0.04432130594639153, 1.7068753225051707],
            "ardell_6_start-1": [0.020604360234527684, -0.20334300194785968, 2.4792636645743573],
            "ardell_6_start-2": [0.027486994236058732, -0.05477198714923579, 2.2354806868475836],
            "ardell_6_start-3": [0.013645339775015527, 0.11488224948470467, 2.013748045055838],
            "ardell_6_start-4": [0.025486682464064013, 0.162901876247143, 1.8877493303637087],
            "ardell_6_start-5": [0.023120428244867854, 0.18676565989741595, 1.8309620849113861],
            "ardell_6_start-6": [0.060516873814746415, 0.2125953394675951, 1.7594622891372789],
            "ardell_6_start-7": [0.09242145815131268, 0.21428090415931048, 1.6750682080785892],
            "ardell_6_start-8": [0.02664899054007593, 0.15827435297175363, 1.607783831949774],
            "ardell_6_start-9": [0.05233908062408706, 0.12863830591332986, 1.5030123296442794],
            "bladder_n-bladder": [-0.010682800370210993, -0.26780302510284404, 0.18436987778036062],
            "bolser_14-1": [-0.030191359837631868, 0.13803182983187548, 2.2866794057980924],
            "bolser_1_branching_point": [0.04393354616171909, 0.07052583787592064, 1.8633896395517144],
            "bolser_6_branching_point": [-0.036429213528043104, 0.058714414593065964, 2.056515250768737],
            "bolser_8_branching_point": [0.21701818859519348, 0.0949744897548108, 2.6603315292853185],
            "bolser_9-2": [-0.021603291975977962, 0.23271910299590398, 1.9191895354865052],
            "bolser_9_branching_point": [0.16539586396154224, -0.42688868533737234, 2.8321714521838643],
            "brain_12-1": [0.025820271808667028, 0.4082595084249503, 2.7943919120024647],
            "brain_33-1": [0.06271257092632289, 0.40010739787982447, 2.799475240307364],
            "brain_34-1": [0.05010831497260702, 0.40358054460656634, 2.7973191326425493],
            "brain_44-1": [-0.044728774383474966, 0.4016674212303391, 2.7824382133736556],
            "cardio_4-1": [0.0038770226360206376, 0.25767914299230826, 2.602417781185171],
            "connection_1-1": [0.0006564560395649186, 0.16375013232616326, 2.0074044740088555],
            "digestive_9-1": [0.010669500132744791, 0.16891559155709826, 1.3365923385869418],
            "external_laryngeal_n_branching_point": [0.0018535909234143062, -0.2352791525525333, 2.6144321267063364],
            "ganglion_1-1": [0.3020475906782037, 0.05386859227091203, 2.785763125387399],
            "ganglion_12-1": [-0.049874255219609094, 0.44358590234270673, 2.7527840034694493],
            "ganglion_2-1": [0.2801674720894303, -0.07392493888403924, 2.762244666472658],
            "ganglion_3-1": [-0.1979577569337593, 0.16756019227087676, 2.701122752982803],
            "ganglion_4-1": [-0.015495802306862448, 0.03742302038104417, 0.38155169919575715],
            "ganglion_5-1": [-0.0025508059609870015, -0.007831127462770288, 0.7581835314790677],
            "keast_1-1": [-0.00881671803277357, -0.10254418047101983, 0.27162845873260666],
            "keast_2-1": [-0.037305588962563356, 0.16154733274538596, 0.2753955993836467],
            "keast_4-1": [0.011321346036656355, 0.17466781302740916, 0.38411452624839415],
            "keast_5-1": [-0.00887712564221177, 0.0658759394985676, 0.5671446078032861],
            "keast_7-1": [-0.01702264168332862, 0.008984090750420613, 0.7641858476274449],
            "label_1-1": [0.24233213453853622, -0.04804262865546726, 2.8196900030814707],
            "label_10-1": [0.022394051014532946, -0.39221724730545904, 2.553522826546636],
            "label_11-1": [0.004500326769316632, -0.36988507641040924, 2.544193524974597],
            "label_12-1": [0.0020658231610839006, -0.34963628645858386, 2.5458664170933822],
            "label_13-1": [0.09775766365636465, -0.39112140797026373, 2.508089173593395],
            "label_14-1": [0.05074621175834829, -0.38093308136865733, 2.5323943856224678],
            "label_15-1": [0.11501767385582265, -0.40029981358564515, 2.47669758320286],
            "label_16-1": [0.12213294103766527, -0.4015620926588465, 2.4715304332076653],
            "label_17-1": [0.13429173986739315, -0.41956600214864165, 2.4844122276747966],
            "label_18-1": [0.13118107193037185, -0.4222180631879017, 2.486686567467264],
            "label_19-1": [0.14382915805118782, -0.41752239965803817, 2.47119868535046],
            "label_2-1": [-0.48346617263352126, -0.025510964732225982, 3.0004268932924205],
            "label_20-1": [0.13883759809750487, -0.42441862943850733, 2.48062607956546],
            "label_21-1": [0.1356183086782928, -0.42084061267668416, 2.476280267916309],
            "label_22-1": [0.14274212856165702, -0.41974267331361736, 2.4821904384494875],
            "label_3-1": [-0.4621861848712916, -0.03713557697482017, 2.996330897308272],
            "label_4-1": [0.45862685592014296, 0.10832561962746678, 2.9968409221178],
            "label_5-1": [0.04749952034604189, -0.4731496701036091, 2.7354285380908574],
            "label_6-1": [-0.07029502684828838, -0.46761336074289117, 2.7428948073608295],
            "label_7-1": [-0.0035482676439523393, -0.451027107463609, 2.670617571761484],
            "label_8-1": [0.06097579573751892, -0.44747025425450593, 2.666349328168109],
            "label_9-1": [0.09335895472040188, -0.3973609184550194, 2.502480936788414],
            "lumbar_splanchnic_n_start": [0.06214928110012197, 0.047795739727471975, 0.7706944263492905],
            "pelvic_splanchnic_n_start": [0.05714659450802255, 0.2724987232954604, 0.4180897394941961],
            "phrenic_n_branching_point": [0.014868625851324674, -0.2519158820205103, 0.9618865499199392],
            "plexus_1-1": [-0.0972171011004839, -0.08591492928330748, 2.6637216043658536],
            "plexus_2-1": [0.06245736813616862, -0.15086108592493655, 2.7107130975329503],
            "plexus_3-1": [0.007353573769681592, 0.0299969365901009, 2.182311533186655],
            "plexus_4-1": [0.007898763945867622, 0.08570871711106318, 2.0448666409377836],
            "plexus_5-1": [0.001366998014844518, 0.11090846662459104, 1.5185139782549681],
            "plexus_6-1": [0.007413189041136012, 0.2996132541214259, 1.25953666413364],
            "plexus_7-1": [-0.007590818809016026, 0.2264225727403843, 2.672684068494318],
            "point_1": [0.006688377104042991, 0.36130614946820705, 2.7692703495009576],
            "point_10": [-0.11489735089858619, 0.39310376770945, 2.7559044885420714],
            "point_11": [0.36452067845039055, -0.1568622972829942, 2.747737663625703],
            "point_12": [0.07121085922231549, -0.3117292589864084, 2.74017039920701],
            "point_13": [0.07329181585594936, -0.4116912341115187, 2.738855750258359],
            "point_14": [0.3685486356450809, -0.18313252602318847, 2.7150843778022886],
            "point_15": [0.21932489686677853, 0.006102080178131754, 2.69484214699987],
            "point_16": [0.36300046622076065, -0.19455385695238217, 2.68695514695242],
            "point_17": [0.30086875262134366, -0.23827354832677328, 2.6758404237429265],
            "point_18": [0.21374122184525918, -0.26369552583726136, 2.667880910742195],
            "point_19": [0.1353412161732525, -0.3380656918358264, 2.6647053250832364],
            "point_2": [0.29739822540454786, 0.061709356440369206, 2.7884411795325077],
            "point_20": [0.05067585034789162, -0.4056171263376282, 2.657691665695792],
            "point_21": [0.07222410740745247, -0.35846792443982944, 2.6345437739823585],
            "point_22": [0.3523824589385048, -0.20269092729744742, 2.659212622546793],
            "point_23": [0.34559202075013856, -0.23744354124282926, 2.6229723985115423],
            "point_24": [0.27535423033568923, -0.22490727360358811, 2.57106784894404],
            "point_25": [0.24275313161888792, -0.21227528809519675, 2.5389170229704394],
            "point_26": [0.26194940315852294, -0.22878849290520317, 2.5077463472424513],
            "point_27": [0.14163316307099644, -0.2454225427743763, 2.508727106601158],
            "point_28": [0.23878083238825953, -0.31199954437195726, 2.5041708811967704],
            "point_29": [0.23228574016095455, -0.3473953047512029, 2.4983348791928774],
            "point_3": [0.25821940647637814, 0.03806082450853832, 2.8046123750973733],
            "point_30": [0.08426315548064099, -0.18366350193835393, 2.18875861255995],
            "point_31": [0.07329625165271166, -0.141374069772622, 2.1235132678885944],
            "point_32": [0.015448450304965746, -0.03504268121155852, 1.9330136062796626],
            "point_33": [0.013405251310978262, 0.09210711286066171, 1.6208546636649368],
            "point_34": [0.010381322135109313, 0.031837478123001625, 1.4033596867558122],
            "point_35": [-0.01420988442545449, 0.14978880647209547, 1.2960058605681164],
            "point_36": [0.11899966493490469, -0.06911533230264504, 0.9119789179345877],
            "point_37": [0.11331594667626851, -0.3176961116041969, 2.6001074027150963],
            "point_38": [0.23412117328909252, -0.3483128685462542, 2.5025418766285346],
            "point_4": [0.30522596313857103, 0.07002006324936161, 2.7827050800615436],
            "point_40": [0.019680145651425986, 0.22891519923578704, 1.78667042625321],
            "point_41": [0.4121777811037062, 0.05411916938289741, 0.6463471039639483],
            "point_42": [0.0024862193913405355, 0.17841731833254715, 1.7845643248178],
            "point_5": [0.2932914423449599, -0.0568511065171585, 2.7557790680702556],
            "point_6": [0.32706409944674325, 0.16514312072794968, 2.7469922860685134],
            "point_7": [0.20466736510213046, 0.005858947748217496, 2.7460993219652514],
            "point_8": [0.12234930029341887, 0.33944518082699426, 2.7339834631697495],
            "point_9": [0.2518730999999884, 0.31992836440027583, 2.7618921805296233],
            "pudendal_n_branching_point": [-0.08520427451684401, 0.002990991433412445, 0.16560735905438223],
            "pudendal_n_start": [0.004438377411185904, 0.3086014494484442, 0.41184835218722776],
            "respiratory_11-1": [0.00870862926395022, -0.13956573447950937, 2.611399713924401],
            "respiratory_17-1": [0.005654647173971988, 0.11172165560925487, 1.6417798627784514],
            "respiratory_5-1": [0.00011767513857698363, -0.30968194470446886, 2.5687619383107485],
            "respiratory_8-1": [0.0014562397246949343, -0.3886006862840372, 2.522291840802067],
            "spinal_47-1": [0.016579958907853496, 0.0854346321548327, 0.9015909183902494],
            "splanchnic_n_branching_point": [0.1213971971853084, -0.2776717225479557, 0.8428440106135134],
            "splanchnic_n_start-1": [0.031010624777538583, 0.14322586852701125, 2.3280524672688845],
            "splanchnic_n_start-10": [0.007700355764729637, 0.31004384273947644, 1.3455285618363262],
            "splanchnic_n_start-11": [0.017724971494985482, 0.3205419428768455, 1.2440209171960177],
            "splanchnic_n_start-12": [-0.0006668443078337933, 0.0764094625747302, 0.9113009497835204],
            "splanchnic_n_start-13": [0.013885685293007963, 0.021302301311731785, 0.8644364898601069],
            "splanchnic_n_start-14": [0.027224841871729168, 0.4833468730079142, 0.9643500831593407],
            "splanchnic_n_start-2": [0.017663436343349662, 0.18945548016747038, 2.2378398878887715],
            "splanchnic_n_start-3": [0.02156211553274803, 0.26001008273665577, 2.1428915211384014],
            "splanchnic_n_start-4": [0.01789715785594272, 0.2849706024166771, 2.0333663289729422],
            "splanchnic_n_start-5": [0.020139339688008417, 0.2997526715941531, 1.9120670061667429],
            "splanchnic_n_start-6": [0.02388214943875268, 0.3005604438383162, 1.7998435619214619],
            "splanchnic_n_start-7": [0.019396432925925522, 0.31932243341565475, 1.683245215274538],
            "splanchnic_n_start-8": [0.02339339954198851, 0.3108591375333426, 1.5714295613650182],
            "splanchnic_n_start-9": [0.02556707992638514, 0.3084249870224194, 1.4718989375969862],
            "urinary_13-1": [0.26381665729734477, 0.010418541089591584, 0.9343432877685937],
        }

        for name, nerveCoordinatesValues in markerTermNameNerveCoordinatesMap.items():
            annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, (name, ''), isMarker=True)
            annotationGroup.createMarkerNode(nodeIdentifier, bodyCoordinates, nerveCoordinatesValues)
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
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong,
                                                       refineElementsCountAcrossMajor)

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
