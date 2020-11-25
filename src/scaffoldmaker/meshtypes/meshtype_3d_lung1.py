'''
Generates 3D lung surface mesh.
'''

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_lung1(Scaffold_base):
    '''
    3D lung scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Lung 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        if parameterSetName == 'Default':
            parameterSetName = 'Mouse 1'
        options['Base parameter set'] = parameterSetName
        options['Refine'] = False
        options['Refine number of elements'] = 4
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
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
        for key in [
            'Refine number of elements']:
            if options[key] < 1:
                options[key] = 1
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
        isMouse = 'Mouse' in parameterSetName
        isHuman = 'Human' in parameterSetName

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

        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        #rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        annotationGroups = [leftLungGroup, lungGroup]

        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        #rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)

        if isHuman:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            annotationGroups.append(upperLeftLungGroup)

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

        # common element field templates
        eftWedgeCollapseXi1_15 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
        eftWedgeCollapseXi1_26 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
        eftWedgeCollapseXi1_37 = eftfactory.createEftWedgeCollapseXi1Quadrant([3, 7])
        eftWedgeCollapseXi1_48 = eftfactory.createEftWedgeCollapseXi1Quadrant([4, 8])
        eftWedgeCollapseXi1_57 = eftfactory.createEftWedgeCollapseXi1Quadrant([5, 7])
        eftWedgeCollapseXi1_68 = eftfactory.createEftWedgeCollapseXi1Quadrant([6, 8])
        eftWedgeCollapseXi2_56 = eftfactory.createEftWedgeCollapseXi2Quadrant([5, 6])
        eftWedgeCollapseXi2_78 = eftfactory.createEftWedgeCollapseXi2Quadrant([7, 8])
        eftTetCollapseXi1Xi2_71 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(7, 1)
        eftTetCollapseXi1Xi2_82 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(8, 2)
        eftTetCollapseXi1Xi2_63 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(6, 3)
        eftTetCollapseXi1Xi2_53 = eftfactory.createEftTetrahedronCollapseXi1Xi2Quadrant(5, 3)

        if isHuman:
            lelementsCount1 = 2
            lelementsCount2 = 4
            lelementsCount3 = 3

            uelementsCount1 = 2
            uelementsCount2 = 4
            uelementsCount3 = 4

            # Create nodes
            nodeIdentifier = 1
            lNodeIds = []
            uNodeIds = []
            d1 = [1.0, 0.0, 0.0]
            d2 = [0.0, 1.0, 0.0]
            d3 = [0.0, 0.0, 1.0]

            # Lower lobe nodes
            for n3 in range(lelementsCount3 + 1):
                lNodeIds.append([])
                for n2 in range(lelementsCount2 + 1):
                    lNodeIds[n3].append([])
                    for n1 in range(lelementsCount1 + 1):
                        lNodeIds[n3][n2].append(None)
                        if (n1 == 0 or n1 == lelementsCount1) and (n2 == 0):
                            continue
                        if (n3 > (lelementsCount3 - 2)) and (n2 > (lelementsCount2 - 2)):
                            continue

                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        x = [1.0 * (n1 - 1), 1.0 * (n2 - 1), 1.0 * n3]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        lNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Upper lobe nodes
            for n3 in range(uelementsCount3 + 1):
                reset = 0
                uNodeIds.append([])
                for n2 in range(uelementsCount2 + 1):
                    uNodeIds[n3].append([])
                    for n1 in range(uelementsCount1 + 1):
                        uNodeIds[n3][n2].append(None)
                        if (n1 == 0 or n1 == uelementsCount1) and (n2 == 0 or n2 == uelementsCount2):
                            continue
                        if (n2 < uelementsCount2 - 2) and (n3 < uelementsCount3 - 2):
                            continue
                        if (n2 == 0 or n2 == uelementsCount2) and (n3 == uelementsCount3):
                            continue
                        if (n1 == 0 or n1 == uelementsCount1) and (n3 == uelementsCount3):
                            continue

                        # Oblique fissure nodes
                        if reset < 3 and n3 < uelementsCount3 - 2:
                            uNodeIds[n3][n2][n1] = lNodeIds[n3][uelementsCount2][n1]
                            reset += 1
                            continue
                        elif (reset < 7) and (n3 == uelementsCount3 - 2):
                            uNodeIds[n3][n2][n1] = lNodeIds[lelementsCount3][n2][n1]
                            reset += 1
                            continue

                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        x = [1.0 * (n1 - 1), 1.0 * (n2 - 1) + 1.5, 1.0 * n3 + 1.0]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        uNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Create elements
            # Lower lobe elements
            elementIdentifier = 1
            for e3 in range(lelementsCount3):
                for e2 in range(lelementsCount2):
                    for e1 in range(lelementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            lNodeIds[e3    ][e2][e1], lNodeIds[e3    ][e2][e1 + 1], lNodeIds[e3    ][e2 + 1][e1], lNodeIds[e3    ][e2 + 1][e1 + 1],
                            lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                        if (e2 == 0) and (e1 == 0):
                            # Back wedge elements
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eftWedgeCollapseXi1_15
                        elif (e2 == 0) and (e1 == lelementsCount1-1):
                            # Back wedge elements
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(1)
                            eft = eftWedgeCollapseXi1_26
                        elif e3 == 1 and e2 == (lelementsCount2 - 2):
                            # Middle wedge
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(6)
                            eft = eftWedgeCollapseXi2_78
                        elif (e3 == (lelementsCount3 - 1)) and (e2 == (lelementsCount2 - 3)):
                            # Remapped cube element 1
                            eft = eftfactory.createEftBasic()
                            setEftScaleFactorIds(eft, [1], [])
                            remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif (e3 == (lelementsCount3 - 1)) and (e2 == (lelementsCount2 - 2)):
                            # Remapped cube element 2
                            nodeIdentifiers[2] = lNodeIds[e3 - 1][e2 + 1][e1    ]
                            nodeIdentifiers[3] = lNodeIds[e3 - 1][e2 + 1][e1 + 1]
                            nodeIdentifiers[6] = lNodeIds[e3 - 1][e2 + 2][e1    ]
                            nodeIdentifiers[7] = lNodeIds[e3 - 1][e2 + 2][e1 + 1]
                            eft = eftfactory.createEftBasic()
                            setEftScaleFactorIds(eft, [1], [])
                            remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif None in nodeIdentifiers:
                            continue

                        # print('element ', elementIdentifier, '|| ', nodeIdentifiers)
                        if eft is eftRegular:
                            element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                        else:
                            elementtemplateCustom.defineField(coordinates, -1, eft)
                            element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                        element.setNodesByIdentifier(eft, nodeIdentifiers)
                        if eft.getNumberOfLocalScaleFactors() == 1:
                            element.setScaleFactors(eft, [-1.0])
                        elementIdentifier += 1
                        leftLungMeshGroup.addElement(element)
                        lungMeshGroup.addElement(element)
                        lowerLeftLungMeshGroup.addElement(element)

            # Upper lobe elements
            for e3 in range(uelementsCount3):
                for e2 in range(uelementsCount2):
                    for e1 in range(uelementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            uNodeIds[e3][e2][e1], uNodeIds[e3][e2][e1 + 1], uNodeIds[e3][e2 + 1][e1],
                            uNodeIds[e3][e2 + 1][e1 + 1],
                            uNodeIds[e3 + 1][e2][e1], uNodeIds[e3 + 1][e2][e1 + 1], uNodeIds[e3 + 1][e2 + 1][e1],
                            uNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                        if (e3 < uelementsCount3 - 1) and (e2 == uelementsCount2-1) and (e1 == 0):
                            # Distal-front wedge elements
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(2)
                            eft = eftWedgeCollapseXi1_37
                        elif (e3 < uelementsCount3 - 1) and (e2 == uelementsCount2-1) and (e1 == uelementsCount1 - 1):
                            # Distal-back wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(3)
                            eft = eftWedgeCollapseXi1_48
                        elif (e3 == uelementsCount3 - 2) and (e2 == 0) and (e1 == 0):
                            # Medial-front wedge elements
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eftWedgeCollapseXi1_15
                        elif (e3 == uelementsCount3 - 2) and (e2 == 0) and (e1 == uelementsCount1 - 1):
                            # Medial-back wedge elements
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(1)
                            eft = eftWedgeCollapseXi1_26
                        elif (e3 == uelementsCount3 - 1) and (0 < e2 < uelementsCount2-1) and (e1 == 0):
                            # Top-front wedge elements
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(4)
                            eft = eftWedgeCollapseXi1_57
                        elif (e3 == uelementsCount3 - 1) and (0 < e2 < uelementsCount2-1) and (e1 == uelementsCount1 - 1):
                            # Top-back wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(5)
                            eft = eftWedgeCollapseXi1_68
                        elif (e3 == uelementsCount3 - 1) and (e2 == 0) and (e1 == 0):
                            # Top-front-medial tetrahedron wedge elements
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eftTetCollapseXi1Xi2_82
                        elif (e3 == uelementsCount3 - 1) and (e2 == 0) and (e1 == uelementsCount1 - 1):
                            # Top-back-medial tetrahedron wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(1)
                            eft = eftTetCollapseXi1Xi2_71
                        elif (e3 == uelementsCount3 - 1) and (e2 == uelementsCount2 - 1) and (e1 == 0):
                            # Top-front-distal tetrahedron wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(2)
                            eft = eftTetCollapseXi1Xi2_63
                        elif (e3 == uelementsCount3 - 1) and (e2 == uelementsCount2 - 1) and (e1 == uelementsCount1 - 1):
                            # Top-front-distal tetrahedron wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(3)
                            eft = eftTetCollapseXi1Xi2_53
                        elif (e3 == uelementsCount3 - 2) and (e2 == uelementsCount2 - 3):
                            # Remapped cube element 1
                            eft = eftfactory.createEftBasic()
                            setEftScaleFactorIds(eft, [1], [])
                            remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                        elif (e3 == uelementsCount3 - 2) and (e2 == uelementsCount2 - 2):
                            # Remapped cube element 2
                            eft = eftfactory.createEftBasic()
                            setEftScaleFactorIds(eft, [1], [])
                            remapEftNodeValueLabel(eft, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                        elif None in nodeIdentifiers:
                            continue

                        # print('element ', elementIdentifier, '|| ', nodeIdentifiers)
                        if eft is eftRegular:
                            element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                        else:
                            elementtemplateCustom.defineField(coordinates, -1, eft)
                            element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                        element.setNodesByIdentifier(eft, nodeIdentifiers)
                        if eft.getNumberOfLocalScaleFactors() == 1:
                            element.setScaleFactors(eft, [-1.0])
                        elementIdentifier += 1
                        leftLungMeshGroup.addElement(element)
                        lungMeshGroup.addElement(element)
                        upperLeftLungMeshGroup.addElement(element)

        elif isMouse:
            elementsCount1 = 2
            elementsCount2 = 4
            elementsCount3 = 4

            # Create nodes
            nodeIdentifier = 1
            lNodeIds = []
            d1 = [0.5, 0.0, 0.0]
            d2 = [0.0, 0.5, 0.0]
            d3 = [0.0, 0.0, 1.0]
            for n3 in range(elementsCount3 + 1):
                lNodeIds.append([])
                for n2 in range(elementsCount2 + 1):
                    lNodeIds[n3].append([])
                    for n1 in range(elementsCount1 + 1):
                        lNodeIds[n3][n2].append(None)
                        if n3 < elementsCount3:
                            if (n1 == 0) and ((n2 == 0) or (n2 == elementsCount2)):
                                continue
                        else:
                            if (n2 == 0) or (n2 == elementsCount2) or (n1 == 0):
                                continue
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        x = [0.5 * (n1 - 1), 0.5 * (n2 - 1), 1.0 * n3]
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        lNodeIds[n3][n2][n1] = nodeIdentifier
                        nodeIdentifier += 1

            # Create elements
            elementIdentifier = 1
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1):
                        eft = eftRegular
                        nodeIdentifiers = [
                            lNodeIds[e3    ][e2][e1], lNodeIds[e3    ][e2][e1 + 1], lNodeIds[e3    ][e2 + 1][e1], lNodeIds[e3    ][e2 + 1][e1 + 1],
                            lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1]]

                        if (e3 < elementsCount3 - 1):
                            if (e2 == 0) and (e1 == 0):
                                # Back wedge elements
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(0)
                                eft = eftWedgeCollapseXi1_15
                            elif (e2 == elementsCount2 - 1) and (e1 == 0):
                                # Front wedge elements
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(2)
                                eft = eftWedgeCollapseXi1_37
                        else:
                            if (e2 == 0) and (e1 == 1):
                                # Top back wedge elements
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(4)
                                eft = eftWedgeCollapseXi2_56
                            elif (e2 == elementsCount2 - 1) and (e1 == 1):
                                # Top front wedge elements
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                eft = eftWedgeCollapseXi2_78
                            elif ((0 < e2 < (elementsCount2 - 1))) and (e1 == 0):
                                # Top middle back wedge element
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                eft = eftWedgeCollapseXi1_57
                            elif (e2 == 0) and (e1 == 0):
                                # Top back tetrahedron element
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(5)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(0)
                                eft = eftTetCollapseXi1Xi2_82
                            elif (e2 == elementsCount2 - 1) and (e1 == 0):
                                # Top front tetrahedron element
                                nodeIdentifiers.pop(7)
                                nodeIdentifiers.pop(6)
                                nodeIdentifiers.pop(4)
                                nodeIdentifiers.pop(2)
                                eft = eftTetCollapseXi1Xi2_63

                        if eft is eftRegular:
                            element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                        else:
                            elementtemplateCustom.defineField(coordinates, -1, eft)
                            element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                        element.setNodesByIdentifier(eft, nodeIdentifiers)
                        if eft.getNumberOfLocalScaleFactors() == 1:
                            element.setScaleFactors(eft, [-1.0])
                        elementIdentifier += 1
                        leftLungMeshGroup.addElement(element)
                        lungMeshGroup.addElement(element)

            # Apex annotation point
            idx = elementsCount1 * elementsCount2 * (elementsCount3 - 1) + elementsCount1 * (elementsCount2 // 2)
            element1 = mesh.findElementByIdentifier(idx)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, 'apex of left lung')
            markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 1.0])

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
