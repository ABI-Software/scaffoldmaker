'''
Generates 3D lung surface mesh.
'''

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
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
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        options['Base parameter set'] = parameterSetName

        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = []
        return optionNames

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''
        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        elementsCount1 = 2
        elementsCount2 = 4
        elementsCount3 = 4

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

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
                    lNodeIds[n3][n2].append([])
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
        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        annotationGroups = [leftLungGroup, rightLungGroup, lungGroup]

        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)

        eft1 = eftfactory.createEftWedgeCollapseXi1AtXi2Zero()
        eft2 = eftfactory.createEftWedgeCollapseXi1AtXi2One()
        eft3 = eftfactory.createEftWedgeCollapseXi2RightAtXi3One()
        eft4 = eftfactory.createEftWedgeCollapseXi2LeftAtXi3One()
        eft5 = eftfactory.createEftWedgeCollapseXi1AtXi3One()
        eft6 = eftfactory.createEftTetrahedronCollapseXi1Xi2AtXi3OneXi1RightAtXi2Zero()
        eft7 = eftfactory.createEftTetrahedronCollapseXi1Xi2AtXi3OneXi1LeftAtXi2Zero()

        elementIdentifier = 1
        for e3 in range(elementsCount3):
            for e2 in range(elementsCount2):
                for e1 in range(elementsCount1):
                    eft = eftRegular
                    nodeIdentifiers = [
                        lNodeIds[e3    ][e2][e1], lNodeIds[e3    ][e2][e1 + 1], lNodeIds[e3    ][e2 + 1][e1], lNodeIds[e3    ][e2 + 1][e1 + 1],
                        lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1]]
                    scalefactors = None
                    if (e3 < elementsCount3 - 1):
                        if (e2 == 0) and (e1 == 0):
                            # Back wedge elements
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eft1
                            scalefactors = [-1.0]
                        elif (e2 == elementsCount2 - 1) and (e1 == 0):
                            # Front wedge elements
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(2)
                            eft = eft2
                    else:
                        if (e2 == 0) and (e1 == 1):
                            # Top back wedge elements
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(4)
                            eft = eft3
                        elif (e2 == elementsCount2 - 1) and (e1 == 1):
                            # Top front wedge elements
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(6)
                            eft = eft4
                            scalefactors = [-1.0]
                        elif (e2 == 1) and (e1 == 0):
                            # Top middle back wedge element
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(4)
                            eft = eft5
                        elif (e2 == 2) and (e1 == 0):
                            # Top middle front wedge element
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(4)
                            eft = eft5
                        if (e2 == 0) and (e1 == 0):
                            # Top back tetrahedron element
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eft6
                            scalefactors = [-1.0]
                        if (e2 == elementsCount2 - 1) and (e1 == 0):
                            # Top front tetrahedron element
                            nodeIdentifiers.pop(7)
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(2)
                            eft = eft7
                            scalefactors = [-1.0]

                    if eft is eftRegular:
                        element = mesh.createElement(elementIdentifier, elementtemplateRegular)
                    else:
                        elementtemplateCustom.defineField(coordinates, -1, eft)
                        element = mesh.createElement(elementIdentifier, elementtemplateCustom)
                    element.setNodesByIdentifier(eft, nodeIdentifiers)
                    if scalefactors:
                        element.setScaleFactors(eft, scalefactors)
                    elementIdentifier += 1
                    leftLungMeshGroup.addElement(element)
                    lungMeshGroup.addElement(element)

        # Apex annotation point
        idx = elementsCount1 * elementsCount2 * (elementsCount3 - 1) + elementsCount1 * (elementsCount2 // 2)
        element1 = mesh.findElementByIdentifier(idx)
        markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
        nodeIdentifier += 1
        cache.setNode(markerPoint)
        markerName.assignString(cache, 'APEX')
        markerLocation.assignMeshLocation(cache, element1, [1.0, 1.0, 1.0])

        fm.endChange()
        return annotationGroups

