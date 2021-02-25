"""
Generates a 2-D tube bifurcation mesh.
"""

# from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
# from scaffoldmaker.annotation.lung_terms import get_lung_term
from __future__ import division
import math
from opencmiss.utils.maths.vectorops import cross, normalize
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import get_curve_circle_points, \
    make_tube_bifurcation_points, make_tube_bifurcation_elements_2d
from scaffoldmaker.utils.interpolation import getCubicHermiteArcLength


class MeshType_2d_tubebifurcationtree1(Scaffold_base):
    '''
    Generates a 2-D tube bifurcation tree mesh.
    '''

    @staticmethod
    def getName():
        return '2D Tube Bifurcation Tree 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {}
        options['Bifurcation tree'] = ScaffoldPackage(MeshType_1d_bifurcationtree1)
        options['Number of elements around root'] = 6
        options['Maximum element length'] = 0.25
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Bifurcation tree',
            'Number of elements around root',
            'Maximum element length']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Bifurcation tree':
            return [ MeshType_1d_bifurcationtree1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
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
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Bifurcation tree':
            if not parameterSetName:
                parameterSetName = MeshType_1d_bifurcationtree1.getParameterSetNames()[0]
            return ScaffoldPackage(MeshType_1d_bifurcationtree1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = False
        for key in [
            'Number of elements around root']:
            if options[key] < 4:
                options[key] = 4
        if (options['Maximum element length'] <= 0.0) or (options['Maximum element length'] >= 0.4):
            options['Maximum element length'] = 0.25
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base bicubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        bifurcationTreeScaffold = options['Bifurcation tree']
        maxElementLength = options['Maximum element length']
        elementsCountAroundRoot = options['Number of elements around root']
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        mesh = fm.findMeshByDimension(2)
        elementtemplateStd = mesh.createElementtemplate()
        elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eftStd = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eftStd.setFunctionNumberOfTerms(n * 4 + 4, 0)
        elementtemplateStd.defineField(coordinates, -1, eftStd)

        # airwayGroup = AnnotationGroup(region, get_lung_term("airway"))
        # airwayMeshGroup = airwayGroup.getMeshGroup(mesh)
        # annotationGroups = [airwayGroup]
        # tracheaGroup = AnnotationGroup(region, get_lung_term("trachea"))
        # tracheaMeshGroup = tracheaGroup.getMeshGroup(mesh)
        # annotationGroups.append(tracheaGroup)
        # rightBronchusGroup = AnnotationGroup(region, get_lung_term("right bronchus"))
        # rightBronchusMeshGroup = rightBronchusGroup.getMeshGroup(mesh)
        # annotationGroups.append(rightBronchusGroup)
        # leftBronchusGroup = AnnotationGroup(region, get_lung_term("left bronchus"))
        # leftBronchusMeshGroup = leftBronchusGroup.getMeshGroup(mesh)
        # annotationGroups.append(leftBronchusGroup)

        ##############
        # Create nodes
        ##############
        child1 = []
        cx = []
        cxd1 = []
        cxd2 = []
        x = []
        d1 = []
        d2 = []
        nextRootNodeId = []

        # Create bifurcation tree
        bifurcationTree = bifurcationTreeScaffold.getScaffoldType().generateBifurcationTree(bifurcationTreeScaffold.getScaffoldSettings())
        Generation = bifurcationTreeScaffold.getScaffoldSettings()['Number of generations']

        # Generate bifurcation tree nodes and elements and return a next nodeIdentifier
        # nodeIdentifier = bifurcationTree.generateZincModel(region)[0]
        nodeIdentifier = 1
        elementIdentifier = 1

        # Create a node template for the thickness - circle nodes around the bifurcation tree
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        # Get a bifurcation root node
        rootNode = bifurcationTree.getRootNode()

        # Create a branch at each generation
        for j in range(0, Generation - 1):
            child1.append([])
            # Create a bifurcation at each branch
            for k in range(2 ** j):
                child1[j].append(None)
                xStem = []
                d1Stem = []
                d2Stem = []
                index = 0
                if j == 0:
                    child1[0][0] = rootNode.getChild(0)
                    xParent1, xParentd1, rParent1, xParent2, xParentd2, rParent2 = rootNode.getChildCurve(0)
                    xRight1, xRightd1, rRight1, xRight2, xRightd2, rRight2 = child1[0][0].getChildCurve(0)
                    xLeft1, xLeftd1, rLeft1, xLeft2, xLeftd2, rLeft2 = child1[0][0].getChildCurve(1)
                else:
                    ind = math.floor(k/2)
                    if (k % 2) == 0:
                        child1[j][k] = child1[j-1][ind].getChild(0)
                    else:
                        child1[j][k] = child1[j-1][ind].getChild(1)
                    # Parental branch
                    xParent1, xParentd1, rParent1, xParent2, xParentd2, rParent2 = child1[j-1][ind].getChildCurve(k % 2)
                    # Get the right(0) [left(1)] children curve from the bifurcation tree
                    xRight1, xRightd1, rRight1, xRight2, xRightd2, rRight2 = child1[j][k].getChildCurve(0)
                    xLeft1, xLeftd1, rLeft1, xLeft2, xLeftd2, rLeft2 = child1[j][k].getChildCurve(1)

                # Set up the curve axis for bifurcation
                axis1 = normalize(xRightd2)
                axis2 = normalize(xLeftd2)
                side = cross(axis1, axis2)
                side[0] = side[0] * -1
                side[1] = side[1] * -1
                side[2] = side[2] * -1

                # Look at each bifurcation (parent and children branch)
                for l in range(3):
                    cx.append([])
                    cxd1.append([])
                    cxd2.append([])
                    x.append([])
                    d1.append([])
                    d2.append([])
                    if l == 0:
                        x1  = xParent1
                        xd1 = xParentd1
                        x2  = xParent2
                        xd2 = xParentd2
                        r1  = rParent1
                        r2  = rParent2
                        rd1 = rd2 = (rParent2 - rParent1)
                    else:
                        x1  = xRight1   if l == 1 else xLeft1
                        xd1 = xRightd1  if l == 1 else xLeftd1
                        x2  = xRight2   if l == 1 else xLeft2
                        xd2 = xRightd2  if l == 1 else xLeftd2
                        r1  = rRight1   if l == 1 else rLeft1
                        rd1 = rd2 = (rRight2 - rRight1) if l == 1 else (rLeft2 - rLeft1) # The radius of the circle between the parental nodes
                        r2  = rRight1   if l == 1 else rLeft1

                    # Curve length between the parental nodes and divided into 2 bifurcation points
                    curveLength = getCubicHermiteArcLength(x1, xd1, x2, xd2)
                    elementsCount = max(2, math.ceil(curveLength / maxElementLength))
                    elementLength = curveLength / elementsCount
                    xi = 1 / elementsCount

                    # Get ring nodes along the curve
                    if l == 0:
                        # Parent branch
                        cx[l], cxd1[l], cxd2[l], x[l], d1[l], d2[l] = get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, 1-xi, elementLength, side, elementsCountAroundRoot)
                        # Add ring nodes along the root
                        if j == 0:
                            rootNodeId = []
                            rootNodeId.append([])
                            xRoot = []
                            d1Root = []
                            d2Root = []
                            for m in range(0, elementsCount-1):
                                xRoot.append([])
                                d1Root.append([])
                                d2Root.append([])
                                xi = m / elementsCount
                                cxStem, cxd1Stem, cxd2Stem, xRoot[index], d1Root[index], d2Root[index] = \
                                    get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side,
                                                            elementsCountAroundRoot)
                                for n in range(elementsCountAroundRoot):
                                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                                    cache.setNode(node)
                                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xRoot[index][n])
                                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Root[index][n])
                                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Root[index][n])
                                    rootNodeId[0].append(nodeIdentifier)
                                    nodeIdentifier = nodeIdentifier + 1
                    else:
                        # Children branch
                        cx[l], cxd1[l], cxd2[l], x[l], d1[l], d2[l] = get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side, elementsCountAroundRoot)
                        # Add ring nodes along the children branch
                        elementsCountBranch = elementsCount + 1 if j == Generation - 2 else elementsCount - 1
                        for m in range(2, elementsCountBranch):
                            xStem.append([])
                            d1Stem.append([])
                            d2Stem.append([])
                            xi = m / elementsCount
                            cxStem, cxd1Stem, cxd2Stem, xStem[index], d1Stem[index], d2Stem[index] = \
                                get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side, elementsCountAroundRoot)
                            print(index)
                            index += 1

                rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
                    make_tube_bifurcation_points(cx[0], x[0], d2[0], cx[2], x[2], d2[2], cx[1], x[1], d2[1])

                paNodeId = []
                for n in range(elementsCountAroundRoot):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x [0][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[0][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[0][n])
                    paNodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1

                roNodeId = []
                for n in range(len(rox)):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox [n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
                    roNodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1

                coNodeId = []
                for n in range(len(cox)):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox [n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
                    coNodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1

                c1NodeId = []
                for n in range(elementsCountAroundRoot):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x [1][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[1][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[1][n])
                    c1NodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1

                c2NodeId = []
                for n in range(elementsCountAroundRoot):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x [2][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[2][n])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[2][n])
                    c2NodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1

                stemNodeId = []
                for m in range(len(xStem)):
                    for n in range(elementsCountAroundRoot):
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xStem[m][n])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Stem[m][n])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Stem[m][n])
                        stemNodeId.append(nodeIdentifier)
                        nodeIdentifier = nodeIdentifier + 1

            ################
            # Create elements
            ################
                # Create an element from the bottom to the top
                elementsCountAlongRoot = math.ceil(len(rootNodeId[k])/elementsCountAroundRoot - 1)
                for e2 in range(elementsCountAlongRoot+1):
                    for e1 in range(elementsCountAroundRoot):
                        if e2 == (elementsCountAlongRoot):
                            # Connect root and bifurcation parents
                            bni1 = e2 * elementsCountAroundRoot + e1  if len(rootNodeId[k]) > elementsCountAroundRoot else e1
                            bni2 = e2 * elementsCountAroundRoot + (e1 + 1) % elementsCountAroundRoot  if len(rootNodeId[k]) > elementsCountAroundRoot else (e1 + 1) % elementsCountAroundRoot
                            bni3 = e1 % elementsCountAroundRoot if j == 0 else (e1 + 2) % elementsCountAroundRoot
                            bni4 = (e1 + 1) % elementsCountAroundRoot if j == 0 else (e1 + 3) % elementsCountAroundRoot
                            nodeIdentifiers = [rootNodeId[k][bni1], rootNodeId[k][bni2], paNodeId[bni3], paNodeId[bni4]]
                        elif len(rootNodeId[k]) != elementsCountAroundRoot:
                            # Connect root
                            bni1 = e2 * elementsCountAroundRoot + e1
                            bni2 = e2 * elementsCountAroundRoot + (e1 + 1) % elementsCountAroundRoot
                            nodeIdentifiers = [rootNodeId[k][bni1], rootNodeId[k][bni2], rootNodeId[k][bni1 + elementsCountAroundRoot], rootNodeId[k][bni2 + elementsCountAroundRoot]]

                        element = mesh.createElement(elementIdentifier, elementtemplateStd)
                        result = element.setNodesByIdentifier(eftStd, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
                        nodeIdentifiers = []

                elementIdentifier = make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier,
                                                                      paNodeId, paStartIndex, c2NodeId,
                                                                      c2StartIndex, c1NodeId, c1StartIndex,
                                                                      roNodeId, coNodeId,
                                                                      useCrossDerivatives)
                # Store the root nodes for the next generation
                for i in range(2):
                    nextRootNodeId.append([])
                    ind = elementsCountAroundRoot * (elementsCount - 1) if j == Generation - 2 else elementsCountAroundRoot * (elementsCount - 3)
                    if i == 0:
                        nextRootNodeId[k*2] = c1NodeId
                        nextRootNodeId[k*2].extend(stemNodeId[:ind])
                    else:
                        nextRootNodeId[(k*2)+1] = c2NodeId
                        nextRootNodeId[(k*2)+1].extend(stemNodeId[ind:])

            # Connect children branch (last branch)
            if j == Generation - 2:
                for e3 in range(len(nextRootNodeId)):
                    for e2 in range(elementsCount - 1):
                        for e1 in range(elementsCountAroundRoot):
                            bni1 = e2 * elementsCountAroundRoot + e1
                            bni2 = e2 * elementsCountAroundRoot + (e1 + 1) % elementsCountAroundRoot
                            nodeIdentifiers = [nextRootNodeId[e3][bni1], nextRootNodeId[e3][bni2], nextRootNodeId[e3][bni1 + elementsCountAroundRoot], nextRootNodeId[e3][bni2 + elementsCountAroundRoot]]
                            element = mesh.createElement(elementIdentifier, elementtemplateStd)
                            result = element.setNodesByIdentifier(eftStd, nodeIdentifiers)
                            elementIdentifier = elementIdentifier + 1

            rootNodeId = nextRootNodeId
            nextRootNodeId = []

        return []
