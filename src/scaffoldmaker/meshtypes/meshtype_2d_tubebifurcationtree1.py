"""
Generates a 2-D tube bifurcation mesh.
"""

from __future__ import division
import math

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.lung_terms import get_lung_term
from opencmiss.maths.vectorops import cross, normalize, magnitude
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1, TreeNode, BifurcationTree
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import get_curve_circle_points, \
    make_tube_bifurcation_points, make_tube_bifurcation_elements_2d, move_1D_central_path
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
        options['Maximum element length'] = 0.5
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
        if options['Maximum element length'] <= 0.0:
            options['Maximum element length'] = 1.0
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

        # Bifurcation tree extraction
        fieldParameters = move_1D_central_path(region, bifurcationTreeScaffold)
        bifurcationTree = bifurcationTreeScaffold.getScaffoldType(). \
            generateBifurcationTree(bifurcationTreeScaffold.getScaffoldSettings())
        generationCount = bifurcationTreeScaffold.getScaffoldSettings()['Number of generations']

        # reconstruct 1D tree from the extracted field parameters
        rootNode = fieldParameters.pop(0)
        x_rootNode = rootNode[1][0][0]
        d1_rootNode = rootNode[1][1][0]
        r_rootNode = rootNode[1][2][0][0]
        bifurcationTree._rootNode = []
        bifurcationTree._rootNode = TreeNode(x_rootNode, d1_rootNode, r_rootNode)
        nextChildNode = fieldParameters.pop(0)
        x_nextChildNode = nextChildNode[1][0][0]
        d1_nextChildNode = nextChildNode[1][1]
        r_nextChildNode = nextChildNode[1][2][0]
        forkNormal = normalize(cross([0.0, 1.0, 0.0], d1_nextChildNode[0]))
        bifurcationTree._rootNode.addChild(createNodeTreeFromParameters(fieldParameters, 1, generationCount, x_nextChildNode,
                                d1_nextChildNode, r_nextChildNode, forkNormal))
        # [nodeIdentifier, elementIdentifier] = bifurcationTree.generateZincModel(region)  # 1D bifurcation tree
        airway2D = TubeAlongTree(region, bifurcationTree, maxElementLength, elementsCountAroundRoot)

        ##############
        # Create nodes
        ##############
        nodeIdentifier = 1
        nodeIdentifier, lastRing = airway2D.createNodesAround1DTree(bifurcationTree._rootNode, [], nodeIdentifier)

        #################
        # Create elements
        #################
        elementIdentifier = 1
        elementIdentifier = airway2D.create2DElementAlong1DTree(elementIdentifier)

        return []

class TubeAlongTree:
    '''
    Create a 2D tube along 1D tree curve.
    '''

    def __init__(self, region, bifurcationTree, maxElementLength, elementsCountAroundRoot):
        """
        :param bifurcationTree:
        :param maxElementLength:
        :param elementsCountAroundRoot:
        :param nodeIdentifier:
        """
        self._maxElementLength = maxElementLength
        self._elementsCountAround = elementsCountAroundRoot
        self._bifurcationTree = bifurcationTree
        self._rootNode = bifurcationTree._rootNode

        self._region = region
        self._fm = region.getFieldmodule()
        self._coordinates = findOrCreateFieldCoordinates(self._fm)
        self._cache = self._fm.createFieldcache()

        # Create nodes around the 1D bifurcation tree
        self._nodes = self._fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._nodetemplate = self._nodes.createNodetemplate()
        self._nodetemplate.defineField(self._coordinates)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        # Create elements around the 1D bifurcation tree
        self._mesh = self._fm.findMeshByDimension(2)
        self._elementtemplateStd = self._mesh.createElementtemplate()
        self._elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        bicubicHermiteBasis = self._fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        self._eftStd = self._mesh.createElementfieldtemplate(bicubicHermiteBasis)
        for n in range(4):
            self._eftStd.setFunctionNumberOfTerms(n * 4 + 4, 0)
        self._elementtemplateStd.defineField(self._coordinates, -1, self._eftStd)

        # Index for bifurcation
        self._countTubeElements = []
        self._paStartIndex = []
        self._c1StartIndex = []
        self._c2StartIndex = []
        self._bifurcationSide = []
        self._tubeNodeIdentifersList = []
        self._paNodeIdentifersList = []
        self._c1NodeIdentifersList = []
        self._c2NodeIdentifersList = []
        self._c1NodeIdentifersList_temp = []
        self._c2NodeIdentifersList_temp = []
        self._coNodeIdentifersList = []
        self._roNodeIdentifersList = []

    def createNodesAround1DTree(self, currentNode, previousRingNodes, nodeIdentifier):
        """
        Create nodes around 1D tree network using depth-first traversal (recursive method)
        :param currentNode:
        :param nodeIdentifier:
        :return:
        """
        lastRingNodes = previousRingNodes
        if currentNode != -1:
            currentNode_iter = iter(currentNode._children)
            currentChild = next(currentNode_iter, -1)
            # Left branch
            if currentChild != -1:
                if len(currentNode._children) > 1:
                    nodeIdentifier, c1RingNodes, c2RingNodes = self.createBifurcationNodes(lastRingNodes, currentNode, nodeIdentifier)
                    c1Node_temp = self._c1NodeIdentifersList_temp.pop(0)
                    c2Node_temp = self._c2NodeIdentifersList_temp.pop(-1)
                    self._tubeNodeIdentifersList.append(c1Node_temp)
                nodeIdentifier, firstRingNodes, lastRingNodes = self.create2DtubeNodesAlong1Dtree(currentNode.getChild(0), currentNode.getChildCurve(0), nodeIdentifier)
                # update nodeIdentifiers at bifurcation-to-bifurcation-elements
                if lastRingNodes['nodeIdentifiers'] == []:
                    lastRingNodes['nodeIdentifiers'] = c1RingNodes['nodeIdentifiers']
                nodeIdentifier, lastRingNodes = self.createNodesAround1DTree(currentChild, lastRingNodes, nodeIdentifier)
            # Right branch
            currentChild = next(currentNode_iter, -1)
            if currentChild != -1:
                self._tubeNodeIdentifersList.append(c2Node_temp)
                nodeIdentifier, _, lastRingNodes = self.create2DtubeNodesAlong1Dtree(currentNode.getChild(1), currentNode.getChildCurve(1), nodeIdentifier)
                # update nodeIdentifiers at bifurcation-to-bifurcation-elements
                if lastRingNodes['nodeIdentifiers'] == []:
                    lastRingNodes['nodeIdentifiers'] = c2RingNodes['nodeIdentifiers']
                nodeIdentifier, lastRingNodes = self.createNodesAround1DTree(currentChild, lastRingNodes, nodeIdentifier)

        return nodeIdentifier, lastRingNodes

    def create2DtubeNodesAlong1Dtree(self, childNode, p2cNodeCurve, nodeIdentifier, isTube=True):
        """
        Generate 2D tube nodes along 1D tree network for both tube element and bifurcation children
        :param childNode:
        :param p2cNodeCurve:
        :param nodeIdentifier:
        :param isTube:
        :return:
        """
        firstRing = True
        firstRingNodes = {}
        lastRingNodes = {}

        parentx, parentd1, parentr, childx, childd1, childr = p2cNodeCurve

        side = self.findSide(childNode)

        curveLength = getCubicHermiteArcLength(parentx, parentd1, childx, childd1)
        elementsCount_temp = max(2, math.ceil(curveLength / self._maxElementLength))
        elementLength = curveLength / elementsCount_temp
        rd1 = rd2 = childr-parentr

        startIndex = 0 if parentx == self._rootNode._x else 2
        endIndex = elementsCount_temp + 1 if childNode._children == [] else elementsCount_temp
        if isTube == False:
            startIndex = 1
            endIndex = startIndex + 1

        # Find a parent ring at bifurcation-to-bifurcation element
        if startIndex == endIndex:
            xi = 1 / elementsCount_temp
            x, d1, d2 = get_curve_circle_points(parentx, parentd1, childx, childd1, parentr, rd1, childr, rd2,
                                        xi, elementLength, side, self._elementsCountAround)
            firstRingNodes = lastRingNodes = {
                'x': x,
                'd1': d1,
                'd2': d2,
                'side': side,
                'nodeIdentifiers': []
            }

        for m in range(startIndex, endIndex):
            xi = m / elementsCount_temp
            x, d1, d2 = get_curve_circle_points(parentx, parentd1, childx, childd1, parentr, rd1, childr, rd2,
                                        xi, elementLength, side, self._elementsCountAround)
            nodeIdentifier_temp = []
            for n in range(self._elementsCountAround):
                node = self._nodes.createNode(nodeIdentifier, self._nodetemplate)
                self._cache.setNode(node)
                self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_VALUE, 1, x[n])
                self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[n])
                self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[n])
                nodeIdentifier_temp.append(nodeIdentifier)
                nodeIdentifier += 1
            # Store variables
            if isTube == True:
                self._tubeNodeIdentifersList.append(nodeIdentifier_temp)
            if firstRing == True:
                firstRing = False
                firstRingNodes = {
                    'x': x,
                    'd1': d1,
                    'd2': d2,
                    'side': side,
                    'nodeIdentifiers': nodeIdentifier_temp
                }
            if m == (endIndex - 1):
                lastRingNodes = {
                    'x': x,
                    'd1': d1,
                    'd2': d2,
                    'side': side,
                    'nodeIdentifiers': nodeIdentifier_temp
                }
        # Store number of tube elements else None in bifurcation nodes
        if isTube == True:
            countTubeElement = (endIndex-startIndex) if (parentx == self._rootNode._x) else (endIndex - startIndex + 1)
            self._countTubeElements.append(countTubeElement)

        return nodeIdentifier, firstRingNodes, lastRingNodes

    def createBifurcationNodes(self, parentRingNode, divergingNode, nodeIdentifier):
        """
        Generate bifurcation nodes and it children-branch nodes
        :param childNode:
        :param p2cNodeCurve:
        :param nodeIdentifier:
        :return:
        """
        child1NodeCurve = divergingNode.getChildCurve(0)
        nodeIdentifier, child1_firstRingNodes, child1_lastRingNodes = self.create2DtubeNodesAlong1Dtree(divergingNode, child1NodeCurve, nodeIdentifier, isTube=False)

        child2NodeCurve = divergingNode.getChildCurve(1)
        nodeIdentifier, child2_firstRingNodes, child2_lastRingNodes = self.create2DtubeNodesAlong1Dtree(divergingNode, child2NodeCurve, nodeIdentifier, isTube=False)

        rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
            make_tube_bifurcation_points(None, parentRingNode['x'], parentRingNode['d2'],
                                         None, child2_firstRingNodes['x'], child2_firstRingNodes['d2'],
                                         None, child1_firstRingNodes['x'], child1_firstRingNodes['d2'])
        # Store identifiers index for generating elements
        self._bifurcationSide.append(parentRingNode['side'])
        self._paNodeIdentifersList.append(parentRingNode['nodeIdentifiers'])
        self._c1NodeIdentifersList.append(child1_firstRingNodes['nodeIdentifiers'])
        self._c2NodeIdentifersList.append(child2_firstRingNodes['nodeIdentifiers'])
        self._c1NodeIdentifersList_temp.append(child1_firstRingNodes['nodeIdentifiers'])
        self._c2NodeIdentifersList_temp.append(child2_firstRingNodes['nodeIdentifiers'])
        self._paStartIndex.append(paStartIndex)
        self._c1StartIndex.append(c1StartIndex)
        self._c2StartIndex.append(c2StartIndex)
        self._countTubeElements.append(None)

        roNodeIdentifersList_temp = []
        for n in range(len(rox)):
            node = self._nodes.createNode(nodeIdentifier, self._nodetemplate)
            self._cache.setNode(node)
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_VALUE, 1, rox[n])
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
            roNodeIdentifersList_temp.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        self._roNodeIdentifersList.append(roNodeIdentifersList_temp)

        coNodeIdentifersList_temp = []
        for n in range(len(cox)):
            node = self._nodes.createNode(nodeIdentifier, self._nodetemplate)
            self._cache.setNode(node)
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_VALUE, 1, cox[n])
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
            self._coordinates.setNodeParameters(self._cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
            coNodeIdentifersList_temp.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
        self._coNodeIdentifersList.append(coNodeIdentifersList_temp)

        return nodeIdentifier, child1_firstRingNodes, child2_firstRingNodes

    def create2DElementAlong1DTree(self, elementIdentifier):
        """
        Generate elements for both bifurcation and 2D tube elements using depth-first travsersal
        :param elementIdentifier:
        :return:
        """
        useCrossDerivatives = False
        tubeNodeIdentifersList = self._tubeNodeIdentifersList
        paNodeId_temp = self._paNodeIdentifersList
        c1NodeId_temp = self._c1NodeIdentifersList
        c2NodeId_temp = self._c2NodeIdentifersList
        roNodeId_temp = self._roNodeIdentifersList
        coNodeId_temp = self._coNodeIdentifersList
        paStartIndex_temp = self._paStartIndex
        c2StartIndex_temp = self._c1StartIndex
        c1StartIndex_temp = self._c2StartIndex
        bifurcationSide_temp = self._bifurcationSide
        for e3 in range(len(self._countTubeElements)):
            # Bifurcation elements
            if None == self._countTubeElements[e3]:
                paNodeId = paNodeId_temp.pop(0)
                c1NodeId = c1NodeId_temp.pop(0)
                c2NodeId = c2NodeId_temp.pop(0)
                roNodeId = roNodeId_temp.pop(0)
                coNodeId = coNodeId_temp.pop(0)
                paStartIndex = paStartIndex_temp.pop(0)
                c2StartIndex = c1StartIndex_temp.pop(0)
                c1StartIndex = c2StartIndex_temp.pop(0)
                bifurcationSide = bifurcationSide_temp.pop(0)

                # re-indexing the tube elements to the bifurcation elements
                pre_diff = 10  # initialise preset
                for i in range(self._elementsCountAround):
                    select_node = paNodeId[i]
                    root_side = self.FindACrossDerivativesAtNode(select_node, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2)
                    # Calculate smallest differences using squared difference
                    diff = magnitude([bifurcationSide[j] - root_side[j] for j in range(3)])
                    if abs(diff) < abs(pre_diff):
                        startIndex = i
                        pre_diff = diff
                paStartIndex = startIndex if startIndex != 0 else paStartIndex
                elementIdentifier = make_tube_bifurcation_elements_2d(self._region, self._coordinates,
                                                                      elementIdentifier,
                                                                      paNodeId, paStartIndex,
                                                                      c2NodeId, c2StartIndex,
                                                                      c1NodeId, c1StartIndex,
                                                                      roNodeId, coNodeId,
                                                                      useCrossDerivatives)
                continue
            # Tube elements
            for e2 in range(self._countTubeElements[e3] - 1):
                for e1 in range(self._elementsCountAround):
                    # re-indexing the tube elements to the bifurcation elements
                    if e1 == 0:
                        pre_diff = 10
                        select_node = tubeNodeIdentifersList[1][0]
                        bni3_side = self.FindACrossDerivativesAtNode(select_node,
                                                                     Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2)
                        for i in range(self._elementsCountAround):
                            select_node = tubeNodeIdentifersList[0][i]
                            bni1_side = self.FindACrossDerivativesAtNode(select_node,
                                                                         Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2)
                            # Calculate smallest differences using squared difference
                            diff = magnitude([bni3_side[j] - bni1_side[j] for j in range(3)])
                            if abs(diff) < abs(pre_diff):
                                startIndex = i
                                pre_diff = diff
                    if startIndex == 0:
                        bni1 = e1
                        bni2 = (bni1 + 1) % self._elementsCountAround
                    else:
                        bni1 = (e1 + startIndex) % self._elementsCountAround
                        bni2 = (bni1 + 1) % self._elementsCountAround
                    bni3 = e1
                    bni4 = (e1 + 1) % self._elementsCountAround

                    nodeIdentifiers = [tubeNodeIdentifersList[0][bni1], tubeNodeIdentifersList[0][bni2],
                                       tubeNodeIdentifersList[1][bni3], tubeNodeIdentifersList[1][bni4]]
                    element = self._mesh.createElement(elementIdentifier, self._elementtemplateStd)
                    element.setNodesByIdentifier(self._eftStd, nodeIdentifiers)
                    elementIdentifier += 1
                tubeNodeIdentifersList.pop(0)
            tubeNodeIdentifersList.pop(0)

        return elementIdentifier

    def findSide(self, Node):
        """
        Calculate a cross derivatives at each ring during node creation for alignment
        :param Node:
        :return:
        """
        side = [1.0, 0.0, 0.0]
        if (Node._children == []) or (len(Node._children) == 1):
            return side
        _, parent1_d1, _, _, _, _ = Node.getChildCurve(0)
        _, parent2_d1, _, _, _, _ = Node.getChildCurve(1)
        side = normalize(cross(parent1_d1, parent2_d1))
        side = [element * -1 for element in side]

        return side

    def FindACrossDerivativesAtNode(self, nodeIdentifer, d1, d2):
        """
        Find a cross derivatives in the first existing nodes for alignment
        :param nodeIdentifer:
        :param VALUE_LABEL_D_DS_1:
        :param VALUE_LABEL_D_DS_2:
        :return:
        """
        fm = self._fm
        coordinates = self._coordinates
        cache = self._cache
        nodes = self._nodes
        nodetemplate = self._nodetemplate
        nodetemplate.defineField(coordinates)
        node = nodes.findNodeByIdentifier(nodeIdentifer)

        versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE)
        if versionsCount == 1:
            cache.setNode(node)
            result0, d1 = coordinates.getNodeParameters(cache, -1, d1, 1, 3)
            result0, d2 = coordinates.getNodeParameters(cache, -1, d2, 1, 3)
            crossProduct = normalize(cross(d1, d2))

        return crossProduct

def createNodeTreeFromParameters(fieldParameters, generation, generationCount, x1, d1, r, forkNormal):
    """
    reconstructing and restoring a 1D bifurcation tree from the extracted 1D bifurcation scaffold.
    :param fieldParameters:
    :param generation:
    :param generationCount:
    :param x1:
    :param d1:
    :param r:
    :param forkNormal:
    :return:
    """
    node = TreeNode(x1, d1[0], r[0])
    d1_1 = d1[0] if len(d1) == 1 else d1[1]
    r_1 = r[0] if len(r) == 1 else r[1]
    d1_2 = d1[0] if len(d1) == 1 else d1[2]
    r_2 = r[0] if len(r) == 1 else r[2]
    if generation < generationCount:
        node_branch1 = fieldParameters.pop(0)
        x_branch1 = node_branch1[1][0][0]
        d1_branch1 = node_branch1[1][1]
        r_branch1 = node_branch1[1][2][0]
        branch1Normal = normalize(cross(forkNormal, d1_branch1[0]))
        node.addChild(
            createNodeTreeFromParameters(fieldParameters, generation + 1, generationCount, x_branch1, d1_branch1, r_branch1, branch1Normal), d1_1, r_1)

        node_branch2 = fieldParameters.pop(0)
        x_branch2 = node_branch2[1][0][0]
        d1_branch2 = node_branch2[1][1]
        r_branch2 = node_branch2[1][2][0]
        branch2Normal = normalize(cross(forkNormal, d1_branch2[0]))
        node.addChild(
            createNodeTreeFromParameters(fieldParameters, generation + 1, generationCount, x_branch2, d1_branch2, r_branch2, branch2Normal), d1_2, r_2)

    return node
