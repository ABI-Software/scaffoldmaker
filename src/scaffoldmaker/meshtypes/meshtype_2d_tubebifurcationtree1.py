"""
Generates a 2-D tube bifurcation tree mesh.
"""

from __future__ import division
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.lung_terms import get_lung_term
from opencmiss.utils.maths.vectorops import cross, normalize
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import get_curve_circle_points, \
    make_tube_bifurcation_points, make_tube_bifurcation_elements_2d, move_1D_central_path, \
    readRegionFromFile
from scaffoldmaker.utils.interpolation import getCubicHermiteArcLength
import math

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
        if options['Number of elements around root'] < 4:
            options['Number of elements around root'] = 4
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

        # Bifurcation tree extraction
        Readfile = True
        if Readfile:
            filename = "D:\\Python\\venv_sparc\\mapclient_workflow\\scaffold\\mesh_6.exf"
            generationCount, bifurcationTree = readRegionFromFile(filename)
            fieldParameters = None
        else:
            fieldParameters = move_1D_central_path(region, bifurcationTreeScaffold)
            bifurcationTree = bifurcationTreeScaffold.getScaffoldType(). \
                generateBifurcationTree(bifurcationTreeScaffold.getScaffoldSettings())
            generationCount = bifurcationTreeScaffold.getScaffoldSettings()['Number of generations']

        # Tube element template
        mesh = fm.findMeshByDimension(2)
        elementtemplateStd = mesh.createElementtemplate()
        elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eftStd = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eftStd.setFunctionNumberOfTerms(n * 4 + 4, 0)
        elementtemplateStd.defineField(coordinates, -1, eftStd)

        ##############
        # Create nodes
        ##############
        child1 = []
        cx = [[], [], []]
        cd1 = [[], [], []]
        cd2 = [[], [], []]
        paStartIndex = []
        c1StartIndex = []
        c2StartIndex = []
        paNodeId = []
        roNodeId = []
        coNodeId = []
        c1NodeId = []
        c2NodeId = []
        stemNodeId = []
        lastNodeId = []
        missingRight = []
        missingLeft = []
        noRootNode = []

        # Generate bifurcation tree nodes and elements and return a next nodeIdentifier
        # nodeIdentifier = bifurcationTree.generateZincModel(region)[0]
        nodeIdentifier = 1

        # Create circle nodes around the 1D bifurcation tree
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        # Get a bifurcation root node
        rootNode = bifurcationTree._rootNode
        if fieldParameters != None:
            # Update central path - root node
            rootNode._x = fieldParameters[0][1][0][0]
            rootNode._d1 = fieldParameters[0][1][1]
        centralPathIndex = []

        # Create a branch at each generation
        for j in range(generationCount - 1):
            x = []
            d1 = []
            d2 = []
            paStartIndex.append([])
            c1StartIndex.append([])
            c2StartIndex.append([])
            child1.append([])
            centralPathIndex.append([])
            paNodeId.append([])
            roNodeId.append([])
            coNodeId.append([])
            c1NodeId.append([])
            c2NodeId.append([])
            stemNodeId.append([])
            missingRight.append([])
            missingLeft.append([])

            # Create a bifurcation at each branch
            for k in range(2 ** j):
                paStartIndex[j].append([])
                c1StartIndex[j].append([])
                c2StartIndex[j].append([])
                x.append([[], [], []])
                d1.append([[], [], []])
                d2.append([[], [], []])
                child1[j].append(None)
                curveLength = [[], [], []]
                centralPathIndex[j].append([])
                lastNodeId.append([]) if j == (generationCount - 2) else None
                missingRight[j].append(None)
                missingLeft[j].append(None)
                stemNodeId[j].append([])
                paNodeId[j].append([])
                c1NodeId[j].append([])
                c2NodeId[j].append([])

                if j == 0:
                    try:
                        child1[j][k] = rootNode.getChild(1)
                        child1[j][k] = rootNode
                        noRootNode = True
                    except:
                        child1[j][k] = rootNode.getChild(0)
                        noRootNode = False
                    if fieldParameters != None:
                        # Update central path - Parent nodes
                        centralPathIndex[j][k] = 1
                        child1[j][k]._x = fieldParameters[centralPathIndex[j][k]][1][0][0]
                        child1[j][k]._d1 = fieldParameters[centralPathIndex[j][k]][1][1]
                        # Update central path - Children nodes
                        nextCentralPathIndex = [centralPathIndex[j][k] + 1, centralPathIndex[j][k] + (2 ** (generationCount - 1))]
                        child1[j][k].getChild(0)._x = fieldParameters[nextCentralPathIndex[0]][1][0][0]
                        child1[j][k].getChild(0)._d1 = fieldParameters[nextCentralPathIndex[0]][1][1]
                        child1[j][k].getChild(1)._x = fieldParameters[nextCentralPathIndex[1]][1][0][0]
                        child1[j][k].getChild(1)._d1 = fieldParameters[nextCentralPathIndex[1]][1][1]
                    # Parental branch
                    xParent1, xParentd1, rParent1, xParent2, xParentd2, rParent2 = rootNode.getChildCurve(0)
                    if noRootNode:
                        xParent2 = xParent1
                        xParentd2 = xParentd1
                        rParent2 = rParent1
                    # Get the right(0) [left(1)] children curve from the bifurcation tree
                    try:
                        xRight1, xRightd1, rRight1, xRight2, xRightd2, rRight2 = child1[j][k].getChildCurve(0)
                    except:
                        missingRight[j][k] = True
                    try:
                        xLeft1, xLeftd1, rLeft1, xLeft2, xLeftd2, rLeft2 = child1[j][k].getChildCurve(1)
                    except:
                        missingLeft[j][k] = True
                else:
                    if (k % 2) == 0:
                        try:
                            child1[j][k] = child1[j-1][k//2].getChild(0)
                        except:
                            continue
                        if fieldParameters != None:
                            # Update central path - Children nodes
                            centralPathIndex[j][k] = centralPathIndex[j - 1][k // 2] + 1
                            nextCentralPathIndex = [centralPathIndex[j][k] + 1,
                                                    centralPathIndex[j][k] + (2 ** (generationCount - j - 1))]
                            child1[j][k].getChild(0)._x = fieldParameters[nextCentralPathIndex[0]][1][0][0]
                            child1[j][k].getChild(0)._d1 = fieldParameters[nextCentralPathIndex[0]][1][1]
                            child1[j][k].getChild(1)._x = fieldParameters[nextCentralPathIndex[1]][1][0][0]
                            child1[j][k].getChild(1)._d1 = fieldParameters[nextCentralPathIndex[1]][1][1]
                    else:
                        try:
                            child1[j][k] = child1[j-1][k//2].getChild(1)
                        except:
                            continue
                        if fieldParameters != None:
                            # Update central path - Children nodes
                            centralPathIndex[j][k] = centralPathIndex[j - 1][k // 2] + (2 ** (generationCount - j))
                            nextCentralPathIndex = [centralPathIndex[j][k] + 1,
                                                    centralPathIndex[j][k] + (2 ** (generationCount - j - 1))]
                            child1[j][k].getChild(0)._x = fieldParameters[nextCentralPathIndex[0]][1][0][0]
                            child1[j][k].getChild(0)._d1 = fieldParameters[nextCentralPathIndex[0]][1][1]
                            child1[j][k].getChild(1)._x = fieldParameters[nextCentralPathIndex[1]][1][0][0]
                            child1[j][k].getChild(1)._d1 = fieldParameters[nextCentralPathIndex[1]][1][1]

                    # Parental branch
                    xParent1, xParentd1, rParent1, xParent2, xParentd2, rParent2 = child1[j-1][k//2].getChildCurve(k % 2)
                    # Get the right(0) [left(1)] children curve from the bifurcation tree
                    try:
                        xRight1, xRightd1, rRight1, xRight2, xRightd2, rRight2 = child1[j][k].getChildCurve(0)
                    except:
                        missingRight[j][k] = True
                    try:
                        xLeft1, xLeftd1, rLeft1, xLeft2, xLeftd2, rLeft2 = child1[j][k].getChildCurve(1)
                    except:
                        missingLeft[j][k] = True

                # Set up the curve axis for bifurcation
                if missingRight[j][k] or missingLeft[j][k]:
                    side = [1.0, 0.0, 0.0]
                else:
                    axis1 = normalize(xRightd2)
                    axis2 = normalize(xLeftd2)
                    side = cross(axis1, axis2)
                    side = [element * -1 for element in side]

                # Check if there are children branches
                for l in range(3):
                    if l == 0:
                        x1  = xParent1
                        xd1 = xParentd1
                        x2  = xParent2
                        xd2 = xParentd2
                        r1  = rParent1
                        r2  = rParent2
                        rd1 = rd2 = (rParent2 - rParent1)
                    else:
                        if (l == 1) and missingRight[j][k]:
                            continue
                        elif (l == 2) and missingLeft[j][k]:
                            continue
                        x1  = xRight1   if l == 1 else xLeft1
                        xd1 = xRightd1  if l == 1 else xLeftd1
                        r1  = rRight1   if l == 1 else rLeft1
                        x2  = xRight2   if l == 1 else xLeft2
                        xd2 = xRightd2  if l == 1 else xLeftd2
                        r2  = rRight1   if l == 1 else rLeft1
                        rd1 = rd2 = (rRight2 - rRight1) if l == 1 else (rLeft2 - rLeft1) # The radius of the circle between the parental nodes
                        if missingRight[j][k] or missingLeft[j][k]:
                            x1Junction = xParent2
                            xd1Junction = xParentd2
                            r1Junction = rParent2
                            x2Junction = x2
                            xd2Junction = xd2
                            r2Junction = r2
                            rd1Junction = rd2Junction = (r2Junction - r1Junction)

                    # Curve length between the parental nodes and divided into 2 bifurcation points
                    curveLength[l] = getCubicHermiteArcLength(x1, xd1, x2, xd2)
                    elementsCount = max(2, math.ceil(curveLength[l] / maxElementLength))
                    elementLength = curveLength[l] / elementsCount

                    # Get ring nodes along the curve
                    if l == 0:
                        skip = False
                        elementsCountEnd = elementsCount if (j == 0) or (missingRight[j][k] and missingLeft[j][k]) else elementsCount - 2
                        elementsCountEnd = 1 if elementsCountEnd == 0 else elementsCountEnd
                        for m in range(elementsCountEnd):
                            x[k][l].append(None)
                            d1[k][l].append(None)
                            d2[k][l].append(None)
                            xi = m/elementsCount if (j == 0) else (m + 2)/elementsCount
                            if elementsCountEnd == 1:
                                xi = (m + 1)/elementsCount
                                skip = True
                            # Parent branch
                            cx[l], cd1[l], cd2[l], x[k][l][m], d1[k][l][m], d2[k][l][m] = get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side, elementsCountAroundRoot)
                    else:
                        elementsCountEnd = 1 if j != (generationCount - 2) else elementsCount
                        for m in range(elementsCountEnd):
                            x[k][l].append(None)
                            d1[k][l].append(None)
                            d2[k][l].append(None)
                            xi = (m + 1)/elementsCount
                            # Children branch
                            cx[l], cd1[l], cd2[l], x[k][l][m], d1[k][l][m], d2[k][l][m] = get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side, elementsCountAroundRoot)

                if missingRight[j][k] and missingLeft[j][k]:
                    pass
                elif missingRight[j][k] or missingLeft[j][k]:
                    for m in range(1):
                        xi = m / elementsCount
                        cxTemp, cd1Temp, cd2Temp, xTemp, d1Temp, d2Temp = get_curve_circle_points(x1Junction, xd1Junction, x2Junction, xd2Junction, r1Junction, rd1Junction,
                                                                                                  r2Junction, rd2Junction, xi, elementLength, side, elementsCountAroundRoot)
                        x[k][0].append(xTemp)
                        d1[k][0].append(d1Temp)
                        d2[k][0].append(d2Temp)
                else:
                    rox, rod1, rod2, cox, cod1, cod2, paStartIndex[j][k], c1StartIndex[j][k], c2StartIndex[j][k] = \
                        make_tube_bifurcation_points(cx[0], x[k][0][-1], d2[k][0][-1], cx[2], x[k][2][0], d2[k][2][0], cx[1], x[k][1][0], d2[k][1][0])

                for l in range(3):
                    # Add children index to stem index
                    if (l == 0) and (j > 0) and (skip is not True):
                        stemNodeId[j][k] = c1NodeId[j-1][k // 2].copy() if (k % 2) == 0 else c2NodeId[j-1][k // 2].copy()
                    for m in range(len(x[k][l])):
                        # Skip replicated nodes at each child branch
                        if skip:
                            paNodeId[j][k] = c1NodeId[j-1][k // 2].copy() if (k % 2) == 0 else c2NodeId[j-1][k // 2].copy()
                            skip = False
                            continue
                        # Skip replicated nodes at each parent (Only for noRootNodes)
                        if noRootNode and (m <= (len(x[k][l]) - 2)):
                            if (x[k][0][m][0] == x[k][0][m-1][0]):
                                x[k][0][m+1] = x[k][0][m]
                                d1[k][0][m + 1] = d1[k][0][m]
                                d2[k][0][m + 1] = d2[k][0][m]
                                continue
                        for n in range(elementsCountAroundRoot):
                            node = nodes.createNode(nodeIdentifier, nodetemplate)
                            cache.setNode(node)
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[k][l][m][n])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[k][l][m][n])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[k][l][m][n])
                            if (l == 0):
                                # Store parent nodes in the last ring at each branch
                                if (m == len(x[k][l]) - 1):
                                    paNodeId[j][k].append(nodeIdentifier)
                                stemNodeId[j][k].append(nodeIdentifier)
                            else:
                                if missingRight[j][k] and missingLeft[j][k]:
                                    continue
                                # Store nodes in the first ring at each child branch
                                if (m == 0):
                                    c1NodeId[j][k].append(nodeIdentifier) if l == 1 else c2NodeId[j][k].append(nodeIdentifier)
                                # Store stemNodeIndex at each end
                                if (j == generationCount - 2):
                                    stemNodeId[j][k].append(nodeIdentifier)
                                    # Store lastNodeIndex at each end
                                    if (m == len(x[k][l]) - 1):
                                        lastNodeId[k].append(nodeIdentifier)
                            nodeIdentifier = nodeIdentifier + 1

                    # Bifurcation nodes
                    roNodeId[j].append([])
                    coNodeId[j].append([])
                    if l == 0 and (missingLeft[j][k] is not True) and (missingRight[j][k] is not True):
                        # Skip rox nodes at each bifurcation (Only for noRootNodes)
                        if noRootNode and (x[k][0][0][0] == x[k][0][1][0]):
                            pass
                        else:
                            for n in range(len(rox)):
                                node = nodes.createNode(nodeIdentifier, nodetemplate)
                                cache.setNode(node)
                                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox [n])
                                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
                                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
                                roNodeId[j][k].append(nodeIdentifier)
                                nodeIdentifier = nodeIdentifier + 1

                        for n in range(len(cox)):
                            node = nodes.createNode(nodeIdentifier, nodetemplate)
                            cache.setNode(node)
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox [n])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
                            coNodeId[j][k].append(nodeIdentifier)
                            nodeIdentifier = nodeIdentifier + 1

        ################
        # Create elements
        ################
        elementIdentifier = 1
        for e4 in range(generationCount-1):
            for e3 in range(2 ** e4):
                elementsCountAlong = int(len(stemNodeId[e4][e3])/elementsCountAroundRoot) - 1
                for e2 in range(elementsCountAlong):
                    for e1 in range(elementsCountAroundRoot):
                        # Connect root
                        bni1 = e2 * elementsCountAroundRoot + e1
                        if (bni1 + 1) % (elementsCountAroundRoot * (e2 + 1)) != 0:
                            bni2 = (bni1 + 1) % (elementsCountAroundRoot * (e2 + 1))
                        else:
                            bni2 = ((bni1 + 1) % (elementsCountAroundRoot * (e2 + 1))) + (elementsCountAroundRoot * e2)
                        bni3 = bni1 + elementsCountAroundRoot
                        bni4 = bni2 + elementsCountAroundRoot
                        nodeIdentifiers = [stemNodeId[e4][e3][bni1], stemNodeId[e4][e3][bni2],
                                           stemNodeId[e4][e3][bni3], stemNodeId[e4][e3][bni4]]
                        # print(nodeIdentifiers)
                        if ((nodeIdentifiers[0] and nodeIdentifiers[1]) in paNodeId[e4][e3]) or \
                                ((nodeIdentifiers[0] and nodeIdentifiers[1]) in lastNodeId[e3]):
                            continue

                        element = mesh.createElement(elementIdentifier, elementtemplateStd)
                        element.setNodesByIdentifier(eftStd, nodeIdentifiers)
                        elementIdentifier += 1

                # Bifurcation/Junction elements
                try:
                    if noRootNode:
                        noRootNode = False
                        elementIdentifier = make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier,
                                                                              [], paStartIndex[e4][e3],
                                                                              c2NodeId[e4][e3],
                                                                              c2StartIndex[e4][e3], c1NodeId[e4][e3],
                                                                              c1StartIndex[e4][e3],
                                                                              paNodeId[e4][e3], coNodeId[e4][e3],
                                                                              useCrossDerivatives)
                        print('success')
                    else:
                        # Bifurcation for no root node
                        elementIdentifier = make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier,
                                                                          paNodeId[e4][e3], paStartIndex[e4][e3], c2NodeId[e4][e3],
                                                                          c2StartIndex[e4][e3], c1NodeId[e4][e3], c1StartIndex[e4][e3],
                                                                          roNodeId[e4][e3], coNodeId[e4][e3],
                                                                          useCrossDerivatives)
                except:
                    if missingRight[e4][e3] and missingLeft[e4][e3]:
                        continue
                    else:
                        cNodeId = paNodeId[e4][e3].copy() + c1NodeId[e4][e3].copy() if c1NodeId != [] \
                            else paNodeId[e4][e3].copy() + c2NodeId[e4][e3].copy()
                        elementsCountAlong = int(len(cNodeId) / elementsCountAroundRoot) - 1
                        for e2 in range(elementsCountAlong):
                            for e1 in range(elementsCountAroundRoot):
                                bni1 = e2 * elementsCountAroundRoot + e1
                                if (bni1 + 1) % (elementsCountAroundRoot * (e2 + 1)) != 0:
                                    bni2 = (bni1 + 1) % (elementsCountAroundRoot * (e2 + 1))
                                else:
                                    bni2 = ((bni1 + 1) % (elementsCountAroundRoot * (e2 + 1))) + (
                                                elementsCountAroundRoot * e2)
                                bni3 = bni1 + elementsCountAroundRoot
                                bni4 = bni2 + elementsCountAroundRoot
                                if (e2 == 0) and (elementsCountAlong > 1):
                                    bni3 = ((bni1 + elementsCountAroundRoot + 2) % elementsCountAroundRoot) + elementsCountAroundRoot
                                    bni4 = ((bni2 + elementsCountAroundRoot + 2) % elementsCountAroundRoot) + elementsCountAroundRoot
                                nodeIdentifiers = [cNodeId[bni1], cNodeId[bni2],
                                                   cNodeId[bni3], cNodeId[bni4]]
                                element = mesh.createElement(elementIdentifier, elementtemplateStd)
                                element.setNodesByIdentifier(eftStd, nodeIdentifiers)
                                elementIdentifier += 1

        return []
