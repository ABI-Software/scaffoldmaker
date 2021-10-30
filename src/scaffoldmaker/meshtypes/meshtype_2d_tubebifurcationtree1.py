"""
Generates a 2-D tube bifurcation mesh.
"""

from __future__ import division
import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import get_curve_circle_points
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

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        ##############
        # Create nodes
        ##############

        bifurcationTree = bifurcationTreeScaffold.getScaffoldType().generateBifurcationTree(bifurcationTreeScaffold.getScaffoldSettings())
        nodeIdentifier = bifurcationTree.generateZincModel(region)[0]

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        rootNode = bifurcationTree.getRootNode()
        child1 = rootNode.getChild(0)

        x1, xd1, r1, x2, xd2, r2 = child1.getChildCurve(0)
        rd1 = rd2 = (r2 - r1)

        curveLength = getCubicHermiteArcLength(x1, xd1, x2, xd2)
        elementsCount = max(2, math.ceil(curveLength/maxElementLength))
        elementLength = curveLength/elementsCount
        side = [ 1.0, 0.0, 0.0 ]
        for i in range(elementsCount + 1):
            xi = i/elementsCount
            x, d1, d2 = get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, elementLength, side, elementsCountAroundRoot)
            for n in range(elementsCountAroundRoot):
                #print('node',nodeIdentifier,i,n,x[n],d1[n],d2[n])
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[n])
                nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################

        elementIdentifier = 1

        return []
