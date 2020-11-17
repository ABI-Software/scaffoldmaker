'''
Generates 3D lung surface mesh, with variable numbers of elements around, along.
'''

import copy
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.geometry import createEllipsePoints, createEllipsoidPoints
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues, mesh_destroy_elements_and_nodes_by_identifiers
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.node import Node


class MeshType_3d_lungs1(Scaffold_base):
    '''
    3D lung scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Lungs 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1',
            'Human 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            'Number of elements around': 8,
            'Number of elements along': 4,
            'Number of elements through wall': 1,
            'Major diameter': 30.0,
            'Minor diameter': 15.0,
            'Height': 70,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4
        }
        if 'Mouse' in parameterSetName:
            options['Major diameter'] = 20.0
            options['Minor diameter'] = 10.0
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Number of elements around',
            'Number of elements along',
            'Major diameter',
            'Minor diameter',
            'Height',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along']
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
    def checkOptions(cls, options):
        for key in [
            'Number of elements around',
            'Number of elements along',
            'Refine number of elements along',
            'Refine number of elements around']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''

        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        else:
            nodetemplate = nodetemplateApex


        cache = fm.createFieldcache()

        elementsCount1 = 2
        elementsCount2 = 4
        elementsCount3 = 4

        # Create nodes
        nodesList = [[[-1.0, -1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, -0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 1.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],

                     [[-1.0, -1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, -0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 1.0, 1.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],

                     [[-1.0, -1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, -0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.5, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],

                     [[-1.0, -1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, -0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.5, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 1.0, 3.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],

                     [[-1.0, -1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, -0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, -0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, -0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 0.5, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-1.0, 1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[-0.5, 1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     [[0.0, 1.0, 4.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]]]
        print('len(nodesList)', len(nodesList))

        nodeIdentifier = 1
        lNodeIds = []
        nix = 0
        # d1 = [ 1.0, 0.0, 0.0 ]
        for n3 in range(elementsCount3 + 1):
            lNodeIds.append([])
            for n2 in range(elementsCount2 + 1):
                lNodeIds[n3].append([])
                for n1 in range(elementsCount1 + 1):
                    lNodeIds[n3][n2].append(None)
                    # if n3 < elementsCount3:
                    #     if (n1 == 0) and ((n2 == 0) or (n2 == elementsCount2)):
                    #         continue
                    # else:
                    #     if (n2 == 0) or (n2 == elementsCount2) or (n1 == 0):
                    #         continue
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    # x = [ n1*1.0, n2*1.0, n3*1.0 ]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nodesList[nix][0])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nodesList[nix][1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nodesList[nix][2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nodesList[nix][3])
                    lNodeIds[n3][n2][n1] = nodeIdentifier
                    nix += 1
                    nodeIdentifier += 1

        # cx = []
        # cd1 = []
        # cd2 = []
        # cd3 = []
        # for n1 in range(len(nodesList)):
        #     x = nodesList[n1][0]
        #     dx_ds1 = nodesList[n1][1]
        #     dx_ds2 = nodesList[n1][2]
        #     dx_ds3 = nodesList[n1][3]
        #     cx.append(x)
        #     cd1.append(dx_ds1)
        #     cd2.append(dx_ds2)
        #     cd3.append(dx_ds3)
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3[n1])
        #     nodeIdentifier += 1
        #

        # Create elements
        from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite

        mesh = fm.findMeshByDimension(3)
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eftRegular = eftfactory.createEftBasic()

        elementtemplateRegular = mesh.createElementtemplate()
        elementtemplateRegular.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateRegular.defineField(coordinates, -1, eftRegular)

        elementtemplateCustom = mesh.createElementtemplate()
        elementtemplateCustom.setElementShapeType(Element.SHAPE_TYPE_CUBE)

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
                        lNodeIds[e3 + 1][e2][e1], lNodeIds[e3 + 1][e2][e1 + 1], lNodeIds[e3 + 1][e2 + 1][e1], lNodeIds[e3 + 1][e2 + 1][e1 + 1] ]
                    scalefactors = None
                    if (e3 < elementsCount3 - 1):
                        if (e2 == 0) and (e1 == 0):
                            # Back wedge elements
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eft1
                            scalefactors = [ -1.0 ]
                        elif (e2 == (elementsCount2 - 1)) and (e1 == 0):
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
                            scalefactors = [ -1.0 ]
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
                            # Top back tetrahedron elements
                            nodeIdentifiers.pop(6)
                            nodeIdentifiers.pop(5)
                            nodeIdentifiers.pop(4)
                            nodeIdentifiers.pop(0)
                            eft = eft6
                            scalefactors = [ -1.0 ]
                        if (e2 == elementsCount2 - 1) and (e1 == 0):
                            # Top front tetrahedron elements
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

        annotationGroups = []
        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        '''
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        '''
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountAround = options['Refine number of elements around']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong)
        return



