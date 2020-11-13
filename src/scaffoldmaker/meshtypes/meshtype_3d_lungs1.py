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


    # @classmethod
    # def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
    #     if optionName == 'Central path LUT':
    #         return list(cls.centralPathDefaultScaffoldPackages_LUT.keys())
    #     assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
    #         'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
    #     return scaffoldType.getParameterSetNames()

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
                     [[0.0, 1.0, 2.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]]]
                     #
                     # [[0.0, -1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, -1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, -1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, -0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, -0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, -0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 0.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 0.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 0.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 0.5, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 1.0, 2.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     #
                     # [[0.0, -1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, -1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, -1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, -0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, -0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, -0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 0.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 0.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 0.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 0.5, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[0.0, 1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-0.5, 1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]],
                     # [[-1.0, 1.0, 3.0], [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]]]

                    # [[0.0, 0.0, 4.7], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.5, 0.0]]]

        print('len(nodesList)', len(nodesList))

        cx = []
        cd1 = []
        cd2 = []
        cd3 = []
        nodeIdentifier = 1
        for n1 in range(len(nodesList)):
            x = nodesList[n1][0]
            dx_ds1 = nodesList[n1][1]
            dx_ds2 = nodesList[n1][2]
            dx_ds3 = nodesList[n1][3]
            cx.append(x)
            cd1.append(dx_ds1)
            cd2.append(dx_ds2)
            cd3.append(dx_ds3)
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3[n1])
            nodeIdentifier += 1

        print('len(cx)', len(cx))
        print('len(cd1)', len(cd1))
        print('len(cd2)', len(cd2))


        # Create elements
        from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite

        mesh = fm.findMeshByDimension(3)
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        # elementtemplate2 = mesh.createElementtemplate()
        # elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        # eft2 = eftfactory.createEftWedgeXi1Zero()
        # elementtemplate2.defineField(coordinates, -1, eft2)

        elementIdentifier = 1
        elementsCountAlongMajorDiameter = 4
        elementsCountAlongMinorDiameter = 2
        elementsCountAlong = 2
        #
        nel = (elementsCountAlongMinorDiameter + 1) * (elementsCountAlongMajorDiameter + 1)
        for e3 in range(elementsCountAlong):
            for e2 in range(elementsCountAlongMajorDiameter):
                # First row of elements #4
                for e1 in range(elementsCountAlongMinorDiameter):
                    if (e2 == 0 or e2 == elementsCountAlongMajorDiameter - 1) and e1 == 0:
                        # nodeIdentifiers = [10, 11, 14, 25, 26, 29]
                        # element = mesh.createElement(elementIdentifier, elementtemplate2)
                        # element.setNodesByIdentifier(eft2, nodeIdentifiers)
                        pass
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni1 = e3 * nel + e2 * (elementsCountAlongMinorDiameter + 1) + e1 + 1
                        bni2 = e3 * nel + e2 * (elementsCountAlongMinorDiameter + 1) + (e1 + 1) % (elementsCountAlongMinorDiameter + 1) + 1
                        bni3 = bni1 + elementsCountAlongMinorDiameter + 1
                        bni4 = bni3 + 1
                        bni5 = bni1 + nel
                        bni6 = bni5 + 1
                        bni7 = bni5 + elementsCountAlongMinorDiameter + 1
                        bni8 = bni7 + 1
                        nodeIdentifiers = [bni1, bni2, bni3, bni4,
                                           bni5, bni6, bni7, bni8]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
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

    # def createEftWedgeXi1Zero(self):
    #     '''
    #     Create a basic tricubic hermite element template for elements
    #     along boundary of tenia coli where nodes on xi1 = 0 are collapsed.
    #     :return: Element field template
    #     '''
    #     eft = self.createEftBasic()
    #
    #     # remap parameters on xi1 = 1 before collapsing nodes
    #     remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [])
    #     remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D2_DS1DS2, [])
    #     remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D2_DS1DS2, [])
    #     remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D3_DS1DS2DS3, [])
    #
    #     ln_map = [1, 2, 1, 3, 4, 5, 4, 6]
    #     remapEftLocalNodes(eft, 6, ln_map)
    #     assert eft.validate(), 'eftfactory_tricubichermite.createEftWedgeXi1Zero:  Failed to validate eft'
    #     return eft

