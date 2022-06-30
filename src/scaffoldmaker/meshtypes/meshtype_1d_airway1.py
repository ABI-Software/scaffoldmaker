'''
Generates 3D lung surface mesh.
'''

import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_bifurcationtree1 import MeshType_1d_bifurcationtree1
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_1d_airway1(Scaffold_base):
    '''
    Generates a 1-D airway from 1D_bifurcationtree1
    '''

    @staticmethod
    def getName():
        return '1D Airway 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        bifurcationTreeParameters = ScaffoldPackage(MeshType_1d_bifurcationtree1).getScaffoldSettings()
        options = {}

        for key, value in bifurcationTreeParameters.items():
            options[key] = value

        return options


    @staticmethod
    def getOrderedOptionNames():
        bifurcationTreeParameters = ScaffoldPackage(MeshType_1d_bifurcationtree1).getScaffoldSettings()
        options = list(bifurcationTreeParameters.keys())

        return options

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''

        bifurcationTreeScaffold = ScaffoldPackage(MeshType_1d_bifurcationtree1).getScaffoldType()
        bifurcationTreeParameters = options

        bifurcationTree = bifurcationTreeScaffold.generateBifurcationTree(bifurcationTreeParameters)
        nodeIdentifier = bifurcationTree.generateZincModel(region)[0]

        fm = region.getFieldmodule()
        mesh = fm.findMeshByDimension(1)

        # Annotation groups
        airwayGroup = AnnotationGroup(region, get_lung_term("respiratory airway"))
        airwayMeshGroup = airwayGroup.getMeshGroup(mesh)
        annotationGroups = [airwayGroup]

        # Annotating groups
        for elementIdentifier in range(1, mesh.getSize() + 1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            airwayMeshGroup.addElement(element)

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
