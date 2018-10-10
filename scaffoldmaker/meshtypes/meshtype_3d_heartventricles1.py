"""
Generates 3-D mesh of left and right ventricles below base plane.
Variant using collapsed/wedge elements at septum junction.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
import scaffoldmaker.utils.vector as vector
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_heartventricles1(object):
    '''
    Generates 3-D mesh of left and right ventricles below base plane.
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around LV free wall' : 5,
            'Number of elements around RV free wall' : 7,
            'Number of elements up LV apex' : 1,
            'Number of elements up RV' : 4,
            'Interventricular sulcus derivative factor' : 0.5,
            'LV outer height' : 1.0,
            'LV outer radius' : 0.5,
            'LV free wall thickness' : 0.12,
            'LV apex thickness' : 0.06,
            'RV height fraction' : 0.9,
            'RV arc around degrees' : 200.0,
            'RV arc apex fraction' : 0.5,
            'RV free wall thickness' : 0.04,
            'RV width' : 0.4,
            'RV extra cross radius base' : 0.15,
            'Ventricular septum thickness' : 0.1,
            'Ventricular septum base radial displacement' : 0.15,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Refine number of elements through LV wall' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            'Number of elements up LV apex',
            'Number of elements up RV',
            'Interventricular sulcus derivative factor',
            'LV outer height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV height fraction',
            'RV arc around degrees',
            'RV arc apex fraction',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Ventricular septum thickness',
            'Ventricular septum base radial displacement',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall',
            'Number of elements up LV apex']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around LV free wall',
            'Number of elements up RV']:
            if options[key] < 2:
                options[key] = 2
        for key in [
            'Number of elements around RV free wall']:
            if options[key] < 3:
                options[key] = 3
        for key in [
            'LV outer height',
            'LV outer radius',
            'LV free wall thickness',
            'LV apex thickness',
            'RV height fraction',
            'RV free wall thickness',
            'RV width',
            'RV extra cross radius base',
            'Ventricular septum thickness',
            'Ventricular septum base radial displacement']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['RV height fraction'] > 0.99:
            options['RV height fraction'] = 0.99
        if options['RV arc around degrees'] < 45.0:
            options['RV arc around degrees'] = 45.0
        elif options['RV arc around degrees'] > 270.0:
            options['RV arc around degrees'] = 270.0
        for key in [
            'Interventricular sulcus derivative factor',
            'RV arc apex fraction']:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 1.0:
                options[key] = 1.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        pointsCountAroundOuter = elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        ivSulcusDerivativeFactor = options['Interventricular sulcus derivative factor']
        lvOuterHeight = options['LV outer height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        lvApexThickness = options['LV apex thickness']
        rvHeightFraction = options['RV height fraction']
        rvArcAroundBaseRadians = math.radians(options['RV arc around degrees'])
        rvArcApexFraction = options['RV arc apex fraction']
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        rvExtraCrossRadiusBase = options['RV extra cross radius base']
        vSeptumThickness = options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = options['Ventricular septum base radial displacement']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        lvGroup = AnnotationGroup(region, 'left ventricle', FMANumber = 7101, lyphID = 'Lyph ID unknown')
        rvGroup = AnnotationGroup(region, 'right ventricle', FMANumber = 7098, lyphID = 'Lyph ID unknown')
        vSeptumGroup = AnnotationGroup(region, 'interventricular septum', FMANumber = 7133, lyphID = 'Lyph ID unknown')
        annotationGroups = [ lvGroup, rvGroup, vSeptumGroup ]

        #################
        # Create nodes
        #################

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplateApex = nodetemplate

        nodeIdentifier = 1

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        Stops at end of ventricles, hence can be called from ventriclesbase.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through wall']

        startRvElementIdentifier = elementsCountAroundLV*elementsCountUpLV + 1
        limitRvElementIdentifier = startRvElementIdentifier + elementsCountUpRV*elementsCountAroundRVFreeWall - 2

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            elementId = element.getIdentifier()
            if elementId < startRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < limitRvElementIdentifier:
                numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitRvElementIdentifier - 1):
                return  # finish on last so can continue in ventriclesbase
            element = meshrefinement._sourceElementiterator.next()

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()
