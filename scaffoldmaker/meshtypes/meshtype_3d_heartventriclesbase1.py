"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valve regions.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_heartventricles1 import MeshType_3d_heartventricles1
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
import scaffoldmaker.utils.vector as vector
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventriclesbase1(object):
    '''
    Generates a 3-D heart ventricles with base plane model, ready to attach the
    atria, mitral and tricuspid valves, with LV + RV outlets ready to attach
    aorta and pulmonary trunk and their valve regions.
    '''

    @staticmethod
    def getName():
        return '3D Heart Ventricles with Base 1'

    @staticmethod
    def getDefaultOptions():
        options = MeshType_3d_heartventricles1.getDefaultOptions()
        # only works with particular numbers of elements around
        options['Number of elements around LV free wall'] = 5
        options['Number of elements around RV free wall'] = 7
        options['Number of elements around atrial free wall'] = 8
        options['Number of elements around atrial septum'] = 2
        # works best with particular numbers of elements up
        options['Number of elements up LV apex'] = 1
        options['Number of elements up RV'] = 4
        # additional options
        options['Atria base inner major axis length'] = 0.55
        options['Atria base inner minor axis length'] = 0.42
        options['Atria major axis rotation degrees'] = 40.0
        options['Atrial septum thickness'] = 0.06
        options['Atrial base wall thickness'] = 0.05
        options['Atrial base slope degrees'] = 30.0
        options['Base height'] = 0.12
        options['Base thickness'] = 0.06
        options['Fibrous ring thickness'] = 0.01
        options['LV outlet front incline degrees'] = 30.0
        options['LV outlet inner diameter'] = 0.3
        options['LV outlet wall thickness'] = 0.025
        options['RV outlet left incline degrees'] = 30.0
        options['RV outlet inner diameter'] = 0.27
        options['RV outlet wall thickness'] = 0.025
        options['Ventricles outlet element length'] = 0.1
        options['Ventricles outlet spacing'] = 0.04
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventricles1.getOrderedOptionNames()
        optionNames.insert(4, 'Number of elements around atrial free wall')
        optionNames.insert(5, 'Number of elements around atrial septum')
        optionNames += [
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atrial septum thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Base height',
            'Base thickness',
            'Fibrous ring thickness',
            'LV outlet front incline degrees',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet left incline degrees',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Ventricles outlet element length',
            'Ventricles outlet spacing']
        # want refinement options last
        for optionName in [
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_heartventricles1.checkOptions(options)
        # only works with particular numbers of elements around
        #options['Number of elements around LV free wall'] = 5
        #options['Number of elements around RV free wall'] = 7
        #options['Number of elements around atrial free wall'] = 8
        #options['Number of elements around atrial septum'] = 2
        # while editing, restrict to limitations of atria:
        if options['Number of elements around atrial free wall'] < 6:
            options['Number of elements around atrial free wall'] = 6
        # need even number of elements around free wall
        if (options['Number of elements around atrial free wall'] % 2) == 1:
            options['Number of elements around atrial free wall'] += 1
        if options['Number of elements around atrial septum'] < 2:
            options['Number of elements around atrial septum'] = 2
        for key in [
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atrial septum thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Base height',
            'Base thickness',
            'Fibrous ring thickness',
            'LV outlet inner diameter',
            'LV outlet wall thickness',
            'RV outlet inner diameter',
            'RV outlet wall thickness',
            'Ventricles outlet element length',
            'Ventricles outlet spacing']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Atria major axis rotation degrees'] < -75.0:
            options['Atria major axis rotation degrees'] = -75.0
        elif options['Atria major axis rotation degrees'] > 75.0:
            options['Atria major axis rotation degrees'] = 75.0

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
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        elementsCountAroundAtrialFreeWall = options['Number of elements around atrial free wall']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundAtria = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
        # from heartventricles1:
        ivSulcusDerivativeFactor = options['Interventricular sulcus derivative factor']
        lvOuterHeight = options['LV outer height']
        lvOuterRadius = options['LV outer radius']
        lvFreeWallThickness = options['LV free wall thickness']
        #lvApexThickness = options['LV apex thickness']
        #lvInnerHeight = lvOuterHeight - lvApexThickness
        lvInnerRadius = lvOuterRadius - lvFreeWallThickness
        #rvInnerHeightFraction = options['RV inner height fraction']
        rvArcAroundBaseRadians = math.radians(options['RV arc around degrees'])
        #rvArcApexFraction = options['RV arc apex fraction']
        rvFreeWallThickness = options['RV free wall thickness']
        rvWidth = options['RV width']
        #rvWidthGrowthFactor = options['RV width growth factor']
        rvSideExtension = options['RV side extension']
        #rvSideExtensionGrowthFactor = options['RV side extension growth factor']
        vSeptumThickness = options['Ventricular septum thickness']
        vSeptumBaseRadialDisplacement = options['Ventricular septum base radial displacement']
        # from heartatria1:
        aBaseInnerMajorMag = 0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = 0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aSeptumThickness = options['Atrial septum thickness']
        aBaseWallThickness = options['Atrial base wall thickness']
        aBaseSlopeRadians = math.radians(options['Atrial base slope degrees'])
        # new:
        baseHeight = options['Base height']
        baseThickness = options['Base thickness']
        fibrousRingThickness = options['Fibrous ring thickness']
        lvOutletFrontInclineRadians = math.radians(options['LV outlet front incline degrees'])
        lvOutletInnerRadius = 0.5*options['LV outlet inner diameter']
        lvOutletWallThickness = options['LV outlet wall thickness']
        lvOutletOuterRadius = lvOutletInnerRadius + lvOutletWallThickness
        rvOutletLeftInclineRadians = math.radians(options['RV outlet left incline degrees'])
        rvOutletInnerRadius = 0.5*options['RV outlet inner diameter']
        rvOutletWallThickness = options['RV outlet wall thickness']
        rvOutletOuterRadius = rvOutletInnerRadius + rvOutletWallThickness
        vOutletElementLength = options['Ventricles outlet element length']
        vOutletSpacing = options['Ventricles outlet spacing']
        useCrossDerivatives = False

        # generate heartventricles1 model to add base plane to
        annotationGroups = MeshType_3d_heartventricles1.generateBaseMesh(region, options)
        lvGroup = findAnnotationGroupByName(annotationGroups, 'left ventricle')
        rvGroup = findAnnotationGroupByName(annotationGroups, 'right ventricle')
        vSeptumGroup = findAnnotationGroupByName(annotationGroups, 'interventricular septum')
        conusArteriosusGroup = AnnotationGroup(region, 'conus arteriosus', FMANumber = 0, lyphID = 'Lyph ID unknown')
        lFibrousRingGroup = AnnotationGroup(region, 'left fibrous ring', FMANumber = 77124, lyphID = 'Lyph ID unknown')
        rFibrousRingGroup = AnnotationGroup(region, 'right fibrous ring', FMANumber = 77125, lyphID = 'Lyph ID unknown')
        annotationGroups += [ conusArteriosusGroup, lFibrousRingGroup, rFibrousRingGroup ]

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

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
        # LV/RV outlet elements are linear through the wall, hence their nodes do not have D_DS3 parameters
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = startNodeIdentifier = getMaximumNodeIdentifier(nodes) + 1

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)
        conusArteriosusMeshGroup = conusArteriosusGroup.getMeshGroup(mesh)
        lFibrousRingMeshGroup = lFibrousRingGroup.getMeshGroup(mesh)
        rFibrousRingMeshGroup = rFibrousRingGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        MeshType_3d_heartventricles1.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startBaseElementIdentifier = element.getIdentifier()

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
