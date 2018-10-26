"""
Generates a 3-D heart ventricles with base plane model, ready to attach the
atria, mitral and tricuspid valves, with LV + RV outlets ready to attach aorta and
pulmonary trunk and their valve regions.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import getLeftAtriumBasePoints
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
        options['Atrial base slope degrees'] = 15.0
        options['Base height'] = 0.14
        options['Base thickness'] = 0.06
        options['Fibrous ring thickness'] = 0.01
        options['LV outlet front incline degrees'] = 15.0
        options['LV outlet inner diameter'] = 0.3
        options['LV outlet wall thickness'] = 0.025
        options['RV outlet left incline degrees'] = 25.0
        options['RV outlet inner diameter'] = 0.27
        options['RV outlet wall thickness'] = 0.025
        options['Ventricles outlet element length'] = 0.1
        options['Ventricles outlet spacing y'] = 0.04
        options['Ventricles outlet spacing z'] = 0.14
        options['Ventricles rotation degrees'] = 20.0
        options['Ventricles translation x'] = -0.22
        options['Ventricles translation y'] = -0.2
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
            'Ventricles outlet spacing y',
            'Ventricles outlet spacing z',
            'Ventricles rotation degrees',
            'Ventricles translation x',
            'Ventricles translation y']
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
            'Ventricles outlet spacing y',
            'Ventricles outlet spacing z']:
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
        vOutletSpacingy = options['Ventricles outlet spacing y']
        vOutletSpacingz = options['Ventricles outlet spacing z']
        vRotationRadians = math.radians(options['Ventricles rotation degrees'])
        vTranslationx = options['Ventricles translation x']
        vTranslationy = options['Ventricles translation y']
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

        # move ventricles to fit atria centred around aorta
        cosVRotationRadians = math.cos(-vRotationRadians)
        sinVRotationRadians = math.sin(-vRotationRadians)
        rotationMatrix = fm.createFieldConstant([ cosVRotationRadians, sinVRotationRadians, 0.0, -sinVRotationRadians, cosVRotationRadians, 0.0, 0.0, 0.0, 1.0 ])
        # offset makes top of fibrous ring (excepting cfb rotation) on z == 0
        ventriclesOffset = fm.createFieldConstant([ vTranslationx, vTranslationy, -(fibrousRingThickness + baseHeight + baseThickness)])
        newCoordinates = fm.createFieldAdd(fm.createFieldMatrixMultiply(3, rotationMatrix, coordinates), ventriclesOffset)
        fieldassignment = coordinates.createFieldassignment(newCoordinates)
        fieldassignment.assign()
        fieldassignment = None
        newCoordinates = None
        ventriclesOffset = None

        # LV outlet points
        elementsCountAroundOutlet = 6
        cosLvOutletFrontInclineRadians = math.cos(lvOutletFrontInclineRadians)
        sinLvOutletFrontInclineRadians = math.sin(lvOutletFrontInclineRadians)
        lvOutletCentre = [ 0.0, 0.0, -fibrousRingThickness ]
        axis1 = [ 0.0, -cosLvOutletFrontInclineRadians, sinLvOutletFrontInclineRadians ]
        axis2 = [ 1.0, 0.0, 0.0 ]
        lvOutletInnerx, lvOutletInnerd1 = createCirclePoints(lvOutletCentre,
            vector.setMagnitude(axis1, lvOutletInnerRadius), vector.setMagnitude(axis2, lvOutletInnerRadius), elementsCountAroundOutlet)
        lvOutletOuterx, lvOutletOuterd1 = createCirclePoints(lvOutletCentre,
            vector.setMagnitude(axis1, lvOutletOuterRadius), vector.setMagnitude(axis2, lvOutletOuterRadius), elementsCountAroundOutlet)
        lvOutletd2 = [ 0.0, vOutletElementLength*sinLvOutletFrontInclineRadians, vOutletElementLength*cosLvOutletFrontInclineRadians ]

        # RV outlet points
        cosRvOutletLeftInclineRadians = math.cos(rvOutletLeftInclineRadians)
        sinRvOutletLeftInclineRadians = math.sin(rvOutletLeftInclineRadians)
        dy = vOutletSpacingy + rvOutletOuterRadius
        dz = vOutletSpacingz
        axis1 = [ 0.0, -cosLvOutletFrontInclineRadians, sinLvOutletFrontInclineRadians ]
        axis2 = [ cosRvOutletLeftInclineRadians, sinRvOutletLeftInclineRadians*sinLvOutletFrontInclineRadians, sinRvOutletLeftInclineRadians*cosLvOutletFrontInclineRadians ]
        axis3 = vector.crossproduct3(axis1, axis2)
        rvOutletCentre = [ (lvOutletOuterx[3][c] - dy*axis1[c] + dz*axis3[c]) for c in range(3) ]
        rvOutletInnerx, rvOutletInnerd1 = createCirclePoints(rvOutletCentre,
            vector.setMagnitude(axis1, rvOutletInnerRadius), vector.setMagnitude(axis2, rvOutletInnerRadius), elementsCountAroundOutlet)
        rvOutletOuterx, rvOutletOuterd1 = createCirclePoints(rvOutletCentre,
            vector.setMagnitude(axis1, rvOutletOuterRadius), vector.setMagnitude(axis2, rvOutletOuterRadius), elementsCountAroundOutlet)
        rvOutletd2 = [ vOutletElementLength*axis3[c] for c in range(3) ]

        # Left av fibrous ring points
        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)
        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)
        aortaOuterRadius = lvOutletOuterRadius
        aBaseFrontInclineRadians = lvOutletFrontInclineRadians
        laCentre, laSeptumRadians, laBaseInnerx, laBaseInnerd1, laBaseInnerd2, laBaseOuterx, laBaseOuterd1, laBaseOuterd2 = \
            getLeftAtriumBasePoints(elementsCountAroundAtrialFreeWall, elementsCountAroundAtrialSeptum,
                aBaseInnerMajorMag, aBaseInnerMinorMag, aMajorAxisRadians,
                aBaseWallThickness, aBaseSlopeHeight, aBaseSlopeLength, aSeptumThickness,
                aortaOuterRadius, aBaseFrontInclineRadians, aBaseSideInclineRadians = 0.0, aBaseBackInclineRadians = 0.0)
        # get d3 from inner-outer difference
        laBaseInnerd3 = []
        laBaseOuterd3 = []
        for n1 in range(elementsCountAroundAtria):
            if laBaseOuterx[n1]:
                innerd3 = outerd3 = [ (laBaseOuterx[n1][c] - laBaseInnerx[n1][c]) for c in range(3) ]
            else:  # septum
                innerd3 = [ -2.0*laBaseInnerx[n1][0], 0.0, 0.0 ]
                outerd3 = None
            laBaseInnerd3.append(innerd3)
            laBaseOuterd3.append(outerd3)
        # fix cfb, crux centre derivative 3, code duplicated from atria1:
        for n1 in [ 0, elementsCountAroundAtrialFreeWall ]:
            laBaseOuterd3[n1] = [ 0.0, laBaseOuterx[n1][1] - laBaseInnerx[n1][1], laBaseOuterx[n1][2] - laBaseInnerx[n1][2] ]

        # displace to get bottom points
        zero = [ 0.0, 0.0, 0.0 ]
        lavInnerx  = [ [], laBaseInnerx ]
        for n1 in range(elementsCountAroundAtria):
            lavInnerx[0].append([ laBaseInnerx[n1][0], laBaseInnerx[n1][1], laBaseInnerx[n1][2] - fibrousRingThickness ])
        lavInnerd1 = [ laBaseInnerd1, laBaseInnerd1 ]
        lavInnerd2 = [ [ zero ]*elementsCountAroundAtria, [ zero ]*elementsCountAroundAtria ]
        lavInnerd3 = [ laBaseInnerd3, laBaseInnerd3 ]
        lavOuterx  = [ [], laBaseOuterx ]
        for n1 in range(elementsCountAroundAtria):
            lavOuterx[0].append([ laBaseOuterx[n1][0], laBaseOuterx[n1][1], laBaseOuterx[n1][2] - fibrousRingThickness ] if (laBaseOuterx[n1]) else None)
        lavOuterd1 = [ laBaseOuterd1, laBaseOuterd1 ]
        lavOuterd2 = [ [ zero ]*elementsCountAroundAtria, [ zero ]*elementsCountAroundAtria ]
        lavOuterd3 = [ laBaseOuterd3, laBaseOuterd3 ]
        # Right av fibrous ring points = mirror
        ravInnerx  = [ [], [] ]
        ravInnerd1 = [ [], [] ]
        ravInnerd2 = [ [], [] ]
        ravInnerd3 = [ [], [] ]
        ravOuterx  = [ [], [] ]
        ravOuterd1 = [ [], [] ]
        ravOuterd2 = [ [], [] ]
        ravOuterd3 = [ [], [] ]
        for n2 in range(2):
            for n1 in range(elementsCountAroundAtria):
                n1l = elementsCountAroundAtria - 1 - n1
                lavx  = lavInnerx [n2][n1l]
                lavd1 = lavInnerd1[n2][n1l]
                lavd2 = lavInnerd2[n2][n1l]
                lavd3 = lavInnerd3[n2][n1l]
                ravInnerx [n2].append([ -lavx [0],  lavx [1],   lavx [2] ])
                ravInnerd1[n2].append([  lavd1[0], -lavd1[1],  -lavd1[2] ])
                ravInnerd2[n2].append([ -lavd2[0],  lavd2[1],   lavd2[2] ])
                ravInnerd3[n2].append([ -lavd3[0],  lavd3[1],   lavd3[2] ])
                if lavOuterx[n2][n1l]:
                    lavx  = lavOuterx [n2][n1l]
                    lavd1 = lavOuterd1[n2][n1l]
                    lavd2 = lavOuterd2[n2][n1l]
                    lavd3 = lavOuterd3[n2][n1l]
                    ravOuterx [n2].append([ -lavx [0],  lavx [1],   lavx [2] ])
                    ravOuterd1[n2].append([  lavd1[0], -lavd1[1],  -lavd1[2] ])
                    ravOuterd2[n2].append([ -lavd2[0],  lavd2[1],   lavd2[2] ])
                    ravOuterd3[n2].append([ -lavd3[0],  lavd3[1],   lavd3[2] ])
                else:
                    # no outer points around septum
                    ravOuterx [n2].append(None)
                    ravOuterd1[n2].append(None)
                    ravOuterd2[n2].append(None)
                    ravOuterd3[n2].append(None)

        # LV outlet nodes
        lvOutletNodeId = [ [], [] ]
        for n3 in range(2):
            if n3 == 0:
                lvOutletx  = lvOutletInnerx
                lvOutletd1 = lvOutletInnerd1
            else:
                lvOutletx  = lvOutletOuterx
                lvOutletd1 = lvOutletOuterd1
            outletNodeId = []
            for n1 in range(elementsCountAroundOutlet):
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)  # if (n3 == 0) else nodetemplate)
                lvOutletNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lvOutletx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lvOutletd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lvOutletd2)
                nodeIdentifier += 1

        # RV outlet nodes
        rvOutletNodeId = [ [], [] ]
        for n3 in range(2):
            if n3 == 0:
                rvOutletx  = rvOutletInnerx
                rvOutletd1 = rvOutletInnerd1
            else:
                rvOutletx  = rvOutletOuterx
                rvOutletd1 = rvOutletOuterd1
            outletNodeId = []
            for n1 in range(elementsCountAroundOutlet):
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)  # if (n3 == 0) else nodetemplate)
                rvOutletNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rvOutletx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rvOutletd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rvOutletd2)
                nodeIdentifier += 1

        # AV fibrous ring nodes
        lavInnerNodeId = [ [], [] ]
        lavOuterNodeId = [ [], [] ]
        ravInnerNodeId = [ [], [] ]
        ravOuterNodeId = [ [], [] ]
        for n3 in range(2):
            for n2 in range(2):
                for i in range(2):
                    if n3 == 0:
                        if i == 0:
                            avx, avd1, avd2, avd3 = lavInnerx[n2], lavInnerd1[n2], lavInnerd2[n2], lavInnerd3[n2]
                            avNodeId = lavInnerNodeId[n2]
                        else:
                            avx, avd1, avd2, avd3 = ravInnerx[n2], ravInnerd1[n2], ravInnerd2[n2], ravInnerd3[n2]
                            avNodeId = ravInnerNodeId[n2]
                    else:
                        if i == 0:
                            avx, avd1, avd2, avd3 = lavOuterx[n2], lavOuterd1[n2], lavOuterd2[n2], lavOuterd3[n2]
                            avNodeId = lavOuterNodeId[n2]
                        else:
                            avx, avd1, avd2, avd3 = ravOuterx[n2], ravOuterd1[n2], ravOuterd2[n2], ravOuterd3[n2]
                            avNodeId = ravOuterNodeId[n2]
                    for n1 in range(elementsCountAroundAtria):
                        if avx[n1] is None:
                            avNodeId.append(None)
                            continue
                        if n3 == 1:
                            if (n2 == 0) and (((i == 0) and (n1 <= 1)) or ((i == 1) and (n1 >= (elementsCountAroundAtria - 2)))):
                                # substitute LV outlet node around cfb / fibrous trigones
                                #avNodeId.append(lvOutletNodeId[1][0] if (n1 == 0) else (lvOutletNodeId[1][-1] if (n1 == 1) else lvOutletNodeId[1][1]))
                                #continue
                                pass
                            if (i == 1) and ((n1 == (elementsCountAroundAtrialSeptum - 1)) or (n1 == (elementsCountAroundAtria - 1))):
                                # find common nodes on right at cfb and crux
                                avNodeId.append(lavOuterNodeId[n2][elementsCountAroundAtria - 1 - n1])
                                continue
                        avNodeId.append(nodeIdentifier)
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, avx [n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, avd1[n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, avd2[n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, avd3[n1])
                        nodeIdentifier += 1

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
