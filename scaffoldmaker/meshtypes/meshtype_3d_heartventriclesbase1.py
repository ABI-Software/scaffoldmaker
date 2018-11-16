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
        options['Number of elements around atrial free wall'] = 6
        options['Number of elements around atrial septum'] = 3
        # works best with particular numbers of elements up
        options['Number of elements up LV apex'] = 1
        options['Number of elements up RV'] = 4
        # additional options
        options['Atria base inner major axis length'] = 0.55
        options['Atria base inner minor axis length'] = 0.45
        options['Atria major axis rotation degrees'] = 40.0
        options['Atrial septum thickness'] = 0.06
        options['Atrial base wall thickness'] = 0.05
        options['Atrial base slope degrees'] = 30.0
        options['Base height'] = 0.15
        options['Base thickness'] = 0.06
        options['Fibrous ring thickness'] = 0.005
        options['LV outlet front incline degrees'] = 15.0
        options['LV outlet inner diameter'] = 0.3
        options['LV outlet wall thickness'] = 0.025
        options['RV outlet left incline degrees'] = 30.0
        options['RV outlet inner diameter'] = 0.27
        options['RV outlet wall thickness'] = 0.025
        options['Ventricles outlet element length'] = 0.1
        options['Ventricles outlet spacing y'] = 0.02
        options['Ventricles outlet spacing z'] = 0.1
        options['Ventricles rotation degrees'] = 16.0
        options['Ventricles translation x'] = -0.19
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
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        Here the number of elements around LV free wall is dependent on
        number of elements around atrial free wall.
        '''
        dependentChanges = MeshType_3d_heartventricles1.checkOptions(options)
        # only works with particular numbers of elements around
        options['Number of elements around RV free wall'] = 7
        # Supports only 6 or 8 elements around atrial free wall:
        if options['Number of elements around atrial free wall'] <= 6:
            options['Number of elements around atrial free wall'] = 6
            requiredElementsCountAroundLVFreeWall = 5
        else:
            options['Number of elements around atrial free wall'] = 8
            requiredElementsCountAroundLVFreeWall = 7
        if options['Number of elements around LV free wall'] != requiredElementsCountAroundLVFreeWall:
            options['Number of elements around LV free wall'] = requiredElementsCountAroundLVFreeWall
            dependentChanges = True
        #if options['Number of elements around atrial septum'] < 2:
        #    options['Number of elements around atrial septum'] = 2
        options['Number of elements around atrial septum'] = 3
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
        return dependentChanges

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
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountAroundRV = 2*elementsCountAroundRVFreeWall
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

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        # generate heartventricles1 model to add base plane to
        annotationGroups = MeshType_3d_heartventricles1.generateBaseMesh(region, options)

        # find/add annotation groups
        lvGroup = findAnnotationGroupByName(annotationGroups, 'left ventricle')
        rvGroup = findAnnotationGroupByName(annotationGroups, 'right ventricle')
        vSeptumGroup = findAnnotationGroupByName(annotationGroups, 'interventricular septum')
        conusArteriosusGroup = AnnotationGroup(region, 'conus arteriosus', FMANumber = 0, lyphID = 'Lyph ID unknown')
        annotationGroups += [ conusArteriosusGroup ]
        # av boundary nodes are put in left and right fibrous ring groups only so they can be found by heart1
        lFibrousRingGroup = AnnotationGroup(region, 'left fibrous ring', FMANumber = 77124, lyphID = 'Lyph ID unknown')
        rFibrousRingGroup = AnnotationGroup(region, 'right fibrous ring', FMANumber = 77125, lyphID = 'Lyph ID unknown')

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

        # discover ventricles top LV inner, RV inner, V Outer nodes, coordinates and derivatives
        startLVInnerNodeId = 2 + (elementsCountUpLV - 1)*elementsCountAroundLV
        lvInnerNodeId = [ (startLVInnerNodeId + n1) for n1 in range(elementsCountAroundLV) ]
        startRVInnerNodeId = startLVInnerNodeId + elementsCountAroundLV + elementsCountAroundRVFreeWall + 1 + elementsCountAroundRV*(elementsCountUpRV - 1)
        rvInnerNodeId = [ (startRVInnerNodeId + n1) for n1 in range(elementsCountAroundRV) ]
        startVOuterNodeId = startRVInnerNodeId + elementsCountAroundRV + 1 + (elementsCountUpLV - 1)*elementsCountAroundLV
        vOuterNodeId = [ (startVOuterNodeId + n1) for n1 in range(elementsCountAroundLV) ]
        for nodeId in [ lvInnerNodeId, rvInnerNodeId, vOuterNodeId ]:
            vx  = []
            vd1 = [] if (nodeId is rvInnerNodeId) else None
            vd2 = []
            #vd3 = [] if False else None
            for n1 in range(len(nodeId)):
                node = nodes.findNodeByIdentifier(nodeId[n1])
                cache.setNode(node)
                result, x  = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                vx.append(x)
                if vd1 is not None:
                    result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                    vd1.append(d1)
                result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                vd2.append(d2)
                #if vd3:
                #    result, d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                #    vd3.append(d3)
            if nodeId is lvInnerNodeId:
                lvInnerx, lvInnerd2 = vx, vd2
            elif nodeId is rvInnerNodeId:
                rvInnerx, rvInnerd1, rvInnerd2 = vx, vd1, vd2
            else:
                vOuterx, vOuterd2 = vx, vd2

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
        zero = [ 0.0, 0.0, 0.0 ]
        lvOutletOuterd3 = [ None ]*elementsCountAroundOutlet

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

        # fix derivative 3 on lv outlet adjacent to rv outlet
        n1 = elementsCountAroundOutlet//2
        lvOutletOuterd3[n1] = interpolateLagrangeHermiteDerivative(lvOutletOuterx[n1], rvOutletOuterx[0], rvOutletd2, 0.0)

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
        # get av bottom derivative 2 from Hermite-Lagrange interpolation from top row of ventricles
        elementsCountLVFreeWallRegular = elementsCountAroundLVFreeWall - 1
        for n1 in range(elementsCountLVFreeWallRegular + 1):
            noa = elementsCountAroundAtrialFreeWall - elementsCountLVFreeWallRegular + n1
            nov = elementsCountAroundLVFreeWall - elementsCountLVFreeWallRegular + n1
            lavInnerd2[0][noa] = interpolateHermiteLagrangeDerivative(lvInnerx[nov], lvInnerd2[nov], lavInnerx[0][noa], 1.0)
            lavOuterd2[0][noa] = interpolateHermiteLagrangeDerivative(vOuterx[nov], vOuterd2[nov], lavOuterx[0][noa], 1.0)
            if n1 == 0:
                # add d1 to d2 since subtracted in use:
                for c in range(3):
                    lavInnerd2[0][noa][c] += lavInnerd1[0][noa][c]
                    lavOuterd2[0][noa][c] += lavOuterd1[0][noa][c]
        elementsCountRVHanging = 2 if (elementsCountAroundAtrialFreeWall == 8) else 0
        elementsCountRVFreeWallRegular = 3 + elementsCountRVHanging
        for n1 in range(elementsCountRVFreeWallRegular + 1):
            noa = elementsCountAroundAtrialSeptum + n1 - 1
            niv = n1
            nov = elementsCountAroundLVFreeWall + n1
            if elementsCountRVHanging:
                if n1 > 1:
                    niv -= 1
                    nov -= 1
                if n1 > 3:
                    niv -= 1
                    nov -= 1
            if elementsCountRVHanging and ((n1 == 2) or (n1 == 4)):
                six  = interpolateCubicHermite(rvInnerx[niv], rvInnerd1[niv], rvInnerx[niv + 1], rvInnerd1[niv + 1], 0.5)
                sox  = interpolateCubicHermite( vOuterx[nov],  vOuterd2[nov],  vOuterx[nov + 1],  vOuterd2[nov + 1], 0.5)
                sid2 = [ 0.5*(rvInnerd2[niv][c] + rvInnerd2[niv + 1][c]) for c in range(3) ]
                sod2 = [ 0.5*( vOuterd2[nov][c] +  vOuterd2[nov + 1][c]) for c in range(3) ]
            else:
                six, sid2 = rvInnerx[niv], rvInnerd2[niv]
                sox, sod2 = vOuterx[nov], vOuterd2[nov]
            ravInnerd2[0][noa] = interpolateHermiteLagrangeDerivative(six, sid2, ravInnerx[0][noa], 1.0)
            ravOuterd2[0][noa] = interpolateHermiteLagrangeDerivative(sox, sod2, ravOuterx[0][noa], 1.0)
        if elementsCountRVHanging:
            # subtract d1 on last point, later map to d1 + d2
            noa = elementsCountAroundAtrialSeptum + 4
            for c in range(3):
                ravInnerd2[0][noa][c] -= ravInnerd1[0][noa][c]
                ravOuterd2[0][noa][c] -= ravOuterd1[0][noa][c]
        for n1 in range(elementsCountAroundAtrialSeptum + 1):
            noa = (elementsCountAroundAtrialFreeWall + n1)%elementsCountAroundAtria
            nov = elementsCountAroundLVFreeWall + n1
            lavInnerd2[0][noa] = interpolateHermiteLagrangeDerivative(lvInnerx[nov], lvInnerd2[nov], lavInnerx[0][noa], 1.0)
            noa = (elementsCountAroundAtrialSeptum - 1 + elementsCountAroundAtria - n1)%elementsCountAroundAtria
            nov = -n1
            ravInnerd2[0][noa] = interpolateHermiteLagrangeDerivative(rvInnerx[nov], rvInnerd2[nov], ravInnerx[0][noa], 1.0)
        # special fix for left cfb
        lavInnerd2[0][1] = interpolateHermiteLagrangeDerivative(lvOutletInnerx[-1], [ -d for d in lvOutletd2 ], lavInnerx[0][1], 1.0)
        # special fix for left-central cfb (fix from above)
        sd2 = [ (-lvOutletInnerd1[0][c] - lvOutletd2[c]) for c in range(3) ]
        pd2 = smoothCubicHermiteDerivativesLine([ lvOutletInnerx[0], lavInnerx[0][0] ], [ sd2, lavInnerd2[0][0] ], fixStartDerivative=True, fixEndDirection=True)
        lavInnerd2[0][0] = pd2[1]
        # special fix for right cfb
        noa = elementsCountAroundAtria - 2
        nov = -elementsCountAroundAtrialSeptum - 1
        mag = baseHeight + baseThickness
        d2 = vector.setMagnitude(vector.crossproduct3(ravInnerd3[0][noa], ravInnerd1[0][noa]), mag)
        pd2 = smoothCubicHermiteDerivativesLine([ rvInnerx[nov], ravInnerx[0][noa]], [ rvInnerd2[nov], d2 ], fixStartDerivative=True, fixEndDirection=True)
        ravInnerd2[0][noa] = pd2[1]

        # set d2 at ra node mid supraventricular crest to be normal to surface; smooth to get final magnitude later
        ravsvcn1 = elementsCountAroundAtria - 3
        mag = baseHeight + baseThickness
        ravInnerd2[0][ravsvcn1] = vector.setMagnitude(vector.crossproduct3(ravInnerd3[0][ravsvcn1], ravInnerd1[0][ravsvcn1]), mag)
        ravOuterd2[0][ravsvcn1] = vector.setMagnitude(vector.crossproduct3(ravOuterd3[0][ravsvcn1], ravOuterd1[0][ravsvcn1]), mag)
        ravsvcn2 = elementsCountAroundAtria - 4

        # copy derivative 3 from av points to LV outlet at centre, left and right cfb; negate as d1 is reversed:
        lvOutletOuterd3[0] = [ -d for d in lavOuterd3[0][0] ]
        lvOutletOuterd3[1] = [ -d for d in ravOuterd3[0][-2] ]
        lvOutletOuterd3[-1] = [ -d for d in lavOuterd3[0][1] ]

        # create point above anterior ventricular septum end
        xi = 0.5
        fd2 = [ (lvOutletOuterx[4][c] - vOuterx[1][c]) for c in range(3) ]
        mag = 1.0*vector.magnitude(fd2)
        fd2 = vector.setMagnitude([ fd2[0], fd2[1], 0.0 ], mag)
        x = interpolateCubicHermite(vOuterx[0], vOuterd2[0], lvOutletOuterx[4], fd2, xi)
        d2 = interpolateCubicHermiteDerivative(vOuterx[0], vOuterd2[0], lvOutletOuterx[4], fd2, xi)
        pd2 = smoothCubicHermiteDerivativesLine([ vOuterx[0], x, lvOutletOuterx[4] ], [ vOuterd2[0], d2, lvOutletd2 ], fixAllDirections=True, fixStartDerivative=True, fixEndDerivative=True)
        pd1 = smoothCubicHermiteDerivativesLine([ rvOutletOuterx[-1], x, lavOuterx[0][2] ], [ [ -d for d in rvOutletd2 ], zero, lavOuterd2[0][2] ], fixStartDerivative=True, fixEndDirection=True, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
        avsx = x
        avsd1 = pd1[1]
        avsd2 = pd2[1]
        avsd3 = interpolateHermiteLagrangeDerivative(lvInnerx[0], lvInnerd2[0], avsx, 1.0)
        lavOuterd2[0][2] = pd1[2]

        # create points on bottom and top of RV supraventricular crest
        ns = elementsCountAroundRVFreeWall//2 + 1
        nf = elementsCountAroundRVFreeWall + 2
        xis = 0.667
        xif = 1.0 - xis
        mx = [ xis*rvInnerx[ns][0] + xif*rvInnerx[nf][0], xis*rvInnerx[ns][1] + xif*rvInnerx[nf][1], -(fibrousRingThickness + baseThickness) ]
        md2 = [ (rvInnerx[nf][c] - rvInnerx[ns][c]) for c in range(3) ]
        sd1 = [ -d for d in ravInnerd2[0][ravsvcn2] ]
        if elementsCountRVHanging == 0:
            for c in range(3):
                sd1[c] += ravInnerd1[0][ravsvcn2][c]
        fd1 = [ (rvOutletInnerd1[2][c] + rvOutletd2[c]) for c in range(3) ]
        pd1 = smoothCubicHermiteDerivativesLine([ ravInnerx[0][ravsvcn2], mx, rvOutletInnerx[1] ], [ sd1, zero, fd1 ],
            fixStartDerivative=True, fixEndDerivative = True)
        pd2 = smoothCubicHermiteDerivativesLine([ rvInnerx[ns], mx, rvInnerx[nf] ], [ rvInnerd2[ns], md2, [ -d for d in rvInnerd2[nf] ] ],
            fixStartDerivative = True, fixEndDerivative = True)
        svcix = [ mx[0], mx[1], mx[2] ]  # list components to avoid reference bug
        svcid1 = pd1[1]
        svcid2 = pd2[1]
        svcid3 = vector.setMagnitude(vector.crossproduct3(svcid1, svcid2), baseThickness)
        sd2 = [ (svcid2[c] - svcid1[c]) for c in range(3) ]
        pd2 = smoothCubicHermiteDerivativesLine([ mx, ravInnerx[0][ravsvcn1] ], [ sd2, ravInnerd2[0][ravsvcn1] ],
            fixStartDerivative=True, fixEndDirection = True)
        ravInnerd2[0][ravsvcn1] = pd2[1]

        #ravInnerd2[0][ravsvcn1] = [ -d for d in pd1[0] ]

        mx = [ (mx[c] + svcid3[c]) for c in range(3) ]
        md2 = svcid2
        nf = 2
        sd1 = [ -d for d in ravOuterd2[0][ravsvcn2] ]
        if elementsCountRVHanging == 0:
            for c in range(3):
                sd1[c] += ravOuterd1[0][ravsvcn2][c]
        fd1 = [ (rvOutletOuterd1[2][c] + rvOutletd2[c]) for c in range(3) ]
        pd1 = smoothCubicHermiteDerivativesLine([ ravOuterx[0][ravsvcn2], mx, rvOutletOuterx[1] ], [ sd1, zero, fd1 ],
            fixStartDerivative=True, fixEndDerivative = True)
        pd2 = smoothCubicHermiteDerivativesLine([ mx, lvOutletOuterx[nf] ], [ md2, svcid2 ], fixStartDirection=True)
        svcox = [ mx[0], mx[1], mx[2] ]  # list components to avoid reference bug
        svcod1 = pd1[1]
        svcod2 = pd2[0]
        svcod3 = svcid3
        lvOutletOuterd3[nf] = [ -d for d in pd2[1] ]
        sd2 = [ (svcod2[c] - svcod1[c]) for c in range(3) ]
        pd2 = smoothCubicHermiteDerivativesLine([ mx, ravOuterx[0][ravsvcn1] ], [ sd2, ravOuterd2[0][ravsvcn1] ],
            fixStartDerivative=True, fixEndDirection = True)
        ravOuterd2[0][ravsvcn1] = pd2[1]

        # LV outlet nodes
        lvOutletNodeId = [ [], [] ]
        for n3 in range(2):
            if n3 == 0:
                lvOutletx  = lvOutletInnerx
                lvOutletd1 = lvOutletInnerd1
                lvOutletd3 = None
            else:
                lvOutletx  = lvOutletOuterx
                lvOutletd1 = lvOutletOuterd1
                lvOutletd3 = lvOutletOuterd3
            outletNodeId = []
            for n1 in range(elementsCountAroundOutlet):
                node = nodes.createNode(nodeIdentifier, nodetemplate if (lvOutletd3 and lvOutletd3[n1]) else nodetemplateLinearS3)
                lvOutletNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lvOutletx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lvOutletd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lvOutletd2)
                if (lvOutletd3 and lvOutletd3[n1]):
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lvOutletd3[n1])
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
                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                rvOutletNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rvOutletx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rvOutletd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rvOutletd2)
                nodeIdentifier += 1

        # Node above above anterior ventricular septum end
        avsNodeId = nodeIdentifier
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, avsx)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, avsd1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, avsd2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, avsd3)
        nodeIdentifier += 1

        # AV fibrous ring nodes
        lavInnerNodeId = [ [], [] ]
        lavOuterNodeId = [ [], [] ]
        ravInnerNodeId = [ [], [] ]
        ravOuterNodeId = [ [], [] ]
        for n3 in range(2):
            for n2 in range(1): # 2):
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
                            if n2 == 0:
                                # substitute LV outlet node around cfb / fibrous trigones
                                if (i == 0) and (n1 <= 1):
                                    avNodeId.append(lvOutletNodeId[1][0] if (n1 == 0) else lvOutletNodeId[1][-1])
                                    continue
                                elif (i == 1) and (n1 >= (elementsCountAroundAtria - 2)):
                                    avNodeId.append(lvOutletNodeId[1][0] if (n1 == (elementsCountAroundAtria - 1)) else lvOutletNodeId[1][1])
                                    continue
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

        # add nodes to left/right fibrous ring groups so heart1 can find them
        for i in range(2):
            fibrousRingGroup = lFibrousRingGroup if (i == 0) else rFibrousRingGroup
            fibrousRingNodesetGroup = fibrousRingGroup.getNodesetGroup(nodes)
            for n3 in range(2):
                if n3 == 0:
                    avNodeId = lavInnerNodeId[0] if (i == 0) else ravInnerNodeId[0]
                else:
                    avNodeId = lavOuterNodeId[0] if (i == 0) else ravOuterNodeId[0]
                for n1 in range(elementsCountAroundAtria):
                    if avNodeId[n1]:
                        node = nodes.findNodeByIdentifier(avNodeId[n1])
                        fibrousRingNodesetGroup.addNode(node)

        # nodes on bottom and top of RV supraventricular crest
        svciNodeId = nodeIdentifier
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, svcix)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, svcid1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, svcid2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, svcid3)
        nodeIdentifier += 1
        svcoNodeId = nodeIdentifier
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, svcox)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, svcod1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, svcod2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, svcod3)
        nodeIdentifier += 1

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)
        conusArteriosusMeshGroup = conusArteriosusGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()

        elementIdentifier = startElementIdentifier = getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        scalefactors5hanging = [ -1.0, 0.5, 0.25, 0.125, 0.75 ]

        # LV base elements row 1, starting at anterior interventricular sulcus
        for e in range(-1, elementsCountAroundLVFreeWall):
            eft1 = eft
            nids = None
            scalefactors = None
            meshGroups = [ lvMeshGroup ]

            if e == -1:
                # 4 node collapsed tetrahedral element on anterior interventricular sulcus
                nids = [ rvInnerNodeId[elementsCountAroundRVFreeWall], lvInnerNodeId[0], vOuterNodeId[0], avsNodeId ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, []), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 1, 2, 3, 3, 4, 4 ]
                remapEftLocalNodes(eft1, 4, ln_map)
                meshGroups += [ rvMeshGroup ]
            else: # regular elements between LV free wall and left atrium
                noa = elementsCountAroundAtrialFreeWall - elementsCountAroundLVFreeWall + e
                nov = e
                nids = [ lvInnerNodeId[nov], lvInnerNodeId[nov + 1], lavInnerNodeId[0][noa], lavInnerNodeId[0][noa + 1],
                          vOuterNodeId[nov],  vOuterNodeId[nov + 1], lavOuterNodeId[0][noa], lavOuterNodeId[0][noa + 1] ]
                if e == 0:
                    # 7 node collapsed hex-wedge "hedge" element
                    nids.pop(2)
                    nids[5] = avsNodeId
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS2, [])
                    remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, []) ])
                    ln_map = [ 1, 2, 1, 3, 4, 5, 6, 7 ]
                    remapEftLocalNodes(eft1, 7, ln_map)
                elif e == 1:
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                elif e == (elementsCountAroundLVFreeWall - 1):
                    # general linear map d3 adjacent to collapsed crux, transition to atria
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element lv base r1', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # RV base elements row 1, starting at crux / posterior interventricular sulcus
        eInfundibulumStart = elementsCountAroundRVFreeWall + elementsCountRVHanging - 3
        for e in range(-1, elementsCountAroundRVFreeWall + elementsCountRVHanging + 2):
            eft1 = eft
            nids = None
            scalefactors = None
            meshGroups = [ rvMeshGroup ]

            noa = elementsCountAroundAtrialSeptum - 1 + e
            niv = e
            if elementsCountRVHanging:
                if e > 1:
                    niv -= 1
                if e > 3:
                    niv -= 1
            nivp = niv + 1
            nov = elementsCountAroundLVFreeWall + niv
            novp = (nov + 1)%elementsCountAroundLV
            if e == -1:
                # crux / posterior interventricular sulcus, collapsed to 6 node wedge
                nids = [ lvInnerNodeId[elementsCountAroundLVFreeWall], rvInnerNodeId[nivp], lavInnerNodeId[0][elementsCountAroundAtrialFreeWall], ravInnerNodeId[0][noa + 1],
                          vOuterNodeId[novp], ravOuterNodeId[0][noa + 1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ lvMeshGroup ]
            elif e < elementsCountRVFreeWallRegular:
                nids = [ rvInnerNodeId[niv], rvInnerNodeId[nivp], ravInnerNodeId[0][noa], ravInnerNodeId[0][noa + 1],
                          vOuterNodeId[nov],  vOuterNodeId[novp], ravOuterNodeId[0][noa], ravOuterNodeId[0][noa + 1] ]
                if e == 0:
                    # general linear map d3 adjacent to collapsed crux, transition to atria
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                elif elementsCountRVHanging:
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1, 102, 104, 108, 304], [])
                    scalefactors = scalefactors5hanging
                    if e in [ 1, 3 ]:
                        # 1st of pair of elements with hanging nodes at xi1=0.5 on xi2 == 0 plane
                        tricubichermite.setEftMidsideXi1HangingNode(eft1, 2, 1, 1, 2, [1, 2, 3, 4, 5])
                        tricubichermite.setEftMidsideXi1HangingNode(eft1, 6, 5, 5, 6, [1, 2, 3, 4, 5])
                    else:  # e in [ 2, 4 ]:
                        # 2nd of pair of elements with hanging nodes at xi1=0.5 on xi2 == 0 plane
                        tricubichermite.setEftMidsideXi1HangingNode(eft1, 1, 2, 1, 2, [1, 2, 3, 4, 5])
                        tricubichermite.setEftMidsideXi1HangingNode(eft1, 5, 6, 5, 6, [1, 2, 3, 4, 5])
                        if e == 4:
                            remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
            elif e == elementsCountRVFreeWallRegular:
                # supraventricular crest outer 1
                nids = [ rvInnerNodeId[niv], rvInnerNodeId[nivp], ravInnerNodeId[0][noa], svciNodeId,
                          vOuterNodeId[nov],  vOuterNodeId[novp], ravOuterNodeId[0][noa], svcoNodeId ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                if elementsCountRVHanging:
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
            elif e == eInfundibulumStart:
                # supraventricular crest outer 2, outer infundibulum 1
                nids = [ rvInnerNodeId[niv], rvInnerNodeId[nivp], svciNodeId, rvOutletNodeId[0][2],
                          vOuterNodeId[nov],  vOuterNodeId[novp], svcoNodeId, rvOutletNodeId[1][2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == (eInfundibulumStart + 1):
                # supraventricular crest outer 3, outer infundibulum 2
                nids = [ rvInnerNodeId[niv], rvInnerNodeId[nivp], rvOutletNodeId[0][2], rvOutletNodeId[0][3],
                          vOuterNodeId[nov],  vOuterNodeId[novp], rvOutletNodeId[1][2], rvOutletNodeId[1][3] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == (eInfundibulumStart + 2):
                # outer infundibulum 3
                nids = [ rvInnerNodeId[niv], rvInnerNodeId[nivp], rvOutletNodeId[0][3], rvOutletNodeId[0][4],
                          vOuterNodeId[nov],  vOuterNodeId[novp], rvOutletNodeId[1][3], rvOutletNodeId[1][4] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == (eInfundibulumStart + 3):
                # outer infundibulum 4, above septum end
                # 7 node collapsed element above tetrahedral septum end
                nids = [ rvInnerNodeId[niv],            rvOutletNodeId[0][4], rvOutletNodeId[0][5],
                            vOuterNodeId[0], avsNodeId, rvOutletNodeId[1][4], rvOutletNodeId[1][5] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, []), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 1, 2, 3, 4, 5, 6, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == (eInfundibulumStart + 4):
                # outer infundibulum 5, above septum
                nids = [ rvInnerNodeId[niv - 1],   rvInnerNodeId[niv], rvOutletNodeId[0][5], rvOutletNodeId[0][0],
                                      avsNodeId, lvOutletNodeId[1][3], rvOutletNodeId[1][5], rvOutletNodeId[1][0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1]), ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, []), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D2_DS1DS2, []) ])  # swap begin
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS3, []) ])  # swap end
                meshGroups += [ conusArteriosusMeshGroup ]

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element rv base r1', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # interventricular septum elements, row 1
        for e in range(elementsCountAroundVSeptum + 1):
            eft1 = eft
            nids = None
            scalefactors = None
            meshGroups = [ lvMeshGroup, rvMeshGroup, vSeptumMeshGroup ]

            lv1 = elementsCountAroundLVFreeWall + e
            lv2 = (lv1 + 1)%elementsCountAroundLV
            rv1 = -e
            rv2 = rv1 - 1

            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            if e < elementsCountAroundAtrialSeptum:
                la1 = (elementsCountAroundAtrialFreeWall + e)%elementsCountAroundAtria
                la2 = (la1 + 1)%elementsCountAroundAtria
                ra1 = elementsCountAroundAtrialSeptum - 1 - e
                ra2 = ra1 - 1
                nids = [ lvInnerNodeId[lv1], lvInnerNodeId[lv2], lavInnerNodeId[0][la1], lavInnerNodeId[0][la2],
                         rvInnerNodeId[rv1], rvInnerNodeId[rv2], ravInnerNodeId[0][ra1], ravInnerNodeId[0][ra2] ]
                if e == 0:
                    # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                    scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                    scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                elif e == (elementsCountAroundAtrialSeptum - 1):
                    # general linear map d3 adjacent to cfb
                    scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                    scaleEftNodeValueLabels(eft1, [ 5, 6, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                else:
                    scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
            elif e == elementsCountAroundAtrialSeptum:
                # cfb: 5 node inclined pyramid element
                la1 = (elementsCountAroundAtrialFreeWall + e)%elementsCountAroundAtria
                ra1 = elementsCountAroundAtrialSeptum - 1 - e
                nids = [ lvInnerNodeId[lv1], lvOutletNodeId[1][0], lavInnerNodeId[0][la1], rvInnerNodeId[rv1], ravInnerNodeId[0][ra1] ]
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS3, [])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                ln_map = [ 1, 2, 3, 2, 4, 2, 5, 2 ]
                remapEftLocalNodes(eft1, 5, ln_map)
            else:
                # wedge elements along remainder of interventricular septum
                lv1 -= 1
                lv2 -= 1
                rv1 += 1
                rv2 += 1
                lo1 = e - elementsCountAroundAtrialSeptum - 1
                lo2 = lo1 + 1
                nids = [ lvInnerNodeId[lv1], lvInnerNodeId[lv2], lvOutletNodeId[1][lo1], lvOutletNodeId[1][lo2],
                         rvInnerNodeId[rv1], rvInnerNodeId[rv2] ]
                scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS3, [])
                if lo1 == 0:
                    scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                elif lv2 == 0:  # final element on v septum
                    nids[3] = avsNodeId
                    scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                    remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [] ), (Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [] ), (Node.VALUE_LABEL_D_DS3, [] ) ])
                else:
                    scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 3, 4 ]
                remapEftLocalNodes(eft1, 6, ln_map)

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element sp base', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)


        # LV base elements row 2
        for e in range(8):
            eft1 = eft
            nids = None
            scalefactors = None
            meshGroups = [ lvMeshGroup ]

            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            if e == 0:
                # 7 node collapsed element where lv outlet ring expands into
                nids = [    lvInnerNodeId[-1], lvInnerNodeId[0], lvOutletNodeId[0][3], lvOutletNodeId[0][4],
                         lvOutletNodeId[1][3],        avsNodeId, lvOutletNodeId[1][4] ]
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1]) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 5, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 1:
                # 6 node collapsed wedge element mid lv 'crest'
                nids = [ lvInnerNodeId[0], lavInnerNodeId[0][2], lvOutletNodeId[0][4],
                                avsNodeId, lavOuterNodeId[0][2], lvOutletNodeId[1][4] ]
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 3, 3, 4, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 2:
                # 7 node collapsed element where lv outlet ring expands into
                nids = [ lavInnerNodeId[0][2], lavInnerNodeId[0][1], lvOutletNodeId[0][4], lvOutletNodeId[0][5],
                         lavOuterNodeId[0][2], lavOuterNodeId[0][1], lvOutletNodeId[1][4] ]
                remapEftNodeValueLabel(eft1, [ 1, 2, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 7, 6 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 3:
                # 6 node wedge element bridge/curtain between mitral and aortic valve orifices
                no = e - 4
                nids = [ lavInnerNodeId[0][1], lavInnerNodeId[0][0], lvOutletNodeId[0][no], lvOutletNodeId[0][no + 1], lvOutletNodeId[1][no], lvOutletNodeId[1][no + 1] ]
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                # GRC cfb d3 should be reversed in future
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 4:
                # tetrahedral cfb shim-bridge connector element
                ni = elementsCountAroundLVFreeWall + elementsCountAroundAtrialSeptum
                nids = [ lavInnerNodeId[0][0], lvInnerNodeId[ni], lvOutletNodeId[0][0], lvOutletNodeId[1][0] ]
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 3, 3, 4, 4, 4, 4 ]
                remapEftLocalNodes(eft1, 4, ln_map)
            else:
                # 6 node shim wedge elements around LV outlet
                no = e - 5
                ni = elementsCountAroundLVFreeWall + elementsCountAroundAtrialSeptum + no
                nids = [ lvInnerNodeId[ni], lvInnerNodeId[ni + 1], lvOutletNodeId[0][no], lvOutletNodeId[0][no + 1], lvOutletNodeId[1][no], lvOutletNodeId[1][no + 1] ]
                if no == 0:
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                tricubichermite.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element lv base r2', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)


        # RV base elements row 2
        for e in range(5):
            eft1 = eft
            nids = None
            scalefactors = None
            meshGroups = [ rvMeshGroup ]

            if e == 0:
                # 6 node collapsed vs-ra shim element
                rvin1 = -elementsCountAroundAtrialSeptum
                nids = [ lvOutletNodeId[1][0], lvOutletNodeId[1][1], rvInnerNodeId[rvin1], rvInnerNodeId[rvin1 - 1], ravInnerNodeId[0][-1], ravInnerNodeId[0][-2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                scalefactors = [ -1.0 ]
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS2, [])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                ln_map = [ 1, 2, 1, 2, 3, 4, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
            elif e == 1:
                # 7-node collapsed rv crest inner 1, by RA-LV outlet junction
                rvin1 = -elementsCountAroundAtrialSeptum - 1
                nids = [ ravInnerNodeId[0][-3], rvInnerNodeId[rvin1 - 1], ravInnerNodeId[0][-2], rvInnerNodeId[rvin1], ravOuterNodeId[0][-3], lvOutletNodeId[1][2], lvOutletNodeId[1][1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                ln_map = [ 1, 2, 3, 4, 5, 6, 7, 7 ]
                remapEftLocalNodes(eft1, 7, ln_map)
            elif e == 2:
                # 8-node rv crest row 2 element 1
                rvin1 = -elementsCountAroundAtrialSeptum - 2
                nids = [ ravInnerNodeId[0][-4], svciNodeId, ravInnerNodeId[0][-3], rvInnerNodeId[rvin1], ravOuterNodeId[0][-4], svcoNodeId, ravOuterNodeId[0][-3], lvOutletNodeId[1][2] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                if elementsCountRVHanging:
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] )])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
            elif e == 3:
                # 8-node wedge rv crest row 2 element 2
                rvin1 = -elementsCountAroundAtrialSeptum - 2
                nids = [ svciNodeId, rvOutletNodeId[0][2], rvInnerNodeId[rvin1], rvOutletNodeId[0][1], svcoNodeId, rvOutletNodeId[1][2], lvOutletNodeId[1][2], rvOutletNodeId[1][1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            elif e == 4:
                # 8-node rv crest inner 4 by rv outlet
                rvin1 = -elementsCountAroundAtrialSeptum - 2
                nids = [ rvInnerNodeId[rvin1], rvOutletNodeId[0][1], rvInnerNodeId[rvin1 - 1], rvOutletNodeId[0][0], lvOutletNodeId[1][2], rvOutletNodeId[1][1], lvOutletNodeId[1][3], rvOutletNodeId[1][0] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                tricubichermite.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                tricubichermite.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])  # must do before following
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary, to swap with D_DS2
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                meshGroups += [ conusArteriosusMeshGroup ]
            else:
                continue

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element rv base r2', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)


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
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        elementsCountAroundLV = elementsCountAroundLVFreeWall + elementsCountAroundRVFreeWall
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountAroundRV = 2*elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        elementsCountAroundAtrialFreeWall = options['Number of elements around atrial free wall']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundAtria = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        MeshType_3d_heartventricles1.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startBaseLvElementIdentifier = element.getIdentifier()
        startBaseRvElementIdentifier = startBaseLvElementIdentifier + elementsCountAroundLVFreeWall + 1
        elementsCountRVHanging = 2 if (elementsCountAroundAtrialFreeWall == 8) else 0
        startBaseSeptumElementIdentifier = startBaseRvElementIdentifier + elementsCountAroundRVFreeWall + elementsCountRVHanging + 3
        startBaseLv2ElementIdentifier = startBaseSeptumElementIdentifier + elementsCountAroundVSeptum + 1
        startBaseRv2ElementIdentifier = startBaseLv2ElementIdentifier + 8
        limitBaseElementIdentifier = startBaseRv2ElementIdentifier + 5
        #print(startBaseLvElementIdentifier, startBaseRvElementIdentifier, startBaseSeptumElementIdentifier, limitBaseElementIdentifier)
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            numberInXi3 = None
            elementId = element.getIdentifier()
            if elementId < startBaseRvElementIdentifier:
                # LV row 1
                if elementId == startBaseLvElementIdentifier:
                    # collapsed element on anterior interventricular sulcus:
                    numberInXi1 = numberInXi3 = refineElementsCountThroughLVWall
                else:
                    numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < startBaseSeptumElementIdentifier:
                # RV row 1
                if elementId == startBaseRvElementIdentifier:
                    # collapsed element on posterior interventricular sulcus:
                    numberInXi1 = numberInXi3 = refineElementsCountThroughLVWall
                else:
                    numberInXi3 = refineElementsCountThroughWall
            elif elementId < startBaseLv2ElementIdentifier:
                # V septum
                numberInXi3 = refineElementsCountThroughLVWall
            elif elementId < startBaseRv2ElementIdentifier:
                # LV row 2
                numberInXi3 = refineElementsCountThroughLVWall
            else:
                # RV row 2
                numberInXi3 = refineElementsCountThroughWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementId == (limitBaseElementIdentifier - 1):
                return  # finish on last so can continue in full heart mesh
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
