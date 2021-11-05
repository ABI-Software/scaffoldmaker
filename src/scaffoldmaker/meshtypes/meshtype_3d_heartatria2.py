"""
Generates a 3-D heart atria model, suitable for attachment to the
3-D Heart Ventricles with Base 2.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getEllipseArcLength, getEllipseRadiansToX, updateEllipseAngleByArcLength
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import computeNodeDerivativeHermiteLagrange, interpolateNodesCubicHermite


class MeshType_3d_heartatria2(Scaffold_base):
    '''
    3-D heart atria model, suitable for attachment to the 3-D Heart Ventricles with Base 2.
    '''

    @staticmethod
    def getName():
        return '3D Heart Atria 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around atria' : 8,
            'Number of elements around atrial septum' : 2,
            'Number of elements up atria' : 4,
            'Atria base inner major axis length' : 0.55,
            'Atria base inner minor axis length' : 0.42,
            'Atria major axis rotation degrees' : 40.0,
            'Atria outer height' : 0.4,
            'Atria outer major arc up degrees' : 100.0,
            'Atria outer minor arc up degrees' : 125.0,
            'Atrial septum thickness' : 0.06,
            'Atrial free wall thickness' : 0.02,
            'Atrial base wall thickness' : 0.05,
            'Atrial base slope degrees' : 30.0,
            'LV outlet outer diameter' : 0.35,
            'Left pulmonary vein inner diameter' : 0.11,
            'Left pulmonary vein wall thickness' : 0.009,
            'Right pulmonary vein inner diameter' : 0.12,
            'Right pulmonary vein wall thickness' : 0.009,
            'Inferior vena cava inner diameter' : 0.22,
            'Inferior vena cava wall thickness' : 0.015,
            'Superior vena cava inner diameter' : 0.2,
            'Superior vena cava wall thickness' : 0.015,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Refine number of elements through atrial wall' : 1,
            'Use cross derivatives' : False,
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around atria',
            'Number of elements around atrial septum',
            'Number of elements up atria',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atria outer height',
            'Atria outer major arc up degrees',
            'Atria outer minor arc up degrees',
            'Atrial septum thickness',
            'Atrial free wall thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'LV outlet outer diameter',
            'Left pulmonary vein inner diameter',
            'Left pulmonary vein wall thickness',
            'Right pulmonary vein inner diameter',
            'Right pulmonary vein wall thickness', 
            'Inferior vena cava inner diameter',
            'Inferior vena cava wall thickness',
            'Superior vena cava inner diameter',
            'Superior vena cava wall thickness',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through atrial wall',
            #,'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        Here the number of elements around atrial free wall must be an even number,
        but the parameter is number of elements around atria, which depends on
        number of elements around atrial septum.
        '''
        dependentChanges = False
        if options['Number of elements around atria'] < 6:
            options['Number of elements around atria'] = 6
        if options['Number of elements around atrial septum'] < 1:
            options['Number of elements around atrial septum'] = 1
        elif options['Number of elements around atrial septum'] > (options['Number of elements around atria'] - 4):
            options['Number of elements around atrial septum'] = options['Number of elements around atria'] - 4
        # need even number of elements around free wall
        if ((options['Number of elements around atria'] - options['Number of elements around atrial septum']) % 2) == 1:
            options['Number of elements around atria'] += 1
            dependentChanges = True
        if options['Number of elements up atria'] < 3:
            options['Number of elements up atria'] = 3
        for key in [
            'Atria outer major arc up degrees',
            'Atria outer minor arc up degrees']:
            if options[key] < 10.0:
                options[key] = 10.0
            elif options[key] > 170.0:
                options[key] = 170.0
        for key in [
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria outer height',
            'Atrial septum thickness',
            'Atrial free wall thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Left pulmonary vein inner diameter',
            'Left pulmonary vein wall thickness',
            'Right pulmonary vein inner diameter',
            'Right pulmonary vein wall thickness',
            'Inferior vena cava inner diameter',
            'Inferior vena cava wall thickness',
            'Superior vena cava inner diameter',
            'Superior vena cava wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['LV outlet outer diameter'] < options['Atrial septum thickness']:
            options['LV outlet outer diameter'] =options['Atrial septum thickness']
        for key in [
            'Atria major axis rotation degrees']:
            if options[key] < -75.0:
                options[key] = -75.0
            elif options[key] > 75.0:
                options[key] = 75.0
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through atrial wall']:
            if options[key] < 1:
                options[key] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundAtria = options['Number of elements around atria']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountUpAtria = options['Number of elements up atria']
        aBaseInnerMajorMag = 0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = 0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aOuterMajorArcUpRadians = math.radians(options['Atria outer major arc up degrees'])
        aOuterMinorArcUpRadians = math.radians(options['Atria outer minor arc up degrees'])
        aOuterHeight = options['Atria outer height']
        lvOutletOuterRadius = 0.5*options['LV outlet outer diameter']
        aSeptumThickness = options['Atrial septum thickness']
        aFreeWallThickness = options['Atrial free wall thickness']
        aBaseWallThickness = options['Atrial base wall thickness']
        aBaseSlopeRadians = math.radians(options['Atrial base slope degrees'])
        lpvInnerRadius = 0.5*options['Left pulmonary vein inner diameter']
        lpvWallThickness = options['Left pulmonary vein wall thickness']
        rpvInnerRadius = 0.5*options['Right pulmonary vein inner diameter']
        rpvWallThickness = options['Right pulmonary vein wall thickness']
        ivcInnerRadius = 0.5*options['Inferior vena cava inner diameter']
        ivcWallThickness = options['Inferior vena cava wall thickness']
        svcInnerRadius = 0.5*options['Superior vena cava inner diameter']
        svcWallThickness = options['Superior vena cava wall thickness']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        laGroup = AnnotationGroup(region, get_heart_term("left atrium myocardium"))
        raGroup = AnnotationGroup(region, get_heart_term("right atrium myocardium"))
        aSeptumGroup = AnnotationGroup(region, get_heart_term("interatrial septum"))
        fossaGroup = AnnotationGroup(region, get_heart_term("fossa ovalis"))
        lipvGroup = AnnotationGroup(region, get_heart_term("left inferior pulmonary vein"))
        lspvGroup = AnnotationGroup(region, get_heart_term("left superior pulmonary vein"))
        ripvGroup = AnnotationGroup(region, get_heart_term("right inferior pulmonary vein"))
        rspvGroup = AnnotationGroup(region, get_heart_term("right superior pulmonary vein"))
        ivcInletGroup = AnnotationGroup(region, get_heart_term("inferior vena cava inlet"))
        svcInletGroup = AnnotationGroup(region, get_heart_term("superior vena cava inlet"))
        annotationGroups = [ laGroup, raGroup, aSeptumGroup, fossaGroup, lipvGroup, lspvGroup, ripvGroup, rspvGroup, ivcInletGroup, svcInletGroup ]

        ##############
        # Create nodes
        ##############

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # LA/RA inlet elements are linear through the wall, hence their nodes do not have D_DS3 parameters
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)
        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)

        aBaseOuterMajorMag = aBaseInnerMajorMag + aBaseSlopeLength
        aBaseOuterMinorMag = aBaseInnerMinorMag + aBaseSlopeLength

        # GRC fudge factor:
        aOuterSeptumHeight = 0.8*aOuterHeight
        # GRC fudge factor:
        aSeptumBaseRowHeight = 0.2*aOuterHeight

        # convert base major/minor sizes to equator using arc up

        aEquatorOuterMajorMag = aBaseOuterMajorMag/math.sin(aOuterMajorArcUpRadians)
        aEquatorOuterMinorMag = aBaseOuterMinorMag/math.sin(aOuterMinorArcUpRadians)
        aEquatorInnerMajorMag = aEquatorOuterMajorMag - aFreeWallThickness
        aEquatorInnerMinorMag = aEquatorOuterMinorMag - aFreeWallThickness

        # following are angles in radians around LA ellipse from major axis
        axInner = aBaseInnerMajorMag*math.cos(aMajorAxisRadians)
        bxInner = aBaseInnerMinorMag*math.sin(aMajorAxisRadians)
        laSeptumRadians = math.atan2(bxInner, axInner)
        laCentreX = -0.5*aSeptumThickness - axInner*math.cos(laSeptumRadians) - bxInner*math.sin(laSeptumRadians)
        #laCruxLeftRadians = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, laSeptumRadians, \
        #    (aSeptumBaseLength/elementsCountAroundAtrialSeptum)*(0.5*elementsCountAroundAtrialSeptum + 1.0))
        axOuter = aBaseOuterMajorMag*math.cos(aMajorAxisRadians)
        bxOuter = aBaseOuterMinorMag*math.sin(aMajorAxisRadians)
        aCruxSideOffset = lvOutletOuterRadius*math.sin(math.pi/3.0)
        laCruxLeftRadians = getEllipseRadiansToX(axOuter, bxOuter, -aCruxSideOffset - laCentreX, math.pi*0.5)
        #print('axInner',axInner,'bxInner',bxInner,'laCentreX',laCentreX)
        #print('laSeptumRadians',laSeptumRadians,'laCruxLeftRadians',laCruxLeftRadians)
        laCentreY = 0.0
        laCentreZ = 0.0
        raSeptumRadians = math.pi*2.0 - laSeptumRadians
        raCruxRightRadians = math.pi*2.0 - laCruxLeftRadians
        raCentreX = -laCentreX
        raCentreY = laCentreY
        raCentreZ = 0.0

        aOuterScaleZ = aOuterHeight/(1.0 - math.cos(aOuterMinorArcUpRadians))
        aInnerScaleZ = (aOuterHeight - aFreeWallThickness)/(1.0 - math.cos(aOuterMinorArcUpRadians))
        aOuterMajorScaleZ = aOuterHeight/(1.0 - math.cos(aOuterMajorArcUpRadians))
        aInnerMajorScaleZ = (aOuterHeight - aFreeWallThickness)/(1.0 - math.cos(aOuterMajorArcUpRadians))
        zOffset = aOuterHeight - aOuterScaleZ  # so base outer is at z = 0

        cruxLeftX = -aCruxSideOffset
        #cruxLeftX = laCentreX + math.cos(laCruxLeftRadians)*aBaseOuterMajorMag*math.cos(-aMajorAxisRadians) \
        #                      + math.sin(laCruxLeftRadians)*aBaseOuterMinorMag*-math.sin(-aMajorAxisRadians)
        cruxLeftY = laCentreY + math.cos(laCruxLeftRadians)*aBaseOuterMajorMag*math.sin(-aMajorAxisRadians) \
                              + math.sin(laCruxLeftRadians)*aBaseOuterMinorMag*math.cos(-aMajorAxisRadians)
        #print('crux left', cruxLeftX, cruxLeftY)

        cruxCentreX = 0.0
        cruxCentreY = cruxLeftY - lvOutletOuterRadius*(1.0 - math.cos(math.pi/3.0))
        #node = nodes.createNode(nodeIdentifier, nodetemplate)
        #cache.setNode(node)
        #result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ cruxCentreX, cruxCentreY, 0.0 ])
        #nodeIdentifier += 1

        aMinorToMajorRadians = aOuterMajorArcUpRadians/aOuterMinorArcUpRadians

        laNodeId = [ [], [] ]
        raNodeId = [ [], [] ]
        laApexNodeId = []
        raApexNodeId = []

        # crux to be collapsed
        n1CruxL = elementsCountAroundAtrialSeptum//2
        n1CruxR = elementsCountAroundAtria - ((elementsCountAroundAtrialSeptum + 1)//2)
        n1SeptumRange = []
        for n1 in range(elementsCountAroundAtria):
            if (n1 < n1CruxL) or (n1 > n1CruxR):
                n1SeptumRange.append(n1)
        #print('n1SeptumRange', n1SeptumRange)

        # compute radians around based on base outer major and minor axis sizes
        atrialPerimeterLength = getApproximateEllipsePerimeter(aBaseOuterMajorMag, aBaseOuterMinorMag)
        atrialSeptumElementLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumRadians, laCruxLeftRadians) \
            /(0.5*elementsCountAroundAtrialSeptum + 1.0)
        atrialFreeWallElementLength = (atrialPerimeterLength - atrialSeptumElementLength*(elementsCountAroundAtrialSeptum + 2)) \
            / (elementsCountAroundAtria - elementsCountAroundAtrialSeptum - 2)
        atrialTransitionElementLength = 0.5*(atrialSeptumElementLength + atrialFreeWallElementLength)
        #atrialPerimeterLengthTmp = atrialSeptumElementLength*(elementsCountAroundAtrialSeptum + 1) + 2.0*atrialTransitionElementLength \
        #    + (elementsCountAroundAtria - elementsCountAroundAtrialSeptum - 3)*atrialFreeWallElementLength
        #print('lengths:',(elementsCountAroundAtrialSeptum + 1),'*',atrialSeptumElementLength, '+ 2 *', \
        #    atrialTransitionElementLength,'+',(elementsCountAroundAtria - elementsCountAroundAtrialSeptum - 3),'*',atrialFreeWallElementLength,
        #    '=', atrialPerimeterLengthTmp, ' VS ', atrialPerimeterLength)
        laRadians = []
        radiansAround = laSeptumRadians
        oddSeptumElementsCount = (elementsCountAroundAtrialSeptum % 2) == 1
        if oddSeptumElementsCount:
            radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, 0.5*atrialSeptumElementLength)
            n1la2ra = [ (elementsCountAroundAtria - n1 - 1) for n1 in range(elementsCountAroundAtria) ]
        else:
            n1la2ra = [ (-n1 % elementsCountAroundAtria) for n1 in range(elementsCountAroundAtria) ]
        #print('n1la2ra', n1la2ra)
        lan1CruxLimit = elementsCountAroundAtrialSeptum//2 + 1
        lan1SeptumLimit = elementsCountAroundAtria - (elementsCountAroundAtrialSeptum + 1)//2 - 1
        #print('lan1CruxLimit', lan1CruxLimit, 'lan1SeptumLimit', lan1SeptumLimit)
        for n1 in range(elementsCountAroundAtria):
            laRadians.append(radiansAround)
            if (n1 < lan1CruxLimit) or (n1 > lan1SeptumLimit):
                elementLength = atrialSeptumElementLength
            elif (n1 == lan1CruxLimit) or (n1 == lan1SeptumLimit):
                elementLength = atrialTransitionElementLength
            else:
                elementLength = atrialFreeWallElementLength
            radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)
        raRadians = [ 2.0*math.pi - laRadians[n1la2ra[n1]] for n1 in range(elementsCountAroundAtria) ]
        ran1SeptumLimit = elementsCountAroundAtrialSeptum//2
        ran1CruxLimit = elementsCountAroundAtria - ran1SeptumLimit - 1

        # detect fibrous ring, means merge with ventriclesbase2
        laBaseNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        raBaseNodeId = [ [-1]*elementsCountAroundAtria, [-1]*elementsCountAroundAtria ]
        # nodes 3 and 7 of the 3rd fibrous ring element are on base centre of atrial septum
        leftFibrousRingGroup = fm.findFieldByName('left fibrous ring').castGroup()
        mergeWithBase = False
        if leftFibrousRingGroup.isValid():
            leftFibrousRingMeshGroup = leftFibrousRingGroup.getFieldElementGroup(fm.findMeshByDimension(3)).getMeshGroup()
            elementiter = leftFibrousRingMeshGroup.createElementiterator()
            element = elementiter.next()
            element = elementiter.next()
            element = elementiter.next()
            eft1 = element.getElementfieldtemplate(coordinates, -1)
            aSeptumNode1 = element.getNode(eft1, 3)
            aSeptumNode2 = element.getNode(eft1, 7)
            if aSeptumNode1.isValid() and aSeptumNode2.isValid():
                # transform ventricles and base to align with atria (easier than other way around)
                cache.setNode(aSeptumNode1)
                result, sx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                cache.setNode(aSeptumNode2)
                result, sx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                aSeptumCentreY = laCentreY + math.cos(laSeptumRadians)*aBaseInnerMajorMag*math.sin(-aMajorAxisRadians) \
                                    + math.sin(laSeptumRadians)*aBaseInnerMinorMag*math.cos(-aMajorAxisRadians)
                tx1 = [ -0.5*aSeptumThickness, aSeptumCentreY, -aBaseSlopeHeight ]
                ax1 = vector.normalise([ (sx2[c] - sx1[c]) for c in range(3) ])
                ax2 = [ -ax1[1], ax1[0], ax1[2] ]
                ax3 = [ 0.0, 0.0, 1.0 ]

                transmat = fm.createFieldConstant(ax1 + ax2 + ax3)
                newCoordinates = fm.createFieldAdd(fm.createFieldMatrixMultiply(3, transmat, fm.createFieldSubtract(coordinates, fm.createFieldConstant(sx1))), fm.createFieldConstant(tx1))
                fieldassignment = coordinates.createFieldassignment(newCoordinates)
                result = fieldassignment.assign()
                # get node ids on base to use for first row atria
                # happen to know last nodes are on fibrous ring, and in which order
                baseNodeIdentifier = nodeIdentifier
                for n3 in range(1, -1, -1):
                    for i in range(1, -1, -1):
                        aBaseNodeId = laBaseNodeId[n3] if (i == 0) else raBaseNodeId[n3]
                        for n1 in range(elementsCountAroundAtria - 1, -1, -1):
                            if (n3 == 1) and \
                                (((i == 0) and ((n1 < (lan1CruxLimit - 1)) or (n1 > (lan1SeptumLimit + 2)))) or \
                                 ((i == 1) and ((n1 < ran1SeptumLimit) or (n1 > ran1CruxLimit)))):
                                # get node from other side below
                                continue
                            baseNodeIdentifier -= 1
                            aBaseNodeId[n1] = baseNodeIdentifier
                # tie up overlaps - only implemented for elementsCountAroundAtrialSeptum == 2:
                laBaseNodeId[1][ 0] = raBaseNodeId[0][0]
                raBaseNodeId[1][-1] = laBaseNodeId[1][1]
                raBaseNodeId[1][ 0] = laBaseNodeId[0][0]
                #print('laBaseNodeId', laBaseNodeId)
                #print('raBaseNodeId', raBaseNodeId)
                mergeWithBase = True

        for n3 in range(2):

            aScaleZ = aInnerScaleZ if (n3 == 0) else aOuterScaleZ
            aMajorScaleZ = aInnerMajorScaleZ if (n3 == 0) else aOuterMajorScaleZ
            aEquatorMajorMag = aEquatorInnerMajorMag if (n3 == 0) else aEquatorOuterMajorMag
            aEquatorMinorMag = aEquatorInnerMinorMag if (n3 == 0) else aEquatorOuterMinorMag
            distanceUpMinor = getEllipseArcLength(aScaleZ, aEquatorMinorMag, 0, aOuterMinorArcUpRadians)
            elementSizeUpMinor = distanceUpMinor/elementsCountUpAtria
            radiansUp = -aOuterMinorArcUpRadians

            radiansUpMinor = []
            radiansUpMajor = []
            derivativesUpMajor = []
            radiansUpPrev = updateEllipseAngleByArcLength(aScaleZ, aEquatorMinorMag, radiansUp, -elementSizeUpMinor)
            zPrev = aScaleZ*math.cos(radiansUpPrev)
            aMajorOffsetZ = aMajorScaleZ - aScaleZ
            radiansUpMajorPrev = getEllipseRadiansToX(aMajorScaleZ, 0.0, aMajorOffsetZ + zPrev, radiansUpPrev*aMinorToMajorRadians)
            for n2 in range(elementsCountUpAtria):
                radiansUpMinor.append(radiansUp)
                z = aScaleZ*math.cos(radiansUp)
                if n2 == 0:
                    radiansUpMajorCurr = getEllipseRadiansToX(aMajorScaleZ, 0.0, aMajorOffsetZ + z, radiansUp*aMinorToMajorRadians)
                radiansUpMajor.append(radiansUpMajorCurr)
                radiansUp = updateEllipseAngleByArcLength(aScaleZ, aEquatorMinorMag, radiansUp, elementSizeUpMinor)
                zNext = aScaleZ*math.cos(radiansUp)
                radiansUpMajorNext = getEllipseRadiansToX(aMajorScaleZ, 0.0, aMajorOffsetZ + zNext, radiansUp*aMinorToMajorRadians)
                derivativeUpMajor = 0.5*getEllipseArcLength(aMajorScaleZ, aEquatorMajorMag, radiansUpMajorPrev, radiansUpMajorNext)
                derivativesUpMajor.append(derivativeUpMajor)
                radiansUpMajorPrev = radiansUpMajorCurr
                radiansUpMajorCurr = radiansUpMajorNext
            apexDerivativeUpMajor = getEllipseArcLength(aMajorScaleZ, aEquatorMajorMag, radiansUpMajorPrev, radiansUpMajorCurr)

            #if n3 == 1:
            #    print('\nradiansUpMinor', radiansUpMinor)
            #    print('radiansUpMajor', radiansUpMajor)
            #    print('derivativesUpMajor', derivativesUpMajor)

            # regular nodes up atria
            for n2 in range(elementsCountUpAtria):

                radiansUp = radiansUpMinor[n2]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)

                # radiansUpMajor = radiansUpMajor[ # radiansUp*aMinorToMajorRadians
                cosRadiansUpMajor = math.cos(radiansUpMajor[n2])
                sinRadiansUpMajor = math.sin(radiansUpMajor[n2])

                if (n3 == 0) and (n2 == 0):
                    aMajorMag = aBaseInnerMajorMag
                    aMinorMag = aBaseInnerMinorMag
                else:
                    aMajorMag = -aEquatorMajorMag*sinRadiansUpMajor
                    aMinorMag = -aEquatorMinorMag*sinRadiansUp

                #print('n', n3, n2, 'radiansUp', radiansUp, 'radiansUpMajor', radiansUpMajor[n2])
                #print('aMajorMag',aMajorMag,'aMinorMag',aMinorMag)

                laDerivatives = []
                finalArcLength = prevArcLength = getEllipseArcLength(aMajorMag, aMinorMag, laRadians[-1] - 2.0*math.pi, laRadians[0])
                for n1 in range(elementsCountAroundAtria):
                    if n1 == (elementsCountAroundAtria - 1):
                        nextArcLength = finalArcLength
                    else:
                        nextArcLength = getEllipseArcLength(aMajorMag, aMinorMag, laRadians[n1], laRadians[n1 + 1])
                    if (n1 <= lan1CruxLimit) or (n1 > lan1SeptumLimit):
                        arcLength = min(prevArcLength, nextArcLength)
                    else:
                        arcLength = max(prevArcLength, nextArcLength)
                    laDerivatives.append(arcLength)
                    prevArcLength = nextArcLength
                #print('Radians:',[math.degrees(r) for r in laRadians])
                raDerivatives = [ laDerivatives[n1la2ra[n1]] for n1 in range(elementsCountAroundAtria)]

                if (n3 == 0) and (n2 == 0):
                    z = zOffset + aOuterScaleZ*math.cos(aOuterMinorArcUpRadians) - aBaseSlopeHeight
                else:
                    z = zOffset + aScaleZ*cosRadiansUp

                laLayerNodeId = [-1]*elementsCountAroundAtria
                laNodeId[n3].append(laLayerNodeId)
                raLayerNodeId = [-1]*elementsCountAroundAtria
                raNodeId[n3].append(raLayerNodeId)

                for i in range(2):
                    if i == 0:
                        # left
                        centreX = laCentreX
                        centreY = laCentreY
                        majorAxisRadians = -aMajorAxisRadians
                        aRadians = laRadians
                        aDerivatives = laDerivatives
                        aLayerNodeId = laLayerNodeId
                        n1Start = lan1CruxLimit - 1
                        n1Stop = lan1SeptumLimit + 1
                        aBaseNodeId = laBaseNodeId[n3]
                    else:
                        # right
                        centreX = raCentreX
                        centreY = raCentreY
                        majorAxisRadians = math.pi + aMajorAxisRadians
                        aRadians = raRadians
                        aDerivatives = raDerivatives
                        aLayerNodeId = raLayerNodeId
                        n1Start = n1la2ra[lan1SeptumLimit + 1]
                        n1Stop = n1la2ra[lan1CruxLimit - 1]
                        aBaseNodeId = raBaseNodeId[n3]

                    sinMajorAxisRadians = math.sin(majorAxisRadians)
                    cosMajorAxisRadians = math.cos(majorAxisRadians)
                    aMajorX =  aMajorMag*cosMajorAxisRadians
                    aMajorY =  aMajorMag*sinMajorAxisRadians
                    aMinorX = -aMinorMag*sinMajorAxisRadians
                    aMinorY =  aMinorMag*cosMajorAxisRadians
                    #print('aMajor', aMajorX, aMajorY,'aMinor',aMinorX,aMinorY)
                    #print('n1 range: ', n1Start, n1Stop)

                    for n1 in range(elementsCountAroundAtria):
                        mergeNode = False
                        if mergeWithBase and (n2 == 0):
                            aLayerNodeId[n1] = aBaseNodeId[n1]
                            node = nodes.findNodeByIdentifier(aLayerNodeId[n1])
                            mergeNode = True
                        if (n2 > 1) and ((n1 <= n1Start) or (n1 >= n1Stop)):
                            continue
                        radiansAround = aRadians[n1]
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        if (n3 == 1) and ((n1 in n1SeptumRange) or ((n2 == 0) and (i == 1) and (n1 == n1CruxR))):
                            continue  # outer septum node created by other side
                        if not mergeNode:
                            node = nodes.createNode(nodeIdentifier, nodetemplate)
                            aLayerNodeId[n1] = nodeIdentifier
                        cache.setNode(node)
                        if (n3 == 1) and (i == 0) and (n2 == 0) and (n1 == n1CruxL):
                            x = [ cruxCentreX, cruxCentreY, z ]
                            dx_ds1 = [ math.pi*lvOutletOuterRadius/-3.0, 0.0, 0.0 ]
                            dx_ds2 = [ 0.0, 0.0, aSeptumBaseRowHeight + 2.0*aBaseSlopeHeight ]
                        else:
                            x = [
                                centreX + cosRadiansAround*aMajorX + sinRadiansAround*aMinorX,
                                centreY + cosRadiansAround*aMajorY + sinRadiansAround*aMinorY,
                                z ]
                            d1x = -sinRadiansAround*aMajorX + cosRadiansAround*aMinorX
                            d1y = -sinRadiansAround*aMajorY + cosRadiansAround*aMinorY
                            scale1 = aDerivatives[n1]/math.sqrt(d1x*d1x + d1y*d1y)
                            dx_ds1 = [ d1x*scale1, d1y*scale1, 0.0 ]

                            d2Minor = [
                                 cosRadiansUp*aEquatorMinorMag*sinMajorAxisRadians,
                                -cosRadiansUp*aEquatorMinorMag*cosMajorAxisRadians,
                                -sinRadiansUp*aScaleZ ]
                            minorScale = elementSizeUpMinor/vector.magnitude(d2Minor)
                            d2Minor = [ d*minorScale for d in d2Minor ]
                            d2Major = [
                                -cosRadiansUpMajor*aEquatorMajorMag*cosMajorAxisRadians,
                                -cosRadiansUpMajor*aEquatorMajorMag*sinMajorAxisRadians,
                                -sinRadiansUpMajor*aMajorScaleZ ]
                            majorScale = derivativesUpMajor[n2]/vector.magnitude(d2Major)
                            d2Major = [ d*majorScale for d in d2Major ]

                            dx_ds2 = [
                                cosRadiansAround*d2Major[0] + sinRadiansAround*d2Minor[0],
                                cosRadiansAround*d2Major[1] + sinRadiansAround*d2Minor[1],
                                cosRadiansAround*cosRadiansAround*d2Major[2] + sinRadiansAround*sinRadiansAround*d2Minor[2] ]

                            if n1 in n1SeptumRange:
                                if n2 == 1:
                                    x[2] = aSeptumBaseRowHeight
                                    dx_ds1[2] = 0.0
                                    if i == 0:
                                        nx = min(x[0], -0.5*aSeptumThickness)
                                    else:
                                        nx = max(x[0],  0.5*aSeptumThickness)
                                    if nx != x[0]:
                                        x[0] = nx
                                        dx_ds1[0] = 0.0
                                dx_ds2 = [ 0.0, 0.0, (aSeptumBaseRowHeight + 2.0*aBaseSlopeHeight) if (n2 == 0) else aSeptumBaseRowHeight ]

                        if not mergeNode:
                            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                        # derivative 3 is set later
                        nodeIdentifier += 1

            # apex nodes
            n2 = elementsCountUpAtria - 1
            for i in range(2):
                if i == 0:
                    # left
                    nid1 = laNodeId[n3][n2][lan1SeptumLimit]
                    nid2 = laNodeId[n3][n2][lan1CruxLimit]
                    apexNodeId = laApexNodeId
                else:
                    # right
                    nid1 = raNodeId[n3][n2][n1la2ra[lan1SeptumLimit]]
                    nid2 = raNodeId[n3][n2][n1la2ra[lan1CruxLimit]]
                    apexNodeId = raApexNodeId
                nid3 = (nid1 + nid2)//2

                node1 = nodes.findNodeByIdentifier(nid1)
                node2 = nodes.findNodeByIdentifier(nid2)
                x, dx_ds2, dx_ds1, dx_ds3 = interpolateNodesCubicHermite(cache, coordinates, 0.5, aFreeWallThickness, \
                    node1, Node.VALUE_LABEL_D_DS2,  2.0, Node.VALUE_LABEL_D_DS1, 1.0, \
                    node2, Node.VALUE_LABEL_D_DS2, -2.0, Node.VALUE_LABEL_D_DS1, -1.0)
                node3 = nodes.findNodeByIdentifier(nid3)
                cache.setNode(node3)
                result, x3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                result, d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                d3 = [ -d for d in d3 ]
                arcLength = interp.computeCubicHermiteArcLength(x, dx_ds1, x3, d3, True)
                scale1 = (2.0*arcLength - vector.magnitude(d3))/vector.magnitude(dx_ds1)
                dx_ds1 = [ d*scale1 for d in dx_ds1 ]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                apexNodeId.append(nodeIdentifier)
                cache.setNode(node)
                #print(n3, i, 'project apex', nid1, nid2, nid3, vector.magnitude(dx_ds1))
                x[2] = aOuterHeight - aFreeWallThickness if (n3 == 0) else aOuterHeight
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # fix inner base derivative 2 to fit incline
        for i in range(2):
            aNodeId = laNodeId[0] if (i == 0) else raNodeId[0]
            for n1 in range(elementsCountAroundAtria):
                if n1 in n1SeptumRange:
                    continue
                node1 = nodes.findNodeByIdentifier(aNodeId[0][n1])
                node2 = nodes.findNodeByIdentifier(aNodeId[1][n1])
                dx_ds2 = computeNodeDerivativeHermiteLagrange(cache, coordinates, node2, Node.VALUE_LABEL_D_DS2, -1.0, node1, -1.0)
                cache.setNode(node1)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)

        # transfer inner septum nodes to outer on opposite side, set derivative 3 to be node difference
        raNodeId[1][0][n1CruxR] = laNodeId[1][0][n1CruxL]
        for n2 in range(elementsCountUpAtria):
            for n1 in n1SeptumRange:
                na1 = n1
                nb1 = n1la2ra[na1]
                laNodeId[1][n2][na1] = raNodeId[0][n2][nb1]
                raNodeId[1][n2][nb1] = laNodeId[0][n2][na1]
            for i in range(2):
                aNodeId = laNodeId if (i == 0) else raNodeId
                for n1 in range(elementsCountAroundAtria):
                    node_i = nodes.findNodeByIdentifier(aNodeId[0][n2][n1])
                    node_o = nodes.findNodeByIdentifier(aNodeId[1][n2][n1])
                    cache.setNode(node_o)
                    result, x_o = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                    cache.setNode(node_i)
                    result, x_i = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                    dx_ds3 = [ (x_o[i] - x_i[i]) for i in range(3) ]
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if n1 in n1SeptumRange:
                        dx_ds3 = [ -d for d in dx_ds3 ]
                    cache.setNode(node_o)
                    result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # fix crux centre dx_ds3:
        cache.setNode(nodes.findNodeByIdentifier(laNodeId[0][0][n1CruxL]))
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        cache.setNode(nodes.findNodeByIdentifier(raNodeId[0][0][n1CruxR]))
        result, x2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        cache.setNode(nodes.findNodeByIdentifier(laNodeId[1][0][n1CruxL]))
        result, xc = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
        d1 = [ (xc[c] - x1[c]) for c in range(3) ]
        d2 = [ (xc[c] - x2[c]) for c in range(3) ]
        dx_ds3 = [ d1[0] + d2[0], d1[1] + d2[1], d1[2] + d2[2] ]
        scale = vector.magnitude(d1)/vector.magnitude(dx_ds3)
        dx_ds3 = [ d*scale for d in dx_ds3 ]
        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3 )

        # left atrial septum posterior/anterior nodeId[n3][n2], also mid (top) septum
        laspNodeId = [ [ None ]*(elementsCountUpAtria + 1), [ None ]*(elementsCountUpAtria + 1) ]
        lasaNodeId = [ [ None ]*(elementsCountUpAtria + 1), [ None ]*(elementsCountUpAtria + 1) ]
        # right atrial septum posterior/anterior nodeId[n3][n2], also mid (top) septum
        raspNodeId = [ [ None ]*(elementsCountUpAtria + 1), [ None ]*(elementsCountUpAtria + 1) ]
        rasaNodeId = [ [ None ]*(elementsCountUpAtria + 1), [ None ]*(elementsCountUpAtria + 1) ]
        na = elementsCountAroundAtrialSeptum//2
        np = elementsCountAroundAtria - ((elementsCountAroundAtrialSeptum + 1)//2)

        # get atrial septum peak
        node1 = nodes.findNodeByIdentifier(laApexNodeId[1])
        cache.setNode(node1)
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        aSeptumBaseCentreY = laCentreY + aBaseInnerMajorMag*math.sin(-aMajorAxisRadians)*math.cos(laSeptumRadians) + aBaseInnerMinorMag*math.cos(-aMajorAxisRadians)*math.sin(laSeptumRadians)
        x2 = [ 0.0, aSeptumBaseCentreY, aOuterSeptumHeight - aFreeWallThickness ]
        d2 = [ 1.0, 0.0, 0.0 ]
        #print('septum top centre ', x2)
        arcLength = interp.computeCubicHermiteArcLength(x1, d1, x2, d2, True)
        scale1 = arcLength/vector.magnitude(d1)
        d1 = [ d*scale1 for d in d1 ]
        scale2 = arcLength/vector.magnitude(d2)
        d2 = [ d*arcLength for d in d2 ]
        # GRC fudge factor:
        xi = 0.7
        lasmx = interp.interpolateCubicHermite(x1, d1, x2, d2, xi)
        lasmd1 = interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        scale1 = (1.0 - xi)*2.0
        lasmd1 = [ d*scale1 for d in lasmd1 ]
        lasmd2 = vector.normalise(vector.crossproduct3([0.0, 0.0, 1.0], lasmd1))
        lasmd3 = vector.crossproduct3(lasmd1, lasmd2)
        scale3 = aFreeWallThickness/vector.magnitude(lasmd3)
        lasmd3 = [ d*scale3 for d in lasmd3 ]
        lasmCurvature1 = -interp.getCubicHermiteCurvature(x1, d1, x2, d2, vector.normalise(lasmd3), xi)

        for n3 in range(2):
            for n2 in range(2):
                laspNodeId[n3][n2] = laNodeId[n3][n2][np]
                lasaNodeId[n3][n2] = laNodeId[n3][n2][na]
                raspNodeId[n3][n2] = raNodeId[n3][n2][n1la2ra[np]]
                rasaNodeId[n3][n2] = raNodeId[n3][n2][n1la2ra[na]]
        n2 = 1
        cache.setNode(nodes.findNodeByIdentifier(laspNodeId[0][n2]))
        result, laspx =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, laspd1 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, laspd2 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        result, laspd3 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
        cache.setNode(nodes.findNodeByIdentifier(laspNodeId[1][n2]))
        result, laspd1_outer =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        laspCurvature1 = (vector.magnitude(laspd1_outer)/vector.magnitude(laspd1) - 1.0)/aFreeWallThickness
        cache.setNode(nodes.findNodeByIdentifier(lasaNodeId[0][n2]))
        result, lasax =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, lasad1 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, lasad2 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        result, lasad3 =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
        cache.setNode(nodes.findNodeByIdentifier(lasaNodeId[1][n2]))
        result, lasad1_outer =  coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        lasaCurvature1 = (vector.magnitude(lasad1_outer)/vector.magnitude(lasad1) - 1.0)/aFreeWallThickness

        # create points arcing over septum peak
        for n3 in range(2):
            #print(str(n3) + '. la septum nodes p', laspNodeId[n3][1],'a', lasaNodeId[n3][1])
            #print(str(n3) + '. ra septum nodes p', raspNodeId[n3][1],'a', rasaNodeId[n3][1])
            arcLength = interp.computeCubicHermiteArcLength(laspx, laspd2, lasmx, lasmd2, True)
            x1 = laspx
            scale1 = arcLength/vector.magnitude(laspd2)
            d1 = [ d*scale1 for d in laspd2 ]
            x2 = lasmx
            scale2 = arcLength/vector.magnitude(lasmd2)
            d2 = [ d*scale2 for d in lasmd2 ]
            derivativeScale = arcLength/(elementsCountUpAtria - 1.0)
            for n2 in range(2, elementsCountUpAtria + 1):
                xi = (n2 - 1.0)/(elementsCountUpAtria - 1.0)
                xr = 1.0 - xi
                x = interp.interpolateCubicHermite(x1, d1, x2, d2, xi)
                dx_ds1 = [ (xi*lasmd1[c] + xr*laspd1[c]) for c in range(3) ]
                dx_ds2 = interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
                scale2 = derivativeScale/vector.magnitude(dx_ds2)
                dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                dx_ds3 = vector.crossproduct3(dx_ds1, dx_ds2)
                scale3 = aFreeWallThickness/vector.magnitude(dx_ds3)
                dx_ds3 = [ d*scale3 for d in dx_ds3 ]
                if n3 == 1:
                    curvature1 = xi*lasmCurvature1 + xr*laspCurvature1
                    curvatureScale1 = 1.0 + aFreeWallThickness*curvature1
                    dx_ds1 = [ d*curvatureScale1 for d in dx_ds1 ]
                    radialVector = vector.normalise(dx_ds3)
                    curvature2 = -interp.getCubicHermiteCurvature(x1, d1, x2, d2, radialVector, xi)
                    curvatureScale2 = 1.0 + aFreeWallThickness*curvature2
                    dx_ds2 = [ d*curvatureScale2 for d in dx_ds2 ]
                    x = [ (x[c] + dx_ds3[c]) for c in range(3) ]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                laspNodeId[n3][n2] = nodeIdentifier
                nodeIdentifier += 1

                # mirror RA about x = 0, then reverse dx_ds1
                x[0] = -x[0]
                dx_ds1 = [ dx_ds1[0], -dx_ds1[1], -dx_ds1[2] ]
                dx_ds2[0] = -dx_ds2[0]
                dx_ds3[0] = -dx_ds3[0]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                raspNodeId[n3][n2] = nodeIdentifier
                nodeIdentifier += 1

            arcLength = interp.computeCubicHermiteArcLength(lasax, lasad2, lasmx, lasmd2, True)
            x1 = lasax
            scale1 = arcLength/vector.magnitude(lasad2)
            d1 = [ d*scale1 for d in lasad2 ]
            x2 = lasmx
            scale2 = -arcLength/vector.magnitude(lasmd2)
            d2 = [ d*scale2 for d in lasmd2 ]
            derivativeScale = arcLength/(elementsCountUpAtria - 1.0)
            for n2 in range(2, elementsCountUpAtria):
                xi = (n2 - 1.0)/(elementsCountUpAtria - 1.0)
                xr = 1.0 - xi
                x = interp.interpolateCubicHermite(x1, d1, x2, d2, xi)
                dx_ds1 = [ (xi*-lasmd1[c] + xr*lasad1[c]) for c in range(3) ]
                dx_ds2 = interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
                scale2 = derivativeScale/vector.magnitude(dx_ds2)
                dx_ds2 = [ d*scale2 for d in dx_ds2 ]
                dx_ds3 = vector.crossproduct3(dx_ds1, dx_ds2)
                scale3 = aFreeWallThickness/vector.magnitude(dx_ds3)
                dx_ds3 = [ d*scale3 for d in dx_ds3 ]
                if n3 == 1:
                    curvature1 = xi*lasmCurvature1 + xr*lasaCurvature1
                    curvatureScale1 = 1.0 + aFreeWallThickness*curvature1
                    dx_ds1 = [ d*curvatureScale1 for d in dx_ds1 ]
                    radialVector = vector.normalise(dx_ds3)
                    curvature2 = -interp.getCubicHermiteCurvature(x1, d1, x2, d2, radialVector, xi)
                    curvatureScale2 = 1.0 + aFreeWallThickness*curvature2
                    dx_ds2 = [ d*curvatureScale2 for d in dx_ds2 ]
                    x = [ (x[c] + dx_ds3[c]) for c in range(3) ]

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                lasaNodeId[n3][n2] = nodeIdentifier
                nodeIdentifier += 1

                # mirror RA about x = 0, then reverse dx_ds1
                x[0] = -x[0]
                dx_ds1 = [ dx_ds1[0], -dx_ds1[1], -dx_ds1[2] ]
                dx_ds2[0] = -dx_ds2[0]
                dx_ds3[0] = -dx_ds3[0]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                rasaNodeId[n3][n2] = nodeIdentifier
                nodeIdentifier += 1

            lasaNodeId[n3][-1] = laspNodeId[n3][-1]
            rasaNodeId[n3][-1] = raspNodeId[n3][-1]

        # atrial septum centre nodes and nodes around fossa ovalis

        aSeptumCentreY = lasmx[1]
        aSeptumCentreZ = 0.5*(lasmx[2] + aSeptumBaseRowHeight)
        fossaMagY = 0.5*(aSeptumCentreY - laspx[1])
        fossaMagZ = 0.25*(lasmx[2] - aSeptumBaseRowHeight)
        elementsCountAroundFossa = elementsCountAroundAtrialSeptum + 2*(elementsCountUpAtria - 1)
        fossaPerimeterLength = getApproximateEllipsePerimeter(fossaMagY, fossaMagZ)
        elementSizeAroundFossa = fossaPerimeterLength/elementsCountAroundFossa
        radiansAround = -0.5*math.pi
        if (elementsCountAroundFossa%2) == 1:
            radiansAround = updateEllipseAngleByArcLength(fossaMagY, fossaMagZ, radiansAround, 0.5*elementSizeAroundFossa)
        fossaRadiansAround = []
        for n1 in range(elementsCountAroundFossa):
            fossaRadiansAround.append(radiansAround)
            radiansAround = updateEllipseAngleByArcLength(fossaMagY, fossaMagZ, radiansAround, elementSizeAroundFossa)

        fossaCentreNodeId = [ -1, -1 ]
        fossaNodeId = [ [ -1 ]*elementsCountAroundFossa, [ -1 ]*elementsCountAroundFossa ]
        dx_ds3 = [ aSeptumThickness, 0.0, 0.0 ]
        for n3 in range(2):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            x = [ aSeptumThickness*(-0.5 if (n3 == 0) else 0.5), aSeptumCentreY, aSeptumCentreZ ]
            dx_ds1 = [ 0.0, fossaMagY, 0.0 ]
            dx_ds2 = [ 0.0, 0.0, fossaMagZ ]
            result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            fossaCentreNodeId[n3] = nodeIdentifier
            nodeIdentifier += 1

            for n1 in range(elementsCountAroundFossa):
                radiansAround = fossaRadiansAround[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[1] = aSeptumCentreY + fossaMagY*cosRadiansAround
                x[2] = aSeptumCentreZ + fossaMagZ*sinRadiansAround
                dx_ds1[1] = -fossaMagY*sinRadiansAround
                dx_ds1[2] =  fossaMagZ*cosRadiansAround
                scale1 = elementSizeAroundFossa/vector.magnitude(dx_ds1)
                dx_ds1[1] *= scale1
                dx_ds1[2] *= scale1
                dx_ds2[1] = -fossaMagY*cosRadiansAround
                dx_ds2[2] = -fossaMagZ*sinRadiansAround
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                fossaNodeId[n3][n1] = nodeIdentifier
                nodeIdentifier += 1

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)

        laMeshGroup = laGroup.getMeshGroup(mesh)
        raMeshGroup = raGroup.getMeshGroup(mesh)
        aSeptumMeshGroup = aSeptumGroup.getMeshGroup(mesh)
        fossaMeshGroup = fossaGroup.getMeshGroup(mesh)
        lipvMeshGroup = lipvGroup.getMeshGroup(mesh)
        lspvMeshGroup = lspvGroup.getMeshGroup(mesh)
        ripvMeshGroup = ripvGroup.getMeshGroup(mesh)
        rspvMeshGroup = rspvGroup.getMeshGroup(mesh)
        ivcInletMeshGroup = ivcInletGroup.getMeshGroup(mesh)
        svcInletMeshGroup = svcInletGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        e1FreeWallStart = elementsCountAroundAtrialSeptum + 2
        e1lspv = e1FreeWallStart + (elementsCountAroundAtria - elementsCountAroundAtrialSeptum - 1)//2
        e1lipv = e1FreeWallStart + (elementsCountAroundAtria - elementsCountAroundAtrialSeptum + 1)//2
        e1n1FreeWallStart = (elementsCountAroundAtrialSeptum//2) - e1FreeWallStart
        for e2 in range(elementsCountUpAtria - 1):

            for i in range(2):

                if i == 0:
                    # left
                    aNodeId = laNodeId
                    bNodeId = raNodeId
                    apexNodeId = laApexNodeId
                else:
                    # right
                    aNodeId = raNodeId
                    bNodeId = laNodeId
                    apexNodeId = raApexNodeId
                e1Start = e1FreeWallStart + 1
                e1Limit = e1Start + elementsCountAroundAtria - elementsCountAroundAtrialSeptum - 2

                for e1 in range(elementsCountAroundAtria + 2):
                    eft1 = eft
                    elementtemplate1 = elementtemplate
                    nids = None
                    meshGroups = []

                    if e1 < e1FreeWallStart:
                        if i == 1:
                            continue  # already made in left atrium loop
                        if e2 > 0:
                            continue
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if (e1 == 0) or (e1 == (e1FreeWallStart - 1)):
                            meshGroups = [ laMeshGroup, raMeshGroup ]
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            if e1 == 0:
                                na1 = e1 + e1n1FreeWallStart + 2
                                nb1 = n1la2ra[na1]
                                nids = [ aNodeId[0][e2][na1], bNodeId[0][e2][nb1], aNodeId[0][e2 + 1][na1], bNodeId[0][e2 + 1][nb1], \
                                         aNodeId[1][e2][na1], bNodeId[1][e2][nb1], aNodeId[1][e2 + 1][na1], bNodeId[1][e2 + 1][nb1] ]
                            else:
                                na1 = e1 + e1n1FreeWallStart + 1
                                nb1 = n1la2ra[na1]
                                nids = [ bNodeId[0][e2][nb1], aNodeId[0][e2][na1], bNodeId[0][e2 + 1][nb1], aNodeId[0][e2 + 1][na1], \
                                         aNodeId[1][e2][na1], bNodeId[1][e2 + 1][nb1], aNodeId[1][e2 + 1][na1] ]
                                # 7 node crux collapsed in xi1 on xi2 == 0, xi3 == 1
                                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 7 ]
                                remapEftLocalNodes(eft1, 7, ln_map)
                        else:
                            meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup ]
                            na1 = e1 + e1n1FreeWallStart + 1
                            na2 = na1 + 1
                            nb1 = n1la2ra[na1]
                            nb2 = n1la2ra[na2]
                            nids = [ aNodeId[0][e2][na1], aNodeId[0][e2][na2], aNodeId[0][e2 + 1][na1], aNodeId[0][e2 + 1][na2], \
                                     bNodeId[0][e2][nb1], bNodeId[0][e2][nb2], bNodeId[0][e2 + 1][nb1], bNodeId[0][e2 + 1][nb2] ]
                            scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ] )
                            if e1 == 1:
                                # septum end 1
                                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            elif elementsCountAroundAtrialSeptum > 1:
                                scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            if e1 == (e1FreeWallStart - 2):
                                # septum end 2
                                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            elif elementsCountAroundAtrialSeptum > 1:
                                scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX
                    else:
                        # atrial free wall
                        if (e2 > 0) and ((e1 < e1Start) or (e1 >= e1Limit)):
                            continue
                        na1 = e1 + e1n1FreeWallStart
                        na2 = na1 + 1
                        nids = [ aNodeId[0][e2][na1], aNodeId[0][e2][na2], aNodeId[0][e2 + 1][na1], aNodeId[0][e2 + 1][na2], \
                                 aNodeId[1][e2][na1], aNodeId[1][e2][na2], aNodeId[1][e2 + 1][na1], aNodeId[1][e2 + 1][na2] ]
                        meshGroups = [ laMeshGroup if (i == 0) else raMeshGroup ]
                        if (e2 == 0) and (((i == 0) and (e1 == e1FreeWallStart)) or ((i == 1) and (e1 == elementsCountAroundAtria + 1))):
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            if i == 0:
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            else:
                                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if eft1.getNumberOfLocalScaleFactors() == 1:
                        result3 = element.setScaleFactors(eft1, [ -1.0 ])
                    else:
                        result3 = ' '
                    #print('create element', 'la' if i == 0 else 'ra', element.isValid(), elementIdentifier, result2, result3, nids)
                    if i == 0:
                        if e2 == (elementsCountUpAtria - 2):
                            if e1 == e1lspv:
                                lspvElementId = elementIdentifier
                            elif e1 == e1lipv:
                                lipvElementId = elementIdentifier
                    elif (e2 == 0) and (e1 == (e1Start - 1)):
                        ivcElementId1 = elementIdentifier
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        # create atria roof / apex elements
        n2 = elementsCountUpAtria - 1
        radiansPerElementAroundApex = math.pi/(lan1SeptumLimit - lan1CruxLimit)
        for i in range(2):
            if i == 0:
                # left
                aNodeId = laNodeId
                apexNodeId = laApexNodeId
                e1Start = lan1CruxLimit
                e1Limit = lan1SeptumLimit
                meshGroups = [ laMeshGroup ]
            else:
                # right
                aNodeId = raNodeId
                apexNodeId = raApexNodeId
                e1Start = n1la2ra[lan1SeptumLimit]
                e1Limit = n1la2ra[lan1CruxLimit]
                meshGroups = [ raMeshGroup ]
            aDerivativeMid = 0.5*(min(aDerivatives) + max(aDerivatives))

            # scale factor identifiers follow convention of offsetting by 100 for each 'version'
            radiansAroundApex = -0.5*math.pi if (i == 0) else 0.5*math.pi
            s = 0
            for e1 in range(e1Start, e1Limit):
                radiansAroundNext = radiansAroundApex + radiansPerElementAroundApex
                va = e1
                vb = e1 + 1
                eft1 = tricubichermite.createEftShellPoleTop(s*100, (s + 1)*100)
                elementtemplateX.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplateX)
                nodeIdentifiers = [ aNodeId[0][n2][va], aNodeId[0][n2][vb], apexNodeId[0], aNodeId[1][n2][va], aNodeId[1][n2][vb], apexNodeId[1] ]
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                scalefactors = [
                    -1.0,
                    math.cos(radiansAroundApex), math.sin(radiansAroundApex), radiansPerElementAroundApex,
                    math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAroundApex,
                    math.cos(radiansAroundApex), math.sin(radiansAroundApex), radiansPerElementAroundApex,
                    math.cos(radiansAroundNext), math.sin(radiansAroundNext), radiansPerElementAroundApex
                ]
                result = element.setScaleFactors(eft1, scalefactors)
                elementIdentifier = elementIdentifier + 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

                radiansAroundApex += radiansPerElementAroundApex
                s += 1

        # create septum arc side elements
        lna = elementsCountAroundAtrialSeptum//2 + 1
        rna = n1la2ra[lna]
        lnp = elementsCountAroundAtria - ((elementsCountAroundAtrialSeptum + 1)//2) - 1
        rnp = n1la2ra[lnp]
        for j in range(3):  # left, centre, right

            if j == 0:
                splNodeId = []
                salNodeId = []
                for n3 in range(2):
                    splNodeId.append([ laNodeId[n3][n2][lnp] for n2 in range(elementsCountUpAtria) ] + [ laApexNodeId[n3] ])
                    salNodeId.append([ laNodeId[n3][n2][lna] for n2 in range(elementsCountUpAtria) ] + [ laApexNodeId[n3] ])
                sprNodeId = laspNodeId
                sarNodeId = lasaNodeId
                meshGroups = [ laMeshGroup ]
            elif j == 1:
                splNodeId = laspNodeId
                salNodeId = lasaNodeId
                sprNodeId = raspNodeId
                sarNodeId = rasaNodeId
                meshGroups = [ laMeshGroup, raMeshGroup ]
            else:
                splNodeId = raspNodeId
                salNodeId = rasaNodeId
                sprNodeId = []
                sarNodeId = []
                for n3 in range(2):
                    sprNodeId.append([ raNodeId[n3][n2][rnp] for n2 in range(elementsCountUpAtria) ] + [ raApexNodeId[n3] ])
                    sarNodeId.append([ raNodeId[n3][n2][rna] for n2 in range(elementsCountUpAtria) ] + [ raApexNodeId[n3] ])
                meshGroups = [ raMeshGroup ]

            # posterior
            for e2 in range(1, elementsCountUpAtria):
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [ splNodeId[0][e2], sprNodeId[0][e2], splNodeId[0][e2 + 1], sprNodeId[0][e2 + 1],
                         splNodeId[1][e2], sprNodeId[1][e2], splNodeId[1][e2 + 1], sprNodeId[1][e2 + 1] ]

                if j == 1: # septum
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = ' '
                #print('create element sp arc', element.isValid(), elementIdentifier, result2, result3, nids)
                if (j == 0) and (e2 == (elementsCountUpAtria - 1)):
                    ripvElementId = elementIdentifier
                elif (j == 2) and (e2 == 1):
                    ivcElementId2 = elementIdentifier
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

            # anterior
            for e2 in range(1, elementsCountUpAtria):
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [ sarNodeId[0][e2], salNodeId[0][e2], sarNodeId[0][e2 + 1], salNodeId[0][e2 + 1],
                         sarNodeId[1][e2], salNodeId[1][e2], sarNodeId[1][e2 + 1], salNodeId[1][e2 + 1] ]

                if j == 1: # septum
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    if e2 == (elementsCountUpAtria - 1):
                        remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                        scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                    else:
                        remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                elif e2 == (elementsCountUpAtria - 1):
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = ' '
                #print('create element sa arc', element.isValid(), elementIdentifier, result2, result3, nids)
                if (j == 0) and (e2 == (elementsCountUpAtria - 1)):
                    rspvElementId = elementIdentifier
                elif (j == 2):
                    if e2 == (elementsCountUpAtria - 2):
                        svcElementId1 = elementIdentifier
                    elif e2 == (elementsCountUpAtria - 1):
                        svcElementId2 = elementIdentifier
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        # fossa ovalis elements

        meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup, fossaMeshGroup ]
        radiansAround0 = fossaRadiansAround[-1] - 2.0*math.pi
        for e1 in range(elementsCountAroundFossa):
            va = e1
            vb = (e1 + 1)%elementsCountAroundFossa
            eft1 = tricubichermite.createEftShellPoleTop(va*100, vb*100)
            elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            nids = [ fossaNodeId[0][va], fossaNodeId[0][vb], fossaCentreNodeId[0], fossaNodeId[1][va], fossaNodeId[1][vb], fossaCentreNodeId[1] ]
            result2 = element.setNodesByIdentifier(eft1, nids)
            radiansAround1 = fossaRadiansAround[va]
            radiansAround2 = fossaRadiansAround[vb]
            if radiansAround2 < radiansAround1:
                radiansAround2 += 2.0*math.pi
            radiansAround3 = fossaRadiansAround[vb + 1 - elementsCountAroundFossa]
            if radiansAround3 < radiansAround2:
                radiansAround3 += 2.0*math.pi
            dRadiansAround1 = 0.5*(radiansAround2 - radiansAround0)
            dRadiansAround2 = 0.5*(radiansAround3 - radiansAround1)
            #print('dRadians',dRadiansAround1,dRadiansAround2)
            scalefactors = [
                -1.0,
                -math.cos(radiansAround1), -math.sin(radiansAround1), dRadiansAround1,
                -math.cos(radiansAround2), -math.sin(radiansAround2), dRadiansAround2,
                -math.cos(radiansAround1), -math.sin(radiansAround1), dRadiansAround1,
                -math.cos(radiansAround2), -math.sin(radiansAround2), dRadiansAround2
            ]
            result3 = element.setScaleFactors(eft1, scalefactors)
            #print('create element fossa', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier = elementIdentifier + 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            radiansAround0 = radiansAround1

        # septum ring outside fossa ovalis

        meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup ]
        n1BaseOffset = -((elementsCountAroundAtrialSeptum + 1)//2)
        for e1 in range(elementsCountAroundFossa):
            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            nids = None

            na1 = (e1 + n1BaseOffset)%elementsCountAroundAtria
            na2 = (e1 + n1BaseOffset + 1)%elementsCountAroundAtria
            nb1 = n1la2ra[na1]
            nb2 = n1la2ra[na2]
            nf1 = (e1 + n1BaseOffset)%elementsCountAroundFossa
            nf2 = (e1 + n1BaseOffset + 1)%elementsCountAroundFossa

            if e1 < elementsCountAroundAtrialSeptum:
                nids = [ laNodeId[0][1][na1], laNodeId[0][1][na2], fossaNodeId[0][nf1], fossaNodeId[0][nf2], \
                         raNodeId[0][1][nb1], raNodeId[0][1][nb2], fossaNodeId[1][nf1], fossaNodeId[1][nf2] ]
                if e1 == 0:
                    # septum end 1
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                elif elementsCountAroundAtrialSeptum > 1:
                    scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                if e1 == (elementsCountAroundAtrialSeptum - 1):
                    # septum end 2
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                elif elementsCountAroundAtrialSeptum > 1:
                    scaleEftNodeValueLabels(eft1, [ 6 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
            elif e1 < (elementsCountAroundAtrialSeptum + elementsCountUpAtria - 1):
                ns1 = e1 - elementsCountAroundAtrialSeptum + 1
                ns2 = ns1 + 1
                nids = [ lasaNodeId[0][ns1], lasaNodeId[0][ns2], fossaNodeId[0][nf1], fossaNodeId[0][nf2], \
                         rasaNodeId[0][ns1], rasaNodeId[0][ns2], fossaNodeId[1][nf1], fossaNodeId[1][nf2] ]
                if ns2 < elementsCountUpAtria:
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                else:
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # finish swap
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            else:
                ns1 = elementsCountUpAtria - (e1 - elementsCountAroundAtrialSeptum - elementsCountUpAtria + 1)
                ns2 = ns1 - 1
                nids = [ laspNodeId[0][ns1], laspNodeId[0][ns2], fossaNodeId[0][nf1], fossaNodeId[0][nf2], \
                         raspNodeId[0][ns1], raspNodeId[0][ns2], fossaNodeId[1][nf1], fossaNodeId[1][nf2] ]
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # finish swap
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # finish swap
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])

            elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, [ -1.0 ])
            #print('create element septum', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier = elementIdentifier + 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # add left atria inlets (pulmonary veins)

        pvElementIds = [ lipvElementId, lspvElementId, ripvElementId, rspvElementId ]
        pvMeshGroups = [ lipvMeshGroup, lspvMeshGroup, ripvMeshGroup, rspvMeshGroup ]
        pvInnerRadii = [ lpvInnerRadius, lpvInnerRadius, rpvInnerRadius, rpvInnerRadius ]
        pvWallThickesses = [ lpvWallThickness, lpvWallThickness, rpvWallThickness, rpvWallThickness ]
        pvInclineRadians = [ -0.25*math.pi, -0.25*math.pi, 0.25*math.pi, -0.25*math.pi ]
        pvInclineAxis = [ 2, 2, 1, 1 ]

        for i in range(len(pvElementIds)):
            element = mesh.findElementByIdentifier(pvElementIds[i])
            pvLength = pvInnerRadii[i]*2.0
            tricubichermite.replaceElementWithInlet4(element, elementIdentifier, nodetemplate, nodeIdentifier, \
                pvLength, pvInnerRadii[i], pvWallThickesses[i], meshGroups = [ laMeshGroup, pvMeshGroups[i] ], \
                revCorners = ([ 3, 4 ] if (i == 3) else []), inclineRadians = pvInclineRadians[i], \
                inclineAxis = pvInclineAxis[i])
            elementIdentifier += 4
            nodeIdentifier += 8

        # add right atria inlets: inferior and superior vena cavae
        #print('ivc elements', ivcElementId1, ivcElementId2, '\nsvc elements', svcElementId1, svcElementId2)
        ivcElement1 = mesh.findElementByIdentifier(ivcElementId1)
        ivcElement2 = mesh.findElementByIdentifier(ivcElementId2)
        svcElement1 = mesh.findElementByIdentifier(svcElementId1)
        svcElement2 = mesh.findElementByIdentifier(svcElementId2)
        diff1 = mesh.getChartDifferentialoperator(1, 1)
        cache.setMeshLocation(ivcElement1, [ 0.5, 1.0, 1.0 ])
        resultx, ix = coordinates.evaluateReal(cache, 3)
        resultd, id = coordinates.evaluateDerivative(diff1, cache, 3)
        cache.setMeshLocation(svcElement1, [ 0.5, 1.0, 1.0 ])
        resultx, sx = coordinates.evaluateReal(cache, 3)
        resultd, sd = coordinates.evaluateDerivative(diff1, cache, 3)
        direction = vector.normalise([ (ix[c] - sx[c]) for c in range(3) ])

        ivcInletScale = -ivcInnerRadius
        ivcInletAxis = [ d*ivcInletScale for d in direction ]
        n = vector.crossproduct3(id, ivcInletAxis)
        ivcInletSide = vector.crossproduct3(n, ivcInletAxis)
        ivcInletDistance = ivcInnerRadius
        ivcInletCentre = [ (ix[c] + direction[c]*ivcInletDistance) for c in range(3) ]
        tricubichermite.replaceTwoElementWithInlet6(ivcElement1, ivcElement2, elementIdentifier, nodetemplate, nodeIdentifier, \
            ivcInletCentre, ivcInletAxis, ivcInletSide, ivcInnerRadius, ivcWallThickness, meshGroups = [ raMeshGroup, ivcInletMeshGroup ])
        elementIdentifier += 6
        nodeIdentifier += 12

        svcInletScale = -svcInnerRadius
        svcInletAxis = [ -d*svcInletScale for d in direction ]
        n = vector.crossproduct3(sd, svcInletAxis)
        svcInletSide = vector.crossproduct3(n, svcInletAxis)
        svcInletDistance = svcInnerRadius
        svcInletCentre = [ (sx[c] - direction[c]*svcInletDistance) for c in range(3) ]
        tricubichermite.replaceTwoElementWithInlet6(svcElement1, svcElement2, elementIdentifier, nodetemplate, nodeIdentifier, \
            svcInletCentre, svcInletAxis, svcInletSide, svcInnerRadius, svcWallThickness, meshGroups = [ raMeshGroup, svcInletMeshGroup ], revCorners = [ 5, 6 ])
        elementIdentifier += 6
        nodeIdentifier += 12

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
        refineElementsCountThroughAtrialWall = options['Refine number of elements through atrial wall']
        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            meshrefinement.refineElementCubeStandard3d(element, refineElementsCountSurface, refineElementsCountSurface, refineElementsCountThroughAtrialWall)
            element = meshrefinement._sourceElementiterator.next()
