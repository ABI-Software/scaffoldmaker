"""
Generates a 3-D heart atria model, suitable for attachment to the
3-D Heart Ventricles with Base 2.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartatria2(object):
    '''
    3-D heart atria model, suitable for attachment to the 3-D Heart Ventricles with Base 2.
    '''

    @staticmethod
    def getName():
        return '3D Heart Atria 2'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around atria' : 8,
            'Number of elements around atrial septum' : 2,
            'Number of elements up atria' : 3,
            'Number of elements radial atrial septum' : 2,
            'Atria base inner major axis length' : 0.5,  # 0.549
            'Atria base inner minor axis length' : 0.35,  # 0.37
            'Atria major axis rotation degrees' : 40.0,
            'Atria outer height' : 0.4,
            'Atria outer major arc up degrees' : 100.0,
            'Atria outer minor arc up degrees' : 120.0,
            'Atrial crux side offset' : 0.195*math.sin(math.pi/3.0), # from ventriclesbase2 LV outlet inner radius + wall thickness
            'Atrial septum thickness' : 0.06,
            'Atrial free wall thickness' : 0.02,
            'Atrial base wall thickness' : 0.05,
            'Atrial base slope degrees' : 30.0,
            'Left pulmonary vein inner diameter' : 0.07,
            'Left pulmonary vein wall thickness' : 0.007,
            'Right pulmonary vein inner diameter' : 0.085,
            'Right pulmonary vein wall thickness' : 0.0085,
            'Inferior vena cava inner diameter' : 0.17,
            'Inferior vena cava wall thickness' : 0.017,
            'Superior vena cava inner diameter' : 0.15,
            'Superior vena cava wall thickness' : 0.015,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around atria',
            'Number of elements around atrial septum',
            'Number of elements up atria',
            'Number of elements radial atrial septum',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atria outer height',
            'Atria outer major arc up degrees',
            'Atria outer minor arc up degrees',
            'Atrial crux side offset',
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
            'Superior vena cava wall thickness',
            #,'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if options['Number of elements around atria'] < 6:
            options['Number of elements around atria'] = 6
        # need even number of elements around for topology
        if (options['Number of elements around atria'] % 2) == 1:
            options['Number of elements around atria'] += 1
        if options['Number of elements around atrial septum'] < 1:
            options['Number of elements around atrial septum'] = 1
        elif options['Number of elements around atrial septum'] > (options['Number of elements around atria'] - 4):
            options['Number of elements around atrial septum'] = options['Number of elements around atria'] - 4
        if options['Number of elements up atria'] < 2:
            options['Number of elements up atria'] = 2
        if options['Number of elements radial atrial septum'] < 1:
            options['Number of elements radial atrial septum'] = 1
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
        if options['Atrial crux side offset'] < 0.55*options['Atrial septum thickness']:
            options['Atrial crux side offset'] = 0.55*options['Atrial septum thickness']
        for key in [
            'Atria major axis rotation degrees']:
            if options[key] < -75.0:
                options[key] = -75.0
            elif options[key] > 75.0:
                options[key] = 75.0

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundAtria = options['Number of elements around atria']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountUpAtria = options['Number of elements up atria']
        elementsCountRadialAtrialSeptum = options['Number of elements radial atrial septum']
        aBaseInnerMajorMag = 0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = 0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aOuterMajorArcUpRadians = math.radians(options['Atria outer major arc up degrees'])
        aOuterMinorArcUpRadians = math.radians(options['Atria outer minor arc up degrees'])
        aOuterHeight = options['Atria outer height']
        aCruxSideOffset = options['Atrial crux side offset']
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
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        laGroup = AnnotationGroup(region, 'left atrium', FMANumber = 7097, lyphID = 'Lyph ID unknown')
        raGroup = AnnotationGroup(region, 'right atrium', FMANumber = 7096, lyphID = 'Lyph ID unknown')
        aSeptumGroup = AnnotationGroup(region, 'interatrial septum', FMANumber = 7108, lyphID = 'Lyph ID unknown')
        fossaGroup = AnnotationGroup(region, 'fossa ovalis', FMANumber = 9246, lyphID = 'Lyph ID unknown')
        lipvGroup = AnnotationGroup(region, 'left inferior pulmonary vein', FMANumber = 49913, lyphID = 'Lyph ID unknown')
        lspvGroup = AnnotationGroup(region, 'left superior pulmonary vein', FMANumber = 49916, lyphID = 'Lyph ID unknown')
        ripvGroup = AnnotationGroup(region, 'right inferior pulmonary vein', FMANumber = 49911, lyphID = 'Lyph ID unknown')
        rspvGroup = AnnotationGroup(region, 'right superior pulmonary vein', FMANumber = 49914, lyphID = 'Lyph ID unknown')
        ivcInletGroup = AnnotationGroup(region, 'inferior vena cava inlet', FMANumber = 10951, lyphID = 'Lyph ID unknown')
        svcInletGroup = AnnotationGroup(region, 'superior vena cava inlet', FMANumber = 4720, lyphID = 'Lyph ID unknown')
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

        nodeIdentifier = 1

        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)
        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)

        # convert base major/minor sizes to equator using arc up

        aEquatorOuterMajorMag = (aBaseInnerMajorMag + aBaseSlopeLength)/math.sin(aOuterMajorArcUpRadians)
        aEquatorOuterMinorMag = (aBaseInnerMinorMag + aBaseSlopeLength)/math.sin(aOuterMinorArcUpRadians)
        aEquatorInnerMajorMag = aEquatorOuterMajorMag - aFreeWallThickness
        aEquatorInnerMinorMag = aEquatorOuterMinorMag - aFreeWallThickness

        # following are angles in radians around LA ellipse from major axis
        axInner = aBaseInnerMajorMag*math.cos(aMajorAxisRadians)
        bxInner = aBaseInnerMinorMag*math.sin(aMajorAxisRadians)
        laSeptumRadians = math.atan2(bxInner, axInner)
        laCentreX = -0.5*aSeptumThickness - axInner*math.cos(laSeptumRadians) - bxInner*math.sin(laSeptumRadians)
        #laCruxLeftRadians = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, laSeptumRadians, \
        #    (aSeptumBaseLength/elementsCountAroundAtrialSeptum)*(0.5*elementsCountAroundAtrialSeptum + 1.0))
        # calculate crux left from aCruxSideOffset = (lvOutletInnerRadius + lvOutletWallThickness)*math.sin(math.pi/3.0)
        axOuter = (aBaseInnerMajorMag + aBaseSlopeLength)*math.cos(aMajorAxisRadians)
        bxOuter = (aBaseInnerMinorMag + aBaseSlopeLength)*math.sin(aMajorAxisRadians)
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

        aMinorToMajorRadians = aOuterMajorArcUpRadians/aOuterMinorArcUpRadians

        laNodeId = [ [], [] ]
        raNodeId = [ [], [] ]
        laApexNodeId = []
        raApexNodeId = []

        n1SeptumRange = []
        for n1 in range(elementsCountAroundAtria):
            if (n1 < (elementsCountAroundAtrialSeptum//2)) or (n1 > (elementsCountAroundAtria - ((elementsCountAroundAtrialSeptum + 1)//2))):
                n1SeptumRange.append(n1)
        #print('n1SeptumRange', n1SeptumRange)

        # compute radians around based on base outer major and minor axis sizes
        aBaseOuterMajorMag = aBaseInnerMajorMag + aBaseSlopeLength
        aBaseOuterMinorMag = aBaseInnerMinorMag + aBaseSlopeLength
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
            n1la2ra = [ -n1 for n1 in range(elementsCountAroundAtria) ]
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
                    else:
                        # right
                        centreX = raCentreX
                        centreY = raCentreY
                        majorAxisRadians = math.pi + aMajorAxisRadians
                        aRadians = raRadians
                        aDerivatives = raDerivatives
                        aLayerNodeId = raLayerNodeId

                    sinMajorAxisRadians = math.sin(majorAxisRadians)
                    cosMajorAxisRadians = math.cos(majorAxisRadians)
                    aMajorX =  aMajorMag*cosMajorAxisRadians
                    aMajorY =  aMajorMag*sinMajorAxisRadians
                    aMinorX = -aMinorMag*sinMajorAxisRadians
                    aMinorY =  aMinorMag*cosMajorAxisRadians
                    #print('aMajor', aMajorX, aMajorY,'aMinor',aMinorX,aMinorY)

                    for n1 in range(elementsCountAroundAtria):
                        radiansAround = aRadians[n1]
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        if (n3 == 1) and (n1 in n1SeptumRange):
                            continue  # outer septum node created by other side
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        aLayerNodeId[n1] = nodeIdentifier
                        cache.setNode(node)
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

                        result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                        # derivative 3 is set later
                        nodeIdentifier += 1

            # apex nodes
            for i in range(2):
                if i == 0:
                    # left
                    centreX = laCentreX
                    centreY = laCentreY
                    apexNodeId = laApexNodeId
                    majorAxisRadians = -aMajorAxisRadians
                else:
                    # right
                    centreX = raCentreX
                    centreY = raCentreY
                    apexNodeId = raApexNodeId
                    majorAxisRadians = math.pi + aMajorAxisRadians
                sinMajorAxisRadians = math.sin(majorAxisRadians)
                cosMajorAxisRadians = math.cos(majorAxisRadians)

                x = [ centreX, centreY, aOuterHeight - aFreeWallThickness if (n3 == 0) else aOuterHeight ]
                dx_ds1 = [ apexDerivativeUpMajor*cosMajorAxisRadians, apexDerivativeUpMajor*sinMajorAxisRadians, 0.0 ]
                dx_ds2 = [ -elementSizeUpMinor*sinMajorAxisRadians, elementSizeUpMinor*cosMajorAxisRadians, 0.0 ]
                dx_ds3 = [ 0.0, 0.0, aFreeWallThickness ]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                apexNodeId.append(nodeIdentifier)
                cache.setNode(node)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                nodeIdentifier += 1

        # fix inner base derivative 2 to fit incline
        for i in range(2):
            aNodeId = laNodeId[0] if (i == 0) else raNodeId[0]
            for n1 in range(elementsCountAroundAtria):
                node1 = nodes.findNodeByIdentifier(aNodeId[0][n1])
                node2 = nodes.findNodeByIdentifier(aNodeId[1][n1])
                dx_ds2 = computeNodeDerivativeHermiteLagrange(cache, coordinates, node2, Node.VALUE_LABEL_D_DS2, -1.0, node1, -1.0)
                cache.setNode(node1)
                result = coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)

        # transfer inner septum nodes to outer on opposite side, set derivative 3 to be node difference
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

        #################
        # Create elements
        #################

        elementIdentifier = 1

        mesh = fm.findMeshByDimension(3)

        laMeshGroup = laGroup.getMeshGroup(mesh)
        raMeshGroup = raGroup.getMeshGroup(mesh)
        aSeptumMeshGroup = aSeptumGroup.getMeshGroup(mesh)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        e1FreeWallStart = elementsCountAroundAtrialSeptum + 2
        e1n1FreeWallStart = (elementsCountAroundAtrialSeptum//2) - e1FreeWallStart
        for e2 in range(elementsCountUpAtria - 1):

            for i in range(2):

                aNodeId, bNodeId = ( laNodeId, raNodeId ) if (i == 0) else ( raNodeId, laNodeId )

                for e1 in range(elementsCountAroundAtria + 2):
                    eft1 = eft
                    elementtemplate1 = elementtemplate
                    nids = None

                    if e1 < e1FreeWallStart:
                        if i == 1:
                            continue  # already made in left atrium loop
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        if (e1 == 0) or (e1 == (e1FreeWallStart - 1)):
                            meshGroups = [ laMeshGroup, raMeshGroup ]
                            if e1 == 0:
                                na1 = e1 + e1n1FreeWallStart + 2
                                nb1 = n1la2ra[na1]
                                nids = [ aNodeId[0][e2][na1], bNodeId[0][e2][nb1], aNodeId[0][e2 + 1][na1], bNodeId[0][e2 + 1][nb1], \
                                         aNodeId[1][e2][na1], bNodeId[1][e2][nb1], aNodeId[1][e2 + 1][na1], bNodeId[1][e2 + 1][nb1] ]
                            else:
                                na1 = e1 + e1n1FreeWallStart + 1
                                nb1 = n1la2ra[na1]
                                nids = [ bNodeId[0][e2][nb1], aNodeId[0][e2][na1], bNodeId[0][e2 + 1][nb1], aNodeId[0][e2 + 1][na1], \
                                         bNodeId[1][e2][nb1], aNodeId[1][e2][na1], bNodeId[1][e2 + 1][nb1], aNodeId[1][e2 + 1][na1] ]
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        else:
                            meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup ]
                            na1 = e1 + e1n1FreeWallStart + 1
                            na2 = na1 + 1
                            nb1 = n1la2ra[na1]
                            nb2 = n1la2ra[na2]
                            nids = [ aNodeId[0][e2][na1], aNodeId[0][e2][na2], aNodeId[0][e2 + 1][na1], aNodeId[0][e2 + 1][na2], \
                                     bNodeId[0][e2][nb1], bNodeId[0][e2][nb2], bNodeId[0][e2 + 1][nb1], bNodeId[0][e2 + 1][nb2] ]
                            scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ] )
                            if e1 == 1:
                                # septum end 1
                                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            elif elementsCountAroundAtrialSeptum > 1:
                                scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                            if e1 == (e1FreeWallStart - 2):
                                # septum end 2
                                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            elif elementsCountAroundAtrialSeptum > 1:
                                scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ] )
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX
                    else:
                        # atrial free wall
                        meshGroups = [ laMeshGroup ] if (i == 0) else [ raMeshGroup ]
                        na1 = e1 + e1n1FreeWallStart
                        na2 = na1 + 1
                        nids = [ aNodeId[0][e2][na1], aNodeId[0][e2][na2], aNodeId[0][e2 + 1][na1], aNodeId[0][e2 + 1][na2], \
                                 aNodeId[1][e2][na1], aNodeId[1][e2][na2], aNodeId[1][e2 + 1][na1], aNodeId[1][e2 + 1][na2] ]

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if eft1.getNumberOfLocalScaleFactors() == 1:
                        result3 = element.setScaleFactors(eft1, [ -1.0 ])
                    else:
                        result3 = ' '
                    #print('create element', 'la' if i == 0 else 'ra', element.isValid(), elementIdentifier, result2, result3, nids)
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        # create atria roof / apex elements
        n2 = elementsCountUpAtria - 1
        for i in range(2):
            if i == 0:
                # left
                nodeId = laNodeId
                apexNodeId = laApexNodeId
                aRadians = laRadians
                aDerivatives = laDerivatives
            else:
                # right
                nodeId = raNodeId
                apexNodeId = raApexNodeId
                aRadians = raRadians
                aDerivatives = raDerivatives
            aDerivativeMid = 0.5*(min(aDerivatives) + max(aDerivatives))

            # scale factor identifiers follow convention of offsetting by 100 for each 'version'
            for e1 in range(elementsCountAroundAtria):
                va = e1
                vb = (e1 + 1)%elementsCountAroundAtria
                eft1 = tricubichermite.createEftShellApexTop(va*100, vb*100)
                elementtemplateX.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplateX)
                nodeIdentifiers = [ nodeId[0][n2][va], nodeId[0][n2][vb], apexNodeId[0], nodeId[1][n2][va], nodeId[1][n2][vb], apexNodeId[1] ]
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                deltaRadiansPrev = aRadians[va] - aRadians[va - 1]
                deltaRadiansCurr = aRadians[vb] - aRadians[va]
                deltaRadiansNext = aRadians[(e1 + 2)%elementsCountAroundAtria] - aRadians[vb]
                if aDerivatives[va] < aDerivativeMid:
                    aDeltaRadians = min(deltaRadiansPrev, deltaRadiansCurr)
                else:
                    aDeltaRadians = max(deltaRadiansPrev, deltaRadiansCurr)
                if aDerivatives[vb] < aDerivativeMid:
                    bDeltaRadians = min(deltaRadiansCurr, deltaRadiansNext)
                else:
                    bDeltaRadians = max(deltaRadiansCurr, deltaRadiansNext)
                scalefactors = [
                    -1.0,
                    -math.cos(aRadians[va]), -math.sin(aRadians[va]), aDeltaRadians,
                    -math.cos(aRadians[vb]), -math.sin(aRadians[vb]), bDeltaRadians,
                    -math.cos(aRadians[va]), -math.sin(aRadians[va]), aDeltaRadians,
                    -math.cos(aRadians[vb]), -math.sin(aRadians[vb]), bDeltaRadians
                ]
                result = element.setScaleFactors(eft1, scalefactors)
                elementIdentifier = elementIdentifier + 1


        fm.endChange()
        return annotationGroups

