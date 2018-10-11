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

        # LV nodes

        # get distance from outside of round outer LV to outside of RV base centre
        rvAddWidthBase = rvWidth - vSeptumBaseRadialDisplacement + vSeptumThickness - lvFreeWallThickness + rvFreeWallThickness

        # get element spacing up
        rvThetaUp = math.acos(rvHeightFraction)
        rvOuterHeight = 0.5*lvOuterHeight*(1.0 + rvHeightFraction)
        rvCrossHeight = rvHeightFraction*lvOuterHeight
        radialDisplacementStartRadiansUp = math.acos(rvOuterHeight/lvOuterHeight)

        # get LV inner points

        lvInnerx  = []
        lvInnerd1 = []
        lvInnerHeight = lvOuterHeight - lvApexThickness
        lvInnerRadius = lvOuterRadius - lvFreeWallThickness
        lengthUp = 0.25*getApproximateEllipsePerimeter(lvInnerHeight, lvInnerRadius)
        lengthUpApex = getEllipseArcLength(lvInnerHeight, lvInnerRadius, 0.0, rvThetaUp)
        lengthUpRV = lengthUp - lengthUpApex
        elementSizeUpRV = lengthUpRV/elementsCountUpRV
        elementSizeUpApex = (lengthUpApex - 0.5*elementSizeUpRV)/(elementsCountUpLVApex - 0.5)
        elementSizeUpApexTransition = 0.5*(elementSizeUpApex + elementSizeUpRV)
        # apex points, noting s1, s2 is x, -y to get out outward s3
        vApexInnerx = [ 0.0, 0.0, -lvInnerHeight ]
        vApexInnerd1 = [ elementSizeUpApex, 0.0, 0.0 ]
        vApexInnerd2 = [ 0.0, -elementSizeUpApex, 0.0, 0.0 ]
        vApexInnerd3 = [ 0.0, 0.0, -lvApexThickness ]

        lvInnerRadiansUp = []
        radiansUp = 0.0
        for n2 in range(elementsCountUpLV):
            if n2 < elementsCountUpLVApex:
                arcLengthUp = elementSizeUpApex if (n2 < (elementsCountUpLVApex - 1)) else elementSizeUpApexTransition
            else:
                arcLengthUp = elementSizeUpRV
            radiansUp = updateEllipseAngleByArcLength(lvInnerHeight, lvInnerRadius, radiansUp, arcLengthUp)
            if n2 == (elementsCountUpLV - 1):
                radiansUp = 0.5*math.pi
            lvInnerRadiansUp.append(radiansUp)

        for n2 in range(elementsCountUpLV):
            radiansUp = lvInnerRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            z = -lvInnerHeight*cosRadiansUp
            lvRadius = lvInnerRadius*sinRadiansUp
            rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvInnerHeight)

            # get radial displacement of centre of septum, a function of radiansUp
            xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
            midSeptumDisplacement = interpolateCubicHermite([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]

            nx, nd1 = getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, rvArcAroundZRadians, z, \
                elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, ivSulcusDerivativeFactor)

            lvInnerx.append(nx)
            lvInnerd1.append(nd1)

        # get ventricles outer points

        vOuterx  = []
        vOuterd1 = []
        lengthUp = 0.25*getApproximateEllipsePerimeter(lvOuterHeight, lvOuterRadius)
        lengthUpApex = getEllipseArcLength(lvOuterHeight, lvOuterRadius, 0.0, rvThetaUp)
        lengthUpRV = lengthUp - lengthUpApex
        elementSizeUpRV = lengthUpRV/elementsCountUpRV
        elementSizeUpApex = (lengthUpApex - 0.5*elementSizeUpRV)/(elementsCountUpLVApex - 0.5)
        elementSizeUpApexTransition = 0.5*(elementSizeUpApex + elementSizeUpRV)
        # apex points, noting s1, s2 is x, -y to get out outward s3
        vApexOuterx = [ 0.0, 0.0, -lvOuterHeight ]
        vApexOuterd1 = [ elementSizeUpApex, 0.0, 0.0 ]
        vApexOuterd2 = [ 0.0, -elementSizeUpApex, 0.0, 0.0 ]
        vApexOuterd3 = [ 0.0, 0.0, -lvApexThickness ]

        lvOuterRadiansUp = []
        radiansUp = 0.0
        for n2 in range(elementsCountUpLV):
            if n2 < elementsCountUpLVApex:
                arcLengthUp = elementSizeUpApex if (n2 < (elementsCountUpLVApex - 1)) else elementSizeUpApexTransition
            else:
                arcLengthUp = elementSizeUpRV
            radiansUp = updateEllipseAngleByArcLength(lvOuterHeight, lvOuterRadius, radiansUp, arcLengthUp)
            if n2 == (elementsCountUpLV - 1):
                radiansUp = 0.5*math.pi
            lvOuterRadiansUp.append(radiansUp)

        for n2 in range(elementsCountUpLV):
            radiansUp = lvOuterRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            z = -lvOuterHeight*cosRadiansUp
            lvRadius = lvOuterRadius*sinRadiansUp

            rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvOuterHeight)

            # get size of RV at z
            xiUpZ = 1.0 + z/rvOuterHeight
            xiUpCross = 1.0 + z/rvCrossHeight
            rvAddWidthRadius, rvAddCrossRadius = getRVOuterSize(xiUpZ, xiUpCross, rvAddWidthBase, rvExtraCrossRadiusBase)

            nx, nd1 = getVentriclesOuterPoints(lvRadius, rvAddWidthRadius, rvAddCrossRadius, rvArcAroundZRadians, z, \
                elementsCountAroundLVFreeWall, elementsCountAroundRVFreeWall, ivSulcusDerivativeFactor)

            vOuterx.append(nx)
            vOuterd1.append(nd1)

        # get radians around apex for each radial line on outer surface
        rvApexArcAroundRadians = rvArcAroundBaseRadians*rvArcApexFraction
        lvApexArcAroundRadians = 2.0*math.pi - rvApexArcAroundRadians
        lvApexRadiansPerElement = lvApexArcAroundRadians/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
        lvApexRadiansPerElementTransition = 0.5*lvApexRadiansPerElement*(1.0 + ivSulcusDerivativeFactor)
        ivSulcusRadiansPerElement = 2.0*lvApexRadiansPerElementTransition - lvApexRadiansPerElement
        rvApexRadiansPerElement = (rvApexArcAroundRadians - ivSulcusRadiansPerElement)/(elementsCountAroundRVFreeWall - 1)
        rvApexRadiansPerElementTransition = 0.5*(rvApexRadiansPerElement + ivSulcusRadiansPerElement)
        apexRadiansAround = []
        radiansAround = -0.5*rvApexArcAroundRadians
        for n1 in range(elementsCountAroundLVFreeWall):
            apexRadiansAround.append(radiansAround)
            if (n1 == 0) or (n1 == (elementsCountAroundLVFreeWall - 1)):
                radiansAround -= lvApexRadiansPerElementTransition
            else:
                radiansAround -= lvApexRadiansPerElement
        for n1 in range(elementsCountAroundRVFreeWall):
            apexRadiansAround.append(radiansAround)
            if (n1 == 0) or (n1 == (elementsCountAroundRVFreeWall - 1)):
                radiansAround -= rvApexRadiansPerElementTransition
            else:
                radiansAround -= rvApexRadiansPerElement

        # calculate derivative 2 on inner surfaces

        lvInnerd2 = []
        for n2 in range(elementsCountUpLV):
            lvInnerd2.append([])
        for n1 in range(pointsCountAroundOuter):
            cosRadiansAround = math.cos(apexRadiansAround[n1])
            sinRadiansAround = math.sin(apexRadiansAround[n1])
            apexd2 = [ (vApexOuterd1[c]*cosRadiansAround + vApexOuterd2[c]*sinRadiansAround) for c in range(3) ]
            nx = [ vApexOuterx ]
            nd2 = [ apexd2 ]
            for n2 in range(elementsCountUpLV):
                nx.append(lvInnerx[n2][n1])
                nd2.append([ 0.0, 0.0, 0.0 ])
            nd2[-1] = [ 0.0, 0.0, elementSizeUpRV ]
            nd2 = smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True)
            for n2 in range(elementsCountUpLV):
                lvInnerd2[n2].append(nd2[n2 + 1])

        # calculate derivative 2 on outer surface

        vOuterd2 = []
        for n2 in range(elementsCountUpLV):
            vOuterd2.append([])
        for n1 in range(pointsCountAroundOuter):
            cosRadiansAround = math.cos(apexRadiansAround[n1])
            sinRadiansAround = math.sin(apexRadiansAround[n1])
            apexd2 = [ (vApexOuterd1[c]*cosRadiansAround + vApexOuterd2[c]*sinRadiansAround) for c in range(3) ]
            nx = [ vApexOuterx ]
            nd2 = [ apexd2 ]
            for n2 in range(elementsCountUpLV):
                nx.append(vOuterx[n2][n1])
                nd2.append([ 0.0, 0.0, 0.0 ])
            nd2[-1] = [ 0.0, 0.0, elementSizeUpRV ]
            nd2 = smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True)
            for n2 in range(elementsCountUpLV):
                vOuterd2[n2].append(nd2[n2 + 1])

        # calculate derivative 3 on LV free wall and sub-RV apex from difference between inner and outer surfaces

        lvInnerd3 = []
        vOuterd3 = []
        for n2 in range(elementsCountUpLV):
            layerInnerd3 = []
            layerOuterd3 = []
            for n1 in range(pointsCountAroundOuter):
                if (n2 < elementsCountUpLVApex) or (n1 <= elementsCountAroundLVFreeWall):
                    ix = lvInnerx[n2][n1]
                    ox = vOuterx[n2][n1]
                    innerd3 = outerd3 = [ (ox[c] - ix[c]) for c in range(3) ]
                else:
                    innerd3 = outerd3 = [ 0.0, 0.0, 0.0 ]
                layerInnerd3.append(innerd3)
                layerOuterd3.append(outerd3)
            lvInnerd3.append(layerInnerd3)
            vOuterd3.append(layerOuterd3)

        # get points on inside of RV, anticlockwise starting on posterior free wall

        septumOuterRadius = lvInnerRadius + vSeptumThickness
        rvInnerx = []
        rvInnerd1 = []
        rvInnerd2 = []
        rvInnerd3 = []
        for n2 in range(elementsCountUpLV):
            layerInnerx = []
            layerInnerd1 = []
            layerInnerd2 = []
            layerInnerd3 = []
            if n2 >= elementsCountUpLVApex:
                for n1 in range(elementsCountAroundLVFreeWall + 1, pointsCountAroundOuter):
                    outerx  = vOuterx [n2][n1]
                    outerd1 = vOuterd1[n2][n1]
                    outerd2 = vOuterd2[n2][n1]
                    unitNormal = vector.normalise(vector.crossproduct3(outerd1, outerd2))
                    innerd3 = vOuterd3[n2][n1] = vector.setMagnitude(unitNormal, rvFreeWallThickness)
                    innerx = [ (outerx[c] - innerd3[c]) for c in range(3) ]
                    # calculate inner d1 from curvature around
                    n1m = n1 - 1
                    n1p = (n1 + 1) % pointsCountAroundOuter
                    curvature = -0.5*(
                        getCubicHermiteCurvature(vOuterx[n2][n1m], vOuterd1[n2][n1m], outerx, outerd1, unitNormal, 1.0) +
                        getCubicHermiteCurvature(outerx, outerd1, vOuterx[n2][n1p], vOuterd1[n2][n1p], unitNormal, 0.0))
                    factor = 1.0 - curvature*rvFreeWallThickness
                    innerd1 = [ factor*c for c in outerd1 ]
                    # calculate inner d2 from curvature up
                    n2m = n2 - 1
                    curvature = -getCubicHermiteCurvature(vOuterx[n2m][n1], vOuterd2[n2m][n1], outerx, outerd2, unitNormal, 1.0)
                    if n2 < (elementsCountUpLVApex - 1):
                        n2p = n2 + 1
                        curvature = 0.5*(curvature - getCubicHermiteCurvature(outerx, outerd2, vOuterx[n2p][n1], vOuterd2[n2p][n1], unitNormal, 0.0))
                    factor = 1.0 - curvature*rvFreeWallThickness
                    innerd2 = [ factor*c for c in outerd2 ]
                    layerInnerx .append(innerx)
                    layerInnerd1.append(innerd1)
                    layerInnerd2.append(innerd2)
                    layerInnerd3.append(innerd3)
                # sample points onto outer septum
                innerRadiansUp = lvInnerRadiansUp[n2]
                z = -lvInnerHeight*math.cos(innerRadiansUp)
                radiansUp = math.acos(-z/lvOuterHeight)
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                lvRadius = septumOuterRadius*sinRadiansUp
                rvArcAroundZRadians = rvArcAroundBaseRadians*(1.0 + (1.0 - rvArcApexFraction)*z/lvOuterHeight)

                # get radial displacement of centre of septum, a function of radiansUp
                xiUp = max(0.0, (radiansUp - radialDisplacementStartRadiansUp)/(0.5*math.pi - radialDisplacementStartRadiansUp))
                midSeptumDisplacement = interpolateCubicHermite([0.0], [0.0], [vSeptumBaseRadialDisplacement], [0.0], xiUp)[0]

                numberAround = 12
                nx, nd1 = getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, rvArcAroundZRadians, z, numberAround, numberAround)
                # extract septum points, reversing order and derivative
                nx = [ layerInnerx [-1] ] + nx [-1:-numberAround:-1] + [ layerInnerx [1 - elementsCountAroundRVFreeWall] ]
                nd1 = [ layerInnerd1[-1] ] + [ ([-d for d in nd1[n1] ]) for n1 in range(-1, -numberAround, -1) ] + [ layerInnerd1[1 - elementsCountAroundRVFreeWall] ]
                px, pd1, _ = sampleCubicHermiteCurves(nx, nd1, [], elementsCountAroundVSeptum + 2,
                    addLengthStart = 0.5*vector.magnitude(nd1[ 0]), lengthFractionStart = 0.5,
                    addLengthEnd = 0.5*vector.magnitude(nd1[-1]), lengthFractionEnd = 0.5,
                    arcLengthDerivatives = False)
                for n1 in range(elementsCountAroundVSeptum + 1):
                    lvx = lvInnerx[n2][-n1]
                    rvx = px[n1 + 1]
                    innerd3 = [ (lvx[c] - rvx[c]) for c in range(3) ]
                    layerInnerx .append(rvx)
                    layerInnerd1.append(pd1[n1 + 1])
                    layerInnerd2.append([ 0.0, 0.0, 0.0 ])
                    layerInnerd3.append(innerd3)
                    if (n1 > 0) and (n1 < elementsCountAroundVSeptum):
                        lvInnerd3[n2][-n1] = [- d for d in innerd3 ]
            rvInnerx.append(layerInnerx)
            rvInnerd1.append(layerInnerd1)
            rvInnerd2.append(layerInnerd2)
            rvInnerd3.append(layerInnerd3)

        # calculate derivative 2 on inner RV septum
        n2Range = range(elementsCountUpLVApex, elementsCountUpLV)
        for n1 in range(elementsCountAroundRVFreeWall - 1, elementsCountAroundRVFreeWall + elementsCountAroundVSeptum):
            nx = [ rvInnerx[n2][n1] for n2 in n2Range ]
            nd2 = [ rvInnerd2[n2][n1] for n2 in n2Range ]
            nd2 = smoothCubicHermiteDerivativesLine([ rvInnerx[n2][n1] for n2 in n2Range ], [ rvInnerd2[n2][n1] for n2 in n2Range ])
            for n2 in n2Range:
                rvInnerd2[n2][n1] = nd2[n2 - elementsCountUpLVApex]

        # get points on RV apex

        n2 = elementsCountUpLVApex
        dFactor = 2.0
        sx = []
        sd1 = []
        sd2 = []
        for n1 in range(elementsCountAroundRVFreeWall - 1):
            r1 = 2*elementsCountAroundRVFreeWall - n1 - 2
            ax = rvInnerx[n2][r1]
            ad1 = [ -d for d in rvInnerd1[n2][r1] ]
            ad2 = [ -dFactor*d for d in rvInnerd2[n2][r1] ]
            bx = rvInnerx[n2][n1]
            bd1 = rvInnerd1[n2][n1]
            bd2 = [ dFactor*d for d in rvInnerd2[n2][n1] ]
            px, pd2, ( pd1, ) = sampleCubicHermiteCurves([ ax, bx ], [ ad2, bd2 ], [ [ ad1, bd1 ] ], elementsCountOut = 2,
                addLengthStart = 0.5*vector.magnitude(ad2)/dFactor, lengthFractionStart = 0.5,
                addLengthEnd = 0.5*vector.magnitude(bd2)/dFactor, lengthFractionEnd = 0.5, arcLengthDerivatives = False)
            sx .append(px[1])
            sd1.append(pd1[1])
            sd2.append(pd2[1])
        n1 = 2*elementsCountAroundRVFreeWall - 1
        ax  = rvInnerx [n2][n1]
        ad1 = rvInnerd1[n2][n1]
        ad2 = [ -d for d in rvInnerd2[n2][n1] ]
        n1 = elementsCountAroundRVFreeWall - 1
        bx  = rvInnerx [n2][n1]
        bd1 = [ -d for d in rvInnerd1[n2][n1] ]
        bd2 = rvInnerd2[n2][n1]
        px, pd1, ( pd2, ) = sampleCubicHermiteCurves([ ax ] + sx + [ bx ], [ ad2 ] + sd1 + [ bd2 ], [ [ ad1 ] + sd2 + [ bd1 ] ], elementsCountAroundRVFreeWall + 2,
            addLengthStart = 0.5*vector.magnitude(ad2)/dFactor, lengthFractionStart = 0.5,
            addLengthEnd = 0.5*vector.magnitude(bd2)/dFactor, lengthFractionEnd = 0.5, arcLengthDerivatives = True)
        n2 = elementsCountUpLVApex - 1
        for n1 in range(1, elementsCountAroundRVFreeWall + 2):
            rvInnerx [n2].append(px [n1])
            rvInnerd1[n2].append(pd1[n1])
            o1 = (elementsCountAroundLVFreeWall + n1 - 1) % elementsCountAroundLV
            d1Factor = 1.0 if (n1 == 1) else (-1.0 if (n1 == (elementsCountAroundRVFreeWall + 1)) else 0.0)
            if (n1 == 1) or (n1 == elementsCountAroundRVFreeWall + 1):
                # collapsed RV corner uses triangle derivative d/dx{1|3} = +/-d2; outside d/dxi2 = +/-d1
                # compute derivative 2 to fit with glm vector on inside row:
                od2 = [ (d1Factor*lvInnerd1[n2][o1][c] + lvInnerd2[n2][o1][c] + lvInnerd3[n2][o1][c]) for c in range(3) ]
                id2 = getHermiteLagrangeEndDerivative(lvInnerx[n2][o1], od2, px[n1])
                rvInnerd2[n2].append(id2)
            else:
                rvInnerd2[n2].append(pd2[n1])
            # compute derivative 3 to fit with glm vector on outside row:
            od3 = [ (d1Factor*vOuterd1[n2][o1][c] + vOuterd2[n2][o1][c] - vOuterd3[n2][o1][c]) for c in range(3) ]
            id3 = getHermiteLagrangeEndDerivative(vOuterx[n2][o1], od3, px[n1])
            rvInnerd3[n2].append([ -d for d in id3 ])
        for li in [ rvInnerx[n2], rvInnerd1[n2], rvInnerd2[n2], rvInnerd3[n2] ]:
            li.append(li.pop(0))

        # create nodes on inner left ventricle

        apexNodeId = [ None, None ]

        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vApexInnerx)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vApexInnerd1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vApexInnerd2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vApexInnerd3)
        apexNodeId[1] = nodeIdentifier
        nodeIdentifier = nodeIdentifier + 1

        lvInnerNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = lvInnerx[n2]
            nd1 = lvInnerd1[n2]
            nd2 = lvInnerd2[n2]
            nd3 = lvInnerd3[n2]
            layerNodeId = []
            for n1 in range(pointsCountAroundOuter):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            lvInnerNodeId.append(layerNodeId)

        # create nodes on inner right ventricle

        rvInnerNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = rvInnerx[n2]
            nd1 = rvInnerd1[n2]
            nd2 = rvInnerd2[n2]
            nd3 = rvInnerd3[n2]
            layerNodeId = []
            for n1 in range(len(nx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            rvInnerNodeId.append(layerNodeId)
        n2 = elementsCountUpLVApex - 1
        for n1 in range(elementsCountAroundRVFreeWall - 1):
            rvInnerNodeId[n2].insert(-1, rvInnerNodeId[n2][elementsCountAroundRVFreeWall - n1 - 2])

        # create nodes on outer ventricles

        apexNodeId = [ None, None ]

        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, vApexOuterx)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, vApexOuterd1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vApexOuterd2)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vApexOuterd3)
        apexNodeId[1] = nodeIdentifier
        nodeIdentifier = nodeIdentifier + 1

        vOuterNodeId = []
        for n2 in range(elementsCountUpLV):
            nx  = vOuterx[n2]
            nd1 = vOuterd1[n2]
            nd2 = vOuterd2[n2]
            nd3 = vOuterd3[n2]
            layerNodeId = []
            for n1 in range(pointsCountAroundOuter):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, nd2[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, nd3[n1])
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
            vOuterNodeId.append(layerNodeId)

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

        norl = elementsCountAroundLV
        nowl = 2 + elementsCountAroundRVFreeWall + elementsCountAroundLV*elementsCountUpLV + elementsCountAroundVSeptum*2*elementsCountUpRV
        norr = 2*elementsCountAroundRVFreeWall

        for e2 in range(elementsCountUpLV):

            if e2 < elementsCountUpLVApex:

                # LV apex

                meshGroups = [ lvMeshGroup ]
                for e1 in range(elementsCountAroundLV):

                    eft1 = eft
                    scalefactors = None

                    va = e1
                    vb = (va + 1)%elementsCountAroundLV
                    if e2 == 0:
                        # create apex elements, varying eft scale factor identifiers around apex
                        # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                        nids = [ 1       , 2 + va       , 2 + vb,
                                 1 + nowl, 2 + va + nowl, 2 + vb + nowl ]
                        eft1 = tricubichermite.createEftShellPoleBottom(va*100, vb*100)
                        # calculate general linear map coefficients
                        radiansa = apexRadiansAround[va]
                        radiansb = apexRadiansAround[vb]
                        dRadiansa = 0.5*(apexRadiansAround[va - 1] - apexRadiansAround[vb])
                        if dRadiansa < 0.0:
                            dRadiansa += math.pi
                        dRadiansb = 0.5*(apexRadiansAround[va] - apexRadiansAround[va + 2 - elementsCountAroundLV])
                        if dRadiansb < 0.0:
                            dRadiansb += math.pi
                        scalefactors = [ -1.0,
                            math.cos(radiansa), math.sin(radiansa), dRadiansa,
                            math.cos(radiansb), math.sin(radiansb), dRadiansb,
                            math.cos(radiansa), math.sin(radiansa), dRadiansa,
                            math.cos(radiansb), math.sin(radiansb), dRadiansb
                        ]
                    else:
                        bni1 = 2 + norl*(e2 - 1)
                        nids = [ bni1        + va, bni1        + vb, bni1        + norl + va, bni1        + norl + vb,
                                 bni1 + nowl + va, bni1 + nowl + vb, bni1 + nowl + norl + va, bni1 + nowl + norl + vb ]

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    #print('create element lv apex', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            if e2 >= elementsCountUpLVApex:

                # LV free wall elements starting from -1 = anterior interventricular sulcus collapsed element

                for e1 in range(-1, elementsCountAroundLVFreeWall):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ lvMeshGroup ]

                    va = e1 if (e1 >= 0) else (elementsCountAroundLV - 1)
                    vb = e1 + 1
                    nids = [ lvInnerNodeId[e2 - 1][va], lvInnerNodeId[e2 - 1][vb], lvInnerNodeId[e2][va], lvInnerNodeId[e2][vb],
                             vOuterNodeId [e2 - 1][va], vOuterNodeId [e2 - 1][vb], vOuterNodeId [e2][va], vOuterNodeId [e2][vb] ]
                    if e1 == -1:
                        # anterior interventricular sulcus: collapsed to 6 element wedge
                        nids[0] = rvInnerNodeId[e2 - 1][elementsCountAroundRVFreeWall - 1]
                        nids[2] = rvInnerNodeId[e2    ][elementsCountAroundRVFreeWall - 1]
                        nids.pop(6)
                        nids.pop(4)
                        meshGroups += [ rvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses triangle derivative d/dx1 = -d2; outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    #print('create element lv', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            if e2 >= (elementsCountUpLVApex - 1):

                # RV free wall elements starting from -1 = posterior interventricular sulcus collapsed element

                for e1 in range(-1, elementsCountAroundRVFreeWall):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ rvMeshGroup ]
                    ua = (e1 - 1) % (2*elementsCountAroundRVFreeWall)
                    ub = e1 % (2*elementsCountAroundRVFreeWall)
                    va = elementsCountAroundLVFreeWall + e1
                    vb = (va + 1)%elementsCountAroundLV
                    e2m = max(e2 - 1, elementsCountUpLVApex - 1)
                    nids = [ rvInnerNodeId[e2m][ua], rvInnerNodeId[e2m][ub], rvInnerNodeId[e2][ua], rvInnerNodeId[e2][ub],
                             vOuterNodeId [e2m][va], vOuterNodeId [e2m][vb], vOuterNodeId [e2][va], vOuterNodeId [e2][vb] ]
                    if e2 == (elementsCountUpLVApex - 1):
                        if e1 == -1:
                            continue
                        nids[0] = lvInnerNodeId[e2m][va]
                        nids[1] = lvInnerNodeId[e2m][vb]
                        nids.pop(7)
                        nids.pop(6)
                        meshGroups += [ rvMeshGroup ]
                        # collapsed elements at RV apex
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                        if e1 == 0:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        elif e1 == (elementsCountAroundVSeptum - 1):
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                            remapEftNodeValueLabel(eft1, [ 7, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                    elif e1 == -1:
                        # posterior interventricular sulcus: collapsed to 6 element wedge
                        nids[0] = lvInnerNodeId[e2 - 1][elementsCountAroundLVFreeWall]
                        nids[2] = lvInnerNodeId[e2    ][elementsCountAroundLVFreeWall]
                        nids.pop(6)
                        nids.pop(4)
                        meshGroups += [ lvMeshGroup ]
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        if e2 == elementsCountUpLVApex:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # collapsed RV corner uses triangle derivative d/dx1 = d2; outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                    elif e1 == 0:
                        # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    elif e1 == (elementsCountAroundRVFreeWall - 1):
                        # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e2 == elementsCountUpLVApex:
                            # collapsed RV corner uses outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                            remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    elif e2 == elementsCountUpLVApex:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    #print('create element rv', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

            # interventricular septum

            if e2 >= elementsCountUpLVApex:

                for e1 in range(elementsCountAroundVSeptum):

                    eft1 = eft
                    scalefactors = None
                    meshGroups = [ lvMeshGroup, rvMeshGroup, vSeptumMeshGroup ]

                    ua = 2*elementsCountAroundRVFreeWall - 1 - e1
                    ub = ua - 1
                    va = elementsCountAroundLVFreeWall + e1
                    vb = (va + 1)%elementsCountAroundLV
                    nids = [ lvInnerNodeId[e2 - 1][va], lvInnerNodeId[e2 - 1][vb], lvInnerNodeId[e2][va], lvInnerNodeId[e2][vb],
                             rvInnerNodeId[e2 - 1][ua], rvInnerNodeId[e2 - 1][ub], rvInnerNodeId[e2][ua], rvInnerNodeId[e2][ub] ]
                    if e2 == elementsCountUpLVApex:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e1 == 0:
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            # collapsed RV corner uses triangle derivative d/dx3 = d2; outside d/dxi2 = -d1
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        elif e1 == (elementsCountAroundVSeptum - 1):
                            remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            # collapsed RV corner uses triangle derivative d/dx3 = d2; outside d/dxi2 = d1
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        else:
                            remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS2, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                        scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                    else:
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [ -1.0 ]
                        if e1 == 0:
                            # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        elif e1 == (elementsCountAroundRVFreeWall - 1):
                            # general linear map d3 adjacent to collapsed anterior interventricular sulcus
                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                        scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])

                    result1 = elementtemplate1.defineField(coordinates, -1, eft1)

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    #print('create element septum', elementIdentifier, result1, result2, result3, nids)
                    elementIdentifier = elementIdentifier + 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

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
        elementsCountAroundVSeptum = elementsCountAroundRVFreeWall
        elementsCountUpLVApex = options['Number of elements up LV apex']
        elementsCountUpRV = options['Number of elements up RV']
        elementsCountUpLV = elementsCountUpLVApex + elementsCountUpRV
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        refineElementsCountThroughRVWall = options['Refine number of elements through wall']

        startPostApexElementIdentifier = elementsCountAroundLV*elementsCountUpLVApex + 1
        elementsCountPostApexLayer = elementsCountAroundLV + 2 + elementsCountAroundVSeptum
        lastVentriclesElementIdentifier = startPostApexElementIdentifier + elementsCountAroundVSeptum + elementsCountUpRV*elementsCountPostApexLayer - 1

        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = refineElementsCountSurface
            numberInXi3 = refineElementsCountThroughLVWall
            elementIdentifier = element.getIdentifier()
            if elementIdentifier >= startPostApexElementIdentifier:
                n1 = (elementIdentifier - startPostApexElementIdentifier - elementsCountAroundVSeptum)
                if n1 < 0:
                    # collapsed elements at RV apex
                    numberInXi2 = refineElementsCountThroughLVWall
                else:
                    n1 = n1 % elementsCountPostApexLayer
                    if (n1 == 0) or (n1 == (elementsCountAroundLVFreeWall + 1)):
                        # collapsed elements on posterior or anterior interventricular sulcus:
                        numberInXi1 = refineElementsCountThroughLVWall
                        #print(n1,'refine collapse element', elementIdentifier, numberInXi1, numberInXi2, numberInXi3)
                    elif (n1 > (elementsCountAroundLVFreeWall + 1)) and (n1 <= (elementsCountAroundLV + 1)):
                        numberInXi3 = refineElementsCountThroughRVWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementIdentifier == lastVentriclesElementIdentifier:
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


def getSeptumPoints(septumArcRadians, lvRadius, radialDisplacement, elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, z, n3):
    '''
    Symmetric 2-cubic interpolation of n septum elements around arc.
    :return: sx[], sd1[]
    '''
    radiansPerElementAroundLVFreeWall = (2.0*math.pi - septumArcRadians)/elementsCountAroundLVFreeWall
    # get cubic curve with arc length scaling across to centre of septum
    radiansAround = 0.5*septumArcRadians
    circleArcLength = lvRadius*radiansAround
    x1 = [ lvRadius - radialDisplacement, 0.0, z ]
    d1 = [ 0.0, circleArcLength, 0.0 ]
    x2 = [ lvRadius*math.cos(radiansAround), lvRadius*math.sin(radiansAround), z ]
    d2 = [ -circleArcLength*math.sin(radiansAround), circleArcLength*math.cos(radiansAround), 0.0 ]
    cubicLength = computeCubicHermiteArcLength(x1, d1, x2, d2, False)
    elementLengthMid = (2.0*cubicLength - lvRadius*radiansPerElementAroundLVFreeWall)/(elementsCountAroundVSeptum + n3 - 1)
    odd = elementsCountAroundVSeptum % 2
    sx, sd1, _ = sampleCubicHermiteCurves([ x1, x2 ], [ d1, d2 ], [], elementsCountAroundVSeptum//2 + odd,
        lengthFractionStart = 0.5 if odd else 1.0,
        addLengthEnd = cubicLength - 0.5*elementsCountAroundVSeptum*elementLengthMid)
    if odd:
        sx = sx[1:]
        sd1 = sd1[1:]
    return sx, sd1


def getLeftVentricleInnerPoints(lvRadius, midSeptumDisplacement, septumArcAroundRadians, z,
    elementsCountAroundLVFreeWall, elementsCountAroundVSeptum, ivSulcusDerivativeFactor = 1.0):
    '''
    Get array of points and derivatives around inside of leftventricle anticlockwise starting from
    anterior interventricular junction.
    LV is assumed to be circular, centred at x, y = (0,0), but optionally displaced inward over septum.
    :param lvRadius: Radius of the LV.
    :param midSeptumDisplacement: Inward displacement of mid septum.
    :param septumArcAroundRadians: Angle in radians around septum.
    :param z: Coordinate to give to all x[2]. All d1[2] are zero.
    :param elementsCountAroundLVFreeWall: Number of elements around LV part of outer wall.
    :param elementsCountAroundVSeptum: Number of elements around Ventricular septum.
    :param ivSulcusDerivativeFactor: Factor, typically less than 1.0 giving derivative magnitude
    at interventricular sulcus as a fraction of LV derivative, to allow a sharper bend.
    :return: Arrays nx[], nd[].
    '''
    lvArcAroundRadians = 2.0*math.pi - septumArcAroundRadians
    lvLength = lvArcAroundRadians*lvRadius
    elementSizeLV = lvLength/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
    elementSizeLVTransition = 0.5*elementSizeLV*(1.0 + ivSulcusDerivativeFactor)
    a = lvRadius - midSeptumDisplacement
    b = lvRadius
    circleLimit = 0.5*math.pi*lvRadius
    ellipseLength = 0.25*getApproximateEllipsePerimeter(a, b)
    sepLength = 2.0*(circleLimit + ellipseLength) - lvLength
    ivSulcusDerivative = ivSulcusDerivativeFactor*elementSizeLV
    elementSizeSep = (sepLength - ivSulcusDerivative)/(elementsCountAroundVSeptum - 1)
    elementSizeSepTransition = 0.5*(elementSizeSep + ivSulcusDerivative)
    # get 2-D points half way around, from mid LV free wall to mid septum
    hx = []
    hd = []
    lvOdd = (elementsCountAroundLVFreeWall % 2)
    sepOdd = (elementsCountAroundVSeptum % 2)
    nIvSulcus = elementsCountAroundLVFreeWall//2
    distance = 0.5*elementSizeLV if lvOdd else 0.0
    halfPointsCount = nIvSulcus + 1 + elementsCountAroundVSeptum//2
    for n in range(halfPointsCount):
        # get coordinate x and derivative direction d (magnitude is set later)
        if distance < circleLimit:
            angleRadians = -math.pi + distance/lvRadius
            x = [ lvRadius*math.cos(angleRadians), lvRadius*math.sin(angleRadians) ]
            d = [ -math.sin(angleRadians), math.cos(angleRadians) ]
        else:
            angleRadians = updateEllipseAngleByArcLength(a, b, -0.5*math.pi, distance - circleLimit)
            x = [ a*math.cos(angleRadians), b*math.sin(angleRadians) ]
            d = [ -a*math.sin(angleRadians), b*math.cos(angleRadians) ]
        if n < nIvSulcus:
            magnitude = elementSizeLV
            distance += elementSizeLV if (n < (nIvSulcus - 1)) else elementSizeLVTransition
        elif n == nIvSulcus:
            magnitude = ivSulcusDerivative
            distance += elementSizeSepTransition
        else:
            magnitude = elementSizeSep
            distance += elementSizeSep
        hx.append(x)
        hd.append(vector.setMagnitude(d, magnitude))
    #distance -= (0.5*elementSizeSep if sepOdd else elementSizeSep)

    nx = []
    nd = []
    lvMirrorPointsLimit = nIvSulcus + lvOdd
    for p in range(lvMirrorPointsLimit):
        n = nIvSulcus - p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    for n in range(halfPointsCount):
        x = hx[n]
        nx.append([ x[0], x[1], z ])
        d = hd[n]
        nd.append([ d[0], d[1], 0.0 ])
    sepMirrorPointsStart = -1 if sepOdd else -2
    sepMirrorPointsLimit = -(elementsCountAroundVSeptum//2 + 1)
    for p in range(sepMirrorPointsStart, sepMirrorPointsLimit, -1):
        n = halfPointsCount + p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    return nx, nd


def getVentriclesOuterPoints(lvRadius, rvAddWidthRadius, rvAddCrossRadius, rvArcAroundRadians, z,
        elementsCountAroundLVFreeWall, elementsCountAroundRVFreeWall, ivSulcusDerivativeFactor = 1.0):
    '''
    Get array of points and derivatives around outside of ventricles anticlockwise starting from
    anterior interventricular sulcus.
    LV is assumed to be circular, centred at x, y = (0,0)
    :param lvRadius: Radius of the LV.
    :param rvAddWidthRadius: Radius to add to LV to get RV width at centre.
    :param rvAddCrossRadius: Additional radius to add only laterally around septum.
    :param rvArcAroundRadians: Angle in radians around outside of RV to tangent nodes where it meets LV.
    :param z: Coordinate to give to all x[2]. All d1[2] are zero.
    :param elementsCountAroundLVFreeWall: Number of elements around LV part of outer wall.
    :param elementsCountAroundRVFreeWall: Number of elements around RV part of outer wall.
    :param ivSulcusDerivativeFactor: Factor, typically less than 1.0 giving derivative magnitude
    at interventricular sulcus as a fraction of LV derivative, to allow a sharper bend.
    :return: Arrays nx[], nd[].
    '''
    lvArcAroundRadians = 2.0*math.pi - rvArcAroundRadians
    lvLength = lvArcAroundRadians*lvRadius
    elementSizeLV = lvLength/(elementsCountAroundLVFreeWall + ivSulcusDerivativeFactor - 1.0)
    elementSizeLVTransition = 0.5*elementSizeLV*(1.0 + ivSulcusDerivativeFactor)
    a = lvRadius
    b = lvRadius + rvAddCrossRadius  # cross radius
    ellipseCentrex = lvRadius + rvAddWidthRadius - a
    joinRadians = -0.5*max(math.pi, rvArcAroundRadians)
    x1 = [ lvRadius*math.cos(joinRadians), lvRadius*math.sin(joinRadians) ]
    d1 = [ -math.sin(joinRadians), math.cos(joinRadians) ]
    x2 = [ ellipseCentrex, -b ]
    d2 = [ 1.0, 0.0 ]
    cubicLength = computeCubicHermiteArcLength(x1, d1, x2, d2, True)
    d1 = [ d1[0]*cubicLength, d1[1]*cubicLength ]
    d2 = [ d2[0]*cubicLength, d2[1]*cubicLength ]
    circleLimit = lvRadius*(math.pi + joinRadians)
    cubicLimit = circleLimit + cubicLength
    ellipseLength = 0.25*getApproximateEllipsePerimeter(a, b)
    rvLength = 2.0*(circleLimit + cubicLength + ellipseLength) - lvLength
    ivSulcusDerivative = ivSulcusDerivativeFactor*elementSizeLV
    elementSizeRV = (rvLength - ivSulcusDerivative)/(elementsCountAroundRVFreeWall - 1)
    elementSizeRVTransition = 0.5*(elementSizeRV + ivSulcusDerivative)
    # get 2-D points half way around, from mid LV free wall to mid RV free wall
    hx = []
    hd = []
    lvOdd = (elementsCountAroundLVFreeWall % 2)
    rvOdd = (elementsCountAroundRVFreeWall % 2)
    nIvSulcus = elementsCountAroundLVFreeWall//2
    distance = 0.5*elementSizeLV if lvOdd else 0.0
    halfPointsCount = nIvSulcus + 1 + elementsCountAroundRVFreeWall//2
    for n in range(halfPointsCount):
        # get coordinate x and derivative direction d (magnitude is set later)
        if distance < circleLimit:
            angleRadians = -math.pi + distance/lvRadius
            x = [ lvRadius*math.cos(angleRadians), lvRadius*math.sin(angleRadians) ]
            d = [ -math.sin(angleRadians), math.cos(angleRadians) ]
        elif distance < cubicLimit:
            xi = (distance - circleLimit)/cubicLength
            x = list(interpolateCubicHermite(x1, d1, x2, d2, xi))
            d = interpolateCubicHermiteDerivative(x1, d1, x2, d2, xi)
        else:
            angleRadians = updateEllipseAngleByArcLength(a, b, -0.5*math.pi, distance - cubicLimit)
            x = [ ellipseCentrex + a*math.cos(angleRadians), b*math.sin(angleRadians) ]
            d = [ -a*math.sin(angleRadians), b*math.cos(angleRadians) ]
        if n < nIvSulcus:
            magnitude = elementSizeLV
            distance += elementSizeLV if (n < (nIvSulcus - 1)) else elementSizeLVTransition
        elif n == nIvSulcus:
            magnitude = ivSulcusDerivative
            distance += elementSizeRVTransition
        else:
            magnitude = elementSizeRV
            distance += elementSizeRV
        hx.append(x)
        hd.append(vector.setMagnitude(d, magnitude))
    #distance -= (0.5*elementSizeRV if rvOdd else elementSizeRV)

    nx = []
    nd = []
    lvMirrorPointsLimit = nIvSulcus + lvOdd
    for p in range(lvMirrorPointsLimit):
        n = nIvSulcus - p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    for n in range(halfPointsCount):
        x = hx[n]
        nx.append([ x[0], x[1], z ])
        d = hd[n]
        nd.append([ d[0], d[1], 0.0 ])
    rvMirrorPointsStart = -1 if rvOdd else -2
    rvMirrorPointsLimit = -(elementsCountAroundRVFreeWall//2 + 1)
    for p in range(rvMirrorPointsStart, rvMirrorPointsLimit, -1):
        n = halfPointsCount + p
        x = hx[n]
        nx.append([ x[0], -x[1], z ])
        d = hd[n]
        nd.append([ -d[0], d[1], 0.0 ])
    return nx, nd


def getRVOuterSize(xiUpWidth, xiUpCross, rvWidthRadius, rvExtraCrossRadiusBase):
    if xiUpWidth < 0.0:
        #print('getRVOuterSize NEGATIVE xiUpWidth', xiUpWidth)
        return 0.0, 0.0
    if xiUpCross < 0.0:
        xiUpCross = 0.0
    xiUpFast = 1.0 - (1.0 - xiUpWidth)*(1.0 - xiUpWidth)
    xiUpSlow = xiUpCross
    rvAddWidthRadius = interpolateCubicHermite([0.0], [0.0], [rvWidthRadius], [0.0], xiUpFast)[0]
    rvAddCrossRadius = interpolateCubicHermite([0.0], [0.0], [rvExtraCrossRadiusBase], [0.0], xiUpSlow)[0]
    #print('getRVOuterSize xiUpWidth', xiUpWidth, ', addWidth', rvAddWidthRadius, ', addCross', rvAddCrossRadius)
    return rvAddWidthRadius, rvAddCrossRadius
