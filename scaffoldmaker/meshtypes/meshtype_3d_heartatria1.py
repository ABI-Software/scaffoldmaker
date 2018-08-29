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
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartatria1(object):
    '''
    3-D heart atria model, suitable for attachment to the 3-D Heart Ventricles with Base 2.
    '''

    @staticmethod
    def getName():
        return '3D Heart Atria 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around atrial free wall' : 6,
            'Number of elements around atrial septum' : 2,
            'Number of elements up atria' : 4,
            'Atria base inner major axis length' : 0.55,
            'Atria base inner minor axis length' : 0.42,
            'Atria major axis rotation degrees' : 40.0,
            'Atria outer height' : 0.4,
            'Atrial septum thickness' : 0.06,
            'Atrial free wall thickness' : 0.02,
            'Atrial base wall thickness' : 0.05,
            'Atrial base slope degrees' : 30.0,
            'Aorta outer diameter' : 0.35,
            'Atrial base front incline degrees' : 30.0,
            'Atrial base back incline degrees' : 30.0,
            'Atrial base side incline degrees' : 30.0,
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
            'Number of elements around atrial free wall',
            'Number of elements around atrial septum',
            'Number of elements up atria',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atria outer height',
            'Atrial septum thickness',
            'Atrial free wall thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Aorta outer diameter',
            'Atrial base front incline degrees',
            'Atrial base back incline degrees',
            'Atrial base side incline degrees',
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
        if options['Number of elements around atrial free wall'] < 4:
            options['Number of elements around atrial free wall'] = 4
        # need even number of elements around free wall
        if (options['Number of elements around atrial free wall'] % 2) == 1:
            options['Number of elements around atrial free wall'] += 1
        if options['Number of elements around atrial septum'] < 1:
            options['Number of elements around atrial septum'] = 1
        if options['Number of elements up atria'] < 3:
            options['Number of elements up atria'] = 3
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
        if options['Aorta outer diameter'] < options['Atrial septum thickness']:
            options['Aorta outer diameter'] = options['Atrial septum thickness']
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

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundAtrialFreeWall = options['Number of elements around atrial free wall']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundAtria = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
        elementsCountUpAtria = options['Number of elements up atria']
        aBaseInnerMajorMag = 0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = 0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aOuterHeight = options['Atria outer height']
        aortaOuterRadius = 0.5*options['Aorta outer diameter']
        aBaseFrontInclineRadians = math.radians(options['Atrial base front incline degrees'])
        aBaseSideInclineRadians = math.radians(options['Atrial base side incline degrees'])
        aBaseBackInclineRadians = math.radians(options['Atrial base back incline degrees'])
        #aortaAxis = [ 0.0, math.sin(aortaInclineRadians), math.cos(aortaInclineRadians) ]
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

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)
        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)

        laCentre, laSeptumRadians, laBaseInnerx, laBaseInnerd1, laBaseInnerd2, laBaseOuterx, laBaseOuterd1, laBaseOuterd2 = \
            getLeftAtriumBasePoints(elementsCountAroundAtrialFreeWall, elementsCountAroundAtrialSeptum,
                aBaseInnerMajorMag, aBaseInnerMinorMag, aMajorAxisRadians,
                aBaseWallThickness, aBaseSlopeHeight, aBaseSlopeLength, aSeptumThickness,
                aortaOuterRadius, aBaseFrontInclineRadians, aBaseSideInclineRadians, aBaseBackInclineRadians)

        laInnerx = [ laBaseInnerx ]
        laInnerd1 = [ laBaseInnerd1 ]
        laInnerd2 = [ laBaseInnerd2 ]
        laOuterx = [ laBaseOuterx ]
        laOuterd1 = [ laBaseOuterd1 ]
        laOuterd2 = [ laBaseOuterd2 ]
        laInnerd3 = [ [ None ]*elementsCountAroundAtria ]
        laOuterd3 = [ [ None ]*elementsCountAroundAtria ]
        for n2 in range(elementsCountUpAtria):
            for lav in (laInnerx, laInnerd1, laInnerd2, laInnerd3, laOuterx, laOuterd1, laOuterd2, laOuterd3):
                lav.append([ None ]*elementsCountAroundAtria)

        # GRC fudge factors:
        aOuterSeptumHeight = 0.7*aOuterHeight
        aSeptumBaseRowHeight = 0.2*aOuterHeight
        ridgeSeptumElementLengthFraction = 0.6

        # get la ridge points from cubic functions from ax = septum crest centre, to bx = apex, to cx = mid outer LV base
        aSeptumBaseCentreY = laCentre[2] \
            + aBaseInnerMajorMag*math.sin(-aMajorAxisRadians)*math.cos(laSeptumRadians) \
            + aBaseInnerMinorMag*math.cos(-aMajorAxisRadians)*math.sin(laSeptumRadians)
        #ax = [ 0.0, aSeptumBaseCentreY, aOuterSeptumHeight ]
        ax = [ 0.0, 0.5*(laBaseOuterx[0][1] + laBaseOuterx[elementsCountAroundAtrialFreeWall][1]), aOuterSeptumHeight ]
        ad1 = [ -1.0, 0.0, 0.0 ]
        n1MidFreeWall = elementsCountAroundAtrialFreeWall//2
        cx = laBaseOuterx[n1MidFreeWall]
        cd1 = [ -d for d in laBaseOuterd2[n1MidFreeWall]]
        # 2. apex of left atrium
        bx = [ laCentre[0], laCentre[1], aOuterHeight ]
        bd1 = [ cx[0], 0.5*(cx[1] - ax[1]), 0.0 ]
        rx, rd1, rd2 = sampleCubicHermiteCurves([ ax, bx, cx ], [ ad1, bd1, cd1 ], [ ad1, bd1, cd1 ],
            elementsCountAroundAtrialFreeWall//2 + 1,
            lengthFractionStart = 0.5*ridgeSeptumElementLengthFraction,
            lengthFractionEnd = 0.9*(0.5*elementsCountAroundAtrialFreeWall - 1.0 + 0.5*ridgeSeptumElementLengthFraction))
        rx.pop(0)
        rd1.pop(0)
        rd2.pop(0)

        # get points on outside arch of left atrium, anterior and posterior
        n2CfbCollapseTop = min(elementsCountUpAtria - 2, elementsCountUpAtria//2)
        for na in range(n1MidFreeWall):
            np = elementsCountAroundAtrialFreeWall - na
            # sample arch from double cubic through anterior, ridge and posterior points
            lx, ld2, ld1 = sampleCubicHermiteCurves(
                [ laBaseOuterx[na], rx[na], laBaseOuterx[np] ],
                [ laBaseOuterd2[na], [ -rd1[na][1], rd1[na][0], 0.0 ], [ -d for d in laBaseOuterd2[np]] ],
                [ laBaseOuterd1[na], rd1[na], [ -d for d in laBaseOuterd1[np]] ],
                2*elementsCountUpAtria)
            if na == 0:
                # move collapsed cfb nodes to x = 0, derivative 1 y = 0, derivative 2 x = 0
                for noa in range(n2CfbCollapseTop + 1):
                    ld1[noa][0] += 2.0*lx[noa][0]
                    ld1[noa][1] = 0.0
                    ld2[noa][0] = 0.0
                    lx[noa][0] = 0.0
            for noa in range(1, elementsCountUpAtria*2):
                if noa <= elementsCountUpAtria:
                    laOuterx[noa][na] = lx[noa]
                    laOuterd1[noa][na] = ld1[noa]
                    laOuterd2[noa][na] = ld2[noa]
                else:
                    nop = elementsCountUpAtria*2 - noa
                    laOuterx[nop][np] = lx[noa]
                    laOuterd1[nop][np] = [ -d for d in ld1[noa] ]
                    laOuterd2[nop][np] = [ -d for d in ld2[noa] ]
            # fix scale of base derivative 2 on anterior and posterior
            laBaseOuterd2[na] = ld2[0]
            laBaseOuterd2[np] = [-d for d in ld2[2*elementsCountUpAtria] ]

        # add end points on ridge
        ax = laBaseOuterx[n1MidFreeWall]
        ad1 = laBaseOuterd2[n1MidFreeWall]
        bx = laOuterx[elementsCountUpAtria][n1MidFreeWall - 1]
        bd1 = [ -d for d in laOuterd1[elementsCountUpAtria][n1MidFreeWall - 1] ]
        endDerivative = vector.magnitude(bd1)
        # following transitions to endDerivative
        ex, ed2, ed1 = sampleCubicHermiteCurves([ ax, bx ], [ ad1, bd1 ], [ ad1, bd1 ], elementsCountUpAtria,
            addLengthEnd = 0.5*endDerivative, lengthFractionEnd = 0.5)
        for n2 in range(1, elementsCountUpAtria):
            laOuterx[n2][n1MidFreeWall] = ex[n2]
            laOuterd1[n2][n1MidFreeWall] = getDoubleCubicHermiteCurvesMidDerivative(
                laOuterx[n2][n1MidFreeWall - 1], laOuterd1[n2][n1MidFreeWall - 1],
                laOuterx[n2][n1MidFreeWall],
                laOuterx[n2][n1MidFreeWall + 1], laOuterd1[n2][n1MidFreeWall + 1])
            laOuterd2[n2][n1MidFreeWall] = ed2[n2]
        laOuterd2[0][n1MidFreeWall] = ed2[0]

        # get inner points
        for n2 in range(1, elementsCountUpAtria + 1):
            for n1 in range(elementsCountAroundAtria):
                ox = laOuterx[n2][n1]
                if ox is None:
                    continue
                od1 = laOuterd1[n2][n1]
                od2 = laOuterd2[n2][n1]
                unitRadial = vector.normalise(vector.crossproduct3(od1, od2))
                id3 = od3 = [ aFreeWallThickness*unitRadial[c] for c in range(3) ]
                laInnerd3[n2][n1] = laOuterd3[n2][n1] = od3
                ix = laInnerx[n2][n1]
                if ix is not None:
                    continue
                ix = [ (ox[c] - od3[c]) for c in range(3) ]

                # calculate inner d1 from curvature around
                curvature = 0.0
                count = 0
                if (n1 < elementsCountAroundAtrialFreeWall) and (laOuterx[n2][n1 + 1] is not None):
                    curvature -= getCubicHermiteCurvature(laOuterx[n2][n1], laOuterx[n2][n1], laOuterx[n2][n1 + 1], laOuterx[n2][n1 + 1], unitRadial, 0.0)
                    count += 1
                if (n1 > 0) and (laOuterx[n2][n1 - 1] is not None):
                    curvature -= getCubicHermiteCurvature(laOuterx[n2][n1 - 1], laOuterx[n2][n1 - 1], laOuterx[n2][n1], laOuterx[n2][n1], unitRadial, 1.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*aFreeWallThickness
                id1 = [ factor*c for c in od1 ]

                # calculate inner d2 from curvature up
                curvature = 0.0
                count = 0
                if (n2 < elementsCountUpAtria) and (laOuterx[n2 + 1][n1] is not None):
                    curvature -= getCubicHermiteCurvature(laOuterx[n2][n1], laOuterd2[n2][n1], laOuterx[n2 + 1][n1], laOuterd2[n2 + 1][n1], unitRadial, 0.0)
                    count += 1
                if n2 > 0:
                    curvature -= getCubicHermiteCurvature(laOuterx[n2 - 1][n1], laOuterd2[n2 - 1][n1], laOuterx[n2][n1], laOuterd2[n2][n1], unitRadial, 1.0)
                    count += 1
                curvature /= count
                factor = 1.0 - curvature*aFreeWallThickness
                id2 = [ factor*c for c in od2 ]

                laInnerx[n2][n1] = ix
                laInnerd1[n2][n1] = id1
                laInnerd2[n2][n1] = id2
                laInnerd3[n2][n1] = id3

        # fix inner base derivative 2 to fit incline
        for n1 in range(1, elementsCountAroundAtrialFreeWall + 1):
            d2 = getHermiteLagrangeEndDerivative(laInnerx[1][n1], [ -d for d in laInnerd2[1][n1] ], laInnerx[0][n1])
            laInnerd2[0][n1] = [ -d for d in d2 ]
        # special fix inner base derivative 2 at cfb = slope of cfbLeft inner derivative 2
        laInnerd2[0][0] = laInnerd2[0][1]

        # fix inner nodes near septum to thicken wall and transition to septum
        ax  = laInnerx [0][0]
        ad1 = laInnerd1[0][0]
        ad2 = laInnerd2[0][0]
        # get start distance to account for aBaseSlopeRadians
        scale2 = aBaseSlopeHeight/ad2[2]
        addLengthStart = vector.magnitude([ ad2[0]*scale2, ad2[1]*scale2, aBaseSlopeHeight ])
        mx  = laInnerx [elementsCountUpAtria][0]
        md1 = laInnerd1[elementsCountUpAtria][0]
        md2 = laInnerd2[elementsCountUpAtria][0]
        md3 = laInnerd3[elementsCountUpAtria][0]
        # GRC fudge factors: 1.5x free wall thickness, 2x as steep
        mx = [ (mx[c] - 0.5*md3[c]) for c in range(3) ]
        md1[2] *= 2.0
        px  = laInnerx [0][elementsCountAroundAtrialFreeWall]
        pd1 = [ -d for d in laInnerd1[0][elementsCountAroundAtrialFreeWall] ]
        pd2 = [ -d for d in laInnerd2[0][elementsCountAroundAtrialFreeWall] ]
        # get start distance to account for aBaseSlopeRadians
        scale2 = -aBaseSlopeHeight/pd2[2]
        addLengthEnd = vector.magnitude([ pd2[0]*scale2, pd2[1]*scale2, aBaseSlopeHeight ])
        ix, id2, id1 = sampleCubicHermiteCurves(
            [ ax , mx , px  ],
            [ ad2, md2, pd2 ],
            [ ad1, md1, pd1 ],
            2*elementsCountUpAtria, addLengthStart, addLengthEnd)
        for noa in range(elementsCountUpAtria*2 + 1):
            nop = elementsCountUpAtria*2 - noa
            if noa <= elementsCountUpAtria:
                laInnerx[noa][0] = ix[noa]
                laInnerd1[noa][0] = id1[noa]
                laInnerd2[noa][0] = id2[noa]
            else:
                laInnerx[nop][elementsCountAroundAtrialFreeWall] = ix[noa]
                laInnerd1[nop][elementsCountAroundAtrialFreeWall] = [ -d for d in id1[noa] ]
                laInnerd2[nop][elementsCountAroundAtrialFreeWall] = [ -d for d in id2[noa] ]

        # fix derivative 3 through wall
        n2 = 0
        for n1 in range(elementsCountAroundAtrialFreeWall + 1, elementsCountAroundAtria):
            laInnerd3[n2][n1] = [ -2.0*laBaseInnerx[n1][0], 0.0, 0.0 ]
        for n2 in range(elementsCountUpAtria + 1):
            n1Limit = n1MidFreeWall if (n2 == elementsCountUpAtria) else (elementsCountAroundAtrialFreeWall + 1)
            for n1 in range(n1Limit):
                if laOuterx[n2][n1] is None:
                    continue
                laInnerd3[n2][n1] = laOuterd3[n2][n1] = [ (laOuterx[n2][n1][c] - laInnerx[n2][n1][c]) for c in range(3) ]
        # fix cfb/aorta centre derivative 3:
        for n2 in range(n2CfbCollapseTop + 1):
            laOuterd3[n2][0] = [ 0.0, laOuterx[n2][0][1] - laInnerx[n2][0][1], laOuterx[n2][0][2] - laInnerx[n2][0][2] ]

        laNodeId = [ [], [] ]
        raNodeId = [ [], [] ]

        # Create nodes around atria
        ran1FreeWallStart = elementsCountAroundAtrialSeptum - 1
        ran1MidFreeWall = ran1FreeWallStart + n1MidFreeWall
        for n3 in range(2):
            for n2 in range(elementsCountUpAtria + 1):
                if n3 == 0:
                    lax = laInnerx[n2]
                    lad1 = laInnerd1[n2]
                    lad2 = laInnerd2[n2]
                    lad3 = laInnerd3[n2]
                else:
                    lax =  laOuterx[n2]
                    lad1 = laOuterd1[n2]
                    lad2 = laOuterd2[n2]
                    lad3 = laOuterd3[n2]
                # left atrium
                aNodeId = [ None ]*elementsCountAroundAtria
                for n1 in range(elementsCountAroundAtria):
                    if (n2 == elementsCountUpAtria) and (n1 >= n1MidFreeWall) and (n1 <= elementsCountAroundAtrialFreeWall):
                        if n1 == n1MidFreeWall:
                            aNodeId[n1] = aNodeId[n1 - 1]
                        else:
                            aNodeId[n1] = aNodeId[elementsCountAroundAtrialFreeWall - n1]
                        continue
                    if lax[n1] is None:
                        continue
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    aNodeId[n1] = nodeIdentifier
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lax[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lad1[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lad2[n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lad3[n1])
                    nodeIdentifier += 1
                laNodeId[n3].append(aNodeId)
                # right atrium
                aNodeId = [ None ]*elementsCountAroundAtria
                for n1 in range(elementsCountAroundAtria):
                    n1l = elementsCountAroundAtria - 1 - n1
                    if (n3 == 1) and (n2 <= n2CfbCollapseTop) and (n1l == 0):
                        aNodeId[n1] = laNodeId[n3][n2][0]
                        continue
                    if (n2 == elementsCountUpAtria) and (n1 >= ran1FreeWallStart) and (n1 <= ran1MidFreeWall):
                        continue  # copy from anterior, below
                    if lax[n1l] is None:
                        continue
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    aNodeId[n1] = nodeIdentifier
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [  -lax[n1l][0],   lax[n1l][1],   lax[n1l][2] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [  lad1[n1l][0], -lad1[n1l][1], -lad1[n1l][2] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ -lad2[n1l][0],  lad2[n1l][1],  lad2[n1l][2] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ -lad3[n1l][0],  lad3[n1l][1],  lad3[n1l][2] ])
                    nodeIdentifier += 1
                if n2 == elementsCountUpAtria:
                    # fix up posterior ridge nodes
                    for n1 in range(ran1FreeWallStart, ran1MidFreeWall):
                        n1l = elementsCountAroundAtria - 1 - n1
                        aNodeId[n1] = aNodeId[n1 + 2*(ran1MidFreeWall - n1)]
                    aNodeId[ran1MidFreeWall] = aNodeId[ran1MidFreeWall + 1]
                raNodeId[n3].append(aNodeId)

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

        # create regular element rows

        for e2 in range(elementsCountUpAtria):

            # left atrium, starting at cfb
            for e1 in range(-1, elementsCountAroundAtrialFreeWall):
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [
                    laNodeId[0][e2][e1], laNodeId[0][e2][e1 + 1], laNodeId[0][e2 + 1][e1], laNodeId[0][e2 + 1][e1 + 1],
                    laNodeId[1][e2][e1], laNodeId[1][e2][e1 + 1], laNodeId[1][e2 + 1][e1], laNodeId[1][e2 + 1][e1 + 1]]
                meshGroups = [ laMeshGroup ]

                if e1 == -1:
                    # cfb straddles left and right atria
                    nids[0] = raNodeId[0][e2][-1]
                    nids[2] = raNodeId[0][e2 + 1][-1]
                    nids[4] = raNodeId[1][e2][-1]
                    nids[6] = raNodeId[1][e2 + 1][-1]
                    meshGroups += [ raMeshGroup ]
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    if e2 < n2CfbCollapseTop:
                        # collapsed to 6 element wedge
                        nids.pop(6)
                        nids.pop(4)
                        remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                        remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                        remapEftLocalNodes(eft1, 6, ln_map)
                    elif e2 == n2CfbCollapseTop:
                        # 7 node element, lower face collapsed to triangle
                        nids.pop(4)
                        remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
                        #remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                        remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        #remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                        remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 7 ]
                        remapEftLocalNodes(eft1, 7, ln_map)
                elif (e1 == 0) and (e2 <= n2CfbCollapseTop):
                    # general linear map d3 adjacent to collapsed cfb
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ] if (e2 < n2CfbCollapseTop) else [ 5 ],
                        Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                elif (e2 >= (elementsCountUpAtria - 1)) and ((e1 == (elementsCountAroundAtrialFreeWall//2 - 1)) or (e1 == (elementsCountAroundAtrialFreeWall//2))):
                    # 6 node wedge on outside of atrium
                    nids.pop(6)
                    nids.pop(2)
                    eft1 = tricubichermite.createEftShellPole90(quadrant = 3 if (e1 < elementsCountAroundAtrialFreeWall//2) else 0)
                elif (e2 == (elementsCountUpAtria - 1)) and (e1 > (elementsCountAroundAtrialFreeWall//2)):
                    # reverse d1 and d2 on ridge
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])

                if eft1 is not eft:
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                elif eft1.getNumberOfLocalScaleFactors() == 2:
                    result3 = element.setScaleFactors(eft1, [ -1.0, 0.5*math.pi ])
                else:
                    result3 = ' '
                #print('create element la', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

            # right atrium, starting at crux
            for e1 in range(-1, elementsCountAroundAtrialFreeWall):
                n1 = ran1FreeWallStart + e1
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [
                    raNodeId[0][e2][n1], raNodeId[0][e2][n1 + 1], raNodeId[0][e2 + 1][n1], raNodeId[0][e2 + 1][n1 + 1],
                    raNodeId[1][e2][n1], raNodeId[1][e2][n1 + 1], raNodeId[1][e2 + 1][n1], raNodeId[1][e2 + 1][n1 + 1]]
                meshGroups = [ laMeshGroup ]

                if e1 == -1:
                    # cfb straddles left and right atria
                    nids[0] = laNodeId[0][e2][elementsCountAroundAtrialFreeWall]
                    nids[2] = laNodeId[0][e2 + 1][elementsCountAroundAtrialFreeWall]
                    nids[4] = laNodeId[1][e2][elementsCountAroundAtrialFreeWall]
                    nids[6] = laNodeId[1][e2 + 1][elementsCountAroundAtrialFreeWall]
                    meshGroups += [ raMeshGroup ]
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                    if e2 == (elementsCountUpAtria - 1):
                        # reverse D_DS1, D_DS2 on ridge
                        scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                elif (e1 == (elementsCountAroundAtrialFreeWall - 1)) and (e2 <= n2CfbCollapseTop):
                    # general linear map d3 adjacent to collapsed cfb
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    remapEftNodeValueLabel(eft1, [ 6, 8 ] if (e2 < n2CfbCollapseTop) else [ 6 ],
                        Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                elif (e2 >= (elementsCountUpAtria - 1)) and ((e1 == (elementsCountAroundAtrialFreeWall//2 - 1)) or (e1 == (elementsCountAroundAtrialFreeWall//2))):
                    # 6 node wedge on outside of atrium
                    nids.pop(6)
                    nids.pop(2)
                    eft1 = tricubichermite.createEftShellPole90(quadrant = 1 if (e1 < elementsCountAroundAtrialFreeWall//2) else 2)
                elif (e2 == (elementsCountUpAtria - 1)) and (e1 < (elementsCountAroundAtrialFreeWall//2)):
                    # reverse d1 and d2 on ridge
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])

                if eft1 is not eft:
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1.getNumberOfLocalScaleFactors() == 1:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                elif eft1.getNumberOfLocalScaleFactors() == 2:
                    result3 = element.setScaleFactors(eft1, [ -1.0, 0.5*math.pi ])
                else:
                    result3 = ' '
                #print('create element ra', element.isValid(), elementIdentifier, result2, result3, nids)
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
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughAtrialWall = options['Refine number of elements through atrial wall']
        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            meshrefinement.refineElementCubeStandard3d(element, refineElementsCountSurface, refineElementsCountSurface, refineElementsCountThroughAtrialWall)
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


def getLeftAtriumBasePoints(elementsCountAroundAtrialFreeWall, elementsCountAroundAtrialSeptum,
        aBaseInnerMajorMag, aBaseInnerMinorMag, aMajorAxisRadians,
        aBaseWallThickness, aBaseSlopeHeight, aBaseSlopeLength, aSeptumThickness,
        aortaOuterRadius, aBaseFrontInclineRadians, aBaseSideInclineRadians, aBaseBackInclineRadians):
    """
    Get points around left atrium based on an ellipse. Points start from central fibrous body
    and wind anticlockwise around LA.
    :return: laCentre, laSeptumRadians, laBaseInnerx, laBaseInnerd1, laBaseInnerd2, laBaseOuterx, laBaseOuterd1, laBaseOuterd2.
    """

    elementsCountAroundAtrium = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
    lvOutletFrontInclineRadians = aBaseFrontInclineRadians  # for now

    aBaseOuterMajorMag = aBaseInnerMajorMag + aBaseSlopeLength
    aBaseOuterMinorMag = aBaseInnerMinorMag + aBaseSlopeLength

    # following are angles in radians around LA ellipse from major axis
    axInner = aBaseInnerMajorMag*math.cos(aMajorAxisRadians)
    bxInner = aBaseInnerMinorMag*math.sin(aMajorAxisRadians)
    laSeptumRadians = math.atan2(bxInner, axInner)
    laCentreX = -0.5*aSeptumThickness - axInner*math.cos(laSeptumRadians) - bxInner*math.sin(laSeptumRadians)
    #laCfbLeftRadians = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, laSeptumRadians, \
    #    (aSeptumBaseLength/elementsCountAroundAtrialSeptum)*(0.5*elementsCountAroundAtrialSeptum + 1.0))
    axOuter = aBaseOuterMajorMag*math.cos(aMajorAxisRadians)
    bxOuter = aBaseOuterMinorMag*math.sin(aMajorAxisRadians)
    cfbSideOffset = aortaOuterRadius*math.sin(math.pi/3.0)
    laCfbLeftRadians = getEllipseRadiansToX(axOuter, bxOuter, -cfbSideOffset - laCentreX, math.pi*0.5)
    #print('axInner',axInner,'bxInner',bxInner,'laCentreX',laCentreX)
    #print('laSeptumRadians',laSeptumRadians,'laCfbLeftRadians',laCfbLeftRadians)
    laCentreY = 0.0
    laCentreZ = 0.0

    # get points on central fibrous body centre and cfbLeft (60 degrees clockwise around aorta)
    # incline rotates about cfb-Left-Right axis, so cfb centre goes up
    cfbLeftX = -cfbSideOffset
    cfbLeftY = laCentreY + math.cos(laCfbLeftRadians)*aBaseOuterMajorMag*math.sin(-aMajorAxisRadians) \
                         + math.sin(laCfbLeftRadians)*aBaseOuterMinorMag*math.cos(-aMajorAxisRadians)
    cfbLeftZ = 0.0
    cosFrontInclineRadians = math.cos(aBaseFrontInclineRadians)
    sinFrontInclineRadians = math.sin(aBaseFrontInclineRadians)
    r = aortaOuterRadius*(1.0 - math.cos(math.pi/3.0))
    cfbX = 0.0
    cfbY = cfbLeftY - r*math.cos(lvOutletFrontInclineRadians)
    cfbZ = cfbLeftZ + r*math.sin(lvOutletFrontInclineRadians)

    pi_3 = math.pi/3.0
    lvOutletDerivativeAround = aortaOuterRadius*pi_3
    cfbLeftDerivative1 = [ \
        -lvOutletDerivativeAround*math.cos(pi_3),
        lvOutletDerivativeAround*math.sin(pi_3)*cosFrontInclineRadians,
        -lvOutletDerivativeAround*math.sin(pi_3)*sinFrontInclineRadians ]

    # compute radians around based on base outer major and minor axis sizes
    atrialPerimeterLength = getApproximateEllipsePerimeter(aBaseOuterMajorMag, aBaseOuterMinorMag)
    atrialSeptumElementLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumRadians, laCfbLeftRadians) \
        /(0.5*elementsCountAroundAtrialSeptum + 1.0)
    atrialFreeWallElementLength = (atrialPerimeterLength - atrialSeptumElementLength*(elementsCountAroundAtrialSeptum + 2)) \
        / (elementsCountAroundAtrialFreeWall - 2)
    atrialTransitionElementLength = 0.5*(atrialSeptumElementLength + atrialFreeWallElementLength)
    #atrialPerimeterLengthTmp = atrialSeptumElementLength*(elementsCountAroundAtrialSeptum + 1) + 2.0*atrialTransitionElementLength \
    #    + (elementsCountAroundAtrialFreeWall - 3)*atrialFreeWallElementLength
    #print('lengths:',(elementsCountAroundAtrialSeptum + 1),'*',atrialSeptumElementLength, '+ 2 *', \
    #    atrialTransitionElementLength,'+',(elementsCountAroundAtrialFreeWall - 3),'*',atrialFreeWallElementLength,
    #    '=', atrialPerimeterLengthTmp, ' VS ', atrialPerimeterLength)
    laRadians = []
    radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumRadians, atrialSeptumElementLength*elementsCountAroundAtrialSeptum/2.0)
    for n1 in range(elementsCountAroundAtrium):
        laRadians.append(radiansAround)
        if (n1 == 0) or (n1 >= elementsCountAroundAtrialFreeWall):
            elementLength = atrialSeptumElementLength
        elif (n1 == 1) or (n1 == (elementsCountAroundAtrialFreeWall - 1)):
            elementLength = atrialTransitionElementLength
        else:
            elementLength = atrialFreeWallElementLength
        radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)

    laBaseInnerx = []
    laBaseInnerd1 = []
    laBaseInnerd2 = []
    laBaseOuterx = []
    laBaseOuterd1 = []
    laBaseOuterd2 = []

    baseDerivative2Scale = aortaOuterRadius
    sinMajorAxisRadians = math.sin(-aMajorAxisRadians)
    cosMajorAxisRadians = math.cos(-aMajorAxisRadians)

    # get base points on inside and outside of atria
    for n3 in range(2):
        if n3 == 0:
            aMajorMag = aBaseInnerMajorMag
            aMinorMag = aBaseInnerMinorMag
            z = -aBaseSlopeHeight
        else:
            aMajorMag = aBaseOuterMajorMag
            aMinorMag = aBaseOuterMinorMag
            z = 0.0

        aMajorX =  aMajorMag*cosMajorAxisRadians
        aMajorY =  aMajorMag*sinMajorAxisRadians
        aMinorX = -aMinorMag*sinMajorAxisRadians
        aMinorY =  aMinorMag*cosMajorAxisRadians

        finalArcLength = prevArcLength = getEllipseArcLength(aMajorMag, aMinorMag, laRadians[-1] - 2.0*math.pi, laRadians[0])
        n1Limit = elementsCountAroundAtrium if (n3 == 0) else (elementsCountAroundAtrialFreeWall + 1)
        for n1 in range(n1Limit):
            radiansAround = laRadians[n1]
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)

            # get derivative around
            if n1 == (elementsCountAroundAtrium - 1):
                nextArcLength = finalArcLength
            else:
                nextArcLength = getEllipseArcLength(aMajorMag, aMinorMag, laRadians[n1], laRadians[n1 + 1])
            if (n1 <= 1) or (n1 > elementsCountAroundAtrialFreeWall):
                derivativeAround = min(prevArcLength, nextArcLength)
            else:
                derivativeAround = max(prevArcLength, nextArcLength)
            prevArcLength = nextArcLength

            if (n3 == 1) and (n1 == 0):
                x = [ cfbX, cfbY, cfbZ ]
                d1 = [ -lvOutletDerivativeAround, 0.0, 0.0 ]
                d2 = [ 0.0, baseDerivative2Scale*sinFrontInclineRadians, baseDerivative2Scale*cosFrontInclineRadians ]
            elif (n3 == 1) and (n1 == 1):
                x = [ cfbLeftX, cfbLeftY, cfbLeftZ ]
                d1 = cfbLeftDerivative1
                d2 = [ 0.0, baseDerivative2Scale*sinFrontInclineRadians, baseDerivative2Scale*cosFrontInclineRadians ]
            else:
                x = [
                    laCentreX + cosRadiansAround*aMajorX + sinRadiansAround*aMinorX,
                    laCentreY + cosRadiansAround*aMajorY + sinRadiansAround*aMinorY,
                    z ]
                d1x = -sinRadiansAround*aMajorX + cosRadiansAround*aMinorX
                d1y = -sinRadiansAround*aMajorY + cosRadiansAround*aMinorY
                scale1 = derivativeAround/math.sqrt(d1x*d1x + d1y*d1y)
                d1 = [ d1x*scale1, d1y*scale1, 0.0 ]

                if (n3 == 0) or (n1 > elementsCountAroundAtrialFreeWall):
                    d2 = [ 0.0, 0.0, baseDerivative2Scale ]  # calculated later
                elif (n1 == elementsCountAroundAtrialFreeWall):
                    # slope directly back near septum = no x component
                    d2 = [ 0.0, -baseDerivative2Scale*math.sin(aBaseBackInclineRadians), baseDerivative2Scale*math.cos(aBaseBackInclineRadians) ]
                else:
                    baseInclineRadians = 0.0
                    sideRadians = laSeptumRadians + math.pi
                    backRadians = sideRadians + math.pi*0.5
                    if radiansAround < sideRadians:
                        xi = (radiansAround - laCfbLeftRadians)/sideRadians
                        baseInclineRadians = (1.0 - xi)*aBaseFrontInclineRadians + xi*aBaseSideInclineRadians
                    elif radiansAround < backRadians:
                        xi = (radiansAround - sideRadians)/(math.pi*0.5)
                        baseInclineRadians = (1.0 - xi)*aBaseSideInclineRadians + xi*aBaseBackInclineRadians
                    else:
                        baseInclineRadians = aBaseBackInclineRadians
                    scale2 = baseDerivative2Scale/vector.magnitude(d1)
                    side = [ d1[1]*scale2, -d1[0]*scale2, 0.0 ]
                    up = [ 0.0, 0.0, baseDerivative2Scale ]
                    d2 = [ (up[c]*math.cos(baseInclineRadians) + side[c]*math.sin(baseInclineRadians)) for c in range(3) ]

            if n3 == 0:
                laBaseInnerx.append(x)
                laBaseInnerd1.append(d1)
                laBaseInnerd2.append(d2)
            else:
                laBaseOuterx.append(x)
                laBaseOuterd1.append(d1)
                laBaseOuterd2.append(d2)

    for n1 in range(elementsCountAroundAtrialFreeWall + 1, elementsCountAroundAtrium):
        laBaseOuterx.append(None)
        laBaseOuterd1.append(None)
        laBaseOuterd2.append(None)

    return [ laCentreX, laCentreY, laCentreZ ], laSeptumRadians, laBaseInnerx, laBaseInnerd1, laBaseInnerd2, laBaseOuterx, laBaseOuterd1, laBaseOuterd2
