"""
Generates a 3-D hip socket mesh with variable numbers of elements
around, up and through the thickness.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import matrix

class MeshType_3d_hip1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Hip 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 8,
            'Number of elements up' : 6,
            'Inner radius': 1.0,
            'CE factor': 0.5,
            'Wall thickness': 0.1,
            'Pincer height factor': 0.2,
            'Pincer width factor': 0.2,
            'Alpha factor': 0.5,
            'Cam width factor': 0.2,
            'Cam thickness factor': 0.1,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements up' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements up',
            'Inner radius',
            'CE factor',
            'Wall thickness',
            'Pincer height factor',
            'Pincer width factor',
            'Alpha factor',
            'Cam width factor',
            'Cam thickness factor',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 2:
            options['Number of elements around'] = 2
        if options['Number of elements around'] % 2 > 0:
            options['Number of elements around'] = options['Number of elements around'] + 1
        for key in [
            'Inner radius',
            'CE factor',
            'Wall thickness',
            'Pincer height factor',
            'Pincer width factor',
            'Alpha factor',
            'Cam width factor',
            'Cam thickness factor']:
            if options[key] < 0.0:
                options[key] = 0.0
        if 0.0 < options['Pincer width factor'] < 0.05:
            options['Pincer width factor'] = 0.05
        if 0.0 < options['Cam width factor'] < 0.05:
            options['Cam width factor'] = 0.05

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite or bicubic Hermite linear mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughWall = 1 # options['Number of elements through wall']
        useCrossDerivatives = False # options['Use cross derivatives']
        useCubicHermiteThroughWall = True # not(options['Use linear through wall'])
        innerRadius = options['Inner radius']
        CEFactor = options['CE factor']
        wallThickness = options['Wall thickness']
        pincerHeightFactor = options['Pincer height factor']
        pincerWidthFactor = options['Pincer width factor']
        alphaFactor = options['Alpha factor']
        camWidthFactor = options['Cam width factor']
        camThicknessFactor = options['Cam thickness factor']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if useCubicHermiteThroughWall:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        elementDataCloud = 20
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementDataCloud
        radiansPerElementUp = CEFactor * math.pi / elementDataCloud
        elementsDataAroundPincer = int(pincerWidthFactor * elementDataCloud)

        xList = []
        for n1 in range(elementDataCloud):
            radiansAround = radiansPerElementAround * n1
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            xAlong = []
            for n2 in range(elementDataCloud):
                radiansUp = radiansPerElementUp * n2
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                v = [ innerRadius * cosRadiansAround * sinRadiansUp,
                      innerRadius * sinRadiansAround * sinRadiansUp,
                     -innerRadius * cosRadiansUp]
                xAlong.append(v)
            xList.append(xAlong)

        # Create cam
        if camThicknessFactor > 0.0:
            radiansCamWidth = camWidthFactor * math.pi
            radiusCamHead = (1 + camThicknessFactor) * innerRadius
            radiansAroundCamList = [-radiansCamWidth, 0.0, radiansCamWidth]
            radiusHeadList = [innerRadius, radiusCamHead, innerRadius]
            elementCamStartUpIdx = int(elementDataCloud * alphaFactor)
            elementsDataAroundCam = int(camWidthFactor * elementDataCloud)
            radiansPerElementUpCam = 0.8 * math.pi / elementDataCloud

            trackRadianUp = elementCamStartUpIdx * radiansPerElementUpCam
            while trackRadianUp < CEFactor * math.pi:
                cosRadiansUp = math.cos(trackRadianUp)
                sinRadiansUp = math.sin(trackRadianUp)
                nx = []
                nd1 = []
                for i in range(3):
                    radiusHead = radiusHeadList[i]
                    radiansAround = radiansAroundCamList[i]
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    v = [ radiusHead * cosRadiansAround * sinRadiansUp,
                          radiusHead * sinRadiansAround * sinRadiansUp,
                          -radiusHead * cosRadiansUp]
                    nx.append(v)
                    d1 = [ radiusHead * -sinRadiansAround * radiansPerElementAround,
                           radiusHead * cosRadiansAround * radiansPerElementAround,
                           0.0]
                    nd1.append(d1)
                sx, sd1 = interp.sampleCubicHermiteCurves(nx, nd1, elementsDataAroundCam, arcLengthDerivatives=True)[0:2]

                count = 1
                for n in range(-int(elementsDataAroundCam * 0.5 - 1), int(elementsDataAroundCam * 0.5)):
                    n2 = int(trackRadianUp * elementDataCloud / (CEFactor * math.pi))
                    xList[n][n2] = sx[count]
                    count += 1
                trackRadianUp += radiansPerElementUp

        # Create pincer
        if pincerHeightFactor > 0.0:
            nx = []
            nd1 = []
            if camThicknessFactor:
                radiansUp = trackRadianUp
            radiansUpToPincer = radiansUp
            radiansUpPincer = radiansUp + pincerHeightFactor * CEFactor * math.pi
            radiansPincerWidth = pincerWidthFactor * math.pi
            radiansUpPincerList = [ radiansUp, radiansUpPincer, radiansUp ]
            radiansAroundPincerList = [-radiansPincerWidth, 0.0, radiansPincerWidth]
            for i in range(3):
                radiansUp = radiansUpPincerList[i]
                radiansAround = radiansAroundPincerList[i]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                radius = radiusCamHead if camThicknessFactor else innerRadius
                v = [ radius * cosRadiansAround * sinRadiansUp,
                      radius * sinRadiansAround * sinRadiansUp,
                      -radius * cosRadiansUp]
                nx.append(v)
                d1 = [ radius * -sinRadiansAround * radiansPerElementAround,
                       radius * cosRadiansAround * radiansPerElementAround,
                        0.0]
                nd1.append(d1)

            sx, sd1 = interp.sampleCubicHermiteCurves(nx, nd1, elementsDataAroundPincer, arcLengthDerivatives=True)[0:2]

            count = 1
            for n in range(-int(elementsDataAroundPincer * 0.5 - 1), int(elementsDataAroundPincer * 0.5)):
                xList[n].append(sx[count])
                count += 1

            # Add more points between hemisphere & pincer
            radiansUpToPincer += radiansPerElementUp
            while radiansUpToPincer < radiansUpPincer:
                radiansUp = radiansUpToPincer
                radiansAround = 0.0
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                v = [ radius * cosRadiansAround * sinRadiansUp,
                      radius * sinRadiansAround * sinRadiansUp,
                      -radius * cosRadiansUp]
                xList[0].insert(len(xList[0]) - 1, v)
                radiansUpToPincer += radiansPerElementUp

        # zero = [0.0, 0.0, 0.0]
        # for n1 in range(len(xList)):
        #     for n2 in range(len(xList[n1])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier = nodeIdentifier + 1

        # Sample along to spread elements out along
        dList = []
        for n1 in range(len(xList)):
            dAlong = []
            for n2 in range(len(xList[n1]) - 1):
                v1 = xList[n1][n2]
                v2 = xList[n1][n2 + 1]
                d = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
                d = [c * arcLengthAround for c in vector.normalise(d)]
                dAlong.append(d)
            dAlong.append(d)
            dList.append(dAlong)

        xSampledAlongList = []
        d2SampledAlongList = []
        for n1 in range(elementDataCloud):
            nx = xList[n1]
            nd = dList[n1]
            sx, sd = interp.sampleCubicHermiteCurves(nx, nd, elementsCountUp, arcLengthDerivatives=True)[:2]
            dSmoothed = interp.smoothCubicHermiteDerivativesLine(sx, sd)
            xSampledAlongList.append(sx)
            d2SampledAlongList.append(dSmoothed)

        xSampledAroundList = []
        d1SampledAroundList = []

        for n2 in range(1, elementsCountUp + 1):
            xAround = []
            dAround = []
            for n1 in range(int(elementDataCloud * 0.5) + 2):
                xAround.append(xSampledAlongList[n1][n2])

            for n1 in range(len(xAround) - 1):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
                d = [c * arcLengthAround for c in vector.normalise(d)]
                dAround.append(d)
            del xAround[-1]

            sx = interp.sampleCubicHermiteCurves(xAround, dAround, int(elementsCountAround * 0.5))[0]
            # Reflect to get the other half
            xOpp = sx[1: -1]
            xOppArranged = []
            for i in range(len(xOpp)):
                xOppArranged.insert(0, [xOpp[i][0], -xOpp[i][1], xOpp[i][2]])
            sx += xOppArranged

            dAround = []
            for n1 in range(elementsCountAround):
                v1 = sx[n1]
                v2 = sx[(n1 + 1) % elementsCountAround]
                d = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
                d = [c * arcLengthAround for c in vector.normalise(d)]
                dAround.append(d)
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(sx, dAround, )

            xSampledAroundList.append(sx)
            d1SampledAroundList.append(d1Smoothed)

        xSampledAroundList.insert(0, [xSampledAlongList[0][0]] * elementsCountAround)

        # Calculate d2
        d2AlongAround = []
        for n1 in range(elementsCountAround):
            xUp = []
            d2Up = []
            for n2 in range(elementsCountUp + 1):
                xUp.append(xSampledAroundList[n2][n1])

            for n2 in range(len(xUp) - 1):
                v1 = xUp[n2]
                v2 = xUp[n2 + 1]
                d = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
                d = [c * arcLengthAround for c in vector.normalise(d)]
                d2Up.append(d)
            d2Up.append(d)

            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xUp, d2Up)
            d2AlongAround.append(d2Smoothed)

        d2SampledAroundList = []
        for n2 in range(elementsCountUp + 1):
            d2Around = []
            for n1 in range(elementsCountAround):
                d2 = d2AlongAround[n1][n2]
                d2Around.append(d2)
            d2SampledAroundList.append(d2Around)

        d1AroundApex = []
        for n1 in range(elementsCountAround):
            d1AroundApex.append(matrix.rotateAboutZAxis(d2SampledAroundList[0][n1], math.pi * 0.5))
        d1SampledAroundList.insert(0, d1AroundApex)

        unitD3SampledAroundList = []
        for n2 in range(len(xSampledAroundList)):
            d3Around = []
            for n1 in range(len(xSampledAroundList[n2])):
                d3 = vector.normalise(vector.crossproduct3(d1SampledAroundList[n2][n1], d2SampledAroundList[n2][n1]))
                d3Around.append(d3)
            unitD3SampledAroundList.append(d3Around)

        # Find curvature
        d1Curvature = []
        for n2 in range(1, elementsCountUp + 1):
            curvatureAround = findD1CurvatureAround(xSampledAroundList[n2], d1SampledAroundList[n2], unitD3SampledAroundList[n2])
            d1Curvature.append(curvatureAround)
        d1Curvature.insert(0, [-1.0] * elementsCountAround)

        d2Curvature = []
        for n1 in range(elementsCountAround):
            xUp = []
            d2Up = []
            d3Up = []
            for n2 in range(elementsCountUp + 1):
                xUp.append(xSampledAroundList[n2][n1])
                d2Up.append(d2SampledAroundList[n2][n1])
                d3Up.append(unitD3SampledAroundList[n2][n1])
            curvatureUp = findCurvatureAlongLine(xUp, d2Up, d3Up)
            d2Curvature.append(curvatureUp)

        # Create node list
        xList = []
        d1List = []
        d2List = []
        d3List = []

        del xSampledAroundList[0][1:]
        d2Apex = d2SampledAroundList[0][0]
        d1Apex = matrix.rotateAboutZAxis(d2Apex, math.pi * 0.5)
        d3Apex = vector.normalise(vector.crossproduct3(d1Apex, d2Apex))
        del d2SampledAroundList[0]
        del d1SampledAroundList[0]
        del unitD3SampledAroundList[0]
        d1SampledAroundList.insert(0, [d1Apex])
        d2SampledAroundList.insert(0, [d2Apex])
        unitD3SampledAroundList.insert(0,[d3Apex])

        for n2 in range(elementsCountUp + 1):
            for n3 in range(elementsCountThroughWall + 1):
                xi3 = 1 / elementsCountThroughWall * n3
                for n1 in range(len(xSampledAroundList[n2])):
                    # Coordinates
                    norm = unitD3SampledAroundList[n2][n1]
                    xIn = xSampledAroundList[n2][n1]
                    xOut = [xIn[i] + norm[i] * wallThickness for i in range(3)]
                    dWall = [wallThickness * c for c in norm]
                    x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                    xList.append(x)

                    # d1
                    factor = 1.0 - wallThickness * (1.0 + xi3) * d1Curvature[n2][n1]
                    d1 = [factor * c for c in d1SampledAroundList[n2][n1]]
                    d1List.append(d1)

                    # d2
                    factor = 1.0 - wallThickness * (1.0 + xi3) * d2Curvature[n1][n2]
                    d2 = [factor * c for c in d2SampledAroundList[n2][n1]]
                    d2List.append(d2)

                    # d3
                    d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
                    d3List.append(d3)

        for n in range(len(xList)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n])
            nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        radiansPerElementAround = 2.0 * math.pi / elementsCountAround

        # Create apex
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1) % elementsCountAround
                eft1 = eftfactory.createEftShellPoleBottom(va * 100, vb * 100)
                elementtemplate2.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate2)
                bni1 = e3 + 1
                bni2 = elementsCountThroughWall + 1 + elementsCountAround * e3 + e1 + 1
                bni3 = elementsCountThroughWall + 1 + elementsCountAround * e3 + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni1 + 1, bni2 + elementsCountAround, bni3 + elementsCountAround]
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                # set general linear map coefficients
                radiansAround = e1 * radiansPerElementAround
                radiansAroundNext = ((e1 + 1) % elementsCountAround) * radiansPerElementAround
                scalefactors = [
                    -1.0,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                ]
                result = element.setScaleFactors(eft1, scalefactors)
                elementIdentifier = elementIdentifier + 1

        # Create regular elements
        now = elementsCountAround * (elementsCountThroughWall + 1)
        for e2 in range(1, elementsCountUp):
            for e3 in range(elementsCountThroughWall):
                for e1 in range(elementsCountAround):
                    bni11 = (e2 - 1) * now + e3 * elementsCountAround + e1 + 1 + (elementsCountThroughWall + 1)
                    bni12 = (e2 - 1) * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1)
                    bni21 = (e2 - 1) * now + (e3 + 1) * elementsCountAround + e1 + 1 + (
                                elementsCountThroughWall + 1)
                    bni22 = (e2 - 1) * now + (e3 + 1) * elementsCountAround + (
                                e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1)
                    nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now,
                                       bni22 + now]
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountThroughWall)

def findD1CurvatureAround(xAround, d1Around, normsAround):
    """
    Rearrange points around so that the group of points starts from left side of the annulus to the greater curvature
    and ends on the right side of the annulus. This put points in consecutive order for calculating curvature.
    The calculated curvatures are then re-arranged such that it starts from the greater curvature and goes to the right
    side of the annulus, followed by the left side of the annulus and closing the loop back at the greater curvature.
    :param xAround: points around a loop joining to the annulus.
    :param d1Around: derivative of points.
    :param normsAround: radial normal at each point.
    :return: curvature in the original order.
    """
    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
    normsLoop = normsAround[int(len(d1Around) * 0.5 + 1):] + normsAround[: int(len(d1Around) * 0.5 + 1)]
    curvature = findCurvatureAlongLine(xLoop, d1Loop, normsLoop)
    # Rearrange to correct order
    d1CurvatureAround = curvature[int(len(xAround) * 0.5):] + curvature[: int(len(xAround) * 0.5):]

    return d1CurvatureAround

def findCurvatureAlongLine(nx, nd, radialVectors):
    """
    Calculate curvature for points lying along a line.
    :param nx: points on line
    :param nd: derivative of points on line
    :param radialVectors: radial direction, assumed normal to curve tangent at point
    :return: curvatures along points on line
    """
    curvature = []
    for n in range(len(nx)):
        if n == 0:
            curvature.append(interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0))
        elif n == len(nx) - 1:
            curvature.append(interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0))
        else:
            curvature.append(0.5 * (
                        interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0) +
                        interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0)))

    return curvature