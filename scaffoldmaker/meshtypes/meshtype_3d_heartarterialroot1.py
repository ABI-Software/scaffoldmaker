"""
Generates a 3-D heart arterial root scaffold with semilunar valve,
for attaching to a 6-element-around bicubic-linear orifice.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
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

class MeshType_3d_heartarterialroot1(object):
    '''
    Generates a 3-D heart arterial root scaffold with semilunar valve,
    for attaching to a 6-element-around bicubic-linear orifice.
    '''

    @staticmethod
    def getName():
        return '3D Heart Arterial Root 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Unit scale' : 1.0,
            'Outer height' : 0.5,
            'Inner depth' : 0.2,
            'Cusp height' : 0.6,
            'Inner diameter': 1.0,
            'Sinus radial displacement': 0.1,
            'Wall thickness': 0.1,
            'Cusp thickness' : 0.02,
            'Aortic not pulmonary' : True,
            'Refine' : False,
            'Refine number of elements surface' : 4,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Unit scale',
            'Outer height',
            'Inner depth',
            'Cusp height',
            'Inner diameter',
            'Sinus radial displacement',
            'Wall thickness',
            'Cusp thickness',
            'Aortic not pulmonary',
            'Refine',
            'Refine number of elements surface'
        ]

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        for key in [
            'Unit scale',
            'Outer height',
            'Inner depth',
            'Cusp height',
            'Inner diameter',
            'Sinus radial displacement',
            'Wall thickness',
            'Cusp thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'Refine number of elements surface']:
            if options[key] < 1:
                options[key] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options, baseCentre=[ 0.0, 0.0, 0.0 ], axisSide1=[ 0.0, -1.0, 0.0 ], axisUp=[ 0.0, 0.0, 1.0]):
        """
        Generate the base bicubic-linear Hermite mesh. See also generateMesh().
        Optional extra parameters allow centre and axes to be set.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :param baseCentre: Centre of valve on ventriculo-arterial junction.
        :param axisSide: Unit vector in first side direction where angle around starts.
        :param axisUp: Unit vector in outflow direction of valve.
        :return: list of AnnotationGroup
         """
        unitScale = options['Unit scale']
        outerHeight = unitScale*options['Outer height']
        innerDepth = unitScale*options['Inner depth']
        cuspHeight = unitScale*options['Cusp height']
        innerRadius = unitScale*0.5*options['Inner diameter']
        sinusRadialDisplacement = unitScale*options['Sinus radial displacement']
        wallThickness = unitScale*options['Wall thickness']
        cuspThickness = unitScale*options['Cusp thickness']
        aorticNotPulmonary = options['Aortic not pulmonary']
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        if aorticNotPulmonary:
            arterialRootGroup = AnnotationGroup(region, 'root of aorta', FMANumber = 3740, lyphID = 'Lyph ID unknown')
            cuspGroups = [
                AnnotationGroup(region, 'right cusp of aortic valve',     FMANumber = 7252, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'left cusp of aortic valve',      FMANumber = 7251, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'posterior cusp of aortic valve', FMANumber = 7253, lyphID = 'Lyph ID unknown') ]
        else:
            arterialRootGroup = AnnotationGroup(region, 'root of pulmonary trunk', FMANumber = 8612, lyphID = 'Lyph ID unknown')
            cuspGroups = [
                AnnotationGroup(region, 'right cusp of pulmonary valve',    FMANumber = 7250, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'anterior cusp of pulmonary valve', FMANumber = 7249, lyphID = 'Lyph ID unknown'),
                AnnotationGroup(region, 'left cusp of pulmonary valve',     FMANumber = 7247, lyphID = 'Lyph ID unknown') ]

        allGroups = [ arterialRootGroup ]  # groups that all elements in scaffold will go in
        annotationGroups = allGroups + cuspGroups

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
        # most nodes in this scaffold do not have a DS3 derivative
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        elementsCountAround = 6
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        axisSide2 = vector.crossproduct3(axisUp, axisSide1)
        outerRadius = innerRadius + wallThickness
        cuspOuterLength2 = 0.5*getApproximateEllipsePerimeter(innerRadius, cuspHeight)
        cuspOuterWallArcLength = cuspOuterLength2*innerRadius/(innerRadius + cuspHeight)
        noduleOuterAxialArcLength = cuspOuterLength2 - cuspOuterWallArcLength
        noduleOuterRadialArcLength = innerRadius
        cuspOuterWalld3 = interpolateHermiteLagrangeDerivative([ 0.0, 0.0 ], [ innerRadius, 0.0 ], [ innerRadius, outerHeight + innerDepth - cuspHeight ], 1.0)
        cuspInnerLength2 = 0.5*getApproximateEllipsePerimeter(innerRadius - 2.0*cuspThickness, cuspHeight - cuspThickness)
        noduleInnerAxialArcLength = cuspInnerLength2*(cuspHeight - cuspThickness)/(innerRadius + cuspHeight - 3.0*cuspThickness)
        noduleInnerRadialArcLength = innerRadius - cuspThickness/math.tan(math.pi/3.0)

        # lower points
        ix, id1 = createCirclePoints([ (baseCentre[c] - axisUp[c]*innerDepth) for c in range(3) ],
            [ axisSide1[c]*innerRadius for c in range(3) ], [ axisSide2[c]*innerRadius for c in range(3) ],
            elementsCountAround)
        ox, od1 = getSemilunarValveSinusPoints(baseCentre, axisSide1, axisSide2, outerRadius, sinusRadialDisplacement,
            startMidCusp=aorticNotPulmonary)
        lowerx, lowerd1 = [ ix, ox ], [ id1, od1 ]

        # upper points
        topCentre = [ (baseCentre[c] + axisUp[c]*outerHeight) for c in range(3) ]
        # twice as many on inner:
        ix, id1 = createCirclePoints(topCentre,
            [ axisSide1[c]*innerRadius for c in range(3) ], [ axisSide2[c]*innerRadius for c in range(3) ],
            elementsCountAround*2)
        ox, od1 = createCirclePoints(topCentre,
            [ axisSide1[c]*outerRadius for c in range(3) ], [ axisSide2[c]*outerRadius for c in range(3) ],
            elementsCountAround)
        upperx, upperd1 = [ ix, ox ], [ id1, od1 ]

        # get lower and upper derivative 2
        zero = [ 0.0, 0.0, 0.0 ]
        nMidCusp = 0 if aorticNotPulmonary else 1
        upperd2factor = outerHeight*(2.0/3.0)
        upd2 = [ d*upperd2factor for d in axisUp ]
        lowerOuterd2 = smoothCubicHermiteDerivativesLine([ lowerx[1][nMidCusp], upperx[1][nMidCusp] ], [ upd2, upd2 ],
            fixStartDirection=True, fixEndDerivative=True)[0]
        lowerd2 = [ [ upd2 ]*elementsCountAround, [ lowerOuterd2 ]*elementsCountAround ]  # lowerd2[0] to be fitted below
        upperd2 = [ [ upd2 ]*(elementsCountAround*2), [ upd2 ]*elementsCountAround ]

        # get lower and upper derivative 3 for points which have them
        lowerd3 = [ [ None ]*elementsCountAround, None ]
        upperd3 = [ [ None ]*(elementsCountAround*2), None ]
        for n1 in range(elementsCountAround):
            radiansAround = n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            if (n1 % 2) == nMidCusp:
                lowerd3[0][n1] = [ cuspOuterWallArcLength*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c]) for c in range(3) ]
            else:
                upperd3[0][n1*2] = [ (cuspOuterWalld3[0]*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c]) + cuspOuterWalld3[1]*axisUp[c]) for c in range(3) ]

        # inner-wall mid sinus points
        midDistance = 0.5*(outerHeight - innerDepth)
        outerArcLength2 = getCubicHermiteArcLength(lowerx[1][nMidCusp], lowerd2[1][nMidCusp], upperx[1][nMidCusp], upperd2[1][nMidCusp])
        outerArcLength2Part = outerArcLength2*midDistance/outerHeight
        mx, md2 = getCubicHermiteCurvesPointAtArcDistance([ lowerx[1][nMidCusp], upperx[1][nMidCusp] ], [ lowerd2[1][nMidCusp], upperd2[1][nMidCusp] ], outerArcLength2Part)[0:2]
        mn = vector.setMagnitude(vector.crossproduct3(lowerd1[1][nMidCusp], md2), wallThickness)
        mx = [ (mx[c] - mn[c]) for c in range(3) ]
        id2 = smoothCubicHermiteDerivativesLine([ lowerx[0][nMidCusp], mx, upperx[0][nMidCusp] ], [ lowerd2[0][nMidCusp], md2, upperd2[0][nMidCusp] ],
            fixAllDirections=True, fixEndDerivative=True)
        lowerd2[0] = [ id2[0] ]*elementsCountAround

        startRadians = 0.0 if aorticNotPulmonary else math.pi/3.0
        cosStartRadians = math.cos(startRadians)
        sinStartRadians = math.sin(startRadians)
        axisRadial = [ (cosStartRadians*axisSide1[c] + sinStartRadians*axisSide2[c]) for c in range(3) ]
        mr = vector.dotproduct([ (mx[c] - baseCentre[c]) for c in range(3) ], axisRadial)
        mc = [ (mx[c] - mr*axisRadial[c]) for c in range(3) ]
        md2a = vector.dotproduct(id2[1], axisUp)
        md2r = vector.dotproduct(id2[1], axisRadial)
        ix, id1 = getSemilunarValveSinusPoints(mc, axisSide1, axisSide2, innerRadius, mr - innerRadius,
            startMidCusp=aorticNotPulmonary)
        # resample to double number of points, halving derivative size:
        jx  = []
        jd1 = []
        xi = 0.5
        for n1 in range(elementsCountAround):
            np = (n1 + 1)%elementsCountAround
            jx .append(interpolateCubicHermite          (ix[n1], id1[n1], ix[np], id1[np], xi))
            jd1.append(interpolateCubicHermiteDerivative(ix[n1], id1[n1], ix[np], id1[np], xi))
        sinusx = []
        sinusd1 = []
        for n1 in range(elementsCountAround):
            sinusx .append(ix [n1])
            sinusx .append(jx [n1])
            sinusd1.append([ 0.5*d for d in id1[n1]])
            sinusd1.append([ 0.5*d for d in jd1[n1]])
        sinusd2 = [ None ]*(elementsCountAround*2)
        for n1 in range(elementsCountAround*2):
            radiansAround = 0.5*n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            n1Cusp = (n1 - nMidCusp*2)%4
            rf = 0.0 if (n1Cusp == 2) else (1.0 if (n1Cusp == 0) else 0.5)
            sinusd2[n1] = [ (rf*md2r*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c]) + md2a*axisUp[c]) for c in range(3) ]

        # cusp nodule points
        noduleCentre = [ (baseCentre[c] + axisUp[c]*(cuspHeight - innerDepth)) for c in range(3) ]
        nodulex  = [ [], [] ]
        noduled1 = [ [], [] ]
        noduled2 = [ [], [] ]
        noduled3 = [ [], [] ]
        for i in range(3):
            nodulex[0].append(noduleCentre)
            n1 = i*2 + nMidCusp
            radiansAround = n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            nodulex[1].append([ (noduleCentre[c] + 2.0*cuspThickness*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c])) for c in range(3) ])
            n1 = i*2 - 1 + nMidCusp
            radiansAround = n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            noduled1[0].append([ noduleOuterRadialArcLength*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c]) for c in range(3) ])
            noduled1[1].append(vector.setMagnitude(noduled1[0][i], noduleInnerRadialArcLength))
            n1 = i*2 + 1 + nMidCusp
            radiansAround = n1*radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            noduled2[0].append([ noduleOuterRadialArcLength*(cosRadiansAround*axisSide1[c] + sinRadiansAround*axisSide2[c]) for c in range(3) ])
            noduled2[1].append(vector.setMagnitude(noduled2[0][i], noduleInnerRadialArcLength))
            noduled3[0].append([ -noduleOuterAxialArcLength*axisUp[c] for c in range(3) ])
            noduled3[1].append([ -noduleInnerAxialArcLength*axisUp[c] for c in range(3) ])

        # Create nodes

        lowerNodeId = [ [], [] ]
        for n3 in range(2):
            for n1 in range(elementsCountAround):
                node = nodes.createNode(nodeIdentifier, nodetemplate if (lowerd3[n3] and lowerd3[n3][n1]) else nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lowerx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lowerd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lowerd2[n3][n1])
                if lowerd3[n3] and lowerd3[n3][n1]:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lowerd3[n3][n1])
                lowerNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

        sinusNodeId = []
        for n1 in range(elementsCountAround*2):
            if ((n1 - nMidCusp*2)%4) == 2:
                sinusNodeId.append(None)
                continue
            node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sinusx [n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sinusd1[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sinusd2[n1])
            sinusNodeId.append(nodeIdentifier)
            nodeIdentifier += 1

        noduleNodeId = [ [], [] ]
        for n3 in range(2):
            for n1 in range(3):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nodulex [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, noduled1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, noduled2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, noduled3[n3][n1])
                noduleNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

        upperNodeId = [ [], [] ]
        for n3 in range(2):
            for n1 in range(len(upperx[n3])):
                node = nodes.createNode(nodeIdentifier, nodetemplate if (upperd3[n3] and upperd3[n3][n1]) else nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, upperx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, upperd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, upperd2[n3][n1])
                if upperd3[n3] and upperd3[n3][n1]:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, upperd3[n3][n1])
                upperNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1



        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        allMeshGroups = [ allGroup.getMeshGroup(mesh) for allGroup in allGroups ]
        cuspMeshGroups = [ cuspGroup.getMeshGroup(mesh) for cuspGroup in cuspGroups ]

        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives()

        #tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # wall elements
        for cusp in range(3):
            n1 = cusp*2 - 1 + nMidCusp
            n2 = n1*2
            for e in range(6):
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]

                if (e == 0) or (e == 5):
                    # 6 node collapsed wedge element expanding from zero width on outer wall of root, attaching to vertical part of cusp
                    if e == 0:
                        nids = [ lowerNodeId[0][n1], sinusNodeId[n2 + 1], upperNodeId[0][n2], upperNodeId[0][n2 + 1], lowerNodeId[1][n1], upperNodeId[1][n1] ]
                        remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    else:
                        nids = [ sinusNodeId[n2 + 3], lowerNodeId[0][n1 - 4], upperNodeId[0][n2 + 3], upperNodeId[0][n2 - 8], lowerNodeId[1][n1 - 4], upperNodeId[1][n1 - 4] ]
                        remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                    ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                    remapEftLocalNodes(eft1, 6, ln_map)
                elif (e == 1) or (e == 4):
                    # 6 node collapsed wedge element on lower wall
                    if e == 1:
                        nids = [ lowerNodeId[0][n1], lowerNodeId[0][n1 + 1], sinusNodeId[n2 + 1], sinusNodeId[n2 + 2], lowerNodeId[1][n1], lowerNodeId[1][n1 + 1] ]
                        remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    else:
                        nids = [ lowerNodeId[0][n1 + 1], lowerNodeId[0][n1 - 4], sinusNodeId[n2 + 2], sinusNodeId[n2 + 3], lowerNodeId[1][n1 + 1], lowerNodeId[1][n1 - 4] ]
                        remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS2, [])
                    ln_map = [ 1, 2, 3, 4, 5, 6, 5, 6 ]
                    remapEftLocalNodes(eft1, 6, ln_map)
                else:
                    if e == 2:
                        nids = [ sinusNodeId[n2 + 1], sinusNodeId[n2 + 2], upperNodeId[0][n2 + 1], upperNodeId[0][n2 + 2],
                                 lowerNodeId[1][n1], lowerNodeId[1][n1 + 1], upperNodeId[1][n1], upperNodeId[1][n1 + 1] ]
                    else:
                        nids = [ sinusNodeId[n2 + 2], sinusNodeId[n2 + 3], upperNodeId[0][n2 + 2], upperNodeId[0][n2 + 3],
                                 lowerNodeId[1][n1 + 1], lowerNodeId[1][n1 - 4], upperNodeId[1][n1 + 1], upperNodeId[1][n1 - 4] ]

                result = elementtemplate1.defineField(coordinates, -1, eft1)
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = 7
                print('create arterial root wall', cusp, e, 'element',elementIdentifier, result, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in allMeshGroups:
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
        for cusp in range(3):
            for e in range(6):
                element = meshrefinement._sourceElementiterator.next()
                numberInXi1 = refineElementsCountSurface
                numberInXi2 = refineElementsCountSurface
                numberInXi3 = 1
                if (e == 0) or (e == 5):
                    numberInXi1 = 1
                elif (e == 1) or (e == 4):
                    numberInXi2 = 1
                meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)

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


def getSemilunarValveSinusPoints(centre, axisSide1, axisSide2, radius, sinusRadialDisplacement, startMidCusp, elementsCountAround = 6):
    '''
    Get points around a circle of radius with every second point displaced radially.
    :return: x[], d1[] for 6 points around
    '''
    assert elementsCountAround == 6, 'getArterialRootLowerPoints.  Only supports 6 elements around'
    px = []
    pd1 = []
    # every second outer points is displaced by sinusRadialDisplacement to make sinus dilatation
    nMidCusp = 0 if startMidCusp else 1
    radiusPlus = radius + sinusRadialDisplacement
    radiansPerElementAround = 2.0*math.pi/elementsCountAround
    for n in range(elementsCountAround):
        radiansAround = n*radiansPerElementAround
        midCusp = (n % 2) == nMidCusp
        r = radiusPlus if midCusp else radius
        rcosRadiansAround = r*math.cos(radiansAround)
        rsinRadiansAround = r*math.sin(radiansAround)
        px.append([ (centre[c] + rcosRadiansAround*axisSide1[c] + rsinRadiansAround*axisSide2[c]) for c in range(3) ])
        pd1.append([ radiansPerElementAround*(-rsinRadiansAround*axisSide1[c] + rcosRadiansAround*axisSide2[c]) for c in range(3) ])
    # smooth to get derivative in sinus
    sd1 = smoothCubicHermiteDerivativesLine(px[1 - nMidCusp:3 - nMidCusp], pd1[1 - nMidCusp:3 - nMidCusp], fixStartDerivative=True, fixEndDirection=True)
    magSinus = vector.magnitude(sd1[1])
    for n in range(nMidCusp, elementsCountAround, 2):
        pd1[n] = vector.setMagnitude(pd1[n], magSinus)

    return px, pd1
