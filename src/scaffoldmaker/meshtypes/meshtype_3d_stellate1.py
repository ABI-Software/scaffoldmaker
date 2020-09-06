"""
Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftLocalNodes
from scaffoldmaker.utils.matrix import rotateAboutZAxis
from scaffoldmaker.utils.vector import magnitude, setMagnitude
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine

class MeshType_3d_stellate1(Scaffold_base):
    """
    Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
    """
    @staticmethod
    def getName():
        return '3D Stellate 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements in all arms (e.g. 4,2,2)' : [4,2,2],
            'Element length along arm' : 1.0,
            'Element width across arm' : 0.5,
            'Element thickness' : 0.5,
            'Refine' : False,
            'Refine number of elements along arm' : 1,
            'Refine number of elements across arm' : 1,
            'Refine number of elements through thickness' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements in all arms (e.g. 4,2,2)',
            'Element length along arm',
            'Element width across arm',
            'Element thickness',
            'Refine',
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness']:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Element length along arm',
            'Element width across arm',
            'Element thickness']:
            if options[key] < 0:
                options[key] = 0.1

        key = 'Number of elements in all arms (e.g. 4,2,2)'
        if len(options[key]) < 3:
            x = [2]*(3 - len(options[key]))
            options[key].append(x)
        for i, numberOfElements in enumerate(options[key]):
            if numberOfElements < 2:
                options[key][i] = 2
            if i > 2:
                options[key].pop(i)

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        armCount = 3
        elementLengths = [options['Element length along arm'],
                 options['Element width across arm'],
                 options['Element thickness']]
        elementsCount1 = options['Number of elements in all arms (e.g. 4,2,2)']
        elementsCount2 = 2
        elementsCount3 = 1
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        # Create nodes
        numNodesPerArm = [0]
        halfArmArcAngleRadians = math.pi / armCount
        dipMultiplier = 1 #elementLengths[1] * 1.2 * 1.5

        nodeIdentifier = 1
        minArmAngle = 2 * math.pi / armCount
        xx = []
        xds1 = []
        xds2 = []
        x_in_nodes = []
        for na in range(armCount):
            elementsCount_i = [elementsCount1[na], elementsCount2, elementsCount3]
            x, ds1, ds2, nWheelEdge = createArm(halfArmArcAngleRadians, elementLengths, minArmAngle * na, minArmAngle, elementsCount_i, dipMultiplier, armCount, na)
            for ix in range(len(x)):
                if na == 0 or ix not in nWheelEdge:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[ix])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1[ix])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2[ix])
                    nodeIdentifier += 1
                    x_in_nodes.append(x[ix])
            numNodesPerArm.append(len(x))
            xx.append(x)
            xds1.append(ds1)
            xds2.append(ds2)

        # Create elements
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives() #createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        nCentre = [elementsCount1[0]+1, int(numNodesPerArm[1]/2) + elementsCount1[0]+1]
        for na in range(armCount):
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1[na]):

                        ### NODES ###
                        no2 = (elementsCount1[na] + 1)
                        no3 = (elementsCount2 + 1) * no2 - 2
                        offset = (cumNumNodesPerArm[na])
                        bni = e3 * no3 + e2 * no2 + e1 + 1 + offset
                        if e2 == 0:
                            if e1 == 0 and na > 0: # and na < armCount -1: # wheelSouth
                                nWh = cumNumNodesPerArm[na - 1] + (2 * elementsCount1[na - 1]) + 2
                                nplUq = int(numNodesPerArm[na+1]/2) - elementsCount1[na] # unused nodes at centre and shared edge
                                npl = int(numNodesPerArm[na+1]/2)  #  nodes at centre and shared edge
                                if na < armCount-1:
                                    cn = cumNumNodesPerArm[na] + elementsCount1[na]-2
                                    no2 = cumNumNodesPerArm[na]
                                    em = elementsCount1[na]
                                    nwPrev = [nWh, nWh + int(numNodesPerArm[na]/2)] # previous arm's edge, depends on armCount.
                                    nodeIdentifiers = [nwPrev[0], no2 + 1, nCentre[0], no2 + em,
                                                       nwPrev[1], no2 + em - 1 + nplUq,
                                                       nCentre[1], bni + (4 * em) - 2]
                                else:
                                    # nplPrev = int(numNodesPerArm[na]/2) - elementsCount1[na-1] - 2
                                    nplPrev = int(numNodesPerArm[na]/2) - 2
                                    no2 = elementsCount1[na]-1
                                    no3 = int(numNodesPerArm[na+1]/2) - 3
                                    nwPrev = [cumNumNodesPerArm[na-1] + 2*(elementsCount1[na-1]),
                                              cumNumNodesPerArm[na-1] + 2*(elementsCount1[na-1]) + nplPrev]
                                    start = cumNumNodesPerArm[na] - 3 # -4 + 1
                                    nodeIdentifiers = [nwPrev[0], start,
                                                       nCentre[0], start + no2,
                                                       nwPrev[1],  start + no3,
                                                       nCentre[1], start + no2 +no3]
                            elif e1 == elementsCount1[na] - 1: # armEnd, south
                                if na == 0:
                                    nodeIdentifiers = [bni, bni + no2 - 1, bni + no2, bni + no3,
                                                       bni + no2 + no3 - 1, bni + no2 + no3]
                                else:
                                    no3 = armCount*elementsCount1[na] - 1
                                    no2 = elementsCount1[na]
                                    if na > 1:
                                        bni -= 4
                                        no3 -= 1
                                    nodeIdentifiers = [bni-1, bni + no2-2, bni + no2 - 1,
                                                       bni + no3 - 1, bni + no2 -2 + no3, bni + no2 + no3 - 1]
                            elif na > 0 and e1 > 0: #  [na=1+, e1=1+, e2=0] for len=3+
                                bni -= 1 + ((armCount+1)*(na-1))
                                no2 = elementsCount1[na]
                                no3 = armCount*no2 - (na-1) - 1
                                nodeIdentifiers = [bni, bni + 1, bni + no2 - 1, bni + no2,
                                                   bni + no3, bni + no3 + 1,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                            else:
                                nodeIdentifiers = [bni, bni + 1, bni + no2 - 1, bni + no2, bni + no3, bni + no3 + 1,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                        else:
                            if e1 == 0 and na > 0:  # and na < armCount -1: # wheelNorth
                                if na < armCount - 1:
                                    bni -= armCount
                                    npl = int(numNodesPerArm[na+1] / 2) - 2
                                    no2 = elementsCount1[na]
                                    nodeIdentifiers = [nCentre[0], bni + 1, bni + no2 + 1, bni + no2 + 2,
                                                       nCentre[1], bni + npl + 1, bni + npl + no2 + 1, bni + npl + no2 + 2]
                                else: # last arm
                                    bni = cumNumNodesPerArm[na] - 2 - (armCount - elementsCount1[na])
                                    nodeIdentifiers = [nCentre[0], bni + 1, 1, bni + no2,
                                                       nCentre[1], bni + no3 - 2,
                                                       int(numNodesPerArm[1]/2)+1, bni + no2 + no3 - armCount]
                            elif e1 == elementsCount1[na] - 1: # armEnd north
                                if na > 0:
                                    no2 = elementsCount1[na]
                                    nplUq = int(numNodesPerArm[na + 1] / 2) - 2
                                    if na > 1:
                                        adj = na - 1
                                        bni -= armCount*na + (armCount-elementsCount1[na]) + 1
                                        if elementsCount1[na] < 3:
                                           bni += 1
                                        if elementsCount1[na]>3:
                                            bni -= elementsCount1[na] - 3
                                        no2 += 1 - adj
                                        no3 = nplUq - adj
                                        nodeIdentifiers = [bni, bni +1, bni+no2,
                                                           bni+no3, bni +no3+1, bni+no2 +no3]
                                    else:
                                        bni -= armCount
                                        nodeIdentifiers = [bni, bni +1, bni+no2 + 1,
                                                           bni+nplUq, bni +nplUq+1, bni+no2 +nplUq+ 1]
                                else:
                                    nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no3 - 1, bni + no3,
                                                       bni + no2 + no3 - 1]

                            elif na > 0 and e1 > 0: #  [na=1+, e1=1+, e2=1] for len=3+
                                adj = na - 1
                                bni -= armCount *na + adj
                                no2 -= adj
                                k = armCount*elementsCount1[na] - na
                                nodeIdentifiers = [bni, bni+1, bni + no2, bni + no2 + 1,
                                                   bni + k, bni + k + 1,
                                                   bni + no2 + k, bni + no2 + k + 1]
                            else:
                                nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no2,
                                                   bni + no3 - 1, bni + no3,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]

                        if e1 == 0:  # wheel
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            if armCount == 3:
                                if e2 == 0:
                                    scaleEftNodeValueLabels(eft1, [1, 5], [Node.VALUE_LABEL_D_DS1,
                                                             Node.VALUE_LABEL_D_DS2], [1])
                                    ns = [3, 7]
                                else:
                                    ns = [1, 5]
                                if na == 0:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                    elif e2 == 1:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS2, [1])])
                                elif na == 2:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                    elif e2 == 1:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                            elif armCount == 4:
                                if e2 == 0:
                                    scaleEftNodeValueLabels(eft1, [1, 5], [Node.VALUE_LABEL_D_DS1,
                                                             Node.VALUE_LABEL_D_DS2], [1])
                                    ns = [3, 7]
                                else:
                                    ns = [1, 5]
                                if na == 0:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 2:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                elif na == 3:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])

                        elif e1 < (elementsCount1[na] - 1):
                            eft1 = eft
                            elementtemplate1 = elementtemplate
                            scalefactors = None
                        else:
                            # rounded ends of arms. Collapse xi2 at xi1 = 1
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2, [])
                            # remapEftNodeValueLabel(eft1, [2, 5], Node.VALUE_LABEL_D_DS2, [])
                            if e2 == 0:
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                nodeIdentifiers = [nodeIdentifiers[0], nodeIdentifiers[2], nodeIdentifiers[1],
                                                   nodeIdentifiers[3], nodeIdentifiers[5], nodeIdentifiers[4]]
                            else:  # e2 == 1
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                            ln_map = [1, 2, 3, 2, 4, 5, 6, 5]
                            remapEftLocalNodes(eft1, 6, ln_map)

                        if eft1 is not eft:
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None
                        elementIdentifier += 1
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
        refineElementsCount1 = options['Refine number of elements along arm']
        refineElementsCount2 = options['Refine number of elements across arm']
        refineElementsCount3 = options['Refine number of elements through thickness']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount1, refineElementsCount2, refineElementsCount3)


def createArm(halfArmArcAngleRadians, elementLengths, armAngle, armAngleConst, elementsCount, dipMultiplier, armCount, armIndex):
    """
    Create single arm unit.
    Base length of element is 1.
    Direction: anticlockwise.
    Minimum arm length is 2 elements.
    :param halfArmArcAngleRadians: angle arising from base node (rad)
    :param elementLengths: [Element length along arm, half Element width across arm, Element thickness]
    :param armAngle: angle from x axis of the arm about origin (rad)
    :param armAngleConst: minimum angle from x axis of the first arm about origin (rad)
    :param elementsCount: list of numbers of elements along arm length, across width and through thickness directions
    :param dipMultiplier: factor that wheel nodes protrude by, relative to unit length
    :param armCount: number of arms in body
    :param armIndex: 0-based
    :return: x: positions of nodes in arm
    :return: nodeIdentifier: number of last node in arm +1
    :return: xnodes_ds1: ds1 derivatives of nodes associated with x
    :return: xnodes_ds2: ds2 derivatives of nodes associated with y
    :return: rmVertexNodes: indices of nodes at arm vertices around central node (including central node)
    """

    elementsCount1, elementsCount2, elementsCount3 = elementsCount
    [elementLength, elementWidth, elementHeight] = elementLengths

    xnodes_ds1 = []
    xnodes_ds2 = []
    dx_ds1 = [elementLength, 0.0, 0.0]
    dx_ds2 = [0.0, elementWidth, 0.0]
    wheelDvMult = [0.5, 0.75]

    nodes_per_layer = (elementsCount1 + 1) * (elementsCount2 + 1) - 2 # discount 2 armEnd corners
    x = []

    # nCentre and rmVertexNodes are in pythonic indexing
    nCentre = elementsCount1
    nCentre = [nCentre, nCentre + nodes_per_layer]
    if armIndex == 0:
        rmVertexNodes = []
    elif armIndex != armCount - 1:
        rmVertexNodes = nCentre + [0, nodes_per_layer]
    else:
        rmVertexNodes = nCentre + [0, nodes_per_layer,
                          2* (elementsCount1+1) - 1, 2* (elementsCount1+1) - 1 + nodes_per_layer ]
    dcent = [elementLength * math.cos(armAngleConst / 2), elementLength * math.sin(armAngleConst / 2), 0.0]
    dvertex = [elementLength * dipMultiplier * math.cos(halfArmArcAngleRadians),
                elementLength * dipMultiplier * math.sin(halfArmArcAngleRadians)]

    dipLength = 0.5*(elementLength + elementWidth)*dipMultiplier
    dvertex = [dipLength * math.cos(halfArmArcAngleRadians),
               dipLength * math.sin(halfArmArcAngleRadians)]
    dipMag = 2*dipLength - elementLength

    nid = 0
    for e3 in range(elementsCount3 + 1):
        x3 = e3 * elementHeight
        for e2 in range(elementsCount2 + 1):
            for e1 in range(elementsCount1 + 1):
                # ignore if armEnd corner nodes
                if (e1 == elementsCount1) and ((e2 == 0) or (e2 == elementsCount2)):
                    pass
                else:
                    if e1 == 0 and e2 == 0:
                        x1 = dvertex[0]
                        x2 = -dvertex[1]
                    elif e1 == 0 and e2 == elementsCount2:
                        x1 = dvertex[0]
                        x2 = dvertex[1]
                    else:
                        x1 = (elementLength*e1)
                        x2 = ((e2 - 1) * elementWidth )
                    x.append([x1, x2, x3])

                    # DERIVATIVES
                    ds1 = dx_ds1.copy()
                    ds2 = dx_ds2.copy()
                    if e2 == 1 and e1 == 0:
                        ds2 = dcent
                        ds1 = [dcent[0], -dcent[1], dcent[2]]
                    elif (e1 == 0) and ((e2 == 0) or (e2 == elementsCount2)):
                        ds2 = [dcent[0], -dcent[1], dcent[2]] if (e2 == 0) else dcent
                        ds2 = setMagnitude(ds2, dipMag)
                        ds1 = rotateAboutZAxis(ds2, -math.pi / 2)
                    elif e1 == elementsCount1 and e2 == elementsCount2-1:
                        ds2 = [0, 1 * (elementLength + elementWidth), 0]
                    xnodes_ds1.append(ds1)
                    xnodes_ds2.append(ds2)
                    nid += 1

    # d1 derivatives at dips around central node
    nidDip = [0, elementsCount1 * 2 + 1,
              nodes_per_layer, elementsCount1 * 2 + 1 + nodes_per_layer]
    for nid in nidDip:
        nx = [x[nid], x[nid+1]]
        nd1 = [xnodes_ds1[nid], xnodes_ds1[nid+1]]
        ds1Smooth = smoothCubicHermiteDerivativesLine(nx, nd1,
                                           fixStartDirection=True, fixEndDerivative=True)
        xnodes_ds1[nid] = ds1Smooth[0]


    # rotate entire arm about origin by armAngle except for centre nodes
    tol = 1e-12
    for j, n in enumerate(x):
        xynew = rotateAboutZAxis(n, armAngle)
        [xnew, ynew] = [r for r in xynew][:2]
        if (abs(xnew) < tol):
            xnew = 0
        if (abs(ynew) < tol):
            ynew = 0
        x[j] = [xnew, ynew, x[j][2]]

        xnodes_ds1[j] = rotateAboutZAxis(xnodes_ds1[j], armAngle)
        xnodes_ds2[j] = rotateAboutZAxis(xnodes_ds2[j], armAngle)

    return (x, xnodes_ds1, xnodes_ds2, rmVertexNodes)