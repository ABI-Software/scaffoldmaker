"""
Stellate mesh, iteration 2

All locations (x) and derivatives (ds) are found separately, and fed as a list to zinc.
X elements are built for the central node, and for nodes at arm ends

"""

from __future__ import division
from math import *
import numpy as np
import matplotlib.pyplot as plt
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftLocalNodes


class MeshType_3d_stellate1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Stellate 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of arms' : 3,
            'Number of elements in long arm' : 4,
            'Number of elements in short arms' : 2,
            'Element length x' : 1.1,
            'Element width y' : 0.5,
            'Element height z' : 0.5,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements 1' : 1,
            'Refine number of elements 2' : 1,
            'Refine number of elements 3' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of arms',
            'Number of elements in long arm',
            'Number of elements in short arms',
            'Element length x',
            'Element width y',
            'Element height z',
            'Refine',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements in long arm',
            'Number of elements in short arms',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3']:
            if options[key] < 1:
                options[key] = 1

        if options['Number of arms'] < 3 or options['Number of arms'] > 4:
            options['Number of arms'] = 3

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elens = [0,0,0]
        numArm = options['Number of arms']
        elens[0] = options['Element length x']
        elens[1] = (options['Element width y'])*2
        elens[2] = options['Element height z']
        elementsCount1 = [options['Number of elements in long arm'], options['Number of elements in short arms']]
        elementsCount1 = elementsCount1 + [2]*(numArm-2)
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
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        # create nodes
        numNodesPerArm = [0]*(numArm+1)
        for na in range(numArm):
            numNodesPerArm[na+1] = (elementsCount1[na] + 1) * (elementsCount2 + 1) * (elementsCount3 + 1)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        ecount = [elementsCount1, elementsCount2, elementsCount3]
        armRotationAngles = (2 * pi) / (numArm * 2)
        dipMultiplier = elens[1] * 1.2 * 0.5
        xnodes, xnodes_ds1, xnodes_ds2, wheelNodes, centreNodes, vertexNodes = createBody(elens, numArm, armRotationAngles, ecount, dipMultiplier)

        # Remove duplicate nodes, but keep the node correspondence
        x_in = [ix[:3] for ix in xnodes]
        dupNodes_dbl = findDuplicateNodes(x_in)
        dupNodes_arr = np.array(dupNodes_dbl)
        dupNodes = list(dupNodes_arr[:,1])

        armAngle = 2*pi/numArm
        dx_ds1 = [ elens[0], 0.0, 0.0 ]
        dx_ds2 = [ 0.0, elens[1]/2, 0.0 ]
        dx_ds2_unit = [ 0.0, elens[1]/2, 0.0 ]
        wheelDvMult = [0.2, 0.3]
        nodeIdentifier = 1

        for n2 in range(len(xnodes)):
            for ia, j in enumerate(cumNumNodesPerArm[1:]):
                if ((n2+1) < j) and ((n2+1) > cumNumNodesPerArm[ia-1]):
                    iArm = ia
            if nodeIdentifier not in dupNodes:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xnodes[n2])
                if xnodes_ds1:
                    ds1 = xnodes_ds1[n2]
                    ds2 = xnodes_ds2[n2]
                else:
                    ds1 = rotateByAngle_2D(dx_ds1, xnodes[n2][-1])
                    ds2 = rotateByAngle_2D(dx_ds2, xnodes[n2][-1])
                    if nodeIdentifier in centreNodes:
                        ds1 = [dipMultiplier *cos(armAngle/2), dipMultiplier *-sin(armAngle/2), 0.0]
                        ds2 = [dipMultiplier *cos(armAngle/2), dipMultiplier *sin(armAngle/2), 0.0]
                    elif nodeIdentifier in vertexNodes:
                        [p, q] = xnodes[n2][:2]
                        ds2 = [dx_ds2_unit[0]*p, dx_ds2_unit[1]*q, dx_ds2_unit[2]]
                        ds1 = rotateByAngle_2D(ds2, -pi/2)
                        for i in range(2):
                            ds1[i] *= wheelDvMult[0]
                        ds1 = [d*elens[0] for d in ds1]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2)
            nodeIdentifier += 1

        # create elements
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        for na in range(numArm):
            no2 = (elementsCount1[na] + 1)
            no3 = (elementsCount2 + 1)*no2
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1[na]):

                        ### NODES ###
                        offset = (cumNumNodesPerArm[na])
                        bni = e3 * no3 + e2 * no2 + e1 + 1 + offset
                        nodeIdentifiers = [bni, bni + 1, bni + no2, bni + no2 + 1, bni + no3, bni + no3 + 1,
                                           bni + no2 + no3, bni + no2 + no3 + 1]
                        for ins, node in enumerate(nodeIdentifiers):
                            if node in dupNodes:
                                iReplacementNode = np.where(dupNodes_arr[:, 1] == (node))[0][0]
                                replacementNode = dupNodes_dbl[iReplacementNode][0]
                                nodeIdentifiers[ins] = replacementNode

                        ### ELEMENTS ###
                        if e1 == 0:  # wheel
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            if numArm == 3:
                                if e2 == 0:
                                    scaleEftNodeValueLabels(eft1, [1, 5],
                                                            [Node.VALUE_LABEL_D_DS1,
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
                            elif numArm == 4:
                                if e2 == 0:
                                    scaleEftNodeValueLabels(eft1, [1, 5],
                                                            [Node.VALUE_LABEL_D_DS1,
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
                            if e2 == 0:
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                nodeIdentifiers = [nodeIdentifiers[0], nodeIdentifiers[3], nodeIdentifiers[2],
                                                   nodeIdentifiers[4], nodeIdentifiers[7], nodeIdentifiers[6]]
                            else:  # e2 == 1
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                nodeIdentifiers = nodeIdentifiers[:3] + nodeIdentifiers[4:7]
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

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCount1 = options['Refine number of elements 1']
        refineElementsCount2 = options['Refine number of elements 2']
        refineElementsCount3 = options['Refine number of elements 3']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount1, refineElementsCount2, refineElementsCount3)


def extrude(x, zheight, e_1d = []):
    x = np.array(x)
    x_0 = np.column_stack([x, [0] * len(x)])
    x_1 = np.column_stack([x, [zheight] * len(x)])
    x_ex = np.vstack([x_0, x_1])

    # 1D elements
    if e_1d:
        e_1d_top = e_1d.copy() + len(x_0)
        e_1d = np.vstack([e_1d, e_1d_top])
        e_z = [[a, a + len(x_0)] for a in range(len(x_0))]
        e_1d_ex = np.vstack([e_1d, e_z])
        return x_ex, e_1d_ex
    else:
        return x_ex


def createArm(thAr, elens, armLength, armAngle, armAngleConst, ecount, n0, dipMultiplier, numArm):
    '''
    Create single arm unit.
    Base length of element is 1.
    Direction: anticlockwise
    Minimum arm length is 2 elements.

    inputs:
        thAr: angle arising from base node (rad)
        elens: [element length x, half element width y, element height z]
        armLength: length of element in arm, specific to that arm.
        armAngle: angle from x axis of the arm about origin (rad)
        armAngleConst: minimum angle from x axis of the first arm about origin (rad)
        ecount: element counts in x,y,z directions
        n0: starting node number of arm
        dipMultiplier: factor that wheel nodes protrude by, relative to unit length
        numArm: numbdr of arms in body
    outputs:
        x: positions of nodes in arm
        nodeIdentifier: number of last node in arm +1
        nWheel: nodes forming the 'wheel' about the central junction of arms
        nCentre: nodes forming central junction of arms
        nWheelVertices: nodes forming vertices of central wheel
        xnodes_ds1: ds1 derivatives of nodes associated with x
        xnodes_ds2: ds2 derivatives of nodes associated with x
    '''

    elementsCount1, elementsCount2, elementsCount3 = ecount
    [elen, ewid, zheight] = elens

    xnodes_ds1 = []
    xnodes_ds2 = []
    dx_ds1 = [elens[0], 0.0, 0.0]
    dx_ds2 = [0.0, elens[1] / 2, 0.0]
    dx_ds2_unit = [0.0, elens[1] / 2, 0.0]
    wheelDvMult = [0.5, 0.75]
    curveAdjust = [0.2, 0.1]

    nodes_per_layer = (elementsCount1 + 1) * (elementsCount2 + 1)
    x = []

    nodeIdentifier = n0
    nCentre = n0 + elementsCount1 + 1
    nCentre = [nCentre, nCentre + nodes_per_layer]
    nWheelVertices = [n0, n0 + nodes_per_layer,
                      n0 + 2 + (2*elementsCount1), n0 + 2 + (2*elementsCount1) + nodes_per_layer ]
    nArmCornerSth = [n0 + elementsCount1, n0 + elementsCount1 + nodes_per_layer]
    nArmCornerNth = [n0 + elementsCount1 + ((elementsCount1+1)*2),
                     n0 + elementsCount1 + ((elementsCount1+1)*2) + nodes_per_layer]
    nArmEndMid = [(n0-1) + ((elementsCount1+1)*2), (n0-1) + ((elementsCount1+1)*2) + nodes_per_layer]
    nWheel = [n0, \
              n0 + 1, \
              n0 + elementsCount1 + 2, \
              n0 + (2*(elementsCount1 + 1)) + 1, \
              n0 + (2*(elementsCount1 + 1))]
    nWheel = nWheel + [n+nodes_per_layer for n in nWheel]

    dcent = [elen * cos(armAngleConst / 2), elen * sin(armAngleConst / 2), 0.0]
    dvertex = [round(elen * dipMultiplier * cos(thAr), 12),
                round(elen * dipMultiplier * sin(thAr), 12)]

    for n3 in range(elementsCount3 + 1):
        x3 = n3 * zheight
        for n2 in range(elementsCount2 + 1):
            for n1 in range(elementsCount1 + 1):
                if nodeIdentifier == n0 or nodeIdentifier == n0 + nodes_per_layer:
                    x1 = dvertex[0]
                    x2 = -dvertex[1]
                elif nodeIdentifier == n0 + (2*(armLength+1)) or nodeIdentifier == n0 + (2*(armLength+1)) + nodes_per_layer:
                    x1 = dvertex[0]
                    x2 = dvertex[1]
                else:
                    x1 = (elen*n1)
                    x2 = ((n2 - 1) * ewid / 2)
                x.append([x1, x2, x3])

                # DERIVATIVES
                ds1 = dx_ds1.copy()
                ds2 = dx_ds2.copy()
                if nodeIdentifier in nCentre:
                    ds2 = dcent
                    ds1 = [dcent[0], -dcent[1], dcent[2]]
                elif nodeIdentifier in nWheelVertices:
                    [p, q] = [x1/dipMultiplier, x2/dipMultiplier] # unit vector
                    ds2mag = round(2*np.linalg.norm([p,q]),6) - np.linalg.norm(dcent)
                    ds2 = [ds2mag*p, ds2mag*q]
                    ds1 = rotateByAngle_2D(ds2, -pi / 2)  # -pi/2 # -armAngle
                    for i in range(2):
                        ds1[i] *= wheelDvMult[0]
                    ds1 = [d * elen for d in ds1]
                elif nodeIdentifier in (nArmEndMid+nArmCornerNth+nArmCornerSth):
                    ds1 = [0.5 * ewid, 0, 0]
                    ds2 = [0, 0.5 * (elen + 0.5 * ewid), 0]
                xnodes_ds1.append(ds1)
                xnodes_ds2.append(ds2)

                nodeIdentifier += 1

    # rotate entire arm about origin by armAngle except for centre nodes
    tol = 1e-12
    for j, n in enumerate(x):
        xynew = rotateByAngle_2D(n, armAngle)
        [xnew, ynew] = [round(r,12) for r in xynew][:2]
        if (abs(xnew) < tol):
            xnew = 0
        if (abs(ynew) < tol):
            ynew = 0
        x[j] = [xnew, ynew, x[j][2]]

        if (j+1) not in nCentre:
            xnodes_ds1[j] = rotateByAngle_2D(xnodes_ds1[j], armAngle)
            xnodes_ds2[j] = rotateByAngle_2D(xnodes_ds2[j], armAngle)

    return (x, nodeIdentifier, nWheel, nCentre, nWheelVertices, xnodes_ds1, xnodes_ds2)


def createBody(elens, numArm, thAr, ecount, dipMultiplier):
    '''
    Create amalgamation of arms about central junction.

    inputs:
        elens: [element length x, half element width y, element height z]
        numArm: numbdr of arms in body
        thAr: angle arising from base node (rad)
        ecount: element counts in x,y,z directions
        dipMultiplier: factor that wheel nodes protrude by, relative to unit length
    outputs:
        x: positions of all nodes in body
        xds1: ds1 derivatives of nodes associated with x
        xds2: ds2 derivatives of nodes associated with x
        nWheel: nodes forming the 'wheel' about the central junction of arms
        nCentre: nodes forming central junction of arms
        nWheelVertices: nodes forming vertices of central wheel
    '''

    [elen, ewid, zheight] = elens

    armLength = ecount[0]
    minArmAngle = 2*pi/numArm

    x = []
    xds1 = []
    xds2 = []
    nWheel = []
    nCentre = []
    nWheelVertices = []

    nextNode = 1
    x_out = []
    for i in range(numArm):
        ecount_i = [ecount[0][i], ecount[1], ecount[2]]
        x_out, nextNode, nwhl, ncntr, nvtx, ds1, ds2 = createArm(thAr, elens, armLength[i], minArmAngle*i, minArmAngle, ecount_i, nextNode, dipMultiplier, numArm)

        x.extend([ix+[minArmAngle*i] for ix in x_out])
        xds1.extend(ds1)
        xds2.extend(ds2)
        nWheel.extend(nwhl)
        nCentre.extend(ncntr)
        nWheelVertices.extend(nvtx)

        endNode = 0
        if isinstance(x_out, list):
            x_out_nd = np.array(x_out)

    # extrude in z if only a single 2D layer was made
    if len(x[0]) < 4:
        x = np.array(x)
        x = extrude(x, zheight)

    return x, xds1, xds2, nWheel, nCentre, nWheelVertices


def findDuplicateNodes(x):
    '''
    Find repeated rows in x

    input:
        x: list
    output:
        dupNodes: [first instance of row, index of where duplicate occured]
    '''
    xnew = []
    x_orig = x.copy()
    dupNodes = []
    for ir, row in enumerate(x_orig):
        iOG = x_orig.index(row)
        if row in x and ir != iOG:
            dupNodes.append([iOG+1, ir+1])
        else:
            xnew.extend(row)

    return dupNodes


def rotateByAngle_2D(x, th):
    '''
    Stiff rotation in xy plane about (0,0)

    inputs:
        x: list
        th: angle to rotate anticlockwise by (rad)
    output:
        xRot: transformed list
    '''

    xRot = [x[0] * cos(th) - x[1] * sin(th),
            x[1] * cos(th) + x[0] * sin(th),
            0]
    return xRot
