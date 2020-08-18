"""
Stellate mesh, iteration 2

All locations (x) and derivatives (ds) are found separately, and fed as a list to zinc.
"""

from __future__ import division
from math import *
import numpy as np
# from numpy import *
import matplotlib.pyplot as plt
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.meshrefinement import MeshRefinement

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
            'Number of arms (1 <= 4)' : 3,
            'Number of elements in arm1' : 4,
            'Number of elements in arm2' : 2,
            'Number of elements in arm3' : 2,
            'Number of elements in arm4' : 2,
            'Element length x' : 1,
            'Element length y' : 1,
            'Element length z' : 1,
            'Central dip multiplier' : 0.75,
            'Arm end position multiplier' : 0.75,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements 1' : 1,
            'Refine number of elements 2' : 1,
            'Refine number of elements 3' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of arms (1 <= 4)',
            'Number of elements in arm1',
            'Number of elements in arm2',
            'Number of elements in arm3',
            'Number of elements in arm4',
            'Element length x',
            'Element length y',
            'Element length z',
            'Central dip multiplier',
            'Arm end position multiplier',
            'Refine',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of arms (1 <= 4)',
            'Number of elements in arm1',
            'Number of elements in arm2',
            'Number of elements in arm3',
            'Number of elements in arm4',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3']:
            if options[key] < 1:
                options[key] = 1

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elens = [0,0,0]
        numArm = options['Number of arms (1 <= 4)']
        elens[0] = options['Element length x']
        elens[1] = options['Element length y']
        elens[2] = options['Element length z']
        elementsCount1 = [options['Number of elements in arm1'], options['Number of elements in arm2'], options['Number of elements in arm3'], options['Number of elements in arm4']]
        elementsCount2 = 2 #options['Element count 2']
        elementsCount3 = 1 #options['Element count 3']
        dipMultiplier = options['Central dip multiplier']
        armEndMult = options['Arm end position multiplier']
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

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        numNodesPerArm = [0]*(numArm+1)
        # must be cumulative
        for na in range(numArm):
            numNodesPerArm[na+1] = (elementsCount1[na] + 1) * (elementsCount2 + 1) * (elementsCount3 + 1)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        plot_graph = 0
        ecount = [elementsCount1, elementsCount2, elementsCount3]
        armRotationAngles = (2 * pi) / (numArm * 2)
        xnodes, xnodes_ds1, xnodes_ds2, wheelNodes, centreNodes, vertexNodes = createBody(elens, numArm, armRotationAngles, ecount, dipMultiplier, armEndMult, plot_graph)
        nodeList_sh = range(1, len(xnodes) + 1)

        # Remove duplicate nodes, but keep the node correspondence
        if True:
            x_in = [ix[:3] for ix in xnodes]
            dupNodes_dbl = findDuplicateNodes(x_in, nodeList_sh)
            dupNodes_arr = np.array(dupNodes_dbl)
            dupNodes = list(dupNodes_arr[:,1])

            # repeat for centre nodes - but no need really.
            xc_in = [xnodes[ic-1][:3] for ic in centreNodes]
            dupCntNodes_dbl = findDuplicateNodes(xc_in, centreNodes)
            dupCntNodes_arr = np.array(dupCntNodes_dbl)
            idupCntNodes = list(np.array(dupCntNodes_dbl)[:,1])
            dupCntNodes = [centreNodes[id-1] for id in idupCntNodes]

        armAngle = 2*pi/numArm
        dx_ds1 = [ elens[0], 0.0, 0.0 ]
        dx_ds2 = [ 0.0, elens[1]/2, 0.0 ]
        dx_ds2_unit = [ 0.0, elens[1]/2, 0.0 ]
        wheelDvMult = [0.2, 0.3]
        nodeIdentifier = 1

        if True:
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
                    elif True:
                        ds1 = rotateByAngle_2D(dx_ds1, xnodes[n2][-1])
                        ds2 = rotateByAngle_2D(dx_ds2, xnodes[n2][-1])
                        if nodeIdentifier in centreNodes:
                            ds1 = [dipMultiplier *cos(armAngle/2), dipMultiplier *-sin(armAngle/2), 0.0]
                            ds2 = [dipMultiplier *cos(armAngle/2), dipMultiplier *sin(armAngle/2), 0.0]
                        elif nodeIdentifier in vertexNodes:
                            if False:
                                thNormal = -armAngle + xnodes[n2][-1] - pi
                                ds1 = rotateByAngle_2D(dx_ds1, thNormal)
                                ds2 = rotateByAngle_2D(dx_ds2, thNormal)

                            [p, q] = xnodes[n2][:2]
                            ds2 = [dx_ds2_unit[1]*p, dx_ds2_unit[1]*q, dx_ds2_unit[2]]
                            ds1 = rotateByAngle_2D(ds2, -pi/2) # -pi/2 # -armAngle
                            for i in range(2):
                                ds1[i] *= wheelDvMult[0]
                                ds2[i] *= wheelDvMult[1]
                            ds1 = [d*dx_ds1[0] for d in ds1]
                    else: # default
                        ds1 = rotateByAngle_2D(dx_ds1, xnodes[n2][-1])
                        ds2 = rotateByAngle_2D(dx_ds2, xnodes[n2][-1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1) #ds1
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2) #ds2
                nodeIdentifier += 1

        # create elements
        elementIdentifier = 1
        for na in range(numArm): # for every arm present ##############################
            no2 = (elementsCount1[na] + 1)
            no3 = (elementsCount2 + 1)*no2
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1[na]):
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        offset = (cumNumNodesPerArm[na])
                        bni = e3*no3 + e2*no2 + e1 + 1 + offset
                        nodeIdentifiers = [ bni, bni + 1, bni + no2, bni + no2 + 1, bni + no3, bni + no3 + 1, bni + no2 + no3, bni + no2 + no3 + 1 ]
                        for ins, node in enumerate(nodeIdentifiers):
                            if node in dupNodes:
                                iReplacementNode = np.where(dupNodes_arr[:,1] == (node))[0][0]
                                replacementNode = dupNodes_dbl[iReplacementNode][0]
                                nodeIdentifiers[ins] = replacementNode
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1

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
    x_0 = column_stack([x, [0] * len(x)])
    x_1 = column_stack([x, [zheight] * len(x)])
    x_ex = vstack([x_0, x_1])

    # 1D elements
    if e_1d:
        e_1d_top = e_1d.copy() + len(x_0)
        e_1d = np.vstack([e_1d, e_1d_top])
        e_z = [[a, a + len(x_0)] for a in range(len(x_0))]
        e_1d_ex = np.vstack([e_1d, e_z])
        return x_ex, e_1d_ex
    else:
        return x_ex


def createArm(thAr, elens, armLength, armAngle, armAngleConst, ecount, startingNode, dipMultiplier, armEndMult):
    '''
    Base length of element is 1
    direction: anticlockwise
    minimum arm length is 2 elements

    inputs:
        thAr: angle arising from base node (rad)
        elen: element length
        ewid: half element width
        armLength: number of elements in arm
        armAngle: angle from x axis of the arm about origin (rad)
    outputs:
        x: positions of nodes in arm
        e_1d: list of pairs of node numbers forming 1D line elements
    '''

    def round_to_6sig(x, sig=6):
        if x == 0:
            return 0
        else:
            rounded = round(x, sig-int(floor(log10(abs(x))))-1)
            return rounded

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

    nodeIdentifier = startingNode
    nCentre = startingNode + elementsCount1 + 1
    nCentre = [nCentre, nCentre + nodes_per_layer]
    nWheelVertices = [startingNode, startingNode + nodes_per_layer,
                      startingNode + 2 + (2*elementsCount1), startingNode + 2 + (2*elementsCount1) + nodes_per_layer ]
    nArmCornerSth = [startingNode + elementsCount1, startingNode + elementsCount1 + nodes_per_layer]
    nArmCornerNth = [startingNode + elementsCount1 + ((elementsCount1+1)*2), startingNode + elementsCount1 + ((elementsCount1+1)*2) + nodes_per_layer]
    nWheel = [startingNode, \
              startingNode + 1, \
              startingNode + elementsCount1 + 2, \
              startingNode + (2*(elementsCount1 + 1)) + 1, \
              startingNode + (2*(elementsCount1 + 1))]
    nWheel = nWheel + [n+nodes_per_layer for n in nWheel]

    for n3 in range(elementsCount3 + 1):
        x3 = n3 * zheight
        for n2 in range(elementsCount2 + 1):
            for n1 in range(elementsCount1 + 1):
                if nodeIdentifier == startingNode or nodeIdentifier == startingNode + nodes_per_layer:
                    x1 = round(elen * dipMultiplier * cos(thAr), 12)
                    x2 = -round(elen * dipMultiplier * sin(thAr), 12)
                elif nodeIdentifier == startingNode + (2*(armLength+1)) or nodeIdentifier == startingNode + (2*(armLength+1)) + nodes_per_layer:
                    x1 = round(elen * dipMultiplier * cos(thAr), 12)
                    x2 = round(elen * dipMultiplier * sin(thAr), 12)
                elif nodeIdentifier == startingNode + elementsCount1 or nodeIdentifier == startingNode + elementsCount1 + nodes_per_layer or nodeIdentifier == startingNode + (nodes_per_layer) - 1 or nodeIdentifier == startingNode + (2*nodes_per_layer) - 1:
                    x1 = elen*(n1 - (1-armEndMult))
                    x2 = ((n2 - 1) * armEndMult * ewid / 2)
                else:
                    x1 = (elen*n1)
                    x2 = ((n2 - 1) * ewid / 2)
                x.append([x1, x2, x3])

                # DERIVATIVES
                # ds1 = rotateByAngle_2D(dx_ds1, armAngle)  # probably the wrong angle argument
                # ds2 = rotateByAngle_2D(dx_ds2, armAngle)
                ds1 = dx_ds1.copy()
                ds2 = dx_ds2.copy()
                if nodeIdentifier in nCentre:
                    # these must consider elen
                    ds1 = [elen*dipMultiplier * cos(armAngleConst / 2), elen*dipMultiplier * -sin(armAngleConst / 2), 0.0]
                    ds2 = [elen*dipMultiplier * cos(armAngleConst / 2), elen*dipMultiplier * sin(armAngleConst / 2), 0.0]
                elif nodeIdentifier in nWheelVertices:
                    [p, q] = [x1, x2]
                    ds2 = [dx_ds2_unit[1] * p, dx_ds2_unit[1] * q, dx_ds2_unit[2]]
                    ds1 = rotateByAngle_2D(ds2, -pi / 2)  # -pi/2 # -armAngle
                    for i in range(2):
                        ds1[i] *= wheelDvMult[0]
                        ds2[i] *= wheelDvMult[1]
                    ds1 = [d * dx_ds1[0] for d in ds1]
                elif nodeIdentifier in nArmCornerSth:
                    ds1 = [ds1[0]*curveAdjust[0], curveAdjust[0], 0]#[d for d in ds1]
                    ds2 = [curveAdjust[1], ds2[1]*curveAdjust[1], 0]
                elif nodeIdentifier in nArmCornerNth:
                    ds1 = [ds1[0]*curveAdjust[0], -curveAdjust[0], 0]#[d for d in ds1]
                    ds2 = [-curveAdjust[1], ds2[1]*curveAdjust[1], 0]
                xnodes_ds1.append(ds1)
                xnodes_ds2.append(ds2)

                nodeIdentifier += 1
    e_1d = []


    # rotate entire arm about origin by armAngle EXCEPT for centre nodes
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


def createBody(elens, numArm, thAr, ecount, dipMultiplier, armEndMult, plot_):

    import numpy as np

    [elen, ewid, zheight] = elens

    armLength = ecount[0]
    minArmAngle = 2*pi/numArm

    x = []
    xds1 = []
    xds2 = []
    nWheel = []
    nCentre = []
    nWheelVertices = []

    if plot_:
        plt.figure()

    nextNode = 1
    x_out = []
    for i in range(numArm):
        ecount_i = [ecount[0][i], ecount[1], ecount[2]]
        x_out, nextNode, nwhl, ncntr, nvtx, ds1, ds2 = createArm(thAr, elens, armLength[i], minArmAngle*i, minArmAngle, ecount_i, nextNode, dipMultiplier, armEndMult)

        x.extend([ix+[minArmAngle*i] for ix in x_out])
        xds1.extend(ds1)
        xds2.extend(ds2)
        nWheel.extend(nwhl)
        nCentre.extend(ncntr)
        nWheelVertices.extend(nvtx)

        endNode = 0
        if isinstance(x_out, list):
            x_out_nd = np.array(x_out)
            # x = np.array(x_out)
        if plot_:
            plt.plot(x_out_nd[:, 0], x_out_nd[:, 1], 'x')
            for node in range(int(len(x_out_nd)/1)):
                plt.text(x_out_nd[node][0], x_out_nd[node][1], str(node+1+endNode))
            endNode += len(x_out_nd)
            if False:
                e_1d_out = [] # obsolete
                for np in e_1d_out:
                    plt.plot(x_out[np][:2,0], x_out[np][:2,1], 'k-', linewidth = 1)
            plt.axis('equal')

    # extrude in z
    if len(x[0]) < 4:
        x = np.array(x)
        x = extrude(x, zheight)  # x, e_1d = extrude(x, e_1d)

    # adjust value of x by small random offset
    if False:
        x += random.rand(len(x), len(x[0])) * 1e-2
        x = list([list(ix) for ix in x])

    # TESTING ELEMENT CREATION
    if False:
        numNodesPerArm = [0]*(numArm+1)
        for na in range(numArm):
            numNodesPerArm[na+1] = (ecount[0][na] + 1) * (ecount[1] + 1) * (ecount[2] + 1)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        elementIdentifier = 1
        scr = dict()
        for na in range(numArm): # for every arm present ##############################
            no2 = (ecount[0][na] + 1) #+ (numArm-1)*(na>0)
            no3 = (ecount[1] + 1)*no2
            for e3 in range(ecount[2]):
                for e2 in range(ecount[1]):
                    for e1 in range(ecount[0][na]):
                        bni = e3*no3 + e2*no2 + e1 + 1 + (cumNumNodesPerArm[na])
                        nodeIdentifiers = [ bni, bni + 1, bni + no2, bni + no2 + 1, bni + no3, bni + no3 + 1, bni + no2 + no3, bni + no2 + no3 + 1 ]
                        scr[elementIdentifier] = nodeIdentifiers
                        elementIdentifier = elementIdentifier + 1

    print('Node number (length of x) = ' + str(len(x)))

    if plot_:
        plt.axis('equal')
        plt.show()

    # return ([ix.tolist() for ix in x])
    return x, xds1, xds2, nWheel, nCentre, nWheelVertices

# Remove duplicate nodes, but keep the node correspondence - replace duplicate node in nodelist with OG node
def findDuplicateNodes(x, nodeList):

    xnew = []
    x_orig = x.copy()
    dupNodes = []
    for ir, row in enumerate(x_orig):
        iOG = x_orig.index(row)
        if row in x and ir != iOG:
            dupNodes.append([iOG+1, ir+1])
        else:
            xnew.extend(row)

    j = 10

    if False:
        if i == 0:
            x.extend(x_out)
        else:
            x2 = x.copy()
            xarr = np.array(x2)
            for row in x_out:
                if True:
                    # if ~(all(xarr == np.array(row), 1)).any():
                    if ~(xarr == np.array(row), 1)[1]:
                        x.extend([row])
                    else:
                        j = 10
                else:
                    x.extend([row])

    return dupNodes


def rotateByAngle_2D(x, th):
    xRot = [x[0] * cos(th) - x[1] * sin(th),
            x[1] * cos(th) + x[0] * sin(th),
            0]
    return xRot

if __name__ == "__main__":

    elens = [1,1, 1] # [x, y] element lengths
    numArm = 3
    ecount = [[4,2,2],2,1]
    armRotationAngles = (2*pi) / (numArm*2)
    dipMultiplier = 0.75
    armEndMult = 0.75
    plot_graph = True
    xnodes, xds1, xds2, nWheel, nCentre, nWheelVertices = createBody(elens, numArm, armRotationAngles, ecount, dipMultiplier, armEndMult, plot_graph)
    x_in = [ix[:3] for ix in xnodes]
    nodeList_sh = findDuplicateNodes(x_in, (range(1, len(xnodes) + 1)))
