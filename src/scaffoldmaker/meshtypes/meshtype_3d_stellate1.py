"""
Stellate mesh, iteration 1
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
            'Element count 2' : 2,
            'Element count 3' : 1,
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
            'Element count 2',
            'Element count 3',
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
            'Element count 2',
            'Element count 3',
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
        numArm = options['Number of arms (1 <= 4)']
        zheight = options['Element count 3']
        elementsCount1 = [options['Number of elements in arm1'], options['Number of elements in arm2'], options['Number of elements in arm3'], options['Number of elements in arm4']]
        elementsCount2 = options['Element count 2']
        elementsCount3 = options['Element count 3']
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
        xnodes, wheelNodes, centreNodes = createBody(numArm, armRotationAngles, zheight, ecount, plot_graph)
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
        xnodes_d1 = []
        dx_ds1 = [ 1.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0, 0.0 ]
        nodeIdentifier = 1

        for n2 in range(len(xnodes)):
            if nodeIdentifier not in dupNodes:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xnodes[n2])
                if xnodes_d1:
                    ds1 = xnodes_d1[n2]
                    ds2 = xnodes_d2[n2]
                elif False:
                    ds2w = [0.0, 0.0, 1.0 / elementsCount3]
                    ds1 = rotateByAngle_2D(dx_ds1, xnodes[n2][-1])
                    ds2 = rotateByAngle_2D(dx_ds2, xnodes[n2][-1])
                    if nodeIdentifier in wheelNodes:
                        ds2 = ds2w.copy()
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2w if nodeIdentifier in wheelNodes else ds2)
                else: # default
                    ds1 = rotateByAngle_2D(dx_ds1, xnodes[n2][-1])
                    ds2 = rotateByAngle_2D(dx_ds2, xnodes[n2][-1])
                if nodeIdentifier in centreNodes:
                    ds1 = 1 * [cos(armAngle/2), -sin(armAngle/2), 0.0]
                    ds2 = 1 * [cos(armAngle/2), sin(armAngle/2), 0.0]
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1) #ds1
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2) #ds2

            nodeIdentifier += 1

        # create elements
        elementIdentifier = 1
        if False:
            numArm = 1
            cumNumNodesPerArm = [0]*(numArm+1)
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


def createArm(thAr, elen, ewid, armLength, armAngle, zheight, ecount, startingNode):
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

    elementsCount1 = ecount[0]
    elementsCount2 = ecount[1]
    elementsCount3 = ecount[2]

    nodes_per_layer = (elementsCount1 + 1) * (elementsCount2 + 1)
    x = []
    dipMultiplier = 1

    if True: # 3Dbox1
        nodeIdentifier = startingNode
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
                    else:
                        x1 = n1
                        x2 = (n2 - 1) * ewid / 2
                    x.append([x1, x2, x3])
                    nodeIdentifier += 1
        e_1d = []
        # x = np.array(x)

    else:
        numN = 9 + ((armLength - 2) * 3)
        x = np.zeros([numN, 2])
        # create arms at 0degs (lying on x-axis)
        # initialise size of x given Number of elements in arm
        x[2] = [round(elen * cos(thAr), 12), round(elen * sin(thAr), 12)]
        x[1] = [round(elen * cos(thAr), 12), -round(elen * sin(thAr), 12)]

        # find nodepairs to create 1D line elements
        e_1d = [[0, 1], [2, 3], [0, 2], [1, 3], [4, 5], [0, 4], [1, 5], [1, 6], [3, 7], [1, 3], [6, 7], [5, 8], [1, 5],
                [6, 8]]
        e_1d = np.np.array(e_1d)

        for i in range(armLength): # -1
            xnew = (i + 1) * elen
            if True:
                n = i + 3
                x[n] = [xnew, -ewid / 2]
                x[n+armLength] = [xnew, 0]
                x[n+(2*armLength)] = [xnew, ewid / 2]
                # can existing 3dbox1 code for elements be used for these nodes? keep in mind the shared nodes.
            else:
                x[9+(i*3)] = [xnew, 0]
                x[10+(i*3)] = [xnew, ewid/2]
                x[11+(i*3)] = [xnew, -ewid/2]
            n_end = [6+(3*i),7+(3*i),8+(3*i)] # the nodes at the end of the 2-element arm
            e_1d_new = [[n_end[0], 9+(i*3)], [n_end[1], 10+(i*3)], [9+(i*3),10+(i*3)], [n_end[2], 11+(i*3)], [9+(i*3),11+(i*3)]]
            e_1d = np.vstack([e_1d, (np.np.array(e_1d_new))])


    # rotate entire arm about origin by armAngle
    tol = 1e-12
    for j, n in enumerate(x):
        xnew = round(n[0]*cos(armAngle), 12) - round(n[1]*sin(armAngle), 12)
        ynew = round(n[1]*cos(armAngle), 12) + round(n[0]*sin(armAngle), 12)
        if (abs(xnew) < tol):
            xnew = 0
        if (abs(ynew) < tol):
            ynew = 0

        x[j] = [xnew, ynew, x[j][2]]

    nWheel = [startingNode, \
              startingNode + 1, \
              startingNode + elementsCount1 + 2, \
              startingNode + (2*(elementsCount1 + 1)) + 1, \
              startingNode + (2*(elementsCount1 + 1))]
    nWheel = nWheel + [n+nodes_per_layer for n in nWheel]

    nCentre = startingNode + elementsCount1 + 1
    nCentre = [nCentre, nCentre + nodes_per_layer]

    return (x, nodeIdentifier, nWheel, nCentre)


def createBody(numArm, thAr, zheight, ecount, plot_):

    import numpy as np

    elen = 1
    ewid = 1

    armLength = ecount[0]
    minArmAngle = 2*pi/numArm

    x = []
    nWheel = []
    nCentre = []

    if plot_:
        plt.figure()

    nextNode = 1
    x_out = []
    for i in range(numArm):
        ecount_i = [ecount[0][i], ecount[1], ecount[2]]
        x_out, nextNode, nwhl, ncntr = createArm(thAr, elen, ewid, armLength[i], minArmAngle*i, zheight, ecount_i, nextNode)

        x.extend([ix+[minArmAngle*i] for ix in x_out])
        nWheel.extend(nwhl)
        nCentre.extend(ncntr)

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
    return x, nWheel, nCentre

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

    numArm = 3
    zheight = 0.1
    ecount = [[4,2,2],2,1]
    armRotationAngles = (2*pi) / (numArm*2)
    plot_graph = True
    xnodes, nWheel, nCentre = createBody(numArm, armRotationAngles, zheight, ecount, plot_graph)
    x_in = [ix[:3] for ix in xnodes]
    nodeList_sh = findDuplicateNodes(x_in, (range(1, len(xnodes) + 1)))
