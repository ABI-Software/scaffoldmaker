"""
Generates a 3-D stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
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
from scaffoldmaker.utils.vector import magnitude

class MeshType_3d_stellate2(Scaffold_base):
    """
    classdocs
    """
    @staticmethod
    def getName():
        return '3D Stellate 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of arms' : 3,
            'Number of elements in long arm' : 4,
            'Number of elements in short arms' : 2,
            'Element length along arm' : 1.0,
            'Element width across arm' : 0.5,
            'Element thickness' : 0.5,
            'Extent of indentation around centre' : 0.5,
            'Refine' : False,
            'Refine number of elements along arm' : 1,
            'Refine number of elements across arm' : 1,
            'Refine number of elements through thickness' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of arms',
            'Number of elements in long arm',
            'Number of elements in short arms',
            'Element length along arm',
            'Element width across arm',
            'Element thickness',
            'Extent of indentation around centre',
            'Refine',
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements in long arm',
            'Number of elements in short arms',
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness']:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Element length along arm',
            'Element width across arm',
            'Element thickness',
            'Extent of indentation around centre']:
            if options[key] < 0:
                options[key] = 0.1

        # if options['Number of arms'] < 3 or options['Number of arms'] > 4:
        #     options['Number of arms'] = 3

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        numArm = options['Number of arms']
        elementLengths = [options['Element length along arm'],
                 options['Element width across arm'],
                 options['Element thickness']]
        elementsLongArm = options['Number of elements in long arm']
        elementsShortArm = options['Number of elements in short arms']
        elementsCount1 = [elementsLongArm, elementsShortArm] + [elementsShortArm]*(numArm-2)
        elementsCount2 = options['Number of elements in short arms']
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
        numNodesPerArm = [0]*(numArm+1)
        for na in range(numArm):
            numNodesPerArm[na+1] = (elementsCount1[na] + 1) * (elementsCount2 + 1) * (elementsCount3 + 1)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        # elementsCount = [elementsCount1, elementsCount2, elementsCount3]
        armRotationAngles = (2 * math.pi) / (numArm * 2)
        dipMultiplier = 1/(options['Extent of indentation around centre']+1) #elementLengths[1] * 1.2 * 1.5

        if True:
            nodeIdentifier = 1
            armLength = elementsCount1[na]
            minArmAngle = 2 * math.pi / numArm
            xx = []
            xds1 = []
            xds2 = []
            for na in range(numArm):
                elementsCount_i = [elementsCount1[na], elementsCount2, elementsCount3]
                x, ds1, ds2, nWheelEdge = createArm(armRotationAngles, elementLengths, minArmAngle * na, minArmAngle, elementsCount_i, dipMultiplier, numArm)
                for ix in range(len(x)):
                    if na == 0 or ix not in nWheelEdge:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[ix])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1[ix])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2[ix])
                        nodeIdentifier += 1
                xx.append(x)
                xds1.append(ds1)
                xds2.append(ds2)

        else:
            xnodes, xnodes_ds1, xnodes_ds2 = createBody(elementLengths, numArm, armRotationAngles, elementsCount,dipMultiplier)
            # Remove duplicate nodes, but keep the node correspondence
            x_in = [ix[:3] for ix in xnodes]
            duplicateIndices_dbl = findduplicateIndices(x_in)
            duplicateIndices = [d[1] for d in duplicateIndices_dbl]
            nodeIdentifier = 1
            for n2 in range(len(xnodes)):
                if nodeIdentifier not in duplicateIndices:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xnodes[n2])
                    ds1 = xnodes_ds1[n2]
                    ds2 = xnodes_ds2[n2]
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2)
                nodeIdentifier += 1

        # Create elements
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives() #createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        for na in range(numArm):
            no2 = (elementsCount1[na] + 1)
            no3 = (elementsCount2 + 1)*no2 - 2
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1[na]):

                        ### NODES ###
                        offset = (cumNumNodesPerArm[na])
                        bni = e3 * no3 + e2 * no2 + e1 + 1 + offset
                        if e2 == 0:
                            if e1 == elementsCount1[na] - 1:
                                nodeIdentifiers = [bni, bni + no2 - 1, bni + no2, bni + no3,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                            else:
                                nodeIdentifiers = [bni, bni + 1, bni + no2 - 1, bni + no2, bni + no3, bni + no3 + 1,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                        else:
                            if e1 == elementsCount1[na] - 1:
                                nodeIdentifiers = [bni-1, bni, bni + no2 - 1, bni + no3 - 1, bni + no3,
                                                   bni + no2 + no3 - 1]
                            else:
                                nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no2, bni + no3 - 1, bni + no3,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                        # for ins, node in enumerate(nodeIdentifiers):
                            # if node in duplicateIndices:
                            #     iReplacementNode = duplicateIndices.index(node)
                            #     replacementNode = duplicateIndices_dbl[iReplacementNode][0]
                            #     nodeIdentifiers[ins] = replacementNode

                        if e1 == 0 and True:  # wheel
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

                        elif e1 < (elementsCount1[na] - 1) and True:
                            eft1 = eft
                            elementtemplate1 = elementtemplate
                            scalefactors = None
                        else:
                            if True:
                                # rounded ends of arms. Collapse xi2 at xi1 = 1
                                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2, [])
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


def createArm(thAr, elementLengths, armAngle, armAngleConst, elementsCount, dipMultiplier, numArm):
    """
    Create single arm unit.
    Base length of element is 1.
    Direction: anticlockwise.
    Minimum arm length is 2 elements.
    :param thAr: angle arising from base node (rad)
    :param elementLengths: [Element length along arm, half Element width across arm, Element thickness]
    :param armAngle: angle from x axis of the arm about origin (rad)
    :param armAngleConst: minimum angle from x axis of the first arm about origin (rad)
    :param elementsCount: element counts in x,y,z directions
    :param dipMultiplier: factor that wheel nodes protrude by, relative to unit length
    :param numArm: number of arms in body
    :return: x: positions of nodes in arm
    :return: nodeIdentifier: number of last node in arm +1
    :return: xnodes_ds1: ds1 derivatives of nodes associated with x
    :return: xnodes_ds2: ds2 derivatives of nodes associated with y
    """

    elementsCount1, elementsCount2, elementsCount3 = elementsCount
    [elen, ewid, zheight] = elementLengths

    xnodes_ds1 = []
    xnodes_ds2 = []
    dx_ds1 = [elen, 0.0, 0.0]
    dx_ds2 = [0.0, ewid, 0.0]
    wheelDvMult = [0.5, 0.75]

    nodes_per_layer = (elementsCount1 + 1) * (elementsCount2 + 1) - 2 # discount 2 armEnd corners
    x = []

    # nCentre and nWheelVertices are in pythonic indexing
    nCentre = elementsCount1
    nCentre = [nCentre, nCentre + nodes_per_layer]
    nWheelVertices = [0, nodes_per_layer,
                      2* (elementsCount1+1) - 1, 2* (elementsCount1+1) - 1 + nodes_per_layer ]

    dcent = [elen * math.cos(armAngleConst / 2), elen * math.sin(armAngleConst / 2), 0.0]
    dvertex = [round(elen * dipMultiplier * math.cos(thAr), 12),
                round(elen * dipMultiplier * math.sin(thAr), 12)]

    nid = 0
    for n3 in range(elementsCount3 + 1):
        x3 = n3 * zheight
        for n2 in range(elementsCount2 + 1):
            for n1 in range(elementsCount1 + 1):
                # ignore if armEnd corner nodes
                if n1 == elementsCount1 and (n2 == 0 or n2 == elementsCount2):
                    pass
                else:
                    if n1 == 0 and n2 == 0:
                        x1 = dvertex[0]
                        x2 = -dvertex[1]
                    elif n1 == 0 and n2 == elementsCount2:
                        x1 = dvertex[0]
                        x2 = dvertex[1]
                    else:
                        x1 = (elen*n1)
                        x2 = ((n2 - 1) * ewid )
                    x.append([x1, x2, x3])

                    # DERIVATIVES
                    ds1 = dx_ds1.copy()
                    ds2 = dx_ds2.copy()
                    if n2 == 1 and n1 == 0:
                        ds2 = dcent
                        ds1 = [dcent[0], -dcent[1], dcent[2]]
                    elif n1 == 0 and (n2 == 0 or n2 == elementsCount2):
                        [p, q] = [x1/dipMultiplier, x2/dipMultiplier] # unit vector
                        ds2mag = round(2*magnitude([p,q]),6) - magnitude(dcent)
                        ds2 = [ds2mag*p, ds2mag*q, 0]
                        ds1 = rotateAboutZAxis(ds2, -math.pi / 2)
                        for i in range(2):
                            ds1[i] *= wheelDvMult[0]
                        ds1 = [d * elen for d in ds1]
                    elif n1 == elementsCount1 and n2 == elementsCount2-1:
                        ds2 = [0, 1 * (elen + ewid), 0]
                    xnodes_ds1.append(ds1)
                    xnodes_ds2.append(ds2)
                    nid += 1

    # rotate entire arm about origin by armAngle except for centre nodes
    tol = 1e-12
    for j, n in enumerate(x):
        xynew = rotateAboutZAxis(n, armAngle)
        [xnew, ynew] = [round(r,12) for r in xynew][:2]
        if (abs(xnew) < tol):
            xnew = 0
        if (abs(ynew) < tol):
            ynew = 0
        x[j] = [xnew, ynew, x[j][2]]

        if (j+1) not in nCentre:
            xnodes_ds1[j] = rotateAboutZAxis(xnodes_ds1[j], armAngle)
            xnodes_ds2[j] = rotateAboutZAxis(xnodes_ds2[j], armAngle)

    return (x, xnodes_ds1, xnodes_ds2, nCentre+nWheelVertices)


def createBody(elementLengths, numArm, thAr, elementsCount, dipMultiplier):
    """
    Create amalgamation of arms about central junction.
    :param elementLengths: [Element length along arm, half Element width across arm, Element thickness]
    :param numArm: number of arms in body
    :param thAr: angle arising from base node (rad)
    :param elementsCount: element counts in x,y,z directions
    :param dipMultiplier: factor that wheel nodes protrude by, relative to unit length
    :return: x: positions of nodes in body
    :return: xds1: ds1 derivatives of nodes associated with x
    :return: xds2: ds2 derivatives of nodes associated with x
    """
    armLength = elementsCount[0]
    minArmAngle = 2*math.pi/numArm

    x = []
    xds1 = []
    xds2 = []

    nextNode = 1
    x_out = []
    for i in range(numArm):
        elementsCount_i = [elementsCount[0][i], elementsCount[1], elementsCount[2]]
        x_out, nextNode, ds1, ds2 = createArm(thAr, elementLengths, armLength[i], minArmAngle*i, minArmAngle, elementsCount_i, nextNode, dipMultiplier, numArm)

        x.extend([ix+[minArmAngle*i] for ix in x_out])
        xds1.extend(ds1)
        xds2.extend(ds2)

    return x, xds1, xds2


def findduplicateIndices(x):
    """
    Find repeated rows in x and return their indices.
        :param x: list
        :return: duplicateIndices: [index of first instance of row, index of where duplicate occurred]
    """
    xnew = []
    x_orig = x.copy()
    duplicateIndices = []
    for ir, row in enumerate(x_orig):
        iOG = x_orig.index(row)
        if row in x and ir != iOG:
            duplicateIndices.append([iOG+1, ir+1])
        else:
            xnew.extend(row)
    return duplicateIndices
