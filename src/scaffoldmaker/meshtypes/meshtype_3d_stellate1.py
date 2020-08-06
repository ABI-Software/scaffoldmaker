"""
Generates a 3-D unit box mesh with variable numbers of elements in 3 directions.
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
            'Number of elements 1' : 1,
            'Number of elements 2' : 1,
            'Number of elements 3' : 1,
            'Use cross derivatives' : False,
            'Refine' : False,
            'Refine number of elements 1' : 1,
            'Refine number of elements 2' : 1,
            'Refine number of elements 3' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements 1',
            'Refine number of elements 2',
            'Refine number of elements 3'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements 1',
            'Number of elements 2',
            'Number of elements 3',
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
        elementsCount1 = options['Number of elements 1']
        elementsCount2 = options['Number of elements 2']
        elementsCount3 = options['Number of elements 3']
        useCrossDerivatives = options['Use cross derivatives']

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
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 1.0 / elementsCount1, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0 / elementsCount2, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 1.0 / elementsCount3 ]
        zero = [ 0.0, 0.0, 0.0 ]
        for n3 in range(elementsCount3 + 1):
            x[2] = n3 / elementsCount3
            for n2 in range(elementsCount2 + 1):
                x[1] = n2 / elementsCount2
                for n1 in range(elementsCount1 + 1):
                    x[0] = n1 / elementsCount1
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        no2 = (elementsCount1 + 1)
        no3 = (elementsCount2 + 1)*no2
        for e3 in range(elementsCount3):
            for e2 in range(elementsCount2):
                for e1 in range(elementsCount1):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni = e3*no3 + e2*no2 + e1 + 1
                    nodeIdentifiers = [ bni, bni + 1, bni + no2, bni + no2 + 1, bni + no3, bni + no3 + 1, bni + no2 + no3, bni + no2 + no3 + 1 ]
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


def createArm(thAr, elen, ewid, armLength, armAngle):
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

    numN = 9 + ((armLength-1)*3)
    x = np.zeros([numN,2])

    # create arms at 0degs (lying on x-axis)
    x[0] = [0,0]
    x[1] = [elen, 0]
    x[2] = [elen*cos(thAr), elen*sin(thAr)]
    x[3] = [1,ewid/2]

    x[4] = [elen*cos(thAr), -elen*sin(thAr)]
    x[5] = [1,-ewid/2]

    x[6] = [elen*2, 0]
    x[7] = [elen*2, ewid/2]
    x[8] = [x[7][0], -x[7][1]]

    # find nodepairs to create 1D line elements
    e_1d = [[0,1],[2,3],[0,2],[1,3],[4,5],[0,4],[1,5],[1,6],[3,7],[1,3],[6,7],[5,8],[1,5],[6,8]]
    e_1d = np.array(e_1d)

    for i in range(armLength-1):
        xnew = (i+2)*elen
        x[9+(i*3)] = [xnew, 0]
        x[10+(i*3)] = [xnew, ewid/2]
        x[11+(i*3)] = [xnew, -ewid/2]
        n_end = [6+(3*i),7+(3*i),8+(3*i)] # the nodes at the end of the 2-element arm
        e_1d_new = [[n_end[0], 9+(i*3)], [n_end[1], 10+(i*3)], [9+(i*3),10+(i*3)], [n_end[2], 11+(i*3)], [9+(i*3),11+(i*3)]]
        e_1d = np.vstack([e_1d, (np.array(e_1d_new))])

    # rotate entire arm about origin by armAngle
    for j, n in enumerate(x):
        xnew = n[0]*cos(armAngle) - n[1]*sin(armAngle)
        ynew = n[1]*cos(armAngle) + n[0]*sin(armAngle)
        x[j] = [xnew, ynew]

    if False:
        plt.figure()
        plt.plot(x[:,0], x[:,1],'x')
        plt.plot(x[[0,1,6],0], x[[0,1,6],1], '-') # line
        plt.show()

    return (x, e_1d)

def createShape(plot_):

    numArm = 3
    thAr = (2*pi) / (numArm*2)
    elen = 1
    ewid = 0.8

    armLength = [4,2,2]
    armAngle = [0,(2*pi/3), (4*pi/3)]

    x = []

    if plot_:
        plt.figure()

    for i in range(numArm):

        x_out, e_1d_out = createArm(thAr, elen, ewid, armLength[i], armAngle[i])

        # only append points that do not exist already? What about elements?
        x.append(x_out)

        if plot_:
            plt.plot(x_out[:, 0], x_out[:, 1], 'x')
            plt.plot(x_out[[0, 1, 6], 0], x_out[[0, 1, 6], 1], '-')  # line
            for np in e_1d_out:
                plt.plot(x_out[np][:,0], x_out[np][:,1], 'k-')
            plt.axis('equal')

    if plot_:
        plt.show()


if __name__ == "__main__":
    plot_graph = True
    createShape(plot_graph)