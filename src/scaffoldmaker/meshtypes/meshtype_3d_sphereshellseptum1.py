"""
Generates a 3-D 'sphere shell septum' mesh crossing from one side of a sphere
shell to the other, including elements for the outside, suitable for joining
two halves of a sphere shell.
It is the middle line in (|).
The number of elements up the sphere and across the septum can be varied.
Only one element throught the wall is currently implemented.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.interpolation import interpolateCubicHermite, interpolateCubicHermiteDerivative


class MeshType_3d_sphereshellseptum1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Sphere Shell Septum 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements up' : 4,
            'Number of elements across' : 2,
            'Wall thickness left' : 0.25,
            'Wall thickness right' : 0.25,
            'Flange length' : 0.125,
            'Bulge radius' : 0.0,
            'Apex scale factor identifier offset' : 10000,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements across',
            'Wall thickness left',
            'Wall thickness right',
            'Flange length',
            'Bulge radius',
            'Apex scale factor identifier offset',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements up'] < 2) :
            options['Number of elements up'] = 2
        if (options['Number of elements across'] < 2) :
            options['Number of elements across'] = 2
        if (options['Wall thickness left'] < 0.0) :
            options['Wall thickness left'] = 0.0
        elif (options['Wall thickness left'] > 0.5) :
            options['Wall thickness left'] = 0.5
        if (options['Wall thickness right'] < 0.0) :
            options['Wall thickness right'] = 0.0
        elif (options['Wall thickness right'] > 0.5) :
            options['Wall thickness right'] = 0.5
        if (options['Apex scale factor identifier offset'] < 0) :
            options['Apex scale factor identifier offset'] = 0
        if (options['Flange length'] < 0.0) :
            options['Flange length'] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountUp = options['Number of elements up']
        elementsCountAcross = options['Number of elements across']
        wallThickness = [ options['Wall thickness left'], options['Wall thickness right'] ]
        flangeLength = options['Flange length']
        bulgeRadius = options['Bulge radius']
        apexScaleFactorIdentifierOffset = options['Apex scale factor identifier offset']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        else:
            nodetemplate = nodetemplateApex

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()
        eftOuter = tricubichermite.createEftTubeSeptumOuter()
        eftInner1 = tricubichermite.createEftTubeSeptumInner1()
        eftInner2 = tricubichermite.createEftTubeSeptumInner2()

        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eftOuterApex0 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        # general linear map at 4 nodes for one derivative
        eftOuterApex0.setNumberOfLocalScaleFactors(8)
        for s in range(8):
            eftOuterApex0.setScaleFactorType(s + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            scaleFactorId = (s % 2) + 1
            if s < 4:
                scaleFactorId += apexScaleFactorIdentifierOffset
            eftOuterApex0.setScaleFactorIdentifier(s + 1, scaleFactorId)
        if useCrossDerivatives:
            noCrossRange = range(4)
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eftOuterApex0.setFunctionNumberOfTerms(n*8 + 4, 0)
            eftOuterApex0.setFunctionNumberOfTerms(n*8 + 6, 0)
            eftOuterApex0.setFunctionNumberOfTerms(n*8 + 7, 0)
            eftOuterApex0.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi1 from ds1 + ds3 for first 4 nodes
        for n in range(4):
            ln = n + 1
            eftOuterApex0.setFunctionNumberOfTerms(n*8 + 2, 2)
            eftOuterApex0.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftOuterApex0.setTermScaling(n*8 + 2, 1, [n*2 + 1])
            eftOuterApex0.setTermNodeParameter(n*8 + 2, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftOuterApex0.setTermScaling(n*8 + 2, 2, [n*2 + 2])
        #print('eftOuterApex0', eftOuterApex0.validate())

        eftOuterApex1 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftOuterApex2 = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        i = 0
        for eftOuterApex in [ eftOuterApex1, eftOuterApex2 ]:
            i += 1
            # general linear map at 4 nodes for one derivative
            eftOuterApex.setNumberOfLocalScaleFactors(9)
            # GRC: allow scale factor identifier for global -1.0 to be prescribed
            eftOuterApex.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
            eftOuterApex.setScaleFactorIdentifier(1, 1)
            for s in range(8):
                eftOuterApex.setScaleFactorType(s + 2, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                scaleFactorId = (s % 2) + 1
                if ((i == 1) and (s < 4)) or ((i == 2) and (s >= 4)):
                    scaleFactorId += apexScaleFactorIdentifierOffset
                result = eftOuterApex.setScaleFactorIdentifier(s + 2, scaleFactorId)
            if useCrossDerivatives:
                noCrossRange = range(4)
            else:
                noCrossRange = range(8)
            for n in noCrossRange:
                eftOuterApex.setFunctionNumberOfTerms(n*8 + 4, 0)
                eftOuterApex.setFunctionNumberOfTerms(n*8 + 6, 0)
                eftOuterApex.setFunctionNumberOfTerms(n*8 + 7, 0)
                eftOuterApex.setFunctionNumberOfTerms(n*8 + 8, 0)
            # general linear map dxi1 from ds1 + ds3 for first 4 nodes
            for n in range(4):
                ln = n + 1
                eftOuterApex.setFunctionNumberOfTerms(n*8 + 2, 2)
                eftOuterApex.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                negate = (i == 1) and (n < 2)  # False  # ((i == 1) and (n < 2)) or ((i == 2) and (n >= 2))
                eftOuterApex.setTermScaling(n*8 + 2, 1, [1, n*2 + 2] if negate else [n*2 + 2])
                eftOuterApex.setTermNodeParameter(n*8 + 2, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftOuterApex.setTermScaling(n*8 + 2, 2, [1, n*2 + 3] if negate else [n*2 + 3])
            if i == 1:
                negateNodes1 = [ 4, 5 ]
                negateNodes2 = [ 0, 1, 4, 5 ]
            else:
                negateNodes1 = [ 6, 7 ]
                negateNodes2 = [ 2, 3, 6, 7]
            for n in negateNodes1:
                eftOuterApex.setTermScaling(n*8 + 2, 1, [1])
            for n in negateNodes2:
                eftOuterApex.setTermScaling(n*8 + 3, 1, [1])
        #print('eftOuterApex1', eftOuterApex1.validate())
        #print('eftOuterApex2', eftOuterApex2.validate())

        # Inner Apex 1: collapsed on xi2 = 0
        eftInnerApex1a = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftInnerApex1b = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftInnerApex1c = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        for eftInnerApex in [ eftInnerApex1a, eftInnerApex1b, eftInnerApex1c ]:
            #print('**** eftInnerApex ****')
            eftInnerApex.setNumberOfLocalNodes(6)
            if eftInnerApex is eftInnerApex1b:
                eftInnerApex.setNumberOfLocalScaleFactors(18)
            else:
                eftInnerApex.setNumberOfLocalScaleFactors(22)
            # GRC: allow scale factor identifier for global -1.0 to be prescribed
            eftInnerApex.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
            eftInnerApex.setScaleFactorIdentifier(1, 1)
            # Global scale factor 4.0 used for cross derivative term to correct first derivative
            eftInnerApex.setScaleFactorType(2, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
            eftInnerApex.setScaleFactorIdentifier(2, 2)
            for layer in range(2):
                so = 2 + layer*(8 if (eftInnerApex is eftInnerApex1b) else 10)
                no = layer*3
                fo = layer*32
                for s in range(6):
                    si = so + s + 1
                    # 3 scale factors per node*direction in xi1, xi2: cos(theta), sin(theta), arc angle radians
                    sid = apexScaleFactorIdentifierOffset + (s // 3)*100 + (s % 3) + 101  # add 100 for different 'version'
                    eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                    eftInnerApex.setScaleFactorIdentifier(si, sid)
                    #print('scalefactor ', si, ' identifier', sid)
                for s in range(2):
                    si = so + s + 7
                    # 2 scale factors per node in xi3: cos(phi), sin(phi)
                    sid = apexScaleFactorIdentifierOffset + s + 1
                    eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                    eftInnerApex.setScaleFactorIdentifier(si, sid)
                    #print('scalefactor ', si, ' identifier', sid)
                if eftInnerApex is not eftInnerApex1b:
                    for s in range(2):
                        si = so + s + 9
                        # 2 scale factors per node in xi3: cos(phi), sin(phi)
                        sid = s + 1  # add 100 for different 'version'
                        eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                        eftInnerApex.setScaleFactorIdentifier(si, sid)
                        #print('scalefactor ', si, ' identifier', sid)

                # basis node 1 -> local node 1
                ln = no + 1
                eftInnerApex.setTermNodeParameter(fo + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                # 0 terms = zero parameter for d/dxi1 basis
                eftInnerApex.setFunctionNumberOfTerms(fo + 2, 0)
                # 2 terms for d/dxi2 via general linear map:
                eftInnerApex.setFunctionNumberOfTerms(fo + 3, 2)
                eftInnerApex.setTermNodeParameter(fo + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 3, 1, [so + 1])
                eftInnerApex.setTermNodeParameter(fo + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo + 3, 2, [so + 2])
                # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
                eftInnerApex.setFunctionNumberOfTerms(fo + 4, 2)
                eftInnerApex.setTermNodeParameter(fo + 4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 4, 1, [so + 2, so + 3])
                eftInnerApex.setTermNodeParameter(fo + 4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo + 4, 2, [1, so + 1, so + 3])
                # 2 terms for d/dx3 via general linear map
                eftInnerApex.setFunctionNumberOfTerms(fo + 5, 2)
                eftInnerApex.setTermNodeParameter(fo + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 5, 1, [1, so + 7])
                eftInnerApex.setTermNodeParameter(fo + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo + 5, 2, [1, so + 8])
                # 0 terms = zero parameter for cross derivative 1 3
                eftInnerApex.setFunctionNumberOfTerms(fo + 6, 0)
                # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                # GRC WIP
                #eftInnerApex.setFunctionNumberOfTerms(fo + 7, 2)
                #eftInnerApex.setTermNodeParameter(fo + 7, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                #eftInnerApex.setTermScaling(fo + 7, 1, [1, 2, so + 7] if (layer == 0) else [2, so + 7])
                #eftInnerApex.setTermNodeParameter(fo + 7, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                #eftInnerApex.setTermScaling(fo + 7, 2, [1, 2, so + 8] if (layer == 0) else [2, so + 8])
                eftInnerApex.setFunctionNumberOfTerms(fo + 7, 1)
                eftInnerApex.setTermNodeParameter(fo + 7, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo + 7, 1, [2, so + 1])
                # 0 terms = zero parameter for cross derivative 1 2 3
                eftInnerApex.setFunctionNumberOfTerms(fo + 8, 0)

                # basis node 2 -> local node 1
                eftInnerApex.setTermNodeParameter(fo + 9, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                # 0 terms = zero parameter for d/dxi1 basis
                eftInnerApex.setFunctionNumberOfTerms(fo + 10, 0)
                # 2 terms for d/dxi2 via general linear map:
                eftInnerApex.setFunctionNumberOfTerms(fo + 11, 2)
                eftInnerApex.setTermNodeParameter(fo + 11, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 11, 1, [so + 4])
                eftInnerApex.setTermNodeParameter(fo + 11, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo + 11, 2, [so + 5])
                # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
                eftInnerApex.setFunctionNumberOfTerms(fo + 12, 2)
                eftInnerApex.setTermNodeParameter(fo + 12, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 12, 1, [so + 5, so + 6])
                eftInnerApex.setTermNodeParameter(fo + 12, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo + 12, 2, [1, so + 4, so + 6])
                # 2 terms for d/dx3 via general linear map
                eftInnerApex.setFunctionNumberOfTerms(fo + 13, 2)
                eftInnerApex.setTermNodeParameter(fo + 13, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo + 13, 1, [1, so + 7])
                eftInnerApex.setTermNodeParameter(fo + 13, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                eftInnerApex.setTermScaling(fo + 13, 2, [1, so + 8])
                # 0 terms = zero parameter for cross derivative 1 3
                eftInnerApex.setFunctionNumberOfTerms(fo + 14, 0)
                # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                # GRC WIP
                #eftInnerApex.setFunctionNumberOfTerms(fo + 15, 2)
                #eftInnerApex.setTermNodeParameter(fo + 15, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                #eftInnerApex.setTermScaling(fo + 15, 1, [2, so + 7] if (layer == 0) else [1, 2, so + 7])
                #eftInnerApex.setTermNodeParameter(fo + 15, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                #eftInnerApex.setTermScaling(fo + 15, 2, [2, so + 8] if (layer == 0) else [1, 2, so + 8])
                eftInnerApex.setFunctionNumberOfTerms(fo + 15, 1)
                eftInnerApex.setTermNodeParameter(fo + 15, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo + 15, 1, [2, so + 4])
                # 0 terms = zero parameter for cross derivative 1 2 3
                eftInnerApex.setFunctionNumberOfTerms(fo + 16, 0)

                # basis nodes 3, 4 -> regular local nodes 2, 3
                for bn in range(2,4):
                    fo2 = fo + bn*8
                    ln = no + bn
                    eftInnerApex.setTermNodeParameter(fo2 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 3, 1, ln, Node.VALUE_LABEL_D_DS2, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    if useCrossDerivatives:
                        eftInnerApex.setTermNodeParameter(fo2 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D2_DS1DS3, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 7, 1, ln, Node.VALUE_LABEL_D2_DS2DS3, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 8, 1, ln, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                    else:
                         eftInnerApex.setFunctionNumberOfTerms(fo2 + 4, 0)
                         eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 0)
                         eftInnerApex.setFunctionNumberOfTerms(fo2 + 7, 0)
                         eftInnerApex.setFunctionNumberOfTerms(fo2 + 8, 0)

                # correct inner apex 1a and 1c which general linear map d/dxi3 and reverse d/dxi1 on 2nd layer
                if eftInnerApex is eftInnerApex1a:
                    fo2 = fo + 2*8
                    ln = no + 2
                    if layer == 1:
                        eftInnerApex.setTermScaling(fo2 + 2, 1, [1])
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 5, 2)
                    result1 = eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    result2 = eftInnerApex.setTermScaling(fo2 + 5, 1, [so + 9])
                    result3 = eftInnerApex.setTermNodeParameter(fo2 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                    result4 = eftInnerApex.setTermScaling(fo2 + 5, 2, [so + 10])
                    # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    eftInnerApex.setTermScaling(fo2 + 6, 1, [1, 2] if (layer == 0) else [2])
                elif eftInnerApex is eftInnerApex1c:
                    fo2 = fo + 3*8
                    ln = no + 3
                    if layer == 1:
                        eftInnerApex.setTermScaling(fo2 + 2, 1, [1])
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 5, 2)
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    eftInnerApex.setTermScaling(fo2 + 5, 1, [1, so + 9])
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                    eftInnerApex.setTermScaling(fo2 + 5, 2, [1, so + 10])
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    eftInnerApex.setTermScaling(fo2 + 6, 1, [2] if (layer == 0) else [1, 2])
            #print('eftInnerApex1.validate()', eftInnerApex.validate())

        # Inner Apex 2: collapsed on xi2 = 2
        eftInnerApex2a = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftInnerApex2b = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        eftInnerApex2c = mesh.createElementfieldtemplate(tricubicHermiteBasis)
        for eftInnerApex in [ eftInnerApex2a, eftInnerApex2b, eftInnerApex2c ]:
            #print('**** eftInnerApex ****')
            eftInnerApex.setNumberOfLocalNodes(6)
            if eftInnerApex is eftInnerApex2b:
                eftInnerApex.setNumberOfLocalScaleFactors(18)
            else:
                eftInnerApex.setNumberOfLocalScaleFactors(22)
            # GRC: allow scale factor identifier for global -1.0 to be prescribed
            eftInnerApex.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
            eftInnerApex.setScaleFactorIdentifier(1, 1)
            # Global scale factor 4.0 used for cross derivative term to correct first derivative
            eftInnerApex.setScaleFactorType(2, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
            eftInnerApex.setScaleFactorIdentifier(2, 2)
            for layer in range(2):
                so = 2 + layer*(8 if (eftInnerApex is eftInnerApex2b) else 10)
                no = layer*3
                fo = layer*32
                for s in range(6):
                    si = so + s + 1
                    # 3 scale factors per node*direction in xi1, xi2: cos(theta), sin(theta), arc angle radians
                    sid = apexScaleFactorIdentifierOffset + (s // 3)*100 + (s % 3) + 101  # add 100 for different 'version'
                    eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                    eftInnerApex.setScaleFactorIdentifier(si, sid)
                    #print('scalefactor ', si, ' identifier', sid)
                for s in range(2):
                    si = so + s + 7
                    # 2 scale factors per node in xi3: cos(phi), sin(phi)
                    sid = apexScaleFactorIdentifierOffset + s + 1
                    eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                    eftInnerApex.setScaleFactorIdentifier(si, sid)
                    #print('scalefactor ', si, ' identifier', sid)
                if eftInnerApex is not eftInnerApex2b:
                    for s in range(2):
                        si = so + s + 9
                        # 2 scale factors per node in xi3: cos(phi), sin(phi)
                        sid = s + 1  # add 100 for different 'version'
                        eftInnerApex.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
                        eftInnerApex.setScaleFactorIdentifier(si, sid)
                        #print('scalefactor ', si, ' identifier', sid)

                # basis nodes 1, 2 -> regular local nodes 1, 2 (for each layer)
                for bn in range(2):
                    fo2 = fo + bn*8
                    ln = no + bn + 1
                    eftInnerApex.setTermNodeParameter(fo2 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 3, 1, ln, Node.VALUE_LABEL_D_DS2, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    if useCrossDerivatives:
                        eftInnerApex.setTermNodeParameter(fo2 + 4, 1, ln, Node.VALUE_LABEL_D2_DS1DS2, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D2_DS1DS3, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 7, 1, ln, Node.VALUE_LABEL_D2_DS2DS3, 1)
                        eftInnerApex.setTermNodeParameter(fo2 + 8, 1, ln, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
                    else:
                        eftInnerApex.setFunctionNumberOfTerms(fo2 + 4, 0)
                        eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 0)
                        eftInnerApex.setFunctionNumberOfTerms(fo2 + 7, 0)
                        eftInnerApex.setFunctionNumberOfTerms(fo2 + 8, 0)

                # correct inner apex 2a and 2c which general linear map d/dxi3 and reverse d/dxi1 on 2nd layer
                if eftInnerApex is eftInnerApex2a:
                    fo2 = fo
                    ln = no + 1
                    if layer == 1:
                        eftInnerApex.setTermScaling(fo2 + 2, 1, [1])
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 5, 2)
                    result1 = eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    result2 = eftInnerApex.setTermScaling(fo2 + 5, 1, [so + 9])
                    result3 = eftInnerApex.setTermNodeParameter(fo2 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                    result4 = eftInnerApex.setTermScaling(fo2 + 5, 2, [so + 10])
                    # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                    # GRC WIP
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    eftInnerApex.setTermScaling(fo2 + 6, 1, [1, 2] if (layer == 0) else [2])
                elif eftInnerApex is eftInnerApex2c:
                    fo2 = fo + 8
                    ln = no + 2
                    if layer == 1:
                        eftInnerApex.setTermScaling(fo2 + 2, 1, [1])
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 5, 2)
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                    eftInnerApex.setTermScaling(fo2 + 5, 1, [1, so + 9])
                    eftInnerApex.setTermNodeParameter(fo2 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                    eftInnerApex.setTermScaling(fo2 + 5, 2, [1, so + 10])
                    # GRC WIP
                    eftInnerApex.setFunctionNumberOfTerms(fo2 + 6, 1)
                    eftInnerApex.setTermNodeParameter(fo2 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                    eftInnerApex.setTermScaling(fo2 + 6, 1, [2] if (layer == 0) else [1, 2])

                # basis node 3 -> local node 3
                ln = no + 3
                fo3 = fo + 16
                eftInnerApex.setTermNodeParameter(fo3 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                # 0 terms = zero parameter for d/dxi1 basis
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 2, 0)
                # 2 terms for d/dxi2 via general linear map:
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 3, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 3, 1, [so + 1])
                eftInnerApex.setTermNodeParameter(fo3 + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo3 + 3, 2, [so + 2])
                # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 4, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 4, 1, [1, so + 2, so + 3])
                eftInnerApex.setTermNodeParameter(fo3 + 4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo3 + 4, 2, [so + 1, so + 3])
                # 2 terms for d/dx3 via general linear map
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 5, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 5, 1, [1, so + 7])
                eftInnerApex.setTermNodeParameter(fo3 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo3 + 5, 2, [1, so + 8])
                # 0 terms = zero parameter for cross derivative 1 3
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 6, 0)
                # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                # GRC WIP
                #eftInnerApex.setFunctionNumberOfTerms(fo3 + 7, 2)
                #eftInnerApex.setTermNodeParameter(fo3 + 7, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                #eftInnerApex.setTermScaling(fo3 + 7, 1, [1, 2, so + 7] if (layer == 0) else [2, so + 7])
                #eftInnerApex.setTermNodeParameter(fo3 + 7, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                #eftInnerApex.setTermScaling(fo3 + 7, 2, [1, 2, so + 8] if (layer == 0) else [2, so + 8])
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 7, 1)
                eftInnerApex.setTermNodeParameter(fo3 + 7, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo3 + 7, 1, [1, 2, so + 1])
                # 0 terms = zero parameter for cross derivative 1 2 3
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 8, 0)

                # basis node 4 -> local node 3
                eftInnerApex.setTermNodeParameter(fo3 + 9, 1, ln, Node.VALUE_LABEL_VALUE, 1)
                # 0 terms = zero parameter for d/dxi1 basis
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 10, 0)
                # 2 terms for d/dxi2 via general linear map:
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 11, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 11, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 11, 1, [so + 4])
                eftInnerApex.setTermNodeParameter(fo3 + 11, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo3 + 11, 2, [so + 5])
                # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 12, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 12, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 12, 1, [1, so + 5, so + 6])
                eftInnerApex.setTermNodeParameter(fo3 + 12, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
                eftInnerApex.setTermScaling(fo3 + 12, 2, [so + 4, so + 6])
                # 2 terms for d/dx3 via general linear map
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 13, 2)
                eftInnerApex.setTermNodeParameter(fo3 + 13, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                eftInnerApex.setTermScaling(fo3 + 13, 1, [1, so + 7])
                eftInnerApex.setTermNodeParameter(fo3 + 13, 2, ln, Node.VALUE_LABEL_D_DS3, 1) 
                eftInnerApex.setTermScaling(fo3 + 13, 2, [1, so + 8])
                # 0 terms = zero parameter for cross derivative 1 3
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 14, 0)
                # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
                # GRC WIP
                #eftInnerApex.setFunctionNumberOfTerms(fo3 + 15, 2)
                #eftInnerApex.setTermNodeParameter(fo3 + 15, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
                #eftInnerApex.setTermScaling(fo3 + 15, 1, [2, so + 7] if (layer == 0) else [1, 2, so + 7])
                #eftInnerApex.setTermNodeParameter(fo3 + 15, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
                #eftInnerApex.setTermScaling(fo3 + 15, 2, [2, so + 8] if (layer == 0) else [1, 2, so + 8])
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 15, 1)
                eftInnerApex.setTermNodeParameter(fo3 + 15, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
                eftInnerApex.setTermScaling(fo3 + 15, 1, [1, 2, so + 4])
                # 0 terms = zero parameter for cross derivative 11 2 3
                eftInnerApex.setFunctionNumberOfTerms(fo3 + 16, 0)
            #print('eftInnerApex2.validate()', eftInnerApex.validate())

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateOuter = mesh.createElementtemplate()
        elementtemplateOuter.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuter.defineField(coordinates, -1, eftOuter)

        elementtemplateOuterApex0 = mesh.createElementtemplate()
        elementtemplateOuterApex0.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuterApex0.defineField(coordinates, -1, eftOuterApex0)
        elementtemplateOuterApex1 = mesh.createElementtemplate()
        elementtemplateOuterApex1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuterApex1.defineField(coordinates, -1, eftOuterApex1)
        elementtemplateOuterApex2 = mesh.createElementtemplate()
        elementtemplateOuterApex2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateOuterApex2.defineField(coordinates, -1, eftOuterApex2)
        #print(result, 'elementtemplateOuterApex2.defineField')

        elementtemplateInner1 = mesh.createElementtemplate()
        elementtemplateInner1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInner1.defineField(coordinates, -1, eftInner1)
        elementtemplateInner2 = mesh.createElementtemplate()
        elementtemplateInner2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInner2.defineField(coordinates, -1, eftInner2)

        elementtemplateInnerApex1a = mesh.createElementtemplate()
        elementtemplateInnerApex1a.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex1a.defineField(coordinates, -1, eftInnerApex1a)
        elementtemplateInnerApex1b = mesh.createElementtemplate()
        elementtemplateInnerApex1b.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex1b.defineField(coordinates, -1, eftInnerApex1b)
        elementtemplateInnerApex1c = mesh.createElementtemplate()
        elementtemplateInnerApex1c.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex1c.defineField(coordinates, -1, eftInnerApex1c)

        elementtemplateInnerApex2a = mesh.createElementtemplate()
        elementtemplateInnerApex2a.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex2a.defineField(coordinates, -1, eftInnerApex2a)
        elementtemplateInnerApex2b = mesh.createElementtemplate()
        elementtemplateInnerApex2b.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex2b.defineField(coordinates, -1, eftInnerApex2b)
        elementtemplateInnerApex2c = mesh.createElementtemplate()
        elementtemplateInnerApex2c.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateInnerApex2c.defineField(coordinates, -1, eftInnerApex2c)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementUp = math.pi/elementsCountUp
        radiansPerElementAcross = math.pi/elementsCountAcross
        bevelAngle = math.pi/2
        sinBevelAngle = math.sin(bevelAngle)
        cosBevelAngle = math.cos(bevelAngle)
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        wallThicknessSeptum = wallThickness[0]
        #wallThicknessMin = min(wallThickness)
        wallThicknessMax = max(wallThickness)
        for n3 in range(2):
            sign = -1.0 if (n3 == 0) else 1.0
            radiusY = 0.5*wallThicknessSeptum
            radiusX = 0.5 - wallThicknessMax

            # create two bottom apex nodes
            radius = radiusX + wallThickness[n3]
            x[0] = 0.0
            x[1] = -sign*(radiusY + flangeLength)
            #if wallThickness[0] > wallThickness[1]:
            #    x[1] -= sideBulge
            #elif wallThickness[0] < wallThickness[1]:
            #    x[1] += sideBulge
            x[2] = -radius
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, wallThicknessSeptum, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, -wallThickness[n3] ])
            nodeIdentifier = nodeIdentifier + 1

            radius = radiusX
            x[0] = 0.0
            x[1] = -sign*(radiusY + flangeLength)
            x[2] = -radius
            dx_ds3[0] = 0.0
            if ((n3 == 1) and (wallThickness[0] > wallThickness[1])) or \
                ((n3 == 0) and (wallThickness[0] < wallThickness[1])):
                dx_ds3[1] = 0.0
                dx_ds3[2] = -wallThickness[n3]
            else:
                dx_ds3[1] = sign*wallThickness[n3]*cosBevelAngle
                dx_ds3[2] = -wallThickness[n3]*sinBevelAngle
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            mag1 = 2.0*radius/elementsCountAcross
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, mag1*sinBevelAngle, sign*mag1*cosBevelAngle ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier = nodeIdentifier + 1

            # create regular rows between apexes
            for n2 in range(1, elementsCountUp):
                xi2 = n2/elementsCountUp
                radiansUp = n2*radiansPerElementUp
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                x[2] = 0.0  # GRC: must set
                for n1 in range(elementsCountAcross + 3):
                    if (n1 == 0) or (n1 == (elementsCountAcross + 2)):
                        flip = -1.0 if (n1 == 0) else 1.0
                        radius = radiusX + wallThickness[n3]
                        x = [ radius*sinRadiansUp, -sign*(radiusY + flangeLength), -radius*cosRadiansUp ]
                        if n1 == 0:
                            x[0] = -x[0]
                        #if wallThickness[0] > wallThickness[1]:
                        #    x[1] -= sideBulge
                        #elif wallThickness[0] < wallThickness[1]:
                        #    x[1] += sideBulge
                        dx_ds1 = [ 0.0, flip*wallThicknessSeptum, 0.0 ]
                        mag2 = radius*radiansPerElementUp
                        dx_ds2 = [ mag2*cosRadiansUp*flip, 0.0, mag2*sinRadiansUp ]
                        mag3 = wallThickness[n3]*flip
                        dx_ds3 = [ mag3*sinRadiansUp, 0.0, -mag3*cosRadiansUp*flip ]
                    elif (n1 == 1) or (n1 == (elementsCountAcross + 1)):
                        flip = -1.0 if (n1 == 1) else 1.0
                        radius = radiusX
                        x[0] = flip*sinRadiansUp*radius
                        x[1] = sign*(-radiusY - flangeLength)
                        x[2] = -cosRadiansUp*radius


                        mag1x = radiusX/elementsCountAcross
                        mag1y = flangeLength
                        #mag1min = radiusX/elementsCountAcross + flangeLength + math.sqrt(mag1x*mag1x + mag1y*mag1y)
                        mag1min = 2.0*math.sqrt(mag1x*mag1x + mag1y*mag1y)
                        mag1max = 2.0*(radiusY + flangeLength)
                        mag1 = mag1min if (mag1min < mag1max) else mag1max
                        dx_ds1 = [ -sign*sinRadiansUp*mag1*cosBevelAngle, flip*mag1*sinBevelAngle, flip*sign*cosRadiansUp*mag1*cosBevelAngle ]
                        mag2 = radius*radiansPerElementUp
                        dx_ds2 = [ mag2*cosRadiansUp*flip, 0.0, mag2*sinRadiansUp ]
                        mag3 = wallThickness[n3]
                        #if ((n3 == 1) and (wallThickness[0] > wallThickness[1])) or \
                        #    ((n3 == 0) and (wallThickness[0] < wallThickness[1])):
                        dx_ds3 = [ flip*sinRadiansUp*mag3, 0.0, -cosRadiansUp*mag3 ]
                        #else:
                        #    dx_ds3 = [ flip*sinRadiansUp*mag3*sinBevelAngle, sign*mag3*cosBevelAngle, -cosRadiansUp*mag3*sinBevelAngle ]
                        if n1 == 1:
                            # Prepare for interpolating interior points
                            v1 = ( -sinRadiansUp*radius, -sign*radiusY, -cosRadiansUp*radius )
                            mag = math.sqrt(dx_ds1[0]*dx_ds1[0] + dx_ds1[2]*dx_ds1[2])
                            scale1 = -sign/mag*radius*sinRadiansUp
                            d1 = ( scale1*dx_ds1[0], 0.0, scale1*dx_ds1[2] )
                            v2 = ( 0.0, -sign*radiusY, 2.0*(xi2 - 0.5)*radius )
                            d2 = ( radius*sinRadiansUp, 0.0, 0.0 )
                            v1x = ( mag2*cosRadiansUp*flip, 0.0, mag2*sinRadiansUp )
                            d1x = ( 0.0, 0.0, 0.0 )
                            v2x = ( 0.0, 0.0, 2.0*radiusX/elementsCountUp )
                            d2x = ( 0.0, 0.0, 0.0 )
                    else:
                        xi = (n1 - 1)/elementsCountAcross
                        cxi = xi*2.0
                        flipHalf = (n1 - 1)*2 > elementsCountAcross
                        if flipHalf:
                            cxi = 2.0 - cxi
                        v = interpolateCubicHermite(v1, d1, v2, d2, cxi)
                        d = interpolateCubicHermiteDerivative(v1, d1, v2, d2, cxi)
                        if flipHalf:
                            x = [ -v[0], v[1], v[2] ]
                            dx_ds1 = [ 2.0*d[0]/elementsCountAcross, 0.0, -2.0*d[2]/elementsCountAcross ]
                        else:
                            x = [ v[0], v[1], v[2] ]
                            dx_ds1 = [ 2.0*d[0]/elementsCountAcross, 0.0, 2.0*d[2]/elementsCountAcross ]
                        dx = interpolateCubicHermite(v1x, d1x, v2x, d2x, cxi)
                        if flipHalf:
                            dx_ds2 = [ -dx[0], 0.0, dx[2] ]
                        else:
                            dx_ds2 = [ dx[0], 0.0, dx[2] ]
                        dx_ds3 = [ 0.0, -wallThicknessSeptum, 0.0 ]
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

            # create two top apex nodes

            radius = radiusX
            x[0] = 0.0
            x[1] = -sign*(radiusY + flangeLength)
            x[2] = radius
            dx_ds3[0] = 0.0
            if ((n3 == 1) and (wallThickness[0] > wallThickness[1])) or \
                ((n3 == 0) and (wallThickness[0] < wallThickness[1])):
                dx_ds3[1] = 0.0
                dx_ds3[2] = wallThickness[n3]
            else:
                dx_ds3[1] = sign*wallThickness[n3]*cosBevelAngle
                dx_ds3[2] = wallThickness[n3]*sinBevelAngle
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            mag1 = 2.0*radius/elementsCountAcross
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, -mag1*sinBevelAngle, sign*mag1*cosBevelAngle ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
            nodeIdentifier = nodeIdentifier + 1

            radius = radiusX + wallThickness[n3]
            x[0] = 0.0
            x[1] = -sign*(radiusY + flangeLength)
            #if wallThickness[0] > wallThickness[1]:
            #    x[1] -= sideBulge
            #elif wallThickness[0] < wallThickness[1]:
            #    x[1] += sideBulge
            x[2] = radius
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, -wallThicknessSeptum, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, wallThickness[n3] ])
            nodeIdentifier = nodeIdentifier + 1

        rimCrossAngle = math.pi/4
        sinRimCrossAngle = math.sin(rimCrossAngle)
        cosRimCrossAngle = math.cos(rimCrossAngle)

        scaleFactorsOuter = [
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsOuterApex1 = [ -1.0,
            -cosRimCrossAngle, sinRimCrossAngle, -cosRimCrossAngle, -sinRimCrossAngle,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsOuterApex2 = [ -1.0,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle,
            -cosRimCrossAngle, sinRimCrossAngle, -cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsInner1 = [ -1.0, 4.0,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, sinRimCrossAngle,
            cosRimCrossAngle, -sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle
        ]
        scaleFactorsInner2 = [ -1.0, 4.0,
            cosRimCrossAngle, -sinRimCrossAngle, cosRimCrossAngle, -sinRimCrossAngle,
            cosRimCrossAngle, sinRimCrossAngle, cosRimCrossAngle, sinRimCrossAngle
        ]

        # create elements
        elementIdentifier = 1
        rno = elementsCountAcross + 3
        wno = (elementsCountUp - 1)*rno + 4

        # bottom apex row
        element = mesh.createElement(elementIdentifier, elementtemplateOuterApex1)
        bni11 = 2
        bni12 = 2 + wno
        bni21 = 4
        bni22 = 4 + wno
        nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 - 1, bni12 - 1, bni21 - 1, bni22 - 1 ]
        result = element.setNodesByIdentifier(eftOuterApex1, nodeIdentifiers)
        #print(result, 'bottom outer apex element 1', elementIdentifier, nodeIdentifiers)
        element.setScaleFactors(eftOuterApex1, scaleFactorsOuterApex1)
        elementIdentifier = elementIdentifier + 1

        for e1 in range(0, elementsCountAcross):
            if e1 == 0:
                eftInnerApex = eftInnerApex1a
                elementtemplateInnerApex = elementtemplateInnerApex1a
            elif e1 < (elementsCountAcross - 1):
                eftInnerApex = eftInnerApex1b
                elementtemplateInnerApex = elementtemplateInnerApex1b
            else:
                eftInnerApex = eftInnerApex1c
                elementtemplateInnerApex = elementtemplateInnerApex1c
            va = e1 + 1
            vb = e1 + 2
            for layer in range(2):
                si = layer*(8 if (eftInnerApex is eftInnerApex1b) else 10)
                eftInnerApex.setScaleFactorIdentifier(si + 3, apexScaleFactorIdentifierOffset + va*100 + 1)
                eftInnerApex.setScaleFactorIdentifier(si + 4, apexScaleFactorIdentifierOffset + va*100 + 2)
                eftInnerApex.setScaleFactorIdentifier(si + 5, apexScaleFactorIdentifierOffset + va*100 + 3)
                eftInnerApex.setScaleFactorIdentifier(si + 6, apexScaleFactorIdentifierOffset + vb*100 + 1)
                eftInnerApex.setScaleFactorIdentifier(si + 7, apexScaleFactorIdentifierOffset + vb*100 + 2)
                eftInnerApex.setScaleFactorIdentifier(si + 8, apexScaleFactorIdentifierOffset + vb*100 + 3)
            # redefine field in template for changes to eftInnerApex:
            elementtemplateInnerApex.defineField(coordinates, -1, eftInnerApex)
            element = mesh.createElement(elementIdentifier, elementtemplateInnerApex)
            bni1 = 2
            bni21 = e1 + 4
            bni22 = e1 + 5
            nodeIdentifiers = [ bni1, bni21, bni22, bni1 + wno, bni21 + wno, bni22 + wno ]
            result = element.setNodesByIdentifier(eftInnerApex, nodeIdentifiers)
            #print('eftInnerApex setNodesByIdentifier result ', result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
            # set general linear map coefficients
            radiansAround1 = math.pi + e1*radiansPerElementAcross
            radiansAround1Next = math.pi + (e1 + 1)*radiansPerElementAcross
            radiansAround2 = math.pi - e1*radiansPerElementAcross
            radiansAround2Next = math.pi - (e1 + 1)*radiansPerElementAcross
            scalefactors = [
                -1.0, 4.0,
                math.sin(radiansAround1), math.cos(radiansAround1), radiansPerElementAcross,
                math.sin(radiansAround1Next), math.cos(radiansAround1Next), radiansPerElementAcross,
                -cosRimCrossAngle, -sinRimCrossAngle,
                math.sin(radiansAround2), math.cos(radiansAround2), -radiansPerElementAcross,
                math.sin(radiansAround2Next), math.cos(radiansAround2Next), -radiansPerElementAcross,
                -cosRimCrossAngle, -sinRimCrossAngle
            ]
            if eftInnerApex is not eftInnerApex1b:
                scalefactors = scalefactors[:10] + [ cosRimCrossAngle, sinRimCrossAngle ] + scalefactors[10:] + [cosRimCrossAngle, -sinRimCrossAngle ]
            result = element.setScaleFactors(eftInnerApex, scalefactors)
            #print('eftInnerApex setScaleFactors', result, scalefactors)
            elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplateOuterApex0)
        bni11 = wno + 2
        bni12 = 2
        bni21 = wno + rno + 1
        bni22 = rno + 1
        nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 - 1, bni12 - 1, bni21 + 1, bni22 + 1 ]
        result = element.setNodesByIdentifier(eftOuterApex0, nodeIdentifiers)
        #print(result, 'bottom outer apex element 2', elementIdentifier, nodeIdentifiers)
        element.setScaleFactors(eftOuterApex0, scaleFactorsOuter)
        elementIdentifier = elementIdentifier + 1

        # 'regular' rows between apexes
        for e2 in range(0, elementsCountUp - 2):
            bn = e2*rno + 2
            element = mesh.createElement(elementIdentifier, elementtemplateOuter)
            bni11 = bn + 2
            bni12 = bn + wno + 2
            bni21 = bn + rno + 2
            bni22 = bn + wno + rno + 2
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 - 1, bni12 - 1, bni21 - 1, bni22 - 1 ]
            result = element.setNodesByIdentifier(eftOuter, nodeIdentifiers)
            #print(result, 'outer1 element', elementIdentifier, nodeIdentifiers)
            element.setScaleFactors(eftOuter, scaleFactorsOuter)
            elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateInner1)
            bni11 = bn + 2
            bni12 = bn + 3
            bni21 = bn + rno + 2
            bni22 = bn + rno + 3
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
            result = element.setNodesByIdentifier(eftInner1, nodeIdentifiers)
            #print(result, 'inner1 element', elementIdentifier, 'nodes', nodeIdentifiers)
            result = element.setScaleFactors(eftInner1, scaleFactorsInner1)
            #print(result, 'inner1 element', elementIdentifier, 'scale factors', scaleFactorsInner1)
            elementIdentifier = elementIdentifier + 1

            for e1 in range(elementsCountAcross - 2):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = bn + e1 + 3
                bni12 = bn + e1 + 4
                bni21 = bn + rno + e1 + 3
                bni22 = bn + rno + e1 + 4
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                #print(result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateInner2)
            bni11 = bn + rno - 2
            bni12 = bn + rno - 1
            bni21 = bn + 2*rno -2
            bni22 = bn + 2*rno - 1
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + wno, bni12 + wno, bni21 + wno, bni22 + wno ]
            result = element.setNodesByIdentifier(eftInner2, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
            #result = element.setScaleFactor(eftInner2, 1, -1.0)
            result = element.setScaleFactors(eftInner2, scaleFactorsInner2)
            #print(result, 'inner2 element', elementIdentifier, 'scale factors', scaleFactorsInner1)
            #print(result, 'element', elementIdentifier, 'scale factors', scaleFactorsInner1)
            elementIdentifier = elementIdentifier + 1

            element = mesh.createElement(elementIdentifier, elementtemplateOuter)
            bni11 = bn + wno + rno - 1
            bni12 = bn + rno - 1
            bni21 = bn + wno + 2*rno - 1
            bni22 = bn + 2*rno - 1
            nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + 1, bni12 + 1, bni21 + 1, bni22 + 1 ]
            result = element.setNodesByIdentifier(eftOuter, nodeIdentifiers)
            #print(result, 'element', elementIdentifier, nodeIdentifiers)
            element.setScaleFactors(eftOuter, scaleFactorsOuter)
            elementIdentifier = elementIdentifier + 1

        # top apex row
        element = mesh.createElement(elementIdentifier, elementtemplateOuterApex0)
        bni11 = wno - rno
        bni12 = 2*wno - rno
        bni21 = wno - 1
        bni22 = 2*wno - 1
        nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 - 1, bni12 - 1, bni21 + 1, bni22 + 1 ]
        result = element.setNodesByIdentifier(eftOuterApex0, nodeIdentifiers)
        #print(result, 'top outer apex element 1', elementIdentifier, nodeIdentifiers)
        element.setScaleFactors(eftOuterApex0, scaleFactorsOuter)
        elementIdentifier = elementIdentifier + 1

        for e1 in range(0, elementsCountAcross):
            if e1 == 0:
                eftInnerApex = eftInnerApex2a
                elementtemplateInnerApex = elementtemplateInnerApex2a
            elif e1 < (elementsCountAcross - 1):
                eftInnerApex = eftInnerApex2b
                elementtemplateInnerApex = elementtemplateInnerApex2b
            else:
                eftInnerApex = eftInnerApex2c
                elementtemplateInnerApex = elementtemplateInnerApex2c
            va = e1 + 1
            vb = e1 + 2
            for layer in range(2):
                si = layer*(8 if (eftInnerApex is eftInnerApex2b) else 10)
                eftInnerApex.setScaleFactorIdentifier(si + 3, apexScaleFactorIdentifierOffset + va*100 + 1)
                eftInnerApex.setScaleFactorIdentifier(si + 4, apexScaleFactorIdentifierOffset + va*100 + 2)
                eftInnerApex.setScaleFactorIdentifier(si + 5, apexScaleFactorIdentifierOffset + va*100 + 3)
                eftInnerApex.setScaleFactorIdentifier(si + 6, apexScaleFactorIdentifierOffset + vb*100 + 1)
                eftInnerApex.setScaleFactorIdentifier(si + 7, apexScaleFactorIdentifierOffset + vb*100 + 2)
                eftInnerApex.setScaleFactorIdentifier(si + 8, apexScaleFactorIdentifierOffset + vb*100 + 3)
            # redefine field in template for changes to eftInnerApex:
            elementtemplateInnerApex.defineField(coordinates, -1, eftInnerApex)
            element = mesh.createElement(elementIdentifier, elementtemplateInnerApex)
            bni11 = wno - rno + e1
            bni12 = bni11 + 1
            bni2 = wno - 1
            nodeIdentifiers = [ bni11, bni12, bni2, bni11 + wno, bni12 + wno, bni2 + wno ]
            result = element.setNodesByIdentifier(eftInnerApex, nodeIdentifiers)
            #print('eftInnerApex setNodesByIdentifier result ', result, 'element', elementIdentifier, 'nodes', nodeIdentifiers)
            # set general linear map coefficients
            radiansAround1 = math.pi + e1*radiansPerElementAcross
            radiansAround1Next = math.pi + (e1 + 1)*radiansPerElementAcross
            radiansAround2 = math.pi - e1*radiansPerElementAcross
            radiansAround2Next = math.pi - (e1 + 1)*radiansPerElementAcross
            scalefactors = [
                -1.0, 4.0, 
                math.sin(radiansAround1), -math.cos(radiansAround1), radiansPerElementAcross,
                math.sin(radiansAround1Next), -math.cos(radiansAround1Next), radiansPerElementAcross,
                -cosRimCrossAngle, -sinRimCrossAngle,
                math.sin(radiansAround2), -math.cos(radiansAround2), -radiansPerElementAcross,
                math.sin(radiansAround2Next), -math.cos(radiansAround2Next), -radiansPerElementAcross,
                -cosRimCrossAngle, -sinRimCrossAngle
            ]
            if eftInnerApex is not eftInnerApex2b:
                scalefactors = scalefactors[:10] + [ cosRimCrossAngle, sinRimCrossAngle ] + scalefactors[10:] + [cosRimCrossAngle, -sinRimCrossAngle ]
            result = element.setScaleFactors(eftInnerApex, scalefactors)
            #print('eftInnerApex setScaleFactors', result, scalefactors)
            elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplateOuterApex2)
        bni11 = 2*wno - 3
        bni12 = wno - 3
        bni21 = 2*wno - 1
        bni22 = wno - 1
        nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + 1, bni12 + 1, bni21 + 1, bni22 + 1 ]
        result = element.setNodesByIdentifier(eftOuterApex2, nodeIdentifiers)
        #print(result, 'top outer apex element 2', elementIdentifier, nodeIdentifiers)
        result = element.setScaleFactors(eftOuterApex2, scaleFactorsOuterApex2)
        #print(result, 'top outer apex element 2 scale', elementIdentifier, scaleFactorsOuterApex2)
        elementIdentifier = elementIdentifier + 1


        if bulgeRadius != 0.0:
            # spherical polar coordinates:
            # r = y - bulgeRadius
            # theta = -x / bulgeRadius
            # phi = z / bulgeRadius
            yxzCoordinates = fm.createFieldComponent(coordinates, [2, 1, 3])
            scale = fm.createFieldConstant([1.0, -1.0/bulgeRadius, -1.0/bulgeRadius ])
            scaleCoordinates = fm.createFieldMultiply(yxzCoordinates, scale)
            offset = fm.createFieldConstant([-bulgeRadius, 0.0, 0.0 ])
            polarCoordinates = fm.createFieldAdd(scaleCoordinates, offset)
            polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_SPHERICAL_POLAR)
            rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
            rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
            newyxzCoordinates = fm.createFieldSubtract(rcCoordinates, offset)
            newCoordinates = fm.createFieldComponent(newyxzCoordinates, [2, 1, 3])
            fieldassignment = coordinates.createFieldassignment(newCoordinates)
            result = fieldassignment.assign()

        fm.endChange()
        return []
