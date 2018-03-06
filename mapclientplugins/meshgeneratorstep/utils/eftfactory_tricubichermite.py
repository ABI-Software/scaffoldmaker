'''
Definitions of standard element field templates shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''
from mapclientplugins.meshgeneratorstep.utils.eft_utils import *
from opencmiss.zinc.element import Elementbasis, Elementfieldtemplate
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK

class eftfactory_tricubichermite:
    '''
    Factory class for creating element field templates for a 3-D mesh using tricubic Hermite basis.
    '''

    def __init__(self, mesh, useCrossDerivatives):
        '''
        :param mesh:  Zinc mesh to create element field templates in.
        :param useCrossDerivatives: Set to True if you want cross derivative terms.
        '''
        assert mesh.getDimension() == 3, 'eftfactory_tricubichermite: not a 3-D Zinc mesh'
        self._mesh = mesh
        self._useCrossDerivatives = useCrossDerivatives
        self._fieldmodule = mesh.getFieldmodule()
        self._tricubicHermiteBasis = self._fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

    def createEftBasic(self):
        '''
        Create the basic tricubic hermite element template with 1:1 mappings to
        node derivatives, with or without cross derivatives.
        :return: Element field template
        '''
        if not self._useCrossDerivatives:
            return self.createEftNoCrossDerivatives()
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftBasic:  Failed to validate eft'
        return eft

    def createEftNoCrossDerivatives(self):
        '''
        Create a basic tricubic hermite element template with 1:1 mappings to
        node derivatives, without cross derivatives.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        for n in range(8):
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftNoCrossDerivatives:  Failed to validate eft'
        return eft

    def createEftSplitXi1LeftStraight(self):
        '''
        Create an element field template suitable for the inner elements of the
        join between left and right chambers, with xi1 bifurcating.
        :return: Element field template
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eft, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftSplitXi1LeftStraight:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightStraight(self):
        '''
        Create an element field template suitable for the inner elements of the
        join between left and right chambers, with xi1 bifurcating.
        Straight through version.
        :return: Element field template
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 5, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eft, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, []) ])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftSplitXi1RightStraight:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightIn(self):
        '''
        Create an element field template suitable for the outer elements of the
        join between left and right chambers, with xi1 merging from the right.
        Right in version i.e. xi1 heading in from right. h-shape.
        :return: Element field template
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
        remapEftNodeValueLabel(eft, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eft, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftSplitXi1RightOut:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightOut(self):
        '''
        Create an element field template suitable for the outer elements of the
        join between left and right chambers, with xi1 bifurcating.
        Right out version i.e. xi1 heading to right. h-shape.
        :return: Element field template
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
        remapEftNodeValueLabel(eft, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eft, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, []) ])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftSplitXi1RightOut:  Failed to validate eft'
        return eft

    def createEftTubeSeptumOuter(self):
        '''
        Create an element field template suitable for the outer elements of
        a tube septum e.g. top and bottom of (|), general linear mapping inner
        derivatives to give reduced continuity across the septum, but with
        continuity around the two openings.
        Cross derivatives are not used on the general mapped nodes.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        # general linear map at 4 nodes for one derivative
        eft.setNumberOfLocalScaleFactors(8)
        for s in range(8):
            eft.setScaleFactorType(s + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eft.setScaleFactorIdentifier(s + 1, (s % 2) + 1)
        if self._useCrossDerivatives:
            noCrossRange = range(4)
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi1 from ds1 + ds3 for first 4 nodes
        for n in range(4):
            ln = n + 1
            eft.setFunctionNumberOfTerms(n*8 + 2, 2)
            eft.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eft.setTermScaling(n*8 + 2, 1, [n*2 + 1])
            eft.setTermNodeParameter(n*8 + 2, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 2, 2, [n*2 + 2])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftTubeSeptumOuter:  Failed to validate eft'
        return eft

    def createEftTubeSeptumInner1(self):
        '''
        Create an element field template suitable for the inner bottom elements
        of a tube septum e.g. inside bottom of (|), general linear mapping inner
        derivatives to give reduced continuity across the septum, but with
        continuity around the two openings.
        Cross derivatives are not used on the general mapped nodes.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        # negate dxi1 plus general linear map at 4 nodes for one derivative
        eft.setNumberOfLocalScaleFactors(10)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eft.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eft.setScaleFactorIdentifier(1, 1)
        # Global scale factor 4.0 used for cross derivative term to correct first derivative
        eft.setScaleFactorType(2, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eft.setScaleFactorIdentifier(2, 2)
        for s in range(8):
            eft.setScaleFactorType(s + 3, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eft.setScaleFactorIdentifier(s + 3, (s % 2) + 1)
        if self._useCrossDerivatives:
            noCrossRange = [ 0, 2, 4, 6 ]
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi3 from ds1 + ds3 for 4 odd nodes
        s = 0
        for n in [ 0, 2, 4, 6 ]:
            ln = n + 1
            # 2 terms for d/dx3 via general linear map
            eft.setFunctionNumberOfTerms(n*8 + 5, 2)
            eft.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eft.setTermScaling(n*8 + 5, 1, [s*2 + 3])
            eft.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 5, 2, [s*2 + 4])
            # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
            # GRC WIP
            eft.setFunctionNumberOfTerms(n*8 + 6, 1)
            eft.setTermNodeParameter(n*8 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 6, 1, [1, 2] if (n < 4) else [2])
            s += 1
        # negate d/dxi1 at 2 nodes
        for n in [4, 6]:
            result = eft.setTermScaling(n*8 + 2, 1, [1])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftTubeSeptumInner1:  Failed to validate eft'
        return eft

    def createEftTubeSeptumInner2(self):
        '''
        Create an element field template suitable for the inner top elements
        of a tube septum e.g. inside top of (|), general linear mapping inner
        derivatives to give reduced continuity across the septum, but with
        continuity around the two openings.
        Cross derivatives are not used on the general mapped nodes.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        # negate dxi1 plus general linear map at 4 nodes for one derivative
        eft.setNumberOfLocalScaleFactors(10)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eft.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eft.setScaleFactorIdentifier(1, 1)
        # Global scale factor 4.0 used for cross derivative term to correct first derivative
        eft.setScaleFactorType(2, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eft.setScaleFactorIdentifier(2, 2)
        for s in range(8):
            eft.setScaleFactorType(s + 3, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eft.setScaleFactorIdentifier(s + 3, (s % 2) + 1)
        if self._useCrossDerivatives:
            noCrossRange = [ 1, 3, 5, 7 ]
        else:
            noCrossRange = range(8)
        for n in noCrossRange:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        # general linear map dxi3 from ds1 + ds3 for 4 even nodes
        s = 0
        for n in [ 1, 3, 5, 7 ]:
            ln = n + 1
            # 2 terms for d/dx3 via general linear map
            eft.setFunctionNumberOfTerms(n*8 + 5, 2)
            eft.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eft.setTermScaling(n*8 + 5, 1, [1, s*2 + 3])
            eft.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 5, 2, [1, s*2 + 4])
            # add d2/dxi1dxi3 correction along xi1 == 0 to fit septum outer xi1 derivative better
            # GRC WIP
            eft.setFunctionNumberOfTerms(n*8 + 6, 1)
            eft.setTermNodeParameter(n*8 + 6, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            eft.setTermScaling(n*8 + 6, 1, [2] if (n < 4) else [1, 2])
            s += 1
        # negate d/dxi1 at 2 nodes
        for n in [5, 7]:
            eft.setTermScaling(n*8 + 2, 1, [1])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftTubeSeptumInner2:  Failed to validate eft'
        return eft

    def setEftLinearDerivativeXi1(self, eft, basisNode1, basisNode2, localNode1, localNode2, minus1scaleFactorIndex):
        '''
        Makes the derivative in Xi1 for basisNode1 and basisNode2 equal to the difference in
        the value parameter at localNode2 - localNode1. This makes the basis equivalent to linear Lagrange.
        Note! Cross derivatives are not handled and are currently unmodified.
        :param minus1scaleFactorIndex: Local scale factor index for general value -1.0
        '''
        for n in [basisNode1 - 1, basisNode2 - 1]:
            f = n*8 + 2
            eft.setFunctionNumberOfTerms(f, 2)
            eft.setTermNodeParameter(f, 1, localNode2, Node.VALUE_LABEL_VALUE, 1)
            eft.setTermScaling(f, 1, [])
            eft.setTermNodeParameter(f, 2, localNode1, Node.VALUE_LABEL_VALUE, 1)
            eft.setTermScaling(f, 2, [minus1scaleFactorIndex])

    def setEftLinearDerivativeXi3(self, eft, basisNode1, basisNode2, localNode1, localNode2, minus1scaleFactorIndex):
        '''
        Makes the derivative in Xi3 for basisNode1 and basisNode2 equal to the difference in
        the value parameter at localNode2 - localNode1. This makes the basis equivalent to linear Lagrange.
        Note! Cross derivatives are not handled and are currently unmodified.
        :param minus1scaleFactorIndex: Local scale factor index for general value -1.0
        '''
        for n in [basisNode1 - 1, basisNode2 - 1]:
            f = n*8 + 5
            eft.setFunctionNumberOfTerms(f, 2)
            eft.setTermNodeParameter(f, 1, localNode2, Node.VALUE_LABEL_VALUE, 1)
            eft.setTermScaling(f, 1, [])
            eft.setTermNodeParameter(f, 2, localNode1, Node.VALUE_LABEL_VALUE, 1)
            eft.setTermScaling(f, 2, [minus1scaleFactorIndex])

    def setEftMidsideXi1HangingNode(self, eft, hangingBasisNode, otherBasisNode, localNode1, localNode2, scaleFactorIndexes):
        '''
        Makes the functions for eft at basisNode work as a hanging node interpolating the parameters
        from localNode1 at xi1=0 and localNode2 at xi1=1 to midside xi1 = 0.5.
        Note! Cross derivatives are not handled and are currently unmodified.
        :param otherBasisNode: Other node along xi1 which needs its dxi1 derivative halved
        :param scaleFactorIndexes: Local scale factor indexes for general values -1.0 0.5 0.25 0.125 0.75
        '''
        n = hangingBasisNode - 1
        o = otherBasisNode - 1
        sfneg1 = scaleFactorIndexes[0]
        sf05 = scaleFactorIndexes[1]
        sf025 = scaleFactorIndexes[2]
        sf0125 = scaleFactorIndexes[3]
        sf075 = scaleFactorIndexes[4]
        # otherBasisNode d/dxi1 must be halved
        eft.setTermScaling(o*8 + 2, 1, [sf05])
        # workaround for Zinc limitation where faces are not found due to only first
        # node in general linear map being reported; reversing order of terms fixes this
        otherLocalNode = eft.getTermLocalNodeIndex(o*8 + 2, 1)
        termOrder = [ 3, 4, 1, 2] if (otherLocalNode == localNode1) else [ 1, 2, 3, 4]
        # value = 0.5*x_1 + 0.125*ds1_1 + 0.5*x_2 - 0.125*ds1_2
        eft.setFunctionNumberOfTerms(n*8 + 1, 4)
        eft.setTermNodeParameter(n*8 + 1, termOrder[0], localNode1, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 1, termOrder[0], [sf05])
        eft.setTermNodeParameter(n*8 + 1, termOrder[1], localNode1, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 1, termOrder[1], [sf0125])
        eft.setTermNodeParameter(n*8 + 1, termOrder[2], localNode2, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 1, termOrder[2], [sf05])
        eft.setTermNodeParameter(n*8 + 1, termOrder[3], localNode2, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 1, termOrder[3], [sfneg1, sf0125])
        # d/dxi1 = -0.75*x_1 - 0.125*ds1_1 + 0.75*x_2 - 0.125*ds1_2
        eft.setFunctionNumberOfTerms(n*8 + 2, 4)
        eft.setTermNodeParameter(n*8 + 2, 1, localNode1, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 2, 1, [sfneg1, sf075])
        eft.setTermNodeParameter(n*8 + 2, 2, localNode1, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 2, 2, [sfneg1, sf0125])
        eft.setTermNodeParameter(n*8 + 2, 3, localNode2, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 2, 3, [sf075])
        eft.setTermNodeParameter(n*8 + 2, 4, localNode2, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 2, 4, [sfneg1, sf0125])
        # d/dxi2 = 0.5*ds2_1 + 0.5*ds2_2
        eft.setFunctionNumberOfTerms(n*8 + 3, 2)
        eft.setTermNodeParameter(n*8 + 3, 1, localNode1, Node.VALUE_LABEL_D_DS2, 1)
        eft.setTermScaling(n*8 + 3, 1, [sf05])
        eft.setTermNodeParameter(n*8 + 3, 2, localNode2, Node.VALUE_LABEL_D_DS2, 1)
        eft.setTermScaling(n*8 + 3, 2, [sf05])
        # d/dxi3 = 0.5*ds3_1 + 0.5*ds3_2
        eft.setFunctionNumberOfTerms(n*8 + 5, 2)
        eft.setTermNodeParameter(n*8 + 5, 1, localNode1, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 5, 1, [sf05])
        eft.setTermNodeParameter(n*8 + 5, 2, localNode2, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 5, 2, [sf05])


    def setEftMidsideXi3HangingNode(self, eft, hangingBasisNode, otherBasisNode, localNode1, localNode2, scaleFactorIndexes):
        '''
        Makes the functions for eft at basisNode work as a hanging node interpolating the parameters
        from localNode1 at xi3=0 and localNode2 at xi3=1 to midside xi3 = 0.5.
        Note! Cross derivatives are not handled and are currently unmodified.
        :param otherBasisNode: Other node along xi3 which needs its dxi3 derivative halved
        :param scaleFactorIndexes: Local scale factor indexes for general values -1.0 0.5 0.25 0.125 0.75
        '''
        n = hangingBasisNode - 1
        o = otherBasisNode - 1
        sfneg1 = scaleFactorIndexes[0]
        sf05 = scaleFactorIndexes[1]
        sf025 = scaleFactorIndexes[2]
        sf0125 = scaleFactorIndexes[3]
        sf075 = scaleFactorIndexes[4]
        # otherBasisNode d/dxi3 must be halved
        eft.setTermScaling(o*8 + 5, 1, [sf05])
        # workaround for Zinc limitation where faces are not found due to only first
        # node in general linear map being reported; reversing order of terms fixes this
        otherLocalNode = eft.getTermLocalNodeIndex(o*8 + 5, 1)
        termOrder = [ 3, 4, 1, 2] if (otherLocalNode == localNode1) else [ 1, 2, 3, 4]
        # value = 0.5*x_1 + 0.125*ds3_1 + 0.5*x_2 - 0.125*ds3_2
        eft.setFunctionNumberOfTerms(n*8 + 1, 4)
        eft.setTermNodeParameter(n*8 + 1, termOrder[0], localNode1, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 1, termOrder[0], [sf05])
        eft.setTermNodeParameter(n*8 + 1, termOrder[1], localNode1, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 1, termOrder[1], [sf0125])
        eft.setTermNodeParameter(n*8 + 1, termOrder[2], localNode2, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 1, termOrder[2], [sf05])
        eft.setTermNodeParameter(n*8 + 1, termOrder[3], localNode2, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 1, termOrder[3], [sfneg1, sf0125])
        # d/dxi1 = 0.5*ds1_1 + 0.5*ds1_2
        eft.setFunctionNumberOfTerms(n*8 + 2, 2)
        eft.setTermNodeParameter(n*8 + 2, 1, localNode1, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 2, 1, [sf05])
        eft.setTermNodeParameter(n*8 + 2, 2, localNode2, Node.VALUE_LABEL_D_DS1, 1)
        eft.setTermScaling(n*8 + 2, 2, [sf05])
        # d/dxi2 = 0.5*ds2_1 + 0.5*ds2_2
        eft.setFunctionNumberOfTerms(n*8 + 3, 2)
        eft.setTermNodeParameter(n*8 + 3, 1, localNode1, Node.VALUE_LABEL_D_DS2, 1)
        eft.setTermScaling(n*8 + 3, 1, [sf05])
        eft.setTermNodeParameter(n*8 + 3, 2, localNode2, Node.VALUE_LABEL_D_DS2, 1)
        eft.setTermScaling(n*8 + 3, 2, [sf05])
        # d/dxi3 = -0.75*x_1 - 0.125*ds3_1 + 0.75*x_2 - 0.125*ds3_2
        eft.setFunctionNumberOfTerms(n*8 + 5, 4)
        eft.setTermNodeParameter(n*8 + 5, 1, localNode1, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 5, 1, [sfneg1, sf075])
        eft.setTermNodeParameter(n*8 + 5, 2, localNode1, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 5, 2, [sfneg1, sf0125])
        eft.setTermNodeParameter(n*8 + 5, 3, localNode2, Node.VALUE_LABEL_VALUE, 1)
        eft.setTermScaling(n*8 + 5, 3, [sf075])
        eft.setTermNodeParameter(n*8 + 5, 4, localNode2, Node.VALUE_LABEL_D_DS3, 1)
        eft.setTermScaling(n*8 + 5, 4, [sfneg1, sf0125])

