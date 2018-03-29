'''
Definitions of standard element field templates shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK
import math

def normalise(v):
    '''
    :return: vector v normalised to unit length
    '''
    mag = 0.0
    for s in v:
        mag += s*s
    mag = math.sqrt(mag)
    return [ s/mag for s in v ]

def crossproduct3(a, b):
    '''
    :return: vector 3-D cross product of a and b
    '''
    return [ a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] ]

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
        Create the basic tricubic hermite element field template with 1:1 mappings to
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
        Create a basic tricubic hermite element field template with 1:1 mappings to
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

    def createEftShellApexBottom(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a tricubic hermite element field for closing bottom apex of a shell.
        Element is collapsed in xi1 on xi2 = 0.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients
        of the radial line from global derivatives, plus the arc subtended by
        the element in radians, so the apex can be rounded.
        Need to create a new template for each sector around apex giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around apex.
        :param nodeScaleFactorOffset0: offset of node scale factors at apex on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at apex on xi1=1
        :return: Element field template
        '''
        # start with full tricubic to remap D2_DS1DS2 at apex
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        if not self._useCrossDerivatives:
            for n in [ 2, 3, 6, 7 ]:
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ])
        # remap parameters before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
        for layer in range(2):
            so = layer*6 + 1
            ln = layer*4 + 1
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 1]), (Node.VALUE_LABEL_D_DS2, [so + 2]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 2, so + 3]), (Node.VALUE_LABEL_D_DS2, [1, so + 1, so + 3]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 2
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, so + 4), (Node.VALUE_LABEL_D_DS2, so + 5) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 5, so + 6]), (Node.VALUE_LABEL_D_DS2, [1, so + 4, so + 6]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        ln_map = [ 1, 1, 2, 3, 4, 4, 5, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftShellApexBottom:  Failed to validate eft'
        return eft

    def createEftShellApexTop(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a tricubic hermite element field for closing top apex of a shell.
        Element is collapsed in xi1 on xi2 = 1.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients
        of the radial line from global derivatives, plus the arc subtended by
        the element in radians, so the apex can be rounded.
        Need to create a new template for each sector around apex giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around apex.
        :param nodeScaleFactorOffset0: offset of node scale factors at apex on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at apex on xi1=1
        :return: Element field template
        '''
        # start with full tricubic to remap D2_DS1DS2 at apex
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        if not self._useCrossDerivatives:
            for n in [ 0, 1, 4, 5 ]:
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3 ])
        # remap parameters before collapsing nodes
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
        for layer in range(2):
            so = layer*6 + 1
            ln = layer*4 + 3
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [so + 1]), (Node.VALUE_LABEL_D_DS2, [so + 2]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, so + 2, so + 3]), (Node.VALUE_LABEL_D_DS2, [so + 1, so + 3]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 4
            # 2 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, so + 4), (Node.VALUE_LABEL_D_DS2, so + 5) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, so + 5, so + 6]), (Node.VALUE_LABEL_D_DS2, [so + 4, so + 6]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        ln_map = [ 1, 2, 3, 3, 4, 5, 6, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftShellApexTop:  Failed to validate eft'
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

    def replaceElementWithInlet4(self, element, startElementId, nodetemplate, startNodeId, tubeLength, innerDiameter, wallThickness):
        '''
        Replace element with 4 element X-layout tube inlet.
        Inlet axis is at given length from centre of xi3=0 face, oriented with dx/dxi1.
        8 new nodes are created.
        '''
        fm = self._mesh.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        fm.beginChange()
        cache = fm.createFieldcache()
        diff1 = self._mesh.getChartDifferentialoperator(1, 1)
        diff2 = self._mesh.getChartDifferentialoperator(1, 2)
        coordinates = getOrCreateCoordinateField(fm)
        cache.setMeshLocation(element, [0.5, 0.5, 1.0])
        result, fc = coordinates.evaluateReal(cache, 3)
        resulta, a = coordinates.evaluateDerivative(diff1, cache, 3)
        resultb, b = coordinates.evaluateDerivative(diff2, cache, 3)
        n = normalise(crossproduct3(a, b))
        #print(resulta, 'a =', a, ',', resultb, ' b=', b, ' fc=', fc, ' n=',n)
        ic = [ (fc[i] + tubeLength*n[i]) for i in range(3) ]
        na = normalise(a)
        nb = normalise(b)
        a = normalise([ -(na[i] + nb[i]) for i in range(3) ])
        b = normalise(crossproduct3(a, n))

        zero = [ 0.0, 0.0, 0.0 ]
        nodeIdentifier = startNodeId
        elementsCountAround = 4
        radiansPerElementAround = math.pi*2.0/elementsCountAround
        for n3 in range(2):
            radius = innerDiameter*0.5 + n3*wallThickness
            for n1 in range(elementsCountAround):
                radiansAround = n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x = [ (ic[i] + radius*(cosRadiansAround*a[i] + sinRadiansAround*b[i])) for i in range(3) ]
                dx_ds1 = [ radiansPerElementAround*radius*(-sinRadiansAround*a[i] + cosRadiansAround*b[i]) for i in range(3) ]
                dx_ds2 = [ -tubeLength*c for c in n ]
                dx_ds3 = [ wallThickness*(cosRadiansAround*a[i] + sinRadiansAround*b[i]) for i in range(3) ]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                if self._useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier = nodeIdentifier + 1

        eft0 = element.getElementfieldtemplate(coordinates, -1)
        nids0 = getElementNodeIdentifiers(element, eft0)
        orig_nids = [ nids0[0], nids0[2], nids0[3], nids0[1], nids0[4], nids0[6], nids0[7], nids0[5] ]

        elementIdentifier = startElementId
        elementtemplate1 = self._mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        for e in range(4):
            eft1 = self.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            if e == 0:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, []) ])
            elif e == 1:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
            elif e == 2:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, [1]) ])
            elif e == 3:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1]) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, [1]) ])
            ea = e
            eb = (e + 1) % 4
            ec = ea + 4
            ed = eb + 4
            nids = [
                startNodeId + ea, startNodeId + eb, orig_nids[ea], orig_nids[eb],
                startNodeId + ec, startNodeId + ed, orig_nids[ec], orig_nids[ed]
            ]
            elementtemplate1.defineField(coordinates, -1, eft1)
            element = self._mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 7
            #print('create element in', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1

        fm.endChange()
