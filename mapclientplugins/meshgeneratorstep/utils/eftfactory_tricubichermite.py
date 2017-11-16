'''
Definitions of standard element field templates shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''
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
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        if not self._useCrossDerivatives:
            for n in range(8):
                eft.setFunctionNumberOfTerms(n*8 + 4, 0)
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
                eft.setFunctionNumberOfTerms(n*8 + 7, 0)
                eft.setFunctionNumberOfTerms(n*8 + 8, 0)
        assert eft.validate(), 'eftfactory_tricubichermite.createEftBasic:  Failed to validate eft'
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

