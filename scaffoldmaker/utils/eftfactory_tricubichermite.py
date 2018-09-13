'''
Definitions of standard element field templates shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.zinc_utils import *
import scaffoldmaker.utils.vector as vector
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK
import math

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

    def createEftShellPole90(self, quadrant):
        '''
        Create a 6-node wedge element for around a pole with 90 degrees between sides.
        Xi1 is around, xi2 is toward pole, xi3 is out of surface.
        :param quadrant: quadrant from 0 to 3 from +s1 direction around +s2 in first quadrant
        Element has two global scale factors to set: 1 = -1.0, 90 = math.pi/2.0
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [ 1, 90 ], [])  # global scale factor 90 = pi/2
        remapEftNodeValueLabel(eft, [ 3, 7, 4, 8 ], Node.VALUE_LABEL_D_DS1, [])
        if quadrant == 0:
            remapEftNodeValueLabel(eft, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
            scaleEftNodeValueLabels(eft, [ 4, 8 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            crossFix37 = ( Node.VALUE_LABEL_D_DS2, [ 1, 2 ])
            crossFix48 = ( Node.VALUE_LABEL_D_DS1, [ 2 ])
        elif quadrant == 1:
            scaleEftNodeValueLabels(eft, [ 3, 7 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
            remapEftNodeValueLabel(eft, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
            crossFix37 = ( Node.VALUE_LABEL_D_DS1, [ 2 ])
            crossFix48 = ( Node.VALUE_LABEL_D_DS2, [ 2 ])
        elif quadrant == 2:
            remapEftNodeValueLabel(eft, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            crossFix37 = ( Node.VALUE_LABEL_D_DS2, [ 2 ])
            crossFix48 = ( Node.VALUE_LABEL_D_DS1, [ 1, 2 ])
        elif quadrant == 3:
            remapEftNodeValueLabel(eft, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
            crossFix37 = ( Node.VALUE_LABEL_D_DS1, [ 1, 2 ])
            crossFix48 = ( Node.VALUE_LABEL_D_DS2, [ 1, 2 ])
        else:
            assert False, 'eftfactory_tricubichermite.createEftShellPole90:  Invalid quadrant'

        # 2 terms for cross derivative 1 2 to correct circular apex
        for ln in [ 3, 7 ]:
            mapEftFunction1Node1Term(eft, (ln - 1)*8 + 4, ln, crossFix37[0], 1, crossFix37[1])
        for ln in [ 4, 8 ]:
            mapEftFunction1Node1Term(eft, (ln - 1)*8 + 4, ln, crossFix48[0], 1, crossFix48[1])

        ln_map = [ 1, 2, 3, 3, 4, 5, 6, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftShellPole90:  Failed to validate eft'
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

    def createEftWedgeRadial(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1):
        '''
        Create a tricubic hermite element field for the central axis of a solid
        cylinder, with xi1 collapsed, xi2 up and xi3 out radially.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients
        of the radial line from global derivatives, plus the arc subtended by
        the element in radians, so the circumferential direction is rounded.
        Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :return: Element field template
        '''
        # start with full tricubic
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        if not self._useCrossDerivatives:
            for n in [ 4, 5, 6, 7 ]:
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
        remapEftNodeValueLabel(eft, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS1, [])
        for layer in range(2):
            so = layer*6 + 1
            ln = layer*2 + 1
            # 2 terms for d/dxi3 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [so + 1]), (Node.VALUE_LABEL_D_DS3, [so + 2]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, so + 2, so + 3]), (Node.VALUE_LABEL_D_DS3, [so + 1, so + 3]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*2 + 2
            # 2 terms for d/dxi3 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, so + 4), (Node.VALUE_LABEL_D_DS3, so + 5) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, so + 5, so + 6]), (Node.VALUE_LABEL_D_DS3, [so + 4, so + 6]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        ln_map = [ 1, 1, 2, 2, 3, 4, 5, 6 ]
        remapEftLocalNodes(eft, 6, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftWedgeRadial:  Failed to validate eft'
        return eft

    def createEftTetrahedronBottom(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1, nodeScaleFactorOffsetUp):
        '''
        Create a tricubic hermite element field for a solid tetrahedron for the axis of a
        solid sphere pole, with xi1 and xi3 collapsed on xi2 = 0, and xi1 collapsed on xi3 = 0.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients of
        the radial line from global derivatives, plus the arc subtended by the element in radians,
        so the circumferential direction is rounded. Collapsed nodes on xi2 = 0 have 2 additional
        scale factors cos and sin coefficients of the inclination angle.
        Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :param nodeScaleFactorOffsetUp: offset of first scale factor for inclination at pole,
        increase by 2 in each layer away from axis. Suggest starting at 100000 on axis.
        :return: Element field template
        '''
        # start with full tricubic
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        for n in [ 2, 3, 6, 7 ]:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            if n > 3:
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffsetUp + 1, nodeScaleFactorOffsetUp + 2,
            nodeScaleFactorOffsetUp + 3, nodeScaleFactorOffsetUp + 4])

        # remap parameters on xi2 = 1 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [])

        for layer in range(2):
            soAround = 1
            soUp = layer*2 + 1 + 12
            ln = layer*4 + 1
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 1, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 2, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [1, soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 2, soAround + 3 , soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 1, soAround + 3, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 2
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 4, soUp + 2] ), (Node.VALUE_LABEL_D_DS2, [soAround + 5, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [1, soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 5, soAround + 6, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 4, soAround + 6, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]) ])

        # remap parameters on xi3 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 3, 4 ], Node.VALUE_LABEL_D_DS1, [])

        soAround = 1 + 6
        remapEftNodeValueLabel(eft, [ 3 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [soAround + 1]), (Node.VALUE_LABEL_D_DS3, [soAround + 2]) ])
        remapEftNodeValueLabel(eft, [ 3 ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 2, soAround + 3]), (Node.VALUE_LABEL_D_DS3, [soAround + 1, soAround + 3]) ])
        remapEftNodeValueLabel(eft, [ 4 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [soAround + 4]), (Node.VALUE_LABEL_D_DS3, [soAround + 5]) ])
        remapEftNodeValueLabel(eft, [ 4 ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 5, soAround + 6]), (Node.VALUE_LABEL_D_DS3, [soAround + 4, soAround + 6]) ])

        ln_map = [ 1, 1, 2, 2, 1, 1, 3, 4 ]
        remapEftLocalNodes(eft, 4, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftTetrahedronTop:  Failed to validate eft'
        return eft

    def createEftTetrahedronTop(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1, nodeScaleFactorOffsetUp):
        '''
        Create a tricubic hermite element field for a solid tetrahedron for the axis of a
        solid sphere pole, with xi1 and xi3 collapsed on xi2 = 1, and xi1 collapsed on xi3 = 0.
        Each collapsed node has 3 scale factors giving the cos, sin coefficients of
        the radial line from global derivatives, plus the arc subtended by the element in radians,
        so the circumferential direction is rounded. Collapsed nodes on xi2 = 1 have 2 additional
        scale factors cos and sin coefficients of the inclination angle.
        Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :param nodeScaleFactorOffsetUp: offset of first scale factor for inclination at pole,
        increase by 2 in each layer away from axis. Suggest starting at 100000 on axis.
        :return: Element field template
        '''
        # start with full tricubic
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        for n in [ 0, 1, 4, 5 ]:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            if n > 1:
                eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffsetUp + 1, nodeScaleFactorOffsetUp + 2,
            nodeScaleFactorOffsetUp + 3, nodeScaleFactorOffsetUp + 4])

        # remap parameters on xi2 = 1 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS3, [])

        for layer in range(2):
            soAround = 1 + 6
            soUp = layer*2 + 1 + 12
            ln = layer*4 + 3
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 1, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 2, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 2, soAround + 3 , soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 1, soAround + 3, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 4
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 4, soUp + 2] ), (Node.VALUE_LABEL_D_DS2, [1, soAround + 5, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 5, soAround + 6, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 4, soAround + 6, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS2, [1]) ])

        # remap parameters on xi3 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2 ], Node.VALUE_LABEL_D_DS1, [])

        soAround = 1
        remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [soAround + 1]), (Node.VALUE_LABEL_D_DS3, [soAround + 2]) ])
        remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 2, soAround + 3]), (Node.VALUE_LABEL_D_DS3, [soAround + 1, soAround + 3]) ])
        remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D_DS3, [ (Node.VALUE_LABEL_D_DS1, [soAround + 4]), (Node.VALUE_LABEL_D_DS3, [soAround + 5]) ])
        remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D2_DS1DS3, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 5, soAround + 6]), (Node.VALUE_LABEL_D_DS3, [soAround + 4, soAround + 6]) ])

        ln_map = [ 1, 1, 2, 2, 3, 4, 2, 2 ]
        remapEftLocalNodes(eft, 4, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftTetrahedronTop:  Failed to validate eft'
        return eft


    def createEftPyramidBottom(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1, nodeScaleFactorOffsetUp):
        '''
        Create a tricubic hermite element field for a solid pyramid for the axis of a
        solid sphere pole, with xi1 and xi3 collapsed on xi2 = 0. Each collapsed node
        has 5 scale factors giving the cos, sin coefficients of the radial line from
        global derivatives, plus the arc subtended by the element in radians, so the
        circumferential direction is rounded, cos and sin coefficients of the inclination
        angle. Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :param nodeScaleFactorOffsetUp: offset of first scale factor for inclination at pole,
        increase by 2 in each layer away from axis. Suggest starting at 100000 on axis.
        :return: Element field template
        '''
        # start with full tricubic
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        for n in [ 2, 3, 6, 7 ]:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffsetUp + 1, nodeScaleFactorOffsetUp + 2,
            nodeScaleFactorOffsetUp + 3, nodeScaleFactorOffsetUp + 4])

        # remap parameters on xi2 = 0 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
        remapEftNodeValueLabel(eft, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [])

        for layer in range(2):
            soAround = 1
            soUp = layer*2 + 6 + 1
            ln = layer*4 + 1
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 1, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 2, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [1, soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 2, soAround + 3 , soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 1, soAround + 3, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 2
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 4, soUp + 2] ), (Node.VALUE_LABEL_D_DS2, [soAround + 5, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [1, soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 5, soAround + 6, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [soAround + 4, soAround + 6, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        ln_map = [ 1, 1, 2, 3, 1, 1, 4, 5 ]
        remapEftLocalNodes(eft, 5, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftPyramidTop:  Failed to validate eft'
        return eft

    def createEftPyramidTop(self, nodeScaleFactorOffset0, nodeScaleFactorOffset1, nodeScaleFactorOffsetUp):
        '''
        Create a tricubic hermite element field for a solid pyramid for the axis of a
        solid sphere pole, with xi1 and xi3 collapsed on xi2 = 1. Each collapsed node
        has 5 scale factors giving the cos, sin coefficients of the radial line from
        global derivatives, plus the arc subtended by the element in radians, so the
        circumferential direction is rounded, cos and sin coefficients of the inclination
        angle. Need to create a new template for each sector around axis giving common
        nodeScaleFactorOffset values on common faces. Suggestion is to start at 0 and
        add 100 for each radial line around axis.
        :param nodeScaleFactorOffset0: offset of node scale factors at axis on xi1=0
        :param nodeScaleFactorOffset1: offset of node scale factors at axis on xi1=1
        :param nodeScaleFactorOffsetUp: offset of first scale factor for inclination at pole,
        increase by 2 in each layer away from axis. Suggest starting at 100000 on axis.
        :return: Element field template
        '''
        # start with full tricubic
        eft = self._mesh.createElementfieldtemplate(self._tricubicHermiteBasis)
        for n in [ 0, 1, 4, 5 ]:
            eft.setFunctionNumberOfTerms(n*8 + 4, 0)
            eft.setFunctionNumberOfTerms(n*8 + 6, 0)
            eft.setFunctionNumberOfTerms(n*8 + 7, 0)
            eft.setFunctionNumberOfTerms(n*8 + 8, 0)

        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        setEftScaleFactorIds(eft, [1], [
            nodeScaleFactorOffset0 + 1, nodeScaleFactorOffset0 + 2, nodeScaleFactorOffset0 + 3,
            nodeScaleFactorOffset1 + 1, nodeScaleFactorOffset1 + 2, nodeScaleFactorOffset1 + 3,
            nodeScaleFactorOffsetUp + 1, nodeScaleFactorOffsetUp + 2,
            nodeScaleFactorOffsetUp + 3, nodeScaleFactorOffsetUp + 4])

        # remap parameters on xi2 = 1 before collapsing nodes
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
        remapEftNodeValueLabel(eft, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS3, [])

        for layer in range(2):
            soAround = 1
            soUp = layer*2 + 1 + 6
            ln = layer*4 + 3
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 1, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 2, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 2, soAround + 3 , soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 1, soAround + 3, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

            ln = layer*4 + 4
            # 3 terms for d/dxi2 via general linear map:
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, [1, soAround + 4, soUp + 2] ), (Node.VALUE_LABEL_D_DS2, [1, soAround + 5, soUp + 2]), (Node.VALUE_LABEL_D_DS3, [soUp + 1]) ])
            # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS1, [soAround + 5, soAround + 6, soUp + 2]), (Node.VALUE_LABEL_D_DS2, [1, soAround + 4, soAround + 6, soUp + 2]) ])
            # zero other cross derivative parameters
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS1DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D2_DS2DS3, [])
            remapEftNodeValueLabel(eft, [ ln ], Node.VALUE_LABEL_D3_DS1DS2DS3, [])

        ln_map = [ 1, 2, 3, 3, 4, 5, 3, 3 ]
        remapEftLocalNodes(eft, 5, ln_map)

        assert eft.validate(), 'eftfactory_tricubichermite.createEftPyramidTop:  Failed to validate eft'
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

    def createEftSplitXi2RightStraight(self):
        '''
        Create an element field template suitable for the inner elements of the
        join between left and right chambers at the bottom of the RV, with xi2 bifurcating.
        Straight through version.
        :return: Element field template
        '''
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1]) ])
        remapEftNodeValueLabel(eft, [ 7, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, []) ])
        assert eft.validate(), 'eftfactory_tricubichermite.createEftSplitXi2RightStraight:  Failed to validate eft'
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

    def setEftLinearDerivative(self, eft, basisNodes, derivative, localNode1, localNode2, minus1scaleFactorIndex):
        '''
        Makes the derivative in xiIndex for basisNodes equal to the difference in
        the value parameter at localNode2 - localNode1. This makes the basis equivalent to linear Lagrange.
        Note! Cross derivatives are not handled and are currently unmodified.
        :basisNodes: List of basis nodes to set derivative at, numbering from 1.
        :param xiIndex:  Xi derivative to set: 1, 2 or 3.
        :param minus1scaleFactorIndex: Local scale factor index for general value -1.0
        '''
        assert derivative in [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ], 'setEftLinearDerivative.  Invalid derivative'
        nodeFunctionIndex = derivative - Node.VALUE_LABEL_VALUE + 1
        for n in basisNodes:
            f = (n - 1)*8 + nodeFunctionIndex
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

    def replaceElementWithInlet4(self, origElement, startElementId, nodetemplate, startNodeId, tubeLength, innerRadius, wallThickness, \
            meshGroups = [], revCorners = [], inclineRadians = -0.25*math.pi, inclineAxis = 2):
        '''
        Replace origElement with 4 element X-layout tube inlet.
        Inlet axis is at given length from centre of xi3=0 face, oriented with dx/dxi1.
        8 new nodes are created. Original element is destroyed.
        :param meshGroups:  Optional list of Zinc MeshGroup for adding new elements to.
        :param revCorners: Optional list of corners to reverse xi1 and xi2 from 1 to 4, indicating local node and local node + 4
        :param inclineRadians: Angle from centre normal toward positive inclineAxis derivative for inlet orientation.
        :param inclineAxis: 1 or 2 (surface xi direction)
        '''
        fm = self._mesh.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        fm.beginChange()
        cache = fm.createFieldcache()
        diff1 = self._mesh.getChartDifferentialoperator(1, 1)
        diff2 = self._mesh.getChartDifferentialoperator(1, 2)
        coordinates = getOrCreateCoordinateField(fm)
        cache.setMeshLocation(origElement, [0.5, 0.5, 1.0])
        result, fc = coordinates.evaluateReal(cache, 3)
        resulta, a = coordinates.evaluateDerivative(diff1, cache, 3)
        resultb, b = coordinates.evaluateDerivative(diff2, cache, 3)
        n = vector.normalise(vector.crossproduct3(a, b))
        t = vector.normalise(a if (inclineAxis == 1) else b)
        n = [ (math.cos(inclineRadians)*n[c] + math.sin(inclineRadians)*t[c]) for c in range(3) ]
        #print(resulta, 'a =', a, ',', resultb, ' b=', b, ' fc=', fc, ' n=',n)
        ic = [ (fc[i] + tubeLength*n[i]) for i in range(3) ]
        if inclineAxis == 2:
            na = vector.normalise(a)
            nb = vector.normalise(vector.crossproduct3(n, a))
        else:
            nb = vector.normalise(b)
            na = vector.normalise(vector.crossproduct3(b, n))
        a = vector.normalise([ -(na[i] + nb[i]) for i in range(3) ])
        b = vector.normalise(vector.crossproduct3(a, n))

        zero = [ 0.0, 0.0, 0.0 ]
        nodeIdentifier = startNodeId
        elementsCountAround = 4
        radiansPerElementAround = math.pi*2.0/elementsCountAround
        for n3 in range(2):
            radius = innerRadius + n3*wallThickness
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

        eft0 = origElement.getElementfieldtemplate(coordinates, -1)
        nids0 = getElementNodeIdentifiers(origElement, eft0)
        orig_nids = [ nids0[0], nids0[2], nids0[3], nids0[1], nids0[4], nids0[6], nids0[7], nids0[5] ]
        #print('orig_nids',orig_nids)

        elementIdentifier = startElementId
        elementtemplate1 = self._mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        neg1 = [] if (1 in revCorners) else [1]
        pos1 = [1] if (1 in revCorners) else []
        neg2 = [] if (2 in revCorners) else [1]
        pos2 = [1] if (2 in revCorners) else []
        neg3 = [] if (3 in revCorners) else [1]
        pos3 = [1] if (3 in revCorners) else []
        neg4 = [] if (4 in revCorners) else [1]
        pos4 = [1] if (4 in revCorners) else []
        #print(elementIdentifier, revCorners, '->',neg1,pos1,neg2,pos2,neg3,pos3,neg4,pos4)

        for e in range(4):
            eft1 = self.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            if e == 0:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg1), (Node.VALUE_LABEL_D_DS2, neg1) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, pos1) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg3), (Node.VALUE_LABEL_D_DS2, pos3) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, pos3) ])
            elif e == 1:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg3), (Node.VALUE_LABEL_D_DS2, pos3) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, pos3) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos4), (Node.VALUE_LABEL_D_DS2, pos4) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, pos4) ])
            elif e == 2:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos4), (Node.VALUE_LABEL_D_DS2, pos4) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, neg4) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos2), (Node.VALUE_LABEL_D_DS2, neg2) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, neg2) ])
            elif e == 3:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos2), (Node.VALUE_LABEL_D_DS2, neg2) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, neg2) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg1), (Node.VALUE_LABEL_D_DS2, neg1) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, neg1) ])
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
            #print('inlet element', element.isValid(), elementIdentifier, elementtemplate1.isValid())
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            else:
                result3 = 7
            #print('create element inlet4', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        self._mesh.destroyElement(origElement)
        fm.endChange()

    def replaceTwoElementWithInlet6(self, origElement1, origElement2, startElementId, nodetemplate, startNodeId, inletCentre, inletAxis, inletSide, innerRadius, wallThickness, meshGroups = [], revCorners = []):
        '''
        Replace origElement1 and origElement2 (on xi2=1 side of origElement1) with 6 element tube inlet.
        Inlet is centred at inletCentre, inward direction and derivative magnitude given by inletAxis,
        with inletSide pointing in -ve xi1 direction for both origElements.
        12 new nodes are created. Original elements are destroyed.
        :param meshGroups:  Optional list of Zinc MeshGroup for adding new elements to.
        :param revCorners: Optional list of corners to reverse xi1 and xi2 from 1 to 6, indicating local node and local node + 4
        '''
        fm = self._mesh.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        fm.beginChange()
        cache = fm.createFieldcache()
        coordinates = getOrCreateCoordinateField(fm)
        a = vector.normalise(inletSide)
        b = vector.normalise(vector.crossproduct3(inletAxis, inletSide))

        zero = [ 0.0, 0.0, 0.0 ]
        nodeIdentifier = startNodeId
        elementsCountAround = 6
        radiansPerElementAround = math.pi*2.0/elementsCountAround
        for n3 in range(2):
            radius = innerRadius + n3*wallThickness
            for n1 in range(elementsCountAround):
                radiansAround = (n1 - 1)*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x = [ (inletCentre[c] + radius*(cosRadiansAround*a[c] + sinRadiansAround*b[c])) for c in range(3) ]
                dx_ds1 = [ radiansPerElementAround*radius*(-sinRadiansAround*a[c] + cosRadiansAround*b[c]) for c in range(3) ]
                dx_ds2 = inletAxis
                dx_ds3 = [ wallThickness*(cosRadiansAround*a[c] + sinRadiansAround*b[c]) for c in range(3) ]
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

        eft0 = origElement1.getElementfieldtemplate(coordinates, -1)
        nids1 = getElementNodeIdentifiers(origElement1, eft0)
        eft0 = origElement2.getElementfieldtemplate(coordinates, -1)
        nids2 = getElementNodeIdentifiers(origElement2, eft0)
        orig_nids = [ nids1[0], nids1[2], nids2[2], nids2[3], nids1[3], nids1[1],
                      nids1[4], nids1[6], nids2[6], nids2[7], nids1[7], nids1[5] ]
        #print('orig_nids',orig_nids)

        elementIdentifier = startElementId
        elementtemplate1 = self._mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        neg1 = [] if (1 in revCorners) else [1]
        pos1 = [1] if (1 in revCorners) else []
        neg2 = [] if (2 in revCorners) else [1]
        pos2 = [1] if (2 in revCorners) else []
        neg3 = [] if (3 in revCorners) else [1]
        pos3 = [1] if (3 in revCorners) else []
        neg4 = [] if (4 in revCorners) else [1]
        pos4 = [1] if (4 in revCorners) else []
        neg5 = [] if (5 in revCorners) else [1]
        pos5 = [1] if (5 in revCorners) else []
        neg6 = [] if (6 in revCorners) else [1]
        pos6 = [1] if (6 in revCorners) else []

        for e in range(6):
            eft1 = self.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            if e == 0:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg1), (Node.VALUE_LABEL_D_DS2, neg1) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, pos1) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D2_DS1DS2, []) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg3) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS2, pos3) ])  # finish swap
            elif e == 1:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D2_DS1DS2, []) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg3) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS2, pos3) ])  # finish swap
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg5), (Node.VALUE_LABEL_D_DS2, pos5) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, pos5) ])
            elif e == 2:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg5), (Node.VALUE_LABEL_D_DS2, pos5) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, pos5) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos6), (Node.VALUE_LABEL_D_DS2, pos6) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, pos6) ])
            elif e == 3:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos6), (Node.VALUE_LABEL_D_DS2, pos6) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, neg6) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D2_DS1DS2, []) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos4) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS2, neg4) ])  # finish swap
            elif e == 4:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D2_DS1DS2, []) ])  # temporary to enable swap
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos4) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D2_DS1DS2, [ (Node.VALUE_LABEL_D_DS2, neg4) ])  # finish swap
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos2), (Node.VALUE_LABEL_D_DS2, neg2) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS2, neg2) ])
            elif e == 5:
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, pos2), (Node.VALUE_LABEL_D_DS2, neg2) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, neg2) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ (Node.VALUE_LABEL_D_DS1, neg1), (Node.VALUE_LABEL_D_DS2, neg1) ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ (Node.VALUE_LABEL_D_DS1, neg1) ])
            ea = e
            eb = (e + 1) % elementsCountAround
            ec = ea + elementsCountAround
            ed = eb + elementsCountAround
            nids = [
                startNodeId + ea, startNodeId + eb, orig_nids[ea], orig_nids[eb],
                startNodeId + ec, startNodeId + ed, orig_nids[ec], orig_nids[ed]
            ]
            elementtemplate1.defineField(coordinates, -1, eft1)
            element = self._mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, [ -1.0 ])
            #print('create element inlet6', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        self._mesh.destroyElement(origElement1)
        self._mesh.destroyElement(origElement2)
        fm.endChange()

    def createAnnulusMesh3d(self,
            startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
            endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
            nodetemplate, nodetemplateLinearS3, nextNodeIdentifier, nextElementIdentifier,
            elementsCountRadial = 1, meshGroups = []):
        '''
        Create an annular mesh from a loop of start points/nodes with specified derivative mappings to
        a loop of end points/nodes with specified derivative mappings.
        Derivative d3 is through the wall. Currently limited to single element layer through wall.
        Points/nodes order cycles fastest around the annulus, then through the wall.
        Note doesn't support cross derivatives.
        Arrays are indexed by n3 (node through wall, size 2), n2 (node along/radial), n1 (node around, variable size)
        and coordinate component c.
        :param startPointsx, startPointsd1, startPointsd2, startPointsd3, endPointsx, endPointsd1, endPointsd2, endPointsd3:
            List array[n3][n1][c] or start/point coordinates and derivatives. To linearise through the wall, pass None to d3.
        :param startNodeId, endNodeId: List array [n3][n1] of existing node identifiers to use at start/end. Pass None for
            argument if no nodes are specified at end. These arguments are 'all or nothing'.
        :param startDerivativesMap, endDerivativesMap: List array[n3][n1] of mappings for d/dxi1, d/dxi2, d/dxi3 at start/end of form:
            ( (1, -1, 0), (1, 0, 0), None ) where the first tuple means d/dxi1 = d/ds1 - d/ds2. Only 0, 1 and -1 may be used.
            None means use default e.g. d/dxi2 = d/ds2.
            Pass None for the entire argument to use the defaults d/dxi1 = d/ds1, d/dxi2 = d/ds2, d/dxi3 = d/ds3.
            Pass a 4th mapping to apply to d/dxi1 on other side of node; if not supplied first mapping applies both sides.
        :param nodetemplate: Full tricubic Hermite node template.
        :param nodetemplateLinearS3: Node template to use where linear through-thickness.
        :param nextNodeIdentifier, nextElementIdentifier: Next identifiers to use and increment.
        :param elementsCountRadial: Optional number of elements in radial direction between start and end.
        :param meshGroups:  Optional list of Zinc MeshGroup for adding new elements to.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier
        '''
        assert (elementsCountRadial >= 1), 'eftfactory_tricubichermite.createAnnulusMesh3d:  Invalid number of radial elements'
        elementsCountWall = 1
        nodesCountWall = elementsCountWall + 1
        assert (len(startPointsx) == nodesCountWall) and (len(startPointsd1) == nodesCountWall) and (len(startPointsd2) == nodesCountWall) and \
            ((startPointsd3 is None) or (len(startPointsd3) == nodesCountWall)) and \
            (len(endPointsx) == nodesCountWall) and (len(endPointsd1) == nodesCountWall) and (len(endPointsd2) == nodesCountWall) and \
            ((endPointsd3 is None) or (len(endPointsd3) == nodesCountWall)) and \
            ((startNodeId is None) or (len(startNodeId) == nodesCountWall)) and \
            ((endNodeId is None) or (len(endNodeId) == nodesCountWall)) and \
            ((startDerivativesMap is None) or (len(startDerivativesMap) == nodesCountWall)) and \
            ((endDerivativesMap is None) or (len(endDerivativesMap) == nodesCountWall)), \
            'eftfactory_tricubichermite.createAnnulusMesh3d:  Mismatch in number of layers through wall'
        elementsCountAround = nodesCountAround = len(startPointsx[0])
        assert (nodesCountAround > 1), 'eftfactory_tricubichermite.createAnnulusMesh3d:  Invalid number of points/nodes around annulus'
        for n3 in range(nodesCountWall):
            assert (len(startPointsx[n3]) == nodesCountAround) and (len(startPointsd1[n3]) == nodesCountAround) and (len(startPointsd2[n3]) == nodesCountAround) and \
                ((startPointsd3 is None) or (len(startPointsd3[n3]) == nodesCountAround)) and \
                (len(endPointsx[n3]) == nodesCountAround) and (len(endPointsd1[n3]) == nodesCountAround) and (len(endPointsd2[n3]) == nodesCountAround) and \
                ((endPointsd3 is None) or (len(endPointsd3[n3]) == nodesCountAround)) and \
                ((startNodeId is None) or (len(startNodeId[n3]) == nodesCountAround)) and \
                ((endNodeId is None) or (len(endNodeId[n3]) == nodesCountAround)) and \
                ((startDerivativesMap is None) or (len(startDerivativesMap[n3]) == nodesCountAround)) and \
                ((endDerivativesMap is None) or (len(endDerivativesMap[n3]) == nodesCountAround)), \
                'eftfactory_tricubichermite.createAnnulusMesh3d:  Mismatch in number of points/nodes in layers through wall'

        fm = self._mesh.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        coordinates = getOrCreateCoordinateField(fm)

        # Build arrays of points from start to end
        px  = [ [], [] ]
        pd1 = [ [], [] ]
        pd2 = [ [], [] ]
        pd3 = [ [], [] ]
        for n3 in range(2):
            px [n3] = [ startPointsx [n3], endPointsx [n3] ]
            pd1[n3] = [ startPointsd1[n3], endPointsd1[n3] ]
            pd2[n3] = [ startPointsd2[n3], endPointsd2[n3] ]
            pd3[n3] = [ startPointsd3[n3] if (startPointsd3 is not None) else None, \
                        endPointsd3[n3] if (endPointsd3 is not None) else None ]
        if elementsCountRadial > 1:
            # add in-between points
            startPointsd = [ startPointsd1, startPointsd2, startPointsd3 ]
            startPointsdslimit = 2 if (startPointsd3 is None) else 3
            endPointsd = [ endPointsd1, endPointsd2, endPointsd3 ]
            endPointsdslimit = 2 if (endPointsd3 is None) else 3
            for n3 in range(2):
                for n2 in range(1, elementsCountRadial):
                    for li in (px[n3], pd1[n3], pd2[n3], pd3[n3]):
                        li.insert(n2, [ None ]*nodesCountAround)
                for n1 in range(nodesCountAround):
                    ax  = startPointsx [n3][n1]
                    if (startDerivativesMap is None) or (startDerivativesMap[n3][n1][0] is None):
                        ad1 = startPointsd1[n3][n1]
                    else:
                        derivativesMap = startDerivativesMap[n3][n1][0]
                        ad1 = [ 0.0, 0.0, 0.0 ]
                        for ds in range(startPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    ad1[c] += derivativesMap[ds]*startPointsd[ds][n3][n1][c]
                        if len(startDerivativesMap[n3][n1]) > 3:
                            # average with d1 map for other side
                            derivativesMap = startDerivativesMap[n3][n1][3]
                            ad1 = [ 0.5*d for d in ad1 ]
                            for ds in range(startPointsdslimit):
                                if derivativesMap[ds] != 0.0:
                                    for c in range(3):
                                        ad1[c] += 0.5*derivativesMap[ds]*startPointsd[ds][n3][n1][c]
                    if (startDerivativesMap is None) or (startDerivativesMap[n3][n1][1] is None):
                        ad2 = startPointsd2[n3][n1]
                    else:
                        derivativesMap = startDerivativesMap[n3][n1][1]
                        ad2 = [ 0.0, 0.0, 0.0 ]
                        for ds in range(startPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    ad2[c] += derivativesMap[ds]*startPointsd[ds][n3][n1][c]

                    bx  = endPointsx [n3][n1]
                    if (endDerivativesMap is None) or (endDerivativesMap[n3][n1][0] is None):
                        bd1 = endPointsd1[n3][n1]
                    else:
                        derivativesMap = endDerivativesMap[n3][n1][0]
                        bd1 = [ 0.0, 0.0, 0.0 ]
                        for ds in range(endPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    bd1[c] += derivativesMap[ds]*endPointsd[ds][n3][n1][c]
                        if len(endDerivativesMap[n3][n1]) > 3:
                            # average with d1 map for other side
                            derivativesMap = endDerivativesMap[n3][n1][3]
                            bd1 = [ 0.5*d for d in bd1 ]
                            for ds in range(endPointsdslimit):
                                if derivativesMap[ds] != 0.0:
                                    for c in range(3):
                                        bd1[c] += 0.5*derivativesMap[ds]*endPointsd[ds][n3][n1][c]
                    if (endDerivativesMap is None) or (endDerivativesMap[n3][n1][1] is None):
                        bd2 = endPointsd2[n3][n1]
                    else:
                        derivativesMap = endDerivativesMap[n3][n1][1]
                        bd2 = [ 0.0, 0.0, 0.0 ]
                        for ds in range(endPointsdslimit):
                            if derivativesMap[ds] != 0.0:
                                for c in range(3):
                                    bd2[c] += derivativesMap[ds]*endPointsd[ds][n3][n1][c]

                    mx, md2, md1 = sampleCubicHermiteCurves([ ax, bx ], [ ad2, bd2 ], [ ad1, bd1 ], elementsCountRadial,
                        addLengthStart = 0.5*vector.magnitude(ad2), lengthFractionStart = 0.5,
                        addLengthEnd = 0.5*vector.magnitude(bd2), lengthFractionEnd = 0.5)
                    for n2 in range(1, elementsCountRadial):
                        px [n3][n2][n1] = mx [n2]
                        pd1[n3][n2][n1] = md1[n2]
                        pd2[n3][n2][n1] = md2[n2]
                        # get derivative 3 from difference across wall
                        if n3 == 1:
                            pd3[1][n2][n1] = pd3[0][n2][n1] = [ (px[1][n2][n1][c] - px[0][n2][n1][c]) for c in range(3) ]

        ##############
        # Create nodes
        ##############

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = nextNodeIdentifier
        nodeId = [ [], [] ]
        for n3 in range(2):
            for n2 in range(elementsCountRadial + 1):
                if (n2 == 0) and (startNodeId is not None):
                    rowNodeId = startNodeId[n3]
                elif (n2 == elementsCountRadial) and (endNodeId is not None):
                    rowNodeId = endNodeId[n3]
                else:
                    rowNodeId = []
                    nodetemplate1 = nodetemplate if (pd3[n3][n2] is not None) else nodetemplateLinearS3
                    for n1 in range(nodesCountAround):
                        node = nodes.createNode(nodeIdentifier, nodetemplate1)
                        rowNodeId.append(nodeIdentifier)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, px[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pd2[n3][n2][n1])
                        if nodetemplate1 is nodetemplate:
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, pd3[n3][n2][n1])
                        nodeIdentifier = nodeIdentifier + 1
                nodeId[n3].append(rowNodeId)

        #################
        # Create elements
        #################

        elementIdentifier = nextElementIdentifier

        eftStandard = self.createEftNoCrossDerivatives()

        elementtemplate1 = self._mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        for e2 in range(elementsCountRadial):
            for e1 in range(elementsCountAround):
                en = (e1 + 1)%elementsCountAround
                mapStartDerivatives = (e2 == 0) and (startDerivativesMap is not None)
                mapStartLinearDerivative3 = (e2 == 0) and (startPointsd3 is None)
                mapEndDerivatives = (e2 == (elementsCountRadial - 1)) and (endDerivativesMap is not None)
                mapEndLinearDerivative3 = (e2 == (elementsCountRadial - 1)) and (endPointsd3 is None)

                if mapStartDerivatives or mapStartLinearDerivative3 or mapEndDerivatives or mapEndLinearDerivative3:
                    eft1 = self.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    if mapStartLinearDerivative3:
                        self.setEftLinearDerivative(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS3, 1, 5, 1)
                        self.setEftLinearDerivative(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS3, 2, 6, 1)
                    if mapStartDerivatives:
                        for i in range(2):
                            lns = [ 1, 5 ] if (i == 0) else [ 2, 6 ]
                            for n3 in range(2):
                                derivativesMap = startDerivativesMap[n3][e1] if (i == 0) else startDerivativesMap[n3][en]
                                # handle different d1 on each side of node
                                d1Map = derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else derivativesMap[3]
                                d2Map = derivativesMap[1]
                                d3Map = derivativesMap[2]
                                # use temporary to safely swap DS1 and DS2:
                                ln = [ lns[n3] ]
                                if d1Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])
                                if d3Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS2DS3, [] ) ])
                                if d2Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d2Map))
                                if d1Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d1Map))
                                if d3Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d3Map))
                    if mapEndLinearDerivative3:
                        self.setEftLinearDerivative(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS3, 3, 7, 1)
                        self.setEftLinearDerivative(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS3, 4, 8, 1)
                    if mapEndDerivatives:
                        for i in range(2):
                            lns = [ 3, 7 ] if (i == 0) else [ 4, 8 ]
                            for n3 in range(2):
                                derivativesMap = endDerivativesMap[n3][e1] if (i == 0) else endDerivativesMap[n3][en]
                                # handle different d1 on each side of node
                                d1Map = derivativesMap[0] if ((i == 1) or (len(derivativesMap) < 4)) else derivativesMap[3]
                                d2Map = derivativesMap[1]
                                d3Map = derivativesMap[2]
                                # use temporary to safely swap DS1 and DS2:
                                ln = [ lns[n3] ]
                                if d1Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])
                                if d3Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D2_DS2DS3, [] ) ])
                                if d2Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d2Map))
                                if d1Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS1DS2, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d1Map))
                                if d3Map is not None:
                                    remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D2_DS2DS3, \
                                        derivativeSignsToExpressionTerms( ( Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ), d3Map))
                else:
                    eft1 = eftStandard

                nids = [ nodeId[0][e2][e1], nodeId[0][e2][en], nodeId[0][e2 + 1][e1], nodeId[0][e2 + 1][en],
                         nodeId[1][e2][e1], nodeId[1][e2][en], nodeId[1][e2 + 1][e1], nodeId[1][e2 + 1][en] ]

                elementtemplate1.defineField(coordinates, -1, eft1)
                element = self._mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if eft1 is not eftStandard:
                    result3 = element.setScaleFactors(eft1, [ -1.0 ])
                else:
                    result3 = 7
                #print('create element orifice', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        fm.endChange()

        return nodeIdentifier, elementIdentifier
