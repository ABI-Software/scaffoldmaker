'''
Definitions of standard element field templates using bicubic Hermite x linear Lagrange basis.

@author: Richard Christie
'''
from scaffoldmaker.utils.eft_utils import *
from opencmiss.zinc.element import Elementbasis, Elementfieldtemplate
from opencmiss.zinc.node import Node
from opencmiss.zinc.status import OK as ZINC_OK

class eftfactory_bicubichermitelinear:
    '''
    Factory class for creating element field templates for a 3-D mesh using bicubic Hermite x linear Lagrange basis.
    '''

    def __init__(self, mesh, useCrossDerivatives, linearAxis = 3,
            d_ds1 = Node.VALUE_LABEL_D_DS1, d_ds2 = Node.VALUE_LABEL_D_DS2):
        '''
        :param mesh:  Zinc mesh to create element field templates in.
        :param useCrossDerivatives: Set to True if you want cross derivative terms.
        :param linearAxis: 1, 2, or 3.
        :param d_ds1: Node derivative to use in Hermite axis 1: Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2.
        :param d_ds2: Node derivative to use in Hermite axis 2, > d_ds1: Node.VALUE_LABEL_D_DS2 or Node.VALUE_LABEL_D_DS3.
        '''
        assert mesh.getDimension() == 3, 'eftfactory_bicubichermitelinear: not a 3-D Zinc mesh'
        assert linearAxis in [ 1, 2, 3 ], 'eftfactory_bicubichermitelinear: linearAxis must be 1, 2 or 3'
        assert d_ds1 in [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], 'eftfactory_bicubichermitelinear: invalid d_ds1'
        assert d_ds2 in [ Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ] and (d_ds2 > d_ds1), 'eftfactory_bicubichermitelinear: invalid d_ds2'
        self._mesh = mesh
        self._useCrossDerivatives = useCrossDerivatives
        self._linearAxis = linearAxis
        self._d_ds1 = d_ds1
        self._d_ds2 = d_ds2
        self._d2_ds1ds2 = Node.VALUE_LABEL_D2_DS2DS3 if (d_ds1 == Node.VALUE_LABEL_D_DS2) \
                     else Node.VALUE_LABEL_D2_DS1DS3 if (d_ds2 == Node.VALUE_LABEL_D_DS3) \
                     else Node.VALUE_LABEL_D2_DS1DS2
        self._fieldmodule = mesh.getFieldmodule()
        self._basis = self._fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        self._basis.setFunctionType(linearAxis, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)

    def _remapDefaultNodeDerivatives(self, eft):
        '''
        Remap the Hermite node derivatives to those chosen in __init__.
        Use only on first create.
        :param eft: The element field template to remap.
        '''
        # must do d_ds2 first!
        if self._d_ds2 != Node.VALUE_LABEL_D_DS2:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D_DS2, [ (self._d_ds2, []) ])
        if self._d_ds1 != Node.VALUE_LABEL_D_DS1:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D_DS1, [ (self._d_ds1, []) ])
        if self._d2_ds1ds2 != Node.VALUE_LABEL_D2_DS1DS2:
            remapEftNodeValueLabel(eft, range(1, 9), Node.VALUE_LABEL_D2_DS1DS2, [ (self._d2_ds1ds2, []) ])

    def createEftBasic(self):
        '''
        Create the basic biicubic Hermite x linear Lagrange element template with 1:1 mappings to
        node derivatives ds1 & ds2, with or without cross derivatives accordinate as initialised.
        :return: Element field template
        '''
        if not self._useCrossDerivatives:
            return self.createEftNoCrossDerivatives()
        eft = self._mesh.createElementfieldtemplate(self._basis)
        self._remapDefaultNodeDerivatives(eft)
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftBasic:  Failed to validate eft'
        return eft

    def createEftNoCrossDerivatives(self):
        '''
        Create a basic tricubic hermite element template with 1:1 mappings to
        node derivatives ds1 & ds2, without cross derivatives.
        :return: Element field template
        '''
        eft = self._mesh.createElementfieldtemplate(self._basis)
        for n in range(8):
            eft.setFunctionNumberOfTerms(n*4 + 4, 0)
        self._remapDefaultNodeDerivatives(eft)
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftNoCrossDerivatives:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightStraight(self):
        '''
        Create an element field template suitable for the inner elements of the
        join between left and right chambers, with xi1 bifurcating to right.
        Straight through version.
        Only works with linearAxis 2.
        :return: Element field template
        '''
        assert linearAxis == 2, 'eftfactory_bicubichermitelinear.createEftSplitXi1RightStraight:  Not linearAxis 2'
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 5, 7 ], self._d_ds1, [ (self._d_ds1, []), (self._d_ds2, [1]) ])
        remapEftNodeValueLabel(eft, [ 6, 8 ], self._d_ds1, [ (self._d_ds1, []), (self._d_ds2, []) ])
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftSplitXi1RightStraight:  Failed to validate eft'
        return eft

    def createEftSplitXi1RightOut(self):
        '''
        Create an element field template suitable for the outer elements of the
        join between left and right chambers, with xi1 bifurcating to right.
        Right out version i.e. xi1 heading to right. h-shape.
        Only works with linearAxis 2.
        :return: Element field template
        '''
        assert linearAxis == 2, 'eftfactory_bicubichermitelinear.createEftSplitXi1RightOut:  Not linearAxis 2'
        eft = self.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eft, [1], [])
        remapEftNodeValueLabel(eft, [ 1, 3 ], self._d_ds1, [ (self._d_ds1, [1]) ])
        remapEftNodeValueLabel(eft, [ 1, 3 ], self._d_ds2, [ (self._d_ds1, [1]), (self._d_ds2, [1]) ])
        remapEftNodeValueLabel(eft, [ 5, 7 ], self._d_ds2, [ (self._d_ds1, [1]), (self._d_ds2, []) ])
        assert eft.validate(), 'eftfactory_bicubichermitelinear.createEftSplitXi1RightOut:  Failed to validate eft'
        return eft
