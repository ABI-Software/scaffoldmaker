"""
Generates a 1-D path mesh.
"""

from __future__ import division

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.interpolation import smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.zinc_utils import make_nodeset_derivatives_orthogonal, \
    get_nodeset_path_field_parameters, setPathParameters

class MeshType_1d_path1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '1D Path 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Coordinate dimensions' : 3,
            'D2 derivatives': False,
            'D3 derivatives': False,
            'Length' : 1.0,
            'Number of elements' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'D2 derivatives',
            'D3 derivatives',
            'Length',
            'Number of elements'
        ]

    @staticmethod
    def checkOptions(options):
        dependentChanges = False
        if (options['Coordinate dimensions'] < 1) :
            options['Coordinate dimensions'] = 1
        elif (options['Coordinate dimensions'] > 3) :
            options['Coordinate dimensions'] = 3
        if (options['Number of elements'] < 1) :
            options['Number of elements'] = 1

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, None
        """
        coordinateDimensions = options['Coordinate dimensions']
        d2derivatives = options['D2 derivatives']
        d3derivatives = options['D3 derivatives']
        length = options['Length']
        elementsCount = options['Number of elements']

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule, components_count=coordinateDimensions)
        cache = fieldmodule.createFieldcache()

        #################
        # Create nodes
        #################

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        if d2derivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if d3derivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)

        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ length/elementsCount, 0.0, 0.0 ]
        if d2derivatives:
            dx_ds2 = [ 0.0, 1.0, 0.0 ]
            d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        if d3derivatives:
            dx_ds3 = [ 0.0, 0.0, 1.0 ]
            d2x_ds1ds3 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCount + 1):
            x[0] = length*n/elementsCount
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            if d2derivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            if d3derivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, d2x_ds1ds3)
            nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################

        mesh = fieldmodule.findMeshByDimension(1)
        cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        for e in range(elementsCount):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1, e + 2 ])
            elementIdentifier = elementIdentifier + 1

        return [], None

    @classmethod
    def makeSideDerivativesNormal(cls, region, options, constructionObject, functionOptions, editGroupName):
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        makeD2Normal = options['D2 derivatives'] and functionOptions['Make D2 normal']
        makeD3Normal = options['D3 derivatives'] and functionOptions['Make D3 normal']
        if not (makeD2Normal or makeD3Normal):
            return False, False
        make_nodeset_derivatives_orthogonal(nodeset, coordinates, makeD2Normal, makeD3Normal, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def smoothSideCrossDerivatives(cls, region, options, constructionObject, functionOptions, editGroupName):
        smoothD12 = options['D2 derivatives'] and functionOptions['Smooth D12']
        smoothD13 = options['D3 derivatives'] and functionOptions['Smooth D13']
        if not (smoothD12 or smoothD13):
            return False, False
        valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
        if smoothD12:
            valueLabels += [Node.VALUE_LABEL_D_DS2]
        if smoothD13:
            valueLabels += [Node.VALUE_LABEL_D_DS3]
        fieldmodule = region.getFieldmodule()
        parameters = get_nodeset_path_field_parameters(
            fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES),
            fieldmodule.findFieldByName('coordinates'),
            valueLabels)
        x = parameters[0]
        d1 = parameters[1]
        nsv = []
        if smoothD12:
            nsv.append(parameters[2])
        if smoothD13:
            nsv.append(parameters[-1])
        dnsv = smoothCurveSideCrossDerivatives(x, d1, nsv)
        modifyParameters = []
        modifyValueLabels = []
        if smoothD12:
            modifyParameters.append(dnsv[0])
            modifyValueLabels.append(Node.VALUE_LABEL_D2_DS1DS2)
        if smoothD13:
            modifyParameters.append(dnsv[-1])
            modifyValueLabels.append(Node.VALUE_LABEL_D2_DS1DS3)
        setPathParameters(region, modifyValueLabels, modifyParameters, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Supply client with functions for smoothing path parameters.
        """
        return Scaffold_base.getInteractiveFunctions() + [
            ("Make side derivatives normal...",
                { 'Make D2 normal': True,
                  'Make D3 normal': True },
                lambda region, options, constructionObject, functionOptions, editGroupName:
                    cls.makeSideDerivativesNormal(region, options, constructionObject, functionOptions, editGroupName)),
            ("Smooth side cross derivatives...",
                { 'Smooth D12' : True,
                  'Smooth D13' : True },
                lambda region, options, constructionObject, functionOptions, editGroupName:
                    cls.smoothSideCrossDerivatives(region, options, constructionObject, functionOptions, editGroupName))
        ]
