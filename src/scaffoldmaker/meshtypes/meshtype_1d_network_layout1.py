"""
Constructs a 1-D network layout mesh with specifiable structure.
"""
from cmlibs.utils.zinc.field import findOrCreateFieldGroup
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.interpolation import smoothCubicHermiteCrossDerivativesLine
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.zinc_utils import make_nodeset_derivatives_orthogonal, \
    get_nodeset_path_field_parameters, setPathParameters


class MeshType_1d_network_layout1(Scaffold_base):
    """
    """
    parameterSetStructureStrings = {
        "Default": "1-2",
        "Bifurcation": "1-2,2-3,2.2-4"
    }

    @classmethod
    def getName(cls):
        return "1D Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default", "Bifurcation"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Structure": cls.parameterSetStructureStrings[parameterSetName]
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Structure"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, NetworkMesh
        """
        structure = options["Structure"]

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        return [], networkMesh

    @classmethod
    def makeSideDerivativesNormal(cls, region, options, functionOptions, editGroupName):
        """
        Make side directions normal to d1 and each other. Works for all versions.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        makeD2Normal = functionOptions['Make D2 normal']
        makeD3Normal = functionOptions['Make D3 normal']
        if not (makeD2Normal or makeD3Normal):
            return False, False
        make_nodeset_derivatives_orthogonal(nodeset, coordinates, makeD2Normal, makeD3Normal, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def smoothSideCrossDerivatives(cls, region, options, functionOptions, editGroupName):
        """
        Smooth side cross derivatives giving rate of change of side directions d2, d3 w.r.t. d1.
        Note: only works for a single path with version 1.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param functionOptions: Which side derivatives to smooth.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        smoothD12 = functionOptions["Smooth D12"]
        smoothD13 = functionOptions["Smooth D13"]
        if not (smoothD12 or smoothD13):
            return False, False
        valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
        if smoothD12:
            valueLabels += [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2]
        if smoothD13:
            valueLabels += [Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]
        fieldmodule = region.getFieldmodule()
        parameters = get_nodeset_path_field_parameters(
            fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES),
            fieldmodule.findFieldByName('coordinates'),
            valueLabels)
        x = parameters[0]
        d1 = parameters[1]
        modifyParameters = []
        modifyValueLabels = []
        if smoothD12:
            d12 = smoothCubicHermiteCrossDerivativesLine(x, d1, parameters[2], parameters[3])
            modifyParameters.append(d12)
            modifyValueLabels.append(Node.VALUE_LABEL_D2_DS1DS2)
        if smoothD13:
            d13 = smoothCubicHermiteCrossDerivativesLine(x, d1, parameters[-2], parameters[-1])
            modifyParameters.append(d13)
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
                {"Make D2 normal": True,
                 "Make D3 normal": True},
                lambda region, options, functionOptions, editGroupName:
                    cls.makeSideDerivativesNormal(region, options, functionOptions, editGroupName)),
            ("Smooth side cross derivatives...",
                {"Smooth D12": True,
                 "Smooth D13": True},
                lambda region, options, functionOptions, editGroupName:
                    cls.smoothSideCrossDerivatives(region, options, functionOptions, editGroupName))
        ]
