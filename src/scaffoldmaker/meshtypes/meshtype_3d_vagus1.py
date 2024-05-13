"""
Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches.
"""

import copy

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, get_group_list, find_or_create_field_group, \
    create_field_group, get_managed_field_names
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    findAnnotationGroupByName
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from cmlibs.maths.vectorops import add, sub, mult, div
from scaffoldmaker.utils.vector import magnitude_squared, magnitude
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, interpolateCubicHermite, interpolateLagrangeHermite, \
    interpolateHermiteLagrange, sampleCubicHermiteCurves
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters, print_node_field_parameters
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelWithNodes
from scaffoldfitter.fitter import Fitter
from scaffoldfitter.fitterstepconfig import FitterStepConfig
from scaffoldfitter.fitterstepalign import FitterStepAlign
from scaffoldfitter.fitterstepfit import FitterStepFit


class MeshType_3d_vagus1(Scaffold_base):
    """
    Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches.
    """

    @staticmethod
    def getName():
        return "3D Vagus 1"

    @staticmethod
    def getParameterSetNames():
        return [
            'Human Left Trunk 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk': 45
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along the trunk',
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options['Number of elements along the trunk'] < 10:
            options['Number of elements along the trunk'] = 10
        return dependentChanges


    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        applyFitting = True
        addBranches = True
        xmlData = False

        marker_TermNameVagusLengthList = {
            "level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
            "level of superior border of the jugular foramen on the vagal trunk": 8.6342,
            "level of inferior border of the jugular foramen on the vagal trunk": 16.7227,
            "level of C1 transverse process on the vagal trunk": 32.1129,
            "level of angle of mandible on the vagal trunk": 42.2450,
            "level of tubercles of the greater horn of hyoid bone on the vagal trunk": 45.6122,
            "level of carotid bifurcation on the vagal trunk": 48.3581,
            "level of laryngeal prominence on the vagal trunk": 68.8431,
            "level of superior border of clavicle on the vagal trunk": 117.5627,
            "level of jugular notch on the vagal trunk": 124.6407,
            "level of sternal angle on the vagal trunk": 151.2352,
            "1 cm superior to esophageal plexus on the vagal trunk": 165.5876,
            "level of esophageal hiatus on the vagal trunk": 254.32879,
            "level of aortic hiatus on the vagal trunk": 291.3695,
            "level of end of trunk on Japanese dataset": 312.5 # note this term is also not on the list of annotations
        }

        lengthToDiameterRatio = 312.5 # calculated from total length of nerve/average diameter of nerve
        rescaled_lengthToDiameterRatio = 100.0
        rescaledmarker_TermNameVagusLengthList = {}
        for term in marker_TermNameVagusLengthList:
            rescaledmarker_TermNameVagusLengthList[term] = marker_TermNameVagusLengthList[term] / lengthToDiameterRatio * rescaled_lengthToDiameterRatio

        # load data from file
        if xmlData:
            trunk_group_name, trunk_data, marker_data, branch_data = load_data_xml(region)
        else:
            trunk_group_name, trunk_data, marker_data, branch_data = load_data(region)
        assert len(marker_data) >= 2, f"At least two landmarks are expected in the data."

        # choose markers for building initial scaffold
        use_marker_names = [list(marker_data.keys())[0],
                            list(marker_data.keys())[-1]]
        assert [name.lower() in marker_data for name in use_marker_names] and \
               [name.lower() in rescaledmarker_TermNameVagusLengthList for name in use_marker_names]

        # setup
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # geometric coordinates
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        mesh = fieldmodule.findMeshByDimension(1)
        elementbasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(elementbasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        elementtemplate.defineField(coordinates, -1, eft)

        # material coordinates
        vagusCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="vagus coordinates")
        vagusNodetemplate = nodes.createNodetemplate()
        vagusNodetemplate.defineField(vagusCoordinates)
        vagusNodetemplate.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        vagusNodetemplate.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        eftVagus = mesh.createElementfieldtemplate(elementbasis)
        vagusElementtemplate = mesh.createElementtemplate()
        vagusElementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        vagusElementtemplate.defineField(vagusCoordinates, -1, eftVagus)

        # branch special node and element field template - geometric coordinates
        nodeTemplateNoValue = nodes.createNodetemplate()
        nodeTemplateNoValue.defineField(coordinates)
        nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 0)
        nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eftNV = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eftNV.setNumberOfLocalNodes(4)
        eftNV.setNumberOfLocalScaleFactors(4)
        for i in range(4):
            eftNV.setScaleFactorType(i + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
        eftNV.setFunctionNumberOfTerms(1, 4)  # 4 terms = 4 cubic basis functions
        eftNV.setTermNodeParameter(1, 1, 3, Node.VALUE_LABEL_VALUE, 1)
        eftNV.setTermScaling(1, 1, [1])
        eftNV.setTermNodeParameter(1, 2, 3, Node.VALUE_LABEL_D_DS1, 1)
        eftNV.setTermScaling(1, 2, [2])
        eftNV.setTermNodeParameter(1, 3, 4, Node.VALUE_LABEL_VALUE, 1)
        eftNV.setTermScaling(1, 3, [3])
        eftNV.setTermNodeParameter(1, 4, 4, Node.VALUE_LABEL_D_DS1, 1)
        eftNV.setTermScaling(1, 4, [4])
        elementtemplateBranchRoot = mesh.createElementtemplate()
        elementtemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
        elementtemplateBranchRoot.defineField(coordinates, -1, eftNV)

        # branch special node and element field template - material coordinates
        vagusNodeTemplateNoValue = nodes.createNodetemplate()
        vagusNodeTemplateNoValue.defineField(vagusCoordinates)
        vagusNodeTemplateNoValue.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_VALUE, 0)
        vagusNodeTemplateNoValue.setValueNumberOfVersions(vagusCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        eftNVvagus = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eftNVvagus.setNumberOfLocalNodes(4)
        eftNVvagus.setNumberOfLocalScaleFactors(4)
        for i in range(4):
            eftNVvagus.setScaleFactorType(i + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
        eftNVvagus.setFunctionNumberOfTerms(1, 4)  # 4 terms = 4 cubic basis functions
        eftNVvagus.setTermNodeParameter(1, 1, 3, Node.VALUE_LABEL_VALUE, 1)
        eftNVvagus.setTermScaling(1, 1, [1])
        eftNVvagus.setTermNodeParameter(1, 2, 3, Node.VALUE_LABEL_D_DS1, 1)
        eftNVvagus.setTermScaling(1, 2, [2])
        eftNVvagus.setTermNodeParameter(1, 3, 4, Node.VALUE_LABEL_VALUE, 1)
        eftNVvagus.setTermScaling(1, 3, [3])
        eftNVvagus.setTermNodeParameter(1, 4, 4, Node.VALUE_LABEL_D_DS1, 1)
        eftNVvagus.setTermScaling(1, 4, [4])
        vagusElementtemplateBranchRoot = mesh.createElementtemplate()
        vagusElementtemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
        vagusElementtemplateBranchRoot.defineField(vagusCoordinates, -1, eftNVvagus)

        # build initial 1d trunk line
        elementsAlongTrunk = options['Number of elements along the trunk']
        step = rescaled_lengthToDiameterRatio / (elementsAlongTrunk - 1)
        trunk_nodes, trunk_ld1, vagus_trunk_nodes, vagus_trunk_ld1 = \
            evaluate_trunk_coordinates(elementsAlongTrunk, step, use_marker_names, marker_data,
                                       rescaledmarker_TermNameVagusLengthList)

        # annotations
        annotationGroups = []
        vagusTrunkGroup = AnnotationGroup(region, (trunk_group_name, 'None'))
        annotationGroups.append(vagusTrunkGroup)

        nodeIdentifier = 1
        elementIdentifier = 1

        vagusNodeIdentifier = 1
        vagusElementIdentifier = 1

        # build trunk
        for n in range(elementsAlongTrunk):
            # geometric coordinates
            lx = trunk_nodes[n]
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, trunk_ld1)
            if n > 0:
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [nodeIdentifier - 1, nodeIdentifier])
                meshGroup = vagusTrunkGroup.getMeshGroup(mesh)
                meshGroup.addElement(element)
                elementIdentifier += 1
            nodeIdentifier += 1

            # material coordinates:
            # lx = vagus_trunk_nodes[n]
            # node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
            # node.merge(vagusNodetemplate)
            # fieldcache.setNode(node)
            # vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
            # vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, vagus_trunk_ld1)
            # if n > 0:
            #     element = mesh.findElementByIdentifier(vagusElementIdentifier)
            #     element.merge(vagusElementtemplate)
            #     element.setNodesByIdentifier(eftVagus, [vagusNodeIdentifier - 1, vagusNodeIdentifier])
            #     vagusElementIdentifier += 1
            # vagusNodeIdentifier += 1

        # set markers
        for marker_name, marker_coordinate in marker_data.items():
            annotationGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
            annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
            nodeIdentifier += 1
            vagusNodeIdentifier += 1

        # geometry fitting - trunk
        if applyFitting:
            # create temporary model file
            sir = region.createStreaminformationRegion()
            srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf")
            region.write(sir)

            fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
            fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf"
            fitter = fit_trunk_model(fitter_model_file, fitter_data_file)
            reset_nodes_with_fitted_nodes(region, fitter)

        trunk_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        trunk_fitted_data = extract_values_for_group_by_label(region, "coordinates", trunk_group_name, Node.VALUE_LABEL_VALUE)

        # building initial branches
        if addBranches:
            for branch_name in branch_data.keys():
                print(branch_name)

                branchGroup = AnnotationGroup(region, (branch_name, 'None'))
                annotationGroups.append(branchGroup)

                # determine branch approximate start and closest trunk node index
                branch_coordinates = branch_data[branch_name]
                min_distance_squared = float('inf')
                for i in range(len(trunk_fitted_data)):
                    trunk_point = trunk_fitted_data[i]
                    for branch_point in branch_coordinates:
                        distance_squared = magnitude_squared(sub(branch_point, trunk_point))
                        if distance_squared <= min_distance_squared:
                            min_distance_squared = distance_squared
                            branch_start_coordinate = branch_point
                            trunk_index = i
                # determine trunk segment where branch starts (previous or next index)
                if trunk_index == len(trunk_fitted_data) - 1:
                    trunk_index -= 1
                if trunk_index > 0 and (trunk_index < len(trunk_fitted_data) - 1):
                    distance_sq_previous_node = magnitude_squared(sub(branch_start_coordinate, trunk_fitted_data[trunk_index - 1]))
                    distance_sq_next_node = magnitude_squared(sub(branch_start_coordinate, trunk_fitted_data[trunk_index + 1]))
                    if distance_sq_previous_node < distance_sq_next_node:
                        trunk_index -= 1

                print('\t branch start coordinate ' + str(branch_start_coordinate))

                # annotationGroup = findOrCreateAnnotationGroupForTerm(
                #     annotationGroups, region, ('start of ' + branch_name, None), isMarker=True)
                # annotationGroup.createMarkerNode(nodeIdentifier, coordinates, branch_start_coordinate)
                # nodeIdentifier += 1
                # vagusNodeIdentifier += 1

                # get coordinates and derivatives of trunk nodes in the segment
                trunk_node_id = trunk_index + 1
                print('\t trunk segment %d - %d' % (trunk_node_id, trunk_node_id + 1))
                node = nodes.findNodeByIdentifier(trunk_node_id)
                fieldcache.setNode(node)
                _, v1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                _, v_v1 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, v_d1 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                node = nodes.findNodeByIdentifier(trunk_node_id + 1)
                fieldcache.setNode(node)
                _, v2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                _, v_v2 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)

                # find xi closest to branch_start on a cubic Hermite curve by bisection
                xi_a = 0
                xi_b = 1
                min_distance_sq = float('inf')
                while (xi_b - xi_a) > 0.01:
                    d_sq_a = magnitude_squared(sub(branch_start_coordinate, interpolateCubicHermite(v1, d1, v2, d2, xi_a)))
                    d_sq_b = magnitude_squared(sub(branch_start_coordinate, interpolateCubicHermite(v1, d1, v2, d2, xi_b)))
                    #print('\txi_a=%s \txi_b=%s \td_sq_a=%s \td_sq_b=%s' % (xi_a, xi_b, d_sq_a, d_sq_b))
                    if d_sq_a >= d_sq_b:
                        xi_a = (xi_a + xi_b) / 2
                    else:
                        xi_b = (xi_a + xi_b) / 2

                branch_root_xi = (xi_a + xi_b) / 2
                print('\t xi = %s' % branch_root_xi)

                # determine branch approximate end
                max_distance_squared = 0
                for branch_node in branch_coordinates:
                    distance_squared = magnitude_squared(sub(branch_node, branch_start_coordinate))
                    if distance_squared >= max_distance_squared:
                        max_distance_squared = distance_squared
                        branch_end_coordinate = branch_node

                # calculate branch coordinates
                branch_start_coordinate = interpolateHermiteLagrange(v1, d1, v2, branch_root_xi)
                vagus_branch_start_coordinate = interpolateHermiteLagrange(v_v1, v_d1, v_v2, branch_root_xi)
                branch_length = magnitude(sub(branch_start_coordinate, branch_end_coordinate))
                elementsAlongBranch = int(branch_length/step)
                if elementsAlongBranch < 3:
                    elementsAlongBranch = 3
                if elementsAlongBranch > 7:
                    elementsAlongBranch = 7
                branch_coordinates = []
                vagus_branch_coordinates = []
                dx, dy, dz = div(sub(branch_end_coordinate, branch_start_coordinate), (elementsAlongBranch - 1))

                for i in range(elementsAlongBranch):
                    branch_coordinates.append([branch_start_coordinate[0] + dx * i,
                                               branch_start_coordinate[1] + dy * i,
                                               branch_start_coordinate[2] + dz * i
                                               ])
                    vagus_branch_coordinates.append([vagus_branch_start_coordinate[0] + dx * i,
                                                     vagus_branch_start_coordinate[1] + dy * i,
                                                     vagus_branch_start_coordinate[2] - dz * i
                                                     ])


                # build branches
                for n in range(elementsAlongBranch):
                    # geometric coordinates
                    branch_ld1 = [dx, dy, dz]
                    if n == 0:
                        # create branch special node
                        node = nodes.createNode(nodeIdentifier, nodeTemplateNoValue)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, branch_ld1)
                    else:
                        lx = branch_coordinates[n]
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, branch_ld1)
                        if n == 1:
                            # create branch special element
                            element = mesh.createElement(elementIdentifier, elementtemplateBranchRoot)
                            element.setNodesByIdentifier(eftNV, [nodeIdentifier - 1, nodeIdentifier,
                                                                 trunk_node_id, trunk_node_id + 1])
                            scalefactors = getCubicHermiteBasis(branch_root_xi)
                            element.setScaleFactors(eftNV, list(scalefactors))
                            meshGroup = branchGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)
                            elementIdentifier += 1
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            element.setNodesByIdentifier(eft, [nodeIdentifier - 1, nodeIdentifier])
                            meshGroup = branchGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)
                            elementIdentifier += 1
                    nodeIdentifier += 1

                branchNodesetGroup = branchGroup.getNodesetGroup(nodes)
                if branchNodesetGroup.isValid():
                    branchNodesetGroup.removeNodesConditional(trunk_group)

                    # # material coordinates
                    # branch_ld1 = [dx, dy, -dz]
                    # if n == 0:
                    #     # branch special node in vagus coordinates
                    #     node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
                    #     node.merge(vagusNodeTemplateNoValue)
                    #     fieldcache.setNode(node)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, branch_ld1)
                    # else:
                    #     lx = vagus_branch_coordinates[n]
                    #     node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
                    #     node.merge(vagusNodetemplate)
                    #     fieldcache.setNode(node)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, branch_ld1)
                    #     if n == 1:
                    #         # create branch special element
                    #         element = mesh.findElementByIdentifier(vagusElementIdentifier)
                    #         element.merge(vagusElementtemplateBranchRoot)
                    #         element.setNodesByIdentifier(eftNVvagus, [vagusNodeIdentifier - 1, vagusNodeIdentifier,
                    #                                                   trunk_node_id, trunk_node_id + 1])
                    #         scalefactors = getCubicHermiteBasis(branch_root_xi)
                    #         element.setScaleFactors(eftNVvagus, list(scalefactors))
                    #         vagusElementIdentifier += 1
                    #     else:
                    #         element = mesh.findElementByIdentifier(vagusElementIdentifier)
                    #         element.merge(vagusElementtemplate)
                    #         element.setNodesByIdentifier(eftVagus, [vagusNodeIdentifier - 1, vagusNodeIdentifier])
                    #         vagusElementIdentifier += 1
                    # vagusNodeIdentifier += 1

            sir = region.createStreaminformationRegion()
            srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_model.exf")
            region.write(sir)

            # geometry fitting - branches
            if applyFitting:
                fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
                fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_model.exf"

                for branch_name in branch_data.keys():
                    print('fitting %s' % branch_name)
                    fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
                    reset_nodes_with_fitted_nodes(region, fitter, branch_name)

        sir = region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_final_model.exf")
        region.write(sir)

        return annotationGroups, None


def evaluate_trunk_coordinates(elementsAlongTrunk, step, use_marker_names, marker_data, termNameVagusLengthList):
    """

    """
    trunk_nodes = []
    vagus_trunk_nodes = []
    x1, y1, z1 = marker_data[use_marker_names[0]]
    x2, y2, z2 = marker_data[use_marker_names[1]]
    t1 = termNameVagusLengthList[use_marker_names[0]]
    t2 = termNameVagusLengthList[use_marker_names[1]]
    dx, dy, dz = [(x2 - x1) / (t2 - t1), (y2 - y1) / (t2 - t1), (z2 - z1) / (t2 - t1)]
    trunk_ld1 = [dx * step, dy * step, dz * step]
    vagus_trunk_ld1 = [0, 0, step]
    for i in range(elementsAlongTrunk):
        trunk_nodes.append([
            x1 + dx * (i * step - t1),
            y1 + dy * (i * step - t1),
            z1 + dz * (i * step - t1)])
        vagus_trunk_nodes.append([0, 0, i * step])

    return trunk_nodes, trunk_ld1, vagus_trunk_nodes, vagus_trunk_ld1


def extract_values_for_group_by_label(region, field_name, group_name, valueLabel):
    """

    """

    fieldmodule = region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fieldmodule.findFieldByName(field_name).castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)

    coordinates_list = []
    group = fieldmodule.findFieldByName(group_name).castGroup()
    if group.isValid():
        group_nodes = group.getNodesetGroup(nodes)
        node_iter = group_nodes.createNodeiterator()
        node = node_iter.next()
        while node.isValid():
            fieldcache.setNode(node)
            result, x = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1, 3)
            coordinates_list.append(x)
            node = node_iter.next()
    return coordinates_list


def reset_nodes_with_fitted_nodes(region, fitter, group_name = None):
    """
    """

    fieldmodule = region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()
    coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    fitter_fieldmodule = fitter.getFieldmodule()
    fitter_fieldcache = fitter_fieldmodule.createFieldcache()
    fitter_coordinates = fitter.getModelCoordinatesField().castFiniteElement()
    fitter_nodes = fitter_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    if group_name:
        group = find_or_create_field_group(fitter_fieldmodule, group_name)
        if group.isValid():
            fitter_nodes = group.getNodesetGroup(fitter_nodes)

    # reset trunk nodes with the fitted nodes
    fitter_node_iter = fitter_nodes.createNodeiterator()
    fitter_node = fitter_node_iter.next()
    while fitter_node.isValid():
        fitter_fieldcache.setNode(fitter_node)
        _, lx = fitter_coordinates.getNodeParameters(fitter_fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        _, ld1 = fitter_coordinates.getNodeParameters(fitter_fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

        node = nodes.findNodeByIdentifier(fitter_node.getIdentifier())
        fieldcache.setNode(node)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
        fitter_node = fitter_node_iter.next()


def fit_trunk_model(modelfile, datafile, trunk_group_name = None):
    """

    """
    fitter = Fitter(modelfile, datafile)
    fitter.load()

    # initial configuration
    fitter_fieldmodule = fitter.getFieldmodule()
    fitter.setModelCoordinatesFieldByName('coordinates')
    fitter.setDataCoordinatesFieldByName('coordinates')
    if trunk_group_name:
        fitter.setModelFitGroupByName(trunk_group_name)
    fitter.setFibreField(fitter_fieldmodule.findFieldByName("zero fibres"))
    fitter.setMarkerGroupByName('marker')  # not necessary, it's marker by default
    fitter.setDiagnosticLevel(0)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupDataWeight('marker', 200.0)
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight(None, 5.0)
    fit1.setNumberOfIterations(10)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    fit2 = FitterStepFit()
    fit2.setGroupDataWeight('marker', 400.0)
    fit1.setGroupDataWeight(None, None)
    fit2.setNumberOfIterations(5)
    fit2.setUpdateReferenceState(True)
    fitter.addFitterStep(fit2)


    # align step
    # align = FitterStepAlign()
    # align.setAlignMarkers(True)
    # align.setAlignGroups(True)
    # align.setScaleProportion(1.0)
    # fitter.addFitterStep(align)
    # align.run()

    # fit step 1
    # fit1 = FitterStepFit()
    # fit1.setGroupStrainPenalty(None, [15.0])
    # fit1.setGroupCurvaturePenalty(None, [50.0])
    # fit1.setGroupDataWeight('marker', 200.0)
    # fit1.setNumberOfIterations(5)
    # fit1.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit1)
    #
    # # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionError()
    rmsTrunkError, maxTrunkError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('left vagus X nerve trunk')
    rmsMarkerError, maxMarkerError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('marker')

    print('(all) RMS error: ' + str(rmsError))
    print('(all) Max error: ' + str(maxError))
    print('(trunk) RMS error: ' + str(rmsTrunkError))
    print('(trunk) Max error: ' + str(maxTrunkError))
    print('(marker) RMS error: ' + str(rmsMarkerError))
    print('(marker) Max error: ' + str(maxMarkerError))

    # fitter_nodes = fitter_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    # fitter_coordinates = fitter.getModelCoordinatesField().castFiniteElement()
    # valueLabels, fieldParameters = get_nodeset_field_parameters(fitter_nodes, fitter_coordinates)
    # print_node_field_parameters(valueLabels, fieldParameters, '{: 1.4f}')

    return fitter


def fit_branches_model(modelfile, datafile, branch_name = None):
    """

    """
    fitter = Fitter(modelfile, datafile)
    fitter.load()

    # initial configuration
    fitter_fieldmodule = fitter.getFieldmodule()
    fitter.setModelCoordinatesFieldByName('coordinates')
    if branch_name:
        fitter.setModelFitGroupByName(branch_name)
    fitter.setFibreField(fitter_fieldmodule.findFieldByName("zero fibres"))
    fitter.setDataCoordinatesFieldByName('coordinates')
    fitter.setMarkerGroupByName('marker')  # not necessary, it's marker by default
    fitter.setDiagnosticLevel(0)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight('marker', 200.0)
    fit1.setNumberOfIterations(5)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionError()
    print('(branches) RMS error: ' + str(rmsError))
    print('(branches) Max error: ' + str(maxError))

    # fitter_nodes = fitter_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    # fitter_coordinates = fitter.getModelCoordinatesField().castFiniteElement()
    # valueLabels, fieldParameters = get_nodeset_field_parameters(fitter_nodes, fitter_coordinates)
    # print_node_field_parameters(valueLabels, fieldParameters, '{: 1.4f}')

    return fitter


def load_data(region):
    """
    Extract data from supplied datafile, separate out data related to
    vagus trunk, vagus branches, fascicles, markers (anatomical landmarks)
    """

    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid()
    fieldmodule = data_region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()

    coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)

    markers = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    marker_coordinates = fieldmodule.findFieldByName("marker_data_coordinates").castFiniteElement()
    marker_names = fieldmodule.findFieldByName("marker_data_name")
    assert marker_coordinates.isValid() and (marker_coordinates.getNumberOfComponents() == 3)

    group_list = get_group_list(fieldmodule)
    group_map = {}
    trunk_group_name = None
    branch_groups = []
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group
        if 'trunk' in group_name:
            trunk_group_name = group_name
        if 'branch' and ('left A' or 'right A') in group_name:
            branch_groups.append(group_name)

    # extract data related to vagus markers, trunk and A-branches
    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fieldcache.setNode(marker_node)
            result, x = marker_coordinates.evaluateReal(fieldcache, 3)
            marker_name = marker_names.evaluateString(fieldcache)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    trunk_data_coordinates = extract_values_for_group_by_label(data_region, "coordinates", trunk_group_name, Node.VALUE_LABEL_VALUE)
    branch_data = {}
    for branch_name in branch_groups:
        branch_data_coordinates = extract_values_for_group_by_label(data_region, "coordinates", branch_name, Node.VALUE_LABEL_VALUE)
        branch_data[branch_name] = branch_data_coordinates

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    return trunk_group_name, trunk_data_coordinates, marker_data, branch_data


def load_data_xml(region):
    """
    Extract data from supplied datafile, separate out data related to
    vagus trunk, vagus branches, fascicles, markers (anatomical landmarks)
    """

    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid()
    fieldmodule = data_region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()

    coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)

    markers = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    #marker_coordinates = fieldmodule.findFieldByName("marker_data_coordinates").castFiniteElement()
    marker_names = fieldmodule.findFieldByName("marker_name")
    #print(marker_coordinates.getNumberOfComponents())
    #assert marker_coordinates.isValid() and (marker_coordinates.getNumberOfComponents() == 3)

    group_list = get_group_list(fieldmodule)
    group_map = {}
    trunk_group_name = None
    branch_groups = []
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group
        if 'trunk' in group_name:
            trunk_group_name = group_name
        if 'branch' and ('left A' or 'right A') in group_name:
            branch_groups.append(group_name)

    # extract data related to vagus markers, trunk and A-branches
    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fieldcache.setNode(marker_node)
            result, x = coordinates.evaluateReal(fieldcache, 3)
            marker_name = marker_names.evaluateString(fieldcache)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    trunk_data_coordinates = extract_values_for_group_by_label(data_region, "coordinates", trunk_group_name, Node.VALUE_LABEL_VALUE)
    branch_data = {}
    for branch_name in branch_groups:
        branch_data_coordinates = extract_values_for_group_by_label(data_region, "coordinates", branch_name, Node.VALUE_LABEL_VALUE)
        branch_data[branch_name] = branch_data_coordinates

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    print(marker_data)

    return trunk_group_name, trunk_data_coordinates, marker_data, branch_data


