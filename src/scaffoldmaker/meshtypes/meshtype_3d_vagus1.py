"""
Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches.
"""

import math

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, get_group_list, find_or_create_field_group, \
    create_field_group, get_managed_field_names
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field, FieldGroup
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    findAnnotationGroupByName
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from cmlibs.maths.vectorops import add, sub, mult, div, dot
from scaffoldmaker.utils.vector import magnitude_squared, magnitude
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, interpolateCubicHermite, interpolateLagrangeHermite, \
    interpolateHermiteLagrange, sampleCubicHermiteCurves
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters, print_node_field_parameters
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelWithNodes
from scaffoldfitter.fitter import Fitter
from scaffoldfitter.fitterstepconfig import FitterStepConfig
from scaffoldfitter.fitterstepalign import FitterStepAlign
from scaffoldfitter.fitterstepfit import FitterStepFit

from scaffoldmaker.utils.read_vagus_data import load_exf_data


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
            'Human Trunk 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk': 30,
            'Iterations (fit trunk)': 1,
            'Apply fitting': False,
            'Add branches': True
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along the trunk',
            'Iterations (fit trunk)',
            'Apply fitting',
            'Add branches'
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options['Number of elements along the trunk'] < 10:
            options['Number of elements along the trunk'] = 10
        if options['Iterations (fit trunk)'] < 1:
            options['Iterations (fit trunk)'] = 1
        return dependentChanges


    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        # setup
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]

        # geometric coordinates
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for value_label in value_labels:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, value_label, 1)

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
        for value_label in value_labels:
            vagusNodetemplate.setValueNumberOfVersions(vagusCoordinates, -1, value_label, 1)

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

        elementsAlongTrunk = options['Number of elements along the trunk']
        iterationsNumber = options['Iterations (fit trunk)']
        applyFitting = options['Apply fitting']
        addBranches = options['Add branches']

        # load data from file
        data_region = region.getParent().findChildByName('data')
        if data_region.isValid():
            marker_data, trunk_group_name, trunk_data, _, branch_data, branch_parents, _ = load_exf_data(data_region)
        assert len(marker_data) >= 2, f"At least two landmarks are expected in the data. Incomplete data."

        # field group used for fitting
        trunkFitCentroidGroup = find_or_create_field_group(fieldmodule, trunk_group_name + '-fit')
        trunkFitCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
        trunkFitCentroidMeshGroup = trunkFitCentroidGroup.getOrCreateMeshGroup(mesh)

        for ii in range(iterationsNumber):
            # annotations
            annotationGroups = []
            vagusTrunkGroup = AnnotationGroup(region, (trunk_group_name, 'None'))
            annotationGroups.append(vagusTrunkGroup)
            vagusTrunkMeshGroup = vagusTrunkGroup.getMeshGroup(mesh)

            if ii == 0:
                tx, td1, vx, vd1, elementLength = estimate_trunk_coordinates(elementsAlongTrunk, marker_data)
            else:
                # read tx from fit_coordinates
                _, node_field_parameters = get_nodeset_field_parameters(nodes, coordinates, value_labels)
                tx = [nodeParameter[1][0][0] for nodeParameter in node_field_parameters]
                td1 = [nodeParameter[1][1][0] for nodeParameter in node_field_parameters]

            trunk_nodes_data_bounds = estimate_trunk_data_boundaries(tx, trunk_data, elementsAlongTrunk)

            nodeIdentifier = 1
            elementIdentifier = 1
            vagusNodeIdentifier = 1
            vagusElementIdentifier = 1

            nodes_before = []
            nodes_after = []
            for n in range(elementsAlongTrunk):
                # geometric coordinates
                lx = tx[n]
                ld1 = td1[n]

                if ii == 0:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                else:
                    node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)

                # add to trunk group used for data fitting
                if trunk_nodes_data_bounds[0] <= nodeIdentifier <= trunk_nodes_data_bounds[-1]:
                    pass
                elif nodeIdentifier < trunk_nodes_data_bounds[0]:
                    nodes_before.append(nodeIdentifier)
                else:
                    nodes_after.append(nodeIdentifier)

                if n > 0:
                    nids = [nodeIdentifier - 1, nodeIdentifier]
                    if ii == 0:
                        line = mesh.createElement(elementIdentifier, elementtemplate)
                    else:
                        line = mesh.findElementByIdentifier(elementIdentifier)
                    line.setNodesByIdentifier(eft, nids)
                    vagusTrunkMeshGroup.addElement(line)
                    # add element to trunk group used for data fitting
                    if nodeIdentifier - 1 >= trunk_nodes_data_bounds[0] and nodeIdentifier <= trunk_nodes_data_bounds[-1]:
                        trunkFitCentroidMeshGroup.addElement(line)
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

            if ii == 0:
                # set markers
                for marker_name, marker_coordinate in marker_data.items():
                    #annotationGroup = findOrCreateAnnotationGroupForTerm(
                    #    annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
                    annotationGroup = findOrCreateAnnotationGroupForTerm(
                        annotationGroups, region, (marker_name, ''), isMarker=True)
                    annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
                    nodeIdentifier += 1
                    vagusNodeIdentifier += 1
            else:
                nodeIdentifier += len(marker_data)
                vagusNodeIdentifier += len(marker_data)

        # geometry fitting - trunk
        if applyFitting:
            # create temporary model file
            sir = region.createStreaminformationRegion()
            srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf")
            region.write(sir)

            print('... Fitting trunk, iteration', str(ii + 1))
            fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
            fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf"
            fitter = fit_trunk_model(fitter_model_file, fitter_data_file, trunk_group_name + '-fit')
            set_fitted_group_nodes(region, fitter, trunk_group_name + '-fit')

            if len(nodes_before) > 0:
                # recalculate unfitted nodes by the first fitted node
                node = nodes.findNodeByIdentifier(nodes_before[-1] + 1)
                fieldcache.setNode(node)
                _, lx = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, ld1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

                node_count = 1
                for i in range(len(nodes_before) - 1, -1, -1):
                    node_id = nodes_before[i]
                    x = [lx[j] - node_count * ld1[j] for j in range(3)]

                    node = nodes.findNodeByIdentifier(node_id)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                    node_count += 1

            if len(nodes_after) > 0:
                # recalculate unfitted nodes by the last fitted node
                node = nodes.findNodeByIdentifier(nodes_after[0] - 1)
                fieldcache.setNode(node)
                _, lx = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, ld1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

                node_count = 1
                for i in range(len(nodes_after)):
                    node_id = nodes_after[i]
                    x = [lx[j] + node_count * ld1[j] for j in range(3)]

                    node = nodes.findNodeByIdentifier(node_id)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                    node_count += 1

        trunk_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        #trunk_fitted_data = extract_values_for_group_by_label(region, "coordinates", trunk_group_name, Node.VALUE_LABEL_VALUE)


        # building initial branches
        if addBranches:
            visited_branches = []

            print('... Adding branches')
            queue = [branch for branch in branch_parents.keys() if branch_parents[branch] == trunk_group_name]
            while queue:
                branch_name = queue.pop(0)
                # print(branch_name)

                if branch_name in visited_branches:
                    continue
                visited_branches.append(branch_name)
                queue.extend([branch for branch in branch_parents.keys() if branch_parents[branch] == branch_name])

                branchGroup = AnnotationGroup(region, (branch_name, 'None'))
                annotationGroups.append(branchGroup)
                branchMeshGroup = branchGroup.getMeshGroup(mesh)

                branch_coordinates = [branch_node[0] for branch_node in branch_data[branch_name]]
                branch_parent_name = branch_parents[branch_name]
                #print(branch_name, ' -> ', branch_parent_name)

                # determine branch approximate start and closest trunk node index
                bx, vbx, bd1, parent_s_nid, parent_f_nid, branch_root_xi, elementsAlongBranch = \
                    estimate_branch_coordinates(region, branch_coordinates, elementLength, branch_parent_name)
                #print('  branch between nodes: ', parent_s_nid, parent_f_nid)

                for n in range(elementsAlongBranch):
                    sx = bx[n]
                    sd1 = bd1

                    if n == 0:
                        # create branch special node
                        node = nodes.createNode(nodeIdentifier, nodeTemplateNoValue)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

                        if n == 1:
                            # create branch root element
                            nids = [nodeIdentifier - 1, nodeIdentifier,
                                    parent_s_nid, parent_f_nid]
                            element = mesh.createElement(elementIdentifier, elementtemplateBranchRoot)
                            element.setNodesByIdentifier(eftNV, nids)
                            scalefactorsNV = getCubicHermiteBasis(branch_root_xi)
                            element.setScaleFactors(eftNV, list(scalefactorsNV))
                            branchMeshGroup.addElement(element)
                            elementIdentifier += 1
                        else:
                            nids = [nodeIdentifier - 1, nodeIdentifier]
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            element.setNodesByIdentifier(eft, nids)
                            branchMeshGroup.addElement(element)
                            elementIdentifier += 1
                    nodeIdentifier += 1

                    # material coordinates
                    # sx = vbx[n]
                    # sd1 = [bd1[0], bd1[1], -bd1[2]]
                    # if n == 0:
                    #     # branch special node in vagus coordinates
                    #     node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
                    #     node.merge(vagusNodeTemplateNoValue)
                    #     fieldcache.setNode(node)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
                    # else:
                    #
                    #     node = nodes.findNodeByIdentifier(vagusNodeIdentifier)
                    #     node.merge(vagusNodetemplate)
                    #     fieldcache.setNode(node)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                    #     vagusCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
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

                # remove trunk nodes from branch group
                parentGroup = find_or_create_field_group(fieldmodule, branch_parent_name)
                branchNodesetGroup = branchGroup.getNodesetGroup(nodes)
                if branchNodesetGroup.isValid():
                    branchNodesetGroup.removeNodesConditional(parentGroup)

                sir = region.createStreaminformationRegion()
                srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_model.exf")
                region.write(sir)

                # geometry fitting - branches
                if applyFitting:
                    fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
                    fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_model.exf"

                    #print('fitting %s' % branch_name)
                    fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
                    set_fitted_group_nodes(region, fitter, branch_name)

        sir = region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_final_model.exf")
        region.write(sir)

        return annotationGroups, None


def estimate_trunk_coordinates(elementsAlongTrunk, marker_data):
    """

    """

    # choose markers for building initial scaffold
    # at the moment uses the first and the last markers in the data
    termNameVagusLengthList = {
        "level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
        "level of superior border of jugular foramen": 8.6342,
        "level of inferior border of jugular foramen": 16.7227,
        "level of C1 transverse process": 32.1129,
        "level of angle of mandible": 42.2450,
        "level of greater horn of hyoid": 45.6122,
        "level of carotid bifurcation": 48.3581,
        "level of laryngeal prominence": 68.8431,
        "level of superior border of the clavicle": 117.5627,
        "level of jugular notch": 124.6407,
        "level of carina": 149.5929,  # not on the list of annotations yet!
        "level of sternal angle": 151.2352,
        "level of 1 cm superior to start of esophageal plexus": 165.5876,
        "level of esophageal hiatus": 254.32879,
        "level of aortic hiatus": 291.3695,
        "level of end of trunk": 312.5  # note this term is also not on the list of annotations
    }

    totalVagusLength = 312.5  # calculated from total length of nerve/average diameter of nerve

    use_markers = [list(marker_data.keys())[0],
                   list(marker_data.keys())[-1]]

    pts = []
    params = []

    for marker in use_markers:
        use_marker_name = marker.replace('left ', '', 1).replace('right ', '', 1).replace(' on the vagus nerve', '', 1)
        assert use_marker_name in termNameVagusLengthList

        pts.append(marker_data[marker])
        params.append(termNameVagusLengthList[use_marker_name])

    step = totalVagusLength / (elementsAlongTrunk - 1)
    dx, dy, dz = [(pts[1][dim] - pts[0][dim]) / (params[1] - params[0]) for dim in range(3)]

    trunk_nodes = []
    vagus_trunk_nodes = []
    trunk_d1 = []
    vagus_trunk_d1 = []
    for i in range(elementsAlongTrunk):
        trunk_nodes.append([pts[0][0] + dx * (i * step - params[0]),
                            pts[0][1] + dy * (i * step - params[0]),
                            pts[0][2] + dz * (i * step - params[0])])
        vagus_trunk_nodes.append([0, 0, i * step])
        trunk_d1.append([dx * step, dy * step, dz * step])
        vagus_trunk_d1.append([0, 0, step])


    return trunk_nodes, trunk_d1, vagus_trunk_nodes, vagus_trunk_d1, step


def estimate_trunk_data_boundaries(trunk_nodes, trunk_data, elementsAlongTrunk):
    """
    """

    # finding trunk nodes at the start and end of trunk data - by projections
    # TODO: change brute force to something (takes too long for CASE data),
    # at the moment assumes trunk_data is numbered top to bottom

    #trunk_data_endpoints = find_dataset_endpoints([trunk_pt[0] for trunk_pt in trunk_data])
    trunk_data_endpoints = [trunk_data[0][0], trunk_data[-1][0]]

    trunk_nodes_data_bounds = []
    for ind, endpoint in enumerate(trunk_data_endpoints):
        ap = sub(endpoint, trunk_nodes[0])
        ab = sub(trunk_nodes[-1], trunk_nodes[0])
        param = dot(ap, ab) / dot(ab, ab)  # between 0 and 1
        nearby_node = param * (elementsAlongTrunk - 1) + 1

        if ind == 0:
            # trunk node near the start of data
            trunk_nodes_data_bounds.append(math.floor(nearby_node))
        else:
            # trunk node near the end of data
            if nearby_node >= elementsAlongTrunk:
                trunk_nodes_data_bounds.append(elementsAlongTrunk)
            else:
                trunk_nodes_data_bounds.append(math.ceil(nearby_node))

    return trunk_nodes_data_bounds


def estimate_branch_coordinates(region, branch_coordinates, elementLength, branch_parent_name):
    """

    """

    fm = region.getFieldmodule()
    fieldcache = fm.createFieldcache()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    vagusCoordinates = fm.findFieldByName("vagus coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    # assumes parent_group_name is known
    branch_start_x, branch_end_x, parent_s_nid, parent_f_nid = find_branch_start_segment(
        region, branch_coordinates, branch_parent_name)

    # determine parent hermite curve parameters
    node = nodes.findNodeByIdentifier(parent_s_nid)
    fieldcache.setNode(node)
    _, px_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, pd1_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
    _, vx_1 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, vd1_1 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

    node = nodes.findNodeByIdentifier(parent_f_nid)
    fieldcache.setNode(node)
    _, px_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, pd1_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
    _, vx_2 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, vd1_2 = vagusCoordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

    # find xi closest to branch_start on a cubic Hermite curve by bisection
    xi_a = 0
    xi_b = 1
    eps = 0.005
    while (xi_b - xi_a) > eps:
        dsq_a = magnitude_squared(sub(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_a)))
        dsq_b = magnitude_squared(sub(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_b)))
        if dsq_a >= dsq_b:
            xi_a = (xi_a + xi_b) / 2
        else:
            xi_b = (xi_a + xi_b) / 2
    branch_root_xi = (xi_a + xi_b) / 2
    #print('\t xi = %s' % branch_root_xi)

    # recalculate branch start parameters
    branch_start_x = interpolateHermiteLagrange(px_1, pd1_1, px_2, branch_root_xi)
    vagus_branch_start_x = interpolateHermiteLagrange(vx_1, vd1_1, vx_2, branch_root_xi)

    branch_length = magnitude(sub(branch_end_x, branch_start_x))
    elementsAlongBranch = int(branch_length / elementLength - 1)
    if elementsAlongBranch < 4:
        elementsAlongBranch = 4
    if elementsAlongBranch > 10:
        elementsAlongBranch = 10
    #print('  branch_length: ', branch_length, ', # elements: ', elementsAlongBranch)

    branch_coordinates = []
    vagus_branch_coordinates = []
    dx, dy, dz = div(sub(branch_end_x, branch_start_x), (elementsAlongBranch - 1))
    for i in range(elementsAlongBranch):
        branch_coordinates.append([branch_start_x[0] + dx * i,
                                   branch_start_x[1] + dy * i,
                                   branch_start_x[2] + dz * i])
        vagus_branch_coordinates.append([vagus_branch_start_x[0] + dx * i,
                                         vagus_branch_start_x[1] + dy * i,
                                         vagus_branch_start_x[2] + dz * i])

    return branch_coordinates, vagus_branch_coordinates, [dx, dy, dz], parent_s_nid, parent_f_nid, \
           branch_root_xi, elementsAlongBranch


def find_branch_start_segment(region, branch_coordinates, parent_group_name):
    """

    """
    branch_ends_points = find_dataset_endpoints(branch_coordinates)

    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    parent_group = find_or_create_field_group(fm, parent_group_name)
    parent_nodeset = parent_group.getNodesetGroup(nodes)

    _, group_parameters = get_nodeset_field_parameters(parent_nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
    group_ids = [parameter[0] for parameter in group_parameters]
    group_x = [parameter[1][0][0] for parameter in group_parameters]

    # find closest to branch distance, branch start, group node closest to the branch
    min_dsq = float('inf')
    for i in range(len(group_x)):
        node_x = group_x[i]
        if node_x is None:
            continue
        for branch_point in branch_ends_points:
            distance_squared = magnitude_squared(sub(node_x, branch_point))
            if distance_squared <= min_dsq:
                min_dsq = distance_squared
                branch_start = branch_point
                closest_index = i

    # determine segment closest to branch  (previous or next to the node)
    if closest_index == len(group_x) - 1:
        closest_index -= 1
    if closest_index > 0 and (closest_index < len(group_x) - 1):
        # TODO: should think about the case when closest_index = 0:
        # actual starting point might not have a node,
        # therefore first segment is not always taken into consideration
        dsq_prev_node = magnitude_squared(sub(branch_start, group_x[closest_index - 1]))
        dsq_next_node = magnitude_squared(sub(branch_start, group_x[closest_index + 1]))
        if dsq_prev_node < dsq_next_node:
            closest_index -= 1

    parent_s_node_id = group_ids[closest_index]
    parent_f_node_id = group_ids[closest_index + 1]
    branch_end = branch_ends_points[0] if branch_ends_points[1] == branch_start else branch_ends_points[1]

    return branch_start, branch_end, parent_s_node_id, parent_f_node_id


def find_dataset_endpoints(coordinate_dataset):
    """
    Given list of node coordinates, find two furthest from each other nodes.
    Returns a list of two node coordinates.
    """

    max_distance_squared = 0
    ends_points = []
    for node_x in coordinate_dataset:
        for node_y in coordinate_dataset:
            distance_squared = magnitude_squared(sub(node_x, node_y))
            if distance_squared >= max_distance_squared:
                max_distance_squared = distance_squared
                ends_points = [node_y, node_x]
    return ends_points




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

    # align step
    # align = FitterStepAlign()
    # align.setAlignMarkers(True)
    # align.setAlignGroups(True)
    # align.setScaleProportion(1.0)
    # fitter.addFitterStep(align)
    # align.run()

    # fit step 1
    # fit1 = FitterStepFit()
    # fit1.setGroupDataWeight('marker', 200.0)
    # fit1.setGroupStrainPenalty(None, [15.0])
    # fit1.setGroupCurvaturePenalty(None, [50.0])
    # fit1.setGroupDataWeight(None, 5.0)
    # fit1.setNumberOfIterations(10)
    # fit1.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit1)
    #
    # # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit1.setGroupDataWeight(None, None)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight(None, 10.0)
    fit1.setNumberOfIterations(10)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    fit2 = FitterStepFit()
    fit1.setGroupDataWeight(None, [10.0])
    fit1.setGroupStrainPenalty(None, [5.0])
    fit1.setGroupCurvaturePenalty(None, [100.0])
    fit2.setNumberOfIterations(5)
    fit2.setUpdateReferenceState(True)
    fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionError()
    rmsTrunkError, maxTrunkError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('left vagus X nerve trunk')
    rmsMarkerError, maxMarkerError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('marker')

    # print('(all) RMS error: ' + str(rmsError))
    # print('(all) Max error: ' + str(maxError))
    # print('(trunk) RMS error: ' + str(rmsTrunkError))
    # print('(trunk) Max error: ' + str(maxTrunkError))
    # print('(marker) RMS error: ' + str(rmsMarkerError))
    # print('(marker) Max error: ' + str(maxMarkerError))

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
    fit1.setGroupDataWeight(None, 5.0)
    fit1.setNumberOfIterations(5)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionErrorForGroup(branch_name)
    #print('(%s) RMS error: %f' % (branch_name, rmsError))
    #print('(%s) Max error: %f' % (branch_name, maxError))

    # fitter_nodes = fitter_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    # fitter_coordinates = fitter.getModelCoordinatesField().castFiniteElement()
    # valueLabels, fieldParameters = get_nodeset_field_parameters(fitter_nodes, fitter_coordinates)
    # print_node_field_parameters(valueLabels, fieldParameters, '{: 1.4f}')

    return fitter


def set_fitted_group_nodes(region, fitter, group_name = None):
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

