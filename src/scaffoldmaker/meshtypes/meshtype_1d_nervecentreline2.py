"""
Generates a 1-D mesh for nerve centrelines.
"""

from __future__ import division

import json
import os
import sys

import numpy as np

from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldCoordinates,\
    findOrCreateFieldFiniteElement, findOrCreateFieldStoredString
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from opencmiss.zinc.element import Element, Elementbasis
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation import nerve_terms
from copy import deepcopy


class MeshType_1d_nervecentreline2(Scaffold_base):
    '''
    Generates a 1-D mesh for nerve centrelines
    '''
    @staticmethod
    def getName():
        return '1D Nerve Centreline 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Read centrelines from flatmap directory': False,
            'Manifest file': ''
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Read centrelines from flatmap directory',
            'Manifest file'
        ]

    @staticmethod
    def checkOptions(options):
        pass

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base cubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: AnnotationGroups list
        """
        read_files = options['Read centrelines from flatmap directory']
        input_file = options['Manifest file']
        manifest_file = r'C:\Users\egha355\Desktop\sparc3\codes\python_packages\rat-flatmap\manifest.json'
        working_dir = sys.argv[2]

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        radius = findOrCreateFieldFiniteElement(fm, "radius", components_count=1, managed=True)
        markerName = findOrCreateFieldStoredString(fm, name="marker_name", managed=True)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.defineField(radius)
        nodetemplate.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.defineField(markerName)

        nodetemplate1 = nodes.createNodetemplate()
        nodetemplate1.defineField(coordinates)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate1.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 3)
        nodetemplate1.defineField(radius)
        nodetemplate1.setValueNumberOfVersions(radius, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # if useCrossDerivatives:
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        #     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        linearBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft1.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 1)
        eft2 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft2.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 2)
        eft3 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft3.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, 3)
        eft32 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        eft32.setTermNodeParameter(2, 1, 2, Node.VALUE_LABEL_D_DS1, 3)
        eftRadius = mesh.createElementfieldtemplate(linearBasis)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate1.defineField(coordinates, -1, eft1)
        elementtemplate1.defineField(radius, -1, eftRadius)

        elementtemplate2 = mesh.createElementtemplate()
        elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate2.defineField(coordinates, -1, eft2)
        elementtemplate3 = mesh.createElementtemplate()
        elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate3.defineField(coordinates, -1, eft3)
        elementtemplate2.defineField(radius, -1, eftRadius)
        elementtemplate3.defineField(radius, -1, eftRadius)

        cache = fm.createFieldcache()

        # names_internal_anatomical_map = {
        #     "brain_44-1": "trigeminal nucleus",
        #     "brain_12-1": "nucleus of solitary tract",
        #     "brain_34-1": "dorsal motor nucleus of vagus nerve",
        #     "brain_33-1": "nucleus ambiguus",
        #     "point_1": "point_1",
        #     "ardell_4_branching_point": "ardell_4_branching_point",
        #     "ganglion_1-1": "superior vagus X ganglion",
        #     "bolser_X-1": "bolser_X-1",
        #     "ardell_3_branching_point": "ardell_3_branching_point",
        #     "point_3": "point_3",
        #     "label_1-1": "dura mater in the posterior cranial fossa",
        #     "ganglion_2-1": "inferior vagus X ganglion",
        #     "point_6": "point_6",
        #     "point_7": "point_7",
        #     "point_8": "point_8",
        #     "point_9": "point_9",
        #     "point_10": "point_10",
        #     "label_4-1": "tympanic membrane",
        #     "label_3-1": "floor of external acoustic meatus",
        #     "label_2-1": "posterior wall of external acoustic meatus",
        #     "point_11": "point_11",
        #     "point_13": "point_13",
        #     "label_6-1": "pharynx)",
        #     "label_5-1": "back part of the tongue*",
        #     "point_14": "point_14",
        #     "cardio_4-1": "carotid body",
        #     "point_16": "point_16",
        #     "point_18": "point_18",
        #     "point_19": "point_19",
        #     "point_20": "point_20",
        #     "label_7-1": "mucosa of pharynx",
        #     "label_8-1": "epiglottic vallecula",
        #     "respiratory_11-1": "epiglottis",
        #     "respiratory_5-1": "Wall of larynx",
        #     "point_21": "point_21",
        #     "label_10-1": "aryepiglottic fold",
        #     "label_11-1": "arytenoideus",
        #     "label_12-1": "mucosa of arytenoid cartilage",
        #     "point_37": "point_37",
        #     "bolser_X-5": "bolser_X-5",
        #     "bolser_X-4": "bolser_X-4",
        #     "external_laryngeal_n_branching_point": "external_laryngeal_n_branching_point",
        #     "plexus_2-1": "pharyngeal nerve plexus",
        #     "bolser_X-14": "bolser_X-14",
        #     "label_15-1": "inferior pharyngeal constrictor",
        #     "label_14-1": "cricothyroid muscle",
        #     "bolser_X-15": "bolser_X-15",
        #     "point_41": "point_41",
        #     "point_22": "point_22",
        #     "point_24": "point_24",
        #     "point_25": "point_25",
        #     "plexus_3-1": "cardiac nerve plexus",
        #     "point_26": "point_26",
        #     "point_28": "point_28",
        #     "label_9-1": "mucosa of larynx",
        #     "point_29": "point_29",
        #     "label_16-1": "ceratocricoid",
        #     "label_17-1": "lateral crico-arytenoid",
        #     "label_18-1": "oblique arytenoid",
        #     "label_19-1": "posterior crico-arytenoid",
        #     "label_20-1": "hyro-arytenoid",
        #     "label_21-1": "transverse arytenoid",
        #     "label_22-1": "vocalis muscle",
        #     "label_13-1": "vocal_cords",
        #     "respiratory_8-1": "Laryngeal mechanoreceptors",
        #     "point_38": "point_38",
        #     "bolser_X-2": "bolser_X-2",
        #     "bolser_X-3": "bolser_X-3",
        #     "point_30": "point_30",
        #     "plexus_4-1": "pulmonary nerve plexus",
        #     "point_31": "point_31",
        #     "respiratory_17-1": "bronchus smooth muscle",
        #     "point_32": "point_32",
        #     "point_33": "point_33",
        #     "plexus_5-1": "esophageal nerve plexus",
        #     "digestive_12-1": "esophagus",
        #     "point_34": "point_34",
        #     "digestive_9-1": "pancreas",
        #     "point_35": "point_35",
        #     "plexus_6-1": "left gastric nerve plexus",
        #     "point_36": "point_36",
        #     "urinary_13-1": "kidney",
        # }

        # markers group, coordinates and derivatives
        update_terms = False
        nerve_terms_changed = False
        # filename = os.path.join(working_dir, 'coordinates_modified.json')
        coordinates_file = os.path.join(working_dir, 'coordinates.json')
        properties_data, anatomicalMap_data = read_manifest(manifest_file)
        centrelines = get_centrelines_list(properties_data)
        points = get_points_in_centrelines(centrelines)

        # if update_terms:
        #     write_points_to_json(names_internal_anatomical_map, points, filename)
        nodes_data = read_json(coordinates_file)
        node_terms, centreline_terms = get_nerve_terms(points, centrelines, properties_data, anatomicalMap_data)
        if nerve_terms_changed:
            def print_terms(terms):
                terms_list = list(set([values for key, values in terms.items()]))
                terms_list.sort()
                for t in terms_list:
                    print(f"('{t[0]}', '{t[1]}'),")

            print_terms(node_terms)
            print_terms(centreline_terms)

        marker_data = get_marker_data(nodes_data, node_terms)
        centreline_data = get_centrline_data(centrelines, centreline_terms)

        nodes_w_ver = []

        annotationGroups = []
        marker_nodeid_map = {}
        # create nodes
        nodeIdentifier = 1

        skip_names = ['spinal_47-1']
        for name, marker in marker_data.items():
            # if name in skip_names:
            #     continue
            # if n3+1 in nodes_w_ver:
            #     nt = nodetemplate1
            # else:
            if 'L6' in name or 'T13' in name:
                continue
            nt = nodetemplate
            node = nodes.createNode(nodeIdentifier, nt)
            cache.setNode(node)

            marker_nodeid_map[name] = nodeIdentifier
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, marker['x'])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, marker["d1"])
            markerName.assignString(cache, name)
            annotationGroup = AnnotationGroup(region, marker["group"])
            annotationGroups.append(annotationGroup)
            annotationGroup.getNodesetGroup(nodes).addNode(node)
            # if n3+1 in nodes_w_ver:
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, nerveMarkerPoints[n3]["d1"])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, nerveMarkerPoints[n3]["d1"])

            radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, marker["r"])
            nodeIdentifier = nodeIdentifier + 1

        # connectivity
        # centrelines = [
        #     {'id': 'n_1', 'connects': ['brain_44-1', 'point_1']},
        #     {'id': 'n_5', 'connects': ['point_1', 'ganglion_1-1']},
        #     {'id': 'n_5', 'connects': ['ganglion_1-1', 'ganglion_2-1']},
        #     {'id': 'n_5', 'connects': ['ganglion_2-1', 'point_11']},
        #     {'id': 'n_5', 'connects': ['point_11', 'point_14']},
        #     {'id': 'n_5', 'connects': ['point_14', 'point_16']},
        #     {'id': 'n_5', 'connects': ['point_16', 'point_22']},
        #     {'id': 'n_5', 'connects': ['point_22', 'point_26']},
        #     {'id': 'n_5', 'connects': ['point_26', 'point_30']},
        #     {'id': 'n_5', 'connects': ['point_30', 'point_31']},
        #     {'id': 'n_5', 'connects': ['point_31', 'point_32']},
        #     {'id': 'n_5', 'connects': ['point_32', 'point_34']},
        #     {'id': 'n_5', 'connects': ['point_34', 'point_35']},
        #     {'id': 'n_5', 'connects': ['point_35', 'point_36']},
        #     {'id': 'n_68', 'connects': ['ganglion_1-1', 'point_3']},
        #     {'id': 'n_6', 'connects': ['ganglion_1-1', 'label_1-1']},
        #     {'id': 'n_7', 'connects': ['ganglion_2-1', 'point_6']},
        #     {'id': 'n_7', 'connects': ['point_6', 'point_8']},
        #     {'id': 'n_7', 'connects': ['point_8', 'point_10']},
        #     {'id': 'n_12', 'connects': ['point_6', 'point_7']},
        #     {'id': 'n_8', 'connects': ['point_8', 'point_9']},
        #     {'id': 'n_9', 'connects': ['point_10', 'label_3-1']},
        #     {'id': 'n_11', 'connects': ['point_10', 'label_4-1']},
        #     {'id': 'n_10', 'connects': ['point_10', 'label_2-1']},
        #     {'id': 'n_13', 'connects': ['point_11', 'point_13']},
        #     {'id': 'n_2', 'connects': ['point_1', 'brain_12-1']},
        #     {'id': 'n_15', 'connects': ['point_13', 'label_6-1']},
        #     {'id': 'n_14', 'connects': ['point_13', 'label_5-1']},
        #     {'id': 'n_16', 'connects': ['point_14', 'cardio_4-1']},
        #     {'id': 'n_3', 'connects': ['point_1', 'brain_34-1']},
        #     {'id': 'n_17', 'connects': ['point_16', 'point_18']},
        #     {'id': 'n_4', 'connects': ['point_1', 'brain_33-1']},
        #     {'id': 'n_18', 'connects': ['point_18', 'point_19']},
        #     {'id': 'n_29', 'connects': ['point_18', 'external_laryngeal_n_branching_point']},
        #     {'id': 'n_24', 'connects': ['point_19', 'point_21']},
        #     {'id': 'n_19', 'connects': ['point_19', 'point_20']},
        #     {'id': 'n_22', 'connects': ['point_20', 'label_7-1']},
        #     {'id': 'n_20', 'connects': ['point_20', 'label_8-1']},
        #     {'id': 'n_23', 'connects': ['point_20', 'respiratory_5-1']},
        #     {'id': 'n_21', 'connects': ['point_20', 'respiratory_11-1']},
        #     {'id': 'n_27', 'connects': ['point_21', 'label_12-1']},
        #     {'id': 'bolser_5', 'connects': ['point_21', 'bolser_X-5']},
        #     {'id': 'n_26', 'connects': ['point_21', 'label_11-1']},
        #     {'id': 'bolser_4', 'connects': ['point_21', 'bolser_X-4']},
        #     {'id': 'n_28', 'connects': ['point_21', 'point_37']},
        #     {'id': 'n_67', 'connects': ['external_laryngeal_n_branching_point', 'plexus_2-1']},
        #     {'id': 'n_68', 'connects': ['external_laryngeal_n_branching_point', 'point_41']},
        #     {'id': 'bolser_15', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-15']},
        #     {'id': 'bolser_14', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-14']},
        #     {'id': 'n_65', 'connects': ['external_laryngeal_n_branching_point', 'label_14-1']},
        #     {'id': 'n_66', 'connects': ['external_laryngeal_n_branching_point', 'label_15-1']},
        #     {'id': 'n_30', 'connects': ['point_22', 'point_24']},
        #     {'id': 'n_25', 'connects': ['point_21', 'label_10-1']},
        #     {'id': 'n_32', 'connects': ['point_24', 'point_25']},
        #     {'id': 'n_34', 'connects': ['point_25', 'plexus_3-1']},
        #     {'id': 'n_35', 'connects': ['point_24', 'plexus_3-1']},
        #     {'id': 'n_36', 'connects': ['point_26', 'point_28']},
        #     {'id': 'n_50', 'connects': ['point_28', 'point_38']},
        #     {'id': 'n_38', 'connects': ['point_28', 'point_29']},
        #     {'id': 'n_47', 'connects': ['point_28', 'label_13-1']},
        #     {'id': 'bolser_2', 'connects': ['point_28', 'bolser_X-2']},
        #     {'id': 'n_49', 'connects': ['point_28', 'respiratory_11-1']},
        #     {'id': 'n_48', 'connects': ['point_28', 'respiratory_8-1']},
        #     {'id': 'n_37', 'connects': ['point_28', 'label_9-1']},
        #     {'id': 'bolser_3', 'connects': ['point_28', 'bolser_X-3']},
        #     {'id': 'n_42', 'connects': ['point_29', 'label_18-1']},
        #     {'id': 'n_43', 'connects': ['point_29', 'label_19-1']},
        #     {'id': 'n_40', 'connects': ['point_29', 'label_16-1']},
        #     {'id': 'n_46', 'connects': ['point_29', 'label_22-1']},
        #     {'id': 'n_44', 'connects': ['point_29', 'label_20-1']},
        #     {'id': 'n_45', 'connects': ['point_29', 'label_21-1']},
        #     {'id': 'n_41', 'connects': ['point_29', 'label_17-1']},
        #     {'id': 'bolser_1', 'connects': ['ganglion_1-1', 'bolser_X-1']},
        #     {'id': 'ardell_3', 'connects': ['ganglion_1-1', 'ardell_3_branching_point']},
        #     {'id': 'ardell_4', 'connects': ['point_1', 'ardell_4_branching_point']},
        #     {'id': 'n_51', 'connects': ['point_30', 'plexus_4-1']},
        #     {'id': 'n_52', 'connects': ['point_31', 'respiratory_17-1']},
        #     {'id': 'n_53', 'connects': ['point_32', 'point_33']},
        #     {'id': 'n_54', 'connects': ['point_33', 'digestive_12-1']},
        #     {'id': 'n_55', 'connects': ['point_33', 'plexus_5-1']},
        #     {'id': 'n_56', 'connects': ['point_34', 'digestive_9-1']},
        #     {'id': 'n_57', 'connects': ['point_35', 'plexus_6-1']},
        #     {'id': 'n_58', 'connects': ['point_36', 'urinary_13-1']},
        #     {'id': 'n_39', 'connects': ['point_29', 'label_15-1']},
        # ]

        elementIdentifier = 1
        #  create elements of each centreline and repeat for other centrelines
        for name, centreline_values in centreline_data.items():
            el = [centreline_values['connects'][0]]
            try:
                # el.extend(centreline_values['contained-in'])
                for c in centreline_values['contained-in']:
                    if c not in skip_names:
                        el.append(c)
            except KeyError:
                pass
            el.extend(centreline_values['connects'][1:])
            centreline = []
            for i in range(len(el)-1):
                centreline.append({'id': name, 'group': centreline_values['group'],
                                   'connects': [el[i], el[i+1]]})
            #  create elements of the centreline
            for e3 in range(len(centreline)):
                if 'L6' in centreline[e3]['connects'][0] or 'T13' in centreline[e3]['connects'][0]:
                    continue
                if 'L6' in centreline[e3]['connects'][1] or 'T13' in centreline[e3]['connects'][1]:
                    continue
                elementtemplate = elementtemplate1
                eft = eft1
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result = element.setNodesByIdentifier(eft, [marker_nodeid_map[centreline[e3]['connects'][0]], marker_nodeid_map[centreline[e3]['connects'][1]]])
                result = element.setNodesByIdentifier(eftRadius, [marker_nodeid_map[centreline[e3]['connects'][0]], marker_nodeid_map[centreline[e3]['connects'][1]]])
                annotationGroup = AnnotationGroup(region, centreline[e3]['group'])
                annotationGroups.append(annotationGroup)
                annotationMeshGroup = annotationGroup.getMeshGroup(mesh)
                annotationMeshGroup.addElement(element)
                elementIdentifier = elementIdentifier + 1
                if 'ardell' in centreline[e3]['group'][0]:
                    annotationGroup = AnnotationGroup(region, ('ardell_all', ''))
                    annotationGroups.append(annotationGroup)
                    annotationMeshGroup = annotationGroup.getMeshGroup(mesh)
                    annotationMeshGroup.addElement(element)
                for c in ['keast', 'n_58', 'bladder_n', 'pelvic_splanchnic_n', 'hypogastric_n',
                          'lumbar_splanchnic_n', 'pudendal_n']:
                    if c in centreline[e3]['group'][0]:
                        annotationGroup = AnnotationGroup(region, ('keast_all', ''))
                        annotationGroups.append(annotationGroup)
                        annotationMeshGroup = annotationGroup.getMeshGroup(mesh)
                        annotationMeshGroup.addElement(element)


        # Make a group consists of all segments
        is_vagus = fm.createFieldConstant(1)
        markerGroup = findOrCreateFieldGroup(fm, 'centreline')
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerPoints.addNodesConditional(is_vagus)
        elementGroup = findOrCreateFieldGroup(fm, 'centreline')
        vagusMeshGroup = elementGroup.createFieldElementGroup(mesh).getMeshGroup()
        # vagusMeshGroup = vagusGroup.getMeshGroup(mesh)
        vagusMeshGroup.addElementsConditional(is_vagus)

        fm.endChange()
        return annotationGroups


def read_json(filename):
    with open(filename, 'r') as f:
        dicf = json.loads(f.read())
    return dicf


def write_json(filename, dic):
    with open(filename, 'w') as f:
        json.dump(dic, f, indent=4)


def read_manifest(manifest_file):
    # input_file = 'C:\flatmap\rat-flatmap\manifest.json'
    manifest_dict = read_json(os.path.join(manifest_file))
    file_dirname = os.path.dirname(manifest_file)
    properties_data = read_json(os.path.join(file_dirname, manifest_dict['properties']))
    anatomicalMap_data = read_json(os.path.join(file_dirname, manifest_dict['anatomicalMap']))
    return properties_data, anatomicalMap_data


def get_centrelines_list(properties_data):
    return properties_data['networks'][0]['centrelines']


def get_points_in_centrelines(centrelines):  #, properties_data, anatomicalMap_data):
    def _get_point_in_centreline(centreline):
        pl = deepcopy(centreline['connects'])
        try:
            pl.extend(centreline['contained-in'])
        except KeyError:
            pass
        pl = list(set(pl))
        pl.sort()
        return pl

    points = []
    for centreline in centrelines:
        points.extend(_get_point_in_centreline(centreline))
        points = list(set(points))
        points.sort()
    return points


def write_points_to_json(names_internal_anatomical_map, points, input_file='', new_coordinates=[]):
    nodes_data = {}
    rng = np.random.default_rng(seed=10)
    for p in points:
        d = {"x": (300 * np.random.rand(3) + [-150, -150, 700]).tolist(), "d1": [10.0, 10.0, 10.0], "r": 6.0}
        # Get the previously given coordinates if exist
        if input_file:
            previous_points = read_json(input_file)
            if p in names_internal_anatomical_map:
                try:
                    d = previous_points[p]
                except KeyError:
                    pass
        nodes_data[p] = d

    filename = os.path.join(sys.argv[2], 'coordinates_modified2.json')
    write_json(filename, nodes_data)


def get_anatomical_names(properties_data, anatomicalMap_data):
    d = {}
    for key in properties_data['features']:
        try:
            d[key] = anatomicalMap_data[properties_data['features'][key]['class']]['name']
        except KeyError:
            pass
    return d


def get_nerve_terms(points, centrelines, properties_data, anatomicalMap_data):
    def get_name_class(p):
        try:
            return properties_data['features'][p]['class']
        except KeyError:
            return p

    def nerve_term(internal_names):
        name_class = []
        terms = {}
        for name in internal_names:
            name_cls = get_name_class(name)
            # print only unique values
            if name_cls in anatomicalMap_data:
                # if name_cls not in name_class:
                terms[name] = (anatomicalMap_data[name_cls]['name'], anatomicalMap_data[name_cls]['term'])
            else:
                terms[name] = (name, '')
            name_class.append(name_cls)
        return terms

    node_terms = nerve_term(points)
    centreline_names = []
    for c in centrelines:
        centreline_names.append(c['id'])
    centreline_terms = nerve_term(centreline_names)
    return node_terms, centreline_terms


def get_marker_data(nodes_data, node_terms):
    marker_data = {}
    for key in node_terms:
        values = nodes_data[key]
        anatomical_name = node_terms[key][0]
        marker_data[key] = (
            {"group": nerve_terms.get_nerve_term(anatomical_name),
             "x": values['x'], "d1": values['d1'], "r": values['r']}
        )
    return marker_data


def get_centrline_data(centrelines, centreline_terms):
    centrline_data = {}
    for centreline in centrelines:
        name = centreline['id']
        anatomical_name = centreline_terms[name][0]
        try:
            values = centreline['contained-in']
        except KeyError:
            values = []
        centrline_data[name] = (
            {"group": nerve_terms.get_nerve_term(anatomical_name), "connects": centreline['connects'],
             "contained-in": values}
        )
    return centrline_data


def modify_points_in_scaffoldmaker(names_internal_anatomical_map, marker_data, scaffoldmaker_nodes):
    new_coordinates = [{'x': c[1][0], 'd1': c[1][1], 'r': 6.0} for c in scaffoldmaker_nodes]
    nodes_data = {}
    c = 0
    for name, marker in marker_data.items():
        # if name == 'spinal_47-1':
        #     continue
        # marker_data order and new_coordinates should be the same
        d = new_coordinates[c]
        # Get the previously given coordinates if exist
        if name not in names_internal_anatomical_map:
            nodes_data[name] = d
        else:
            nodes_data[name] = {'x': marker_data[name]['x'], 'd1': marker_data[name]['d1'], 'r': marker_data[name]['r']}
        c += 1

    filename = os.path.join(sys.argv[2], 'coordinates_modified9.json')
    write_json(filename, nodes_data)
