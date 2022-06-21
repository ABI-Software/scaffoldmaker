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

        names_internal_anatomical_map = {
            "brain_44-1": "trigeminal nucleus",
            "brain_12-1": "nucleus of solitary tract",
            "brain_34-1": "dorsal motor nucleus of vagus nerve",
            "brain_33-1": "nucleus ambiguus",
            "point_1": "point_1",
            "ardell_4_branching_point": "ardell_4_branching_point",
            "ganglion_1-1": "superior vagus X ganglion",
            "bolser_X-1": "bolser_X-1",
            "ardell_3_branching_point": "ardell_3_branching_point",
            "point_3": "point_3",
            "label_1-1": "dura mater in the posterior cranial fossa",
            "ganglion_2-1": "inferior vagus X ganglion",
            "point_6": "point_6",
            "point_7": "point_7",
            "point_8": "point_8",
            "point_9": "point_9",
            "point_10": "point_10",
            "label_4-1": "tympanic membrane",
            "label_3-1": "floor of external acoustic meatus",
            "label_2-1": "posterior wall of external acoustic meatus",
            "point_11": "point_11",
            "point_13": "point_13",
            "label_6-1": "pharynx)",
            "label_5-1": "back part of the tongue*",
            "point_14": "point_14",
            "cardio_4-1": "carotid body",
            "point_16": "point_16",
            "point_18": "point_18",
            "point_19": "point_19",
            "point_20": "point_20",
            "label_7-1": "mucosa of pharynx",
            "label_8-1": "epiglottic vallecula",
            "respiratory_11-1": "epiglottis",
            "respiratory_5-1": "Wall of larynx",
            "point_21": "point_21",
            "label_10-1": "aryepiglottic fold",
            "label_11-1": "arytenoideus",
            "label_12-1": "mucosa of arytenoid cartilage",
            "point_37": "point_37",
            "bolser_X-5": "bolser_X-5",
            "bolser_X-4": "bolser_X-4",
            "external_laryngeal_n_branching_point": "external_laryngeal_n_branching_point",
            "plexus_2-1": "pharyngeal nerve plexus",
            "bolser_X-14": "bolser_X-14",
            "label_15-1": "inferior pharyngeal constrictor",
            "label_14-1": "cricothyroid muscle",
            "bolser_X-15": "bolser_X-15",
            "point_41": "point_41",
            "point_22": "point_22",
            "point_24": "point_24",
            "point_25": "point_25",
            "plexus_3-1": "cardiac nerve plexus",
            "point_26": "point_26",
            "point_28": "point_28",
            "label_9-1": "mucosa of larynx",
            "point_29": "point_29",
            "label_16-1": "ceratocricoid",
            "label_17-1": "lateral crico-arytenoid",
            "label_18-1": "oblique arytenoid",
            "label_19-1": "posterior crico-arytenoid",
            "label_20-1": "hyro-arytenoid",
            "label_21-1": "transverse arytenoid",
            "label_22-1": "vocalis muscle",
            "label_13-1": "vocal_cords",
            "respiratory_8-1": "Laryngeal mechanoreceptors",
            "point_38": "point_38",
            "bolser_X-2": "bolser_X-2",
            "bolser_X-3": "bolser_X-3",
            "point_30": "point_30",
            "plexus_4-1": "pulmonary nerve plexus",
            "point_31": "point_31",
            "respiratory_17-1": "bronchus smooth muscle",
            "point_32": "point_32",
            "point_33": "point_33",
            "plexus_5-1": "esophageal nerve plexus",
            "digestive_12-1": "esophagus",
            "point_34": "point_34",
            "digestive_9-1": "pancreas",
            "point_35": "point_35",
            "plexus_6-1": "left gastric nerve plexus",
            "point_36": "point_36",
            "urinary_13-1": "kidney",
        }

        # markers group, coordinates and derivatives
        update_terms = False
        nerve_terms_changed = False
        filename = os.path.join(working_dir, 'coordinates_modified.json')
        coordinates_file = os.path.join(working_dir, 'coordinates_modified3.json')
        properties_data, anatomicalMap_data = read_manifest(manifest_file)
        centrelines = get_centrelines_list(properties_data)
        points = get_points_in_centrelines(centrelines)

        if update_terms:
            write_points_to_json(names_internal_anatomical_map, points, filename)
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

        def modify_points_in_scaffoldmaker(names_internal_anatomical_map, marker_data, scaffoldmaker_nodes):
            new_coordinates = [{'x': c[1][0], 'd1': c[1][1], 'r': 6.0} for c in scaffoldmaker_nodes]
            nodes_data = {}
            c = 0
            for name, marker in marker_data.items():
                # marker_data order and new_coordinates should be the same
                d = new_coordinates[c]
                # Get the previously given coordinates if exist
                if name not in names_internal_anatomical_map:
                    nodes_data[name] = d
                else:
                    nodes_data[name] = {'x': marker_data[name]['x'], 'd1': marker_data[name]['d1'], 'r': marker_data[name]['r']}
                c += 1

            filename = os.path.join(sys.argv[2], 'coordinates_modified3.json')
            write_json(filename, nodes_data)

        for name, marker in marker_data.items():
            # if n3+1 in nodes_w_ver:
            #     nt = nodetemplate1
            # else:
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
                el.extend(centreline_values['contained-in'])
            except KeyError:
                pass
            el.extend(centreline_values['connects'][1:])
            centreline = []
            for i in range(len(el)-1):
                centreline.append({'id': name, 'group': centreline_values['group'],
                                   'connects': [el[i], el[i+1]]})
            #  create elements of the centreline
            for e3 in range(len(centreline)):
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

# [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ]
# [
# ( 1, [ [   3.455, -48.725,1373.805], [   7.159, -42.376, -49.194] ] ),
# ( 2, [ [   1.901, -39.784,1413.770], [   6.816,  14.543,  33.326] ] ),
# ( 3, [ [   0.953, -40.783,1402.753], [  21.290, -19.940, -25.742] ] ),
# ( 4, [ [   2.854, -35.733,1403.720], [   8.980,  23.114,  15.060] ] ),
# ( 5, [ [   2.836, -27.054,1390.265], [   4.748,  24.711,  -4.111] ] ),
# ( 6, [ [   0.634, -31.872,1383.777], [   6.298,  -8.085, -10.669] ] ),
# ( 7, [ [   1.636, -25.119,1384.707], [   2.911,  23.509,  -9.258] ] ),
# ( 8, [ [   0.593, -25.557,1369.152], [   0.031,  -0.122,   0.049] ] ),
# ( 9, [ [   0.028, -27.692,1362.915], [  -6.473,  -6.356, -55.570] ] ),
# ( 10, [ [  -4.944, -35.575,1287.670], [ -10.070,  11.760, -79.078] ] ),
# ( 11, [ [  -2.585, -20.015,1345.160], [  48.012,  49.983,  42.410] ] ),
# ( 12, [ [  -2.938, -22.195,1339.684], [  -1.725,  -1.113, -32.115] ] ),
# ( 13, [ [  -5.521, -14.114,1326.391], [ -25.845,  53.045,  39.940] ] ),
# ( 14, [ [  -4.276, -12.900,1317.935], [   3.147, -12.875, -33.054] ] ),
# ( 15, [ [  -8.418,  -6.384,1304.632], [  11.265,  98.719,  31.090] ] ),
# ( 16, [ [  -1.810, -18.396,1285.229], [   4.985, -27.015, -11.020] ] ),
# ( 17, [ [  -7.549, -72.483, 949.168], [ -37.820, -28.122,   9.726] ] ),
# ( 18, [ [  18.345, -64.209, 957.825], [ -16.261,  13.278, 106.361] ] ),
# ( 19, [ [   9.730, -50.250, 842.635], [  31.336, -12.348, -36.113] ] ),
# ( 20, [ [  -1.258, -54.298, 923.134], [  -6.207,  19.293,  -5.432] ] ),
# ( 21, [ [  -4.672, -48.725, 939.781], [  42.309,  18.100, -70.678] ] ),
# ( 22, [ [  11.035, -47.573, 833.262], [  12.555,  16.731, -13.915] ] ),
# ( 23, [ [  -4.008, -39.402, 915.140], [  21.480, -33.279, -49.652] ] ),
# ( 24, [ [  11.931, -43.800, 926.036], [  39.362,  -0.730,  24.867] ] ),
# ( 25, [ [   2.271, -43.232,1026.651], [ -11.764, -18.545, 257.721] ] ),
# ( 26, [ [  -2.907, -52.069, 887.133], [  46.157,   7.277,-106.324] ] ),
# ( 27, [ [  15.200, -37.581, 816.921], [   3.910,  -2.296, -43.440] ] ),
# ( 28, [ [   1.588, -37.797,1030.716], [ -32.940,   9.759,  96.576] ] ),
# ( 29, [ [  -0.105, -56.179, 885.988], [  42.328,   3.471, -76.006] ] ),
# ( 30, [ [  19.267, -34.944, 794.137], [ -13.437,   0.782, -16.575] ] ),
# ( 31, [ [   9.959, -52.538,1373.268], [ -17.780, -23.273,  -9.355] ] ),
# ( 32, [ [  -3.177, -31.965,1380.796], [ -11.517,   7.077, -19.632] ] ),
# ( 33, [ [  -5.173, -25.199,1259.076], [  -3.205,  15.737, -72.210] ] ),
# ( 34, [ [  13.681, -39.517,1237.581], [  21.699,   7.687, 147.931] ] ),
# ( 35, [ [   1.492, -50.482, 817.638], [  42.851, -30.781, 207.969] ] ),
# ( 36, [ [  -0.682, -45.918, 918.367], [  19.809,  -7.831, -10.676] ] ),
# ( 37, [ [   6.952, -50.984, 864.099], [   9.296,   6.825, -24.566] ] ),
# ( 38, [ [  12.414, -54.131, 882.909], [  14.772,  42.409, -30.746] ] ),
# ( 39, [ [   2.756, -37.888,1226.694], [  17.772, -58.644, -29.571] ] ),
# ( 40, [ [ -15.208, -65.938,1039.556], [  15.996, -17.254, -38.828] ] ),
# ( 41, [ [  17.975, -98.393, 997.599], [  -9.408,  -5.118, -54.980] ] ),
# ( 42, [ [   6.392, -40.549,1021.816], [ -18.767,   2.768, -19.495] ] ),
# ( 43, [ [   9.719, -50.506,1048.404], [  11.839, -11.246,  -9.134] ] ),
# ( 44, [ [  -0.683, -28.936,1220.403], [  11.275,  -9.109, -34.126] ] ),
# ( 45, [ [   9.179, -42.523,1099.168], [   4.018, -32.972, -85.401] ] ),
# ( 46, [ [  -8.526, -50.547,1125.811], [  26.270, -28.957, -54.131] ] ),
# ( 47, [ [  -5.298, -92.084,1148.966], [  -5.472, -38.852, -79.123] ] ),
# ( 48, [ [   9.658, -53.026,1117.141], [  15.487, -32.808,   1.700] ] ),
# ( 49, [ [  -2.626, -22.934,1110.372], [  -6.729, -23.803, -17.327] ] ),
# ( 50, [ [   2.634, -70.539,1386.866], [ 104.107, -80.434, -81.714] ] ),
# ( 51, [ [   3.387, -36.733,1073.751], [  11.738, -13.322, -11.814] ] ),
# ( 52, [ [  24.604, -49.926,1438.819], [  38.293,  -7.439,  -3.074] ] ),
# ( 53, [ [  13.305, -31.413,1380.814], [  17.454,   1.948, -14.650] ] ),
# ( 54, [ [  -5.772, -18.583,1340.752], [  -0.710,   8.886, -23.109] ] ),
# ( 55, [ [ -12.311, -28.832,1297.535], [   7.033,  -0.399, -40.080] ] ),
# ( 56, [ [  -3.276, -40.161,1266.565], [  -6.502,  -1.879, -49.531] ] ),
# ( 57, [ [  -4.744, -27.088,1250.960], [ -13.293, -11.248, -24.204] ] ),
# ( 58, [ [  20.414, -59.061, 939.728], [ -27.690, -13.635, -23.475] ] ),
# ( 59, [ [  -1.690, -15.350, 926.330], [   9.728,  -1.829, -71.950] ] ),
# ( 60, [ [   8.265, -37.053, 860.079], [   4.462, -24.081, -41.968] ] ),
# ( 61, [ [  12.305, -45.622, 849.311], [   2.122, -17.456, -41.190] ] ),
# ( 62, [ [   1.645, -53.336,1293.589], [   9.247, -36.255, -15.738] ] ),
# ( 63, [ [   3.420, -42.089,1182.132], [  11.698, -46.306, -37.647] ] ),
# ( 64, [ [ -17.052, -50.759,1178.264], [   9.170, -10.482, -52.206] ] ),
# ( 65, [ [  -8.003, -29.381,1173.898], [  -0.324, -55.700,-106.426] ] ),
# ( 66, [ [   7.117, -36.146,1138.923], [  -3.704, -24.080, -28.794] ] ),
# ( 67, [ [   1.171, -28.825,1108.790], [  -5.129,  15.233,   3.383] ] ),
# ( 68, [ [  21.013, -46.483,1416.194], [ -14.730,  -4.852, -32.613] ] ),
# ( 69, [ [  23.793, -28.112,1354.889], [ -10.743, -27.002, -24.854] ] ),
# ( 70, [ [  -1.803, -12.171,1322.728], [  22.706, -21.326, -13.338] ] ),
# ( 71, [ [   7.168, -23.875,1267.339], [  13.765, -11.259, -43.997] ] ),
# ( 72, [ [ -12.018, -23.376,1241.316], [  34.460, -12.809,  24.367] ] ),
# ( 73, [ [ -20.227, -36.412,1240.365], [ -15.764, -26.870,  -2.197] ] ),
# ( 74, [ [  -2.739, -73.060, 932.590], [ -13.155,  -5.058, -32.893] ] ),
# ( 75, [ [   3.411, -85.988, 870.941], [  10.207, -68.624,  25.931] ] ),
# ( 76, [ [  -2.268, -89.016, 846.419], [  -6.422, -52.462,   5.734] ] ),
# ( 77, [ [   1.629, -83.470, 828.719], [  24.282, -20.716, -29.849] ] ),
# ( 78, [ [   6.746, -64.581,1262.879], [ -38.967,  18.762,  -6.710] ] ),
# ( 79, [ [   1.343, -42.576,1208.391], [   5.656, -26.947,  45.332] ] ),
# ( 80, [ [   4.202, -31.197,1141.803], [  23.713,  31.414,   9.052] ] ),
# ( 81, [ [  -5.901, -87.907, 978.117], [  44.986, -31.290, -56.321] ] ),
# ( 82, [ [  -2.936, -47.936,1128.183], [ -17.184,   0.327,  22.667] ] ),
# ( 83, [ [  -5.867, -18.050,1117.061], [   0.268,  -5.041,  29.619] ] ),
# ( 84, [ [   6.804, -59.999, 982.311], [ -52.580,  87.830, 106.155] ] ),
# ( 85, [ [   4.306, -49.330,1057.172], [ -19.391,   4.141,-271.906] ] ),
# ( 86, [ [   0.008, -39.739, 944.118], [ -15.756,  29.125,-119.097] ] ),
# ( 87, [ [  -0.891, -41.022, 923.098], [ -19.507,  16.443,-117.142] ] ),
# ( 88, [ [  17.231, -44.521, 877.407], [  14.191,  -4.323, -76.508] ] ),
# ( 89, [ [   6.623, -41.300, 831.502], [  -7.302,  -0.403, -43.686] ] ),
# ( 90, [ [  -1.396, -58.172, 792.537], [ -14.362,   9.561,  -2.386] ] ),
# ( 91, [ [  -0.812, -46.231, 800.934], [ -19.016,   4.630, -10.216] ] ),
# ( 92, [ [   5.205, -47.063,1002.431], [   8.646,  -0.992, -12.808] ] ),
# ( 93, [ [  -1.798, -33.884,1015.901], [ -21.134,   4.418, 208.651] ] ),
# ( 94, [ [   3.796, -42.470,1034.805], [  24.371,-120.764,-124.976] ] ),
# ( 95, [ [  -2.535, -22.132,1036.624], [ -13.929,  24.452, -64.542] ] ),
# ( 96, [ [  -3.135, -22.908,1251.136], [ -40.491,  82.112, -60.471] ] ),
# ( 97, [ [   0.363, -41.275,1184.820], [  -1.554,  49.665, 129.944] ] ),
# ( 98, [ [  -3.655, -46.065,1152.646], [ -25.193, -44.686, 125.455] ] ),
# ( 99, [ [ -16.460, -57.749,1163.117], [ -26.727,  94.308, 137.090] ] ),
# ( 100, [ [   1.697, -41.425,1110.449], [  19.549, -42.359, 126.257] ] ),
# ( 101, [ [   0.563, -36.107,1126.778], [   1.342, -27.772,-131.378] ] ),
# ( 102, [ [  -5.330, -60.702,1052.430], [ -43.526, -42.199,-141.567] ] ),
# ( 103, [ [  13.782, -47.932,1079.049], [   6.930,  19.852,   7.798] ] ),
# ( 104, [ [  21.746, -50.031, 917.765], [   5.642,  -8.050, -52.323] ] ),
# ( 105, [ [  22.666, -58.716, 825.397], [   9.521,   5.526,  20.866] ] ),
# ( 106, [ [  -1.787, -53.069,1007.754], [  -1.206,  15.692, -10.545] ] ),
# ( 107, [ [   8.775, -99.877, 972.872], [ -16.678,  42.382, -55.692] ] ),
# ( 108, [ [   2.402, -49.203,1001.011], [  46.666,  -2.380, -35.157] ] ),
# ( 109, [ [   3.189, -26.409,1043.291], [  -6.835,  29.694, -19.606] ] ),
# ( 110, [ [   6.247, -37.598,1220.023], [  -3.622,  25.658, -98.274] ] ),
# ( 111, [ [  -2.758, -64.386,1106.560], [ -32.583,   9.176,-161.296] ] ),
# ( 112, [ [   9.215, -43.834,1110.750], [ -25.815,  19.444,-144.431] ] ),
# ( 113, [ [   3.139, -56.256,1123.526], [ -10.325, -44.251,-142.793] ] ),
# ( 114, [ [  -5.615, -30.338,1103.264], [   2.612,   5.132,-146.767] ] ),
# ( 115, [ [   8.254, -54.839,1095.467], [  51.600,  -2.988, -12.963] ] ),
# ( 116, [ [  11.020, -62.812,1055.156], [ -57.537,  12.981,   7.205] ] ),
# ( 117, [ [   5.591, -44.864,1070.139], [   2.478, -10.959,  -2.839] ] ),
# ( 118, [ [   0.209, -20.161,1362.145], [   3.333,  26.140,-147.627] ] ),
# ( 119, [ [  -1.249, -62.200, 899.203], [   8.290,   7.561,-113.841] ] ),
# ( 120, [ [   3.121, -57.169, 997.867], [  -7.465, -18.880, -55.451] ] ),
# ( 121, [ [  -6.453, -78.752,1021.977], [  20.700,   2.385, -58.504] ] ),
# ( 122, [ [  12.382,-105.827,1024.730], [  -1.364, -18.723, -75.948] ] ),
# ( 123, [ [  17.462, -46.752,1248.722], [ -31.629, -52.245, -86.448] ] ),
# ( 124, [ [  -2.907,  -6.555,1234.286], [   1.143,  15.275, -91.456] ] ),
# ( 125, [ [  -1.376, -34.543,1171.408], [  -4.451,  18.717,-112.562] ] ),
# ( 126, [ [  -8.069, -42.579,1140.594], [  12.620,  24.618,-105.461] ] ),
# ( 127, [ [ -12.946, -37.112,1156.058], [   5.613,  21.397,-113.579] ] ),
# ( 128, [ [   0.863, -36.252,1119.870], [ -10.129,  23.807,-108.979] ] ),
# ( 129, [ [   7.355, -67.146,1112.033], [ -28.421,  30.498, -20.321] ] ),
# ( 130, [ [   1.518, -82.824,1064.515], [   7.128, -14.528, -88.991] ] ),
# ( 131, [ [   5.389, -25.178,1070.145], [ -21.412,   0.996,  -8.502] ] ),
# ( 132, [ [  -1.356, -60.566,1387.978], [  -0.740,  26.090,  16.122] ] ),
# ( 133, [ [  -5.621, -55.836,1389.921], [  -9.524, -11.719,  14.685] ] ),
# ( 134, [ [   1.107, -52.356,1373.738], [  -4.645,   8.435,  40.068] ] ),
# ( 135, [ [ -11.842, -15.988,1350.643], [   2.077, -28.480, -36.110] ] ),
# ( 136, [ [   8.069, -40.880,1347.005], [  -2.859, -16.445,  63.249] ] ),
# ( 137, [ [  -5.746, -12.727,1240.584], [  15.809,   4.687,  75.821] ] ),
# ( 138, [ [   1.045, -51.145,1237.363], [ -25.700, -39.149,  74.120] ] ),
# ( 139, [ [   5.975, -45.508,1275.795], [   7.147,   0.932,  54.618] ] ),
# ( 140, [ [ -19.891, -66.033,1255.339], [   9.182,  -1.316,  60.601] ] ),
# ( 141, [ [   2.030, -28.867,1231.570], [ -21.140,  25.360,  83.264] ] ),
# ( 142, [ [  20.933, -56.876, 921.263], [  -6.517,   6.299,  42.838] ] ),
# ( 143, [ [   5.336, -57.558,1081.998], [  13.034, -11.589, 289.140] ] ),
# ( 144, [ [   5.350, -80.254, 911.167], [  -4.665,  25.334,  31.417] ] ),
# ( 145, [ [  -6.560, -45.433, 892.133], [ -18.785,   4.300,  32.577] ] ),
# ( 146, [ [  18.244, -46.456, 901.466], [  -1.550,   0.364, -30.327] ] ),
# ( 147, [ [   5.242, -56.679, 859.349], [   6.299,  46.191,  66.425] ] ),
# ( 148, [ [  14.126, -40.208, 842.438], [   5.598,   3.038, -33.922] ] ),
# ( 149, [ [  23.410, -62.124, 830.614], [  81.483,  49.137,  37.778] ] ),
# ( 150, [ [   2.726,   7.249, 842.605], [  12.751,   6.554, -85.795] ] ),
# ( 151, [ [   4.425, -52.624,1280.025], [ -16.755,  66.142, -42.270] ] ),
# ( 152, [ [  14.871, -67.821,1014.986], [   2.481,  17.304, -28.566] ] ),
# ( 153, [ [   7.917, -84.692, 996.048], [ -14.051,   4.059, -27.244] ] ),
# ( 154, [ [ -14.942, -44.855,1031.816], [ -18.092, -15.004, -24.324] ] ),
# ( 155, [ [   4.995, -40.207,1051.291], [  -4.128,  15.032,  -3.046] ] ),
# ( 156, [ [ -12.638, -45.216,1295.599], [ -29.454, -13.997,  14.843] ] ),
# ( 157, [ [   1.176, -25.291,1228.877], [   3.680, -31.332,  79.849] ] ),
# ( 158, [ [  15.143, -90.487,1245.293], [  -0.824,   6.012,  84.382] ] ),
# ( 159, [ [   3.418, -47.473,1106.922], [ -24.269, -24.744, 137.544] ] ),
# ( 160, [ [   0.902, -30.725,1165.665], [  -6.890,   6.853,  64.586] ] ),
# ( 161, [ [   9.772, -46.030,1131.599], [  38.160,  34.017, 107.337] ] ),
# ( 162, [ [   3.691, -27.873,1163.989], [ -41.358, -82.294, 117.491] ] ),
# ( 163, [ [  -9.193, -37.250,1129.545], [  29.508,  78.686,  90.799] ] ),
# ( 164, [ [   1.237, -32.048,1154.667], [   2.095,  21.950,  47.667] ] ),
# ( 165, [ [  -1.964, -32.934,1104.007], [ -35.922,  39.139,  97.114] ] ),
# ( 166, [ [  10.252, -48.078,1132.703], [   9.629, -13.165,  63.016] ] ),
# ( 167, [ [ -14.662, -44.337,1115.964], [   5.672, -29.348, -10.214] ] ),
# ( 168, [ [   4.011, -63.654,1048.707], [  24.389,   7.534,-210.111] ] ),
# ( 169, [ [   8.649, -44.189,1074.469], [   0.619,  -7.131,  -4.803] ] ),
# ( 170, [ [  -7.409, -37.858,1078.471], [  -2.267, -25.506, -36.311] ] ),
# ( 171, [ [  -2.145, -44.749, 900.640], [  -2.942,  41.934,-133.326] ] ),
# ( 172, [ [  -3.660, -46.547,1061.045], [   4.283, -36.655, -49.573] ] ),
# ( 173, [ [   2.768, -61.876, 998.441], [  -0.217,   3.008,   1.096] ] ),
# ( 174, [ [  -1.285, -46.343,1043.732], [   6.487,  -8.057, -26.039] ] ),
# ( 175, [ [   2.360, -54.764, 974.828], [  17.361,  25.807,  -0.781] ] ),
# ( 176, [ [  -0.100, -32.980,1054.911], [   5.498, -27.051, -22.632] ] ),
# ( 177, [ [   2.636, -41.291,1031.230], [ -20.503,  37.817, -79.243] ] ),
# ( 178, [ [ -12.274, -60.355,1283.117], [   8.142,  -4.392, -10.738] ] ),
# ( 179, [ [  -3.701, -23.741,1280.966], [   9.571, -16.237, -52.885] ] ),
# ( 180, [ [   1.769, -41.560,1244.864], [  -0.610,  -2.799, -20.401] ] ),
# ( 181, [ [  -6.872, -10.276,1290.300], [ -15.297,  47.886,   3.149] ] ),
# ( 182, [ [  -7.600, -16.894,1259.020], [  10.361, -25.722, -47.577] ] ),
# ( 183, [ [  -5.128,  -9.966,1136.584], [  -5.459, -21.599,-101.603] ] ),
# ( 184, [ [  -7.092, -12.976,1265.248], [   3.740,  26.994,  52.606] ] ),
# ( 185, [ [  -7.515, -13.773,1216.275], [   1.846, -57.221, -26.110] ] ),
# ( 186, [ [  -8.069, -15.546,1116.269], [  -1.921,   4.144,   0.491] ] ),
# ( 187, [ [ -11.410, -68.722,1150.673], [  -1.903,  -0.102,  24.579] ] ),
# ( 188, [ [  -7.481,  -9.250,1198.444], [   5.068, -41.835, -39.101] ] ),
# ( 189, [ [  -4.009, -23.512,1098.220], [   2.532,  -7.603, -11.652] ] ),
# ( 190, [ [  -2.469,  -8.763,1223.131], [  -9.452,  24.480,  88.753] ] ),
# ( 191, [ [  -6.003,  -9.311,1178.039], [   8.793, -49.705, -11.665] ] ),
# ( 192, [ [  -6.690, -33.539,1075.771], [   6.499, -13.425, -44.271] ] ),
# ( 193, [ [  -8.212,  -5.806,1203.326], [ -39.087,  82.031,  65.678] ] ),
# ( 194, [ [   1.413, -37.389,1111.442], [   7.519,  -4.515,  25.336] ] ),
# ( 195, [ [  -4.466, -28.164,1078.385], [   0.522,   7.509, -25.591] ] ),
# ( 196, [ [  -7.741, -15.816,1140.766], [   2.498, -18.863, -20.517] ] ),
# ( 197, [ [  -7.358, -34.534,1115.398], [  -0.697,  24.084,  18.761] ] ),
# ( 198, [ [  -8.957, -23.786,1120.051], [  33.884,-108.377, 151.894] ] ),
# ( 199, [ [  -3.504, -58.032, 953.246], [ -16.276,  60.762,-126.569] ] ),
# ( 200, [ [  -6.901, -29.725,1103.016], [  14.255, -15.179, -16.522] ] ),
# ( 201, [ [  -0.243, -38.584,1071.314], [   3.841, -10.524,   4.102] ] ),
# ( 202, [ [   4.749, -45.630,1191.393], [  10.000,  10.000,  10.000] ] ),
# ( 203, [ [  -0.741, -20.507,1319.756], [  -0.025,   3.569, -16.965] ] ),
# ( 204, [ [  -2.344, -15.896,1306.302], [   0.316,   0.569, -13.411] ] ),
# ( 205, [ [   3.067, -30.152,1298.762], [   1.811,   1.330,  14.466] ] ),
# ( 206, [ [   7.681, -32.740,1260.429], [  -3.378,   0.451,  -3.024] ] ),
# ( 207, [ [   4.341, -27.947,1245.391], [  -5.192,  -1.389,  -7.802] ] ),
# ( 208, [ [   3.545, -30.025,1217.509], [  -6.873,  -8.380,  -9.689] ] ),
# ( 209, [ [  -0.063, -18.974,1184.760], [  -6.071, -13.913,  -7.293] ] ),
# ( 210, [ [  -1.152, -26.119,1165.248], [  -0.268,  -4.741,   4.101] ] ),
# ( 211, [ [   1.414, -32.945,1133.265], [  -6.148,  10.961,  14.976] ] ),
# ( 212, [ [   4.993, -54.083,1256.725], [  12.261,   7.312,   5.792] ] ),
# ( 213, [ [ -19.628, -40.944,1239.340], [  27.327, -39.348,  24.416] ] ),
# ( 214, [ [  -3.084,-102.096,1260.982], [ -24.887,  12.822, -66.630] ] ),
# ( 215, [ [   3.296, -68.465,1338.180], [   4.507,   9.351,-244.960] ] ),
# ( 216, [ [  15.337, -81.890,1209.048], [  10.000,  10.000,  10.000] ] ),
# ( 217, [ [   2.032, -35.215,1168.344], [ -11.391,  -4.786, -49.995] ] ),
# ( 218, [ [ -10.817, -31.445,1148.091], [ -21.184,   2.545, -64.332] ] ),
# ( 219, [ [  -1.017, -39.530,1121.162], [  -0.210,  -8.064, -44.924] ] ),
# ( 220, [ [   5.965, -54.985,1094.610], [   7.044,   4.256, -38.288] ] ),
# ( 221, [ [  -7.113, -65.443, 987.279], [  13.914, -36.684,-117.832] ] ),
# ( 222, [ [   1.980, -43.737,1084.341], [   2.431,  -4.833,  22.860] ] ),
# ( 223, [ [  -0.127, -82.599,1202.228], [  10.000,  10.000,  10.000] ] ),
# ( 224, [ [   2.345, -53.825,1395.294], [  10.672,   6.902, -39.806] ] ),
# ( 225, [ [   8.159, -50.724,1364.578], [  43.808,  50.207,  -1.666] ] ),
# ( 226, [ [   0.920, -26.701,1307.745], [  19.410,  18.849,  56.667] ] ),
# ( 227, [ [   5.841, -37.715,1204.181], [   9.468, -11.657,-100.217] ] ),
# ( 228, [ [  12.733, -60.748,1203.597], [  29.955, -61.265, -62.535] ] ),
# ( 229, [ [ -20.818, -35.407,1205.060], [   3.370,  -1.653,  -0.436] ] ),
# ( 230, [ [  -4.411, -38.644,1134.152], [   1.015, -16.055, -84.740] ] ),
# ( 231, [ [   0.625, -53.557,1168.361], [   4.458, -38.753,  -7.683] ] ),
# ( 232, [ [  10.939, -53.961,1138.047], [  10.987, -27.275,  26.556] ] ),
# ( 233, [ [   6.977,-143.098, 803.603], [  -0.318,   1.086,  -1.479] ] ),
# ( 234, [ [  -9.021, -78.999,1291.725], [  -5.163,  61.814, -71.337] ] ),
# ( 235, [ [  17.504, -86.497,1252.626], [ -40.337,  -0.939,-309.083] ] ),
# ( 236, [ [ -15.705, -53.409,1281.176], [ -17.608,  -5.627,-207.834] ] ),
# ( 237, [ [  12.604,-115.494,1510.482], [  26.349, -65.975, -58.140] ] ),
# ( 238, [ [  -8.915,  -6.574,1211.738], [  -3.550, -62.802, 124.226] ] ),
# ( 239, [ [  27.594,-146.970,1475.180], [   2.251,   1.874,  -7.725] ] ),
# ( 240, [ [   0.526, -56.094,1488.768], [  -6.052,  -3.281,-147.777] ] ),
# ( 241, [ [   2.607, -56.401,1489.594], [   6.595,  -5.741, -13.996] ] ),
# ( 242, [ [   1.878, -56.261,1489.300], [   5.962,  -5.903, -15.305] ] ),
# ( 243, [ [  -3.068, -58.888,1485.941], [  -1.068,  -8.877, -18.460] ] ),
# ( 244, [ [ -15.128, -90.734,1239.663], [   6.966, -62.105, -51.336] ] ),
# ( 245, [ [  -6.516, -11.417,1238.790], [  14.025, -71.971, -73.016] ] ),
# ( 246, [ [   4.051, -80.134, 952.868], [  -0.155,  -1.132,  -1.261] ] ),
# ( 247, [ [  19.268, -95.333,1427.899], [ -16.905, -26.707, -47.062] ] ),
# ( 248, [ [  19.140, -71.172,1471.654], [  38.042, -29.689,-121.480] ] ),
# ( 249, [ [  -3.782, -16.071,1236.370], [  21.082, -84.809, -71.818] ] ),
# ( 250, [ [  35.504, -78.359,1463.385], [ -31.406,  18.032,   4.080] ] ),
# ( 251, [ [ -11.065, -63.248,1435.031], [  39.574, -10.233, -41.578] ] ),
# ( 252, [ [  -4.885, -72.882, 860.710], [ -39.062, -35.820,   8.014] ] ),
# ( 253, [ [  13.524, -85.505, 853.309], [  29.587, -31.858, -19.279] ] ),
# ( 254, [ [   4.853,-128.079, 823.090], [   5.566, -34.542, -32.836] ] ),
# ( 255, [ [   1.702, -63.954, 808.886], [  -0.567,  30.873,  -8.715] ] ),
# ( 256, [ [  19.635, -50.920, 811.273], [  -9.831,   4.611,  49.735] ] ),
# ( 257, [ [  16.025, -77.837, 861.976], [ -19.674,  14.367,   8.814] ] ),
# ( 258, [ [  -5.261, -35.953, 894.691], [  -4.846, -10.428, -64.773] ] ),
# ( 259, [ [  19.059, -76.211,1480.118], [  -2.900,   9.504,  27.931] ] ),
# ( 260, [ [   2.098,-109.068,1405.412], [ -13.903, -11.924, -30.233] ] ),
# ( 261, [ [   0.755,-107.608,1406.119], [ -16.558,  -8.990, -28.768] ] ),
# ( 262, [ [   0.392,-105.992,1406.523], [ -17.229,  -5.748, -27.876] ] ),
# ( 263, [ [  12.808,-107.604,1396.148], [ -19.636,  -2.129,  -2.232] ] ),
# ( 264, [ [   6.623,-108.164,1402.941], [  -5.133,   0.639,  -1.747] ] ),
# ( 265, [ [  15.735,-106.477,1387.002], [ -14.654,  -8.380, -26.468] ] ),
# ( 266, [ [  16.942,-106.518,1385.462], [ -15.075,  -0.833, -15.714] ] ),
# ( 267, [ [  18.305,-110.622,1389.761], [ -12.813,  -9.045,  -7.615] ] ),
# ( 268, [ [  17.812,-111.228,1390.543], [ -13.704, -10.185,  -6.014] ] ),
# ( 269, [ [  20.305,-109.662,1385.525], [  -8.642,  -6.984, -15.735] ] ),
# ( 270, [ [ -65.468, -70.375,1520.633], [ -80.971, -22.695,  37.571] ] ),
# ( 271, [ [  19.183,-111.350,1388.633], [ -11.055, -10.492,  -9.859] ] ),
# ( 272, [ [  18.880,-110.329,1387.243], [ -11.656,  -8.454, -12.627] ] ),
# ( 273, [ [  19.598,-110.678,1389.008], [ -10.237,  -9.160,  -9.120] ] ),
# ( 274, [ [ -57.275, -70.875,1519.911], [ -65.695, -23.698,  36.302] ] ),
# ( 275, [ [  60.913, -52.711,1519.712], [ 113.397,   7.643,  25.844] ] ),
# ( 276, [ [   9.582,-154.877,1443.576], [  10.727,  -4.838,  -5.242] ] ),
# ( 277, [ [ -14.000,-152.355,1446.593], [ -21.628,  -2.840,  -2.635] ] ),
# ( 278, [ [  -0.539,-131.251,1434.710], [  -8.340, -26.663,  10.038] ] ),
# ( 279, [ [   9.361,-130.529,1433.801], [   8.391, -24.395,   8.129] ] ),
# ( 280, [ [  12.293,-108.004,1394.644], [ -20.525,  -2.858,  -4.989] ] ),
# ( 281, [ [  17.129, -66.221, 938.958], [   5.828,   1.212, -49.422] ] ),
# ( 282, [ [  11.498, -67.907, 806.052], [  19.146,   3.945, -66.035] ] ),
# ( 283, [ [   3.453, -80.187,1208.950], [  -6.956, -19.955,-272.404] ] ),
# ( 284, [ [  -8.381, -58.164,1429.721], [   3.055,  15.166, -88.114] ] ),
# ( 285, [ [   4.832,-136.024,1448.806], [  -9.797, -29.003,  -1.031] ] ),
# ( 286, [ [  25.584,-101.952,1409.817], [  11.866, -11.879, -16.532] ] ),
# ( 287, [ [   5.488, -47.520,1218.815], [ -13.011, 110.974,-117.481] ] ),
# ( 288, [ [  -0.497, -88.117,1055.124], [  -5.512,   0.035,-117.313] ] ),
# ( 289, [ [   2.251, -72.056, 934.224], [  -0.434,   4.350,  -6.767] ] ),
# ( 290, [ [  -3.681, -51.274,1437.257], [  29.986, -21.817,  62.139] ] ),
# ( 291, [ [  -0.628, -62.693,1481.347], [   3.263, -21.056, -37.149] ] ),
# ( 292, [ [  -6.224, -58.028,1472.469], [ -32.570,  -0.612,  56.461] ] ),
# ( 293, [ [  36.380, -80.588,1457.154], [  -0.613, -18.668, -20.229] ] ),
# ( 294, [ [   8.015,-102.434,1436.256], [  -2.443, -38.545,   4.918] ] ),
# ( 295, [ [  11.232,-116.557,1452.669], [  -2.323, -21.230,   7.172] ] ),
# ( 296, [ [  35.703, -82.095,1448.654], [  -6.075,   3.203, -19.179] ] ),
# ( 297, [ [  15.757, -54.707,1436.215], [  -1.959,  46.295, -86.490] ] ),
# ( 298, [ [  34.972, -83.394,1441.867], [   0.378,  -3.862, -25.655] ] ),
# ( 299, [ [  47.342, -89.596,1384.372], [  35.766, -30.305, -53.949] ] ),
# ( 300, [ [  20.440, -92.164,1438.379], [ -27.350, -19.027,   2.437] ] ),
# ( 301, [ [  14.327, -97.349,1437.414], [  -6.064,  -6.427,  -4.386] ] ),
# ( 302, [ [  14.027, -67.328,1469.790], [  -7.892, -15.789,   9.051] ] ),
# ( 303, [ [   6.677,-105.285,1436.500], [  -4.765, -21.028, -15.215] ] ),
# ( 304, [ [  13.106, -99.879,1430.559], [  -8.104,  -6.446, -20.042] ] ),
# ( 305, [ [  34.265, -84.451,1435.143], [  -1.385,  -7.670, -39.261] ] ),
# ( 306, [ [  36.305, -98.299,1367.884], [  29.613, -35.026, -60.996] ] ),
# ( 307, [ [  29.099, -94.640,1421.734], [ -38.123,  -9.292,   2.806] ] ),
# ( 308, [ [   4.377, -94.270,1418.331], [  -4.038, -14.153, -22.622] ] ),
# ( 309, [ [  32.290, -92.737,1396.093], [ -14.248,   2.936, -50.723] ] ),
# ( 310, [ [  17.234, -73.100,1379.748], [  16.620, -11.719, -41.967] ] ),
# ( 311, [ [  30.962, -99.435,1395.981], [ -15.002, -14.028,   2.755] ] ),
# ( 312, [ [  30.621,-102.116,1398.184], [ -11.786,  -7.944,  -9.212] ] ),
# ( 313, [ [  17.463, -71.807,1476.875], [   1.299,   0.601,   0.450] ] ),
# ( 314, [ [  15.985,-109.344,1318.216], [  -7.353,   7.298, -75.630] ] ),
# ( 315, [ [  17.739,-111.935,1303.563], [  -2.316,  19.120, -92.511] ] ),
# ( 316, [ [   6.058, -96.193,1168.221], [  -6.653,   6.353,-118.141] ] ),
# ( 317, [ [   4.126, -91.315,1143.098], [  -3.726,   6.361, -58.476] ] ),
# ( 318, [ [   4.490,-120.333, 982.412], [  -5.288,  48.187, -94.960] ] ),
# ( 319, [ [  -4.485, -76.724, 919.081], [   9.185,  52.286, -36.530] ] ),
# ( 320, [ [   5.830, -19.109, 858.554], [  44.392,   7.136, -47.818] ] ),
# ( 321, [ [   8.089, -99.441,1425.022], [   0.924,  -3.507,  -4.296] ] ),
# ( 322, [ [  30.626,-102.783,1394.211], [  -6.061,  -3.101,   2.663] ] ),
# ( 323, [ [  16.284, -83.248,1459.373], [  11.763,  -6.290, -18.355] ] ),
# ( 324, [ [   5.571, -54.485,1228.819], [   0.985,  -3.311, -23.523] ] ),
# ( 325, [ [  11.611, -66.853,1242.190], [  -1.580,  59.097,-202.489] ] ),
# ( 326, [ [   0.939, -37.534,1237.464], [   5.748, -20.681, -36.544] ] ),
# ( 327, [ [  17.284, -76.799,1444.923], [  14.501,  14.331,  10.708] ] ),
# ( 328, [ [  15.862, -66.876,1461.220], [ -11.043,  -0.118,   9.971] ] ),
# ( 329, [ [  14.086, -75.333,1461.325], [   2.338,  -5.241,  -3.046] ] ),
# ( 330, [ [   4.653, -62.902,1463.956], [  -7.282,   8.493,   9.078] ] ),
# ( 331, [ [  12.273, -59.511,1466.982], [   7.606,  -0.578,  -1.022] ] ),
# ( 332, [ [  17.160, -34.378, 806.448], [  27.543,  24.741,   3.359] ] ),
# ( 333, [ [  15.319, -62.224, 813.439], [   3.774, -11.411, -71.256] ] ),
# ( 334, [ [  12.660,-104.603,1398.421], [  -9.935,  -6.725, -48.078] ] ),
# ( 335, [ [  13.336, -56.091,1200.123], [  -6.295,  89.789,-110.937] ] ),
# ( 336, [ [   5.960,-114.822,1414.415], [  -7.322, -15.601, -18.775] ] ),
# ( 337, [ [  11.744,-108.913,1398.904], [ -22.373,  -4.705,   2.951] ] ),
# ( 338, [ [   2.887, -73.508,1401.050], [   4.801,   8.319,-333.069] ] ),
# ( 339, [ [  -0.076,-123.109,1042.699], [  10.000,  10.000,  10.000] ] ),
# ( 340, [ [ -11.729, -35.572,1218.308], [  -9.681, -10.065, -37.881] ] ),
# ( 341, [ [  -5.146, -60.037,1046.049], [   6.713, -18.634, -28.201] ] ),
# ( 342, [ [  -8.073, -66.896,1030.617], [  -9.770,  -3.013,  -8.410] ] ),
# ( 343, [ [   1.206, -66.887,1023.603], [  -1.322, -29.040, -12.503] ] ),
# ( 344, [ [  -8.962, -58.967,1014.734], [ -21.929, -23.536, -54.510] ] ),
# ( 345, [ [   7.763, -64.499,1014.489], [  -2.400,  -6.782,   3.410] ] ),
# ( 346, [ [ -17.307, -40.506,1204.246], [  -9.891, -13.565, -34.650] ] ),
# ( 347, [ [  -2.922, -33.223,1210.302], [  -5.806, -14.492, -11.204] ] ),
# ( 348, [ [  -2.755, -47.993,1204.295], [   2.207, -17.950,  25.580] ] ),
# ( 349, [ [  -7.127, -38.497,1187.670], [  -4.807,  -3.773,  13.463] ] ),
# ( 350, [ [   6.596, -34.950,1150.273], [   2.792,   9.205,  51.282] ] ),
# ( 351, [ [  -7.912, -36.723,1133.987], [  -2.150, -17.376,   5.268] ] ),
# ( 352, [ [  -1.400, -46.133,1113.427], [   7.303, -24.780,  64.259] ] ),
# ( 353, [ [   5.091, -45.029,1100.265], [   6.857, -10.876,   7.767] ] ),
# ( 354, [ [  13.820, -46.607, 875.366], [  -6.292, -13.759,  18.035] ] )
# ]
