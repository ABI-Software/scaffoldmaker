"""
Generates a 1-D mesh for nerve centrelines.
"""

from __future__ import division

import json
import os
import sys

from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldCoordinates,\
    findOrCreateFieldFiniteElement, findOrCreateFieldStoredString
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from opencmiss.zinc.element import Element, Elementbasis
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation import nerve_terms


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

        # markers group, coordinates and derivatives
        def read_centrelines(manifest_file):
            def read_json(filename):
                with open(filename, 'r') as f:
                    dicf = json.loads(f.read())
                return dicf

            # input_file = 'C:\flatmap\rat-flatmap\manifest.json'
            manifest_dict = read_json(os.path.join(manifest_file))
            file_dirname = os.path.dirname(manifest_file)
            properties_data = read_json(os.path.join(file_dirname, manifest_dict['properties']))
            anatomicalMap_data = read_json(os.path.join(file_dirname, manifest_dict['anatomicalMap']))
            return properties_data, anatomicalMap_data

        def get_anatomical_names(properties_data, anatomicalMap_data):
            d = {}
            for key in properties_data['features']:
                try:
                    d[key] = anatomicalMap_data[properties_data['features'][key]['class']]['name']
                except KeyError:
                    pass
            return d

        def get_anatomical_name(internal_name, terms_dict):
            try:
                return terms_dict[internal_name]
            except KeyError:
                return internal_name

        if input_file:
            properties_data, anatomicalMap_data = read_centrelines(input_file)

        def marker_record(name, x, d1, radius_size):
            anatomical_name = get_anatomical_name(name, d)
            return {"name": name, "group": nerve_terms.get_nerve_term(anatomical_name), "x": x, "d1": d1, "r": radius_size}

        def get_centrelines_list(properties_data):
            return properties_data['networks'][0]['centrelines']

        def get_points_in_centrelines(centrelines):
            import numpy as np

            def _get_point_in_centreline(centreline):
                pl = centreline['connects']
                try:
                    pl.extend(centreline['contained-in'])
                except KeyError:
                    pass
                pl = list(set(pl))
                pl.sort()
                return pl

            def _get_anatomical_name(points):
                # ("trigeminal nucleus", "UBERON:0002925")
                # ('point_1', ''),
                terms = []
                for key in points:
                    try:
                        terms.append((anatomicalMap_data[properties_data['features'][key]['class']]['name'],
                                      anatomicalMap_data[properties_data['features'][key]['class']]['term']))
                    except KeyError:
                        terms.append((key, ""))

                terms = list(set(terms))
                terms.sort()
                for term in terms:
                    print('{},'.format(term))
                return terms

            points = []
            for centreline in centrelines:
                points.extend(_get_point_in_centreline(centreline))
                points = list(set(points))
                points.sort()

            d = {}
            rng = np.random.default_rng(10)
            for p in points:
                d[p] = {"x": (1000*np.random.rand(3)).tolist(), "d1": [10.0, 10.0, 10.0], "r": 6.0}

            filename = os.path.join(sys.argv[2], 'coordintes.json')
            with open(filename, 'w') as f:
                json.dump(d, f, indent=4)

            with open(filename, 'r') as f:
                dic = json.load(f)

            _get_anatomical_name(points)
            a=1


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

        if input_file:
            d = get_anatomical_names(properties_data, anatomicalMap_data)
        else:
            d = names_internal_anatomical_map

        nerveMarkerPoints = [
            marker_record("brain_44-1", [-3.068, -58.888, 1485.941], [-5.129, -2.099, -5.383], 6.0),
            marker_record("brain_12-1", [0.526, -56.094, 1488.768], [-3.962, 3.013, 3.485], 6.0),
            marker_record("brain_34-1", [1.878, -56.261, 1489.3], [-3.416, 3.366, 3.014], 6.0),
            marker_record("brain_33-1", [2.607, -56.401, 1489.594], [-3.918, 3.78, 3.375], 6.0),
            marker_record("point_1", [-0.628, -62.693, 1481.347], [25.147, 6.953, 2.966], 6.0),
            marker_record("ardell_4_branching_point", [3.296, -68.465, 1338.18], [-0.341, -1.935, -12.395], 6.0),
            marker_record("ganglion_1-1", [19.14, -71.172, 1471.654], [3.23, -2.449, -7.906], 6.0),
            marker_record("bolser_X-1", [2.263, -68.435, 1345.906], [-55.339, 23.836, -191.394], 6.0),
            marker_record("ardell_3_branching_point", [0.501, -67.23, 1351.551], [-58.22, 25.91, -179.434], 6.0),
            marker_record("point_3", [17.463, -71.807, 1476.875], [-6.734, 0.139, 6.414], 6.0),
            marker_record("label_1-1", [19.059, -76.211, 1480.118], [-6.259, 0.228, 3.865], 6.0),
            marker_record("ganglion_2-1", [35.504, -78.359, 1463.385], [-1.057, -4.299, -8.305], 6.0),
            marker_record("point_6", [15.862, -66.876, 1461.22], [-16.113, -0.11, 1.188], 6.0),
            marker_record("point_7", [14.086, -75.333, 1461.325], [3.242, -4.337, -0.252], 6.0),
            marker_record("point_8", [4.653, -62.902, 1463.956], [-7.282, 8.493, 9.077], 6.0),
            marker_record("point_9", [12.273, -59.511, 1466.982], [7.606, -0.578, -1.022], 6.0),
            marker_record("point_10", [-6.224, -58.028, 1472.469], [-26.479, -6.517, 43.508], 6.0),
            marker_record("label_4-1", [60.913, -52.711, 1519.712], [113.397, 7.643, 25.844], 6.0),
            marker_record("label_3-1", [-57.275, -70.875, 1519.911], [-65.695, -23.698, 36.302], 6.0),
            marker_record("label_2-1", [-65.468, -70.375, 1520.633], [-80.971, -22.695, 37.571], 6.0),
            marker_record("point_11", [36.38, -80.588, 1457.154], [-0.118, -2.112, -6.641], 6.0),
            marker_record("point_13", [11.232, -116.557, 1452.669], [-17.74, -38.752, -6.999], 6.0),
            marker_record("label_6-1", [-14.0, -152.355, 1446.593], [-32.062, -32.179, -5.049], 6.0),
            marker_record("label_5-1", [9.582, -154.877, 1443.576], [12.923, -33.908, -10.012], 6.0),
            marker_record("point_14", [35.703, -82.095, 1448.654], [-6.773, -2.898, -16.891], 6.0),
            marker_record("cardio_4-1", [-3.468, -73.077, 1344.169], [-68.816, 21.122, -164.179], 6.0),
            marker_record("point_16", [34.972, -83.394, 1441.867], [0.902, 1.073, -8.921], 6.0),
            marker_record("point_18", [20.44, -92.164, 1438.379], [-5.057, -3.352, -3.807], 6.0),
            marker_record("point_19", [14.327, -97.349, 1437.414], [-5.856, -6.144, -4.017], 6.0),
            marker_record("point_20", [6.677, -105.285, 1436.5], [-6.059, -19.197, -11.926], 6.0),
            marker_record("label_7-1", [-0.539, -131.251, 1434.71], [-7.68, -30.025, 7.655], 6.0),
            marker_record("label_8-1", [9.361, -130.529, 1433.801], [9.952, -27.252, 5.685], 6.0),
            marker_record("respiratory_11-1", [1.823, -112.755, 1413.751], [-3.0, 3.499, -27.6], 6.0),
            marker_record("respiratory_5-1", [5.96, -114.822, 1414.415], [3.841, 0.102, -26.783], 6.0),
            marker_record("point_21", [13.106, -99.879, 1430.559], [-7.534, -5.737, -20.761], 6.0),
            marker_record("label_10-1", [2.098, -109.068, 1405.412], [-14.437, -12.602, -29.442], 6.0),
            marker_record("label_11-1", [0.755, -107.608, 1406.119], [-17.08, -9.671, -27.975], 6.0),
            marker_record("label_12-1", [0.392, -105.992, 1406.523], [-17.761, -6.441, -27.108], 6.0),
            marker_record("point_37", [12.367, -99.729, 1421.523], [1.399, 2.283, -5.055], 6.0),
            marker_record("bolser_X-5", [4.033, -107.383, 1404.8], [-10.611, -9.27, -30.754], 6.0),
            marker_record("bolser_X-4", [1.376, -109.143, 1406.395], [-15.831, -12.715, -27.404], 6.0),
            marker_record("external_laryngeal_n_branching_point", [19.268, -95.333, 1427.899], [-2.347, -4.675, -8.715],
                          6.0),
            marker_record("plexus_2-1", [4.832, -136.024, 1448.806], [-15.246, -52.924, 50.184], 6.0),
            marker_record("bolser_X-14", [7.465, -106.606, 1402.082], [-14.226, -9.153, -28.678], 6.0),
            marker_record("label_15-1", [6.818, -105.579, 1401.325], [-15.423, -7.075, -30.01], 6.0),
            marker_record("label_14-1", [6.623, -108.164, 1402.941], [-15.945, -12.281, -27.056], 6.0),
            marker_record("bolser_X-15", [8.188, -105.664, 1400.628], [-12.738, -7.253, -31.429], 6.0),
            marker_record("point_41", [19.826, -96.819, 1421.789], [0.141, 1.614, -5.649], 6.0),
            marker_record("point_22", [34.265, -84.451, 1435.143], [-3.332, -7.049, -20.25], 6.0),
            marker_record("point_24", [29.099, -94.64, 1421.734], [-3.226, -6.683, -10.8], 6.0),
            marker_record("point_25", [27.685, -97.76, 1415.52], [-1.761, -3.668, -5.999], 6.0),
            marker_record("plexus_3-1", [25.584, -101.952, 1409.817], [-2.883, -5.866, -8.726], 6.0),
            marker_record("point_26", [32.29, -92.737, 1396.093], [0.321, -5.982, -10.694], 6.0),
            marker_record("point_28", [30.962, -99.435, 1395.981], [-4.261, -2.824, 1.454], 6.0),
            marker_record("label_9-1", [12.293, -108.004, 1394.644], [-22.842, -7.903, -1.896], 6.0),
            marker_record("point_29", [30.621, -102.116, 1398.184], [-5.481, -5.016, 2.894], 6.0),
            marker_record("label_15-1", [15.735, -106.477, 1387.002], [-17.503, -0.767, -12.83], 6.0),
            marker_record("label_16-1", [16.942, -106.518, 1385.462], [-15.129, -0.845, -15.793], 6.0),
            marker_record("label_17-1", [18.305, -110.622, 1389.761], [-12.874, -9.06, -7.701], 6.0),
            marker_record("label_18-1", [17.812, -111.228, 1390.543], [-13.768, -10.203, -6.101], 6.0),
            marker_record("label_19-1", [20.305, -109.662, 1385.525], [-8.697, -6.995, -15.815], 6.0),
            marker_record("label_20-1", [19.183, -111.35, 1388.633], [-11.115, -10.505, -9.945], 6.0),
            marker_record("label_21-1", [18.88, -110.329, 1387.243], [-11.712, -8.465, -12.709], 6.0),
            marker_record("label_22-1", [19.598, -110.678, 1389.008], [-10.296, -9.173, -9.205], 6.0),
            marker_record("label_13-1", [12.808, -107.604, 1396.148], [-21.763, -7.091, 1.085], 6.0),
            marker_record("respiratory_8-1", [11.744, -108.913, 1398.904], [-23.743, -9.631, 6.505], 6.0),
            marker_record("respiratory_11-1", [12.66, -104.603, 1398.421], [-21.174, -1.106, 5.361], 6.0),
            marker_record("point_38", [30.626, -102.783, 1394.211], [-0.863, -1.635, -4.163], 6.0),
            marker_record("bolser_X-2", [13.545, -109.239, 1392.497], [-20.354, -10.352, -6.156], 6.0),
            marker_record("bolser_X-3", [13.402, -106.487, 1391.751], [-20.209, -4.791, -7.477], 6.0),
            marker_record("point_30", [28.921, -109.629, 1316.751], [-3.914, -3.308, -13.084], 6.0),
            marker_record("plexus_4-1", [20.571, -47.974, 1214.357], [-11.614, 110.759, -121.987], 6.0),
            marker_record("point_31", [28.321, -112.239, 1305.054], [0.365, -9.088, -15.532], 6.0),
            marker_record("respiratory_17-1", [13.336, -56.091, 1200.123], [-4.754, 7.795, -8.567], 6.0),
            marker_record("point_32", [6.058, -96.193, 1168.221], [-7.486, 1.904, -38.19], 6.0),
            marker_record("point_33", [4.126, -91.315, 1143.098], [-6.426, 7.025, -34.377], 6.0),
            marker_record("plexus_5-1", [-0.497, -88.117, 1055.124], [-5.067, -1.88, -108.785], 6.0),
            marker_record("digestive_12-1", [-0.849, -78.626, 1059.708], [-5.782, 17.092, -99.868], 6.0),
            marker_record("point_34", [4.49, -120.333, 982.412], [-5.29, 48.18, -94.925], 6.0),
            marker_record("digestive_9-1", [4.051, -80.134, 952.868], [-0.119, 8.718, 0.525], 6.0),
            marker_record("point_35", [-4.485, -76.724, 919.081], [2.729, 53.763, -37.471], 6.0),
            marker_record("plexus_6-1", [2.251, -72.056, 934.224], [1.907, -1.596, 6.514], 6.0),
            marker_record("point_36", [-9.096, -18.384, 858.438], [73.681, 23.211, -65.919], 6.0),
            marker_record("urinary_13-1", [84.559, -64.284, 856.328], [82.082, -83.081, 44.569], 6.0),
        ]
        nodes_w_ver = []

        def coordinates_modified(nerveMarkerPoints):

            filename = os.path.join(sys.argv[2], 'coordintes.json')
            filename2 = os.path.join(sys.argv[2], 'coordintes_modified.json')
            with open(filename, 'r') as f:
                dic = json.load(f)

            for pd in nerveMarkerPoints:
                if pd['name'] in dic:
                    dic[pd['name']]['x'] = pd['x']
                    dic[pd['name']]['d1'] = pd['d1']
                    dic[pd['name']]['r'] = pd['r']
            with open(filename2, 'w') as f:
                json.dump(dic, f, indent=4)
            a=1

        cl = get_centrelines_list(properties_data)
        get_points_in_centrelines(cl)
        coordinates_modified(nerveMarkerPoints)

        annotationGroups = []
        marker_node_dict = {}
        # create nodes
        nodeIdentifier = 1
        for n3 in range(len(nerveMarkerPoints)):
            if n3+1 in nodes_w_ver:
                nt = nodetemplate1
            else:
                nt = nodetemplate
            node = nodes.createNode(nodeIdentifier, nt)
            cache.setNode(node)

            marker_node_dict[nerveMarkerPoints[n3]["name"]] = nodeIdentifier
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nerveMarkerPoints[n3]["x"])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nerveMarkerPoints[n3]["d1"])
            markerName.assignString(cache, nerveMarkerPoints[n3]["group"][0])
            annotationGroup = AnnotationGroup(region, nerveMarkerPoints[n3]["group"])
            annotationGroups.append(annotationGroup)
            annotationGroup.getNodesetGroup(nodes).addNode(node)
            if n3+1 in nodes_w_ver:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, nerveMarkerPoints[n3]["d1"])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 3, nerveMarkerPoints[n3]["d1"])

            radius.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nerveMarkerPoints[n3]["r"])
            nodeIdentifier = nodeIdentifier + 1

        # connectivity
        centrelines = [
            {'id': 'n_1', 'connects': ['brain_44-1', 'point_1']},
            {'id': 'n_5', 'connects': ['point_1', 'ganglion_1-1']},
            {'id': 'n_5', 'connects': ['ganglion_1-1', 'ganglion_2-1']},
            {'id': 'n_5', 'connects': ['ganglion_2-1', 'point_11']},
            {'id': 'n_5', 'connects': ['point_11', 'point_14']},
            {'id': 'n_5', 'connects': ['point_14', 'point_16']},
            {'id': 'n_5', 'connects': ['point_16', 'point_22']},
            {'id': 'n_5', 'connects': ['point_22', 'point_26']},
            {'id': 'n_5', 'connects': ['point_26', 'point_30']},
            {'id': 'n_5', 'connects': ['point_30', 'point_31']},
            {'id': 'n_5', 'connects': ['point_31', 'point_32']},
            {'id': 'n_5', 'connects': ['point_32', 'point_34']},
            {'id': 'n_5', 'connects': ['point_34', 'point_35']},
            {'id': 'n_5', 'connects': ['point_35', 'point_36']},
            {'id': 'n_68', 'connects': ['ganglion_1-1', 'point_3']},
            {'id': 'n_6', 'connects': ['ganglion_1-1', 'label_1-1']},
            {'id': 'n_7', 'connects': ['ganglion_2-1', 'point_6']},
            {'id': 'n_7', 'connects': ['point_6', 'point_8']},
            {'id': 'n_7', 'connects': ['point_8', 'point_10']},
            {'id': 'n_12', 'connects': ['point_6', 'point_7']},
            {'id': 'n_8', 'connects': ['point_8', 'point_9']},
            {'id': 'n_9', 'connects': ['point_10', 'label_3-1']},
            {'id': 'n_11', 'connects': ['point_10', 'label_4-1']},
            {'id': 'n_10', 'connects': ['point_10', 'label_2-1']},
            {'id': 'n_13', 'connects': ['point_11', 'point_13']},
            {'id': 'n_2', 'connects': ['point_1', 'brain_12-1']},
            {'id': 'n_15', 'connects': ['point_13', 'label_6-1']},
            {'id': 'n_14', 'connects': ['point_13', 'label_5-1']},
            {'id': 'n_16', 'connects': ['point_14', 'cardio_4-1']},
            {'id': 'n_3', 'connects': ['point_1', 'brain_34-1']},
            {'id': 'n_17', 'connects': ['point_16', 'point_18']},
            {'id': 'n_4', 'connects': ['point_1', 'brain_33-1']},
            {'id': 'n_18', 'connects': ['point_18', 'point_19']},
            {'id': 'n_29', 'connects': ['point_18', 'external_laryngeal_n_branching_point']},
            {'id': 'n_24', 'connects': ['point_19', 'point_21']},
            {'id': 'n_19', 'connects': ['point_19', 'point_20']},
            {'id': 'n_22', 'connects': ['point_20', 'label_7-1']},
            {'id': 'n_20', 'connects': ['point_20', 'label_8-1']},
            {'id': 'n_23', 'connects': ['point_20', 'respiratory_5-1']},
            {'id': 'n_21', 'connects': ['point_20', 'respiratory_11-1']},
            {'id': 'n_27', 'connects': ['point_21', 'label_12-1']},
            {'id': 'bolser_5', 'connects': ['point_21', 'bolser_X-5']},
            {'id': 'n_26', 'connects': ['point_21', 'label_11-1']},
            {'id': 'bolser_4', 'connects': ['point_21', 'bolser_X-4']},
            {'id': 'n_28', 'connects': ['point_21', 'point_37']},
            {'id': 'n_67', 'connects': ['external_laryngeal_n_branching_point', 'plexus_2-1']},
            {'id': 'n_68', 'connects': ['external_laryngeal_n_branching_point', 'point_41']},
            {'id': 'bolser_15', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-15']},
            {'id': 'bolser_14', 'connects': ['external_laryngeal_n_branching_point', 'bolser_X-14']},
            {'id': 'n_65', 'connects': ['external_laryngeal_n_branching_point', 'label_14-1']},
            {'id': 'n_66', 'connects': ['external_laryngeal_n_branching_point', 'label_15-1']},
            {'id': 'n_30', 'connects': ['point_22', 'point_24']},
            {'id': 'n_25', 'connects': ['point_21', 'label_10-1']},
            {'id': 'n_32', 'connects': ['point_24', 'point_25']},
            {'id': 'n_34', 'connects': ['point_25', 'plexus_3-1']},
            {'id': 'n_35', 'connects': ['point_24', 'plexus_3-1']},
            {'id': 'n_36', 'connects': ['point_26', 'point_28']},
            {'id': 'n_50', 'connects': ['point_28', 'point_38']},
            {'id': 'n_38', 'connects': ['point_28', 'point_29']},
            {'id': 'n_47', 'connects': ['point_28', 'label_13-1']},
            {'id': 'bolser_2', 'connects': ['point_28', 'bolser_X-2']},
            {'id': 'n_49', 'connects': ['point_28', 'respiratory_11-1']},
            {'id': 'n_48', 'connects': ['point_28', 'respiratory_8-1']},
            {'id': 'n_37', 'connects': ['point_28', 'label_9-1']},
            {'id': 'bolser_3', 'connects': ['point_28', 'bolser_X-3']},
            {'id': 'n_42', 'connects': ['point_29', 'label_18-1']},
            {'id': 'n_43', 'connects': ['point_29', 'label_19-1']},
            {'id': 'n_40', 'connects': ['point_29', 'label_16-1']},
            {'id': 'n_46', 'connects': ['point_29', 'label_22-1']},
            {'id': 'n_44', 'connects': ['point_29', 'label_20-1']},
            {'id': 'n_45', 'connects': ['point_29', 'label_21-1']},
            {'id': 'n_41', 'connects': ['point_29', 'label_17-1']},
            {'id': 'bolser_1', 'connects': ['ganglion_1-1', 'bolser_X-1']},
            {'id': 'ardell_3', 'connects': ['ganglion_1-1', 'ardell_3_branching_point']},
            {'id': 'ardell_4', 'connects': ['point_1', 'ardell_4_branching_point']},
            {'id': 'n_51', 'connects': ['point_30', 'plexus_4-1']},
            {'id': 'n_52', 'connects': ['point_31', 'respiratory_17-1']},
            {'id': 'n_53', 'connects': ['point_32', 'point_33']},
            {'id': 'n_54', 'connects': ['point_33', 'digestive_12-1']},
            {'id': 'n_55', 'connects': ['point_33', 'plexus_5-1']},
            {'id': 'n_56', 'connects': ['point_34', 'digestive_9-1']},
            {'id': 'n_57', 'connects': ['point_35', 'plexus_6-1']},
            {'id': 'n_58', 'connects': ['point_36', 'urinary_13-1']},
            {'id': 'n_39', 'connects': ['point_29', 'label_15-1']},
        ]

        elementIdentifier = 1
        for e3 in range(len(centrelines)):
            elementtemplate = elementtemplate1
            eft = eft1
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result = element.setNodesByIdentifier(eft, [marker_node_dict[centrelines[e3]['connects'][0]], marker_node_dict[centrelines[e3]['connects'][1]]])
            result = element.setNodesByIdentifier(eftRadius, [marker_node_dict[centrelines[e3]['connects'][0]], marker_node_dict[centrelines[e3]['connects'][1]]])
            annotationGroup = AnnotationGroup(region, nerve_terms.get_nerve_term(centrelines[e3]['id']))
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
