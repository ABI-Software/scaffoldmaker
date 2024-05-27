"""
Constructs a 1-D network layout mesh with specifiable structure.
"""
from cmlibs.maths.vectorops import magnitude, mult
from cmlibs.utils.zinc.field import Field, findOrCreateFieldCoordinates, findOrCreateFieldStoredString,\
    findOrCreateFieldGroup
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.scene import scene_get_selection_group
from cmlibs.zinc.field import Field, FieldGroup
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.interpolation import smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.zinc_utils import clearRegion, get_nodeset_field_parameters, \
    get_nodeset_path_ordered_field_parameters, make_nodeset_derivatives_orthogonal, \
    set_nodeset_field_parameters
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    findAnnotationGroupByName, getAnnotationGroupForTerm
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from enum import Enum


class MeshType_1d_vagus_path1(Scaffold_base):
    """
    Defines branching network layout with side dimensions.
    """

    parameterSetStructureStrings = {
        "Default": "1-2",
        "Left vagus nerve 1": "1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19.1-20-21-22-23-24-25-26-27-28-29-30-31-32-33-34.1, 34.2-35-36, 19.2-37-38.1, 38.2-39-40, 38.3-41-42, 34.3-43.1-44, 43.2-45-46",
        "Left vagus nerve 2":
            "1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19.1-20-21-22-23-24-25-26-27-28-29-30-31-32-33-34.1, 34.2-35-36, 19.2-37-38.1, 38.2-39-40, 38.3-41-42, 34.3-43.1-44, 43.2-45-46, 8.2-47-48.1-49-50.1-53.1-56-57, 48.2-51-52, 50.2-54-55, 53.2-58-59"
    }

    @classmethod
    def getName(cls):
        return "1D Vagus Path 1"

    @classmethod
    def getParameterSetNames(cls):
        return list(cls.parameterSetStructureStrings.keys())

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = parameterSetName
        options["Structure"] = cls.parameterSetStructureStrings[parameterSetName]
        options["Define inner coordinates"] = False  # can be overridden by parent scaffold
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            #  "Base parameter set"  # Hidden.
            #  "Structure"  # Hidden so must edit via interactive function.
            #  "Define inner coordinates"  # Hidden as enabled by parent scaffold.
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
        parameterSetName = options['Base parameter set']
        structure = options["Structure"]
        defineInnerCoordinates = options["Define inner coordinates"]
        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        fieldcache = fieldmodule.createFieldcache()

        outerPoints = []
        bifurcating_nodes = []
        nodesByGroup = []

        if "Default" in parameterSetName:
            return [], []

        elif "Left vagus nerve 1" in parameterSetName:
            outerPoints = [[[10.900,-70.490,1499.510], [3.058,1.607,-2.565], [0.350,1.950,0.280], [1.530,-0.360,-0.370], [0.890,-0.410,1.750], [-2.080,1.850,1.180]],
                            [[15.210,-68.890,1497.320], [5.470,1.545,-1.738], [1.490,1.310,0.190], [0.760,-0.910,0.200], [-0.350,0.110,1.970], [-0.400,-0.820,-0.740]],
                            [[21.640,-67.580,1496.360], [7.078,-1.968,-2.119], [1.490,1.310,0.190], [0.760,-0.910,0.200], [-0.350,0.110,1.970], [-0.400,-0.820,-0.740]],
                            [[27.230,-73.040,1493.470], [4.913,-3.320,-5.768], [1.880,0.140,0.670], [0.210,-0.590,0.150], [0.040,-1.970,0.320], [0.170,-1.010,-0.750]],
                            [[30.430,-73.620,1486.260], [3.272,-1.026,-8.446], [1.880,0.140,0.670], [0.210,-0.590,0.150], [0.040,-1.970,0.320], [0.170,-1.010,-0.750]],
                            [[33.620,-75.160,1476.690], [2.810,-2.377,-10.701], [1.930,0.110,0.500], [0.060,-0.040,-0.220], [0.000,-1.960,0.420], [-0.020,0.020,0.100]],
                            [[35.880,-78.480,1464.950], [1.174,-2.874,-10.412], [1.990,0.060,0.230], [0.030,-0.060,-0.240], [-0.000,-1.930,0.530], [-0.000,0.010,0.030]],
                            [[36.190,-80.900,1455.970], [0.077,-2.012,-7.865], [2.000,0.000,0.020], [0.000,-0.040,-0.150], [-0.010,-1.940,0.490], [0.000,-0.010,-0.030]],
                            [[36.090,-82.530,1449.240], [-0.273,-1.468,-5.977], [2.000,-0.020,-0.090], [-0.010,-0.030,-0.140], [0.000,-1.940,0.480], [0.000,-0.010,-0.050]],
                            [[35.690,-83.830,1444.020], [-0.459,-0.978,-4.847], [1.980,-0.050,-0.260], [-0.010,-0.010,-0.130], [0.000,-1.960,0.410], [-0.000,-0.020,-0.110]],
                            [[35.190,-84.520,1439.570], [-1.535,-1.266,-9.722], [1.970,-0.050,-0.360], [-0.010,0.000,-0.090], [-0.000,-1.980,0.270], [-0.000,-0.020,-0.140]],
                            [[32.150,-86.090,1424.650], [-3.898,-1.195,-17.097], [1.960,-0.030,-0.410], [-0.010,0.010,-0.040], [-0.000,-2.000,0.140], [0.000,-0.010,-0.040]],
                            [[27.300,-86.750,1405.430], [-4.991,-2.075,-21.658], [1.950,-0.040,-0.440], [0.010,-0.000,0.100], [0.000,-1.990,0.190], [0.000,0.020,0.130]],
                            [[22.310,-90.520,1381.490], [-1.732,-6.353,-25.914], [2.000,-0.030,-0.130], [0.000,0.100,0.430], [0.000,-1.950,0.460], [-0.010,0.070,0.280]],
                            [[24.420,-99.580,1354.430], [3.271,-8.211,-23.148], [1.950,0.160,0.430], [-0.010,0.050,0.170], [-0.010,-1.850,0.760], [-0.000,0.020,0.060]],
                            [[28.380,-106.780,1335.240], [2.320,-5.265,-15.766], [1.980,0.090,0.280], [0.020,-0.080,-0.210], [-0.000,-1.900,0.630], [0.010,-0.050,-0.170]],
                            [[29.470,-110.340,1323.140], [0.149,-3.408,-15.340], [2.000,0.010,0.040], [0.010,-0.030,-0.220], [0.000,-1.950,0.440], [0.020,-0.050,-0.250]],
                            [[28.210,-113.130,1304.780], [-3.637,2.138,-20.620], [2.000,0.050,-0.140], [-0.030,0.120,-0.210], [0.060,-2.000,0.100], [-0.040,0.250,-0.880]],
                            [[22.080,-105.130,1283.990], [[-6.244,13.100,-18.015],[-0.718,13.905,-13.121]],
                             [[1.940,0.280,-0.400],[1.940,0.280,-0.400]],
                             [[-0.020,0.110,-0.080],[-0.020,0.110,-0.080]],
                             [[-0.110,-1.360,-1.460],[-0.110,-1.360,-1.460]],
                             [[-0.030,0.250,-0.700],[-0.030,0.250,-0.700]]],
                            [[16.200,-88.020,1270.070], [-4.734,14.710,-14.740], [1.960,0.270,-0.300], [0.030,-0.100,0.110], [-0.000,-1.480,-1.340], [0.060,-0.250,0.390]],
                            [[12.590,-75.800,1255.040], [-2.455,9.070,-17.173], [1.960,0.270,-0.300], [0.030,-0.100,0.110], [-0.000,-1.480,-1.340], [0.060,-0.250,0.390]],
                            [[11.390,-70.270,1236.640], [-0.968,2.075,-17.928], [1.960,0.270,-0.300], [0.030,-0.100,0.110], [-0.000,-1.480,-1.340], [0.060,-0.250,0.390]],
                            [[10.670,-71.290,1219.990], [-0.754,-4.029,-18.499], [2.000,-0.020,-0.100], [0.010,-0.010,0.120], [-0.000,-1.970,0.370], [0.010,0.030,0.750]],
                            [[9.920,-78.850,1200.330], [-0.604,-7.415,-17.005], [2.000,0.050,0.060], [-0.030,-0.040,-0.200], [0.020,-1.810,0.850], [-0.000,-0.010,-0.050]],
                            [[9.450,-85.880,1185.970], [-1.682,-5.454,-14.731], [2.000,0.050,0.060], [-0.030,-0.040,-0.200], [0.020,-1.810,0.850], [-0.000,-0.010,-0.050]],
                            [[6.620,-89.700,1171.150], [-3.136,-3.360,-15.171], [1.940,-0.080,-0.460], [-0.000,0.040,0.080], [-0.010,-1.970,0.320], [-0.010,0.040,0.080]],
                            [[3.180,-92.570,1155.680], [-1.953,-4.067,-14.129], [1.940,-0.080,-0.460], [-0.000,0.040,0.080], [-0.010,-1.970,0.320], [-0.010,0.040,0.080]],
                            [[2.540,-97.510,1143.140], [0.861,-6.135,-11.183], [1.990,0.120,0.210], [-0.140,0.070,0.790], [0.000,-1.740,0.990], [0.000,-0.020,-0.100]],
                            [[4.700,-104.460,1133.610], [3.004,-6.082,-9.175], [1.990,0.120,0.210], [-0.140,0.070,0.790], [0.000,-1.740,0.990], [0.000,-0.020,-0.100]],
                            [[8.400,-109.640,1124.940], [5.234,-3.868,-9.301], [1.990,0.120,0.210], [-0.140,0.070,0.790], [0.000,-1.740,0.990], [0.000,-0.020,-0.100]],
                            [[15.040,-111.850,1115.450], [7.442,1.216,-8.984], [1.660,0.070,1.110], [-0.230,-0.500,0.360], [-0.000,-2.000,0.130], [0.000,0.100,-1.090]],
                            [[22.470,-107.440,1107.910], [6.570,4.713,-6.419], [1.490,-0.740,1.100], [0.010,-0.390,-0.160], [0.000,-1.660,-1.120], [-0.000,0.280,-0.710]],
                            [[28.110,-102.670,1102.620], [5.393,3.779,-5.902], [1.490,-0.740,1.100], [0.010,-0.390,-0.160], [0.000,-1.660,-1.120], [-0.000,0.280,-0.710]],
                            [[33.150,-99.910,1096.270], [[4.630,1.720,-6.715],[9.276,20.368,-11.397],[8.794,-0.387,-13.708]],
                             [[1.630,-0.810,0.820],[1.630,-0.810,0.820],[1.630,-0.810,0.820]],
                             [[0.010,0.610,0.020],[0.010,0.610,0.020],[0.010,0.610,0.020]],
                             [[-0.000,-1.430,-1.400],[-0.000,-1.430,-1.400],[-0.000,-1.430,-1.400]],
                             [[-0.000,-0.020,0.890],[-0.000,-0.020,0.890],[-0.000,-0.020,0.890]]],
                            [[45.239,-93.595,1085.478], [7.044,-0.868,-7.462], [2.026,-0.129,1.495], [0.100,0.770,-0.270], [-0.000,-1.760,0.940], [0.000,0.350,1.570]],
                            [[54.360,-102.450,1077.350], [1.669,-10.067,-4.216], [1.790,0.810,0.350], [0.550,-0.430,-1.430], [0.000,-0.800,1.830], [0.000,1.570,0.210]],
                            [[23.330,-88.250,1275.900], [3.463,14.329,-4.449], [1.980,-0.050,-0.260], [-0.010,-0.010,-0.130], [0.000,-1.960,0.410], [-0.000,-0.020,-0.110]],
                            [[32.287,-79.285,1268.902], [[11.678,3.730,-9.609],[5.310,-2.835,-11.763],[7.040,5.114,-10.458]],
                             [[1.970,-0.050,-0.360],[1.970,-0.050,-0.360],[1.970,-0.050,-0.360]],
                             [[-0.010,0.000,-0.090],[-0.010,0.000,-0.090],[-0.010,0.000,-0.090]],
                             [[-0.000,-1.980,0.270],[-0.000,-1.980,0.270],[-0.000,-1.980,0.270]],
                             [[-0.000,-0.020,-0.140],[-0.000,-0.020,-0.140],[-0.000,-0.020,-0.140]]],
                            [[39.200,-80.300,1257.430], [9.220,1.186,-9.464], [1.960,-0.030,-0.410], [-0.010,0.010,-0.040], [-0.000,-2.000,0.140], [0.000,-0.010,-0.040]],
                            [[49.150,-77.930,1247.120], [11.111,4.117,-9.268], [1.960,-0.030,-0.410], [-0.010,0.010,-0.040], [-0.000,-2.000,0.140], [0.000,-0.010,-0.040]],
                            [[42.582,-77.775,1260.806], [5.396,-1.511,-1.137], [1.950,-0.040,-0.440], [0.010,-0.000,0.100], [0.000,-1.990,0.190], [0.000,0.020,0.130]],
                            [[52.020,-83.330,1261.290], [8.778,-8.775,5.301], [1.950,-0.040,-0.440], [0.010,-0.000,0.100], [0.000,-1.990,0.190], [0.000,0.020,0.130]],
                            [[39.570, -106.360, 1081.580], [[1.676, -6.540, -9.622], [-4.285, -3.707, -2.507]],
                             [[1.470, 0.640, 1.200], [1.470, 0.640, 1.200]],
                             [[0.100, 0.770, -0.270], [0.100, 0.770, -0.270]],
                             [[-0.000, -1.760, 0.940], [-0.000, -1.760, 0.940]],
                             [[0.000, 0.350, 1.570], [0.000, 0.350, 1.570]]],
                            [[45.981, -110.526, 1067.815], [8.047, 0.871, -12.742], [1.790, 0.810, 0.350], [0.550, -0.430, -1.430], [0.000, -0.800, 1.830], [0.000, 1.570, 0.210]],
                            [[39.021, -111.664, 1076.034], [1.344, -5.047, -5.083], [1.790, 0.810, 0.350], [0.550, -0.430, -1.430], [0.000, -0.800, 1.830], [0.000, 1.570, 0.210]],
                            [[38.728, -119.318, 1072.496], [-4.649, -7.186, 1.582], [1.790, 0.810, 0.350], [0.550, -0.430, -1.430], [0.000, -0.800, 1.830], [0.000, 1.570, 0.210]]]

            bifurcating_nodes = [6, 19, 34, 38, 43]
            nodesByGroup = [
                ['left vagus X nerve trunk', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 43, 44]],
                ['left A lateral pulmonary branch of vagus nerve', [43, 45, 46]],
                ['left A medial pulmonary branch of vagus nerve', [34, 35, 36]],
                ['left A thoracic cardiac branch of vagus nerve', [19, 37, 38, 41, 42]],
                ['left B thoracic cardiac branch of vagus nerve', [38, 39, 40]]
            ]

        elif "Left vagus nerve 2" in parameterSetName:
            outerPoints = [
                [[10.900, -70.490, 1499.510], [3.058, 1.607, -2.565], [-0.247, 1.802, 0.834], [1.873, 0.320, 0.299],
                 [1.389, -0.447, 1.376], [-0.225, -1.001, -0.156]],
                [[15.210, -68.890, 1497.320], [5.470, 1.545, -1.738], [-0.061, 1.581, 1.212], [-1.184, -0.578, 0.262],
                 [0.782, -1.104, 1.479], [-0.545, -0.190, 0.187]],
                [[21.640, -67.580, 1496.360], [7.078, -1.968, -2.119], [0.661, 1.803, 0.535], [3.760, -1.579, -0.395],
                 [0.364, -0.682, 1.849], [-0.737, 0.312, 0.033]],
                [[27.230, -73.040, 1493.470], [4.913, -3.320, -5.768], [1.600, 0.396, 1.135], [0.641, -0.863, -0.326],
                 [-0.179, -1.785, 0.875], [1.514, -1.903, -3.041]],
                [[30.430, -73.620, 1486.260], [3.265, -1.024, -8.429], [1.867, 0.144, 0.706], [-0.049, -0.210, 0.247],
                 [0.054, -1.978, 0.261], [-0.135, -0.059, 0.837]],
                [[33.608, -75.101, 1476.718], [2.811, -2.378, -10.704], [1.934, 0.107, 0.484], [0.083, -0.063, -0.324],
                 [-0.000, -1.957, 0.435], [-0.105, 0.097, 0.428]],
                [[35.880, -78.480, 1464.950], [1.176, -2.880, -10.433], [1.992, 0.054, 0.210], [0.052, -0.111, -0.429],
                 [-0.004, -1.929, 0.532], [0.019, -0.040, -0.136]],
                [[36.190, -80.900, 1455.970], [[0.077, -2.012, -7.865], [3.280, -5.664, 2.587]],
                 [[2.000, -0.000, 0.020], [3.409, 1.290, -1.497]], [[0.002, 0.018, 0.060], [-3.804, -3.798, -3.086]],
                 [[-0.005, -1.939, 0.496], [0.371, 0.990, 1.698]], [[0.003, -0.021, -0.087], [2.308, -0.640, 0.865]]],
                [[36.090, -82.530, 1449.240], [-0.273, -1.468, -5.977], [2.000, -0.019, -0.087],
                 [-0.008, -0.027, -0.138], [0.002, -1.941, 0.477], [0.012, 0.033, 0.131]],
                [[35.690, -83.830, 1444.020], [-0.459, -0.978, -4.847], [1.989, -0.034, -0.181], [0.029, 0.082, 0.365],
                 [0.002, -1.963, 0.396], [-0.006, -0.034, -0.126]],
                [[35.190, -84.520, 1439.570], [-1.535, -1.266, -9.722], [1.979, -0.043, -0.307],
                 [-0.028, 0.017, -0.106], [-0.003, -1.982, 0.259], [0.023, 0.013, 0.079]],
                [[32.150, -86.090, 1424.650], [-3.898, -1.195, -17.096], [1.953, -0.032, -0.443],
                 [-0.049, 0.002, -0.203], [-0.001, -2.000, 0.140], [-0.090, -0.017, -0.393]],
                [[27.300, -86.750, 1405.430], [-4.991, -2.075, -21.658], [1.949, -0.040, -0.445],
                 [-0.035, -0.023, -0.097], [0.003, -1.990, 0.190], [0.037, 0.026, 0.326]],
                [[22.310, -90.520, 1381.490], [-1.732, -6.353, -25.914], [2.000, -0.029, -0.127], [0.063, 0.229, 0.956],
                 [0.002, -1.946, 0.477], [0.011, 0.134, 0.524]],
                [[24.420, -99.580, 1354.430], [3.271, -8.211, -23.148], [1.986, 0.095, 0.247], [-0.038, 0.083, 0.273],
                 [0.007, -1.885, 0.670], [0.003, 0.013, 0.042]],
                [[28.380, -106.780, 1335.240], [2.320, -5.265, -15.766], [1.983, 0.084, 0.264], [0.065, -0.125, -0.302],
                 [-0.004, -1.899, 0.634], [0.012, -0.041, -0.110]],
                [[29.470, -110.340, 1323.140], [0.149, -3.408, -15.340], [2.000, 0.005, 0.018], [-0.005, 0.043, -0.169],
                 [0.001, -1.951, 0.433], [-0.008, 0.087, 0.200]],
                [[28.210, -113.130, 1304.780], [-3.637, 2.138, -20.620], [1.975, 0.071, -0.341],
                 [-0.118, 0.196, -0.597], [0.035, -1.992, -0.213], [-0.119, 0.144, -1.463]],
                [[22.080, -105.130, 1283.990], [[-6.244, 13.100, -18.015], [-0.787, 15.239, -14.380]],
                 [[1.926, 0.310, -0.442], [1.997, -0.015, -0.125]], [[0.019, 0.061, 0.122], [-0.115, 1.147, -1.434]],
                 [[-0.009, -1.618, -1.173], [-0.101, -1.372, -1.449]],
                 [[-0.356, 1.096, -1.336], [-0.325, 1.866, -1.521]]],
                [[16.200, -88.020, 1270.070], [-4.734, 14.710, -14.740], [1.951, 0.299, -0.329], [0.042, -0.099, 0.151],
                 [-0.020, -1.416, -1.407], [0.238, -0.695, 0.535]],
                [[12.590, -75.800, 1255.040], [-2.455, 9.070, -17.173], [1.983, 0.205, -0.175], [0.037, -0.121, 0.312],
                 [0.098, -1.757, -0.942], [0.173, -0.477, 0.907]],
                [[11.390, -70.270, 1236.640], [-0.968, 2.075, -17.928], [1.984, 0.246, -0.079], [0.009, -0.087, -0.156],
                 [0.234, -1.968, -0.240], [-0.049, -0.106, 0.616]],
                [[10.670, -71.290, 1219.990], [-0.754, -4.029, -18.499], [2.001, -0.015, -0.078],
                 [0.000, -0.119, 0.029], [0.002, -1.959, 0.427], [-0.098, 0.191, 0.724]],
                [[9.920, -78.850, 1200.330], [-0.604, -7.415, -17.005], [2.000, -0.006, -0.068], [-0.002, 0.094, 0.136],
                 [0.022, -1.833, 0.799], [-0.005, 0.129, 0.270]],
                [[9.450, -85.880, 1185.970], [-1.682, -5.454, -14.731], [1.990, -0.049, -0.209],
                 [-0.061, -0.188, -0.541], [0.026, -1.876, 0.692], [-0.053, -0.316, -0.801]],
                [[6.620, -89.700, 1171.150], [-3.136, -3.360, -15.171], [1.956, -0.064, -0.390],
                 [-0.031, 0.007, -0.113], [0.021, -1.949, 0.427], [-0.031, -0.035, -0.216]],
                [[3.180, -92.570, 1155.680], [-1.953, -4.067, -14.129], [1.977, -0.024, -0.266], [0.079, 0.101, 0.466],
                 [0.050, -1.919, 0.545], [0.062, 0.145, 0.585]],
                [[2.540, -97.510, 1143.140], [0.861, -6.135, -11.183], [2.000, 0.069, 0.116], [-0.032, 0.241, 0.504],
                 [0.005, -1.755, 0.963], [-0.078, 0.437, 0.751]],
                [[4.700, -104.460, 1133.610], [3.004, -6.082, -9.175], [1.934, 0.278, 0.449], [-0.034, -0.004, 0.155],
                 [-0.016, -1.671, 1.103], [0.100, -0.212, -0.274]],
                [[8.400, -109.640, 1124.940], [5.234, -3.868, -9.301], [1.779, 0.385, 0.841], [-0.254, -0.049, 0.552],
                 [0.029, -1.843, 0.783], [0.033, -0.182, -0.465]],
                [[15.040, -111.850, 1115.450], [7.442, 1.216, -8.984], [1.535, 0.048, 1.278], [-0.276, -0.663, 0.331],
                 [0.170, -1.993, -0.129], [0.790, 0.202, -1.871]],
                [[22.470, -107.440, 1107.910], [6.570, 4.713, -6.419], [1.538, -0.707, 1.055], [0.011, -0.268, -0.210],
                 [0.042, -1.634, -1.157], [-0.130, 0.248, -0.355]],
                [[28.110, -102.670, 1102.620], [5.393, 3.779, -5.902], [1.580, -0.681, 1.008], [0.043, -0.051, -0.100],
                 [-0.024, -1.676, -1.095], [-0.672, -0.430, 0.579]],
                [[33.150, -99.910, 1096.270],
                 [[4.630, 1.720, -6.715], [7.150, 15.699, -8.784], [9.921, -0.437, -15.464]],
                 [[1.588, -0.827, 0.883], [1.830, -0.475, 0.641], [1.525, -0.809, 1.001]],
                 [[-0.089, -0.260, -0.074], [-0.349, -1.022, 2.920], [0.458, 1.220, 0.579]],
                 [[-0.485, -1.774, -0.789], [0.305, -1.070, -1.663], [-0.707, -1.829, -0.402]],
                 [[0.029, 0.253, -0.446], [1.121, -0.173, 0.219], [-0.086, -0.562, 3.343]]],
                [[45.239, -93.595, 1085.478], [11.651, -1.436, -12.342], [1.826, -0.103, 1.736], [-0.213, 0.468, 0.272],
                 [-0.175, -1.987, 0.066], [-0.770, 0.034, 2.213]],
                [[54.360, -102.450, 1077.350], [2.223, -13.411, -5.616], [1.973, 0.274, 0.126], [0.356, -2.278, -4.988],
                 [-0.010, -0.773, 1.842], [2.928, -1.202, -0.482]],
                [[23.330, -88.250, 1275.900], [3.831, 15.852, -4.922], [1.928, -0.505, -0.127], [-0.418, -1.698, 1.100],
                 [-0.265, -0.530, -1.913], [-0.059, -0.464, 0.056]],
                [[32.287, -79.285, 1268.902],
                 [[9.823, 3.137, -8.082], [5.248, -2.802, -11.626], [7.739, 5.622, -11.497]],
                 [[1.308, -0.842, 1.262], [1.834, 0.231, 0.772], [1.699, -0.688, 0.807]],
                 [[2.146, 1.856, -0.419], [0.109, -0.816, -0.169], [-5.176, 10.836, 7.363]],
                 [[-0.217, -1.749, -0.942], [0.040, -1.938, 0.485], [-0.225, -1.720, -0.992]],
                 [[0.825, -1.334, 2.130], [-0.423, -0.029, -0.107], [11.375, -0.925, 5.860]]],
                [[39.200, -80.300, 1257.430], [9.770, 1.257, -10.029], [1.436, -0.314, 1.360], [-0.294, -0.568, 0.179],
                 [-0.103, -1.972, -0.347], [-0.113, 0.141, -0.811]],
                [[49.150, -77.930, 1247.120], [11.120, 4.120, -9.275], [1.328, -0.922, 1.182], [-0.170, -0.710, -0.365],
                 [-0.245, -1.693, -1.045], [-0.112, 0.433, -0.646]],
                [[42.582, -77.775, 1260.806], [11.777, -3.298, -2.481], [0.492, 1.925, -0.222], [2.006, -2.249, -1.967],
                 [0.441, 0.112, 1.946], [-2.413, 1.643, -1.600]],
                [[52.020, -83.330, 1261.290], [6.505, -6.503, 3.928], [1.474, 0.832, -1.064], [0.283, 1.656, -1.072],
                 [0.365, 1.271, 1.500], [2.315, 0.920, 1.295]],
                [[39.570, -106.360, 1081.580], [[2.405, -9.386, -13.810], [-5.193, -4.493, -3.038]],
                 [[1.975, 0.014, 0.334], [0.493, -1.446, 1.295]], [[-0.027, 0.596, -0.178], [6.542, -0.777, -3.568]],
                 [[-0.174, -1.658, 1.097], [-1.355, 0.694, 1.290]], [[0.913, -0.299, -0.389], [0.943, -11.695, 0.663]]],
                [[45.981, -110.526, 1067.815], [8.189, 0.886, -12.967], [1.512, 0.822, 1.011], [-0.494, 0.955, 0.159],
                 [0.753, -1.816, 0.351], [0.298, -0.005, -0.334]],
                [[39.021, -111.664, 1076.034], [1.568, -5.887, -5.929], [1.933, 0.495, 0.019], [-0.560, -0.081, -0.344],
                 [0.332, -1.353, 1.431], [-0.460, 1.100, 0.178]],
                [[38.728, -119.318, 1072.496], [-5.071, -7.838, 1.726], [1.535, -0.769, 1.017], [-1.382, -0.287, 2.451],
                 [-0.700, 0.823, 1.679], [-1.088, 5.059, -1.024]],
                [[39.676, -84.700, 1457.627], [-0.784, -3.974, 0.816], [1.955, -0.330, 0.271], [-2.648, 2.244, 0.223],
                 [-0.195, 0.436, 1.938], [0.117, 0.333, 0.132]],
                [[38.747, -87.062, 1457.995], [[-1.667, -1.652, -0.090], [7.134, -4.236, -3.417]],
                 [[1.138, -1.132, -0.286], [1.004, 1.730, -0.049]], [[-5.394, -3.462, -3.050], [-1.676, 0.627, -0.026]],
                 [[0.193, -0.302, 1.964], [0.545, -0.275, 1.479]], [[-1.008, 0.289, 0.085], [0.020, 0.430, -0.864]]],
                [[35.160, -89.575, 1457.802], [0.012, -1.294, 1.216], [1.869, 0.499, 0.513], [3.047, 3.780, -4.307],
                 [-0.714, 1.273, 1.362], [-3.092, 0.209, -0.952]],
                [[34.749, -91.950, 1459.539], [[-0.232, -1.383, -1.290], [-0.610, -1.701, -0.102]],
                 [[1.467, -1.052, 0.863], [1.291, -0.461, -0.032]], [[-1.435, -2.186, 0.640], [0.655, 3.019, 1.826]],
                 [[-1.336, -0.886, 1.190], [0.006, -0.122, 1.992]], [[0.100, -2.332, 0.555], [1.177, 2.923, -0.668]]],
                [[43.337, -88.990, 1455.261], [1.677, -0.783, -1.601], [0.825, 1.824, -0.028], [1.180, -0.786, -0.171],
                 [0.953, -0.413, 1.201], [-0.157, -0.176, 0.101]],
                [[45.707, -90.450, 1452.310], [-0.261, -1.599, -4.023], [1.977, -0.309, -0.005],
                 [-0.402, -2.390, 1.432], [-0.210, -1.353, 0.551], [-2.467, 0.644, 0.679]],
                [[34.325, -94.875, 1457.084], [[-0.408, -1.566, -1.654], [-0.693, 1.011, -2.864]],
                 [[1.663, -0.984, 0.522], [1.463, 1.360, 0.126]], [[0.382, 0.049, -0.551], [-3.299, -0.725, -1.584]],
                 [[-1.053, -1.094, 1.295], [1.287, -1.313, -0.775]], [[0.865, 0.469, 0.331], [-5.491, -1.302, -0.513]]],
                [[34.504, -98.267, 1459.802], [0.240, -2.029, 0.244], [1.899, 0.292, 0.562], [0.340, -2.130, -1.581],
                 [-0.654, 0.177, 2.119], [-0.219, -2.534, -0.279]],
                [[34.900, -104.258, 1459.373], [1.031, -2.469, -3.084], [1.926, 0.473, 0.265], [-1.919, 9.009, 3.535],
                 [0.193, -1.490, 1.257], [0.732, 6.172, 1.238]],
                [[33.310, -98.196, 1453.203], [-0.481, -2.748, -4.844], [1.957, 0.248, -0.335], [0.521, 1.579, 0.075],
                 [0.378, -1.720, 0.938], [-1.282, -0.623, -1.583]],
                [[32.696, -101.165, 1448.080], [-0.537, -1.773, -3.276], [0.515, -1.801, 0.890],
                 [-0.519, -2.257, 0.799], [-1.912, -0.309, 0.481], [-0.082, 0.981, 1.355]],
                [[33.160, -93.310, 1451.928], [0.666, 0.600, -2.406], [1.932, -0.086, 0.513], [0.833, -1.668, 3.491],
                 [0.039, -1.938, -0.472], [-0.117, 0.684, 1.633]],
                [[35.042, -92.561, 1445.463], [-1.162, -0.857, -5.000], [1.951, -0.091, -0.438],
                 [-0.831, 2.215, -2.002], [-0.015, -1.966, 0.341], [-0.067, -0.775, -1.309]]
            ]

            bifurcating_nodes = [8, 19, 34, 38, 43, 48, 50, 53]

            nodesByGroup = [
                ['left vagus X nerve trunk', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 43, 44]],
                ['left A lateral pulmonary branch of vagus nerve', [43, 45, 46]],
                ['left A medial pulmonary branch of vagus nerve', [34, 35, 36]],
                ['left A thoracic cardiac branch of vagus nerve', [19, 37, 38, 41, 42]],
                ['left B thoracic cardiac branch of vagus nerve', [38, 39, 40]],
                # additional branches
                ['left A superior laryngeal nerve', [8, 47, 48]],
                ['left B internal branch of superior laryngeal nerve', [48, 49, 50, 53]],
                ['left B external branch of superior laryngeal nerve', [48, 51, 52]],
                # not in the interlex termset
                ['left superior/upper ramus of Internal branch of L Superior laryngeal nerve', [50, 54, 55]],
                ['left middle ramus of Internal branch of L Superior laryngeal nerve', [53, 56, 57]],
                ['left inferior ramus of Internal branch of L Superior laryngeal nerve', [53, 58, 59]]
            ]

        for n in range(len(outerPoints)):
            node = nodes.findNodeByIdentifier(n + 1)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, outerPoints[n][0])
            if (n + 1) in bifurcating_nodes:
                for v in range(len(outerPoints[n][1])):
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, v + 1, outerPoints[n][1][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, v + 1, outerPoints[n][2][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, v + 1, outerPoints[n][3][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, v + 1, outerPoints[n][4][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, v + 1, outerPoints[n][5][v])
            else:
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, outerPoints[n][1])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, outerPoints[n][2])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, outerPoints[n][3])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, outerPoints[n][4])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, outerPoints[n][5])

        annotationGroups = []
        for group_n in range(len(nodesByGroup)):
            branch_name = nodesByGroup[group_n][0]

            branch_field_group = findOrCreateFieldGroup(fieldmodule, branch_name)
            branch_nodeset = branch_field_group.createNodesetGroup(nodes)

            for n in range(len(nodesByGroup[group_n][1])):
                node = nodes.findNodeByIdentifier(nodesByGroup[group_n][1][n])
                branch_nodeset.addNode(node)

        # create markers
        markers = [
            ["level of aortic hiatus on the vagus nerve", [33.06685565736445, -99.79675650489413, 1096.533330802349]],
            ["level of esophageal hiatus on the vagus nerve", [3.146209400855563, -100.0347608405455, 1139.307993871152]],
            ["level of sternal angle on the vagus nerve", [22.05489076832703, -104.8332280606236, 1283.890550463673]],
            ["level of jugular notch on the vagus nerve", [29.38481050304174, -110.2391429293425, 1323.108962286499]],
            ["level of superior border of the clavicle on the vagus nerve", [28.41567111857362, -107.0795092601865, 1333.898096789418]],
            ["level of laryngeal prominence on the vagus nerve", [27.86043189405682, -86.62430668141346, 1407.793206838608]],
            ["level of greater horn of hyoid on the vagus nerve", [35.62824151311525, -83.82816000000001, 1443.993933756558]],
            ["level of carotid bifurcation on the vagus nerve", [35.11637594753900, -84.55354200000002, 1439.799758114304]],
            ["level of angle of the mandible on the vagus nerve", [36.11596477766252, -82.68291542919793, 1449.201631366949]],
            ["level of C1 transverse process on the vagus nerve", [35.61590793674574, -78.40798854214731, 1465.004007224821]],
            ["level of inferior border of jugular foramen on the vagus nerve", [29.00379183524906, -73.43913993329343, 1489.753246601004]],
            ["level of superior border of jugular foramen on the vagus nerve", [22.73942104901961, -68.45420763548567, 1496.786534505575]]
        ]

        # set markers
        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        data_coordinates = findOrCreateFieldCoordinates(fieldmodule, "marker_data_coordinates")
        marker_names = findOrCreateFieldStoredString(fieldmodule, name="marker_data_name")

        dnodetemplate = datapoints.createNodetemplate()
        dnodetemplate.defineField(data_coordinates)
        dnodetemplate.setValueNumberOfVersions(data_coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        dnodetemplate.defineField(marker_names)

        marker_fieldgroup = findOrCreateFieldGroup(fieldmodule, 'marker')
        marker_nodesetgroup = marker_fieldgroup.createNodesetGroup(datapoints)
        for n in range(len(markers)):
            node = datapoints.createNode(n + 1, dnodetemplate)
            fieldcache.setNode(node)
            data_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, markers[n][1])
            marker_names.assignString(fieldcache, markers[n][0])
            marker_nodesetgroup.addNode(node)

        if defineInnerCoordinates:
            cls._defineInnerCoordinates(region, coordinates, options, networkMesh)

        return [], networkMesh

    @classmethod
    def _defineInnerCoordinates(cls, region, coordinates, options, networkMesh):
        """
        Copy coordinates to inner coordinates via in-memory model file.
        Assign using the interactive function.
        :param region: Region to define field in.
        :param coordinates: Standard/outer coordinate field.
        :param options: Options used to generate scaffold.
        :param networkMesh: Network mesh object used to generate scaffold.
        """
        assert options["Define inner coordinates"]
        coordinates.setName("inner coordinates")  # temporarily rename
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        region.write(sir)
        result, buffer = srm.getBuffer()
        coordinates.setName("coordinates")  # restore name before reading inner coordinates back in
        sir = region.createStreaminformationRegion()
        sir.createStreamresourceMemoryBuffer(buffer)
        region.read(sir)
        functionOptions = {
            "To field": {"coordinates": False, "inner coordinates": True},
            "From field": {"coordinates": True, "inner coordinates": False},
            "Mode": {"Scale": True, "Offset": False},
            "D2 value": 0.5,
            "D3 value": 0.5}
        cls.assignCoordinates(region, options, networkMesh, functionOptions, editGroupName=None)

    @classmethod
    def editStructure(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Edit structure safely, to prevent accidental changes.
        Copies functionOptions["Structure"] to options["Structure"] and regenerates with
        default geometric coordinates.
        :param region: Region containing model to clear and re-generate.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Contents replaced.
        :param functionOptions: functionOptions["Structure"] contains new structure string.
        :param editGroupName: Name of Zinc group to put edited nodes in. Cleared.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            clearRegion(region)
            structure = options["Structure"] = functionOptions["Structure"]
            networkMesh.build(structure)
            networkMesh.create1DLayoutMesh(region)
            coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
            coordinates.setManaged(True)  # since cleared by clearRegion
            defineInnerCoordinates = options["Define inner coordinates"]
            if defineInnerCoordinates:
                cls._defineInnerCoordinates(region, coordinates, options, networkMesh)

        return True, False  # settings changed, nodes not changed (since reset to original coordinates)

    class AssignCoordinatesMode(Enum):
        SCALE = 1,  # scale side derivative magnitude by value
        OFFSET = 2   # offset side derivative by absolute distance

    @classmethod
    def assignCoordinates(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Assign coordinates or inner coordinates by scaling or offsetting from either original
        values or values from other field.
        If elements are selected, applied only to nodes used by the element and versions used in the enclosing segment.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Unused.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        innerCoordinates = fieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if (functionOptions["To field"]["inner coordinates"] or functionOptions["From field"]["inner coordinates"]) \
                and not innerCoordinates.isValid():
            print("Assign coordinates:  inner coordinates field not defined")
            return False, False
        mode = None
        if functionOptions["Mode"]["Scale"]:
            mode = cls.AssignCoordinatesMode.SCALE
        elif functionOptions["Mode"]["Offset"]:
            mode = cls.AssignCoordinatesMode.OFFSET
        else:
            print("Assign coordinates:  Invalid mode")
            return False, False
        toCoordinates = coordinates if functionOptions["To field"]["coordinates"] else innerCoordinates
        fromCoordinates = coordinates if functionOptions["From field"]["coordinates"] else innerCoordinates
        d2Value = functionOptions["D2 value"]
        d3Value = functionOptions["D3 value"]
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        selectionGroup = scene_get_selection_group(region.getScene(), inherit_root_region=region.getRoot())
        selectionMeshGroup = None
        mesh1d = fieldmodule.findMeshByDimension(1)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        if selectionGroup:
            selectionMeshGroup = selectionGroup.getMeshGroup(mesh1d)
            if not selectionMeshGroup.isValid():
                print("Assign coordinates:  Selection contains no elements. Clear it to assign globally.")
                return False, False
        valueLabels = [
            Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        with ChangeManager(fieldmodule):
            # get all node parameters (from selection if any)
            editNodeset = nodeset
            originalNodeParameters = None
            if selectionGroup:
                # make group of only nodes being edited
                tmpGroup = fieldmodule.createFieldGroup()
                tmpGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
                tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
                tmpMeshGroup.addElementsConditional(selectionGroup)
                editNodeset = tmpGroup.getNodesetGroup(nodes)
                _, originalNodeParameters = get_nodeset_field_parameters(editNodeset, toCoordinates, valueLabels)
                del tmpMeshGroup
                del tmpGroup
            _, nodeParameters = get_nodeset_field_parameters(editNodeset, fromCoordinates, valueLabels)

            modifyVersions = None  # default is to modify all versions
            if selectionGroup:
                nodeIdentifierIndexes = {}
                modifyVersions = []
                for n in range(len(nodeParameters)):
                    nodeIdentifierIndexes[nodeParameters[n][0]] = n
                    versionsCount = len(nodeParameters[n][1][1])
                    modifyVersions.append([False] * versionsCount if selectionGroup else [True] * versionsCount)
                networkSegments = networkMesh.getNetworkSegments()
                for networkSegment in networkSegments:
                    nodeIdentifiers = networkSegment.getNodeIdentifiers()
                    nodeVersions = networkSegment.getNodeVersions()
                    elementIdentifiers = networkSegment.getElementIdentifiers()
                    for e in range(len(elementIdentifiers)):
                        elementIdentifier = elementIdentifiers[e]
                        element = selectionMeshGroup.findElementByIdentifier(elementIdentifier)
                        if element.isValid():
                            for n in [e, e + 1]:
                                nodeIndex = nodeIdentifierIndexes.get(nodeIdentifiers[n])
                                # print("Node identifier", nodeIdentifiers[n], "index", nodeIndex, "version", nodeVersions[n])
                                if nodeIndex is not None:
                                    modifyVersions[nodeIndex][nodeVersions[n] - 1] = True

            for n in range(len(nodeParameters)):
                modifyVersion = modifyVersions[n] if modifyVersions else None
                nNodeParameters = nodeParameters[n][1]
                oNodeParameters = originalNodeParameters[n][1] if modifyVersions else None
                versionsCount = len(nNodeParameters[1])
                for v in range(versionsCount):
                    if (not modifyVersions) or modifyVersion[v]:
                        for dd in range(2):
                            scale = d2Value if (dd == 0) else d3Value
                            if mode == cls.AssignCoordinatesMode.OFFSET:
                                mag = magnitude(nNodeParameters[2 + 2 * dd][v])
                                scale = (mag + scale) / mag if (abs(mag) > 0.0) else 1.0
                            for d in [2 + 2 * dd, 3 + 2 * dd]:
                                nNodeParameters[d][v] = mult(nNodeParameters[d][v], scale)
                    else:
                        # copy original derivative versions
                        for d in range(2, 6):
                            nNodeParameters[d][v] = oNodeParameters[d][v]

            set_nodeset_field_parameters(editNodeset, toCoordinates, valueLabels, nodeParameters, editGroupName)
            del editNodeset

        return False, True  # settings not changed, nodes changed

    @classmethod
    def makeSideDerivativesNormal(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Make side directions normal to d1 and each other. Works for all versions.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Unused.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        innerCoordinates = fieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if functionOptions["Field"]["inner coordinates"] and not innerCoordinates.isValid():
            print("Make side derivatives normal:  inner coordinates field not defined")
            return False, False
        useCoordinates = coordinates if functionOptions["Field"]["coordinates"] else innerCoordinates
        makeD2Normal = functionOptions['Make D2 normal']
        makeD3Normal = functionOptions['Make D3 normal']
        if not (makeD2Normal or makeD3Normal):
            return False, False
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        make_nodeset_derivatives_orthogonal(nodeset, useCoordinates, makeD2Normal, makeD3Normal, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def smoothSideCrossDerivatives(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Smooth side cross derivatives giving rate of change of side directions d2, d3 w.r.t. d1.
        If a single element in a segment is selected, the whole segment is smoothed, and if the segment
        connects to others with the same version, they are also smoothed with it.
        Also detects loops back to the start of a segment.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from.
        Used to determine connected paths for smoothing.
        :param functionOptions: Which side derivatives to smooth.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        innerCoordinates = fieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if functionOptions["Field"]["inner coordinates"] and not innerCoordinates.isValid():
            print("Make side derivatives normal:  inner coordinates field not defined")
            return False, False
        useCoordinates = coordinates if functionOptions["Field"]["coordinates"] else innerCoordinates
        smoothD12 = functionOptions["Smooth D12"]
        smoothD13 = functionOptions["Smooth D13"]
        if not (smoothD12 or smoothD13):
            return False, False

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        selectionGroup = scene_get_selection_group(region.getScene(), inherit_root_region=region.getRoot())
        selectionMeshGroup = None
        mesh1d = fieldmodule.findMeshByDimension(1)
        if selectionGroup:
            selectionMeshGroup = selectionGroup.getMeshGroup(mesh1d)
            if not selectionMeshGroup.isValid():
                print("Smooth side cross derivatives:  Selection must contain elements to smooth segment chains "
                    "containing them, or be clear it to smooth all segment chains.")
                return False, False
        getValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
        setValueLabels = []
        if smoothD12:
            getValueLabels.append(Node.VALUE_LABEL_D_DS2)
            setValueLabels.append(Node.VALUE_LABEL_D2_DS1DS2)
        if smoothD13:
            getValueLabels.append(Node.VALUE_LABEL_D_DS3)
            setValueLabels.append(Node.VALUE_LABEL_D2_DS1DS3)

        # determine segment chains which must be smoothed together:
        # currently only links segments aligned in the same direction
        networkSegments = networkMesh.getNetworkSegments()
        segmentChains = []
        segmentChainsLoop = []  # True if same index segment chain is a loop
        processedSegments = set()
        for segment in networkSegments:
            if segment in processedSegments:
                continue
            segmentChain = [segment]
            processedSegments.add(segment)
            # add other non-processed segments attached and using the same derivative version
            startNetworkNode = segment.getNetworkNodes()[0]
            startNodeVersion = segment.getNodeVersions()[0]
            endNetworkNode = segment.getNetworkNodes()[-1]
            endNodeVersion = segment.getNodeVersions()[-1]
            while True:
                for outSegment in endNetworkNode.getOutSegments():
                    if ((outSegment.getNodeVersions()[0] == endNodeVersion) and
                            outSegment not in processedSegments):
                        segmentChain.append(outSegment)
                        processedSegments.add(outSegment)
                        endNetworkNode = outSegment.getNetworkNodes()[-1]
                        endNodeVersion = outSegment.getNodeVersions()[-1]
                        break
                else:
                    break
            while True:
                for inSegment in startNetworkNode.getInSegments():
                    if ((inSegment.getNodeVersions()[-1] == startNodeVersion) and
                            inSegment not in processedSegments):
                        segmentChain.insert(0, inSegment)
                        processedSegments.add(inSegment)
                        startNetworkNode = inSegment.getNetworkNodes()[0]
                        startNodeVersion = inSegment.getNodeVersions()[0]
                        break
                else:
                    break
            segmentChains.append(segmentChain)
            segmentChainsLoop.append((startNetworkNode == endNetworkNode) and (startNodeVersion == endNodeVersion))

        with ChangeManager(fieldmodule):

            editNodeset = nodes
            if selectionGroup:
                # include only chains containing a selected element
                chainIndex = 0
                while chainIndex < len(segmentChains):
                    segmentChain = segmentChains[chainIndex]
                    for segment in segmentChain:
                        if segment.hasLayoutElementsInMeshGroup(selectionMeshGroup):
                            break
                    else:
                        segmentChains.pop(chainIndex)
                        segmentChainsLoop.pop(chainIndex)
                        continue
                    chainIndex += 1
                # make group of only nodes being edited
                tmpGroup = fieldmodule.createFieldGroup()
                editNodeset = tmpGroup.createNodesetGroup(nodes)
                for segmentChain in segmentChains:
                    for segment in segmentChain:
                        for nodeIdentifier in segment.getNodeIdentifiers():
                            editNodeset.addNode(nodes.findNodeByIdentifier(nodeIdentifier))
                del tmpGroup

            _, nodeParameters = get_nodeset_field_parameters(editNodeset, useCoordinates, setValueLabels)
            nodeIdentifierIndexes = {}
            for n in range(len(nodeParameters)):
                nodeIdentifierIndexes[nodeParameters[n][0]] = n

            for chainIndex in range(len(segmentChains)):
                # get parameters for chain
                segmentChain = segmentChains[chainIndex]
                loop = segmentChainsLoop[chainIndex]
                segmentsCount = len(segmentChain)
                nx = []
                nd1 = []
                sideVectorsCount = 2 if smoothD12 and smoothD13 else 1
                nsv = [[] for s in range(sideVectorsCount)]
                nodeIdentifiers = []
                nodeVersions = []
                for segmentIndex in range(segmentsCount):
                    segment = segmentChain[segmentIndex]
                    segmentNodeIdentifiers = segment.getNodeIdentifiers()
                    segmentNodeVersions = segment.getNodeVersions()
                    segmentParameters = get_nodeset_path_ordered_field_parameters(
                        nodes, useCoordinates, getValueLabels, segmentNodeIdentifiers, segmentNodeVersions)
                    nodesCount = len(segmentNodeIdentifiers)
                    if loop or (segmentIndex < (segmentsCount - 1)):
                        nodesCount -= 1
                    nx += segmentParameters[0][:nodesCount]
                    nd1 += segmentParameters[1][:nodesCount]
                    for s in range(sideVectorsCount):
                        nsv[s] += segmentParameters[2 + s][:nodesCount]
                    nodeIdentifiers += segmentNodeIdentifiers[:nodesCount]
                    nodeVersions += segmentNodeVersions[:nodesCount]
                dnsv = smoothCurveSideCrossDerivatives(nx, nd1, nsv, loop=loop)
                for n in range(len(nodeIdentifiers)):
                    nodeIndex = nodeIdentifierIndexes.get(nodeIdentifiers[n])
                    nodeVersion = nodeVersions[n] - 1
                    assert nodeIndex is not None
                    for s in range(sideVectorsCount):
                        nodeParameters[nodeIndex][1][s][nodeVersion] = dnsv[s][n]

            set_nodeset_field_parameters(editNodeset, useCoordinates, setValueLabels, nodeParameters, editGroupName)
            del editNodeset

        return False, True  # settings not changed, nodes changed

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Supply client with functions for smoothing path parameters.
        """
        # add choice of field to base functions
        modifiedBaseInteractiveFunctions = []
        for interactiveFunction in Scaffold_base.getInteractiveFunctions():
            dct = {"Field": {"coordinates": True, "inner coordinates": False}}
            dct.update(interactiveFunction[1])
            modifiedBaseInteractiveFunctions.append((
               interactiveFunction[0], dct,
               interactiveFunction[2]))
        return modifiedBaseInteractiveFunctions + [
            ("Edit structure...", {
                "Structure": None},  # None = take value from options
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.editStructure(region, options, networkMesh, functionOptions, editGroupName)),
            ("Assign coordinates...", {
                "To field": {"coordinates": True, "inner coordinates": False},
                "From field": {"coordinates": True, "inner coordinates": False},
                "Mode": {"Scale": True, "Offset": False},
                "D2 value": 1.0,
                "D3 value": 1.0},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.assignCoordinates(region, options, networkMesh, functionOptions, editGroupName)),
            ("Make side derivatives normal...", {
                "Field": {"coordinates": True, "inner coordinates": False},
                "Make D2 normal": True,
                "Make D3 normal": True},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.makeSideDerivativesNormal(region, options, networkMesh, functionOptions, editGroupName)),
            ("Smooth side cross derivatives...", {
                "Field": {"coordinates": True, "inner coordinates": False},
                "Smooth D12": True,
                "Smooth D13": True},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.smoothSideCrossDerivatives(region, options, networkMesh, functionOptions, editGroupName))
        ]
