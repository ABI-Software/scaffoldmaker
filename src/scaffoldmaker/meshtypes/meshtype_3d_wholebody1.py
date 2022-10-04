"""
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
"""

from __future__ import division

import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from opencmiss.zinc.node import Node
from opencmiss.zinc.element import Element
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.annotation.nerve_terms import get_nerve_term
from scaffoldmaker.annotation import heart_terms, bladder_terms, lung_terms, stomach_terms, brainstem_terms
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field, FieldFindMeshLocation
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.vector import setMagnitude
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues

from opencmiss.zinc.context import Context


class MeshType_3d_wholebody1(Scaffold_base):
    """
Generates body coordinates using a solid cylinder of all cube elements,
 with variable numbers of elements in major, minor, shell and for each section of abdomen, thorax, neck and head.
    """
    cylinder1Settings = {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.5,
                'Number of elements': 1
            }

    axis1 = [0, 0, 1]
    axis2 = [1, 0, 0]
    axis3 = [0, 1, 0]
    centralPathDefaultScaffoldPackages = {
        'Default': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': cylinder1Settings,
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, 0.0, 0.0], setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]],
                    [setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis1, cylinder1Settings['Length']), setMagnitude(axis2, 0.5), [0.0, 0.0, 0.0], setMagnitude(axis3, 0.5), [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Whole Body 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human Coarse',
            'Human Fine',
            'Rat Coarse',
            'Rat Fine'
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Human Coarse'
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Default']
        options = {}
        options['Base parameter set'] = parameterSetName
        options['Central path'] = copy.deepcopy(centralPathOption)
        options['Number of elements across major'] = 6
        options['Number of elements across minor'] = 6
        options['Number of elements across shell'] = 1
        options['Number of elements across transition'] = 1
        options['Number of elements in abdomen'] = 5
        options['Number of elements in thorax'] = 3
        options['Number of elements in neck'] = 1
        options['Number of elements in head'] = 2
        options['Shell thickness proportion'] = 0.2
        options['Discontinuity on the core boundary'] = True
        options['Lower half'] = False
        options['Use cross derivatives'] = False
        options['Refine'] = False
        options['Refine number of elements across major'] = 1
        options['Refine number of elements along'] = 1

        if 'Coarse' in parameterSetName:
            pass
        if 'Fine' in parameterSetName:
            options['Number of elements across major'] = 10
            options['Number of elements across minor'] = 10
            options['Number of elements across shell'] = 1
            options['Number of elements across transition'] = 1
            options['Number of elements in abdomen'] = 10
            options['Number of elements in thorax'] = 6
            options['Number of elements in neck'] = 2
            options['Number of elements in head'] = 4
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements across transition',
            'Number of elements in abdomen',
            'Number of elements in thorax',
            'Number of elements in neck',
            'Number of elements in head',
            'Shell thickness proportion',
            'Discontinuity on the core boundary',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] != options['Number of elements across major']:
            options['Number of elements across minor'] = options['Number of elements across major']
            dependentChanges = True

        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1

        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1

        if options['Shell thickness proportion'] < 0.07 or options['Shell thickness proportion'] > 0.7:
            options['Shell thickness proportion'] = 2*options['Number of elements across shell']/options['Number of elements across major']

        if options['Number of elements in abdomen'] < 1:
            options['Number of elements in abdomen'] = 1
        if options['Number of elements in head'] < 1:
            options['Number of elements in head'] = 1
        if options['Number of elements in neck'] < 1:
            options['Number of elements in neck'] = 1
        if options['Number of elements in thorax'] < 1:
            options['Number of elements in thorax'] = 1

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: List of AnnotationGroup
        """

        baseParameterSetName = options['Base parameter set']
        isHuman = 'Human' in baseParameterSetName
        isRat = 'Rat' in baseParameterSetName

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        elementsCountAlongAbdomen = options['Number of elements in abdomen']
        elementsCountAlongHead = options['Number of elements in head']
        elementsCountAlongNeck = options['Number of elements in neck']
        elementsCountAlongThorax = options['Number of elements in thorax']
        shellRadiusProportion = options['Shell thickness proportion']
        shellProportion = 1/(1/shellRadiusProportion-1)*(elementsCountAcrossMajor/2/elementsCountAcrossShell - 1)
        discontinuity = options['Discontinuity on the core boundary']
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAlong = elementsCountAlongAbdomen + elementsCountAlongThorax + elementsCountAlongNeck + elementsCountAlongHead

        # region.readFile(
        #     fileName=r'C:\Users\egha355\Desktop\sparc3\codes\mapclient_workflows\workflowandfilesnewer2\whole-body.exf')
        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        mesh = fieldmodule.findMeshByDimension(3)

        # read whole body and get the nerves coordinates.
        markerGroup = findOrCreateFieldGroup(fieldmodule, "nerve_marker")
        markerName = findOrCreateFieldStoredString(fieldmodule, name="nerve_marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fieldmodule, mesh, name="nerve_marker_location")
        markerBodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="nerve_marker_body_coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        markerTemplateInternal.defineField(markerBodyCoordinates)
        fieldcache = fieldmodule.createFieldcache()

        # context = Context('nerves_centrelines')
        # region_centreline = context.createRegion()
        # region_centreline.readFile(
        #     filename=r'C:\Users\egha355\Desktop\sparc3\codes\mapclient_workflows\workflowandfilesnewer2\nerves_centrelines_anatomical_terms.exf')
        bodyMarkerPoints = [
            {"name": "C3_C5_spinal_n-phrenic_n", "x": [-1.452000000000079, -79.02119241690319, 1419.237697325146]},
            {"name": "C3_dorsal_root_end", "x": [1.396362832762293, -55.00539191135986, 1443.926916470878]},
            {"name": "C3_spinal", "x": [-0.626774973288383, -63.22616471288262, 1440.230044789433]},
            {"name": "C3_ventral_root_end", "x": [0.4858310577370967, -72.5919787158736, 1443.788470489924]},
            {"name": "C4_dorsal_root_end", "x": [1.164187341105573, -57.19028747488394, 1429.93350827451]},
            {"name": "C4_spinal", "x": [-0.5069999999999472, -65.71002921159031, 1426.59836086468]},
            {"name": "C4_ventral_root_end", "x": [0.5667239873419503, -74.44667178877705, 1430.294861316902]},
            {"name": "C5_dorsal_root_end", "x": [2.173598614138331, -60.42114565823348, 1415.807665476658]},
            {"name": "C5_spinal", "x": [0.02800000000000881, -68.24970100036765, 1410.454521889242]},
            {"name": "C5_ventral_root_end", "x": [1.676395914893624, -75.10201344553963, 1415.732585272465]},
            {"name": "C6_dorsal_root_end", "x": [4.882350184702168, -59.1440776310951, 1401.346897212325]},
            {"name": "C6_spinal", "x": [-1.30117143690536, -67.5507086835333, 1395.86269025303]},
            {"name": "C7_dorsal_root_end", "x": [4.008350343127929, -57.78605272356008, 1384.104049712696]},
            {"name": "C7_spinal", "x": [-0.7650000000000441, -64.58619464257978, 1379.770361458746]},
            {"name": "C8_dorsal_root_end", "x": [1.798709380245413, -53.97641867054194, 1368.958443066533]},
            {"name": "C8_spinal", "x": [-0.9549857257449329, -60.39046326827202, 1365.549452962976]},
            {"name": "L1_dorsal_root_end", "x": [4.685701126153761, -51.95555777868519, 1040.965084959929]},
            {"name": "L1_spinal", "x": [-1.385258279328277, -68.28496698781258, 1034.176082973243]},
            {"name": "L1_ventral_root_end", "x": [4.575050287268543, -81.85704170115409, 1040.630206539322]},
            {"name": "L2_dorsal_root_end", "x": [4.633723554346656, -54.56666724789302, 1006.823502398816]},
            {"name": "L2_spinal", "x": [-0.75, -71.05631277336954, 1000.505816674577]},
            {"name": "L2_ventral_root_end", "x": [4.839532071705589, -88.85690616504627, 1006.462797567221]},
            {"name": "L5_spinal", "x": [-0.9763829652823178, -43.57788014910619, 908.6694912598865]},
            {"name": "L5_ventral_root_end", "x": [5.102627813026095, -63.38712510346181, 911.3513731278158]},
            {"name": "S1_dorsal_root_end", "x": [5.956829968687891, -21.08340052016078, 890.7335602659736]},
            {"name": "S1_spinal", "x": [-0.1050000000000881, -29.0606498699599, 886.5936709471108]},
            {"name": "S1_ventral_root_end", "x": [5.618702313593648, -47.59252996810586, 889.0519291628217]},
            {"name": "S43_C3", "x": [1.425063060609067, -71.71810165621304, 1442.050240627686]},
            {"name": "S43_C4", "x": [0.2390704999228729, -73.36366794719446, 1428.269653359314]},
            {"name": "S43_C5", "x": [1.852252439519412, -73.16334165676288, 1413.24187389275]},
            {"name": "S43_L1", "x": [4.411206711081232, -84.77015985932856, 1039.334066158587]},
            {"name": "S43_L2", "x": [4.884659726799778, -84.87726932853623, 1002.793911898404]},
            {"name": "S43_L5", "x": [4.726074605855652, -58.48225629508627, 909.7502183074077]},
            {"name": "S43_S1", "x": [5.249478917179068, -41.97307471334122, 886.0499300842902]},
            {"name": "S43_T1", "x": [1.777893687727883, -60.18363072973604, 1348.206592824946]},
            {"name": "S43_T10", "x": [3.666875093710453, -53.959472797532, 1126.244854124646]},
            {"name": "S43_T11", "x": [4.349602376174059, -60.03004392000878, 1101.343332732087]},
            {"name": "S43_T12", "x": [3.510778497324851, -73.81028620243222, 1071.821164183794]},
            {"name": "S43_T2", "x": [3.356488708167721, -52.76245739764429, 1326.923395206977]},
            {"name": "S43_T3", "x": [3.075556887551878, -44.42742009271725, 1307.951920178592]},
            {"name": "S43_T4", "x": [2.433051214861331, -39.10305279194664, 1284.555130418246]},
            {"name": "S43_T5", "x": [1.906625066936073, -36.31555641658424, 1257.394833466355]},
            {"name": "S43_T6", "x": [2.548062560242515, -38.01459769742712, 1234.110017705623]},
            {"name": "S43_T7", "x": [3.672863748329545, -38.31072918942786, 1208.495424568588]},
            {"name": "S43_T8", "x": [3.333448945785247, -41.85006644003461, 1182.342939805176]},
            {"name": "S43_T9", "x": [3.384215990631239, -46.4028293411147, 1157.51572660028]},
            {"name": "S44_C3", "x": [1.746230242904172, -55.27287813589524, 1442.892470303966]},
            {"name": "S44_C4", "x": [1.163205514396105, -58.46499039604345, 1428.916549501825]},
            {"name": "S44_C5", "x": [1.917007776697941, -61.1926621848203, 1414.763083156512]},
            {"name": "S44_C6", "x": [5.014580398486289, -60.31452036816938, 1400.498144330297]},
            {"name": "S44_C7", "x": [4.426339163537312, -58.7061049274637, 1383.127783957719]},
            {"name": "S44_C8", "x": [1.882307476737945, -54.755597581118, 1368.021061739523]},
            {"name": "S44_L1", "x": [4.757472075362164, -53.31450304408124, 1039.122687882434]},
            {"name": "S44_L2", "x": [4.905212933970203, -56.17478975055089, 1005.079274253941]},
            {"name": "S44_S1", "x": [5.183053365102638, -22.11270026008045, 889.0540880372521]},
            {"name": "S44_T1", "x": [2.456693922907685, -48.39494831575267, 1350.606027969544]},
            {"name": "S44_T2", "x": [3.069789842706019, -39.09663555166344, 1329.295928107802]},
            {"name": "S44_T3", "x": [3.006613721554929, -30.49803480934947, 1309.175322738125]},
            {"name": "S44_T4", "x": [2.159488708167728, -23.65364626221666, 1286.078194388809]},
            {"name": "S44_T5", "x": [1.916090977243976, -18.23450003671882, 1258.786294498959]},
            {"name": "S44_T6", "x": [2.508323932398037, -14.13736253031412, 1235.106919769094]},
            {"name": "S45_C3", "x": [1.623884426753343, -55.91189744234539, 1442.018079694941]},
            {"name": "S45_C4", "x": [1.398722478337561, -59.63947194866488, 1428.064742071309]},
            {"name": "S45_C5", "x": [2.14011999452222, -62.34055150055152, 1414.089193019728]},
            {"name": "S45_C6", "x": [5.018164348819896, -61.35778055225428, 1399.546248129108]},
            {"name": "S45_C7", "x": [4.496225716894767, -59.74478409786778, 1382.493087424291]},
            {"name": "S45_C8", "x": [2.060261622545009, -55.79038115932958, 1367.503726639785]},
            {"name": "S45_L1", "x": [4.365435761872759, -57.37208949047061, 1037.835890799218]},
            {"name": "S45_L2", "x": [4.682563985479157, -58.02526238324887, 1003.269183665057]},
            {"name": "S45_S1", "x": [5.892149106423208, -22.84233193461059, 888.863283206919]},
            {"name": "S45_T1", "x": [1.961866904579874, -49.52594831575269, 1348.91312773834]},
            {"name": "S45_T2", "x": [3.220721663321772, -40.79610728653967, 1328.414292611359]},
            {"name": "S45_T3", "x": [4.202000000000079, -31.91082177681067, 1307.968903645272]},
            {"name": "S45_T4", "x": [2.160801188086915, -25.23255233699244, 1285.116664048388]},
            {"name": "S45_T5", "x": [1.687363694780525, -19.84345130301514, 1257.619038770136]},
            {"name": "S45_T6", "x": [2.547284170015381, -16.39566995492389, 1233.157007265408]},
            {"name": "S46_L1", "x": [4.679787381893746, -70.9586033673817, 1031.133998458225]},
            {"name": "S46_L2", "x": [5.11205336510258, -73.11445182699745, 997.3921513331876]},
            {"name": "S46_L3", "x": [5.723851209441035, -76.51866317166098, 971.1794385712689]},
            {"name": "S46_L4", "x": [5.30850015793199, -63.24391816699476, 932.3618358187399]},
            {"name": "S46_L5", "x": [4.432266088496743, -46.70766198300391, 900.2523734445531]},
            {"name": "S46_S1", "x": [5.506989537555432, -31.1793028908065, 882.2705681266476]},
            {"name": "S46_S2", "x": [5.431468244158385, -22.71088030029158, 868.3838711854864]},
            {"name": "S46_T10", "x": [3.299789842705889, -45.13452086587362, 1120.093112890011]},
            {"name": "S46_T11", "x": [3.610806914326476, -53.71774159845062, 1093.695147780607]},
            {"name": "S46_T12", "x": [3.294255753013744, -64.79126449747476, 1064.073329896692]},
            {"name": "S46_T2", "x": [-1.880000000000088, -45.16650124140453, 1324.157659309933]},
            {"name": "S46_T3", "x": [2.143113668006121, -38.42710189406393, 1304.497266657217]},
            {"name": "S46_T4", "x": [2.062556887551913, -32.75476247924082, 1281.366888264051]},
            {"name": "S46_T5", "x": [1.886022797859931, -29.72766718766461, 1254.093375886867]},
            {"name": "S46_T6", "x": [2.665733008702827, -27.9774378890066, 1228.895536678729]},
            {"name": "S46_T7", "x": [2.685443219545767, -28.88590384550561, 1201.321310652545]},
            {"name": "S46_T8", "x": [3.673267205492586, -32.83587231921514, 1175.918284093885]},
            {"name": "S46_T9", "x": [3.238227336012321, -37.73926544971265, 1149.966858964759]},
            {"name": "S48_L1", "x": [4.828106519629028, -72.5090704337647, 1033.440878242578]},
            {"name": "S48_L2", "x": [4.760308675290831, -75.05770175769652, 999.1057557468655]},
            {"name": "S48_T10", "x": [3.758738734942221, -46.01045631270791, 1121.584441237124]},
            {"name": "S48_T11", "x": [3.758011452478777, -54.61099048574193, 1095.701402734419]},
            {"name": "S48_T12", "x": [3.319255753013748, -65.89377321281734, 1065.772034629317]},
            {"name": "S48_T2", "x": [-2.327000000000079, -46.93673438061327, 1325.346373555896]},
            {"name": "S48_T3", "x": [2.540721663321684, -39.40088033995292, 1305.483141882336]},
            {"name": "S48_T4", "x": [2.227625066936134, -33.60868282409854, 1282.36504346514]},
            {"name": "S48_T5", "x": [1.745363694780622, -30.92374684670274, 1255.328615631373]},
            {"name": "S48_T6", "x": [2.463977416335586, -29.23682629793244, 1229.787638196456]},
            {"name": "S48_T7", "x": [2.946488708167764, -31.48503342775221, 1202.706764256303]},
            {"name": "S48_T8", "x": [3.506215990631232, -34.55867694564652, 1177.776135009935]},
            {"name": "S48_T9", "x": [3.371488708167762, -38.27331306563198, 1151.451382778509]},
            {"name": "S49_L1", "x": [4.913585226232137, -70.09367508406223, 1029.57914856364]},
            {"name": "S49_L2", "x": [4.686691640573151, -72.60703912219866, 995.3723250569293]},
            {"name": "S49_T10", "x": [3.368215990631266, -44.49114468922804, 1117.580819951241]},
            {"name": "S49_T11", "x": [3.340409183402647, -53.7364776381595, 1090.969701393491]},
            {"name": "S49_T12", "x": [3.445051214861363, -65.1064044461588, 1062.501927792385]},
            {"name": "S49_T2", "x": [-1.861000000000106, -44.35330839820738, 1322.811140633451]},
            {"name": "S49_T3", "x": [1.910363694780666, -37.71228767058076, 1303.509429567448]},
            {"name": "S49_T4", "x": [2.177858022090072, -31.76683312870471, 1279.76935344474]},
            {"name": "S49_T5", "x": [1.926789842705966, -28.86236171231769, 1252.656359961557]},
            {"name": "S49_T6", "x": [3.103926201474179, -27.06967431988249, 1226.935802610672]},
            {"name": "S49_T7", "x": [3.05781825970739, -28.44793348467504, 1198.705071446394]},
            {"name": "S49_T8", "x": [3.381284170015407, -31.97490777641049, 1174.037749429856]},
            {"name": "S49_T9", "x": [3.13525575301394, -37.75141587415965, 1148.339644282405]},
            {"name": "S50_C3_B", "x": [1.417704584671582, -59.36980576095892, 1441.258897124201]},
            {"name": "S50_C3_T", "x": [1.083057904245435, -67.5233635899951, 1440.905515218172]},
            {"name": "S50_C4_B", "x": [1.172437405436072, -61.47353482843, 1427.305712331369]},
            {"name": "S50_C4_T", "x": [0.1482963413665204, -70.32908276288579, 1427.340255684711]},
            {"name": "S50_C5_B", "x": [2.351904656797745, -63.77099615841755, 1412.453988261626]},
            {"name": "S50_C5_T", "x": [2.073023798132263, -70.73283304295593, 1412.616567448561]},
            {"name": "S50_C6_B", "x": [4.323099919785121, -62.85938786440641, 1398.725139446278]},
            {"name": "S50_C7_B", "x": [3.923085645530189, -61.10329905038709, 1381.826426534944]},
            {"name": "S50_C8_B", "x": [1.959076047815395, -57.03538115932949, 1366.538734667652]},
            {"name": "S50_L1", "x": [4.513953918617396, -75.40793832586787, 1034.856167367138]},
            {"name": "S50_L1_B", "x": [4.410821236477061, -59.56620764864508, 1036.462656427869]},
            {"name": "S50_L2", "x": [5.079340589064327, -78.28967127077404, 1000.187803696305]},
            {"name": "S50_L2_B", "x": [5.061787381893915, -60.0983805414235, 1001.456324644062]},
            {"name": "S50_L5_B", "x": [4.317436330385054, -34.10878106158436, 909.2473794645805]},
            {"name": "S50_L5_T", "x": [5.531053365102559, -50.43876140340036, 909.2570130104522]},
            {"name": "S50_S1_B", "x": [5.461734227367216, -25.01946828821469, 887.6687429461322]},
            {"name": "S50_S1_T", "x": [5.649606572273125, -35.95197287607826, 886.6892755422007]},
            {"name": "S50_T1", "x": [2.11000164064043, -57.25201154144072, 1347.160938274842]},
            {"name": "S50_T10", "x": [3.602215990631292, -48.32660844424751, 1123.4652691734]},
            {"name": "S50_T11", "x": [4.015079631862904, -56.50344872268458, 1097.7358901097]},
            {"name": "S50_T12", "x": [3.151721663321792, -67.3033730202298, 1067.525785053932]},
            {"name": "S50_T1_B", "x": [1.958704565051177, -50.6962169414447, 1347.732357700654]},
            {"name": "S50_T2", "x": [3.077147811247113, -48.96690093313754, 1326.408787674909]},
            {"name": "S50_T2_B", "x": [3.218721663321694, -42.19914560238276, 1327.737230492698]},
            {"name": "S50_T3", "x": [2.796227336012417, -41.02679619555145, 1306.653099925999]},
            {"name": "S50_T3_B", "x": [2.884448945785177, -33.09491311640767, 1307.714562886412]},
            {"name": "S50_T4", "x": [2.290323932397984, -35.22378115063803, 1283.921849083613]},
            {"name": "S50_T4_B", "x": [2.45421599063115, -26.80078154287644, 1284.093259997403]},
            {"name": "S50_T5", "x": [1.756460291166312, -32.70055962317758, 1256.673267033398]},
            {"name": "S50_T5_B", "x": [1.675937546855313, -21.65437680219524, 1257.706967409594]},
            {"name": "S50_T6", "x": [2.657448945785319, -31.2305959856563, 1230.631495892052]},
            {"name": "S50_T6_B", "x": [2.513721663321824, -18.94289791332757, 1231.737852259781]},
            {"name": "S50_T7", "x": [3.309409183402628, -33.3666201916177, 1206.003615915446]},
            {"name": "S50_T8", "x": [3.488517125169293, -36.86307759049078, 1179.459390927572]},
            {"name": "S50_T9", "x": [3.527272824634327, -41.15121692766776, 1154.42360802674]},
            {"name": "T10_spinal", "x": [-0.8606022690761863, -41.560010432646, 1121.392568308061]},
            {"name": "T10_ventral_root_end", "x": [3.982698972559741, -56.85107211253722, 1129.450597547725]},
            {"name": "T11_spinal", "x": [-1.906999999999947, -50.48922937106168, 1095.445867064696]},
            {"name": "T11_ventral_root_end", "x": [4.502681900939168, -62.98244893803199, 1104.258328329053]},
            {"name": "T12_spinal", "x": [-1.284999999999956, -62.28033449578611, 1066.939107725524]},
            {"name": "T12_ventral_root_end", "x": [3.534545542170981, -75.19721133692644, 1074.410545704433]},
            {"name": "T1_dorsal_root_end", "x": [2.604655831037832, -47.53106032973373, 1351.839001603161]},
            {"name": "T1_spinal", "x": [-1.146999999999991, -54.91167120057669, 1345.3053323167]},
            {"name": "T1_ventral_root_end", "x": [2.124504005773964, -62.31045670141069, 1351.469890493544]},
            {"name": "T2_dorsal_root_end", "x": [3.359914856093053, -37.70927222458091, 1330.212092912388]},
            {"name": "T2_spinal", "x": [-1.518477255689072, -45.38357170691999, 1326.974685142999]},
            {"name": "T2_ventral_root_end", "x": [3.315108048864532, -55.11921498097823, 1329.600683956361]},
            {"name": "T3_dorsal_root_end", "x": [2.995448945785284, -29.39015650229117, 1310.052747793607]},
            {"name": "T3_spinal", "x": [-0.7746817938414013, -36.86770965150439, 1306.321088164786]},
            {"name": "T3_ventral_root_end", "x": [2.942613721555035, -47.10857805174037, 1309.844650416601]},
            {"name": "T4_dorsal_root_end", "x": [2.29018757362963, -21.8017126171507, 1286.683231248634]},
            {"name": "T4_spinal", "x": [-1.168039762382603, -30.88270820651022, 1284.174214701056]},
            {"name": "T4_ventral_root_end", "x": [2.213090977243966, -41.30390350398378, 1285.970400918695]},
            {"name": "T5_dorsal_root_end", "x": [1.841460291166267, -17.28644368834746, 1260.266901816675]},
            {"name": "T5_spinal", "x": [-2.181999999999947, -26.84341199456139, 1257.022325894247]},
            {"name": "T5_ventral_root_end", "x": [1.863130739626737, -39.52651268509298, 1259.78646978069]},
            {"name": "T6_dorsal_root_end", "x": [2.370488708167833, -12.45405858577584, 1236.755913479022]},
            {"name": "T6_spinal", "x": [-1.431000000000018, -25.39477069450544, 1231.041671203913]},
            {"name": "T6_ventral_root_end", "x": [2.356051214861469, -40.34938542283608, 1236.507002316686]},
            {"name": "T7_spinal", "x": [-0.9104090763047987, -26.61616160300789, 1205.270885648818]},
            {"name": "T7_ventral_root_end", "x": [3.599244407632693, -41.37995472125812, 1210.94121080918]},
            {"name": "T8_spinal", "x": [-1.120725825158564, -29.58954377503989, 1178.000316659255]},
            {"name": "T8_ventral_root_end", "x": [3.287448945785187, -42.86577822265541, 1184.702432817903]},
            {"name": "T9_spinal", "x": [-0.6351136680062436, -34.51006516724997, 1153.002957761351]},
            {"name": "T9_ventral_root_end", "x": [3.625323932397939, -49.34088568227134, 1160.409522611124]},
            {"name": "ardell_1_branching_point", "x": [59.89, -103.627, 1236.608]},
            {"name": "ardell_1_start-1", "x": [-0.744, -70.846, 1324.135]},
            {"name": "ardell_1_start-2", "x": [-2.344, -63.781, 1302.899]},
            {"name": "ardell_1_start-3", "x": [7.209000000000026, -59.37635073463456, 1270.988380639674]},
            {"name": "ardell_1_start-4", "x": [18.135, -55.372, 1249.411]},
            {"name": "ardell_1_start-5", "x": [2.194, -56.007, 1226.082]},
            {"name": "ardell_1_start-6", "x": [-6.74, -58.959, 1209.862]},
            {"name": "ardell_1_start-7", "x": [-8.395, -63.535, 1194.158]},
            {"name": "ardell_1_start-8", "x": [-3.28, -68.288, 1181.049]},
            {"name": "ardell_1_start-9", "x": [7.783, -75.385, 1161.554]},
            {"name": "ardell_2_branching_point", "x": [26.884, -83.413, 1247.799]},
            {"name": "ardell_2_start", "x": [-19.68108649796534, -82.19347448846338, 1301.463981661421]},
            {"name": "ardell_3_branching_point", "x": [37.934, -86.854, 1232.911]},
            {"name": "ardell_4_branching_point", "x": [14.9517344655693, -105.2010680107312, 1249.297516612408]},
            {"name": "ardell_5_branching_point", "x": [70.968, -122.424, 1222.93]},
            {"name": "ardell_5_start-1", "x": [11.364, -70.44, 1218.599]},
            {"name": "ardell_5_start-2", "x": [10.687, -87.712, 1199.718]},
            {"name": "ardell_5_start-3", "x": [19.007, -90.483, 1195.472]},
            {"name": "ardell_5_start-4", "x": [35.768, -97.093, 1186.35]},
            {"name": "ardell_5_start-5", "x": [36.889, -68.395, 1191.55]},
            {"name": "ardell_5_start-6", "x": [65.575, -68.869, 1182.642]},
            {"name": "ardell_6_branching_point", "x": [53.673, -91.566, 1219.179]},
            {"name": "ardell_6_start-1", "x": [2.345, -93.33, 1386.558]},
            {"name": "ardell_6_start-2", "x": [8.159, -89.298, 1324.808]},
            {"name": "ardell_6_start-3", "x": [4.396, -74.224, 1280.529]},
            {"name": "ardell_6_start-4", "x": [8.427, -66.618, 1253.756]},
            {"name": "ardell_6_start-5", "x": [7.187, -62.572, 1241.094]},
            {"name": "ardell_6_start-6", "x": [20.414, -55.322, 1222.845]},
            {"name": "ardell_6_start-7", "x": [31.316, -52.699, 1202.567]},
            {"name": "ardell_6_start-8", "x": [8.619, -71.316, 1195.269]},
            {"name": "ardell_6_start-9", "x": [19.577, -78.843, 1173.876]},
            {"name": "bladder_n-bladder", "x": [-2.552386253899746, -130.4537400391762, 820.582949907351]},
            {"name": "bolser_14-1", "x": [-9.021, -67.073, 1332.194]},
            {"name": "bolser_1_branching_point", "x": [17.504, -86.497, 1252.626]},
            {"name": "bolser_6_branching_point", "x": [-15.705, -81.961, 1290.933]},
            {"name": "bolser_8_branching_point", "x": [12.604, -73.595, 1437.216]},
            {"name": "bolser_9-2", "x": [-8.915, -53.794, 1257.483]},
            {"name": "bolser_9_branching_point", "x": [27.594, -137.354, 1472.661]},
            {"name": "brain_12-1", "x": [0.5259999999999295, -56.01677892707869, 1488.845221072921]},
            {"name": "brain_33-1", "x": [2.607, -56.401, 1489.594]},
            {"name": "brain_34-1", "x": [1.878, -56.261, 1489.3]},
            {"name": "brain_44-1", "x": [-3.068, -58.888, 1485.941]},
            {"name": "cardio_4-1", "x": [-0.4535302471967846, -71.31111190790045, 1420.283035707635]},
            {"name": "connection_1-1", "x": [-0.7699999999999119, -65.75869116574363, 1277.609317426094]},
            {"name": "digestive_9-1", "x": [4.050999999999929, -88.41149897649373, 1128.013983486292]},
            {"name": "external_laryngeal_n_branching_point",
             "x": [0.01456026895135378, -94.55048029904339, 1425.031631366861]},
            {"name": "ganglion_1-1", "x": [19.03884588693559, -71.172, 1471.286580044339]},
            {"name": "ganglion_12-1", "x": [-3.781999999999947, -44.49158106445484, 1468.081013610422]},
            {"name": "ganglion_2-1", "x": [22.94426222963473, -78.35900000000002, 1463.789969692503]},
            {"name": "ganglion_3-1", "x": [-10.70395637423707, -69.98163604107317, 1449.718866793385]},
            {"name": "ganglion_4-1", "x": [-4.026772846292005, -94.34455037786651, 878.2933909732724]},
            {"name": "ganglion_5-1", "x": [-1.041999999999903, -96.97343596359889, 1011.520471864983]},
            {"name": "keast_1-1", "x": [-1.980000000000088, -112.2919311293448, 847.5093359887412]},
            {"name": "keast_2-1", "x": [-8.908000000000097, -54.5411333598874, 846.0150700655889]},
            {"name": "keast_4-1", "x": [2.331999999999947, -62.09159901786794, 875.0749854330861]},
            {"name": "keast_5-1", "x": [-2.583000000000097, -96.94482837478135, 935.3696184131826]},
            {"name": "keast_7-1", "x": [-5.255637692019879, -89.20502365869832, 1013.450699297463]},
            {"name": "label_1-1", "x": [19.059, -76.211, 1480.118]},
            {"name": "label_10-1", "x": [3.085186791877135, -109.0680000000001, 1408.86715377157]},
            {"name": "label_11-1", "x": [0.755, -107.608, 1406.119]},
            {"name": "label_12-1", "x": [0.392, -105.992, 1406.523]},
            {"name": "label_13-1", "x": [12.80657789560075, -107.6021924208387, 1396.149541186971]},
            {"name": "label_14-1", "x": [6.623, -108.164, 1402.941]},
            {"name": "label_15-1", "x": [15.735, -106.477, 1387.002]},
            {"name": "label_16-1", "x": [16.942, -106.518, 1385.462]},
            {"name": "label_17-1", "x": [18.305, -110.622, 1389.761]},
            {"name": "label_18-1", "x": [17.812, -111.228, 1390.543]},
            {"name": "label_19-1", "x": [20.305, -109.662, 1385.525]},
            {"name": "label_2-1", "x": [-65.468, -70.375, 1520.633]},
            {"name": "label_20-1", "x": [19.183, -111.35, 1388.633]},
            {"name": "label_21-1", "x": [18.88, -110.329, 1387.243]},
            {"name": "label_22-1", "x": [19.598, -110.678, 1389.008]},
            {"name": "label_3-1", "x": [-57.275, -70.875, 1519.911]},
            {"name": "label_4-1", "x": [60.913, -52.711, 1519.712]},
            {"name": "label_5-1", "x": [9.582, -154.877, 1443.576]},
            {"name": "label_6-1", "x": [-14.0, -152.355, 1446.593]},
            {"name": "label_7-1", "x": [-0.539, -131.251, 1434.71]},
            {"name": "label_8-1", "x": [9.361, -130.529, 1433.801]},
            {"name": "label_9-1", "x": [12.293, -108.004, 1394.644]},
            {"name": "lumbar_splanchnic_n_start", "x": [17.12899999999994, -84.43558340690586, 1011.200213818207]},
            {"name": "pelvic_splanchnic_n_start", "x": [11.49799999999992, -42.82028621808127, 878.3340180993188]},
            {"name": "phrenic_n_branching_point", "x": [3.405000000000088, -127.2422887569768, 1151.854221194133]},
            {"name": "plexus_1-1", "x": [-8.381000000000018, -82.82003422210069, 1438.486785263734]},
            {"name": "plexus_2-1", "x": [5.250738332356949, -85.69587432205755, 1451.118102341404]},
            {"name": "plexus_3-1", "x": [1.754222381762545, -82.3427353211713, 1313.597783594231]},
            {"name": "plexus_4-1", "x": [2.219394702024564, -78.51945223830805, 1287.385118563105]},
            {"name": "plexus_5-1", "x": [-0.4969999999999912, -82.22631192321117, 1182.07848441355]},
            {"name": "plexus_6-1", "x": [2.250999999999929, -78.17836681127925, 1099.527903904539]},
            {"name": "plexus_7-1", "x": [-1.112, -68.774, 1442.445]},
            {"name": "point_1", "x": [-0.628, -62.693, 1481.347]},
            {"name": "point_10", "x": [-6.224, -58.028, 1472.469]},
            {"name": "point_11", "x": [34.55985130071155, -80.58799999999997, 1457.495277881116]},
            {"name": "point_12", "x": [8.015, -99.059, 1456.929]},
            {"name": "point_13", "x": [9.69009742357528, -116.5569999999999, 1453.666701667099]},
            {"name": "point_14", "x": [35.70300000000001, -82.40857087649809, 1448.810785438249]},
            {"name": "point_15", "x": [14.45262938065343, -75.50319122623176, 1446.741138247483]},
            {"name": "point_16", "x": [34.972, -83.394, 1441.867]},
            {"name": "point_17", "x": [29.51560153559797, -87.31299999999996, 1439.653814690327]},
            {"name": "point_18", "x": [20.44, -92.16399999999989, 1438.234069931184]},
            {"name": "point_19", "x": [14.327, -97.34900000000007, 1437.758038621684]},
            {"name": "point_2", "x": [18.58492862694576, -70.79974298662209, 1472.165829671721]},
            {"name": "point_20", "x": [6.253132108443852, -105.5393132326375, 1435.677300212148]},
            {"name": "point_21", "x": [8.351303122764456, -99.87899999999993, 1430.510350682337]},
            {"name": "point_22", "x": [34.03838372377241, -84.45099999999994, 1435.256308138114]},
            {"name": "point_23", "x": [36.305, -89.725, 1426.744]},
            {"name": "point_24", "x": [29.77884882868279, -91.62321601409755, 1413.811364373984]},
            {"name": "point_25", "x": [27.68499999999996, -92.26284527408687, 1404.556090717768]},
            {"name": "point_26", "x": [32.29000000000004, -92.57137916848482, 1395.424652344061]},
            {"name": "point_27", "x": [17.23400000000003, -95.99695396187428, 1395.634502139794]},
            {"name": "point_28", "x": [30.96200000000004, -99.1679757611211, 1394.806093348933]},
            {"name": "point_29", "x": [30.62100000000006, -102.4364290866547, 1392.843515222423]},
            {"name": "point_3", "x": [17.463, -71.807, 1476.875]},
            {"name": "point_30", "x": [28.921, -109.629, 1316.751]},
            {"name": "point_31", "x": [28.21282515516866, -112.3540307459024, 1304.945825155169]},
            {"name": "point_32", "x": [5.604583350861463, -109.3295179037709, 1270.878821384835]},
            {"name": "point_33", "x": [4.12624009902753, -85.33741296632425, 1205.089241105896]},
            {"name": "point_34", "x": [4.490000000000044, -101.2477922228036, 1168.426855785408]},
            {"name": "point_35", "x": [-4.484999999999956, -92.54011426247118, 1123.925189560716]},
            {"name": "point_36", "x": [37.35727812716359, -83.31845214668414, 1093.807605465439]},
            {"name": "point_37", "x": [12.367, -99.729, 1421.523]},
            {"name": "point_38", "x": [30.626, -102.783, 1394.211]},
            {"name": "point_4", "x": [18.55859023261447, -70.4983231698189, 1470.56637824997]},
            {"name": "point_40", "x": [5.571000000000062, -55.17470764842699, 1229.232824589056]},
            {"name": "point_41", "x": [87.6319508924071, -70.72918834780654, 929.7169821909588]},
            {"name": "point_42", "x": [-0.5430284170014598, -64.84088305937085, 1232.798940819478]},
            {"name": "point_5", "x": [22.97196468527639, -76.79900000000006, 1461.98689405583]},
            {"name": "point_6", "x": [15.862, -66.876, 1461.22]},
            {"name": "point_7", "x": [14.086, -75.333, 1461.325]},
            {"name": "point_8", "x": [4.653, -62.902, 1463.956]},
            {"name": "point_9", "x": [12.273, -59.511, 1466.982]},
            {"name": "pudendal_n_branching_point", "x": [-21.196, -81.676, 818.357]},
            {"name": "pudendal_n_start", "x": [0.2599999999999559, -35.95073918975689, 874.720608643658]},
            {"name": "respiratory_11-1", "x": [0.6035592810405546, -88.41365342813772, 1424.11677450975]},
            {"name": "respiratory_17-1", "x": [0.7606351837905805, -80.50398602753384, 1208.232427391316]},
            {"name": "respiratory_5-1", "x": [0.02602457179082696, -102.1063383681233, 1412.804349526629]},
            {"name": "respiratory_8-1", "x": [0.3785760961309229, -108.913, 1399.902314261826]},
            {"name": "spinal_47-1", "x": [4.893000000000053, -80.41561418791244, 1058.76329151085]},
            {"name": "splanchnic_n_branching_point", "x": [40.723, -164.527, 1031.848]},
            {"name": "splanchnic_n_start-1", "x": [6.72599999999993, -67.46669873058676, 1341.139422819567]},
            {"name": "splanchnic_n_start-10", "x": [2.013943165996861, -67.5331029166776, 1118.477820719743]},
            {"name": "splanchnic_n_start-11", "x": [4.436000000000106, -75.36208889012215, 1092.620020009017]},
            {"name": "splanchnic_n_start-12", "x": [0.06399999999989422, -84.55919291500041, 1064.352927792975]},
            {"name": "splanchnic_n_start-13", "x": [4.291999999999903, -87.87369105074335, 1053.10832944315]},
            {"name": "splanchnic_n_start-14", "x": [5.22400000000007, -20.60425391745015, 1014.692611259957]},
            {"name": "splanchnic_n_start-2", "x": [4.118000000000053, -61.43837636946041, 1320.788821409292]},
            {"name": "splanchnic_n_start-3", "x": [5.425, -52.72849084736327, 1300.256409723436]},
            {"name": "splanchnic_n_start-4", "x": [4.659999999999956, -48.35448841843786, 1277.491314174933]},
            {"name": "splanchnic_n_start-5", "x": [4.990000000000044, -44.97291851826787, 1252.039098360451]},
            {"name": "splanchnic_n_start-6", "x": [5.759000000000026, -43.60629103275866, 1227.926782398807]},
            {"name": "splanchnic_n_start-7", "x": [4.10099999999993, -42.60459214145332, 1200.490711095613]},
            {"name": "splanchnic_n_start-8", "x": [5.184000000000027, -48.08531162914964, 1173.843307306872]},
            {"name": "splanchnic_n_start-9", "x": [5.912999999999965, -54.36318694445057, 1148.254579674199]},
            {"name": "urinary_13-1", "x": [79.83627524840969, -51.66973364795776, 1070.676874605029]},
        ]
        # annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_nerve_term(termName))
        # annotationGroup.createMarkerNode(nodeIdentifier, bodyCoordinates, nerveCoordinatesValues)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        findMarkerLocation = fieldmodule.createFieldFindMeshLocation(markerBodyCoordinates, coordinates, mesh)
        findMarkerLocation.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_EXACT)
        for bodyMarkerPoint in bodyMarkerPoints:
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            fieldcache.setNode(markerPoint)
            markerBodyCoordinates.assignReal(fieldcache, bodyMarkerPoint["x"])
            markerName.assignString(fieldcache, bodyMarkerPoint["name"])

            element, xi = findMarkerLocation.evaluateMeshLocation(fieldcache, 3)
            print(element.getIdentifier(), xi)
            markerLocation.assignMeshLocation(fieldcache, element, xi)
            nodeIdentifier += 1
        sir = region.createStreaminformationRegion()
        # srm = sir.createStreamresourceMemory()
        # sir.setResourceGroupName(srm, organ_name)
        # sir.setResourceFieldNames(srm, fieldNames)
        with ChangeManager(fieldmodule):
            region.write(sir)
            region.writeFile(r'C:\Users\egha355\Desktop\sparc3\codes\mapclient_workflows\workflowandfilesnewer2\whole-body222222211111.exf')




        bodyGroup = AnnotationGroup(region, get_body_term("body"))
        coreGroup = AnnotationGroup(region, get_body_term("core"))
        non_coreGroup = AnnotationGroup(region, get_body_term("non core"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        neckGroup = AnnotationGroup(region, get_body_term("neck core"))
        headGroup = AnnotationGroup(region, get_body_term("head core"))
        annotationGroups = [bodyGroup, coreGroup, non_coreGroup, abdomenGroup, thoraxGroup, neckGroup, headGroup]

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            elementsCountAcrossTransition, shellProportion,
                            [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fieldmodule, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # body coordinates
        bodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="body coordinates")
        tmp_region = region.createRegion()
        tmp_fieldmodule = tmp_region.getFieldmodule()
        tmp_body_coordinates = findOrCreateFieldCoordinates(tmp_fieldmodule, name="body coordinates")
        tmp_cylinder = CylinderMesh(tmp_fieldmodule, tmp_body_coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)
        sir = tmp_region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        tmp_region.write(sir)
        result, buffer = srm.getBuffer()
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemoryBuffer(buffer)
        region.read(sir)

        del srm
        del sir
        del tmp_body_coordinates
        del tmp_fieldmodule
        del tmp_region

        # Groups of different parts of the body
        is_body = fieldmodule.createFieldConstant(1)
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        bodyMeshGroup.addElementsConditional(is_body)

        coreMeshGroup = coreGroup.getMeshGroup(mesh)

        # core group
        e1a = elementsCountAcrossShell
        e1z = elementsCountAcrossMinor - elementsCountAcrossShell - 1
        e2a = elementsCountAcrossShell
        e2b = e2a + elementsCountAcrossTransition
        e2z = elementsCountAcrossMajor - elementsCountAcrossShell - 1
        e2y = e2z - elementsCountAcrossTransition
        e1oc = elementsCountAcrossMinor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oc = elementsCountAcrossMajor - 2*elementsCountAcrossShell - 2*elementsCountAcrossTransition
        e2oCore = e2oc * e1oc + 2 * elementsCountAcrossTransition * (e2oc + e1oc)
        elementsCountAround = cylinder1.getElementsCountAround()
        e2oShell = elementsCountAround * elementsCountAcrossShell
        e2o = e2oCore + e2oShell
        elementId = cylinder1.getElementIdentifiers()
        for e3 in range(elementsCountAlong):
            for e2 in range(elementsCountAcrossMajor):
                for e1 in range(elementsCountAcrossMinor):
                    coreElement = ((e2 >= e2a) and (e2 <= e2z)) and ((e1 >= e1a) and (e1 <= e1z))
                    if coreElement:
                        elementIdentifier = elementId[e3][e2][e1]
                        if elementIdentifier:
                            element = mesh.findElementByIdentifier(elementIdentifier)
                            coreMeshGroup.addElement(element)

        is_non_core = fieldmodule.createFieldNot(coreGroup.getGroup())
        non_coreMeshGroup = non_coreGroup.getMeshGroup(mesh)
        non_coreMeshGroup.addElementsConditional(is_non_core)

        abdomenMeshGroup = abdomenGroup.getMeshGroup(mesh)
        thoraxMeshGroup = thoraxGroup.getMeshGroup(mesh)
        neckMeshGroup = neckGroup.getMeshGroup(mesh)
        headMeshGroup = headGroup.getMeshGroup(mesh)
        meshGroups = [abdomenMeshGroup, thoraxMeshGroup, neckMeshGroup, headMeshGroup]

        abdomenRange = [1, elementsCountAlongAbdomen*e2o]
        thoraxRange = [abdomenRange[1]+1, abdomenRange[1]+elementsCountAlongThorax*e2o]
        neckRange = [thoraxRange[1]+1, thoraxRange[1] + elementsCountAlongNeck*e2o]
        headRange = [neckRange[1]+1, elementsCountAlong*e2o]
        groupsRanges = [abdomenRange, thoraxRange, neckRange, headRange]

        totalElements = e2o*elementsCountAlong
        for elementIdentifier in range(1, totalElements+1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            if coreMeshGroup.containsElement(element):
                ri = 0
                for groupRange in groupsRanges:
                    if (elementIdentifier >= groupRange[0]) and (elementIdentifier <= groupRange[1]):
                        meshGroups[ri].addElement(element)
                        break
                    ri += 1

        if discontinuity:
            # create discontinuity in d3 on the core boundary
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            elementtemplate = mesh.createElementtemplate()
            undefineNodetemplate = nodes.createNodetemplate()
            undefineNodetemplate.undefineField(coordinates)
            nodetemplate = nodes.createNodetemplate()
            fieldcache = fieldmodule.createFieldcache()
            with ChangeManager(fieldmodule):
                localNodeIndexes = [1, 2, 3, 4]
                valueLabel = Node.VALUE_LABEL_D_DS3
                for e3 in range(elementsCountAlong):
                    for e2 in range(elementsCountAcrossMajor):
                        for e1 in range(elementsCountAcrossMinor):
                            regularRowElement = (((e2 >= e2b) and (e2 <= e2y)) and ((e1 == e1a - 1) or (e1 == e1z + 1)))
                            non_coreFirstLayerElement = (e2 == e2a - 1) or regularRowElement or (e2 == e2z + 1)
                            elementIdentifier = elementId[e3][e2][e1]
                            if elementIdentifier and non_coreFirstLayerElement:
                                element = mesh.findElementByIdentifier(elementIdentifier)
                                eft = element.getElementfieldtemplate(coordinates, -1)
                                nodeIds = get_element_node_identifiers(element, eft)
                                for localNodeIndex in localNodeIndexes:
                                    node = element.getNode(eft, localNodeIndex)
                                    nodetemplate.defineFieldFromNode(coordinates, node)
                                    versionsCount = nodetemplate.getValueNumberOfVersions(coordinates, -1, valueLabel)
                                    if versionsCount == 1:
                                        fieldcache.setNode(node)
                                        result0, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                                        result0, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                                        result0, d2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                                        result0, d3 = coordinates.getNodeParameters(fieldcache, -1, valueLabel, 1, 3)
                                        result1 = node.merge(undefineNodetemplate)
                                        result2 = nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 2)
                                        result3 = node.merge(nodetemplate)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                                        result4 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 1, d3)
                                        result5 = coordinates.setNodeParameters(fieldcache, -1, valueLabel, 2, d3)
                                remapEftNodeValueLabelsVersion(eft, localNodeIndexes, [valueLabel], 2)
                                result1 = elementtemplate.defineField(coordinates, -1, eft)
                                result2 = element.merge(elementtemplate)
                                result3 = element.setNodesByIdentifier(eft, nodeIds)
        else:
            fieldcache = fieldmodule.createFieldcache()

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fieldmodule, "marker")
        markerName = findOrCreateFieldStoredString(fieldmodule, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fieldmodule, mesh, name="marker_location")
        markerBodyCoordinates = findOrCreateFieldCoordinates(fieldmodule, name="marker_body_coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        markerTemplateInternal.defineField(markerBodyCoordinates)
        #
        middleLeft = elementsCountAcrossMinor//2
        topElem = elementsCountAcrossMajor - 1
        middleRight = middleLeft - 1
        neckFirstElem = elementsCountAlongAbdomen+elementsCountAlongThorax
        thoraxFirstElem = elementsCountAlongAbdomen
        middleDown = elementsCountAcrossMajor//2 - 1

        # organ landmarks groups
        apexOfHeart = heart_terms.get_heart_term('apex of heart')
        leftAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('left atrium epicardium venous midpoint')
        rightAtriumEpicardiumVenousMidpoint = heart_terms.get_heart_term('right atrium epicardium venous midpoint')
        apexOfUrinaryBladder = bladder_terms.get_bladder_term('apex of urinary bladder')
        leftUreterJunctionWithBladder = bladder_terms.get_bladder_term('left ureter junction with bladder')
        rightUreterJunctionWithBladder = bladder_terms.get_bladder_term('right ureter junction with bladder')
        urethraJunctionWithBladderDorsal = bladder_terms.get_bladder_term('urethra junction of dorsal bladder neck')
        urethraJunctionWithBladderVentral = bladder_terms.get_bladder_term('urethra junction of ventral bladder neck')
        gastroesophagalJunctionOnLesserCurvature = stomach_terms.get_stomach_term('esophagogastric junction along the lesser curvature on serosa')
        limitingRidgeOnGreaterCurvature = stomach_terms.get_stomach_term('limiting ridge at the greater curvature on serosa')
        pylorusOnGreaterCurvature = stomach_terms.get_stomach_term('gastroduodenal junction along the greater curvature on serosa')
        junctionBetweenFundusAndBodyOnGreaterCurvature = stomach_terms.get_stomach_term("fundus-body junction along the greater curvature on serosa")
        apexOfLeftLung = lung_terms.get_lung_term('apex of left lung')
        ventralBaseOfLeftLung = lung_terms.get_lung_term('ventral base of left lung')
        dorsalBaseOfLeftLung = lung_terms.get_lung_term('dorsal base of left lung')
        apexOfRightLung = lung_terms.get_lung_term('apex of right lung')
        ventralBaseOfRightLung = lung_terms.get_lung_term('ventral base of right lung')
        dorsalBaseOfRightLung = lung_terms.get_lung_term('dorsal base of right lung')
        laterodorsalTipOfMiddleLobeOfRightLung = lung_terms.get_lung_term('laterodorsal tip of middle lobe of right lung')
        apexOfRightLungAccessoryLobe = lung_terms.get_lung_term('apex of right lung accessory lobe')
        ventralBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('ventral base of right lung accessory lobe')
        dorsalBaseOfRightLungAccessoryLobe = lung_terms.get_lung_term('dorsal base of right lung accessory lobe')
        medialBaseOfLeftLung = lung_terms.get_lung_term("medial base of left lung")
        medialBaseOfRightLung = lung_terms.get_lung_term("medial base of right lung")
        brainstemDorsalMidlineCaudalPoint = brainstem_terms.get_brainstem_term('brainstem dorsal midline caudal point')
        brainstemDorsalMidlineCranialPoint = brainstem_terms.get_brainstem_term('brainstem dorsal midline cranial point')
        brainstemVentralMidlineCaudalPoint = brainstem_terms.get_brainstem_term('brainstem ventral midline caudal point')
        brainstemVentralMidlineCranialPoint = brainstem_terms.get_brainstem_term('brainstem ventral midline cranial point')

        # marker coordinates. In future we want to have only one table for all species.
        if isRat:
            bodyMarkerPoints = [
                {"group": ("left hip joint", ''), "x": [0.367, 0.266, 0.477]},
                {"group": ("right hip joint", ''), "x": [-0.367, 0.266, 0.477]},
                {"group": ("left shoulder joint", ''), "x": [0.456, -0.071, 2.705]},
                {"group": ("right shoulder joint", ''), "x": [-0.456, -0.071, 2.705]},
                {"group": ("along left femur", ''), "x": [0.456, 0.07, 0.633]},
                {"group": ("along right femur", ''), "x": [-0.456, 0.07, 0.633]},
                {"group": ("along left humerus", ''), "x": [0.423, -0.173, 2.545]},
                {"group": ("along right humerus", ''), "x": [-0.423, -0.173, 2.545]},
                {"group": apexOfUrinaryBladder, "x": [-0.124, -0.383, 0.434]},
                {"group": leftUreterJunctionWithBladder, "x": [-0.111, -0.172, 0.354]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.03, -0.196, 0.363]},
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.03, -0.26, 0.209]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.037, -0.304, 0.203]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [-0.032, 0.418, 2.713]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [-0.017, 0.203, 2.941]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [-0.028, 0.388, 2.72]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [-0.019, 0.167, 2.95]},
                {"group": apexOfHeart, "x": [0.096, -0.128, 1.601]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.127, -0.083, 2.079]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [0.039, -0.082, 2.075]},
                {"group": apexOfLeftLung, "x": [0.172, -0.175, 2.337]},
                {"group": ventralBaseOfLeftLung, "x": [0.274, -0.285, 1.602]},
                {"group": dorsalBaseOfLeftLung, "x": [0.037, 0.31, 1.649]},
                {"group": apexOfRightLung, "x": [-0.086, -0.096, 2.311]},
                {"group": ventralBaseOfRightLung, "x": [0.14, -0.357, 1.662]},
                {"group": dorsalBaseOfRightLung, "x": [-0.054, 0.304, 1.667]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.258, -0.173, 2.013]},
                {"group": apexOfRightLungAccessoryLobe, "x": [0.041, -0.063, 1.965]},
                {"group": ventralBaseOfRightLungAccessoryLobe, "x": [0.143, -0.356, 1.66]},
                {"group": dorsalBaseOfRightLungAccessoryLobe, "x": [0.121, -0.067, 1.627]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [0.12, 0.009, 1.446]},
                {"group": limitingRidgeOnGreaterCurvature, "x": [0.318, 0.097, 1.406]},
                {"group": pylorusOnGreaterCurvature, "x": [0.08, -0.111, 1.443]},
            ]
        elif isHuman:
            bodyMarkerPoints = [
                {"group": urethraJunctionWithBladderDorsal, "x": [-0.0071, -0.2439, 0.1798]},
                {"group": urethraJunctionWithBladderVentral, "x": [-0.007, -0.2528, 0.1732]},
                {"group": leftUreterJunctionWithBladder, "x": [0.1074, 0.045, 0.1728]},
                {"group": rightUreterJunctionWithBladder, "x": [-0.1058, 0.0533, 0.1701]},
                {"group": apexOfUrinaryBladder, "x": [0.005, 0.1286, 0.1264]},
                {"group": brainstemDorsalMidlineCaudalPoint, "x": [0.0068, 0.427, 2.7389]},
                {"group": brainstemDorsalMidlineCranialPoint, "x": [0.008, -0.0231, 3.0778]},
                {"group": brainstemVentralMidlineCaudalPoint, "x": [0.0054, 0.3041, 2.7374]},
                {"group": brainstemVentralMidlineCranialPoint, "x": [0.0025, -0.2308, 3.091]},
                {"group": apexOfHeart, "x": [0.1373, -0.1855, 1.421]},
                {"group": leftAtriumEpicardiumVenousMidpoint, "x": [0.0024, 0.1452, 1.8022]},
                {"group": rightAtriumEpicardiumVenousMidpoint, "x": [-0.0464, 0.0373, 1.7491]},
                {"group": apexOfLeftLung, "x": [0.0655, -0.0873, 2.3564]},
                {"group": apexOfRightLung, "x": [-0.088, -0.0363, 2.3518]},
                {"group": laterodorsalTipOfMiddleLobeOfRightLung, "x": [-0.2838, -0.0933, 1.9962]},
                {"group": ventralBaseOfLeftLung, "x": [0.219, -0.2866, 1.4602]},
                {"group": medialBaseOfLeftLung, "x": [0.0426, -0.0201, 1.4109]},
                {"group": ventralBaseOfRightLung, "x": [-0.2302, -0.2356, 1.3926]},
                {"group": medialBaseOfRightLung, "x": [-0.0363, 0.0589, 1.3984]},
                {"group": dorsalBaseOfLeftLung, "x": [0.1544, 0.2603, 1.3691]},
                {"group": dorsalBaseOfRightLung, "x": [0.0369, -0.2524, 0.912]},
                {"group": gastroesophagalJunctionOnLesserCurvature, "x": [-0.0062, -0.3259, 0.8586]},
                {"group": pylorusOnGreaterCurvature, "x": [-0.0761, -0.3189, 0.8663]},
                {"group": junctionBetweenFundusAndBodyOnGreaterCurvature, "x": [0.1884, -0.1839, 0.9639]},
            ]

        nodeIdentifier = cylinder1._endNodeIdentifier
        findMarkerLocation = fieldmodule.createFieldFindMeshLocation(markerBodyCoordinates, bodyCoordinates, mesh)
        findMarkerLocation.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_EXACT)
        for bodyMarkerPoint in bodyMarkerPoints:
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            fieldcache.setNode(markerPoint)
            markerBodyCoordinates.assignReal(fieldcache, bodyMarkerPoint["x"])
            markerName.assignString(fieldcache, bodyMarkerPoint["group"][0])

            element, xi = findMarkerLocation.evaluateMeshLocation(fieldcache, 3)
            markerLocation.assignMeshLocation(fieldcache, element, xi)
            nodeIdentifier += 1

        # innnn  k,jjjd
        markerTermNameNerveCoordinatesMap = {
            'Gray communicating ramus': [-0.0071, -0.2439, 0.1798]
        }
        # bodyCoordinates = findOrCreateFieldCoordinates(fieldmodule)
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        for termName, nerveCoordinatesValues in markerTermNameNerveCoordinatesMap.items():
            annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_nerve_term(termName))
            annotationGroup.createMarkerNode(nodeIdentifier, bodyCoordinates, nerveCoordinatesValues)
            nodeIdentifier += 1

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """

        # create 2d surface mesh groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("body"))
        coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))
        non_coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("non core"))
        abdomenGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdomen"))
        thoraxGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thorax"))
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("neck core"))

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))
        coreBoundaryGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("core boundary"))
        diaphragmGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("diaphragm"))
        spinalCordGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("spinal cord"))

        is_exterior = fm.createFieldIsExterior()
        is_on_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_skin = fm.createFieldAnd(is_exterior, is_on_face_xi3_1)

        skinMeshGroup = skinGroup.getMeshGroup(mesh2d)
        skinMeshGroup.addElementsConditional(is_skin)

        is_core_boundary = fm.createFieldAnd(coreGroup.getGroup(), non_coreGroup.getGroup())
        coreBoundaryMeshGroup = coreBoundaryGroup.getMeshGroup(mesh2d)
        coreBoundaryMeshGroup.addElementsConditional(is_core_boundary)

        is_diaphragm = fm.createFieldAnd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        diaphragmMeshGroup = diaphragmGroup.getMeshGroup(mesh2d)
        diaphragmMeshGroup.addElementsConditional(is_diaphragm)

        # spinal cord
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()
        zero = fm.createFieldConstant(0)
        zero_m = fm.createFieldConstant(-0.01)
        zero_p = fm.createFieldConstant(0.01)
        comp2 = cls.axis2.index(max(cls.axis2)) + 1
        ax2_comp = fm.createFieldComponent(coordinates, comp2)
        ax2_gt_zero_m = fm.createFieldGreaterThan(ax2_comp, zero_m)
        ax2_lt_zero_p = fm.createFieldLessThan(ax2_comp, zero_p)
        ax2_gt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_gt_zero_m)
        ax2_lt_zero_xi10 = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0), ax2_lt_zero_p)
        is_ax2_zero = fm.createFieldAnd(ax2_lt_zero_xi10, ax2_gt_zero_xi10)
        comp3 = cls.axis3.index(max(cls.axis3)) + 1
        ax3_comp = fm.createFieldComponent(coordinates, comp3)
        ax3_positive = fm.createFieldGreaterThan(ax3_comp, zero)
        is_ax2_zero_ax3_positive = fm.createFieldAnd(is_ax2_zero, ax3_positive)
        is_abdomen_thorax = fm.createFieldAdd(abdomenGroup.getGroup(), thoraxGroup.getGroup())
        is_abdomen_thorax_neck = fm.createFieldAdd(is_abdomen_thorax, neckGroup.getGroup())
        is_abdomen_thorax_neck_boundary = fm.createFieldAnd(is_core_boundary, is_abdomen_thorax_neck)
        is_spinal_cord = fm.createFieldAnd(is_ax2_zero_ax3_positive, is_abdomen_thorax_neck_boundary)

        mesh1d = fm.findMeshByDimension(1)
        spinalCordMeshGroup = spinalCordGroup.getMeshGroup(mesh1d)
        spinalCordMeshGroup.addElementsConditional(is_spinal_cord)
