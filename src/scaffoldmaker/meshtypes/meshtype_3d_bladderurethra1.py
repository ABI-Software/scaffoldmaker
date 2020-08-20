'''
Generates 3D bladder and urethra meshes along the central path,
with variable numbers of elements around, along and through
wall.
'''

import copy
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues, mesh_destroy_elements_and_nodes_by_identifiers
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_bladderurethra1(Scaffold_base):
    '''
    Generates 3D bladder and urethra meshes with variable numbers
    of elements around, along the central path, and through the wall.
    '''
    centralPathDefaultScaffoldPackages_LUT = {
        'Cat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 8
                },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 15.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 15.0], [0.0, 0.0, 15.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 30.0], [0.0, 0.0, 15.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 45.0], [0.0, 0.0, 15.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.07, 1.2, 60.4], [0.4, 1.4, 13.5], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.3, 3.2, 75.2], [0.5, 1.5, 12.7], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.6, 5.6, 90.2], [0.5, 1.8, 13.2], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[1.0, 8.7, 105.0], [0.7, 3.7, 12.9], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[1.7, 13.9, 119.2], [0.9, 5.5, 13.3], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]]])
            }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 10.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 10.0], [0.0, 0.0, 10.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 20.0], [0.0, 0.0, 10], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 30.0], [0.0, 0.0, 10], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 40.0], [-0.4, 0.3, 7.7], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[-0.5, -2.1, 49.6], [-1.0, -1.9, 7.2], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[-0.9, -3.8, 59.3], [-1.4, -0.6, 8.8], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[-1.3, -3.6, 71.5], [-0.1, 2.0, 8.4], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]],
                    [[0.0, 0.0, 80.0], [0.9, 4.0, 5.5], [0.0, 0.5, 0.0], [0.0, 0.0, -0.5]]])
            })
        }
    ostiumDefaultScaffoldPackages = {
        'Ostium Cat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 2,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 2.0,
                'Ostium length': 1.0,
                'Ostium wall thickness': 0.5,
                'Use linear through ostium wall': True,
                'Vessel inner diameter': 0.8,
                'Vessel wall thickness': 0.25,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Ostium Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 1.0,
                'Ostium length': 0.5,
                'Ostium wall thickness': 0.02,
                'Use linear through ostium wall': True,
                'Vessel inner diameter': 0.4,
                'Vessel wall thickness': 0.1,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        })
    }

    @staticmethod
    def getName():
        return '3D Bladder with Urethra 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cat 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Rat 1' in parameterSetName:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Rat 1']
        else:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Cat 1']
        if 'Rat 1' in parameterSetName:
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Ostium Rat 1']
        else:
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Ostium Cat 1']
        options = {
            'Central path LUT': copy.deepcopy(centralPathOption_LUT),
            'Number of elements along': 20,
            'Number of elements around': 12,
            'Number of elements through wall': 1,
            'Major diameter': 30.0,
            'Minor diameter': 25.0,
            'Neck diameter 1': 5.0,
            'Neck diameter 2': 3.0,
            'Bladder wall thickness': 0.5,
            'Neck angle degrees': 45,
            'Include urethra': True,
            'Urethra diameter 1': 1.5,
            'Urethra diameter 2': 1.0,
            'Urethra wall thickness': 0.5,
            'Length factor': 0.5,
            'Include ureter': False,
            'Ureter': copy.deepcopy(ostiumOption),
            'Ostium position around': 0.65,  # should be on the dorsal part (> 0.5)
            'Ostium position down': 0.75,
            'Number of elements radially on annulus': 1,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements along': 4,
            'Refine number of elements through wall': 1
            }
        if 'Rat' in parameterSetName:
            options['Major diameter'] = 20.0
            options['Minor diameter'] = 10.0
            options['Bladder wall thickness'] = 0.2
            options['Neck diameter 1'] = 3.5
            options['Neck diameter 2'] = 2.0
            options['Ostium position around'] = 0.65  # should be on the dorsal part (> 0.5)
            options['Ostium position down'] = 0.75
            options['Urethra diameter 1'] = 0.75
            options['Urethra diameter 2'] = 0.65
            options['Urethra wall thickness'] = 0.25
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Central path LUT',
            'Number of elements along',
            'Number of elements around',
            'Number of elements through wall',
            'Major diameter',
            'Minor diameter',
            'Neck diameter 1',
            'Neck diameter 2',
            'Bladder wall thickness',
            'Neck angle degrees',
            'Include urethra',
            'Urethra diameter 1',
            'Urethra diameter 2',
            'Urethra wall thickness',
            'Length factor',
            'Include ureter',
            'Ureter',
            'Ostium position around',
            'Ostium position down',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path LUT':
            return [MeshType_1d_path1]
        if optionName == 'Ureter':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path LUT':
            return list(cls.centralPathDefaultScaffoldPackages_LUT.keys())
        if optionName == 'Ureter':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
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
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path LUT':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages_LUT.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages_LUT[parameterSetName])
        if optionName == 'Ureter':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path LUT'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path LUT'):
            options['Central path LUT'] = cls.getOptionScaffoldPackage('Central path LUT', MeshType_1d_path1)
        for key in [
            'Number of elements along',
            'Number of elements around',
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements through wall'] > 1:
            if options['Include ureter']:
                options['Number of elements through wall'] = 1

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        '''
        Update ostium sub-scaffold options which depend on parent options.
        '''
        bladderWallThickness = options['Bladder wall thickness']
        ostiumOptions = options['Ureter']
        ostiumDefaultOptions = ostiumOptions.getScaffoldSettings()
        ostiumDefaultOptions['Ostium wall thickness'] = bladderWallThickness

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        '''
        cls.updateSubScaffoldOptions(options)
        elementsCountAlong = options['Number of elements along']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        majorDiameter = options['Major diameter']
        minorDiameter = options['Minor diameter']
        neckDiameter1 = options['Neck diameter 1']
        neckDiameter2 = options['Neck diameter 2']
        bladderWallThickness = options['Bladder wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        neckAngleRadians = math.radians(options['Neck angle degrees'])

        includeUrethra = options['Include urethra']
        urethraDiameter1 = options['Urethra diameter 1']
        urethraDiameter2 = options['Urethra diameter 2']
        urethraWallThickness = options['Urethra wall thickness']
        lengthFactor = options['Length factor']

        centralPath = options['Central path LUT']

        includeUreter = options['Include ureter']
        ostiumOptions = options['Ureter']
        ostiumDefaultOptions = ostiumOptions.getScaffoldSettings()
        elementsCountAroundOstium = ostiumDefaultOptions['Number of elements around ostium']
        elementsCountAnnulusRadially = options['Number of elements radially on annulus']
        ostiumPositionAround = options['Ostium position around']
        ostiumPositionDown = options['Ostium position down']

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        fm = region.getFieldmodule()
        fm.beginChange()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fm.findMeshByDimension(3)

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], '],')
        del tmpRegion

        # Find arcLength
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            length += arcLength
        segmentLength = length / elementsCountAlong

        if includeUrethra:
            urethraLength = lengthFactor * length
            elementsCountAlongUrethra = int(urethraLength / segmentLength)
            elementsCountAlongBladder = elementsCountAlong - elementsCountAlongUrethra
        else:
            elementsCountAlongBladder = elementsCountAlong

        # Sample central path
        # sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        # sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        if includeUrethra:
            sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
            sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        else:
            m = int((1 - lengthFactor) * len(cx))
            cx = cx[:m+1]
            cd1 = cd1[:m+1]
            cd2 = cd2[:m+1]
            cd12 = cd12[:m+1]
            sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongBladder)
            sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Create bladder
        # Create line down the bladder with major radius
        R1_max = 0.5 * majorDiameter
        R2_max = 0.5 * neckDiameter1
        # Create the body part
        if includeUrethra:
            a = (1 - lengthFactor) * length / 2
        else:
            a = length / 2
        b = R1_max
        xLoop, d2Loop = createEllipsePoints([0.0, 0.0, a], math.pi - neckAngleRadians, [0.0, 0.0, a], [0.0, b, 0.0],
                                            elementsCountAlongBladder, startRadians=math.pi)
        d1Loop = [0, 0, 0] * len(xLoop)
        nodesAlongMax_x = []
        nodesAlongMax_d1 = []
        nodesAlongMax_d2 = []
        for n in range(len(xLoop)):
            nodesAlongMax_x.append(xLoop[n])
            nodesAlongMax_d1.append(d1Loop[n])
            nodesAlongMax_d2.append(d2Loop[n])
        # Create the neck part
        px = [[nodesAlongMax_x[-1][0], nodesAlongMax_x[-1][1], nodesAlongMax_x[-1][2]], [0, (0 - R2_max), 2 * a]]
        pd2 = [nodesAlongMax_d2[-1], [0, 0, segmentLength]]
        pd2_smoothed = interp.smoothCubicHermiteDerivativesLine(px, pd2, fixAllDirections=False,
                                                                fixStartDerivative=True, fixEndDerivative=True,
                                                                fixStartDirection=False, fixEndDirection=False)
        pd1 = [0, 0, 0] * len(px)
        for n in range(1, len(pd2_smoothed)):
            nodesAlongMax_x.append(px[n])
            nodesAlongMax_d2.append(pd2_smoothed[n])
            nodesAlongMax_d1.append(pd1[n])
        nx_max, nd2_max, ne_max, nxi_max, nnf_max = interp.sampleCubicHermiteCurves(nodesAlongMax_x, nodesAlongMax_d2,
                                                                                    elementsCountAlongBladder)

        # Create line down the bladder with minor radius
        R1_min = 0.5 * minorDiameter
        R2_min = 0.5 * neckDiameter2
        # Create the body part
        if includeUrethra:
            a = (1 - lengthFactor) * length / 2
        else:
            a = length / 2
        b = R1_min
        xLoop, d2Loop = createEllipsePoints([0.0, 0.0, a], math.pi - neckAngleRadians, [0.0, 0.0, a], [0.0, b, 0.0],
                                            elementsCountAlongBladder, startRadians=math.pi)
        d1Loop = [0, 0, 0] * len(xLoop)
        nodesAlongMin_x = []
        nodesAlongMin_d1 = []
        nodesAlongMin_d2 = []
        for n in range(len(xLoop)):
            nodesAlongMin_x.append(xLoop[n])
            nodesAlongMin_d1.append(d1Loop[n])
            nodesAlongMin_d2.append(d2Loop[n])
        # Create the neck part
        px = [[nodesAlongMin_x[-1][0], nodesAlongMin_x[-1][1], nodesAlongMin_x[-1][2]], [0, (0 - R2_min), 2 * a]]
        pd2 = [nodesAlongMin_d2[-1], [0, 0, segmentLength]]
        pd2_smoothed = interp.smoothCubicHermiteDerivativesLine(px, pd2, fixAllDirections=False,
                                                                fixStartDerivative=True, fixEndDerivative=True,
                                                                fixStartDirection=False, fixEndDirection=False)
        pd1 = [0, 0, 0] * len(px)
        for n in range(1, len(pd2_smoothed)):
            nodesAlongMin_x.append(px[n])
            nodesAlongMin_d2.append(pd2_smoothed[n])
            nodesAlongMin_d1.append(pd1[n])
        nx_min, nd2_min, ne_min, nxi_min, nnf_min = interp.sampleCubicHermiteCurves(nodesAlongMin_x, nodesAlongMin_d2,
                                                                                    elementsCountAlongBladder)

        # Create the whole bladder nodes
        # Create ellipses surrounding the center line
        innerNodes_x = []
        innerNodes_d1 = []
        innerNodes_d2 = []
        nd1 = matrix.rotateAboutZAxis(nd2_max[0], 0.5 * math.pi)
        innerNodes_d1 += [nd1] * elementsCountAround
        for n1 in range(0, len(nx_max)):
            xAround, d1Around = createEllipsePoints([0.0, 0.0, nx_max[n1][2]], 2 * math.pi, [nx_min[n1][1], 0.0, 0.0],
                                                    [0.0, nx_max[n1][1], 0.0], elementsCountAround, startRadians=0.0)
            innerNodes_x += xAround
            if n1 >= 1:
                innerNodes_d1 += d1Around
        # Create lines from apex go down the bladder
        nodesInLineToDown_x = []
        nodesInLineToDown_d2 = []
        for n2 in range(elementsCountAround):
            x = [innerNodes_x[n1 * elementsCountAround + n2] for n1 in range(0, elementsCountAlongBladder + 1)]
            nodesInLineToDown_x += x
            for n1 in range(0, elementsCountAlongBladder + 1):
                if n1 == 0:
                    d2 = [segmentLength, 0.0, 0.0]
                elif n1 == elementsCountAlongBladder:
                    if not includeUrethra:
                        badderSegmentLength = (1 - lengthFactor) * length / elementsCountAlongBladder
                        d2 = [0.0, 0.0, badderSegmentLength]
                        # d2 = [0.0, 0.0, segmentLength]
                else:
                    v1 = nodesInLineToDown_x[(elementsCountAlongBladder + 1) * n2 + n1 - 1]
                    v2 = nodesInLineToDown_x[(elementsCountAlongBladder + 1) * n2 + n1]
                    d2 = [v2[c] - v1[c] for c in range(3)]
                nodesInLineToDown_d2.append(d2)
        # Smoothing the derivatives
        smoothed_d2 = []
        for n1 in range(elementsCountAround):
            lineSmoothingNodes = []
            lineSmoothingNodes_d2 = []
            for n2 in range(elementsCountAlongBladder + 1):
                lineSmoothingNodes.append(nodesInLineToDown_x[n1 * (elementsCountAlongBladder + 1) + n2])
                lineSmoothingNodes_d2.append(nodesInLineToDown_d2[n1 * (elementsCountAlongBladder + 1) + n2])
            smd2 = smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2, fixAllDirections=True,
                                                    fixStartDerivative=True, fixEndDerivative=True,
                                                    fixStartDirection=False, fixEndDirection=False)
            smoothed_d2 += smd2
        # Re-arrange the derivatives order
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround):
                rd2 = smoothed_d2[n1 * (elementsCountAlongBladder + 1) + n2]
                innerNodes_d2.append(rd2)

        # Create urethra
        if includeUrethra:
            last_z_bladder = innerNodes_x[-1][2]
            radiansPerElementAround = 2.0 * math.pi / elementsCountAround
            for n2 in range(elementsCountAlongUrethra + 1):
                for n1 in range(elementsCountAround):
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        -urethraDiameter1 * sinRadiansAround,
                        urethraDiameter2 * cosRadiansAround,
                        last_z_bladder + n2 * segmentLength
                    ]
                    dx_ds1 = [
                        -radiansPerElementAround * urethraDiameter1 * cosRadiansAround,
                        radiansPerElementAround * urethraDiameter2 * -sinRadiansAround,
                        0.0
                    ]
                    dx_ds2 = [0, 0, segmentLength]
                    if n2 == 0:
                        pass
                    else:
                        innerNodes_x.append(x)
                        innerNodes_d1.append(dx_ds1)
                        innerNodes_d2.append(dx_ds2)

        # Project reference point for warping onto central path
        sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
             tubemesh.getPlaneProjectionOnCentralPath(innerNodes_x, elementsCountAround, elementsCountAlong,
                                                      length, sx, sd1, sd2, sd12)

        innerRadiusAlong = []
        for n2 in range(elementsCountAlong + 1):
            firstNodeAlong = innerNodes_x[n2 * elementsCountAround]
            radius = vector.magnitude(firstNodeAlong[c] - [0.0, 0.0, firstNodeAlong[2]][c] for c in range(3))
            innerRadiusAlong.append(radius)

        # Warp points
        segmentAxis = [0.0, 0.0, 1.0]
        xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
            tubemesh.warpSegmentPoints(innerNodes_x, innerNodes_d1, innerNodes_d2, segmentAxis, sxRefList, sd1RefList,
                                       sd2ProjectedListRef, elementsCountAround, elementsCountAlong,
                                       zRefList, innerRadiusAlong, closedProximalEnd=True)

        if includeUrethra:
            wallThicknessList = [bladderWallThickness] * (elementsCountAlongBladder + 1) + [urethraWallThickness] * \
                                elementsCountAlongUrethra
        else:
            wallThicknessList = [bladderWallThickness] * (elementsCountAlong + 1)
        transitElementList = [0] * elementsCountAround

        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList,
                                                                                        d2WarpedList, d3WarpedUnitList,
                                                                                        wallThicknessList,
                                                                                        elementsCountAround,
                                                                                        elementsCountAlong,
                                                                                        elementsCountThroughWall,
                                                                                        transitElementList)

        # Deal with multiple nodes at end point for closed proximal end
        xApexInner = xList[0]
        # Arclength between apex point and corresponding point on next face
        mag = interp.getCubicHermiteArcLength(xList[0], d2List[0], xList[2 * elementsCountAround],
                                              d2List[2 * elementsCountAround])
        d2ApexInner = vector.setMagnitude(sd2[0], mag)
        d1ApexInner = vector.crossproduct3(sd1[0], d2ApexInner)
        d1ApexInner = vector.setMagnitude(d1ApexInner, mag)
        d3ApexUnit = vector.normalise(
            vector.crossproduct3(vector.normalise(d1ApexInner), vector.normalise(d2ApexInner)))
        d3ApexInner = [d3ApexUnit[c] * bladderWallThickness / elementsCountThroughWall for c in range(3)]

        xFinal = []
        d1Final = []
        d2Final = []
        d3Final = []
        for n3 in range(elementsCountThroughWall + 1):
            xApex = [xApexInner[c] + d3ApexUnit[c] * bladderWallThickness / elementsCountThroughWall * n3 for c in range(3)]
            xFinal.append(xApex)
            d1Final.append(d1ApexInner)
            d2Final.append(d2ApexInner)
            d3Final.append(d3ApexInner)

        xFinal += xList[(elementsCountThroughWall + 1) * elementsCountAround:]
        d1Final += d1List[(elementsCountThroughWall + 1) * elementsCountAround:]
        d2Final += d2List[(elementsCountThroughWall + 1) * elementsCountAround:]
        d3Final += d3List[(elementsCountThroughWall + 1) * elementsCountAround:]

        xFlat = d1Flat = d2Flat = d3Flat = []
        xTexture = d1Texture = d2Texture = d3Texture = []

        # Obtain elements count along body and neck of the bladder for defining annotation groups
        if includeUrethra:
            bladderLength = length - urethraLength
            bodyLength = ostiumPositionDown * bladderLength
        else:
            bodyLength = ostiumPositionDown * length
        elementsCountAlongBody = int(bodyLength / segmentLength)
        elementsCountAlongNeck = elementsCountAlongBladder - elementsCountAlongBody

        # Create annotation groups for bladder and urethra
        bodyGroup = AnnotationGroup(region, get_bladder_term("Dome of the Bladder"))
        neckGroup = AnnotationGroup(region, get_bladder_term("neck of urinary bladder"))
        bladderGroup = AnnotationGroup(region, get_bladder_term("urinary bladder"))
        urethraGroup = AnnotationGroup(region, get_bladder_term("urethra"))
        if includeUrethra:
            elementsCountAlongGroups = [elementsCountAlongBody, elementsCountAlongNeck, elementsCountAlongUrethra]
            annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup], [urethraGroup]]
        else:
            elementsCountAlongGroups = [elementsCountAlongBody, elementsCountAlongNeck]
            annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        annotationGroupsAround = []
        for i in range(elementsCountAround):
            annotationGroupsAround.append([])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([])

        neckMeshGroup = neckGroup. getMeshGroup(mesh)
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        urinaryBladderMeshGroup = bladderGroup.getMeshGroup(mesh)
        urethraMeshGroup = urethraGroup. getMeshGroup(mesh)

        # Create nodes and elements
        nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xFinal, d1Final, d2Final, d3Final, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier,
            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)

        # Define outer layer of the bladder to create the trackSurface on it
        outerNodes_x = []
        outerNodes_d1 = []
        outerNodes_d2 = []
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround):
                outerNodes_x.append(xList[(2 * n2 + 1) * elementsCountAround + n1])
                outerNodes_d1.append(d1List[(2 * n2 + 1) * elementsCountAround + n1])
                outerNodes_d2.append(d2List[(2 * n2 + 1) * elementsCountAround + n1])

        if includeUreter:
            elementsCount1 = elementsCountAround // 2
            elementsCount2 = elementsCountAlongBladder
            # Create trackSurface at the outer layer of the bladder for ostium 1
            nodesOnTrackSurface_x = []
            nodesOnTrackSurface_d1 = []
            nodesOnTrackSurface_d2 = []
            for n2 in range(elementsCountAlongBladder + 1):
                for n1 in range(elementsCountAround // 2 + 1):
                    nodesOnTrackSurface_x.append(outerNodes_x[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface_d1.append(outerNodes_d1[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface_d2.append(outerNodes_d2[n2 * elementsCountAround + n1])
            trackSurfaceOstium1 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface_x,
                                               nodesOnTrackSurface_d1, nodesOnTrackSurface_d2)
            ostium1Position = trackSurfaceOstium1.createPositionProportion(ostiumPositionAround, ostiumPositionDown)
            ostium1Position.xi1 = 1.0
            ostium1Position.xi2 = 1.0
            ostiumElementPositionAround = ostium1Position.e1
            ostiumElementPositionDown = ostium1Position.e2
            # Create trackSurface at the outer layer of the bladder for ostium 2
            nodesOnTrackSurface2_x = []
            nodesOnTrackSurface2_d1 = []
            nodesOnTrackSurface2_d2 = []
            for n2 in range(elementsCountAlongBladder + 1):
                for n1 in range(elementsCountAround // 2, elementsCountAround):
                    nodesOnTrackSurface2_x.append(outerNodes_x[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface2_d1.append(outerNodes_d1[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface2_d2.append(outerNodes_d2[n2 * elementsCountAround + n1])
                nodesOnTrackSurface2_x.append(outerNodes_x[n2 * elementsCountAround])
                nodesOnTrackSurface2_d1.append(outerNodes_d1[n2 * elementsCountAround])
                nodesOnTrackSurface2_d2.append(outerNodes_d2[n2 * elementsCountAround])
            trackSurfaceOstium2 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface2_x,
                                               nodesOnTrackSurface2_d1, nodesOnTrackSurface2_d2)
            ostium2Position = TrackSurfacePosition(elementsCountAround - ostiumElementPositionAround,
                                                   ostiumElementPositionDown - 1, 0.0, 1.0)
            annulusMeshGroups = [neckMeshGroup, urinaryBladderMeshGroup]
            generateOstiumsAndAnnulusMeshOnBladder(region, fm, nodes, mesh, ostiumDefaultOptions,
                                                  elementsCountAround, elementsCountAroundOstium, trackSurfaceOstium1,
                                                  ostium1Position, trackSurfaceOstium2, ostium2Position,
                                                  ostiumElementPositionDown, ostiumElementPositionAround,
                                                  xFinal, d1Final, d2Final, nextNodeIdentifier,
                                                  nextElementIdentifier, elementsCountAnnulusRadially,
                                                  annulusMeshGroups)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        '''
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        '''
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        '''
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        '''
        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("Dome of the Bladder"))
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("neck of urinary bladder"))
        urinaryBladderGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("urinary bladder"))

        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))

        is_body = bodyGroup.getFieldElementGroup(mesh2d)
        is_body_serosa = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_lumen = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_neck = neckGroup.getFieldElementGroup(mesh2d)
        is_neck_serosa = fm.createFieldAnd(is_neck, is_exterior_face_xi3_1)
        is_neck_lumen = fm.createFieldAnd(is_neck, is_exterior_face_xi3_0)

        serosaOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("serosa of body of urinary bladder"))
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_serosa)
        lumenOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("lumen of body of urinary bladder"))
        lumenOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_lumen)

        serosaOfNeck = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("serosa of neck of urinary bladder"))
        serosaOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_serosa)
        lumenOfNeck = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("lumen of neck of urinary bladder"))
        lumenOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_lumen)

        if options['Include urethra'] == True:
            urethraGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("urethra"))

            is_urethra = urethraGroup.getFieldElementGroup(mesh2d)
            is_urethra_serosa = fm.createFieldAnd(is_urethra, is_exterior_face_xi3_1)
            is_urethra_lumen = fm.createFieldAnd(is_urethra, is_exterior_face_xi3_0)

            serosaOfUrethra = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("serosa of urethra"))
            serosaOfUrethra.getMeshGroup(mesh2d).addElementsConditional(is_urethra_serosa)
            lumenOfUrethra = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("lumen of urethra"))
            lumenOfUrethra.getMeshGroup(mesh2d).addElementsConditional(is_urethra_lumen)

def generateOstiumsAndAnnulusMeshOnBladder(region, fm, nodes, mesh, ostiumDefaultOptions,elementsCountAround,
                                           elementsCountAroundOstium, trackSurfaceOstium1, ostium1Position,
                                           trackSurfaceOstium2, ostium2Position, ostiumElementPositionDown,
                                           ostiumElementPositionAround, xBladder, d1Bladder, d2Bladder,
                                           nextNodeIdentifier, nextElementIdentifier, elementsCountAnnulusRadially,
                                           annulusMeshGroups = []):

    # Create ureters on the surface
    # Ureter 1
    centerUreter1_x, centerUreter1_d1, centerUreter1_d2 = trackSurfaceOstium1.evaluateCoordinates(ostium1Position,
                                                                                          derivatives=True)
    td1, td2, td3 = calculate_surface_axes(centerUreter1_d1, centerUreter1_d2, [1.0, 0.0, 0.0])
    m1 = 2 * elementsCountAround * (ostiumElementPositionDown - 1) + ostiumElementPositionAround + 2
    ureter1StartCornerx = xBladder[m1]
    v1 = [(ureter1StartCornerx[c] - centerUreter1_x[c]) for c in range(3)]
    ostium1Direction = vector.crossproduct3(td3, v1)
    nodeIdentifier, elementIdentifier, (o1_x, o1_d1, o1_d2, _, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, ostiumDefaultOptions, trackSurfaceOstium1, ostium1Position, ostium1Direction,
                           startNodeIdentifier=nextNodeIdentifier, startElementIdentifier=nextElementIdentifier,
                           ostiumMeshGroups=None)
    # Ureter 2
    centerUreter2_x, centerUreter2_d1, centerUreter2_d2 = trackSurfaceOstium2.evaluateCoordinates(ostium2Position,
                                                                                                  derivatives=True)
    ad1, ad2, ad3 = calculate_surface_axes(centerUreter2_d1, centerUreter2_d2, [1.0, 0.0, 0.0])
    m2 = 2 * elementsCountAround * (ostiumElementPositionDown - 1) + elementsCountAround - ostiumElementPositionAround
    ureter2StartCornerx = xBladder[m2]
    v2 = [(ureter2StartCornerx[c] - centerUreter2_x[c]) for c in range(3)]
    ostium2Direction = vector.crossproduct3(ad3, v2)
    nodeIdentifier, elementIdentifier, (o2_x, o2_d1, o2_d2, _, o2_NodeId, o2_Positions) = \
        generateOstiumMesh(region, ostiumDefaultOptions, trackSurfaceOstium2, ostium2Position, ostium2Direction,
                           startNodeIdentifier=nodeIdentifier, startElementIdentifier=elementIdentifier)

    # Create annulus mesh around ostiums
    endPoints1_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endPoints1_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endPoints1_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endNode1_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endDerivativesMap = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endPoints2_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endPoints2_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endPoints2_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
    endNode2_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]

    for n3 in range(2):
        n1 = 0
        endNode1_Id[n3][n1] = (2 * ostiumElementPositionDown - 2 + n3) * elementsCountAround + \
                              ostiumElementPositionAround + 1 + 2
        endNode1_Id[n3][n1 + 1] = endNode1_Id[n3][n1] + 2 * elementsCountAround
        endNode1_Id[n3][n1 + 2] = endNode1_Id[n3][n1 + 1] + 2 * elementsCountAround
        endNode1_Id[n3][n1 + 3] = endNode1_Id[n3][n1 + 2] + 1
        endNode1_Id[n3][n1 + 4] = endNode1_Id[n3][n1 + 3] + 1
        endNode1_Id[n3][n1 + 5] = endNode1_Id[n3][n1 + 1] + 2
        endNode1_Id[n3][n1 + 6] = endNode1_Id[n3][n1] + 2
        endNode1_Id[n3][n1 + 7] = endNode1_Id[n3][n1] + 1

        endNode2_Id[n3][n1] = (2 * ostiumElementPositionDown - 2 + n3) * elementsCountAround + elementsCountAround - \
                              ostiumElementPositionAround - 1 + 2
        endNode2_Id[n3][n1 + 1] = endNode2_Id[n3][n1] + 2 * elementsCountAround
        endNode2_Id[n3][n1 + 2] = endNode2_Id[n3][n1 + 1] + 2 * elementsCountAround
        endNode2_Id[n3][n1 + 3] = endNode2_Id[n3][n1 + 2] + 1
        endNode2_Id[n3][n1 + 7] = endNode2_Id[n3][n1] + 1
        if ostiumElementPositionAround == 0:
            endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] - elementsCountAround + 1
            endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 4] - 2 * elementsCountAround
            endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1 + 5] - 2 * elementsCountAround
        else:
            endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] + 1
            endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 1] + 2
            endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1] + 2

    for n3 in range(2):
        for n1 in range(elementsCountAroundOstium):
            nc1 = endNode1_Id[n3][n1] - 1
            endPoints1_x[n3][n1] = xBladder[nc1]
            endPoints1_d1[n3][n1] = d1Bladder[nc1]
            endPoints1_d2[n3][n1] = [d2Bladder[nc1][c] for c in range(3)]
            nc2 = endNode2_Id[n3][n1] - 1
            endPoints2_x[n3][n1] = xBladder[nc2]
            endPoints2_d1[n3][n1] = d1Bladder[nc2]
            endPoints2_d2[n3][n1] = d2Bladder[nc2]

    for n1 in range(elementsCountAroundOstium):
        if n1 == 0:
            endDerivativesMap[0][n1] = ((-1, 0, 0), (-1, -1, 0), None, (0, 1, 0))
            endDerivativesMap[1][n1] = ((-1, 0, 0), (-1, -1, 0), None, (0, 1, 0))
        elif n1 == 1:
            endDerivativesMap[0][n1] = ((0, 1, 0), (-1, 0, 0), None)
            endDerivativesMap[1][n1] = ((0, 1, 0), (-1, 0, 0), None)
        elif n1 == 2:
            endDerivativesMap[0][n1] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
            endDerivativesMap[1][n1] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
        elif n1 == 3:
            endDerivativesMap[0][n1] = ((1, 0, 0), (0, 1, 0), None)
            endDerivativesMap[1][n1] = ((1, 0, 0), (0, 1, 0), None)
        elif n1 == 4:
            endDerivativesMap[0][n1] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
            endDerivativesMap[1][n1] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
        elif n1 == 5:
            endDerivativesMap[0][n1] = ((0, -1, 0), (1, 0, 0), None)
            endDerivativesMap[1][n1] = ((0, -1, 0), (1, 0, 0), None)
        elif n1 == 6:
            endDerivativesMap[0][n1] = ((0, -1, 0), (1, -1, 0), None, (-1, 0, 0))
            endDerivativesMap[1][n1] = ((0, -1, 0), (1, -1, 0), None, (-1, 0, 0))
        else:
            endDerivativesMap[0][n1] = ((-1, 0, 0), (0, -1, 0), None)
            endDerivativesMap[1][n1] = ((-1, 0, 0), (0, -1, 0), None)

    nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nodeIdentifier, elementIdentifier,
        o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
        endPoints1_x, endPoints1_d1, endPoints1_d2, None, endNode1_Id, endDerivativesMap,
        elementsCountRadial=elementsCountAnnulusRadially, meshGroups=annulusMeshGroups)

    nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nodeIdentifier, elementIdentifier,
        o2_x, o2_d1, o2_d2, None, o2_NodeId, None,
        endPoints2_x, endPoints2_d1, endPoints2_d2, None, endNode2_Id, endDerivativesMap,
        elementsCountRadial=elementsCountAnnulusRadially, meshGroups=annulusMeshGroups)

    # Store elements to be deleted later from bladder mesh
    element_identifiers = []
    for n3 in range(2):
        elementIdxUnderOstium1 = ostiumElementPositionDown * elementsCountAround + ostiumElementPositionAround + \
                                 n3 * elementsCountAround + 1
        element_identifiers.append(elementIdxUnderOstium1)
        element_identifiers.append(elementIdxUnderOstium1 + 1)
        elementIdxUnderOstium2 = ostiumElementPositionDown * elementsCountAround + elementsCountAround - \
                                 ostiumElementPositionAround + n3 * elementsCountAround
        element_identifiers.append(elementIdxUnderOstium2)
        element_identifiers.append(elementIdxUnderOstium2 - 1)

        # Delete elements under annulus mesh
        mesh_destroy_elements_and_nodes_by_identifiers(mesh, element_identifiers)
    return
