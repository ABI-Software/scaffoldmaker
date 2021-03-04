"""
Generates a 3-D stomach mesh along the central line, with variable
numbers of elements around oesophagus and duodenum, along and through
wall, with variable radius and thickness along.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftNodeValueLabel
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector


class MeshType_3d_stomach1(Scaffold_base):
    """
    Generates a 3-D stomach mesh with variable numbers of elements around the oesophagus and duodenum,
    along the central line, and through wall. The stomach is created by a function that generates a bean
    volume defined by a central path as its longitudinal axis. D2 of the central path points to the greater
    curvature of the stomach and magnitude of D2 and D3 are the radii of the stomach in the respective
    direction.
    """
    centralPathDefaultScaffoldPackages = {
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 5
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[17.0, 14.0, 0.0], [-12.4, -22.4, 0.0], [7.9, -4.3, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, 0.0]],
                    [[10.0, 0.0, 0.0], [-19.8, -2.6, 0.0], [1.1, -8.9, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 0.0], [-19.8, -2.6, 0.0], [1.1, -8.9, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, 0.0]],
                    [[-10.0, 7.0, 0.0], [-8.1, 9.4, 0.0], [-1.9, -1.6, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.5],
                     [0.0, 0.0, 0.0]],
                    [[-15.0, 7.0, 0.0], [-8.1, 9.4, 0.0], [-1.9, -1.6, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.5],
                     [0.0, 0.0, 0.0]],
                    [[-17.0, 14.0, 0.0], [0.2, 2.6, 0.0], [-3.5, 0.3, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 3.5],
                     [0.0, 0.0, 0.0]]])
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[17.0, 14.0, 0.0], [-12.4, -22.4, 0.0], [7.9, -4.3, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 0.0], [-19.8, -2.6, 0.0], [1.1, -8.9, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, 0.0]],
                    [[-15.0, 7.0, 0.0], [-6.2, 7.0, 0.0], [-2.0, -1.4, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.5],
                     [0.0, 0.0, 0.0]],
                    [[-19.2, 13.9, 0.0], [-4.1, 9.0, 0.0], [-3.2, -1.4, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 3.5],
                     [0.0, 0.0, 0.0]]])
        }),
    }

    ostiumDefaultScaffoldPackages = {
        'Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 5.0,
                'Ostium length': 5.0,
                'Ostium wall thickness': 0.5,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 1.25,
                'Vessel wall thickness': 0.5,
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
        'Human 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 5.0,
                'Ostium length': 5.0,
                'Ostium wall thickness': 0.5,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 1.25,
                'Vessel wall thickness': 0.5,
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
    }

    @staticmethod
    def getName():
        return '3D Stomach 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Rat 1',
            'Human 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around oesophagus': 8,
            'Number of elements around duodenum': 12,
            'Number of elements along': 8,
            'Number of elements through wall': 1,
            'Wall thickness': 0.5,
            'Gastro-oesophagal junction': copy.deepcopy(ostiumOption),
            'Gastro-oesophagal junction position along factor': 0.3,
            'Use cross derivatives': False,
            'Use linear through wall' : False, # need to deal with wedge not available in bicubichermite
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        # Add default options for Rat and human later

        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around oesophagus',
            'Number of elements around duodenum',
            'Number of elements along',
            'Number of elements through wall',
            'Wall thickness',
            'Gastro-oesophagal junction',
            'Gastro-oesophagal junction position along factor',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        if optionName == 'Gastro-oesophagal junction':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        if optionName == 'Gastro-oesophagal junction':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
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
        if optionName == 'Gastro-oesophagal junction':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Gastro-oesophagal junction'].getScaffoldType() in cls.getOptionValidScaffoldTypes(
                'Gastro-oesophagal junction'):
            options['Gastro-oesophagal junction'] = cls.getOptionScaffoldPackage('Gastro-oesophagal junction',
                                                                                 MeshType_3d_ostium1)
        if options['Number of elements around oesophagus'] < 8:
            options['Number of elements around oesophagus'] = 8
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        for key in [
            'Number of elements around oesophagus',
            'Number of elements around duodenum']:
            if options[key] % 2:
                options[key] += 1
        if options['Number of elements along'] < 8:
            options['Number of elements along'] = 8
        if options['Number of elements through wall'] < 1:
            options['Number of elements through wall'] = 1
        for key in [
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1

        cls.updateSubScaffoldOptions(options)

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        '''
        Update ostium sub-scaffold options which depend on parent options.
        '''
        wallThickness = options['Wall thickness']
        ostiumOptions = options['Gastro-oesophagal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        ostiumSettings['Ostium wall thickness'] = wallThickness

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        cls.updateSubScaffoldOptions(options)
        centralPath = options['Central path']
        elementsCountAroundOesophagus = options['Number of elements around oesophagus']
        elementsCountAroundDuodenum = options['Number of elements around duodenum']
        elementsCountAlong = options['Number of elements along']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not (options['Use linear through wall'])

        GOJPositionAlongFactor = options['Gastro-oesophagal junction position along factor']
        GOJOptions = options['Gastro-oesophagal junction']
        GOJSettings = GOJOptions.getScaffoldSettings()
        oesophagusDiameter = GOJSettings['Ostium diameter']

        elementsCountAlongTrackSurface = 20

        ############################################################################################
        zero = [0.0, 0.0, 0.0]

        fm = region.getFieldmodule()
        fm.beginChange()

        coordinates = findOrCreateFieldCoordinates(fm)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if useCubicHermiteThroughWall:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            if useCrossDerivatives:
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        cache = fm.createFieldcache()
        ###########################################################################################

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd3, cd12, cd13 = extractPathParametersFromRegion(tmpRegion,
                                                                        [Node.VALUE_LABEL_VALUE,
                                                                         Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                                                         Node.VALUE_LABEL_D_DS3,
                                                                         Node.VALUE_LABEL_D2_DS1DS2,
                                                                         Node.VALUE_LABEL_D2_DS1DS3])
        del tmpRegion
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], ',', cd3[i], ',', cd13[i], '],')

        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongTrackSurface)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        # Calculate length of central path
        stomachCentralPathLength = 0.0
        for e in range(len(sx) - 1):
            arcLength = interp.getCubicHermiteArcLength(sx[e], sd1[e], sx[e + 1], sd1[e + 1])
            # print(e+1, arcLength)
            stomachCentralPathLength += arcLength
        lengthElementAlong = stomachCentralPathLength / elementsCountAlongTrackSurface

        nodeIdentifier = firstNodeIdentifier
        elementIdentifier = firstElementIdentifier

        # Fundus diameter
        fundusRadius = vector.magnitude(sd2[0])
        elementsAlongFundus = int(fundusRadius / lengthElementAlong)

        d2Apex = []
        d2 = sd2[0]
        for n1 in range(elementsCountAroundDuodenum):
            rotAngle = n1 * 2.0 * math.pi / elementsCountAroundDuodenum
            rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(sd1[0]), rotAngle)
            d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            d2Apex.append(d2Rot)

        xEllipses = []
        d1Ellipses = []
        for n in range(elementsAlongFundus + 1, len(sx)):
            px, pd1 = createEllipsePoints(sx[n], 2 * math.pi, sd2[n], sd3[n], elementsCountAroundDuodenum,
                                          startRadians=0.0)
            xEllipses.append(px)
            d1Ellipses.append(pd1)

        # Find d2
        d2Raw = []
        for n1 in range(elementsCountAroundDuodenum):
            xAlong = []
            d2Along = []
            for n2 in range(len(xEllipses) - 1):
                v1 = xEllipses[n2][n1]
                v2 = xEllipses[n2 + 1][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xAlong.append(v1)
                d2Along.append(d2)
            xAlong.append(xEllipses[-1][n1])
            d2Along.append(d2)
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
            d2Raw.append(d2Smoothed)

        # Rearrange d2
        d2Ellipses = []
        for n2 in range(len(xEllipses)):
            d2Around = []
            for n1 in range(elementsCountAroundDuodenum):
                d2 = d2Raw[n1][n2]
                d2Around.append(d2)
            d2Ellipses.append(d2Around)

        # for n2 in range(len(xEllipses)):
        #     for n1 in range(elementsCountAroundDuodenum):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xEllipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Ellipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Ellipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Merge fundus and body
        xAll = [[sx[0]] * elementsCountAroundDuodenum] + xEllipses
        d2All = [d2Apex] + d2Ellipses

        # for n2 in range(len(xAll)):
        #     for n1 in range(elementsCountAroundDuodenum):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAll[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2All[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Spread out elements
        xRaw = []
        d2Raw = []
        for n1 in range(elementsCountAroundDuodenum):
            xAlong = []
            d2Along = []
            for n2 in range(len(xAll)):
                xAlong.append(xAll[n2][n1])
                d2Along.append(d2All[n2][n1])
            xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                            elementsCountAlongTrackSurface,
                                                                            arcLengthDerivatives=True)[0:2]
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
            xRaw.append(xSampledAlong)
            d2Raw.append(d2Smoothed)

        # Rearrange x and d2
        xSampledAll = []
        d1SampledAll = []
        d2SampledAll = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            xAround = []
            d1Around = []
            d2Around = []
            for n1 in range(elementsCountAroundDuodenum):
                x = xRaw[n1][n2]
                d2 = d2Raw[n1][n2]
                xAround.append(x)
                d2Around.append(d2)

                # Calculate d1
                if n2 > 0:
                    v1 = xRaw[n1][n2]
                    v2 = xRaw[n1 + 1 if n1 < elementsCountAroundDuodenum - 2 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1Around.append(d1)
                else:
                    d1Around.append(d2Raw[int(elementsCountAroundDuodenum * 0.75)][0])

            if n2 > 0:
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
            else:
                d1Smoothed = d1Around

            xSampledAll.append(xAround)
            d1SampledAll.append(d1Smoothed)
            d2SampledAll.append(d2Around)

        # for n2 in range(elementsCountAlongTrackSurface + 1):
        #     for n1 in range(elementsCountAroundDuodenum):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAll[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAll[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAll[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Create tracksurface
        xTrackSurface = []
        d1TrackSurface = []
        d2TrackSurface = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            for n1 in range(elementsCountAroundDuodenum):
                xTrackSurface.append(xSampledAll[n2][n1])
                d1TrackSurface.append(d1SampledAll[n2][n1])
                d2TrackSurface.append(d2SampledAll[n2][n1])

        trackSurfaceStomach = TrackSurface(elementsCountAroundDuodenum, elementsCountAlongTrackSurface,
                                           xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

        # for n2 in range(len(xTrackSurface)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1TrackSurface[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2TrackSurface[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Set up gastro-oesophagal junction
        GOJSettings['Number of elements around ostium'] = elementsCountAroundOesophagus
        GOJPosition = trackSurfaceStomach.createPositionProportion(0.5, GOJPositionAlongFactor)
        xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(GOJPosition, derivatives=True)
        axis1 = d1Centre

        # fm = region.getFieldmodule()
        mesh = fm.findMeshByDimension(3)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        nextNodeIdentifier = nodeIdentifier
        nextElementIdentifier = elementIdentifier
        nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
            generateOstiumMesh(region, GOJSettings, trackSurfaceStomach, GOJPosition, axis1,
                               nextNodeIdentifier, nextElementIdentifier)
        bodyStartNode = nextNodeIdentifier

        # From oesophagus to duodenum along lesser curvature (LC)
        elementsOesoToDuodLC = elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) # CHECK
        startProportion1, startProportion2 = trackSurfaceStomach.getProportion(
            o1_Positions[int(elementsCountAroundOesophagus * 0.5)])
        d1Start = o1_d2[1][int(elementsCountAroundOesophagus * 0.5)]

        xOesoToDuodLC, d2OesoToDuodLC = getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1,
                                                                               startProportion2, 0.5, 1.0,
                                                                               elementsOesoToDuodLC)

        # From oesophagus to duodenum along greater curvature (GC)
        xAlongGC = []
        d2AlongGC = []
        xAlongGC.append(o1_x[1][0])
        d2AlongGC.append(o1_d2[1][0])

        oesoStartProportion2 = trackSurfaceStomach.getProportion(o1_Positions[0])[1]
        elementsAlongUpstreamOfOeso = int(elementsCountAlongTrackSurface * oesoStartProportion2)

        # From oesophagus to fundus apex
        arcLengthOesoApex = 0.0
        for n2 in range(elementsAlongUpstreamOfOeso):
            nAlong = elementsAlongUpstreamOfOeso - n2
            v1 = xSampledAll[nAlong][int(elementsCountAroundDuodenum * 0.5)]
            v2 = xSampledAll[nAlong - 1][int(elementsCountAroundDuodenum * 0.5)]
            d = [v2[c] - v1[c] for c in range(3)]
            arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
            arcLengthOesoApex += arcLengthAround
            d2 = [c * arcLengthAround for c in vector.normalise(d)]
            xAlongGC.append(v1)
            d2AlongGC.append(d2)

        # From fundus apex to duodenum
        for n2 in range(len(xSampledAll)):
            xAlongGC.append(xSampledAll[n2][0])
            d2AlongGC.append(d2SampledAll[n2][0])

        elementsOesoToDuodGC = elementsCountAlong + int(elementsCountAroundDuodenum * 0.5) - 2 # Check, old - elementsCountAlong + int(elementsCountAroundOesophagus * 0.5)
        xOesoToDuodGC, d2OesoToDuodGC = interp.sampleCubicHermiteCurvesSmooth(xAlongGC, d2AlongGC, elementsOesoToDuodGC,
                                                                              derivativeMagnitudeStart=0.5 * vector.magnitude(
                                                                                  d2AlongGC[0]))[0:2]

        arcLength = 0.0
        for e in range(len(xOesoToDuodGC) - 1):
            arcLength += interp.getCubicHermiteArcLength(xOesoToDuodGC[e], d2OesoToDuodGC[e],
                                                         xOesoToDuodGC[e + 1], d2OesoToDuodGC[e + 1])
            if arcLength > arcLengthOesoApex:
                nodesCountFromOesoToApex = e + 2
                break

        d2OesoToDuodGC = interp.smoothCubicHermiteDerivativesLine(xOesoToDuodGC, d2OesoToDuodGC)

        nodeIdentifier = nextNodeIdentifier

        # for n2 in range(len(xOesoToDuodLC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOesoToDuodLC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2OesoToDuodLC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # for n2 in range(len(xOesoToDuodGC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOesoToDuodGC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2OesoToDuodGC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Spread out elements around for each egroup around
        # Elements around oesophagus
        # Decide where to put extra elements if elementsOesophagus/2 is odd number.
        # Right now, extra elements go downstream of bifurcation
        # For even numbers, we split elements at bifurcation

        ptsOnTrackSurfaceGC = []
        xAlongAround = []
        d1AlongAround = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceGC.append(xSampledAll[n2][0])

        for idx in range(math.ceil(elementsCountAroundOesophagus * 0.25)):
            # First half
            ostiumIdx = int(elementsCountAroundOesophagus * 0.25) + idx
            GCIdx = -(elementsCountAlong - int(elementsCountAroundOesophagus * 0.25)) + idx # 4 - elementsCountAlong - math.ceil(elementsCountAroundOesophagus * 0.25) + idx

            # Find closest position to point on GC as track surface is not a simple geometry near dome
            GCPosition, d1GC = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                trackSurfaceStomach, 0.0, elementsCountAlongTrackSurface)
            GCProportion1, GCProportion2 = trackSurfaceStomach.getProportion(GCPosition)

            endPosition = o1_Positions[ostiumIdx]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(o1_d1[1][ostiumIdx]), math.pi)
            d2 = o1_d2[1][ostiumIdx]
            d1EndOstium = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
            d1EndTrackSurface = trackSurfaceStomach.evaluateCoordinates(endPosition, derivatives=True)[1]
            # remove d1Ave when cleaning up code
            xi = 0 #ostiumIdx / int(elementsCountAroundOesophagus * 0.5)
            d1Ave = [d1EndOstium[c] * (1 - xi) + d1EndTrackSurface[c] * xi for c in range(3)]

            xFirstHalf, d1FirstHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, 0.0, GCProportion2, endProportion1,
                                                       endProportion2, int(0.5 * elementsCountAroundDuodenum + 1),
                                                       startDerivative=d1GC, endDerivative=d1Ave,
                                                       endDerivativeMagnitude= 0.5 * vector.magnitude(d1Ave))

            # Second half
            ostiumIdx2 = elementsCountAroundOesophagus - (int(elementsCountAroundOesophagus * 0.25) + idx)
            startPosition = o1_Positions[ostiumIdx2]
            d1StartOstium = o1_d2[1][ostiumIdx2]
            startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
            d1StartTrackSurface = trackSurfaceStomach.evaluateCoordinates(startPosition, derivatives=True)[1]
            # Clean up
            d1Ave = [d1StartOstium[c] * (1 - xi) + d1StartTrackSurface[c] * xi for c in range(3)]

            xSecondHalf, d1SecondHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1, startProportion2, 1.0,
                                                       GCProportion2, int(0.5 * elementsCountAroundDuodenum + 1),
                                                       startDerivative=d1Ave, endDerivative=d1GC,
                                                       startDerivativeMagnitude= 0.5 * vector.magnitude(d1Ave))

            xAround = xFirstHalf[:-1] + xSecondHalf[1:-1]
            d1Around = d1FirstHalf[:-1] + d1SecondHalf[1:-1]

            xAlongAround.append(xAround)
            d1AlongAround.append(d1Around)

        # Elements downstream of oesophagus
        for idx in range(-(elementsCountAlong - int(elementsCountAroundOesophagus * 0.5)), 0):  #4 - elementsCountAlong, 0):
            startPosition, d1Start = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[idx],
                                                                                    ptsOnTrackSurfaceGC,
                                                                                    trackSurfaceStomach, 0.0,
                                                                                    elementsCountAlongTrackSurface)
            startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)

            endPosition = trackSurfaceStomach.findNearestPosition(xOesoToDuodLC[idx])
            endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
            d1End = trackSurfaceStomach.evaluateCoordinates(endPosition, derivatives=True)[1]

            xFirstHalf, d1FirstHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, 0.0, startProportion2, endProportion1,
                                                       endProportion2, int(0.5 * elementsCountAroundDuodenum),
                                                       startDerivative=d1Start, endDerivative=d1End)

            xSecondHalf, d1SecondHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, endProportion1, endProportion2, 1.0,
                                                       startProportion2, int(0.5 * elementsCountAroundDuodenum),
                                                       startDerivative=d1End, endDerivative=d1Start)

            xAround = xFirstHalf + xSecondHalf[1:-1]
            d1Around = d1FirstHalf + d1SecondHalf[1:-1]
            d1Around = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)

            xAlongAround.append(xAround)
            d1AlongAround.append(d1Around)

        # for n2 in range(len(xAlongAround)):
        #     for n1 in range(len(xAlongAround[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongAround[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1AlongAround[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # # Joining loops at fundus
        ptsOnTrackSurfaceOesoToFundus = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceOesoToFundus.append(xSampledAll[n2][int(elementsCountAroundDuodenum * 0.5)])

        xLoopsRight = []
        d2LoopsRight = []
        xLoopsLeft = []
        d2LoopsLeft = []

        for n in range(int(elementsCountAroundDuodenum * 0.5 - 1)):
            GCIdx = n + 1
            if GCIdx < nodesCountFromOesoToApex:
                ptsOnTrackSurface = ptsOnTrackSurfaceOesoToFundus
                proportion1 = 0.5
            else:
                ptsOnTrackSurface = ptsOnTrackSurfaceGC
                proportion1 = 0.0
            d2GC = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[GCIdx],
                                                                  ptsOnTrackSurface,
                                                                  trackSurfaceStomach, proportion1,
                                                                  elementsCountAlongTrackSurface)[1]

            if GCIdx < nodesCountFromOesoToApex:
                rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d2OesoToDuodGC[GCIdx]), math.pi)
                d2GCRot = [rotFrame[j][0] * d2GC[0] + rotFrame[j][1] * d2GC[1] + rotFrame[j][2] * d2GC[2] for j in range(3)]
                d2GC = d2GCRot

            d2End = [xAlongAround[1][int(elementsCountAroundDuodenum * 0.5) - n][c] -
                     xAlongAround[0][int(elementsCountAroundDuodenum * 0.5) - n][c] for c in range(3)]

            nx = [ xOesoToDuodGC[GCIdx], xAlongAround[0][int(elementsCountAroundDuodenum * 0.5) - n]]
            nd2 = [d2GC, d2End]
            xRight, d2Right = interp.sampleCubicHermiteCurves(nx, nd2, int(elementsCountAroundOesophagus * 0.25) + (1 if n > 0 else 0), arcLengthDerivatives=True)[0:2]

            # Calculate and append d2 for elements downstream of loop
            if n == 0:
                xAnnulusRight = xRight
                d2AnnulusRight = d2Right
                x1 = xAlongAround[1][int(elementsCountAroundDuodenum * 0.5)]
                x2 = xOesoToDuodLC[1]
                d2 = findDerivativeBetweenPoints(x1, x2)
                xAnnulusRight += [x1, x2]
                d2AnnulusRight += [d2, d2OesoToDuodLC[1]]
                d2AnnulusRight = interp.smoothCubicHermiteDerivativesLine(xAnnulusRight, d2AnnulusRight,
                                                                          fixStartDirection=True,
                                                                          fixEndDerivative= True)
                # for n2 in range(len(xAnnulusRight)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusRight[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2AnnulusRight[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                #     nodeIdentifier += 1
            else:
                for n2 in range(1, len(xAlongAround)):
                    x1 = xAlongAround[n2][int(elementsCountAroundDuodenum * 0.5) - n]
                    x2 = xAlongAround[n2 - 1][int(elementsCountAroundDuodenum * 0.5) - n]
                    d2 = findDerivativeBetweenPoints(x2, x1)
                    xRight.append(x1)
                    d2Right.append(d2)

                # for n2 in range(len(xRight)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xRight[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2Right[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                #     nodeIdentifier += 1

                d2Right = interp.smoothCubicHermiteDerivativesLine(xRight, d2Right, fixStartDirection=True)

                xLoopsRight.append(xRight)
                d2LoopsRight.append(d2Right)

            # Repeat for left side
            rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d2OesoToDuodGC[GCIdx]), math.pi)
            d2GCRot = [rotFrame[j][0] * d2GC[0] + rotFrame[j][1] * d2GC[1] + rotFrame[j][2] * d2GC[2] for j in range(3)]
            d2GC = d2GCRot

            xEnd = xAlongAround[0][int(elementsCountAroundDuodenum * 0.5 + 1) + n]
            d2End = [xAlongAround[1][int(elementsCountAroundDuodenum * 0.5 + 1) + n][c] -
                     xAlongAround[0][int(elementsCountAroundDuodenum * 0.5 + 1) + n][c] for c in range(3)]

            nx = [xOesoToDuodGC[GCIdx], xEnd]
            nd2 = [d2GC, d2End]
            xLeft, d2Left = interp.sampleCubicHermiteCurves(nx, nd2, int(elementsCountAroundOesophagus * 0.25) + (1 if n > 0 else 0),
                                                              arcLengthDerivatives=True)[0:2]

            # Calculate and append d2 for elements downstream of loop
            if n == 0:
                # Deal with addition of elementsAlong later
                xAnnulusLeft = xLeft
                d2AnnulusLeft = d2Left
                x1 = xAlongAround[1][int(elementsCountAroundDuodenum * 0.5) + 1]
                x2 = xOesoToDuodLC[1]
                d2 = findDerivativeBetweenPoints(x1, x2)
                xAnnulusLeft += [x1, x2]
                d2AnnulusLeft += [d2, d2OesoToDuodLC[1]]
                d2AnnulusLeft = interp.smoothCubicHermiteDerivativesLine(xAnnulusLeft, d2AnnulusLeft,
                                                                          fixStartDirection=True,
                                                                          fixEndDerivative=True)

                # for n2 in range(len(xAnnulusLeft)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusLeft[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2AnnulusLeft[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                #     nodeIdentifier += 1

            else:
                # Issues when elementsCountAroundOesophagus > 8
                for n2 in range(1, len(xAlongAround)):
                    if n2 == 1:
                        x1 = xAlongAround[n2][int(elementsCountAroundDuodenum * 0.5 + 1) + n]
                        x2 = xAlongAround[n2 - 1][int(elementsCountAroundDuodenum * 0.5 + 1) + n]
                    elif n2 == 2:
                        x1 = xAlongAround[n2][int(elementsCountAroundDuodenum * 0.5) + n]
                        x2 = xAlongAround[n2 - 1][int(elementsCountAroundDuodenum * 0.5 + 1) + n]
                    else:
                        x1 = xAlongAround[n2][int(elementsCountAroundDuodenum * 0.5 ) + n]
                        x2 = xAlongAround[n2 - 1][int(elementsCountAroundDuodenum * 0.5 ) + n]
                    d2 = findDerivativeBetweenPoints(x2, x1)
                    xLeft.append(x1)
                    d2Left.append(d2)

                # for n2 in range(len(xLeft)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLeft[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2Left[n2])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                #     nodeIdentifier += 1

                d2Left = interp.smoothCubicHermiteDerivativesLine(xLeft, d2Left, fixStartDirection=True)

                xLoopsLeft.append(xLeft)
                d2LoopsLeft.append(d2Left)

        # for n1 in range(len(xLoopsRight)):
        #     for n2 in range(len(xLoopsRight[n1])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopsRight[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2LoopsRight[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1
        #
        # for n1 in range(len(xLoopsLeft)):
        #     for n2 in range(len(xLoopsLeft[n1])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopsLeft[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2LoopsLeft[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        xBottom = []
        d2Bottom = []
        for n1 in range(-1, 2):
            xAlong = []
            d2Along = []
            for n2 in range(len(xAlongAround)):
                x1 = xAlongAround[n2][n1]
                x2 = xAlongAround[n2 - 1][n1]
                d2 = findDerivativeBetweenPoints(x2, x1)
                xAlong.append(x1)
                d2Along.append(d2)
            d2Along = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
            xBottom.append(xAlong)
            d2Bottom.append(d2Along)

        # for n1 in range(len(xBottom)):
        #     for n2 in range(len(xBottom[n1])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xBottom[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Bottom[n1][n2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Join row running along next to GC to oesophagus
        # Find point between second row coming down fundus and second row running along the bottom
        nx = [xLoopsLeft[-1][1], xBottom[0][0]]
        nd1 = [[xLoopsLeft[-1][1][c] - xLoopsLeft[-2][1][c] for c in range(3)], d2Bottom[0][0]]
        x, d = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]
        xBifurcationBottomLeft = x[1]

        xAround2Left = []
        d1Around2Left = []
        xAround2Left.append(xAnnulusLeft[1])
        for i in range(len(xLoopsLeft)):
            xAround2Left.append(xLoopsLeft[i][1])
        xAround2Left.append(xBifurcationBottomLeft)
        for i in range(len(xBottom[0])):
            xAround2Left.append(xBottom[0][i])

        for n1 in range(len(xAround2Left) - 1):
            x1 = xAround2Left[n1]
            x2 = xAround2Left[n1 + 1]
            d1 = findDerivativeBetweenPoints(x1, x2)
            d1Around2Left.append(d1)
        d1Around2Left.append(d1)
        d1Around2Left = interp.smoothCubicHermiteDerivativesLine(xAround2Left, d1Around2Left)

        # for n1 in range(len(xAround2Left)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAround2Left[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1Around2Left[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Right
        # Find point between second row coming down fundus and second row running along the bottom
        xAround2Right = []
        d1Around2Right = []

        nx = [xLoopsRight[-1][1], xBottom[2][0]]
        nd1 = [[xLoopsRight[-1][1][c] - xLoopsRight[-2][1][c] for c in range(3)], d2Bottom[2][0]]
        x, d = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]
        xBifurcationBottomRight = x[1]

        xAround = []
        d1Around = []
        xAround.append(xAnnulusRight[1])
        for i in range(len(xLoopsRight)):
            xAround.append(xLoopsRight[i][1])
        xAround.append(xBifurcationBottomRight)
        for i in range(len(xBottom[2])):
            xAround.append(xBottom[2][i])

        for n1 in range(len(xAround) - 1):
            x1 = xAround[n1]
            x2 = xAround[n1 + 1]
            d1 = findDerivativeBetweenPoints(x1, x2)
            d1Around.append(d1)
        d1Around.append(d1)
        d1Around = interp.smoothCubicHermiteDerivativesLine(xAround, d1Around)
        xAround2Right += xAround
        d1Around2Right += d1Around

        # Calculate d in opposite direction to get d1 in correct direction going towards oesophagus
        xAround.reverse()
        d1AroundReverse = []
        for n1 in range(len(xAround) - 1):
            x1 = xAround[n1]
            x2 = xAround[n1 + 1]
            d1 = findDerivativeBetweenPoints(x1, x2)
            d1AroundReverse.append(d1)
        d1AroundReverse.append(d1)
        d1AroundReverse = interp.smoothCubicHermiteDerivativesLine(xAround, d1AroundReverse)

        for i in range(int(elementsCountAroundDuodenum * 0.5) - 1):
            d1Around2Right[i] = d1AroundReverse[-1 - i]

        xAround.reverse()

        # for n2 in range(len(xAround2Right)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAround2Right[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1Around2Right[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Calculate d1 for loop joining to bifurcation
        xAroundBifurcationLoop = []
        xAroundBifurcationLoop.append(xAnnulusLeft[int(elementsCountAroundOesophagus * 0.25)])
        for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2):
            xAroundBifurcationLoop.append(xLoopsLeft[n1][int(elementsCountAroundOesophagus * 0.25)])
        xAroundBifurcationLoop.append(xBifurcationBottomLeft)
        xAroundBifurcationLoop.append(xOesoToDuodGC[int(elementsCountAroundDuodenum * 0.5)])
        xAroundBifurcationLoop.append(xBifurcationBottomRight)
        for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 3, -1, -1):
            xAroundBifurcationLoop.append(xLoopsRight[-int(elementsCountAroundDuodenum * 0.5 - 2) + n1][int(elementsCountAroundOesophagus * 0.25)])
        xAroundBifurcationLoop.append(xAnnulusRight[int(elementsCountAroundOesophagus * 0.25)])

        d1AroundBifurcationLoop = []
        for n1 in range(len(xAroundBifurcationLoop) - 1):
            x1 = xAroundBifurcationLoop[n1]
            x2 = xAroundBifurcationLoop[n1 + 1]
            d = findDerivativeBetweenPoints(x1, x2)
            d1AroundBifurcationLoop.append(d)
        d1AroundBifurcationLoop.append(d)
        d1AroundBifurcationLoop = interp.smoothCubicHermiteDerivativesLine(xAroundBifurcationLoop,
                                                                           d1AroundBifurcationLoop,
                                                                           fixStartDirection=True, fixEndDirection=True)

        # nodeIdentifier = 10000
        # for n1 in range(len(xAroundBifurcationLoop)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAroundBifurcationLoop[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1AroundBifurcationLoop[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Annulus to LC
        # xAnnulusToDuodLC = xOesoToDuodLC #[1:]
        # d2AnnulusToDuodLC = []
        # for n2 in range(len(xAnnulusToDuodLC) - 1):
        #     x1 = xAnnulusToDuodLC[n2]
        #     x2 = xAnnulusToDuodLC[n2 + 1]
        #     d2 = findDerivativeBetweenPoints(x1, x2)
        #     d2AnnulusToDuodLC.append(d2)
        # d2AnnulusToDuodLC.append(d2)
        #
        # d2AnnulusToDuodLC = interp.smoothCubicHermiteDerivativesLine(xAnnulusToDuodLC, d2AnnulusToDuodLC)

        xLC = xOesoToDuodLC  # [1:]
        d2LC = []
        for n2 in range(len(xLC) - 1):
            x1 = xLC[n2]
            x2 = xLC[n2 + 1]
            d2 = findDerivativeBetweenPoints(x1, x2)
            d2LC.append(d2)
        d2LC.append(d2)

        d2LC = interp.smoothCubicHermiteDerivativesLine(xLC, d2LC)
        xAnnulusToDuodLC = xLC[1:]
        d2AnnulusToDuodLC = d2LC[1:]
        # for n2 in range(len(xAnnulusToDuodLC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusToDuodLC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2AnnulusToDuodLC[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Arrange nodes
        xOuter = []
        d1Outer = []
        d2Outer = []

        for n2 in range(elementsCountAlong + 1):
            xAround = []
            d1Around = []
            d2Around = []
            if n2 == 0:
                for i in range(int(elementsCountAroundDuodenum * 0.5) - 2):
                    xAround.append(xOesoToDuodGC[i+1])
                    d1Around.append(d2OesoToDuodGC[i + 1])
                    d2Around.append(d2AnnulusLeft[0] if i == 0 else d2LoopsLeft[i-1][0])

            elif n2 == 1:
                xAround.append(xOesoToDuodGC[int(elementsCountAroundDuodenum * 0.5) - 1])
                d1Around.append(d2OesoToDuodGC[int(elementsCountAroundDuodenum * 0.5) - 1])
                d2Around.append(d2LoopsLeft[int(elementsCountAroundDuodenum * 0.5) - 3][0])
                # Right
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, - 1, -1):
                    xAround.append(xAround2Right[n1])
                    d1Around.append(d1Around2Right[n1])
                    d2Around.append(d2LoopsRight[n1 - 1][n2] if n1 > 0 else d2AnnulusRight[n2])
                # Left
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
                    xAround.append(xAround2Left[n1])
                    d1Around.append(d1Around2Left[n1])
                    d2Around.append(d2LoopsLeft[n1 - 1][n2] if n1 > 0 else d2AnnulusLeft[n2])

            elif n2 > 1 and n2 < int(elementsCountAroundOesophagus * 0.25):
                pass # Handle later

            elif n2 == int(elementsCountAroundOesophagus * 0.25): # Upstream bifurcation row
                xAround = xAroundBifurcationLoop[int(elementsCountAroundDuodenum * 0.5):]
                d1Around = d1AroundBifurcationLoop[int(elementsCountAroundDuodenum * 0.5): - 1]
                d1Around.append(d1AlongAround[0][int(elementsCountAroundDuodenum * 0.5)]) # Annulus right
                d1Around.append(d1AlongAround[0][int(elementsCountAroundDuodenum * 0.5) + 1]) # Annulus left
                xAround += xAroundBifurcationLoop[: int(elementsCountAroundDuodenum * 0.5)]
                d1Around += d1AroundBifurcationLoop[1 : int(elementsCountAroundDuodenum * 0.5)]

                d2Around.append(d2OesoToDuodGC[int(elementsCountAroundDuodenum * 0.5)]) # on GC
                d2Around.append(d1Around2Right[int(elementsCountAroundDuodenum * 0.5) - 1]) # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Around.append(d2LoopsRight[n1 - 1][n2] if n1 > 0 else d2AnnulusRight[int(elementsCountAroundOesophagus * 0.25)])
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Around.append(d2LoopsLeft[n1 - 1][n2] if n1 > 0 else d2AnnulusLeft[int(elementsCountAroundOesophagus * 0.25)])
                d2Around.append(d1Around2Left[int(elementsCountAroundDuodenum * 0.5) - 1])  # next to GC

            elif n2 == int(elementsCountAroundOesophagus * 0.25) + 1: # Downstream bifurcation row:
                xAround = xAlongAround[0][: int(elementsCountAroundDuodenum * 0.5)] + \
                          xAlongAround[0][int(elementsCountAroundDuodenum * 0.5) + 2:]
                d1Around = d1AlongAround[0][: int(elementsCountAroundDuodenum * 0.5)] + \
                          d1AlongAround[0][int(elementsCountAroundDuodenum * 0.5) + 2:]
                d2Around.append(d2OesoToDuodGC[int(elementsCountAroundDuodenum * 0.5) + 1])  # on GC
                d2Around.append(d1Around2Right[int(elementsCountAroundDuodenum * 0.5)])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, 0, -1):
                    d2Around.append(d2LoopsRight[n1 - 1][n2])
                for n1 in range(1, int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Around.append(d2LoopsLeft[n1 - 1][n2])
                d2Around.append(d1Around2Left[int(elementsCountAroundDuodenum * 0.5)])  # next to GC

            elif n2 > int(elementsCountAroundOesophagus * 0.25) + 1 and n2 <=  int(elementsCountAroundOesophagus * 0.5): # n2 == 4
                idx = n2 - (int(elementsCountAroundOesophagus * 0.25) + 1)
                xAround = xAlongAround[idx][: int(elementsCountAroundDuodenum * 0.5) + 1] + \
                          xAlongAround[idx][int(elementsCountAroundDuodenum * 0.5) + 1:]
                d1Around = d1AlongAround[idx][: int(elementsCountAroundDuodenum * 0.5) + 1] + \
                           d1AlongAround[idx][int(elementsCountAroundDuodenum * 0.5) + 1:]
                d2Around.append(d2OesoToDuodGC[int(elementsCountAroundDuodenum * 0.5) + 1 + idx])  # on GC
                d2Around.append(d1Around2Right[int(elementsCountAroundDuodenum * 0.5) + idx])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Around.append(d2LoopsRight[n1 - 1][n2] if n1 > 0 else d2AnnulusRight[int(elementsCountAroundOesophagus * 0.25) + 1])
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Around.append(d2LoopsLeft[n1 - 1][n2] if n1 > 0 else d2AnnulusLeft[int(elementsCountAroundOesophagus * 0.25) + 1])
                d2Around.append(d1Around2Left[int(elementsCountAroundDuodenum * 0.5) + 1])  # next to GC

            elif n2 > int(elementsCountAroundOesophagus * 0.5): # n2 > 4
                idx = n2 - (int(elementsCountAroundOesophagus * 0.25) + 1)
                GCIdx = int(elementsCountAroundDuodenum * 0.5) + n2 - 2
                xAround = xAlongAround[idx]
                d1Around = d1AlongAround[idx]
                d2Around.append(d2OesoToDuodGC[GCIdx])  # on GC
                d2Around.append(d1Around2Right[GCIdx - 1])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Around.append(d2LoopsRight[n1 - 1][n2] if n1 > 0 else d2AnnulusToDuodLC[n2 - 5])
                for n1 in range(1, int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Around.append(d2LoopsLeft[n1 - 1][n2])
                d2Around.append(d1Around2Left[GCIdx - 1])  # next to GC

                # for n1 in range(len(xAround)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAround[n1])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Around[n1])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Around[n1])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                #     nodeIdentifier += 1

            xOuter.append(xAround)
            d1Outer.append(d1Around)
            d2Outer.append(d2Around)

        d3UnitOuter = []
        for n2 in range(elementsCountAlong + 1):
            d3Around = []
            for n1 in range(len(xOuter[n2])):
                d3Around.append(vector.normalise(
                    vector.crossproduct3(vector.normalise(d1Outer[n2][n1]), vector.normalise(d2Outer[n2][n1]))))
            d3UnitOuter.append(d3Around)

        # Calculate curvatures
        # Curvatures along GC
        xGC = []
        dGC = []
        norms = []
        for n1 in range(len(xOuter[0])):
            xGC.append(xOuter[0][n1])
            dGC.append(d1Outer[0][n1])
            norms.append(d3UnitOuter[0][n1])
        for n2 in range(1, elementsCountAlong + 1):
            xGC.append(xOuter[n2][0])
            dGC.append(d1Outer[n2][0] if n2 == 1 else d2Outer[n2][0])
            norms.append(d3UnitOuter[n2][0])
        curvatureAlongGC = findCurvatureAlongLine(xGC, dGC, norms) # 1st len(xOuter[0]) + 1 are for d1, the rest for d2

        # for m1 in range(len(xGC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xGC[m1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dGC[m1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Curvature along left row next to GC and apply same to right side
        norms = []
        xTest = []
        for n in range(int(len(xOuter[1]) * 0.5)): # d1s
            xTest.append(xOuter[1][n + int(len(xOuter[1]) * 0.5) + 1]) # KM
            norms.append(d3UnitOuter[1][n + int(len(xOuter[1]) * 0.5) + 1])
        for n2 in range(2, elementsCountAlong + 1): # d2s
            xTest.append(xOuter[n2][-1]) # KM
            norms.append(d3UnitOuter[n2][-1])
        curvatureAlong2Left = findCurvatureAlongLine(xAround2Left, d1Around2Left, norms)

        # Curvature at GC around loops - excluding annulus points
        curvatureLoopsRight = []
        curvatureLoopsLeft = []
        curvaturesOnGC = []
        for n1 in range(len(xLoopsRight)):
            normsRight = []
            normsLeft = []
            for n2 in range(len(xLoopsRight[n1])):
                if n2 == 0:
                    norm = d3UnitOuter[0][n1 + 1] if n1 < len(xOuter[0]) - 1 else d3UnitOuter[n1 + 2 - len(xOuter[0])][0]
                    normsRight.append(norm)
                    normsLeft.append(norm)
                elif n2 == int(elementsCountAroundOesophagus * 0.25) + 1: # bifurcation downstream
                    normsRight.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) - n1])
                    normsLeft.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) + 1 + n1])
                else:
                    normsRight.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) - 1 - n1])
                    normsLeft.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) + n1 + (2 if n2 <= elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) else 1)])

            curvatureRight = findCurvatureAlongLine(xLoopsRight[n1], d2LoopsRight[n1], normsRight)
            curvatureLoopsRight.append(curvatureRight)
            curvatureLeft = findCurvatureAlongLine(xLoopsLeft[n1], d2LoopsLeft[n1], normsLeft)
            curvatureLoopsLeft.append(curvatureLeft)

            curvatureGC = 0.5 * (curvatureRight[0] + curvatureLeft[0])
            curvaturesOnGC.append(curvatureGC)

        # Curvature around annulus
        normsRight = []
        normsLeft = []
        for n2 in range(len(xAnnulusRight)):
            if n2 == 0:
                normsRight.append(d3UnitOuter[0][0])
                normsLeft.append(d3UnitOuter[0][0])
            elif n2 > int(elementsCountAroundOesophagus * 0.25):
                normsRight.append(d3UnitOuter[n2 + 1][int(len(xOuter[n2 + 1]) * 0.5)])
                normsLeft.append(d3UnitOuter[n2 + 1][int(len(xOuter[n2 + 1]) * 0.5) + (1 if n2 < elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) else 0)])
            else:
                normsRight.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5)])
                normsLeft.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) + 1])
        d2AnnulusRight[-1] = d2Outer[elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) + 1][int(elementsCountAroundDuodenum * 0.5)]
        d2AnnulusLeft[-1] = d2Outer[elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) + 1][int(elementsCountAroundDuodenum * 0.5)]
        curvatureAnnulusRight = findCurvatureAlongLine(xAnnulusRight, d2AnnulusRight, normsRight)
        curvatureAnnulusLeft = findCurvatureAlongLine(xAnnulusLeft, d2AnnulusLeft, normsLeft)
        curvatureAnnulusGC = (curvatureAnnulusRight[0] + curvatureAnnulusLeft[0]) * 0.5

        # Curvature along LC
        norms = []
        for n in range(int(elementsCountAroundOesophagus * 0.5) + 1, elementsCountAlong + 1):
            norms.append(d3UnitOuter[n][int(len(xOuter[n]) * 0.5)])
        curvatureAlongLC = findCurvatureAlongLine(xOesoToDuodLC[1:], d2OesoToDuodLC[1:], norms)

        # Assemble curvatures into matrix
        d1CurvatureOuter = []
        d2CurvatureOuter = []
        xCheck = []
        # GC -> n2 == 0
        d1CurvatureOuter.append(curvatureAlongGC[: len(xOuter[0])])
        d2CurvatureOuter.append([curvatureAnnulusGC] + curvaturesOnGC[:-1])
        xCheck.append(xGC[: len(xOuter[0])]) # KM

        # n2 == 1
        n2 = 1
        curvature = []
        xAround = []
        curvature.append(curvatureAlongGC[len(xOuter[0])])
        xAround.append(xGC[len(xOuter[0])]) # KM
        for n1 in range(int(len(xOuter[1]) * 0.5)):
            xAround.append(xOuter[1][-(1 + n1)]) # KM
            curvature.append(curvatureAlong2Left[-(1 + n1)])
        for n1 in range(int(len(xOuter[1]) * 0.5)):
            xAround.append(xOuter[1][int(len(xOuter[1]) * 0.5) + n1]) # KM
            curvature.append(curvatureAlong2Left[n1])
        xCheck.append(xAround) # KM
        d1CurvatureOuter.append(curvature)

        d2Curvature = []
        d2Curvature.append(curvatureAlongGC[-1])
        for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, - 1, -1):
            d2Curvature.append(curvatureLoopsRight[n1 - 1][n2] if n1 > 0 else curvatureAnnulusRight[n2])
        for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
            d2Curvature.append(curvatureLoopsLeft[n1 - 1][n2] if n1 > 0 else curvatureAnnulusLeft[n2])
        d2CurvatureOuter.append(d2Curvature)

        # Deal with additional elementsAroundOesophagus later -
        # n2 > 1 and n2 < int(elementsCountAroundOesophagus * 0.25)

        # Downstream of second column
        for n2 in range(int(elementsCountAroundOesophagus * 0.25), elementsCountAlong + 1): #int(elementsCountAroundOesophagus * 0.5) + 1):
            d2Curvature = []
            elementsCountAround = len(xOuter[n2])
            GCIdx = len(curvatureAlongGC) - (elementsCountAlong - n2 + 1)
            if n2 == int(elementsCountAroundOesophagus * 0.25): # upstream bifurcation - ignore bifurcation
                # re-arrange nodes to calculate curvature
                xAround = xOuter[n2][int(elementsCountAround * 0.5) + 2:] + xOuter[n2][: int(elementsCountAround * 0.5)]
                d1Around = d1Outer[n2][int(elementsCountAround * 0.5) + 2:] + d1Outer[n2][:int(elementsCountAround * 0.5)]
                norms = d3UnitOuter[n2][int(elementsCountAround * 0.5) + 2:] + d3UnitOuter[n2][:int(elementsCountAround * 0.5)]
                curvature = findCurvatureAlongLine(xAround, d1Around, norms)
                # Arrange curvatures to match original order
                x = xAround[int(elementsCountAround * 0.5) - 1: ] + [[0.0, 0.0, 0.0]] + [[0.0, 0.0, 0.0]] + xAround[:int(elementsCountAround * 0.5) - 1] # KM
                xCheck.append(x) # KM
                curvatureArranged = curvature[int(elementsCountAround * 0.5) - 1: ] + [0.0] + [0.0] + curvature[:int(elementsCountAround * 0.5) - 1]
                d1CurvatureOuter.append(curvatureArranged)

                d2Curvature.append(curvatureAlongGC[GCIdx]) # on GC
                d2Curvature.append(curvatureAlong2Left[GCIdx])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Curvature.append(curvatureLoopsRight[n1 - 1][n2] if n1 > 0 else curvatureAnnulusRight[int(elementsCountAroundOesophagus * 0.25)])
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Curvature.append(curvatureLoopsLeft[n1 - 1][n2] if n1 > 0 else curvatureAnnulusLeft[int(elementsCountAroundOesophagus * 0.25)])
                d2Curvature.append(curvatureAlong2Left[int(elementsCountAroundDuodenum * 0.5) - 1])  # next to GC
                d2CurvatureOuter.append(d2Curvature)

            elif n2 == int(elementsCountAroundOesophagus * 0.25) + 1: # downstream bifurcation - include bifurcation
                elementsCountAroundUpStreamBifurcation = len(xOuter[int(elementsCountAroundOesophagus * 0.25)])
                xAround = [xOuter[n2-1][int(elementsCountAround * 0.5) + 2]] + xOuter[n2][int(elementsCountAround * 0.5) + 1:] + xOuter[n2][: int(elementsCountAround * 0.5) + 1] + [xOuter[n2-1][int(elementsCountAround * 0.5) + 1]]
                d1Around = [d1Outer[n2-1][int(elementsCountAround * 0.5) + 2]] + d1Outer[n2][int(elementsCountAround * 0.5) + 1:] + d1Outer[n2][:int(elementsCountAround * 0.5) + 1] + [d1Outer[n2-1][int(elementsCountAround * 0.5) + 1]]
                norms = [d3UnitOuter[n2-1][int(elementsCountAround * 0.5) + 2]] + d3UnitOuter[n2][int(elementsCountAround * 0.5) + 1:] + d3UnitOuter[n2][:int(elementsCountAround * 0.5) + 1] + [d3UnitOuter[n2-1][int(elementsCountAround * 0.5) + 1]]
                curvature = findCurvatureAlongLine(xAround, d1Around, norms)
                # replace curvature at bifurcation in previous group of elements around
                xCheck[len(xCheck) - 1][int(elementsCountAroundUpStreamBifurcation * 0.5)] = xAround[-1]
                xCheck[len(xCheck) - 1][int(elementsCountAroundUpStreamBifurcation * 0.5) + 1] = xAround[0]
                d1CurvatureOuter[len(d1CurvatureOuter) - 1][int(elementsCountAroundUpStreamBifurcation * 0.5)] = curvature[-1]
                d1CurvatureOuter[len(d1CurvatureOuter) - 1][int(elementsCountAroundUpStreamBifurcation * 0.5) + 1] = curvature[0]
                # rearrange the rest
                x = xAround[int(elementsCountAround * 0.5) + 1: - 1] + xAround[1:int(elementsCountAround * 0.5)+ 1] #KM
                xCheck.append(x) # KM
                curvatureArranged = curvature[int(elementsCountAround * 0.5) + 1: - 1] + curvature[1:int(elementsCountAround * 0.5) + 1]
                d1CurvatureOuter.append(curvatureArranged)

                d2Curvature.append(curvatureAlongGC[GCIdx])  # on GC
                d2Curvature.append(curvatureAlong2Left[GCIdx])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, 0, -1):
                    d2Curvature.append(curvatureLoopsRight[n1 - 1][n2])
                for n1 in range(1, int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Curvature.append(curvatureLoopsLeft[n1 - 1][n2])
                d2Curvature.append(curvatureAlong2Left[int(elementsCountAroundDuodenum * 0.5)])  # next to GC
                d2CurvatureOuter.append(d2Curvature)

            elif n2 > int(elementsCountAroundOesophagus * 0.25) + 1 and n2 <=  int(elementsCountAroundOesophagus * 0.5):
                xAround = xOuter[n2][int(elementsCountAround * 0.5) + 1:] + xOuter[n2][: int(elementsCountAround * 0.5) + 1]
                d1Around = d1Outer[n2][int(elementsCountAround * 0.5) + 1:] + d1Outer[n2][: int(elementsCountAround * 0.5) + 1]
                norms = d3UnitOuter[n2][int(elementsCountAround * 0.5) + 1:] + d3UnitOuter[n2][: int(elementsCountAround * 0.5) + 1]
                curvature = findCurvatureAlongLine(xAround, d1Around, norms)
                x = xAround[int(elementsCountAround * 0.5): ] + xAround[:int(elementsCountAround * 0.5)] # KM
                curvatureArranged = curvature[int(elementsCountAround * 0.5):] + curvature[:int(elementsCountAround * 0.5)]
                d1CurvatureOuter.append(curvatureArranged)
                xCheck.append(x) # KM

                d2Curvature.append(curvatureAlongGC[GCIdx]) # on GC
                d2Curvature.append(curvatureAlong2Left[GCIdx])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Curvature.append(curvatureLoopsRight[n1 - 1][n2] if n1 > 0 else curvatureAnnulusRight[int(elementsCountAroundOesophagus * 0.25) + 1])
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Curvature.append(curvatureLoopsLeft[n1 - 1][n2] if n1 > 0 else curvatureAnnulusLeft[int(elementsCountAroundOesophagus * 0.25) + 1])
                d2Curvature.append(curvatureAlong2Left[int(elementsCountAroundDuodenum * 0.5) + 1])  # next to GC
                d2CurvatureOuter.append(d2Curvature)

            elif n2 > int(elementsCountAroundOesophagus * 0.5): # Downstream of oesophagus
                d1CurvatureOuter.append(findCurvatureAroundLoop(xOuter[n2], d1Outer[n2], d3UnitOuter[n2]))

                d2Curvature.append(curvatureAlongGC[GCIdx])  # on GC
                d2Curvature.append(curvatureAlong2Left[GCIdx])  # next to GC
                for n1 in range(int(elementsCountAroundDuodenum * 0.5) - 2, -1, -1):
                    d2Curvature.append(curvatureLoopsRight[n1 - 1][n2] if n1 > 0 else curvatureAlongLC[n2 - (elementsCountAlong - int(elementsCountAroundOesophagus * 0.5) + 1)])
                for n1 in range(1, int(elementsCountAroundDuodenum * 0.5) - 1):
                    d2Curvature.append(curvatureLoopsLeft[n1 - 1][n2])
                d2Curvature.append(curvatureAlong2Left[GCIdx - 1])  # next to GC
                d2CurvatureOuter.append(d2Curvature)

        # Create inner nodes
        xList = []
        d1List = []
        d2List = []
        d3List = []
        nodeIdx = bodyStartNode

        idxMat = []

        for n2 in range(elementsCountAlong + 1):
            idxThroughWall = []
            for n3 in range(elementsCountThroughWall + 1):
                xi3 = 1 / elementsCountThroughWall * n3
                idxAround = []
                for n1 in range(len(xOuter[n2])):
                    # Coordinates
                    norm = d3UnitOuter[n2][n1]
                    xOut = xOuter[n2][n1]
                    xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
                    dWall = [wallThickness * c for c in norm]
                    x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                    xList.append(x)

                    # d1
                    factor = 1.0 - wallThickness * xi3 * d1CurvatureOuter[n2][n1]
                    d1 = [factor * c for c in d1Outer[n2][n1]]
                    d1List.append(d1)

                    # d2
                    factor = 1.0 - wallThickness * xi3 * d2CurvatureOuter[n2][n1]
                    d2 = [factor * c for c in d2Outer[n2][n1]]
                    d2List.append(d2)

                    # d3
                    d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
                    d3List.append(d3)

                    idxAround.append(nodeIdx)
                    nodeIdx += 1

                idxThroughWall.append(idxAround)
            idxMat.append(idxThroughWall)

        # Assemble endPoints for annulus
        endPoints_x = [[None] * elementsCountAroundOesophagus, [None] * elementsCountAroundOesophagus]
        endPoints_d1 = [[None] * elementsCountAroundOesophagus, [None] * elementsCountAroundOesophagus]
        endPoints_d2 = [[None] * elementsCountAroundOesophagus, [None] * elementsCountAroundOesophagus]
        endNode_Id = [[None] * elementsCountAroundOesophagus, [None] * elementsCountAroundOesophagus]
        endDerivativesMap = [[None] * elementsCountAroundOesophagus, [None] * elementsCountAroundOesophagus]

        thicknessIdx = [0, -1]
        for nAround in range(elementsCountAroundOesophagus):
            for n3 in range(len(thicknessIdx)):
                if nAround == 0:
                    idx = idxMat[nAround][thicknessIdx[n3]][0]
                elif nAround <= int(elementsCountAroundOesophagus * 0.25):
                    idx = idxMat[nAround][thicknessIdx[n3]][int((len(xOuter[nAround]) - 1) * 0.5)]
                elif int(elementsCountAroundOesophagus * 0.25) < nAround < int(elementsCountAroundOesophagus * 0.5):
                    idx = idxMat[nAround + 1][thicknessIdx[n3]][int((len(xOuter[nAround + 1]) - 1) * 0.5)]
                elif nAround == int(elementsCountAroundOesophagus * 0.5):
                    idx = idxMat[nAround + 1][thicknessIdx[n3]][int(len(xOuter[nAround + 1]) * 0.5)]
                elif nAround > int(elementsCountAroundOesophagus * 0.5):
                    idx = endNode_Id[n3][int(elementsCountAroundOesophagus * 0.5) - (nAround - int(elementsCountAroundOesophagus * 0.5))] + 1

                endPoints_x[n3][nAround] = xList[idx]
                endPoints_d1[n3][nAround] = d1List[idx]
                endPoints_d2[n3][nAround] = d2List[idx]
                endNode_Id[n3][nAround] = idx

        for nAround in range(elementsCountAroundOesophagus):
            if nAround == 0:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)
            elif 0 < nAround < int(elementsCountAroundOesophagus * 0.5):
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, 1, 0), (-1, 0, 0), None)
            elif nAround == int(elementsCountAroundOesophagus * 0.5):
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = (None, None, None)
            elif int(elementsCountAroundOesophagus * 0.5) < nAround < elementsCountAroundOesophagus:
                endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)

        for n2 in range(len(xList)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n2])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier += 1

        # for n2 in range(len(xOuter)):
        #     for n1 in range(len(xOuter[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Outer[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Outer[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3UnitOuter[n2][n1])
        #         nodeIdentifier += 1

        # Create element
        mesh = fm.findMeshByDimension(3)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eftStandard = eftfactory.createEftBasic()

        elementtemplateStandard = mesh.createElementtemplate()
        elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateStandard.defineField(coordinates, -1, eftStandard)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementIdentifier = nextElementIdentifier

        # Row 1
        e2 = 0
        startNode = bodyStartNode
        elementsCountAround1 = len(xOuter[e2])
        elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
        elementsCountAround2 = len(xOuter[e2 + 1])

        for e3 in range(elementsCountThroughWall):
            for e1 in range(int(elementsCountAround1) * 2 + 1):
                if e1 != elementsCountAround1:
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e1 < elementsCountAround1:
                        scaleFactors = [-1.0]
                        if e1 == 0:
                            bni11 = startNode + elementsAroundThroughWall + e3 * elementsCountAround2 + e1
                            bni12 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1 - 1
                        else:
                            bni11 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1
                            bni12 = bni11 - 1
                        bni21 = startNode + elementsAroundThroughWall + 1 + e1 + e3 * elementsCountAround2
                        bni22 = bni21 + 1
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + (elementsCountAround2 if e1 == 0 else elementsCountAround1),
                                           bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                        eft1 = eftfactory.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scaleEftNodeValueLabels(eft1, [1, 2, 5, 6], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D2_DS1DS2,
                                                                     Node.VALUE_LABEL_D2_DS1DS3,
                                                                     Node.VALUE_LABEL_D3_DS1DS2DS3], [1])
                        scaleEftNodeValueLabels(eft1, [1, 2, 5, 6], [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                                                                     Node.VALUE_LABEL_D2_DS2DS3,
                                                                     Node.VALUE_LABEL_D3_DS1DS2DS3], [1])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 > elementsCountAround1:
                        if e1 < elementsCountAround1 * 2:
                            bni11 = startNode + e1 - elementsCountAround1 - 1 + elementsCountAround1 * e3
                            bni12 = bni11 + 1
                        else:
                            bni11 = startNode + elementsCountAround1 + e3 * elementsCountAround1 - 1
                            bni12 = startNode + elementsAroundThroughWall + e3 * elementsCountAround2
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + 1
                        bni22 = bni21 + 1
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1,
                                           bni12 + (elementsCountAround1 if e1 < elementsCountAround1 * 2 else elementsCountAround2),
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if scaleFactors:
                        result3 = element.setScaleFactors(eft1, scaleFactors)
                    elementIdentifier += 1

        # Row 2
        e2 = 1
        startNode = bodyStartNode
        for e in range(e2):
            startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

        elementsCountAround1 = len(xOuter[e2])
        elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
        elementsCountAround2 = len(xOuter[e2 + 1])

        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround1 + 2):
                if e1 != int(elementsCountAround1 * 0.5 + 1):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard

                    if e1 < 2:
                        bni11 = startNode + e3 * elementsCountAround1 + e1
                        bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1)  # % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + e1
                        bni22 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + (e1 + 1)  # % elementsCountAround2
                        if e1 == 0: # Remap derivatives of element adjacent to GC
                            scaleFactors = [-1.0]
                            nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                               bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e1 == 1: # Bottom right wedge
                            nodeIdentifiers = [bni11, bni21, bni22,
                                               bni11 + elementsCountAround1,
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                            eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 > 1 and e1 < elementsCountAround1:
                        bni11 = startNode + e3 * elementsCountAround1 + e1 - 1
                        bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                        bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                    elif e1 >= elementsCountAround1:
                        bni11 = startNode + e3 * elementsCountAround1 + e1 - 2
                        bni12 = startNode + e3 * elementsCountAround1 + (e1 - 1) % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                        bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                        if e1 == elementsCountAround1: # Bottom left wedge
                            nodeIdentifiers = [bni12, bni21, bni22,
                                               bni12 + elementsCountAround1,
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                            eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
                        elif e1 == elementsCountAround1 + 1: # Remap derivatives of element adjacent to GC
                            scaleFactors = [-1.0]
                            nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                               bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if scaleFactors:
                        result3 = element.setScaleFactors(eft1, scaleFactors)
                    elementIdentifier += 1

        # Upstream bifurcation
        e2 = int(elementsCountAroundOesophagus * 0.25)
        startNode = bodyStartNode
        for e in range(e2):
            startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

        elementsCountAround1 = len(xOuter[e2])
        elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
        elementsCountAround2 = len(xOuter[e2 + 1])

        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround1):
                if e1 != int(elementsCountAround1 * 0.5):
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    bni11 = startNode + e3 * elementsCountAround1 + e1
                    bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                    bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                    if e1 < int(elementsCountAround1 * 0.5) - 1:
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                    elif e1 == int(elementsCountAround1 * 0.5) - 1: # right wedge
                        nodeIdentifiers = [bni11, bni12, bni21,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2]
                        eft1 = eftfactory.createEftWedgeCollapseXi2([4, 8])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 == int(elementsCountAround1 * 0.5) + 1:  # left wedge
                        bni21 = bni21 - 1
                        nodeIdentifiers = [bni11, bni12, bni21,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2]
                        eft1 = eft1 = eftfactory.createEftWedgeCollapseXi2([3, 7])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 > int(elementsCountAround1 * 0.5) + 1:
                        bni21 = bni21 - 2
                        bni22 = startNode + elementsAroundThroughWall + (e1 - 1) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if scaleFactors:
                        result3 = element.setScaleFactors(eft1, scaleFactors)
                    elementIdentifier += 1

        # Downstream bifurcation
        e2 = int(elementsCountAroundOesophagus * 0.25) + 1
        startNode = bodyStartNode
        for e in range(e2):
            startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

        elementsCountAround1 = len(xOuter[e2])
        elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
        elementsCountAround2 = len(xOuter[e2 + 1])
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround1 + 1):
                if e1 < int(elementsCountAround1 * 0.5) + 1:
                    bni11 = startNode + e3 * elementsCountAround1 + e1
                elif e1 == int(elementsCountAround1 * 0.5) + 1:
                    bni11 = startNode - len(xOuter[e2-1]) * (elementsCountThroughWall + 1) + e3 * len(xOuter[e2 - 1]) + e1 + 1
                elif e1 > int(elementsCountAround1 * 0.5) + 1:
                    bni11 = startNode + e3 * elementsCountAround1 + e1 - 1

                if e1 < int(elementsCountAround1 * 0.5):
                    bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                elif e1 == int(elementsCountAround1 * 0.5):
                    bni12 = startNode - len(xOuter[e2-1]) * (elementsCountThroughWall + 1) + e3 * len(xOuter[e2-1]) + e1 + 1
                elif e1 > int(elementsCountAround1 * 0.5):
                    bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1

                if e1 > int(elementsCountAround1 * 0.5):
                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + 1
                    bni22 = startNode + elementsAroundThroughWall + (e1 + 2) % elementsCountAround2 + elementsCountAround2 * e3
                else:
                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                    bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                   bni11 + (len(xOuter[e2 - 1]) if e1 == int(elementsCountAround1 * 0.5) + 1 else elementsCountAround1),
                                   bni12 + (len(xOuter[e2 - 1]) if e1 == int(elementsCountAround1 * 0.5) else elementsCountAround1),
                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                element = mesh.createElement(elementIdentifier, elementtemplateStandard)
                result = element.setNodesByIdentifier(eftStandard, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        # Penultimate row connecting to annulus and beyond
        for e2 in range(int(elementsCountAroundOesophagus * 0.5), elementsCountAlong):
            startNode = bodyStartNode
            for e in range(e2):
                startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

            elementsCountAround1 = len(xOuter[e2])
            elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
            elementsCountAround2 = len(xOuter[e2 + 1])

            for e3 in range(elementsCountThroughWall):
                for e1 in range(elementsCountAround1 - (1 if e2 == int(elementsCountAroundOesophagus * 0.5) else 0)):
                    scaleFactors = [] # check if needed
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e2 == int(elementsCountAroundOesophagus * 0.5):
                        bni11 = startNode + e3 * elementsCountAround1 + e1 + (0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                        bni12 = startNode + e3 * elementsCountAround1 + (e1 + (1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround1
                        # Remap elements next to annulus
                        if e1 == int(elementsCountAround1 * 0.5) - 1:
                            scaleFactors = [-1.0] # check if needed
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e1 == int(elementsCountAround1 * 0.5):
                            scaleFactors = [-1.0]  # check if needed
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                    else:
                        bni11 = startNode + e3 * elementsCountAround1 + e1
                        bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                    bni22 = startNode + elementsAroundThroughWall + (e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                       bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if scaleFactors:
                        result3 = element.setScaleFactors(eft1, scaleFactors)
                    elementIdentifier += 1

        # Annulus
        nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
            endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap)

        for n2 in range(len(xTrackSurface)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1TrackSurface[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2TrackSurface[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            nodeIdentifier += 1

        fm.endChange()

        annotationGroup = []
        return annotationGroup

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

def findClosestPositionAndDerivativeOnTrackSurface(x, nx, trackSurface, nxProportion1, elementsCountAlongTrackSurface):
    """
    Find the closest position and derivative around the tracksurface of a point sitting near the fundus of stomach.
    Use a startPosition to improve the search for nearest position on the track surface as the fundus has a curved
    and complex track surface.

    :param x: coordinates of point of interest
    :param nx: coordinates of points along curve where point of interest lies.
    :param trackSurface: track surface where point sits
    :param nxProportion1: proportion around track surface of curve
    :param elementsCountAlongTrackSurface: number of elements along track surface
    :return: position and derivative of point around track surface
    """
    closestIdxOnNx = interp.getNearestPointIndex(nx, x)
    closestPositionToPoint = trackSurface.createPositionProportion(nxProportion1, closestIdxOnNx / elementsCountAlongTrackSurface)
    xPosition = trackSurface.findNearestPosition(x, closestPositionToPoint)
    d = trackSurface.evaluateCoordinates(xPosition, derivatives=True)[1]

    return xPosition, d

def getSmoothedSampledPointsOnTrackSurface(trackSurface, startProportion1, startProportion2, endProportion1,
                                           endProportion2, elementsOut, startDerivative = None, endDerivative = None,
                                           startDerivativeMagnitude = None, endDerivativeMagnitude = None):
    """
    Create smoothly spaced out hermite curve points between two points a and b on the surface,
    each defined by their proportions over the surface in directions 1 and 2.
    :param trackSurface: track surface
    :param startProportion1, startProportion2: proportion of start point in direction around and along track surface
    :param endProportion1, endProportion2: proportion of end point in direction around and along track surface
    :param elementsOut: number of elements out
    :param startDerivative, endDerivative: optional derivative vectors in 3-D world coordinates
        to match at the start and end of the curves. If omitted, fits in with other derivative or is
        in a straight line from a to b.
    :param derivativeMagnitudeStart, derivativeMagnitudeEnd: optional magnitude of derivatives to match at the start and
        end of the curves.
    :return: coordinates and derivative of sampled points
    """

    mx, md2, md1, md3, mProportions = \
        trackSurface.createHermiteCurvePoints(startProportion1, startProportion2, endProportion1, endProportion2,
                                              elementsOut, startDerivative, endDerivative)

    xSampled, dSampled = trackSurface.resampleHermiteCurvePointsSmooth(mx, md2, md1, md3, mProportions,
                                                                       startDerivativeMagnitude,
                                                                       endDerivativeMagnitude)[0:2]
    return xSampled, dSampled

def findDerivativeBetweenPoints(v1, v2):
    """

    :param v1:
    :param v2:
    :return:
    """
    d = [v2[c] - v1[c] for c in range(3)]
    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
    d = [c * arcLengthAround for c in vector.normalise(d)]

    return d

def findCurvatureAroundLoop(nx, nd, radialVectors):
    """

    :param nx:
    :param nd:
    :param radialVectors:
    :return:
    """
    curvature = []
    for n in range(len(nx)):
        prevIdx = n - 1 if n > 0 else -1
        nextIdx = n + 1 if n < len(nx) - 1 else 0
        kappam = interp.getCubicHermiteCurvature(nx[prevIdx], nd[prevIdx], nx[n], nd[n], radialVectors[n], 1.0)
        kappap = interp.getCubicHermiteCurvature(nx[n], nd[n], nx[nextIdx], nd[nextIdx], radialVectors[n], 0.0)
        curvature.append(0.5 * (kappam + kappap))

    return curvature

def findCurvatureAlongLine(nx, nd, radialVectors):
    """

    :param nx:
    :param nd:
    :param radialVectors:
    :return:
    """
    curvature = []
    for n in range(len(nx)):
        if n == 0:
            curvature.append(interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0))
        elif n == len(nx) - 1:
            curvature.append(interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0))
        else:
            curvature.append(0.5 * (
                        interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0) +
                        interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0)))

    return curvature

