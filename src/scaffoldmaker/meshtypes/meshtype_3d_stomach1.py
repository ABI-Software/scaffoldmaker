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
from scaffoldmaker.utils.bifurcation import get_bifurcation_triple_point
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
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[70.8, 72.3, 0.0], [-15.1, -43.8, 0.0], [37.3, -18.7, 0.0], [2.2, 7.2, 0.0], [0.0, 0.0, 38.7],
                     [0.0, 0.0, 4.5]],
                    [[51.1, 30.7, 0.0], [-21.4, -33.7, 0.0], [33.9, -17.5, 0.0], [-9.0, -4.8, 0.0], [0.0, 0.0, 41.4],
                     [0.0, 0.0, 0.9]],
                    [[27.8, 2.7, 0.0], [-27.7, -18.4, 0.0], [20.5, -27.0, 0.0], [-17.1, -3.5, 0.0], [0.0, 0.0, 40.9],
                     [0.0, 0.0, -5.0]],
                    [[-0.3, -5.8, 0.0], [-28.0, -2.6, 0.0], [0.4, -25.7, 0.0], [-12.6, 3.7, 0.0], [0.0, 0.0, 32.3],
                     [0.0, 0.0, -9.4]],
                    [[-26.5, -2.9, 0.0], [-20.6, 8.0, 0.0], [-5.4, -19.9, 0.0], [-5.6, 7.7, 0.0], [0.0, 0.0, 22.1],
                     [0.0, 0.0, -6.8]],
                    [[-40.6, 7.4, 0.0], [-11.0, 12.5, 0.0], [-10.9, -10.9, 0.0], [-1.3, 7.9, 0.0], [0.0, 0.0, 17.6],
                     [0.0, 0.0, -5.9]],
                    [[-48.1, 21.0, 0.0], [-6.4, 15.7, 0.0], [-8.4, -4.0, 0.0], [0.1, 4.9, 0.0], [0.0, 0.0, 10.5],
                     [0.0, 0.0, -3.1]],
                    [[-52.9, 38.6, 0.0], [-3.2, 19.3, 0.0], [-11.2, -1.5, 0.0], [-5.7, 0.1, 0.0], [0.0, 0.0, 12.0],
                     [0.0, 0.0, 6.1]]])

        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 6
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[12.4, 11.9, 0.0], [2.3, -12.3, 0.0], [8.7, -1.4, 0.0], [0.5, -3.7, 0.0], [0.0, 0.0, 7.4],
                     [0.0, 0.0, 1.6]],
                    [[9.5, -0.2, 0.0], [-6.0, -9.6, 0.0], [6.6, -5.1, 0.0], [-3.9, -4.3, 0.0], [0.0, 0.0, 8.5],
                     [0.0, 0.0, 0.7]],
                    [[2.3, -5.2, 0.0], [-9.3, -2.2, 0.0], [2.1, -9.1, -0.0], [-6.4, -1.6, 0.0], [0.0, 0.0, 9.0],
                     [0.0, 0.0, -0.1]],
                    [[-7.8, -3.9, 0.0], [-7.5, 4.2, 0.0], [-6.0, -8.0, 0.0], [-3.4, 3.4, 0.0], [0.0, 0.0, 8.1],
                     [0.0, 0.0, -1.5]],
                    [[-11.7, 1.7, 0.0], [-1.7, 7.2, 0.0], [-6.4, -3.5, 0.0], [1.4, 4.0, 0.0], [0.0, 0.0, 6.2],
                     [0.0, 0.0, -2.8]],
                    [[-10.7, 9.4, 0.0], [0.1, 6.6, 0.0], [-2.9, 0.0, 0.0], [1.1, 1.3, 0.0], [0.0, 0.0, 2.4],
                     [0.0, 0.0, -1.0]],
                    [[-11.3, 14.8, 0.0], [-1.2, 4.1, 0.0], [-3.5, -0.3, 0.0], [-2.3, -1.9, 0.0], [0.0, 0.0, 3.4],
                     [0.0, 0.0, 3.0]]])
        }),
        'Stomach coordinates': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[17.9, 0.0, 0.0], [-21.6, 0.0, 0.0], [0.0, -9.0, -0.7], [0.0, -1.5, 1.2], [0.0, -0.7, 9.0],
                     [0.0, 1.2, 1.5]],
                    [[0.0, 0.0, 0.0], [-14.2, 0.0, 0.0], [0.0, -9.0, -0.1], [0.0, 1.4, 0.0], [0.0, -0.1, 9.0],
                     [0.0, 0.0, -1.4]],
                    [[-10.5, 0.0, 0.0], [-9.2, 0.0, 0.0], [0.0, -6.8, -0.4], [0.0, 3.5, 0.0], [0.0, -0.4, 6.8],
                     [0.0, 0.0, -3.5]],
                    [[-18.4, 0.0, 0.0], [-7.8, 0.0, 0.0], [0.0, -2.4, -0.1], [0.0, 1.6, 0.2], [0.0, -0.1, 2.4],
                     [0.0, 0.2, -1.6]],
                    [[-26.2, 0.0, 0.0], [-7.8, 0.0, 0.0], [0.0, -3.6, 0.1], [0.0, -4.0, 0.2], [0.0, 0.1, 3.6],
                     [0.0, 0.2, 4.0]]])
        }),
    }

    ostiumDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 25.0,
                'Ostium length': 15.0,
                'Ostium wall thickness': 5.0,
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 5.0,
                'Vessel wall thickness': 5.0,
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
        'Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 4.0,
                'Ostium length': 3.5,
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
            'Human 1',
            'Rat 1',
            'Stomach coordinates']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        elif 'Generic 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Generic 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        elif 'Stomach coordinates' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Stomach coordinates']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around oesophagus': 8,
            'Number of elements around duodenum': 12,
            'Number of elements between annulus and duodenum': 6,
            'Number of elements through wall': 1,
            'Number of radial elements in annulus': 1,  # KM
            'Wall thickness': 5.0,
            'Limiting ridge': False,
            'Fundus end position along factor': 0.3,
            'Gastro-oesophagal junction': copy.deepcopy(ostiumOption),
            'Gastro-oesophagal junction position along factor': 0.35,
            'Annulus derivative factor': 1.0,
            'Show track surface': False, # KM
            'Make stomach': True,  # KM
            'Show central path': False, # KM
            'Use cross derivatives': False,
            'Use linear through wall' : False, # need to deal with wedge not available in bicubichermite
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        if 'Rat 1' in parameterSetName:
            options['Number of elements around oesophagus'] = 12
            options['Number of elements around duodenum'] = 14
            options['Number of elements between annulus and duodenum'] = 2
            options['Wall thickness'] = 0.5
            options['Gastro-oesophagal junction position along factor'] = 0.55
            options['Annulus derivative factor'] = 0.2
            options['Limiting ridge'] = True
            options['Fundus end position along factor'] = 0.5
        elif 'Stomach coordinates' in parameterSetName:
            options['Wall thickness'] = 0.5
            options['Gastro-oesophagal junction position along factor'] = 0.45
            options['Fundus end position along factor'] = 0.45

        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around oesophagus',
            'Number of elements around duodenum',
            'Number of elements between annulus and duodenum',
            'Number of elements through wall',
            'Number of radial elements in annulus',
            'Wall thickness',
            'Limiting ridge',
            'Fundus end position along factor',
            'Gastro-oesophagal junction',
            'Gastro-oesophagal junction position along factor',
            'Annulus derivative factor',
            'Show track surface',
            'Make stomach',
            'Show central path',
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
        if options['Number of elements around oesophagus'] % 4 > 0:
            options['Number of elements around oesophagus'] = options['Number of elements around oesophagus'] // 4 * 4
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        if options['Number of elements between annulus and duodenum'] < 2:
            options['Number of elements between annulus and duodenum'] = 2
        for key in [
            'Number of elements around oesophagus',
            'Number of elements around duodenum']:
            if options[key] % 2:
                options[key] += 1
        # if options['Number of elements along'] < 8:
        #     options['Number of elements along'] = 8
        for key in [
            'Fundus end position along factor',
            'Annulus derivative factor']:
            if options[key] < 0:
                options[key] = 0.0
        for key in [
            'Number of elements through wall',
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
        elementsCountAroundOeso = options['Number of elements around oesophagus']
        elementsCountAroundDuod = options['Number of elements around duodenum']
        elementsAlongAnnulusToDuod = options['Number of elements between annulus and duodenum']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not (options['Use linear through wall'])

        GOJPositionAlongFactor = options['Gastro-oesophagal junction position along factor']
        GOJOptions = options['Gastro-oesophagal junction']
        GOJSettings = GOJOptions.getScaffoldSettings()
        useLimitingRidge = options['Limiting ridge']
        fundusEndPositionAlongFactor = options['Fundus end position along factor']
        elementsCountAnnulus = options['Number of radial elements in annulus']
        annulusDerivativeFactor = options['Annulus derivative factor']

        showTrackSurface = options['Show track surface']
        makeStomach = options['Make stomach']
        showCentralPath = options['Show central path']

        elementsCountAlongTrackSurface = 20
        elementsAroundHalfOeso = int(elementsCountAroundOeso * 0.5)
        elementsAroundQuarterOeso = int(elementsCountAroundOeso * 0.25)
        elementsAroundHalfDuod = int(elementsCountAroundDuod * 0.5)

        # if useLimitingRidge: # Rat / mouse
        #     annulusDerivativeFactor = 0.2
        # else:
        #     annulusDerivativeFactor = 1.0

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
        #     cd3New = vector.setMagnitude(cd3[i], vector.magnitude(cd2[i]))
        #     print(i, vector.magnitude(cd2[i]), vector.magnitude(cd3[i]), cd3New)
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], ',', cd3[i], ',', cd13[i], '],')

        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongTrackSurface)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        # Make sampled d2 and d3 normal to central path
        # d2Check = []
        # for c in range(len(sx)):
        #     td2 = vector.vectorRejection(sd2[c], sd1[c])
        #     sd2[c] = vector.setMagnitude(td2, vector.magnitude(sd2[c]))
        #     d2Check.append(matrix.rotateAboutZAxis(sd2[c], math.pi))
        #     td3 = vector.vectorRejection(sd3[c], sd1[c])
        #     sd3[c] = vector.setMagnitude(td3, vector.magnitude(sd3[c]))

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
        for n1 in range(elementsCountAroundDuod):
            rotAngle = n1 * 2.0 * math.pi / elementsCountAroundDuod
            rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(sd2[0]), vector.normalise(sd3[0]))) # vector.normalise(sd1[0])
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
            d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            d2Apex.append(d2Rot)

        xEllipses = []
        d1Ellipses = []
        for n in range(elementsAlongFundus + 1, len(sx)):
            px, pd1 = createEllipsePoints(sx[n], 2 * math.pi, sd2[n], sd3[n], elementsCountAroundDuod,
                                          startRadians=0.0)
            xEllipses.append(px)
            d1Ellipses.append(pd1)

        # Find d2
        d2Raw = []
        for n1 in range(elementsCountAroundDuod):
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
            for n1 in range(elementsCountAroundDuod):
                d2 = d2Raw[n1][n2]
                d2Around.append(d2)
            d2Ellipses.append(d2Around)

        # for n2 in range(len(xEllipses)):
        #     for n1 in range(elementsCountAroundDuod):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xEllipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Ellipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Ellipses[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Merge fundus and body
        xAll = [[sx[0]] * elementsCountAroundDuod] + xEllipses
        d2All = [d2Apex] + d2Ellipses

        # for n2 in range(len(xAll)):
        #     for n1 in range(elementsCountAroundDuod):
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
        for n1 in range(elementsCountAroundDuod):
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
            for n1 in range(elementsCountAroundDuod):
                x = xRaw[n1][n2]
                d2 = d2Raw[n1][n2]
                xAround.append(x)
                d2Around.append(d2)

                # Calculate d1
                if n2 > 0:
                    v1 = xRaw[n1][n2]
                    v2 = xRaw[n1 + 1 if n1 < elementsCountAroundDuod - 2 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1Around.append(d1)
                else:
                    d1Around.append(d2Raw[int(elementsCountAroundDuod * 0.75)][0])

            if n2 > 0:
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
            else:
                d1Smoothed = d1Around

            xSampledAll.append(xAround)
            d1SampledAll.append(d1Smoothed)
            d2SampledAll.append(d2Around)

        # for n2 in range(elementsCountAlongTrackSurface + 1):
        #     for n1 in range(elementsCountAroundDuod):
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
            for n1 in range(elementsCountAroundDuod):
                xTrackSurface.append(xSampledAll[n2][n1])
                d1TrackSurface.append(d1SampledAll[n2][n1])
                d2TrackSurface.append(d2SampledAll[n2][n1])

        trackSurfaceStomach = TrackSurface(elementsCountAroundDuod, elementsCountAlongTrackSurface,
                                           xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

        # Set up gastro-oesophagal junction
        GOJSettings['Number of elements around ostium'] = elementsCountAroundOeso
        GOJPosition = trackSurfaceStomach.createPositionProportion(0.5, GOJPositionAlongFactor)
        xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(GOJPosition, derivatives=True)
        axis1 = d1Centre

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        nextNodeIdentifier = nodeIdentifier
        nextElementIdentifier = elementIdentifier
        nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
            generateOstiumMesh(region, GOJSettings, trackSurfaceStomach, GOJPosition, axis1,
                               nextNodeIdentifier, nextElementIdentifier)

        bodyStartNode = nextNodeIdentifier
        nodeIdentifier = nextNodeIdentifier

        fundusEndPosition = trackSurfaceStomach.createPositionProportion(0.0, fundusEndPositionAlongFactor)
        xFundusEnd, d1FundusEnd, d2FundusEnd = trackSurfaceStomach.evaluateCoordinates(fundusEndPosition, derivatives=True)

        # node = nodes.createNode(nodeIdentifier, nodetemplate)
        # cache.setNode(node)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFundusEnd)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2FundusEnd)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1FundusEnd)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        # nodeIdentifier += 1

        elementsAlongFundus = elementsAroundQuarterOeso + (0 if useLimitingRidge else 1)
        elementsCountAlong = elementsAlongAnnulusToDuod + elementsAroundHalfOeso + 1

        # From fundus end to duodenum
        elementsAlongFundusEndToDuod = elementsCountAlong - elementsAlongFundus
        xAlongDownFundus = []
        d2AlongDownFundus = []
        xAlongDownFundus.append(xFundusEnd)
        d2AlongDownFundus.append(d2FundusEnd)
        startIdx = fundusEndPositionAlongFactor*elementsCountAlongTrackSurface
        startIdx = math.ceil(startIdx) + (1 if startIdx - math.ceil(startIdx) == 0 else 0)
        for n2 in range(startIdx, len(xSampledAll)):
            xAlongDownFundus.append(xSampledAll[n2][0])
            d2AlongDownFundus.append(d2SampledAll[n2][0])

        # From oesophagus to fundus end
        elementsAlongGCFromOesoToFundusEnd = elementsAroundHalfDuod - 2 + elementsAlongFundus
        xAlongUpFundus = []
        d2AlongUpFundus = []
        xAlongUpFundus.append(o1_x[1][0])
        d2AlongUpFundus.append([annulusDerivativeFactor * c for c in o1_d2[1][0]])

        # From oesophagus to fundus apex
        oesoStartProportion2 = trackSurfaceStomach.getProportion(o1_Positions[0])[1]
        elementsAlongUpstreamOfOeso = int(elementsCountAlongTrackSurface * oesoStartProportion2)
        arcLengthOesoApex = 0.0
        for n2 in range(elementsAlongUpstreamOfOeso):
            nAlong = elementsAlongUpstreamOfOeso - n2
            v1 = xSampledAll[nAlong][elementsAroundHalfDuod]
            v2 = xSampledAll[nAlong - 1][elementsAroundHalfDuod]
            d = [v2[c] - v1[c] for c in range(3)]
            arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
            arcLengthOesoApex += arcLengthAround
            d2 = [c * arcLengthAround for c in vector.normalise(d)]
            xAlongUpFundus.append(v1)
            d2AlongUpFundus.append(d2)

        # From fundus apex to end
        for n2 in range(startIdx - (0 if useLimitingRidge else 1)):
            xAlongUpFundus.append(xSampledAll[n2][0])
            d2AlongUpFundus.append(d2SampledAll[n2][0])

        if useLimitingRidge:
            # Sample from oesophagus to fundus end
            # xAlongGCOesoToFundusEnd, d2AlongGCOesoToFundusEnd = \
            #     interp.sampleCubicHermiteCurves(xAlongUpFundus, d2AlongUpFundus, elementsAlongGCFromOesoToFundusEnd,
            #                                     addLengthStart=0.5 * vector.magnitude(d2AlongUpFundus[0]),
            #                                     lengthFractionStart=0.5, arcLengthDerivatives=True)[0:2]
            xAlongGCOesoToFundusEnd, d2AlongGCOesoToFundusEnd = \
                interp.sampleCubicHermiteCurvesSmooth(xAlongUpFundus, d2AlongUpFundus,
                                                      elementsAlongGCFromOesoToFundusEnd,
                                                      derivativeMagnitudeStart=vector.magnitude(d2AlongUpFundus[0]))[0:2]

            # Sample from limiting ridge to duodenum
            xAlongDownFundus[0] = xAlongGCOesoToFundusEnd[-1]
            d2AlongDownFundus[0] = d2AlongGCOesoToFundusEnd[-1]
            xAlongGCFundusEndToDuod, d2AlongGCFundusEndToDuod = \
                interp.sampleCubicHermiteCurves(xAlongDownFundus, d2AlongDownFundus, elementsAlongFundusEndToDuod,
                                                addLengthStart=0.5 * vector.magnitude(d2AlongDownFundus[0]),
                                                lengthFractionStart=0.5, arcLengthDerivatives=True)[0:2]

        else:
            # Sample from end of fundus to duodenum
            xAlongGCFundusEndToDuod, d2AlongGCFundusEndToDuod = \
                interp.sampleCubicHermiteCurves(xAlongDownFundus, d2AlongDownFundus, elementsAlongFundusEndToDuod,
                                                arcLengthDerivatives=True)[0:2]

            # Sample from ostium to end of fundus
            # Provide first element size so that we can transit
            xAlongUpFundus[-1] = xAlongGCFundusEndToDuod[0]
            d2AlongUpFundus[-1] = d2AlongGCFundusEndToDuod[0]
            # xAlongGCOesoToFundusEnd, d2AlongGCOesoToFundusEnd = \
            #     interp.sampleCubicHermiteCurves(xAlongUpFundus, d2AlongUpFundus, elementsAlongGCFromOesoToFundusEnd,
            #                                     addLengthStart=0.5* vector.magnitude(d2AlongUpFundus[0]),
            #                                     lengthFractionStart = 0.5, addLengthEnd=0.5 * vector.magnitude(d2AlongUpFundus[-1]),
            #                                     lengthFractionEnd=0.5, arcLengthDerivatives=True)[0:2]
            xAlongGCOesoToFundusEnd, d2AlongGCOesoToFundusEnd = \
                interp.sampleCubicHermiteCurvesSmooth(xAlongUpFundus, d2AlongUpFundus, elementsAlongGCFromOesoToFundusEnd,
                                                derivativeMagnitudeStart= vector.magnitude(d2AlongUpFundus[0]),
                                                derivativeMagnitudeEnd = vector.magnitude(d2AlongUpFundus[-1]))[0:2]

        xOesoToDuodGC = xAlongGCOesoToFundusEnd[:-1] + xAlongGCFundusEndToDuod
        d2OesoToDuodGC = d2AlongGCOesoToFundusEnd[:-1] + d2AlongGCFundusEndToDuod

        # for n2 in range(len(xAlongUpFundus)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongUpFundus[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2AlongUpFundus[n2])
        #     nodeIdentifier += 1
        #
        # for n2 in range(len(xAlongGCFundusEndToDuod)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongGCFundusEndToDuod[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2AlongGCFundusEndToDuod[n2])
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

        arcLength = 0.0
        for e in range(len(xOesoToDuodGC) - 1):
            arcLength += interp.getCubicHermiteArcLength(xOesoToDuodGC[e], d2OesoToDuodGC[e],
                                                         xOesoToDuodGC[e + 1], d2OesoToDuodGC[e + 1])
            if arcLength > arcLengthOesoApex:
                nodesCountFromOesoToApex = e + 2
                break

        ptsOnTrackSurfaceGC = []
        xAlongAround = []
        d1AlongAround = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceGC.append(xSampledAll[n2][0])

        # Rings adjacent to triple point (excluding triple point rings) to 6 point junction
        # First half
        for n2 in range(elementsAroundQuarterOeso + 1, elementsAroundHalfOeso):
            ostiumIdx = n2
            GCIdx = elementsAroundHalfDuod - 1 + n2
            GCPosition, d1GC = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                              trackSurfaceStomach, 0.0,
                                                                              elementsCountAlongTrackSurface)
            GCProportion1, GCProportion2 = trackSurfaceStomach.getProportion(GCPosition)

            endPosition = o1_Positions[ostiumIdx]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(o1_d1[1][ostiumIdx]), math.pi)
            d2 = o1_d2[1][ostiumIdx]
            d1EndOstium = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
            d1EndTrackSurface = trackSurfaceStomach.evaluateCoordinates(endPosition, derivatives=True)[1]

            xFirstHalf, d1FirstHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, 0.0, GCProportion2, endProportion1,
                                                       endProportion2, elementsAroundHalfDuod + 1,
                                                       startDerivative=d1GC, endDerivative=d1EndOstium,
                                                       endDerivativeMagnitude=annulusDerivativeFactor * vector.magnitude(d2))

            # for n2 in range(len(xFirstHalf)):
            #     node = nodes.createNode(nodeIdentifier, nodetemplate)
            #     cache.setNode(node)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFirstHalf[n2])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1FirstHalf[n2])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            #     nodeIdentifier += 1

            # Second half
            ostiumIdx2 = -n2
            startPosition = o1_Positions[ostiumIdx2]
            d1StartOstium = o1_d2[1][ostiumIdx2]
            startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
            d1StartTrackSurface = trackSurfaceStomach.evaluateCoordinates(startPosition, derivatives=True)[1]

            xSecondHalf, d1SecondHalf = \
                getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1, startProportion2, 1.0,
                                                       GCProportion2, elementsAroundHalfDuod + 1, startDerivative=d1StartOstium,
                                                       endDerivative=d1GC,
                                                       startDerivativeMagnitude=annulusDerivativeFactor * vector.magnitude(
                                                           d1StartOstium))

            xAround = xFirstHalf[:-1] + xSecondHalf[1:-1]
            d1Around = d1FirstHalf[:-1] + d1SecondHalf[1:-1]

            xAlongAround.append(xAround)
            d1AlongAround.append(d1Around)

        # Elements downstream of 6 pt junction
        xTests = [] # KM
        for idx in range(-(elementsCountAlong - elementsAroundHalfOeso - 1), 0):
            # Search for point on central path and use that to make ellipse
            xStart = xOesoToDuodGC[idx]
            startPosition, d1Start = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[idx],
                                                                                    ptsOnTrackSurfaceGC,
                                                                                    trackSurfaceStomach, 0.0,
                                                                                    elementsCountAlongTrackSurface)
            startProportion2 = trackSurfaceStomach.getProportion(startPosition)[1]
            if 0.0 < startProportion2 < 1.0:
                closestIdxOnCentralPath = interp.getNearestPointIndex(sx, xStart)
                if 0 < closestIdxOnCentralPath < len(sx) - 1:
                    # Check if xStart is closer to upstream or downstream of closestIdx
                    xOnGCPrevElem = [sx[closestIdxOnCentralPath - 1][c] + sd2[closestIdxOnCentralPath - 1][c] for c in range(3)]
                    distBetweenXOnGCPrevElem = vector.magnitude([xStart[c] - xOnGCPrevElem[c] for c in range(3)])
                    xOnGCNextElem = [sx[closestIdxOnCentralPath + 1][c] + sd2[closestIdxOnCentralPath + 1][c] for c in range(3)]
                    distBetweenXOnGCNextElem = vector.magnitude([xStart[c] - xOnGCNextElem[c] for c in range(3)])
                    eiLowerLimit = closestIdxOnCentralPath - (1 if distBetweenXOnGCNextElem > distBetweenXOnGCPrevElem else 0)

                elif closestIdxOnCentralPath == len(sx) - 1:
                    eiLowerLimit = closestIdxOnCentralPath - 1

                elif closestIdxOnCentralPath == 0:
                    eiLowerLimit = closestIdxOnCentralPath

                xiLowerLimit = 0.0
                xiUpperLimit = 1.0
                xOnGCLowerLimit = [sx[eiLowerLimit][c] + sd2[eiLowerLimit][c] for c in range(3)]
                xOnGCUpperLimit = [sx[eiLowerLimit + 1][c] + sd2[eiLowerLimit + 1][c] for c in range(3)]
                distBetweenXAndXStartPrev = vector.magnitude([xStart[c] - xOnGCLowerLimit[c] for c in range(3)])

                # Search for point on central path which is orthogonal to xStart
                tol = 1e-8
                xLowerLimit = sx[eiLowerLimit]
                d1LowerLimit = sd1[eiLowerLimit]
                d2LowerLimit = sd2[eiLowerLimit]
                d12LowerLimit = sd12[eiLowerLimit]
                xUpperLimit = sx[eiLowerLimit + 1]
                d1UpperLimit = sd1[eiLowerLimit + 1]
                d2UpperLimit = sd2[eiLowerLimit + 1]
                d12UpperLimit = sd12[eiLowerLimit + 1]

                for iter in range(100):
                    xiGuess = 0.5 * (xiLowerLimit + xiUpperLimit)
                    x = interp.interpolateCubicHermite(xLowerLimit, d1LowerLimit, xUpperLimit, d1UpperLimit, xiGuess)
                    d2 = interp.interpolateCubicHermite(d2LowerLimit, d12LowerLimit, d2UpperLimit, d12UpperLimit, xiGuess)
                    xGuess = [x[c] + d2[c] for c in range(3)]
                    distBetweenXAndXStart = vector.magnitude([xStart[c] - xGuess[c] for c in range(3)])
                    distBetweenXOnGCLowerLimitAndXStart = vector.magnitude([xStart[c] - xOnGCLowerLimit[c] for c in range(3)])
                    distBetweenXOnGCUpperLimitAndXStart = vector.magnitude([xStart[c] - xOnGCUpperLimit[c] for c in range(3)])

                    if abs(distBetweenXAndXStart - distBetweenXAndXStartPrev) < tol:
                        xProjection = x
                        xiProjection = xiGuess
                        break
                    elif distBetweenXOnGCLowerLimitAndXStart < distBetweenXOnGCUpperLimitAndXStart:
                        xiUpperLimit = xiGuess
                        xOnGCUpperLimit = xGuess
                        distBetweenXAndXStartPrev = distBetweenXAndXStart
                    else:
                        xiLowerLimit = xiGuess
                        xOnGCLowerLimit = xGuess
                        distBetweenXAndXStartPrev = distBetweenXAndXStart

                if iter > 98:
                    print('Search for projection on central path - Max iters reached:', iter)

                axis1 = interp.interpolateCubicHermite(sd2[eiLowerLimit], sd12[eiLowerLimit],
                                                       sd2[eiLowerLimit + 1], sd12[eiLowerLimit + 1], xiProjection)
                axis2 = interp.interpolateCubicHermite(sd3[eiLowerLimit], sd13[eiLowerLimit],
                                                       sd3[eiLowerLimit + 1], sd13[eiLowerLimit + 1], xiProjection)
                xAround, d1Around = createEllipsePoints(xProjection, 2 * math.pi, axis1, axis2, elementsCountAroundDuod, startRadians=0.0)

            elif startProportion2 == 0:
                xAround, d1Around = createEllipsePoints(sx[0], 2 * math.pi, sd2[0], sd3[0],
                                                        elementsCountAroundDuod, startRadians=0.0)
            elif startProportion2 == 1.0:
                xAround, d1Around = createEllipsePoints(sx[-1], 2 * math.pi, sd2[-1], sd3[-1],
                                                        elementsCountAroundDuod, startRadians=0.0)

            d1Around = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)

            # Average adjacent ring with first downstream ring that is not adjacent to oesophagus
            if idx == -(elementsCountAlong - elementsAroundHalfOeso - 1):
                xAve = []
                dAve = []
                xAve.append(xOesoToDuodGC[idx - 1])

                for n in range(1, elementsCountAroundDuod):
                    startPosition = trackSurfaceStomach.findNearestPosition(xAround[n])
                    startProportion1, startProportion2 = trackSurfaceStomach.getProportion(startPosition)
                    if n == elementsAroundHalfDuod: # Point to search for 6 point junction - Work here if need to move 6pt junction closer to annulus
                        endPosition = o1_Positions[elementsAroundHalfOeso]
                    else:
                        endPosition = trackSurfaceStomach.findNearestPosition(xAlongAround[-1][n + (0 if n < elementsAroundHalfDuod else 1)])
                    endProportion1, endProportion2 = trackSurfaceStomach.getProportion(endPosition)
                    xSampled = \
                        getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, startProportion1,
                                                               startProportion2, endProportion1, endProportion2, 2)[0]
                                                               #endDerivativeMagnitude = annulusDerivativeFactor * vector.magnitude(o1_d2[1][elementsAroundHalfOeso]))[0]

                    xAve.append(xSampled[1])

                for n in range(len(xAve)):
                    v1 = xAve[n]
                    v2 = xAve[(n + 1) % len(xAve)]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    dAve.append(d1)
                dAve = interp.smoothCubicHermiteDerivativesLoop(xAve, dAve)

                xAlongAround.append(xAve)
                d1AlongAround.append(dAve)

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

        # Sample 2 loops next to annulus from point on GC to point on first ring on xAlongAround
        ptsOnTrackSurfaceOesoToFundus = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            ptsOnTrackSurfaceOesoToFundus.append(xSampledAll[n2][elementsAroundHalfDuod])

        xLoopsRight = []
        xLoopsLeft = []
        d2LoopsRight = []  # KM
        d2LoopsLeft = []  # KM

        for nLoop in range(1, elementsAroundHalfDuod - 1): #3):
            GCIdx = nLoop + 1
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

            for nSide in range(2):
                if nSide == 0:
                    xEnd = xAlongAround[0][elementsAroundHalfDuod - nLoop]
                    d2End = [xAlongAround[1][elementsAroundHalfDuod - nLoop][c] -
                             xAlongAround[0][elementsAroundHalfDuod - nLoop][c] for c in range(3)]
                else:
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d2OesoToDuodGC[GCIdx]), math.pi)
                    d2GCRot = [rotFrame[j][0] * d2GC[0] + rotFrame[j][1] * d2GC[1] + rotFrame[j][2] * d2GC[2] for j in range(3)]
                    d2GC = d2GCRot

                    xEnd = xAlongAround[0][elementsAroundHalfDuod + 1 + nLoop]
                    d2End = [xAlongAround[1][elementsAroundHalfDuod + nLoop][c] -
                             xAlongAround[0][elementsAroundHalfDuod + (0 if elementsCountAroundOeso > 8 else 1) + nLoop][c] for c in range(3)]

                nx = [xOesoToDuodGC[GCIdx], xEnd]
                nd2 = [d2GC, d2End]
                x, d2 = interp.sampleCubicHermiteCurves(nx, nd2, elementsAroundQuarterOeso + 2, arcLengthDerivatives=True)[0:2]

                # Find closest sampled points onto track surface
                xProjectedPoints = []
                d2ProjectedPoints = []
                xProjectedPoints.append(xOesoToDuodGC[GCIdx])
                d2ProjectedPoints.append(d2GC)
                for n2 in range(1, len(x)):
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[n2])
                    xProjected = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    xProjectedPoints.append(xProjected)

                for n2 in range(1, len(xProjectedPoints) - 1):
                    d2 = findDerivativeBetweenPoints(xProjectedPoints[n2], xProjectedPoints[n2 + 1])
                    d2ProjectedPoints.append(d2)
                d2ProjectedPoints.append(d2End)

                # Sample points again
                xLoop, d2Loop = interp.sampleCubicHermiteCurves(xProjectedPoints, d2ProjectedPoints,
                                                                elementsAroundQuarterOeso + 2,
                                                                addLengthEnd=0.5 * vector.magnitude(
                                                                   d2ProjectedPoints[-1]),
                                                                lengthFractionEnd=0.5,
                                                                arcLengthDerivatives=True)[0:2]
                (xLoopsRight if nSide == 0 else xLoopsLeft).append(xLoop)
                (d2LoopsRight if nSide == 0 else d2LoopsLeft).append(d2Loop)

        # for n2 in range(len(xLoopsRight)):
        #     for n1 in range(len(xLoopsRight[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopsRight[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2LoopsRight[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1
        #
        # for n2 in range(len(xLoopsLeft)):
        #     for n1 in range(len(xLoopsLeft[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopsLeft[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2LoopsLeft[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Find triple point
        xTriplePts = [[None], [None]]  # Right, left
        d1TriplePts = [[None], [None]]
        d2TriplePts = [[None], [None]]
        d3TriplePtsNorm = [[None], [None]]

        for nSide in range(2):
            ostiumIdx = elementsAroundQuarterOeso if nSide == 0 else -elementsAroundQuarterOeso
            p1x = o1_x[1][ostiumIdx]
            d = o1_d2[1][ostiumIdx]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(o1_d1[1][ostiumIdx], math.pi)
            p1d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
            p1d = [annulusDerivativeFactor * c for c in p1d]

            xLoops = xLoopsRight if nSide == 0 else xLoopsLeft
            p2x = xLoops[0][elementsAroundQuarterOeso + 1]  # downstream bifurcation
            p2d = findDerivativeBetweenPoints(xLoops[0][elementsAroundQuarterOeso + 1],
                                              xLoops[1][elementsAroundQuarterOeso + 1])

            p3x = xLoops[0][elementsAroundQuarterOeso]
            p3d = findDerivativeBetweenPoints(xLoops[0][elementsAroundQuarterOeso],
                                              xLoops[1][elementsAroundQuarterOeso])

            xTriplePts[nSide], d1TriplePts[nSide], d2TriplePts[nSide] = get_bifurcation_triple_point(p1x, p1d,
                                                                                                     p2x, p2d,
                                                                                                     p3x, p3d)
            d3TriplePtsNorm[nSide] = vector.normalise(
                vector.crossproduct3(vector.normalise(d1TriplePts[nSide]),
                                     vector.normalise(d2TriplePts[nSide])))

            # Make sure triple point is on track surface
            triplePointPosition = trackSurfaceStomach.findNearestPosition(xTriplePts[nSide])
            xTriplePts[nSide] = trackSurfaceStomach.evaluateCoordinates(triplePointPosition)

            # node = nodes.createNode(nodeIdentifier, nodetemplate)
            # cache.setNode(node)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, p1x)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, p1d)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            # nodeIdentifier += 1
            #
            # node = nodes.createNode(nodeIdentifier, nodetemplate)
            # cache.setNode(node)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, p2x)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, p2d)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            # nodeIdentifier += 1
            #
            # node = nodes.createNode(nodeIdentifier, nodetemplate)
            # cache.setNode(node)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, p3x)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, p3d)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            # nodeIdentifier += 1

            # node = nodes.createNode(nodeIdentifier, nodetemplate)
            # cache.setNode(node)
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTriplePts[nSide])
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1TriplePts[nSide])
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2TriplePts[nSide])
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3TriplePtsNorm[nSide])
            # nodeIdentifier += 1

        # Plan A - Sample points from GC to bottom of loops to create nodes running on row 2
        xBifurcationRings = []
        d1BifurcationRings = []
        xUp = []
        d1Up = []
        for n2 in range(elementsAroundQuarterOeso):
            xAround = []
            d1Around = []
            xAroundRight = []
            d1AroundRight = []
            xAroundLeft = []
            d1AroundLeft = []
            loopIdx = n2 + 2
            ostiumIdx = loopIdx + (0 if n2 < elementsAroundQuarterOeso - 1 else -1)
            GCIdx = elementsAlongGCFromOesoToFundusEnd - elementsAroundQuarterOeso + n2 + (2 if useLimitingRidge else 1)
            d1GC = findClosestPositionAndDerivativeOnTrackSurface(xOesoToDuodGC[GCIdx], ptsOnTrackSurfaceGC,
                                                                  trackSurfaceStomach, 0.0,
                                                                  elementsCountAlongTrackSurface)[1]
            for nSide in range(2):
                if nSide == 0: # Right side
                    xAroundRight.append(xOesoToDuodGC[GCIdx])
                    d1AroundRight.append(d1GC)
                    xOnLastLoopRight = xLoopsRight[-1][loopIdx]
                    d1OnLastLoopRight = findDerivativeBetweenPoints(xLoopsRight[-1][loopIdx], xLoopsRight[-2][loopIdx])

                    nx = [xOesoToDuodGC[GCIdx], xOnLastLoopRight]
                    nd1 = [d1GC, d1OnLastLoopRight]
                    x, d1 = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]

                    # Find closest sampled points onto track surface
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[1])
                    x[1] = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    d1[1] = findDerivativeBetweenPoints(x[1], x[2])

                    # Sample points again
                    x, d1 = interp.sampleCubicHermiteCurves(x, d1, 2, arcLengthDerivatives=True)[0:2]

                    xAroundRight.append(x[1])
                    d1AroundRight.append(d1[1])

                    for n in range(len(xLoopsRight) - 1):
                        xAroundRight.append(xLoopsRight[-(1 + n)][loopIdx])
                        d1AroundRight.append(findDerivativeBetweenPoints(xLoopsRight[-(1 + n)][loopIdx], xLoopsRight[-(1 + n + 1)][loopIdx]))

                    if loopIdx < elementsAroundQuarterOeso:  # additional elements upstream of triple point
                        xLoop = xLoopsRight[0][loopIdx]
                        xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                        xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                        xOstium = o1_x[1][ostiumIdx]
                        ostiumPosition = o1_Positions[ostiumIdx]
                        ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                        d = findDerivativeBetweenPoints(xLoop, xOstium)
                        endDerivativeMag = vector.magnitude(o1_d2[1][ostiumIdx]) * annulusDerivativeFactor
                        xSampled, dSampled = \
                            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, xLoopProportion1,
                                                                   xLoopProportion2, ostiumProportion1,
                                                                   ostiumProportion2, 2, # startDerivative=d,
                                                                   endDerivativeMagnitude=endDerivativeMag)[0:2]
                        xAroundRight += xSampled[:2]
                        d1AroundRight += dSampled[:2]

                    else: # connected to triple point
                        xAroundRight += [xLoopsRight[0][loopIdx]] + [xTriplePts[0]]
                        d1AroundRight += [findDerivativeBetweenPoints(xLoopsRight[0][loopIdx], xTriplePts[0])] + \
                                         [d1TriplePts[0]]

                    # for n1 in range(len(xAroundRight)):
                    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                    #     cache.setNode(node)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAroundRight[n1])
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d1AroundRight[n1])
                    #     nodeIdentifier += 1

                else: # left side
                    if loopIdx < elementsAroundQuarterOeso:  # additional elements upstream of triple point
                        xLoop = xLoopsLeft[0][loopIdx]
                        xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                        xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                        xOstium = o1_x[1][-ostiumIdx]
                        ostiumPosition = o1_Positions[-ostiumIdx]
                        ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                        d = findDerivativeBetweenPoints(xOstium, xLoop)
                        startDerivativeMag = vector.magnitude(o1_d2[1][-ostiumIdx]) * annulusDerivativeFactor
                        xSampled, dSampled = \
                            getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, ostiumProportion1, ostiumProportion2,
                                                                   xLoopProportion1, xLoopProportion2, 2, # endDerivative=d,
                                                                   startDerivativeMagnitude=startDerivativeMag)[0:2]
                        xAroundLeft.append(xSampled[1])
                        d1AroundLeft.append(dSampled[1])
                    else:
                        xAroundLeft.append(xTriplePts[1])
                        d1AroundLeft.append(d1TriplePts[1])

                    for n in range(len(xLoopsLeft) - 1):
                        xAroundLeft.append(xLoopsLeft[n][loopIdx])
                        d1AroundLeft.append(findDerivativeBetweenPoints(xLoopsLeft[n][loopIdx], xLoopsLeft[n + 1][loopIdx]))

                    xOnLastLoopLeft = xLoopsLeft[-1][loopIdx]
                    d1OnLastLoopLeft = findDerivativeBetweenPoints(xLoopsLeft[-2][loopIdx], xLoopsLeft[-1][loopIdx])

                    nx = [xOnLastLoopLeft, xOesoToDuodGC[GCIdx]]
                    nd1 = [d1OnLastLoopLeft, d1GC]
                    x, d1 = interp.sampleCubicHermiteCurves(nx, nd1, 2, arcLengthDerivatives=True)[0:2]

                    # Find closest sampled points onto track surface
                    projectedPosition = trackSurfaceStomach.findNearestPosition(x[1])
                    x[1] = trackSurfaceStomach.evaluateCoordinates(projectedPosition)
                    d1[1] = findDerivativeBetweenPoints(x[1], x[2])

                    # Sample points again
                    x, d1 = interp.sampleCubicHermiteCurves(x, d1, 2, arcLengthDerivatives=True)[0:2]

                    xAroundLeft += [xOnLastLoopLeft] + [x[1]] + [xOesoToDuodGC[GCIdx]]
                    d1AroundLeft += [findDerivativeBetweenPoints(xOnLastLoopLeft, x[1])] + [d1[1]] + [d1GC]

                    # for n1 in range(len(xAroundLeft)):
                    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                    #     cache.setNode(node)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAroundLeft[n1])
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d1AroundLeft[n1])
                    #     nodeIdentifier += 1

            xAround = xAroundRight + xAroundLeft[:-1]
            d1Around = d1AroundRight + d1AroundLeft[:-1]

            xUp.append(xAround)
            d1Up.append(d1Around)

            if loopIdx >= elementsAroundQuarterOeso:
                xBifurcationRings.append(xAround)
                d1BifurcationRings.append(d1Around)

        # for n2 in range(len(xBifurcationRings)):
        #     for n1 in range(len(xBifurcationRings[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xBifurcationRings[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1BifurcationRings[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # for n2 in range(len(xUp)):
        #     for n1 in range(len(xUp[n2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xUp[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Up[n2][n1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Row 2
        xRow2Right = []
        d1Row2Right = []
        xRow2Left = []
        d1Row2Left = []

        for nSide in range(2):
            loopIdx = 1
            ostiumIdx = 1
            if nSide == 0:
                xRow2Right.append(xUp[0][1])
                d1Row2Right.append(findDerivativeBetweenPoints(xUp[0][1], xLoopsRight[-1][1]))
                # Append rows upwards in loops
                for n in range(len(xLoopsRight) - 1):
                    xRow2Right.append(xLoopsRight[-(1 + n)][1])
                    d1Row2Right.append(
                        findDerivativeBetweenPoints(xLoopsRight[-(1 + n)][1], xLoopsRight[-(1 + n + 1)][1]))

                xLoop = xLoopsRight[0][loopIdx]
                xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                xOstium = o1_x[1][ostiumIdx]
                ostiumPosition = o1_Positions[ostiumIdx]
                ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                d = findDerivativeBetweenPoints(xLoop, xOstium)
                endDerivativeMag = vector.magnitude(o1_d2[1][ostiumIdx]) * annulusDerivativeFactor
                xSampled, dSampled = \
                    getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, xLoopProportion1,
                                                           xLoopProportion2, ostiumProportion1,
                                                           ostiumProportion2, 2, # startDerivative=d,
                                                           endDerivativeMagnitude=endDerivativeMag)[0:2]
                xRow2Right += xSampled[0:2]
                d1Row2Right += [findDerivativeBetweenPoints(xSampled[0], xSampled[1])] + [dSampled[1]] # dSampled[0:2]

            else:
                xLoop = xLoopsLeft[0][loopIdx]
                xLoopPosition = trackSurfaceStomach.findNearestPosition(xLoop)
                xLoopProportion1, xLoopProportion2 = trackSurfaceStomach.getProportion(xLoopPosition)
                xOstium = o1_x[1][-ostiumIdx]
                ostiumPosition = o1_Positions[-ostiumIdx]
                ostiumProportion1, ostiumProportion2 = trackSurfaceStomach.getProportion(ostiumPosition)
                d = findDerivativeBetweenPoints(xOstium, xLoop)
                startDerivativeMag = vector.magnitude(o1_d2[1][-ostiumIdx]) * annulusDerivativeFactor

                xSampled, dSampled = \
                    getSmoothedSampledPointsOnTrackSurface(trackSurfaceStomach, ostiumProportion1, ostiumProportion2,
                                                           xLoopProportion1, xLoopProportion2, 2, # endDerivative=d,
                                                           startDerivativeMagnitude=startDerivativeMag)[0:2]
                xRow2Left += xSampled[1:]
                d1Row2Left += [dSampled[1]] + [findDerivativeBetweenPoints(xSampled[1], xSampled[2])]

                for n in range(1, len(xLoopsLeft)):
                    xRow2Left.append(xLoopsLeft[n][loopIdx])
                    d1Row2Left.append(findDerivativeBetweenPoints(xLoopsLeft[n-1][loopIdx], xLoopsLeft[n][loopIdx]))

                xRow2Left.append(xUp[0][-1])
                d1Row2Left.append(findDerivativeBetweenPoints(xLoopsLeft[-1][1], xUp[0][-1]))

        # for n1 in range(len(xRow2Right)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xRow2Right[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d1Row2Right[n1])
        #     # print(vector.magnitude(d1Row2Right[n1]), vector.magnitude(d1Row2Left[-(1+n1)]))
        #     nodeIdentifier += 1
        #
        # for n1 in range(len(xRow2Left)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xRow2Left[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d1Row2Left[n1])
        #     nodeIdentifier += 1

        # Smooth derivatives from triple point to 6 point junction
        # Start from GC at upstream bifurcation ring to annulus to 6 point junction ring on right then left
        xLoopTripleTo6Pt = []
        dLoopTripleTo6Pt = []

        xLoopTripleTo6Pt += xBifurcationRings[0][0:int(len(xBifurcationRings[0]) * 0.5) + 1]
        for n2 in range(elementsAroundQuarterOeso - 1):
            xLoopTripleTo6Pt.append(xAlongAround[n2][int(len(xAlongAround[n2]) * 0.5)])
            junctionIdx = n2 + 1
        xLoopTripleTo6Pt += xAlongAround[junctionIdx][int(len(xAlongAround[junctionIdx]) * 0.5):] + \
                            xAlongAround[junctionIdx][0: int(len(xAlongAround[junctionIdx]) * 0.5 + 1)]
        for n2 in range(elementsAroundQuarterOeso - 1): # Note order here - going upstream
            idx = junctionIdx - 1 - n2
            xLoopTripleTo6Pt.append(xAlongAround[idx][int(len(xAlongAround[idx]) * 0.5) + 1])
        xLoopTripleTo6Pt += xBifurcationRings[0][int(len(xBifurcationRings[0]) * 0.5 + 1):]

        for n in range(len(xLoopTripleTo6Pt)):
            d = findDerivativeBetweenPoints(xLoopTripleTo6Pt[n], xLoopTripleTo6Pt[(n+1) % len(xLoopTripleTo6Pt)])
            dLoopTripleTo6Pt.append(d)
        dSmoothLoopTripleTo6Pt = interp.smoothCubicHermiteDerivativesLoop(xLoopTripleTo6Pt, dLoopTripleTo6Pt)
        # curvature - DONE

        # for n1 in range(len(xLoopTripleTo6Pt)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopTripleTo6Pt[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dSmoothLoopTripleTo6Pt[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero) #dLoopGCTriplePt[n1])
        #     nodeIdentifier += 1

        # Smooth derivatives around top loop
        # Starts from GC at downstream bifurcation ring to annulus and back
        xLoopGCTriplePt = []
        dLoopGCTriplePt = []

        xLoopGCTriplePt += xBifurcationRings[1][:int(len(xBifurcationRings[1]) * 0.5) + 1]

        for n2 in range(elementsAroundQuarterOeso - 2):
            idx = -(3 + n2)
            xLoopGCTriplePt.append(xUp[idx][int(len(xUp[idx]) * 0.5)])

        xLoopGCTriplePt += [xRow2Right[-1]] + [xOesoToDuodGC[1]] + [xRow2Left[0]]

        for n2 in range(elementsAroundQuarterOeso - 2):
            xLoopGCTriplePt.append(xUp[n2][int(len(xUp[n2]) * 0.5) + 1])

        xLoopGCTriplePt += xBifurcationRings[1][int(len(xBifurcationRings[1]) * 0.5) + 1:]

        for n in range(len(xLoopGCTriplePt)):
            d = findDerivativeBetweenPoints(xLoopGCTriplePt[n], xLoopGCTriplePt[(n+1) % len(xLoopGCTriplePt)])
            dLoopGCTriplePt.append(d)
        dSmoothLoopGCTriplePt = interp.smoothCubicHermiteDerivativesLoop(xLoopGCTriplePt, dLoopGCTriplePt)
        # curvature - DONE

        # for n1 in range(len(xLoopGCTriplePt)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoopGCTriplePt[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dSmoothLoopGCTriplePt[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero) #dLoopGCTriplePt[n1])
        #     nodeIdentifier += 1

        # Resample loop before 6 point junction
        for n in range(1, elementsCountAroundDuod + 1):
            if elementsAroundHalfDuod <= n <= elementsAroundHalfDuod + 1:
                pass
            else:
                if elementsCountAroundOeso > 8:
                    xStart = xAlongAround[-(elementsAlongAnnulusToDuod + 3)][n]
                else:
                    xStart = xUp[-1][n]
                xEnd = xAlongAround[-(elementsAlongAnnulusToDuod + 1)][n - (1 if n > elementsAroundHalfDuod else 0)]

                startPosition = trackSurfaceStomach.findNearestPosition(xStart)
                d2Start = trackSurfaceStomach.evaluateCoordinates(startPosition, derivatives=True)[2]

                endPosition = trackSurfaceStomach.findNearestPosition(xEnd)
                d2End = trackSurfaceStomach.evaluateCoordinates(endPosition, derivatives=True)[2]

                xSampled = interp.sampleCubicHermiteCurves([xStart, xEnd], [d2Start, d2End], 2)[0]
                xNewPosition = trackSurfaceStomach.findNearestPosition(xSampled[1])
                xNew = trackSurfaceStomach.evaluateCoordinates(xNewPosition)
                xAlongAround[-(elementsAlongAnnulusToDuod + 2)][n] = xNew

        # xStarts = []
        # d2Starts = []
        # xEnds = []
        # d2Ends = []
        # for n in range(elementsAroundHalfDuod + 2, elementsCountAroundDuod + 1):
        #     if elementsCountAroundOeso > 8:
        #         xStart = xAlongAround[-(elementsAlongAnnulusToDuod + 3)][n]
        #     else:
        #         xStart = xUp[-1][n]
        #     xEnd = xAlongAround[-(elementsAlongAnnulusToDuod + 1)][n - 1]
        #
        #     startPosition = trackSurfaceStomach.findNearestPosition(xStart)
        #     d2Start = trackSurfaceStomach.evaluateCoordinates(startPosition, derivatives=True)[2]
        #     xStarts.append(xStart)
        #     d2Starts.append(d2Start)
        #
        #     endPosition = trackSurfaceStomach.findNearestPosition(xEnd)
        #     d2End = trackSurfaceStomach.evaluateCoordinates(endPosition, derivatives=True)[2]
        #     xEnds.append(xEnd)
        #     d2Ends.append(d2End)
        #
        #     xSampled = interp.sampleCubicHermiteCurves([xStart, xEnd], [d2Start, d2End], 2)[0]
        #     xNewPosition = trackSurfaceStomach.findNearestPosition(xSampled[1])
        #     xNew = trackSurfaceStomach.evaluateCoordinates(xNewPosition)
        #
        #     xAlongAround[-(elementsAlongAnnulusToDuod + 2)][n] = xNew


        # Assemble nodes and d1
        xOuter = []
        d1Outer = []
        countUp = 0
        countDown = 0
        for n2 in range(elementsCountAlong + 1):
            xAround = []
            d1Around = []
            d2Around = []
            if n2 == 0:
                # print(n2, 'GC')
                for i in range(elementsAroundHalfDuod - 2):
                    xAround.append(xOesoToDuodGC[i + 1])
                    d1Around.append(d2OesoToDuodGC[i + 1])

            elif n2 == 1:
                # print(n2, 'Row 2')
                xAround = [xOesoToDuodGC[i + n2 + 1]] + xRow2Right[1:] + xRow2Left[:-1]
                d1Around = [d2OesoToDuodGC[i + n2 + 1]] + d1Row2Right[1:] + d1Row2Left[:-1] # need to replace d1Row2 after smoothing - DONE

            elif n2 > 1 and n2 < elementsAroundQuarterOeso + 2:
                # print(n2, 'Before triple point + triple point')
                xAround = xUp[countUp]
                if n2 < elementsAroundQuarterOeso: # upstream of triple pt
                    # smooth d1 around - Make into function?
                    d1Around = d1Up[countUp]
                    xLoop = xAround[int(len(xAround) * 0.5 + 1): ] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1): ] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    d1LoopSmooth = interp.smoothCubicHermiteDerivativesLine(xLoop, d1Loop, fixStartDerivative=True,
                                                                            fixEndDerivative=True)
                    # Need to do curvature and rearrange to correct order
                    d1Around = []
                    d1Around = d1LoopSmooth[int(len(xAround) * 0.5) : ] + d1LoopSmooth[: int(len(xAround) * 0.5) : ]

                elif n2 == elementsAroundQuarterOeso: # upstream bifurcation
                    # take smoothed d1 from dSmoothTripleTo6Pt
                    d1Around = dSmoothLoopTripleTo6Pt[: int(len(xBifurcationRings[0]) * 0.5) + 1] + \
                               dSmoothLoopTripleTo6Pt[-int(len(xBifurcationRings[0]) * 0.5) : ]

                elif n2 > elementsAroundQuarterOeso: # downstream bifurcation
                    # take smoothed d1 from dSmoothGCToTriplePt
                    d1Around = dSmoothLoopGCTriplePt[: int(len(xBifurcationRings[1]) * 0.5) + 1] + \
                               dSmoothLoopGCTriplePt[-int(len(xBifurcationRings[1]) * 0.5) : ]
                countUp += 1

            elif n2 > elementsAroundQuarterOeso + 1:
                # print(n2, 'Downstream of triple point')
                xAround = xAlongAround[countDown]
                d1Around = d1AlongAround[countDown]

                # smooth d1 around - Make into function?
                if n2 < elementsAroundHalfOeso + 1:
                    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    d1LoopSmooth = interp.smoothCubicHermiteDerivativesLine(xLoop, d1Loop, fixStartDerivative=True,
                                                                            fixEndDerivative=True)
                    # Need to do curvature and rearrange to correct order
                    d1Around = []
                    d1Around = d1LoopSmooth[int(len(xAround) * 0.5):] + d1LoopSmooth[: int(len(xAround) * 0.5):]

                elif n2 == elementsAroundHalfOeso + 1: # 6 point junction ring
                    # take smoothed d1 from dSmoothedTripleTo6Pt
                    startRightIdx = int(len(xBifurcationRings[0]) * 0.5 + elementsAroundQuarterOeso + len(xAlongAround[junctionIdx]) * 0.5)
                    endRightIdx = startRightIdx + int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    startLeftIdx = startRightIdx - int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    d1Around = dSmoothLoopTripleTo6Pt[startRightIdx: endRightIdx] + \
                               dSmoothLoopTripleTo6Pt[startLeftIdx : startRightIdx]
                countDown += 1

            xOuter.append(xAround)
            d1Outer.append(d1Around)

        # for m2 in range(len(xOuter)):
        #     for m1 in range(len(xOuter[m2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[m2][m1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Outer[m2][m1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #         nodeIdentifier += 1

        # Calculate d2
        xRegularLoops = []
        d2RegularLoops = []
        d2RegularOrderedLoops = []

        for n1 in range(elementsAroundHalfDuod - 2):
            xRegularLoop = []
            d1RegularRightLoop = []
            d2RegularLoop = []
            for n2 in range(elementsCountAlong):
                idx = -(1 + n2)
                xRegularLoop.append(xOuter[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
                d1RegularRightLoop.append(d1Outer[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
            xRegularLoop.append(xOesoToDuodGC[n1 + 2])
            for n2 in range(elementsCountAlong):
                xRegularLoop.append(xOuter[n2 + 1][int(len(xOuter[n2 + 1]) * 0.5 + n1 + (1 if n2 >= elementsAroundHalfOeso else 2))])

            for n in range(len(xRegularLoop) - 1):
                d = findDerivativeBetweenPoints(xRegularLoop[n], xRegularLoop[n+1])
                d2RegularLoop.append(d)
            d2RegularLoop.append(d)

            d2SmoothRegularLoop = interp.smoothCubicHermiteDerivativesLine(xRegularLoop, d2RegularLoop)
            d2SmoothRegularOrderedLoop = copy.deepcopy(d2SmoothRegularLoop)
            # curvature - DONE

            # Switch direction on right side
            for n2 in range(elementsCountAlong):
                rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1RegularRightLoop[n2]), d2SmoothRegularLoop[n2]))
                rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                d = d2SmoothRegularLoop[n2]
                d2SmoothRegularLoop[n2] = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]

            xRegularLoops.append(xRegularLoop)
            d2RegularLoops.append(d2SmoothRegularLoop)
            d2RegularOrderedLoops.append(d2SmoothRegularOrderedLoop)

        # for m1 in range(len(xRegularLoops)):
        #     for m2 in range(len(xRegularLoops[m1])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xRegularLoops[m1][m2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2RegularOrderedLoops[m1][m2])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero) #d2RegularLoops[m1][m2])
        #         nodeIdentifier += 1

        # Smooth d2 along row 2
        xLoop2Right = []
        d1Loop2Right = []
        d2Loop2Right = []

        for n2 in range(len(xAlongAround) + len(xUp) - 1):
            idx = -(1 + n2)
            xLoop2Right.append(xOuter[idx][1])
            d1Loop2Right.append(d1Outer[idx][1])
        xLoop2Right += xRow2Right
        d1Loop2Right += d1Row2Right

        for n in range(len(xLoop2Right) - 1):
            d = findDerivativeBetweenPoints(xLoop2Right[n], xLoop2Right[n + 1])
            d2Loop2Right.append(d)
        d2Loop2Right.append(d1Row2Right[-1])

        d2Loop2Right = interp.smoothCubicHermiteDerivativesLine(xLoop2Right, d2Loop2Right, fixEndDirection=True)
        # curvature - USING LEFT SIDE - DONE

        # Switch direction of d2 for downstream nodes
        for n2 in range(len(xAlongAround) + len(xUp)):
            rotAxis = vector.normalise(
                vector.crossproduct3(vector.normalise(d1Loop2Right[n2]), d2Loop2Right[n2]))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
            d = d2Loop2Right[n2]
            d2Loop2Right[n2] = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
            idxSwitchToD1 = n2

        # for m in range(len(xLoop2Right)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoop2Right[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero) #d2Loop2Right[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2Loop2Right[m])
        #     nodeIdentifier += 1

        # Left
        xLoop2Left = []
        d2Loop2Left = []

        xLoop2Left += xRow2Left
        for n2 in range(3, len(xOuter)):
            xLoop2Left.append(xOuter[n2][-1])

        d2Loop2Left.append(d1Row2Left[0])
        for n in range(1, len(xLoop2Left) - 1):
            d = findDerivativeBetweenPoints(xLoop2Left[n], xLoop2Left[n + 1])
            d2Loop2Left.append(d)
        d2Loop2Left.append(d)

        d2Loop2Left = interp.smoothCubicHermiteDerivativesLine(xLoop2Left, d2Loop2Left, fixStartDirection=True)
        # curvature - DONE

        # for m in range(len(xLoop2Left)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLoop2Left[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Loop2Left[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Smooth lower curvature
        xLC = []
        d2LC = []
        for n2 in range(elementsAroundHalfOeso + 1, elementsCountAlong + 1):
            xLC.append(xOuter[n2][int(len(xOuter[n2]) * 0.5)])

        for n in range(len(xLC) - 1):
            d = findDerivativeBetweenPoints(xLC[n], xLC[n + 1])
            d2LC.append(d)
        d2LC.append(d)

        d2LC = interp.smoothCubicHermiteDerivativesLine(xLC, d2LC, fixStartDirection=True)
        # curvature - DONE

        # for m in range(len(xLC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xLC[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2LC[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Smooth greater curvature
        d2GC = []
        for n in range(len(xOesoToDuodGC) - 1):
            d = findDerivativeBetweenPoints(xOesoToDuodGC[n], xOesoToDuodGC[n + 1])
            d2GC.append(d)
        d2GC.append(d)
        d2GC = interp.smoothCubicHermiteDerivativesLine(xOesoToDuodGC, d2GC, fixStartDirection=True)

        # for m in range(len(xOesoToDuodGC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOesoToDuodGC[m])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2GC[m])
        #     nodeIdentifier += 1

        # Update d1 for upstream nodes
        for n1 in range(1, len(xRow2Right)):
            d1Outer[1][n1] = d2Loop2Right[idxSwitchToD1 + n1]
        for n1 in range(1, len(xRow2Left)):
            d1Outer[1][int(len(d1Outer[1]) * 0.5) + n1] = d2Loop2Left[n1 - 1]

        # Assemble d2
        d2Outer = []
        for n2 in range(elementsCountAlong + 1):
            d2Around = []
            xAround = [] # KM
            if n2 == 0:
                # print(n2, 'GC')
                xAround.append(xLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5)]) #KM
                d2Around.append(dSmoothLoopGCTriplePt[int(len(dSmoothLoopGCTriplePt) * 0.5)])

                for n1 in range(len(xOuter[0]) - 1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5)])
                    d2Around.append(d2RegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5)])
                    nextIdx = n1 + 1

            elif n2 == 1:
                # print(n2, 'Row 2')
                xAround.append(xRegularLoops[nextIdx][int(len(xRegularLoops[nextIdx]) * 0.5)])
                d2Around.append(d2RegularLoops[nextIdx][int(len(xRegularLoops[nextIdx]) * 0.5)])

                for n1 in range(nextIdx, -1, -1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) - n2])
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                # right point on annulus
                xAround.append(xLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2]) # KM
                d2 = dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2]
                rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5)]), vector.normalise(d2)))
                rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                d2Around.append([rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)])

                # left point on annulus
                xAround.append(xLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2])  # KM
                d2Around.append(dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2])

                for n1 in range(nextIdx + 1):
                        xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) + n2])
                        d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

            elif n2 > 1 and n2 < elementsAroundQuarterOeso + 2:
                # print(n2, 'Before triple point + triple point')
                # GC
                xAround.append(xOesoToDuodGC[len(xOuter[0]) + n2])
                d2Around.append(d2GC[len(xOuter[0]) + n2])

                # Row 2 right
                xAround.append(xLoop2Right[-(len(xOuter[0]) + n2)])
                d2Around.append(d2Loop2Right[-(len(xOuter[0]) + n2)])

                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) - n2])
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                # Annulus right
                xAround.append(xLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2 + (1 if n2 > elementsAroundQuarterOeso else 0)])  # KM
                d2 = dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) - n2 + (1 if n2 > elementsAroundQuarterOeso else 0)]
                if n2 <= elementsAroundQuarterOeso: # Rotate to point towards duodenum
                    rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5)]),
                                             vector.normalise(d2)))
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                    d2Around.append(
                        [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)])
                else: # just take d2 as-is cos we are going to remove this point later
                    d2Around.append(d2)

                # Annulus left
                xAround.append(xLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2 - (1 if n2 > elementsAroundQuarterOeso else 0)])  # KM
                d2Around.append(dSmoothLoopGCTriplePt[int(len(xLoopGCTriplePt) * 0.5) + n2 - (1 if n2 > elementsAroundQuarterOeso else 0)])

                # Regular down left
                for n1 in range(nextIdx + 1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) + n2])
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

                # Row 2 left
                xAround.append(xLoop2Left[len(xOuter[0]) + n2 - 1])
                d2Around.append(d2Loop2Left[len(xOuter[0]) + n2 - 1])

                # for m1 in range(len(d2Around)):
                #     node = nodes.createNode(nodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAround[m1])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2Around[m1])
                #     nodeIdentifier += 1

            elif n2 > elementsAroundQuarterOeso + 1:
                # print(n2, 'Downstream of triple point')
                # GC
                xAround.append(xOesoToDuodGC[len(xOuter[0]) + n2])
                d2Around.append(d2GC[len(xOuter[0]) + n2])

                # Row 2 right
                xAround.append(xLoop2Right[-(len(xOuter[0]) + n2)])
                d2Around.append(d2Loop2Right[-(len(xOuter[0]) + n2)])

                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) - n2])
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) - n2])

                if n2 <= elementsAroundHalfOeso + 1:
                    # Annulus right
                    # print('between triple and 6 pt')
                    idx = int(len(xBifurcationRings[0]) * 0.5 + n2 - elementsAroundQuarterOeso - 1)
                    xAround.append(xLoopTripleTo6Pt[idx])
                    if n2 == elementsAroundHalfOeso + 1:
                        d1 = dSmoothLoopTripleTo6Pt[idx]
                        d1Outer[n2][int(len(d1Outer[n2]) * 0.5)] = d1
                    else:
                        d2Around.append(dSmoothLoopTripleTo6Pt[idx])

                    # Annulus left - need to rotate to point towards duodenum
                    xAround.append(xLoopTripleTo6Pt[-idx])
                    # d2Around.append(dSmoothLoopTripleTo6Pt[-idx])
                    d2 = dSmoothLoopTripleTo6Pt[-idx]
                    if n2 < elementsAroundHalfOeso + 1:
                        rotAxis = vector.normalise(
                            vector.crossproduct3(vector.normalise(d1Outer[n2][int(len(d1Outer[n2]) * 0.5 + 1)]),
                                                 vector.normalise(d2)))
                    else: # use d2 on previous overlapping point to rotate
                        rotAxis = vector.normalise(
                            vector.crossproduct3(vector.normalise(d1), vector.normalise(d2)))
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                    d2Around.append(
                        [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in
                         range(3)])

                elif n2 > elementsAroundHalfOeso + 1:
                    # print('beyond 6 pt')
                    # LC
                    xAround.append(xLC[n2 - (elementsAroundHalfOeso + 1)])
                    d2Around.append(d2LC[n2 - (elementsAroundHalfOeso + 1)])

                # Regular down left
                for n1 in range(nextIdx + 1):
                    xAround.append(xRegularLoops[n1][int(len(xRegularLoops[n1]) * 0.5) + n2])
                    d2Around.append(d2RegularLoops[n1][int(len(d2RegularLoops[n1]) * 0.5) + n2])

                # Row 2 left
                xAround.append(xLoop2Left[len(xOuter[0]) + n2 - 1])
                d2Around.append(d2Loop2Left[len(xOuter[0]) + n2 - 1])

            d2Outer.append(d2Around)

        # remove triple point on both sides from downstream ring
        n2Idx = elementsAroundQuarterOeso + 1
        n1Idx = int(len(xOuter[n2Idx]) * 0.5)
        del xOuter[n2Idx][n1Idx: n1Idx + 2], d1Outer[n2Idx][n1Idx: n1Idx + 2], d2Outer[n2Idx][n1Idx: n1Idx + 2]

        d3UnitOuter = []
        for n2 in range(elementsCountAlong + 1):
            d3Around = []
            for n1 in range(len(xOuter[n2])):
                d3Around.append(vector.normalise(
                    vector.crossproduct3(vector.normalise(d1Outer[n2][n1]), vector.normalise(d2Outer[n2][n1]))))
            d3UnitOuter.append(d3Around)

        # for m2 in range(len(xOuter)):
        #     for m1 in range(len(xOuter[m2])):
        #         node = nodes.createNode(nodeIdentifier, nodetemplate)
        #         cache.setNode(node)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[m2][m1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Outer[m2][m1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Outer[m2][m1])
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3UnitOuter[m2][m1])
        #         nodeIdentifier += 1

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
        curvatureAlongGC = findCurvatureAlongLine(xGC, dGC, norms)  # 1st len(xOuter[0]) + 1 are for d1, the rest for d2

        # for m1 in range(len(xGC)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xGC[m1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dGC[m1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1

        # Curvature along rows adjacent to GC - calculate with left and use the same for right
        norms = []
        xTest = []
        for n in range(int(len(xOuter[1]) * 0.5)):  # d1s
            xTest.append(xOuter[1][n + int(len(xOuter[1]) * 0.5) + 1])  # KM
            norms.append(d3UnitOuter[1][n + int(len(xOuter[1]) * 0.5) + 1])
        for n2 in range(2, elementsCountAlong + 1):  # d2s
            xTest.append(xOuter[n2][-1])  # KM
            norms.append(d3UnitOuter[n2][-1])
        curvatureAlong2Left = findCurvatureAlongLine(xLoop2Left, d2Loop2Left, norms)

        # for m1 in range(len(xTest)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTest[m1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, norms[m1])
        #     nodeIdentifier += 1

        # Curvature along LC
        xTest = [] # KM
        norms = []
        for n in range(elementsAroundHalfOeso + 1, elementsCountAlong + 1):
            xTest.append(xOuter[n][int(len(xOuter[n]) * 0.5)])
            norms.append(d3UnitOuter[n][int(len(xOuter[n]) * 0.5)])
        curvatureAlongLC = findCurvatureAlongLine(xLC[1:], d2LC[1:], norms)

        # Curvature along path from triple point to 6 point junction
        xTest = [] # KM
        norms = []
        idxToAnnulus = elementsAroundHalfDuod + 1
        xTest += xOuter[elementsAroundQuarterOeso][: idxToAnnulus]
        norms += d3UnitOuter[elementsAroundQuarterOeso][:idxToAnnulus]

        for n2 in range(elementsAroundQuarterOeso):
            idx = elementsAroundQuarterOeso + 2 + n2
            if idx < elementsAroundHalfOeso + 1:
                xTest.append(xOuter[idx][idxToAnnulus - 1])
                norms.append(d3UnitOuter[idx][idxToAnnulus - 1])
        xTest += xOuter[elementsAroundHalfOeso + 1][idxToAnnulus - 1:] + \
                 xOuter[elementsAroundHalfOeso + 1][: idxToAnnulus]
        norms += d3UnitOuter[elementsAroundHalfOeso + 1][idxToAnnulus - 1:] + \
                 d3UnitOuter[elementsAroundHalfOeso + 1][: idxToAnnulus]

        for n2 in range(elementsAroundQuarterOeso - 1):
            idx = elementsAroundHalfOeso - n2
            xTest.append(xOuter[idx][idxToAnnulus])
            norms.append(d3UnitOuter[idx][idxToAnnulus])

        xTest += xOuter[elementsAroundQuarterOeso][idxToAnnulus:]
        norms += d3UnitOuter[elementsAroundQuarterOeso][idxToAnnulus:]

        curvatureLoopTripleTo6Pt = findCurvatureAroundLoop(xLoopTripleTo6Pt, dSmoothLoopTripleTo6Pt, norms)

        # Curvature along path from GC to triple point
        xTest = []  # KM
        norms = []
        xTest += xOuter[elementsAroundQuarterOeso + 1][:idxToAnnulus - 1]
        norms += d3UnitOuter[elementsAroundQuarterOeso + 1][:idxToAnnulus - 1]

        for n2 in range(elementsAroundQuarterOeso):
            idx = elementsAroundQuarterOeso - n2
            xTest.append(xOuter[idx][int(len(xOuter[idx]) * 0.5)])
            norms.append(d3UnitOuter[idx][int(len(xOuter[idx]) * 0.5)])
        xTest.append(xOuter[0][0])
        norms.append(d3UnitOuter[0][0])
        for n2 in range(1, elementsAroundQuarterOeso + 1):
            xTest.append(xOuter[n2][int(len(xOuter[n2]) * 0.5) + 1])
            norms.append(d3UnitOuter[n2][int(len(xOuter[n2]) * 0.5) + 1])
        xTest += xOuter[elementsAroundQuarterOeso + 1][idxToAnnulus - 1:]
        norms += d3UnitOuter[elementsAroundQuarterOeso + 1][idxToAnnulus - 1:]
        curvatureLoopGCTriplePt = findCurvatureAroundLoop(xLoopGCTriplePt, dSmoothLoopGCTriplePt, norms)

        # Curvature around regular loops
        curvatureRegularLoops = []
        for n1 in range(elementsAroundHalfDuod - 2):
            xTest = []
            norms = []
            for n2 in range(elementsCountAlong):
                idx = -(1 + n2)
                xTest.append(xOuter[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
                norms.append(d3UnitOuter[idx][int(len(xOuter[idx]) * 0.5 - 1 - n1)])
            if n1 < elementsAroundHalfDuod - 3:
                xTest.append(xOuter[0][(n1 + 1)])
                norms.append(d3UnitOuter[0][n1 + 1])
            else:
                xTest.append(xOuter[idx][0])
                norms.append(d3UnitOuter[idx][0])
            for n2 in range(elementsCountAlong):
                xTest.append(xOuter[n2 + 1][int(len(xOuter[n2 + 1]) * 0.5 + n1 + (1 if n2 >= elementsAroundHalfOeso else 2))])
                norms.append(d3UnitOuter[n2 + 1][int(
                    len(xOuter[n2 + 1]) * 0.5 + n1 + (1 if n2 >= elementsAroundHalfOeso else 2))])
            curvatureLoop = findCurvatureAlongLine(xRegularLoops[n1], d2RegularOrderedLoops[n1], norms)
            curvatureRegularLoops.append(curvatureLoop)

        # Assemble curvatures
        d1Curvature = []
        d2Curvature = []
        countUp = 0
        countDown = 0
        for n2 in range(elementsCountAlong + 1):
            d1CurvatureAround = []
            d2CurvatureAround = []
            if n2 == 0:
                # print(n2, 'GC')
                for i in range(elementsAroundHalfDuod - 2):
                    d1CurvatureAround.append(curvatureAlongGC[i])

                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5)])
                for n1 in range(len(xOuter[0]) - 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5)])
                    nextIdx = n1 + 1

            elif n2 == 1:
                # print(n2, 'Row 2')
                d1CurvatureAround.append(curvatureAlongGC[i + n2])
                xTest = []
                dTest = []
                for n in range(int(len(xOuter[1]) * 0.5) - 1, -1, -1):
                    d1CurvatureAround.append(curvatureAlong2Left[n])
                d1CurvatureAround += curvatureAlong2Left[:int(len(xOuter[1]) * 0.5)]
                xTest = xLoop2Left[:int(len(xOuter[1]) * 0.5)]
                dTest = d2Loop2Left[:int(len(xOuter[1]) * 0.5)]

                d2CurvatureAround.append(curvatureRegularLoops[nextIdx][int(len(curvatureRegularLoops[nextIdx]) * 0.5)])
                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                # right point on annulus
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) - n2])
                # left point on annulus
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) + n2])
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])

            elif n2 > 1 and n2 < elementsAroundQuarterOeso + 2:
                # print(n2, 'Before triple point + triple point')
                xAround = xOuter[n2]
                if n2 < elementsAroundQuarterOeso:  # upstream of triple pt
                    # smooth d1 around - Make into function?
                    d1Around = d1Outer[n2]
                    norms = d3UnitOuter[n2]
                    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    normsLoop = norms[int(len(d1Around) * 0.5 + 1):] + norms[: int(len(d1Around) * 0.5 + 1)]
                    curvature = findCurvatureAlongLine(xLoop, d1Loop, normsLoop)
                    # Rearrange to correct order
                    d1CurvatureAround = curvature[int(len(xAround) * 0.5):] + curvature[: int(len(xAround) * 0.5):]

                elif n2 == elementsAroundQuarterOeso:  # upstream bifurcation
                    # take smoothed d1 from dSmoothTripleTo6Pt
                    d1CurvatureAround = curvatureLoopTripleTo6Pt[: int(len(xBifurcationRings[0]) * 0.5) + 1] + \
                                        curvatureLoopTripleTo6Pt[-int(len(xBifurcationRings[0]) * 0.5):]

                elif n2 > elementsAroundQuarterOeso:  # downstream bifurcation
                    # take smoothed d1 from dSmoothGCToTriplePt
                    d1CurvatureAround = curvatureLoopGCTriplePt[: int(len(xBifurcationRings[1]) * 0.5) + 1] + \
                                        curvatureLoopGCTriplePt[-int(len(xBifurcationRings[1]) * 0.5):]

                # GC
                d2CurvatureAround.append(curvatureAlongGC[len(xOuter[0]) + n2 - 1])
                # Row 2 right
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])
                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                # Annulus right
                d2CurvatureAround.append(curvatureLoopGCTriplePt[
                    int(len(curvatureLoopGCTriplePt) * 0.5) - n2 + (1 if n2 > elementsAroundQuarterOeso else 0)])
                # Annulus left
                d2CurvatureAround.append(curvatureLoopGCTriplePt[int(len(curvatureLoopGCTriplePt) * 0.5) + n2 - (
                    1 if n2 > elementsAroundQuarterOeso else 0)])
                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])
                # Row 2 left
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])

            elif n2 > elementsAroundQuarterOeso + 1:
                # print(n2, 'Downstream of triple point')
                xAround = xOuter[n2]
                d1Around = d1Outer[n2]
                normsAround = d3UnitOuter[n2]

                # smooth d1 around - Make into function?
                if n2 < elementsAroundHalfOeso + 1:
                    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    normsLoop = normsAround[int(len(normsAround) * 0.5 + 1):] + normsAround[: int(len(normsAround) * 0.5 + 1)]
                    curvature = findCurvatureAlongLine(xLoop, d1Loop, normsLoop)
                    # Rearrange to correct order
                    d1CurvatureAround = curvature[int(len(xAround) * 0.5):] + curvature[: int(len(xAround) * 0.5):]

                elif n2 == elementsAroundHalfOeso + 1:  # 6 point junction ring
                    # take smoothed d1 from dSmoothedTripleTo6Pt
                    startRightIdx = int(len(xBifurcationRings[0]) * 0.5 + elementsAroundQuarterOeso + len(
                        xAlongAround[junctionIdx]) * 0.5)
                    endRightIdx = startRightIdx + int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    startLeftIdx = startRightIdx - int(len(xAlongAround[junctionIdx]) * 0.5) + 1
                    d1CurvatureAround = curvatureLoopTripleTo6Pt[startRightIdx: endRightIdx] + \
                                        curvatureLoopTripleTo6Pt[startLeftIdx: startRightIdx]

                if n2 > elementsAroundHalfOeso + 1: # closed rings beyond 6 point junction
                    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
                    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
                    normsLoop = normsAround[int(len(normsAround) * 0.5 + 1):] + normsAround[ : int(len(normsAround) * 0.5 + 1)]
                    curvature = findCurvatureAroundLoop(xLoop, d1Loop, normsLoop)
                    # Rearrange to correct order
                    d1CurvatureAround = curvature[int(len(xLoop) * 0.5) - 1:] + curvature[: int(len(xAround) * 0.5) - 1]

                # GC
                d2CurvatureAround.append(curvatureAlongGC[len(xOuter[0]) + n2 - 1])
                # Row 2 right
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])
                # Regular up right
                for n1 in range(nextIdx, -1, -1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) - n2])
                if n2 <= elementsAroundHalfOeso + 1:
                    # Annulus right
                    # print('between triple and 6 pt')
                    idx = int(len(xBifurcationRings[0]) * 0.5 + n2 - elementsAroundQuarterOeso - 1)
                    if n2 == elementsAroundHalfOeso + 1:
                        d1CurvatureAround[int(len(d1Outer[n2]) * 0.5)] = curvatureLoopTripleTo6Pt[idx]
                    else:
                        d2CurvatureAround.append(curvatureLoopTripleTo6Pt[idx])
                    # Annulus left
                    d2CurvatureAround.append(curvatureLoopTripleTo6Pt[-idx])
                elif n2 > elementsAroundHalfOeso + 1:
                    # print('beyond 6 pt')
                    # LC
                    d2CurvatureAround.append(curvatureAlongLC[n2 - (elementsAroundHalfOeso + 1) - 1])
                # Regular down left
                for n1 in range(nextIdx + 1):
                    d2CurvatureAround.append(curvatureRegularLoops[n1][int(len(curvatureRegularLoops[n1]) * 0.5) + n2])
                # Row 2 left
                d2CurvatureAround.append(curvatureAlong2Left[len(xOuter[0]) + n2 - 1])

            d1Curvature.append(d1CurvatureAround)
            d2Curvature.append(d2CurvatureAround)

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
                    factor = 1.0 + wallThickness * (1.0 - xi3) * d1Curvature[n2][n1]
                    d1 = [factor * c for c in d1Outer[n2][n1]]
                    d1List.append(d1)

                    # d2
                    factor = 1.0 + wallThickness * (1.0 - xi3) * d2Curvature[n2][n1]
                    d2 = [factor * c for c in d2Outer[n2][n1]]
                    d2List.append(d2)

                    # d3
                    d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
                    d3List.append(d3)

                    idxAround.append(nodeIdx)
                    nodeIdx += 1

                idxThroughWall.append(idxAround)
            idxMat.append(idxThroughWall)

        if makeStomach:
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

            for e2 in range(elementsAroundQuarterOeso + 2):
                # Row 1
                if e2 == 0:
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
                                    # scaleFactors = [-1.0]
                                    if e1 == 0:
                                        bni11 = startNode + elementsAroundThroughWall + e3 * elementsCountAround2 + e1
                                        bni12 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1 - 1
                                    else:
                                        bni11 = startNode + elementsCountAround1 - e1 + e3 * elementsCountAround1
                                        bni12 = bni11 - 1
                                    bni21 = startNode + elementsAroundThroughWall + 1 + e1 + e3 * elementsCountAround2
                                    bni22 = bni21 + 1
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + (
                                                           elementsCountAround2 if e1 == 0 else elementsCountAround1),
                                                       bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                                    eft1 = eftfactory.createEftNoCrossDerivatives()
                                    scaleFactors = [-1.0]
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scaleEftNodeValueLabels(eft1, [1, 2, 5, 6],
                                                            [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D2_DS1DS2,
                                                             Node.VALUE_LABEL_D2_DS1DS3,
                                                             Node.VALUE_LABEL_D3_DS1DS2DS3], [1])
                                    scaleEftNodeValueLabels(eft1, [1, 2, 5, 6],
                                                            [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
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
                                                       bni12 + (
                                                           elementsCountAround1 if e1 < elementsCountAround1 * 2 else elementsCountAround2),
                                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                                element = mesh.createElement(elementIdentifier, elementtemplate1)
                                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                                if scaleFactors:
                                    result3 = element.setScaleFactors(eft1, scaleFactors)
                                elementIdentifier += 1

                # Row 2
                elif e2 == 1:
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
                                    bni12 = startNode + e3 * elementsCountAround1 + (
                                                e1 + 1)  # % elementsCountAround1
                                    bni21 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + e1
                                    bni22 = startNode + elementsAroundThroughWall + elementsCountAround2 * e3 + (
                                                e1 + 1)  # % elementsCountAround2
                                    if e1 == 0:  # Remap derivatives of element adjacent to GC
                                        scaleFactors = [-1.0]
                                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                           bni11 + elementsCountAround1,
                                                           bni12 + elementsCountAround1,
                                                           bni21 + elementsCountAround2,
                                                           bni22 + elementsCountAround2]
                                        eft1 = eftfactory.createEftNoCrossDerivatives()
                                        setEftScaleFactorIds(eft1, [1], [])
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                    elif e1 == 1:  # Bottom right wedge
                                        nodeIdentifiers = [bni11, bni21, bni22,
                                                           bni11 + elementsCountAround1,
                                                           bni21 + elementsCountAround2,
                                                           bni22 + elementsCountAround2]
                                        eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
                                    elementtemplateX.defineField(coordinates, -1, eft1)
                                    elementtemplate1 = elementtemplateX

                                elif e1 > 1 and e1 < elementsCountAround1:
                                    bni11 = startNode + e3 * elementsCountAround1 + e1 - 1
                                    bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1
                                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                    bni22 = startNode + elementsAroundThroughWall + (
                                                e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                                elif e1 >= elementsCountAround1:
                                    bni11 = startNode + e3 * elementsCountAround1 + e1 - 2
                                    bni12 = startNode + e3 * elementsCountAround1 + (e1 - 1) % elementsCountAround1
                                    bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                    bni22 = startNode + elementsAroundThroughWall + (
                                                e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                                    if e1 == elementsCountAround1:  # Bottom left wedge
                                        nodeIdentifiers = [bni12, bni21, bni22,
                                                           bni12 + elementsCountAround1,
                                                           bni21 + elementsCountAround2,
                                                           bni22 + elementsCountAround2]
                                        eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
                                    elif e1 == elementsCountAround1 + 1:  # Remap derivatives of element adjacent to GC
                                        scaleFactors = [-1.0]
                                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                           bni11 + elementsCountAround1,
                                                           bni12 + elementsCountAround1,
                                                           bni21 + elementsCountAround2,
                                                           bni22 + elementsCountAround2]
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

                # Additional elements between second and upstream bifurcation ring
                elif e2 > 1 and e2 < elementsAroundQuarterOeso:
                    startNode = bodyStartNode
                    for e in range(e2):
                        startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

                    elementsCountAround1 = len(xOuter[e2])
                    elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
                    elementsCountAround2 = len(xOuter[e2 + 1])

                    for e3 in range(elementsCountThroughWall):
                        for e1 in range(elementsCountAround1):
                            if e1 != int(elementsCountAround1 * 0.5):
                                scaleFactors = []
                                eft1 = eftStandard
                                elementtemplate1 = elementtemplateStandard
                                bni11 = startNode + e3 * elementsCountAround1 + e1
                                bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                bni22 = startNode + elementsAroundThroughWall + (
                                        e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                                nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                   bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                   bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                                element = mesh.createElement(elementIdentifier, elementtemplate1)
                                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                                if scaleFactors:
                                    result3 = element.setScaleFactors(eft1, scaleFactors)
                                elementIdentifier += 1

                # Upstream bifurcation
                elif e2 == elementsAroundQuarterOeso:
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
                                bni22 = startNode + elementsAroundThroughWall + (
                                            e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                                if e1 < int(elementsCountAround1 * 0.5) - 1:
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                                elif e1 == int(elementsCountAround1 * 0.5) - 1:  # right wedge
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
                                    eft1 = eftfactory.createEftWedgeCollapseXi2([3, 7])
                                    elementtemplateX.defineField(coordinates, -1, eft1)
                                    elementtemplate1 = elementtemplateX

                                elif e1 > int(elementsCountAround1 * 0.5) + 1:
                                    bni21 = bni21 - 2
                                    bni22 = startNode + elementsAroundThroughWall + (
                                                e1 - 1) % elementsCountAround2 + elementsCountAround2 * e3
                                    nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                                       bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                                       bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                                element = mesh.createElement(elementIdentifier, elementtemplate1)
                                result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                                if scaleFactors:
                                    result3 = element.setScaleFactors(eft1, scaleFactors)
                                elementIdentifier += 1

                # Downstream bifurcation
                elif e2 == elementsAroundQuarterOeso + 1:
                    startNode = bodyStartNode
                    for e in range(e2):
                        startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

                    elementsCountAround1 = len(xOuter[e2])
                    elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
                    elementsCountAround2 = len(xOuter[e2 + 1])
                    for e3 in range(elementsCountThroughWall):
                        for e1 in range(elementsCountAround1 + 1):
                            eft1 = eftStandard
                            elementtemplate1 = elementtemplateStandard
                            if e1 < int(elementsCountAround1 * 0.5) + 1:
                                bni11 = startNode + e3 * elementsCountAround1 + e1
                            elif e1 == int(elementsCountAround1 * 0.5) + 1:
                                bni11 = startNode - len(xOuter[e2 - 1]) * (elementsCountThroughWall + 1) + e3 * len(
                                    xOuter[e2 - 1]) + e1 + 1
                            elif e1 > int(elementsCountAround1 * 0.5) + 1:
                                bni11 = startNode + e3 * elementsCountAround1 + e1 - 1

                            if e1 < int(elementsCountAround1 * 0.5):
                                bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                            elif e1 == int(elementsCountAround1 * 0.5):
                                bni12 = startNode - len(xOuter[e2 - 1]) * (elementsCountThroughWall + 1) + e3 * len(
                                    xOuter[e2 - 1]) + e1 + 1
                            elif e1 > int(elementsCountAround1 * 0.5):
                                bni12 = startNode + e3 * elementsCountAround1 + e1 % elementsCountAround1

                            if e1 > int(elementsCountAround1 * 0.5):
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + 1
                                bni22 = startNode + elementsAroundThroughWall + (
                                            e1 + 2) % elementsCountAround2 + elementsCountAround2 * e3
                            else:
                                bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                                bni22 = startNode + elementsAroundThroughWall + (
                                            e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3

                            nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                               bni11 + (len(xOuter[e2 - 1]) if e1 == int(
                                                   elementsCountAround1 * 0.5) + 1 else elementsCountAround1),
                                               bni12 + (len(xOuter[e2 - 1]) if e1 == int(
                                                   elementsCountAround1 * 0.5) else elementsCountAround1),
                                               bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                            if e1 == int(elementsCountAround1 * 0.5):
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == int(elementsCountAround1 * 0.5) + 1:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [ ])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            element = mesh.createElement(elementIdentifier, elementtemplate1)
                            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                            if scaleFactors:
                                result3 = element.setScaleFactors(eft1, scaleFactors)
                            elementIdentifier += 1

            # Rows between downstream and penultimate ring
            for e2 in range(elementsAroundQuarterOeso + 2, elementsAroundHalfOeso):
                startNode = bodyStartNode
                for e in range(e2):
                    startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)
                elementsCountAround1 = len(xOuter[e2])
                elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
                elementsCountAround2 = len(xOuter[e2 + 1])

                for e3 in range(elementsCountThroughWall):
                    for e1 in range(elementsCountAround1 - 1):
                        bni11 = startNode + e3 * elementsCountAround1 + e1 + (
                            0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                        bni12 = startNode + e3 * elementsCountAround1 + (
                                    e1 + (1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3 + (
                            0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                        bni22 = startNode + elementsAroundThroughWall + (e1 + (1 if e1 < int(
                            elementsCountAround1 * 0.5) else 2)) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]
                        element = mesh.createElement(elementIdentifier, elementtemplateStandard)
                        result = element.setNodesByIdentifier(eftStandard, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1

            # Penultimate row connecting to annulus and beyond
            for e2 in range(elementsAroundHalfOeso, elementsCountAlong):
                startNode = bodyStartNode
                for e in range(e2):
                    startNode += len(xOuter[e]) * (elementsCountThroughWall + 1)

                elementsCountAround1 = len(xOuter[e2])
                elementsAroundThroughWall = elementsCountAround1 * (elementsCountThroughWall + 1)
                elementsCountAround2 = len(xOuter[e2 + 1])

                for e3 in range(elementsCountThroughWall):
                    for e1 in range(
                            elementsCountAround1 - (1 if e2 == elementsAroundHalfOeso else 0)):
                        scaleFactors = []
                        eft1 = eftStandard
                        elementtemplate1 = elementtemplateStandard
                        if e2 == elementsAroundHalfOeso:
                            bni11 = startNode + e3 * elementsCountAround1 + e1 + (
                                0 if e1 < int(elementsCountAround1 * 0.5) else 1)
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + (
                                1 if e1 < int(elementsCountAround1 * 0.5) else 2)) % elementsCountAround1
                            # Remap elements next to annulus
                            if e1 == int(elementsCountAround1 * 0.5) - 1:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX
                        else:
                            bni11 = startNode + e3 * elementsCountAround1 + e1
                            bni12 = startNode + e3 * elementsCountAround1 + (e1 + 1) % elementsCountAround1
                        bni21 = startNode + elementsAroundThroughWall + e1 + elementsCountAround2 * e3
                        bni22 = startNode + elementsAroundThroughWall + (
                                    e1 + 1) % elementsCountAround2 + elementsCountAround2 * e3
                        nodeIdentifiers = [bni11, bni12, bni21, bni22,
                                           bni11 + elementsCountAround1, bni12 + elementsCountAround1,
                                           bni21 + elementsCountAround2, bni22 + elementsCountAround2]

                        if e2 == elementsAroundHalfOeso + 1:
                            if e1 == int(elementsCountAround1 * 0.5) - 1:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == int(elementsCountAround1 * 0.5):
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, ([(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])]))
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        if scaleFactors:
                            result3 = element.setScaleFactors(eft1, scaleFactors)
                        elementIdentifier += 1

            # Annulus
            # Assemble endPoints for annulus
            endPoints_x = [[None] * elementsCountAroundOeso, [None] * elementsCountAroundOeso]
            endPoints_d1 = [[None] * elementsCountAroundOeso, [None] * elementsCountAroundOeso]
            endPoints_d2 = [[None] * elementsCountAroundOeso, [None] * elementsCountAroundOeso]
            endNode_Id = [[None] * elementsCountAroundOeso, [None] * elementsCountAroundOeso]
            endDerivativesMap = [[None] * elementsCountAroundOeso, [None] * elementsCountAroundOeso]
            endProportions = []

            thicknessIdx = [0, -1]
            for nAround in range(elementsCountAroundOeso):
                for n3 in range(len(thicknessIdx)):
                    if nAround == 0:
                        idx = idxMat[nAround][thicknessIdx[n3]][0]
                    elif nAround <= elementsAroundQuarterOeso:
                        idx = idxMat[nAround][thicknessIdx[n3]][int((len(xOuter[nAround]) - 1) * 0.5)]
                    elif elementsAroundQuarterOeso < nAround < elementsAroundHalfOeso:
                        idx = idxMat[nAround + 1][thicknessIdx[n3]][int((len(xOuter[nAround + 1]) - 1) * 0.5)]
                    elif nAround == elementsAroundHalfOeso:
                        idx = idxMat[nAround + 1][thicknessIdx[n3]][int(len(xOuter[nAround + 1]) * 0.5)]
                    elif nAround > elementsAroundHalfOeso:
                        idx = endNode_Id[n3][elementsAroundHalfOeso - (nAround - elementsAroundHalfOeso)] + 1

                    endPoints_x[n3][nAround] = xList[idx - bodyStartNode]
                    endPoints_d1[n3][nAround] = d1List[idx - bodyStartNode]
                    endPoints_d2[n3][nAround] = d2List[idx - bodyStartNode]
                    endNode_Id[n3][nAround] = idx

                    if n3 == len(thicknessIdx) - 1: # outer layer
                        endPosition = trackSurfaceStomach.findNearestPosition(endPoints_x[n3][nAround])
                        # xCheck = trackSurfaceStomach.evaluateCoordinates(endPosition) # KM
                        # print(xCheck, endPoints_x[n3][nAround]) # KM
                        endProportions.append(trackSurfaceStomach.getProportion(endPosition))

            for nAround in range(elementsCountAroundOeso):
                if nAround == 0:
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)
                elif nAround == elementsAroundQuarterOeso:
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
                elif 0 < nAround < elementsAroundHalfOeso:
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, 1, 0), (-1, 0, 0), None)
                elif nAround == elementsAroundHalfOeso:
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
                elif nAround == int(elementsCountAroundOeso * 0.75):
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
                elif elementsAroundHalfOeso < nAround < elementsCountAroundOeso:
                    endDerivativesMap[0][nAround] = endDerivativesMap[1][nAround] = ((0, -1, 0), (1, 0, 0), None)

            startProportions = []
            for n in range(elementsCountAroundOeso):
                startProportions.append(trackSurfaceStomach.getProportion(o1_Positions[n]))

            nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
                nodes, mesh, nodeIdentifier, elementIdentifier,
                o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
                endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap,
                elementsCountRadial = elementsCountAnnulus, tracksurface=trackSurfaceStomach,
                startProportions = startProportions, endProportions = endProportions,
                rescaleStartDerivatives = True, rescaleEndDerivatives = True)

        if showTrackSurface:
            for n2 in range(len(xTrackSurface)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1TrackSurface[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2TrackSurface[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                nodeIdentifier += 1

        if showCentralPath:
            for n2 in range(len(sx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sd2[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sd1[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3[n2])
                nodeIdentifier += 1

        # nodeIdentifier = 10000
        # for n1 in range(len(xStarts)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xStarts[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2Starts[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1
        #
        # for n1 in range(len(xEnds)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xEnds[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2Ends[n1])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier += 1


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
                                           startDerivativeMagnitude = None, endDerivativeMagnitude = None, curveMode = 1):
    """
    EDIT
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
                                              elementsOut, startDerivative, endDerivative, curveMode)


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

