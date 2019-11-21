"""
Generates a 3-D bladder mesh with variable numbers of elements around and up.
"""

from __future__ import division
import math
import copy
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import zinc_utils
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_bladder1(Scaffold_base):
    '''
    3-D bladder scaffold.
    '''

    ostiumDefaultScaffoldPackages = {
        'Ostium Cat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 0.3,
                'Ostium length': 0.2,
                'Ostium wall thickness': 0.05,
                'Ostium inter-vessel distance': 0.8,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': False,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.15,
                'Vessel wall thickness': 0.04,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                # 'Use cross derivatives' : False,
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
                'Number of elements around ostium': 8,
                'Number of elements along': 1,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 0.15,
                'Ostium length': 0.05,
                'Ostium wall thickness': 0.02,
                'Ostium inter-vessel distance': 0.8,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': False,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.04,
                'Vessel wall thickness': 0.01,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                # 'Use cross derivatives' : False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
                },
            })
        }

    @staticmethod
    def getName():
        return '3D Bladder 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cat 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Rat' in parameterSetName:
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Ostium Rat 1']
        else:
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Ostium Cat 1']

        options = {
            'Number of elements up 1': 8,
            'Number of elements up 2': 16,
            'Number of elements around': 8,  # should be even
            'Number of elements through wall': 1,
            'Height': 5.0,
            'Major diameter': 3.0,
            'Minor diameter': 3.0,
            'Urethra radius': 0.5,
            'Bladder wall thickness': 0.05,
            'Number of elements around ostium': 8,# implemented for 8
            'Ureter': copy.deepcopy(ostiumOption),
            'Ostium position around': 0.15,
            'Ostium position up': 0.25,
            'Number of elements radially on annulus': 1,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements up': 4,
            'Refine number of elements through wall': 1
        }

        if 'Rat' in parameterSetName:
            options = {
                'Ureter': copy.deepcopy(ostiumOption),
                'Number of elements up 1': 8,
                'Number of elements up 2': 16,
                'Number of elements around': 16,  # should be even
                'Number of elements through wall': 1,
                'Height': 3.0,
                'Major diameter': 2.5,
                'Minor diameter': 1.5,
                'Urethra radius': 0.35,
                'Bladder wall thickness': 0.05,
                'Number of elements around ostium': 8,
                'Ostium position around': 0.55,
                'Ostium position up': 0.65,
                'Number of elements radially on annulus': 2,
                'Use cross derivatives': False
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up 1',
            'Number of elements up 2',
            'Number of elements around',
            'Number of elements through wall',
            'Height',
            'Major diameter',
            'Minor diameter',
            'Urethra radius',
            'Bladder wall thickness',
            'Ostium position around',
            'Ostium position up',
            'Ureter',
            'Number of elements radially on annulus',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements through wall'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Ureter':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
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
        if optionName == 'Ureter':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Ureter'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Ureter'):
            options['Ureter'] = cls.getOptionScaffoldPackage('Ureter', MeshType_3d_ostium1)
        if (options['Number of elements up 1'] < 4):
            options['Number of elements up 1'] = 4
        if (options['Number of elements up 2'] < 4):
            options['Number of elements up 2'] = 4
        if (options['Number of elements around'] < 8):
            options['Number of elements around'] = 8
        if (options['Number of elements through wall'] != 1):
            options['Number of elements through wall'] = 1
        if (options['Major diameter'] < options['Minor diameter']):
            options['Major diameter'] = options['Minor diameter']

    @staticmethod
    def generateBaseMesh(region, options):
        '''
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        '''

        elementsCountUp1 = options['Number of elements up 1']
        elementsCountUp2 = options['Number of elements up 2']
        elementsCountAround = options['Number of elements around']
        height = options['Height']
        print('height =', height)
        majorDiameter = options['Major diameter']
        minorDiameter = options['Minor diameter']
        radius = options['Urethra radius']
        bladderWallThickness = options['Bladder wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        elementsCountAroundOstium = options['Number of elements around ostium']
        elementsCountAnnulusRadially = options['Number of elements radially on annulus']
        ostiumPositionAround = options['Ostium position around']
        ostiumPositionUp = options['Ostium position up']

        ostiumOptions = options['Ureter']
        ostiumDefaultOptions = ostiumOptions.getScaffoldSettings()

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = zinc_utils.getOrCreateCoordinateField(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        else:
            nodetemplate = nodetemplateApex

        mesh = fm.findMeshByDimension(3)
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        # create lower part of the bladder
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementUp1 = (math.pi/4)/elementsCountUp1
        zero = [0.0, 0.0, 0.0]

        # create lower half of the ellipsoidal
        neckHeight = height - height * math.cos(math.pi / 4)
        listTotalLower_x = []
        listTotalLower_d1 = []
        listTotalLower_d2 = []
        for n3 in range(2):
            inner_x = []
            inner_d1 = []
            inner_d2 = []
            outer_x = []
            outer_d1 = []
            outer_d2 = []
            for n2 in range(0, elementsCountUp1+1):
                radiansUp = n2 * radiansPerElementUp1
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                MajorDiameter = majorDiameter * sinRadiansUp - n3 * bladderWallThickness
                MinorDiameter = minorDiameter * sinRadiansUp - n3 * bladderWallThickness
                if n2 == 0:
                    for n1 in range(elementsCountAround):
                        radiansAround = n1 * radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            -MajorDiameter * sinRadiansAround,
                            MinorDiameter * cosRadiansAround,
                            -height - neckHeight
                        ]
                        dx_ds1 = [
                            -MajorDiameter * cosRadiansAround * radiansPerElementAround,
                            MinorDiameter * -sinRadiansAround * radiansPerElementAround,
                            0.0
                        ]
                        dx_ds2 = [
                            -majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUp1,
                            minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUp1,
                            height * sinRadiansUp * radiansPerElementUp1
                        ]
                        outer_x.append(x)
                        outer_d1.append(dx_ds1)
                        outer_d2.append(dx_ds2)
                else:
                    for n1 in range(elementsCountAround):
                        neckHeight = height - height * math.cos(math.pi/4)
                        radiansAround = n1 * radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            -MajorDiameter * sinRadiansAround,
                            MinorDiameter * cosRadiansAround,
                            -height - neckHeight + n2 * 2 * neckHeight / elementsCountUp1
                        ]
                        dx_ds1 = [
                            -MajorDiameter * cosRadiansAround * radiansPerElementAround,
                            MinorDiameter * -sinRadiansAround * radiansPerElementAround,
                            0.0
                        ]
                        dx_ds2 = [
                            -majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUp1,
                            minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUp1,
                            height * sinRadiansUp * radiansPerElementUp1
                        ]
                        outer_x.append(x)
                        outer_d1.append(dx_ds1)
                        outer_d2.append(dx_ds2)

            # create tube nodes
            radiansPerElementAround = 2.0 * math.pi / elementsCountAround
            for n2 in range(0, elementsCountUp1 + 1):
                radiansUp = n2 * radiansPerElementUp1
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                Radius = radius - n3 * bladderWallThickness
                if n2 == 0:
                    for n1 in range(elementsCountAround):
                        radiansAround = n1 * radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            -Radius * sinRadiansAround,
                            Radius * cosRadiansAround,
                            -height - neckHeight
                        ]
                        dx_ds1 = [
                            -radiansPerElementAround * Radius * cosRadiansAround,
                            radiansPerElementAround * Radius * -sinRadiansAround,
                            0.0
                        ]
                        dx_ds2 = [0, 0, height / elementsCountUp1]
                        inner_x.append(x)
                        inner_d1.append(dx_ds1)
                        inner_d2.append(dx_ds2)
                else:
                    for n1 in range(elementsCountAround):
                        neckHeight = height - height* math.cos(math.pi/4)
                        radiansAround = n1 * radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            -Radius * sinRadiansAround,
                            Radius * cosRadiansAround,
                            -height - neckHeight + n2 * 2 * neckHeight / elementsCountUp1
                        ]
                        dx_ds1 = [
                            -radiansPerElementAround * Radius * cosRadiansAround,
                            radiansPerElementAround * Radius * -sinRadiansAround,
                            0.0
                        ]
                        dx_ds2 = [0, 0, height / elementsCountUp1]
                        inner_x.append(x)
                        inner_d1.append(dx_ds1)
                        inner_d2.append(dx_ds2)

            # interpolation between the lower half of the ellipsoidal and the tube
            m1 = 0
            z_bottom = outer_x[-1][2]
            z_top = outer_x[0][2]
            delta_z = z_top - z_bottom
            interpolatedNodes = []
            interpolatedNodes_d1 = []
            interpolatedNodes_d2 = []
            for n2 in range(elementsCountUp1+1):
                xi = 1.0 - (outer_x[m1][2] - z_bottom) / delta_z
                for n1 in range(elementsCountAround):
                    phi_inner, _, phi_outer, _ = getCubicHermiteBasis(xi)
                    x = [(phi_inner*inner_x[m1][c] + phi_outer*outer_x[m1][c]) for c in range(3)]
                    d1 = [(phi_inner*inner_d1[m1][c] + phi_outer*outer_d1[m1][c]) for c in range(3)]
                    d2 = [(phi_inner*inner_d2[m1][c] + phi_outer*outer_d2[m1][c]) for c in range(3)]

                    interpolatedNodes.append(x)
                    interpolatedNodes_d1.append(d1)
                    interpolatedNodes_d2.append(d2)
                    m1 += 1

            # smoothing the derivatives
            sd2Raw = []
            for n1 in range(elementsCountAround):
                lineSmoothingNodes = []
                lineSmoothingNodes_d2 = []
                for n2 in range(elementsCountUp1+1):
                        lineSmoothingNodes.append(interpolatedNodes[n1 + n2 * elementsCountAround])
                        lineSmoothingNodes_d2.append(interpolatedNodes_d2[n1 + n2 * elementsCountAround])
                sd2 = smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2,
                                                        fixAllDirections=False,
                                                        fixStartDerivative=True, fixEndDerivative=True,
                                                        fixStartDirection=False, fixEndDirection=False)
                sd2Raw.append(sd2)

            # rearrange the derivatives order
            d2RearrangedList = []
            for n2 in range(elementsCountUp1+1):
                for n1 in range(elementsCountAround):
                    d2 = sd2Raw[n1][n2]
                    d2RearrangedList.append(d2)

            # create tracksurface at the lower part of the bladder
            if n3 == 0:
                nodesOnTrackSurface = []
                nodesOnTrackSurface_d1 = []
                nodesOnTrackSurface_d2 = []
                for n2 in range(elementsCountUp1+1):
                    for n1 in range(elementsCountAround):
                        if (n1 <= elementsCountAround / 2):
                            nodesOnTrackSurface.append(interpolatedNodes[n2 * elementsCountAround + n1])
                            nodesOnTrackSurface_d1.append(interpolatedNodes_d1[n2 * elementsCountAround + n1])
                            nodesOnTrackSurface_d2.append(d2RearrangedList[n2 * elementsCountAround + n1])

            # set nodes and derivatives of the lower part of the bladder
            elementsCount1 = elementsCountAround // 2
            elementsCount2 = elementsCountUp1
            tracksurfaceOstium1 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface, nodesOnTrackSurface_d1,
                                        nodesOnTrackSurface_d2)
            ostium1Position = tracksurfaceOstium1.createPositionProportion(ostiumPositionAround, ostiumPositionUp)
            ostium1Position.xi1 = 1.0
            ostium1Position.xi2 = 1.0
            ostiumElementPositionAround = ostium1Position.e1
            ostiumElementPositionUp = ostium1Position.e2
            for n2 in range(len(interpolatedNodes)):
                if n2 == (ostiumElementPositionUp + 1) * elementsCountAround + ostiumElementPositionAround + 1:
                    pass
                elif n2 == (ostiumElementPositionUp + 1) * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 1:
                    pass
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, interpolatedNodes[n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, interpolatedNodes_d1[n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2RearrangedList[n2])
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    listTotalLower_x.append(interpolatedNodes[n2])
                    listTotalLower_d1.append(interpolatedNodes_d1[n2])
                    listTotalLower_d2.append(d2RearrangedList[n2])
                nodeIdentifier += 1

            # create upper part of the bladder
            radiansPerElementAround = 2.0 * math.pi / elementsCountAround
            radiansPerElementUp2 = (3 * math.pi / 4) / elementsCountUp2
            # create regular rows
            for n2 in range(1, elementsCountUp2):
                radiansUp = (math.pi / 4) + n2 * radiansPerElementUp2
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                MajorDiameter = majorDiameter * sinRadiansUp - n3 * bladderWallThickness
                MinorDiameter = minorDiameter * sinRadiansUp - n3 * bladderWallThickness
                for n1 in range(elementsCountAround):
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    et = (2 * height + neckHeight) / (elementsCountUp1 + elementsCountUp2)
                    x = [
                        -MajorDiameter * sinRadiansAround,
                        MinorDiameter * cosRadiansAround,
                        -height * cosRadiansUp
                    ]
                    dx_ds1 = [
                        -MajorDiameter * cosRadiansAround * radiansPerElementAround,
                        MinorDiameter * -sinRadiansAround * radiansPerElementAround,
                        0.0
                    ]
                    dx_ds2 = [
                        -majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUp2,
                        minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUp2,
                        height*sinRadiansUp * radiansPerElementUp2
                    ]
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    nodeIdentifier += 1

            # create apex node
            node = nodes.createNode(nodeIdentifier, nodetemplateApex)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [0.0, 0.0, height - n3 * bladderWallThickness])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [height*radiansPerElementUp2/2, 0.0, 0.0])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, height*radiansPerElementUp2/2, 0.0])
            nodeIdentifier += 1

        # create ureters on the surface
        elementIdentifier = 1

        # ureter 1
        centerUreter1_x, centerUreter1_d1, centerUreter1_d2 = tracksurfaceOstium1.evaluateCoordinates(ostium1Position, derivatives=True)
        td1, td2, td3 = calculate_surface_axes(centerUreter1_d1, centerUreter1_d2, [1.0, 0.0, 0.0])
        m1 = ostiumElementPositionUp * elementsCountAround + ostiumElementPositionAround
        ureter1StartCornerx = listTotalLower_x[m1]
        v1 = [(ureter1StartCornerx[c] - centerUreter1_x[c]) for c in range(3)]
        ostium1Direction = vector.crossproduct3(td3, v1)

        nodeIdentifier, elementIdentifier, (o1_x, o1_d1, o1_d2, _, o1_NodeId, o1_Positions) = \
            generateOstiumMesh(region, ostiumDefaultOptions, tracksurfaceOstium1, ostium1Position, ostium1Direction,
                            startNodeIdentifier=nodeIdentifier,
                            startElementIdentifier=elementIdentifier)

        # ureter 2
        tracksurfaceOstium2 = tracksurfaceOstium1.createMirrorX()
        ostium2Position = TrackSurfacePosition(elementsCountAround - ostiumElementPositionAround, ostiumElementPositionUp - 1, 0.0, 1.0)
        centerUreter2_x, centerUreter2_d1, centerUreter2_d2 = tracksurfaceOstium2.evaluateCoordinates(ostium2Position, derivatives = True)
        ad1, ad2, ad3 = calculate_surface_axes(centerUreter2_d1, centerUreter2_d2, [1.0, 0.0, 0.0])
        if elementsCountAroundOstium == 4:
            m2 = ostiumElementPositionUp * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 1
        else:
            m2 = ostiumElementPositionUp * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 2
        ureter2StartCornerx = listTotalLower_x[m2]
        v2 = [(ureter2StartCornerx[c] - centerUreter2_x[c]) for c in range(3)]
        ostium2Direction = vector.crossproduct3(ad3, v2)

        nodeIdentifier, elementIdentifier, (o2_x, o2_d1, o2_d2, _, o2_NodeId, o2_Positions) = \
            generateOstiumMesh(region, ostiumDefaultOptions, tracksurfaceOstium2, ostium2Position, ostium2Direction,
                           startNodeIdentifier=nodeIdentifier,
                           startElementIdentifier=elementIdentifier)

        # Create annulus mesh around ostium
        endPoints1_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endNode1_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endDerivativesMap = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endNode2_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]

        nodeCountsEachWallLayer = (elementsCountUp1 + elementsCountUp2) * elementsCountAround + 1
        nodeCountsUpperPart = elementsCountAround * (elementsCountUp2 - 1) + 1

        for n3 in range(2):
            n1 = 0
            endNode1_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround) + ostiumElementPositionAround + 1
            endNode1_Id[n3][n1 + 1] = endNode1_Id[n3][n1] + elementsCountAround
            endNode1_Id[n3][n1 + 2] = endNode1_Id[n3][n1 + 1] + elementsCountAround
            endNode1_Id[n3][n1 + 3] = endNode1_Id[n3][n1 + 2] + 1
            endNode1_Id[n3][n1 + 4] = endNode1_Id[n3][n1 + 3] + 1
            endNode1_Id[n3][n1 + 5] = endNode1_Id[n3][n1 + 1] + 2
            endNode1_Id[n3][n1 + 6] = endNode1_Id[n3][n1] + 2
            endNode1_Id[n3][n1 + 7] = endNode1_Id[n3][n1] + 1
            if ostiumElementPositionAround == 0:
                endNode2_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround)\
                                        + elementsCountAround - ostiumElementPositionAround - 1
                endNode2_Id[n3][n1 + 1] = endNode2_Id[n3][n1] + elementsCountAround
                endNode2_Id[n3][n1 + 2] = endNode2_Id[n3][n1 + 1] + elementsCountAround
                endNode2_Id[n3][n1 + 3] = endNode2_Id[n3][n1 + 2] + 1
                endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] - elementsCountAround + 1
                endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 4] - elementsCountAround
                endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1 + 5] - elementsCountAround
                endNode2_Id[n3][n1 + 7] = endNode2_Id[n3][n1] + 1
            else:
                endNode2_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround)\
                                        + elementsCountAround - ostiumElementPositionAround - 1
                endNode2_Id[n3][n1 + 1] = endNode2_Id[n3][n1] + elementsCountAround
                endNode2_Id[n3][n1 + 2] = endNode2_Id[n3][n1 + 1] + elementsCountAround
                endNode2_Id[n3][n1 + 3] = endNode2_Id[n3][n1 + 2] + 1
                endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] + 1
                endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 1] + 2
                endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1] + 2
                endNode2_Id[n3][n1 + 7] = endNode2_Id[n3][n1] + 1

        for n3 in range(2):
            for n1 in range(elementsCountAroundOstium):
                if 1 < n1 < 5:
                    nc1 = endNode1_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 - 2 * (1 - n3)
                elif n1 == 5:
                    nc1 = endNode1_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 1 - 2 * (1 - n3)
                else:
                    nc1 = endNode1_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 * (1 - n3)
                endPoints1_x[n3][n1] = listTotalLower_x[nc1]
                endPoints1_d1[n3][n1] = listTotalLower_d1[nc1]
                endPoints1_d2[n3][n1] = [listTotalLower_d2[nc1][c] for c in range (3)]
                if n1 == 1:
                    nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 1 - 2 * (1 - n3)
                elif 1 < n1 < 5:
                    nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 - 2 * (1 - n3)
                elif n1 == 5:
                    if ostiumElementPositionAround == 0:
                        nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 * (1 - n3)
                    else:
                        nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 - 2 * (1 - n3)
                else:
                    nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsUpperPart - 1 - 2 * (1 - n3)
                endPoints2_x[n3][n1] = listTotalLower_x[nc2]
                endPoints2_d1[n3][n1] = listTotalLower_d1[nc2]
                endPoints2_d2[n3][n1] = listTotalLower_d2[nc2]

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
            elementsCountRadial= elementsCountAnnulusRadially)

        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            o2_x, o2_d1, o2_d2, None, o2_NodeId, None,
            endPoints2_x, endPoints2_d1, endPoints2_d2, None, endNode2_Id, endDerivativesMap,
            elementsCountRadial= elementsCountAnnulusRadially)

        # create elements
        for e3 in range(1):
            ewl = elementsCountEachWallLayer = (e3 + 1) * ((elementsCountUp1 + elementsCountUp2) * elementsCountAround + 1)
            # create lower part of the bladder
            for e2 in range(elementsCountUp1):
                for e1 in range(elementsCountAround):
                    if elementsCountAroundOstium == 4:
                        if e2 == ostiumElementPositionUp and (e1 == ostiumElementPositionAround or e1 == elementsCountAround - ostiumElementPositionAround - 1):
                            pass
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            bni1 = e2 * elementsCountAround + e1 + 1
                            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni2 + elementsCountAround
                            nodeIdentifiers = [bni1 + ewl, bni2 + ewl, bni3 + ewl, bni4 + ewl,
                                               bni1, bni2, bni3, bni4]
                            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                            elementIdentifier += 1
                    if elementsCountAroundOstium == 6:
                        if e2 == ostiumElementPositionUp and (e1 == ostiumElementPositionAround or e1 == ostiumElementPositionAround + 1):
                            pass
                        elif e2 == ostiumElementPositionUp and (e1 == elementsCountAround - ostiumElementPositionAround - 2 or e1 == elementsCountAround - 1 - ostiumElementPositionAround):
                            pass
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            bni1 = e2 * elementsCountAround + e1 + 1
                            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni2 + elementsCountAround
                            nodeIdentifiers = [bni1 + ewl, bni2 + ewl, bni3 + ewl, bni4 + ewl,
                                               bni1, bni2, bni3, bni4]
                            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                            elementIdentifier += 1
                    if elementsCountAroundOstium == 8:
                        if e2 == ostiumElementPositionUp and (e1 == ostiumElementPositionAround or e1 == ostiumElementPositionAround + 1):
                            pass
                        elif e2 == ostiumElementPositionUp + 1 and (e1 == ostiumElementPositionAround or e1 == ostiumElementPositionAround + 1):
                            pass
                        elif e2 == ostiumElementPositionUp and (e1 == elementsCountAround - ostiumElementPositionAround - 2 or e1 == elementsCountAround - 1 - ostiumElementPositionAround):
                            pass
                        elif e2 == ostiumElementPositionUp + 1 and (e1 == elementsCountAround - ostiumElementPositionAround - 2 or e1 == elementsCountAround - 1 - ostiumElementPositionAround):
                            pass
                        else:
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            bni1 = e2 * elementsCountAround + e1 + 1
                            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni2 + elementsCountAround
                            nodeIdentifiers = [bni1 + ewl, bni2 + ewl, bni3 + ewl, bni4 + ewl,
                                               bni1, bni2, bni3, bni4]
                            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                            elementIdentifier += 1

            # create upper part of the bladder
            for e2 in range(elementsCountUp1, (elementsCountUp1 + elementsCountUp2 - 1)):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni1 = e2 * elementsCountAround + e1 + 1
                    bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                    bni3 = bni1 + elementsCountAround
                    bni4 = bni2 + elementsCountAround
                    nodeIdentifiers = [bni1 + ewl, bni2 + ewl, bni3 + ewl, bni4 + ewl,
                                       bni1, bni2, bni3, bni4]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier += 1

            # create Apex elements
            bni3 = 1 + (elementsCountUp1 + elementsCountUp2) * elementsCountAround
            elementtemplateApex = mesh.createElementtemplate()
            elementtemplateApex.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1) % elementsCountAround
                eftApex = eftfactory.createEftShellPoleTop(va, vb)
                elementtemplateApex.defineField(coordinates, -1, eftApex)

                # redefine field in template for changes to eftApex:
                element = mesh.createElement(elementIdentifier, elementtemplateApex)
                bni1 = bni3 - elementsCountAround + e1
                bni2 = bni3 - elementsCountAround + (e1 + 1) % elementsCountAround
                nodeIdentifiers = [bni1 + ewl, bni2 + ewl, bni3 + ewl, bni1, bni2, bni3]
                result = element.setNodesByIdentifier(eftApex, nodeIdentifiers)
                elementIdentifier += 1

        fm.endChange()

    @classmethod
    def generateMesh(cls, region, options):
        '''
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        '''
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountThroughWall)
