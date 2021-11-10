"""
Generates a 3-D bladder mesh with variable numbers of elements around and up.
"""

from __future__ import division

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes


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
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.15,
                'Vessel wall thickness': 0.04,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                # 'Use cross derivatives': False,
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
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.04,
                'Vessel wall thickness': 0.01,
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                # 'Use cross derivatives': False,
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
            'Number of elements up neck': 6,
            'Number of elements up body': 8,
            'Number of elements around': 8,  # should be even
            'Number of elements through wall': 1,
            'Number of elements around ostium': 8,  # implemented for 8
            'Number of elements radially on annulus': 1,
            'Height': 5.0,
            'Major diameter': 6.0,
            'Minor diameter': 6.0,
            'Bladder wall thickness': 0.05,
            'Urethra diameter': 1.0,
            'Ureter': copy.deepcopy(ostiumOption),
            'Ostium position around': 0.55,  # should be on the dorsal part
            'Ostium position up': 0.55,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements up': 4,
            'Refine number of elements through wall': 1
        }
        if 'Rat' in parameterSetName:
            options['Number of elements around'] = 16  # should be even
            options['Height'] = 3.0
            options['Major diameter'] = 5.0
            options['Minor diameter'] = 3.0
            options['Bladder wall thickness'] = 0.02
            options['Urethra diameter'] = 0.7
            options['Ostium position around'] = 0.65  # should be on the dorsal part
            options['Ostium position up'] = 0.75
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up neck',
            'Number of elements up body',
            'Number of elements around',
            'Number of elements through wall',
            'Number of elements radially on annulus',
            'Height',
            'Major diameter',
            'Minor diameter',
            'Bladder wall thickness',
            'Urethra diameter',
            'Ostium position around',
            'Ostium position up',
            'Ureter',
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
        if (options['Number of elements up neck'] < 4):
            options['Number of elements up neck'] = 4
        if (options['Number of elements up body'] < 4):
            options['Number of elements up body'] = 4
        if (options['Number of elements around'] < 8):
            options['Number of elements around'] = 8
        elif (options['Number of elements around'] % 2) == 1:
            options['Number of elements around'] += 1
        if (options['Number of elements through wall'] != 1):
            options['Number of elements through wall'] = 1
        if(options['Number of elements radially on annulus'] > 2):
            options['Number of elements radially on annulus'] = 2
        if (options['Major diameter'] < options['Minor diameter']):
            options['Major diameter'] = options['Minor diameter']

    @classmethod
    def generateBaseMesh(cls, region, options):
        '''
        Generate the base bicubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        '''
        elementsCountUpNeck = options['Number of elements up neck']
        elementsCountUpBody = options['Number of elements up body']
        elementsCountAround = options['Number of elements around']
        height = options['Height']
        majorDiameter = options['Major diameter']
        minorDiameter = options['Minor diameter']
        radius = 0.5 * options['Urethra diameter']
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
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        mesh = fm.findMeshByDimension(3)

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

        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        neckGroup = AnnotationGroup(region, get_bladder_term("neck of urinary bladder"))
        bodyGroup = AnnotationGroup(region, get_bladder_term("dome of the bladder"))
        urinaryBladderGroup = AnnotationGroup(region, get_bladder_term("urinary bladder"))
        annotationGroups = [neckGroup, bodyGroup, urinaryBladderGroup]

        neckMeshGroup = neckGroup.getMeshGroup(mesh)
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        urinaryBladderMeshGroup = urinaryBladderGroup.getMeshGroup(mesh)

        # create nodes
        # create neck of the bladder
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementUpNeck = (math.pi/4)/elementsCountUpNeck

        # create lower part of the ellipsoidal
        neckHeight = height - height * math.cos(math.pi / 4)
        ellipsoidal_x = []
        ellipsoidal_d1 = []
        ellipsoidal_d2 = []
        for n2 in range(0, elementsCountUpNeck+1):
            radiansUp = n2 * radiansPerElementUpNeck
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            majorRadius = 0.5 * majorDiameter * sinRadiansUp
            minorRadius = 0.5 * minorDiameter * sinRadiansUp
            if n2 == 0:
                for n1 in range(elementsCountAround):
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        -majorRadius * sinRadiansAround,
                        minorRadius * cosRadiansAround,
                        -height - neckHeight
                    ]
                    dx_ds1 = [
                        -majorRadius * cosRadiansAround * radiansPerElementAround,
                        minorRadius * -sinRadiansAround * radiansPerElementAround,
                        0.0
                    ]
                    dx_ds2 = [
                        -0.5 * majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUpNeck,
                        0.5 * minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUpNeck,
                        height * sinRadiansUp * radiansPerElementUpNeck
                    ]
                    ellipsoidal_x.append(x)
                    ellipsoidal_d1.append(dx_ds1)
                    ellipsoidal_d2.append(dx_ds2)
            else:
                for n1 in range(elementsCountAround):
                    neckHeight = height - height * math.cos(math.pi/4)
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        -majorRadius * sinRadiansAround,
                        minorRadius * cosRadiansAround,
                        -height - neckHeight + n2 * 2 * neckHeight / elementsCountUpNeck
                    ]
                    dx_ds1 = [
                        -majorRadius * cosRadiansAround * radiansPerElementAround,
                        minorRadius * -sinRadiansAround * radiansPerElementAround,
                        0.0
                    ]
                    dx_ds2 = [
                        -0.5 * majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUpNeck,
                        0.5 * minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUpNeck,
                        height * sinRadiansUp * radiansPerElementUpNeck
                    ]
                    ellipsoidal_x.append(x)
                    ellipsoidal_d1.append(dx_ds1)
                    ellipsoidal_d2.append(dx_ds2)

        # create tube nodes
        radiansPerElementAround = 2.0 * math.pi / elementsCountAround
        tube_x = []
        tube_d1 = []
        tube_d2 = []
        for n2 in range(0, elementsCountUpNeck + 1):
            radiansUp = n2 * radiansPerElementUpNeck
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            if n2 == 0:
                for n1 in range(elementsCountAround):
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        -radius * sinRadiansAround,
                        radius * cosRadiansAround,
                        -height - neckHeight
                    ]
                    dx_ds1 = [
                        -radiansPerElementAround * radius * cosRadiansAround,
                        radiansPerElementAround * radius * -sinRadiansAround,
                        0.0
                    ]
                    dx_ds2 = [0, 0, height / (2 * elementsCountUpNeck)]
                    tube_x.append(x)
                    tube_d1.append(dx_ds1)
                    tube_d2.append(dx_ds2)
            else:
                for n1 in range(elementsCountAround):
                    neckHeight = height - height* math.cos(math.pi/4)
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [
                        -radius * sinRadiansAround,
                        radius * cosRadiansAround,
                        -height - neckHeight + n2 * 2 * neckHeight / elementsCountUpNeck
                    ]
                    dx_ds1 = [
                        -radiansPerElementAround * radius * cosRadiansAround,
                        radiansPerElementAround * radius * -sinRadiansAround,
                        0.0
                    ]
                    dx_ds2 = [0, 0, height / elementsCountUpNeck]
                    tube_x.append(x)
                    tube_d1.append(dx_ds1)
                    tube_d2.append(dx_ds2)

        # interpolation between the lower part of the ellipsoidal and the tube
        m1 = 0
        z_bottom = ellipsoidal_x[-1][2]
        z_top = ellipsoidal_x[0][2]
        delta_z = z_top - z_bottom
        interpolatedNodes = []
        interpolatedNodes_d1 = []
        interpolatedNodes_d2 = []
        for n2 in range(elementsCountUpNeck+1):
            xi = 1.0 - (ellipsoidal_x[m1][2] - z_bottom) / delta_z
            for n1 in range(elementsCountAround):
                phi_inner, _, phi_outer, _ = getCubicHermiteBasis(xi)
                x = [(phi_inner*tube_x[m1][c] + phi_outer*ellipsoidal_x[m1][c]) for c in range(3)]
                d1 = [(phi_inner*tube_d1[m1][c] + phi_outer*ellipsoidal_d1[m1][c]) for c in range(3)]
                d2 = [(phi_inner*tube_d2[m1][c] + phi_outer*ellipsoidal_d2[m1][c]) for c in range(3)]
                interpolatedNodes.append(x)
                interpolatedNodes_d1.append(d1)
                interpolatedNodes_d2.append(d2)
                m1 += 1

        # smoothing the derivatives
        sd2Raw = []
        for n1 in range(elementsCountAround):
            lineSmoothingNodes = []
            lineSmoothingNodes_d2 = []
            for n2 in range(elementsCountUpNeck+1):
                    lineSmoothingNodes.append(interpolatedNodes[n1 + n2 * elementsCountAround])
                    lineSmoothingNodes_d2.append(interpolatedNodes_d2[n1 + n2 * elementsCountAround])
            sd2 = smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2,
                                                    fixAllDirections=False,
                                                    fixStartDerivative=True, fixEndDerivative=True,
                                                    fixStartDirection=False, fixEndDirection=False)
            sd2Raw.append(sd2)

        # re-arrange the derivatives order
        d2RearrangedList = []
        for n2 in range(elementsCountUpNeck+1):
            for n1 in range(elementsCountAround):
                d2 = sd2Raw[n1][n2]
                d2RearrangedList.append(d2)

        # create tracksurface at the outer layer of the neck
        nodesOnTrackSurface = []
        nodesOnTrackSurface_d1 = []
        nodesOnTrackSurface_d2 = []
        for n2 in range(elementsCountUpNeck+1):
            for n1 in range(elementsCountAround):
                if (n1 <= elementsCountAround / 2):
                    nodesOnTrackSurface.append(interpolatedNodes[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface_d1.append(interpolatedNodes_d1[n2 * elementsCountAround + n1])
                    nodesOnTrackSurface_d2.append(d2RearrangedList[n2 * elementsCountAround + n1])

        # nodes and derivatives of the neck of the bladder
        listOuterNeck_x = []
        listOuterNeck_d1 = []
        listOuterNeck_d2 = []
        elementsCount1 = elementsCountAround // 2
        elementsCount2 = elementsCountUpNeck
        tracksurfaceOstium1 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface, nodesOnTrackSurface_d1,
                                    nodesOnTrackSurface_d2)
        ostium1Position = tracksurfaceOstium1.createPositionProportion(ostiumPositionAround, ostiumPositionUp)
        ostium1Position.xi1 = 1.0
        ostium1Position.xi2 = 1.0
        ostiumElementPositionAround = ostium1Position.e1
        ostiumElementPositionUp = ostium1Position.e2
        for n2 in range(len(interpolatedNodes)):
            listOuterNeck_x.append(interpolatedNodes[n2])
            listOuterNeck_d1.append(interpolatedNodes_d1[n2])
            listOuterNeck_d2.append(d2RearrangedList[n2])

        # create body of the bladder
        radiansPerElementAround = 2.0 * math.pi / elementsCountAround
        radiansPerElementUpBody = (3 * math.pi / 4) / elementsCountUpBody
        # create regular rows
        listOuterBody_x = []
        listOuterBody_d1 = []
        listOuterBody_d2 = []
        for n2 in range(1, elementsCountUpBody):
            radiansUp = (math.pi / 4) + n2 * radiansPerElementUpBody
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            majorRadius = 0.5 * majorDiameter * sinRadiansUp
            minorRadius = 0.5 * minorDiameter * sinRadiansUp
            for n1 in range(elementsCountAround):
                radiansAround = n1 * radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x = [
                    -majorRadius * sinRadiansAround,
                    minorRadius * cosRadiansAround,
                    -height * cosRadiansUp
                ]
                dx_ds1 = [
                    -majorRadius * cosRadiansAround * radiansPerElementAround,
                    minorRadius * -sinRadiansAround * radiansPerElementAround,
                    0.0
                ]
                dx_ds2 = [
                    -0.5 * majorDiameter * sinRadiansAround * cosRadiansUp * radiansPerElementUpBody,
                    0.5 * minorDiameter * cosRadiansAround * cosRadiansUp * radiansPerElementUpBody,
                    height*sinRadiansUp * radiansPerElementUpBody
                ]
                listOuterBody_x.append(x)
                listOuterBody_d1.append(dx_ds1)
                listOuterBody_d2.append(dx_ds2)

        # create outer apex node
        outerApexNode_x = []
        outerApexNode_d1 = []
        outerApexNode_d2 = []
        x = [0.0, 0.0, height]
        dx_ds1 = [height*radiansPerElementUpBody/2, 0.0, 0.0]
        dx_ds2 = [0.0, height*radiansPerElementUpBody/2, 0.0]
        outerApexNode_x.append(x)
        outerApexNode_d1.append(dx_ds1)
        outerApexNode_d2.append(dx_ds2)

        # set nodes of outer layer of the bladder
        listTotalOuter_x = listOuterNeck_x + listOuterBody_x + outerApexNode_x
        listTotalOuter_d1 = listOuterNeck_d1 + listOuterBody_d1 + outerApexNode_d1
        listTotalOuter_d2 = listOuterNeck_d2 + listOuterBody_d2 + outerApexNode_d2

        outerLayer_x = []
        outerLayer_d1 = []
        outerLayer_d2 = []
        for n2 in range(len(listTotalOuter_x)):
            if (n2 != (ostiumElementPositionUp + 1) * elementsCountAround + ostiumElementPositionAround + 1) and\
                    (n2 != (ostiumElementPositionUp + 1) * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, listTotalOuter_x[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, listTotalOuter_d1[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, listTotalOuter_d2[n2])
                nodeIdentifier += 1
                outerLayer_x.append(listTotalOuter_x[n2])
                outerLayer_d1.append(listTotalOuter_d1[n2])
                outerLayer_d2.append(listTotalOuter_d2[n2])

        # create and set nodes of inner layer of the bladder
        listTotalInner_x = []
        listTotalInner_d1 = []
        listTotalInner_d2 = []
        for n2 in range(elementsCountUpNeck + elementsCountUpBody):
            loop_x = [listTotalOuter_x[n2 * elementsCountAround + n1] for n1 in range(elementsCountAround)]
            loop_d1 = [listTotalOuter_d1[n2 * elementsCountAround + n1] for n1 in range(elementsCountAround)]
            loop_d2 = [listTotalOuter_d2[n2 * elementsCountAround + n1] for n1 in range(elementsCountAround)]
            for n1 in range(elementsCountAround):
                x, d1, _, _ = interp.projectHermiteCurvesThroughWall(loop_x, loop_d1, loop_d2, n1,
                                                                     -bladderWallThickness, loop=True)
                listTotalInner_x.append(x)
                listTotalInner_d1.append(d1)

        listInner_d2 = []
        for n2 in range(elementsCountAround):
            nx = [listTotalOuter_x[n1 * elementsCountAround + n2] for n1 in range(elementsCountUpNeck + elementsCountUpBody)]
            nd1 = [listTotalOuter_d1[n1 * elementsCountAround + n2] for n1 in range(elementsCountUpNeck + elementsCountUpBody)]
            nd2 = [listTotalOuter_d2[n1 * elementsCountAround + n2] for n1 in range(elementsCountUpNeck + elementsCountUpBody)]
            for n1 in range(elementsCountUpNeck + elementsCountUpBody):
                _, d2, _, _ = interp.projectHermiteCurvesThroughWall(nx, nd2, nd1, n1,
                                                                     bladderWallThickness, loop=False)
                listInner_d2.append(d2)

        # re-arrange the derivatives order
        for n2 in range(elementsCountUpNeck + elementsCountUpBody):
            for n1 in range(elementsCountAround):
                rearranged_d2 = listInner_d2[n1 * (elementsCountUpNeck + elementsCountUpBody) + n2]
                listTotalInner_d2.append(rearranged_d2)

        innerLayer_x = []
        innerLayer_d1 = []
        innerLayer_d2 = []
        for n2 in range(len(listTotalInner_x)):
            if (n2 != (ostiumElementPositionUp + 1) * elementsCountAround + ostiumElementPositionAround + 1) and \
                    (n2 != (ostiumElementPositionUp + 1) * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, listTotalInner_x[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, listTotalInner_d1[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, listTotalInner_d2[n2])
                nodeIdentifier += 1

                innerLayer_x.append(listTotalInner_x[n2])
                innerLayer_d1.append(listTotalInner_d1[n2])
                innerLayer_d2.append(listTotalInner_d2[n2])

        # create inner apex node
        x = [0.0, 0.0, height - bladderWallThickness]
        dx_ds1 = [height*radiansPerElementUpBody/2, 0.0, 0.0]
        dx_ds2 = [0.0, height*radiansPerElementUpBody/2, 0.0]
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
        listTotalInner_x.append(x)
        listTotalInner_d1.append(dx_ds1)
        listTotalInner_d2.append(dx_ds2)
        innerLayer_x.append(x)
        innerLayer_d1.append(dx_ds1)
        innerLayer_d2.append(dx_ds2)
        nodeIdentifier += 1

        # create ureters on the surface
        elementIdentifier = 1
        # ureter 1
        centerUreter1_x, centerUreter1_d1, centerUreter1_d2 = tracksurfaceOstium1.evaluateCoordinates(ostium1Position, derivatives=True)
        td1, td2, td3 = calculate_surface_axes(centerUreter1_d1, centerUreter1_d2, [1.0, 0.0, 0.0])
        m1 = ostiumElementPositionUp * elementsCountAround + ostiumElementPositionAround
        ureter1StartCornerx = listOuterNeck_x[m1]
        v1 = [(ureter1StartCornerx[c] - centerUreter1_x[c]) for c in range(3)]
        ostium1Direction = vector.crossproduct3(td3, v1)
        nodeIdentifier, elementIdentifier, (o1_x, o1_d1, o1_d2, _, o1_NodeId, o1_Positions) = \
            generateOstiumMesh(region, ostiumDefaultOptions, tracksurfaceOstium1, ostium1Position, ostium1Direction,
            startNodeIdentifier=nodeIdentifier, startElementIdentifier=elementIdentifier)

        # ureter 2
        tracksurfaceOstium2 = tracksurfaceOstium1.createMirrorX()
        ostium2Position = TrackSurfacePosition(elementsCountAround - ostiumElementPositionAround, ostiumElementPositionUp - 1, 0.0, 1.0)
        centerUreter2_x, centerUreter2_d1, centerUreter2_d2 = tracksurfaceOstium2.evaluateCoordinates(ostium2Position, derivatives =True)
        ad1, ad2, ad3 = calculate_surface_axes(centerUreter2_d1, centerUreter2_d2, [1.0, 0.0, 0.0])
        if elementsCountAroundOstium == 4:
            m2 = ostiumElementPositionUp * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 1
        else:
            m2 = ostiumElementPositionUp * elementsCountAround + elementsCountAround - ostiumElementPositionAround - 2
        ureter2StartCornerx = listOuterNeck_x[m2]
        v2 = [(ureter2StartCornerx[c] - centerUreter2_x[c]) for c in range(3)]
        ostium2Direction = vector.crossproduct3(ad3, v2)
        nodeIdentifier, elementIdentifier, (o2_x, o2_d1, o2_d2, _, o2_NodeId, o2_Positions) = \
            generateOstiumMesh(region, ostiumDefaultOptions, tracksurfaceOstium2, ostium2Position, ostium2Direction,
            startNodeIdentifier=nodeIdentifier, startElementIdentifier=elementIdentifier)

        # create annulus mesh around ostium
        endPoints1_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endNode1_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endDerivativesMap = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints2_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endNode2_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]

        nodeCountsEachWallLayer = (elementsCountUpNeck + elementsCountUpBody) * elementsCountAround - 1
        for n3 in range(2):
            n1 = 0
            endNode1_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround) + ostiumElementPositionAround + 1
            endNode1_Id[n3][n1 + 1] = endNode1_Id[n3][n1] + elementsCountAround
            endNode1_Id[n3][n1 + 2] = endNode1_Id[n3][n1 + 1] + elementsCountAround - 2
            endNode1_Id[n3][n1 + 3] = endNode1_Id[n3][n1 + 2] + 1
            endNode1_Id[n3][n1 + 4] = endNode1_Id[n3][n1 + 3] + 1
            endNode1_Id[n3][n1 + 5] = endNode1_Id[n3][n1 + 1] + 1
            endNode1_Id[n3][n1 + 6] = endNode1_Id[n3][n1] + 2
            endNode1_Id[n3][n1 + 7] = endNode1_Id[n3][n1] + 1
            if ostiumElementPositionAround == 0:
                endNode2_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround)\
                                        + elementsCountAround - ostiumElementPositionAround - 1
                endNode2_Id[n3][n1 + 1] = endNode2_Id[n3][n1] + elementsCountAround - 1
                endNode2_Id[n3][n1 + 2] = endNode2_Id[n3][n1 + 1] + elementsCountAround - 1
                endNode2_Id[n3][n1 + 3] = endNode2_Id[n3][n1 + 2] + 1
                endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] - elementsCountAround + 1
                endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 4] - elementsCountAround + 2
                endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1 + 5] - elementsCountAround
                endNode2_Id[n3][n1 + 7] = endNode2_Id[n3][n1] + 1
            else:
                endNode2_Id[n3][n1] = ((1 - n3) * nodeCountsEachWallLayer) + (ostiumElementPositionUp * elementsCountAround)\
                                        + elementsCountAround - ostiumElementPositionAround - 1
                endNode2_Id[n3][n1 + 1] = endNode2_Id[n3][n1] + elementsCountAround - 1
                endNode2_Id[n3][n1 + 2] = endNode2_Id[n3][n1 + 1] + elementsCountAround - 1
                endNode2_Id[n3][n1 + 3] = endNode2_Id[n3][n1 + 2] + 1
                endNode2_Id[n3][n1 + 4] = endNode2_Id[n3][n1 + 3] + 1
                endNode2_Id[n3][n1 + 5] = endNode2_Id[n3][n1 + 1] + 1
                endNode2_Id[n3][n1 + 6] = endNode2_Id[n3][n1] + 2
                endNode2_Id[n3][n1 + 7] = endNode2_Id[n3][n1] + 1

        for n3 in range(2):
            for n1 in range(elementsCountAroundOstium):
                nc1 = endNode1_Id[n3][n1] - (1 - n3) * nodeCountsEachWallLayer - 1
                endPoints1_x[n3][n1] = innerLayer_x[nc1]
                endPoints1_d1[n3][n1] = innerLayer_d1[nc1]
                endPoints1_d2[n3][n1] = [innerLayer_d2[nc1][c] for c in range(3)]
                nc2 = endNode2_Id[n3][n1] - (1 - n3) * nodeCountsEachWallLayer - 1
                endPoints2_x[n3][n1] = innerLayer_x[nc2]
                endPoints2_d1[n3][n1] = innerLayer_d1[nc2]
                endPoints2_d2[n3][n1] = innerLayer_d2[nc2]

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
            elementsCountRadial=elementsCountAnnulusRadially, meshGroups=[neckMeshGroup, urinaryBladderMeshGroup])

        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            o2_x, o2_d1, o2_d2, None, o2_NodeId, None,
            endPoints2_x, endPoints2_d1, endPoints2_d2, None, endNode2_Id, endDerivativesMap,
            elementsCountRadial=elementsCountAnnulusRadially, meshGroups=[neckMeshGroup, urinaryBladderMeshGroup])

        # create elements
        for e3 in range(1):
            newl = (e3 + 1) * ((elementsCountUpNeck + elementsCountUpBody) * elementsCountAround - 1)
            # create bladder neck elements
            for e2 in range(elementsCountUpNeck):
                for e1 in range(elementsCountAround):
                    if e2 == ostiumElementPositionUp:
                        if (e1 == ostiumElementPositionAround or e1 == ostiumElementPositionAround + 1):
                            pass
                        elif (e1 == elementsCountAround - ostiumElementPositionAround - 2 or e1 == elementsCountAround - 1 - ostiumElementPositionAround):
                            pass
                        else:
                            bni1 = e2 * elementsCountAround + e1 + 1
                            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                            if e1 < ostiumElementPositionAround:
                                bni3 = bni1 + elementsCountAround
                                bni4 = bni2 + elementsCountAround
                            elif (ostiumElementPositionAround + 1 < e1 < elementsCountAround - ostiumElementPositionAround - 2):
                                bni3 = bni1 + elementsCountAround - 1
                                bni4 = bni2 + elementsCountAround - 1
                            elif e1 > elementsCountAround - ostiumElementPositionAround - 1:
                                bni3 = bni1 + elementsCountAround - 2
                                if e1 == elementsCountAround - 1:
                                    bni4 = bni2 + elementsCountAround
                                else:
                                    bni4 = bni2 + elementsCountAround - 2
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni4 + newl,
                                               bni1, bni2, bni3, bni4]
                            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                            neckMeshGroup.addElement(element)
                            urinaryBladderMeshGroup.addElement(element)
                            elementIdentifier += 1
                    elif e2 == ostiumElementPositionUp + 1:
                        if (e1 == ostiumElementPositionAround or e1 == ostiumElementPositionAround + 1):
                            pass
                        elif (e1 == elementsCountAround - ostiumElementPositionAround - 2 or e1 == elementsCountAround - 1 - ostiumElementPositionAround):
                            pass
                        else:
                            if e1 < ostiumElementPositionAround:
                                bni1 = e2 * elementsCountAround + e1 + 1
                                bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                                bni3 = bni1 + elementsCountAround - 2
                                bni4 = bni2 + elementsCountAround - 2
                            elif (ostiumElementPositionAround + 1 < e1 < elementsCountAround - ostiumElementPositionAround - 2):
                                bni1 = e2 * elementsCountAround + e1
                                bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround
                                bni3 = bni1 + elementsCountAround - 1
                                bni4 = bni2 + elementsCountAround - 1
                            elif e1 > elementsCountAround - ostiumElementPositionAround - 1:
                                bni1 = e2 * elementsCountAround + e1 - 1
                                bni3 = bni1 + elementsCountAround
                                if e1 == elementsCountAround - 1:
                                    bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                                    bni4 = bni2 + elementsCountAround - 2
                                else:
                                    bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround - 1
                                    bni4 = bni2 + elementsCountAround
                            element = mesh.createElement(elementIdentifier, elementtemplate)
                            nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni4 + newl,
                                               bni1, bni2, bni3, bni4]
                            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                            neckMeshGroup.addElement(element)
                            urinaryBladderMeshGroup.addElement(element)
                            elementIdentifier += 1
                    elif e2 > ostiumElementPositionUp + 1:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni1 = e2 * elementsCountAround + e1 - 1
                        bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround - 1
                        bni3 = bni1 + elementsCountAround
                        bni4 = bni2 + elementsCountAround
                        nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni4 + newl,
                                           bni1, bni2, bni3, bni4]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        neckMeshGroup.addElement(element)
                        urinaryBladderMeshGroup.addElement(element)
                        elementIdentifier += 1
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni1 = e2 * elementsCountAround + e1 + 1
                        bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                        bni3 = bni1 + elementsCountAround
                        bni4 = bni2 + elementsCountAround
                        nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni4 + newl,
                                           bni1, bni2, bni3, bni4]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        neckMeshGroup.addElement(element)
                        urinaryBladderMeshGroup.addElement(element)
                        elementIdentifier += 1

            # create bladder body elements
            for e2 in range(elementsCountUpNeck, (elementsCountUpNeck + elementsCountUpBody - 1)):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni1 = e2 * elementsCountAround + e1 - 1
                    bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround - 1
                    bni3 = bni1 + elementsCountAround
                    bni4 = bni2 + elementsCountAround
                    nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni4 + newl,
                                       bni1, bni2, bni3, bni4]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    bodyMeshGroup.addElement(element)
                    urinaryBladderMeshGroup.addElement(element)
                    elementIdentifier += 1
            # create apex elements
            bni3 = (elementsCountUpNeck + elementsCountUpBody) * elementsCountAround - 1
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
                nodeIdentifiers = [bni1 + newl, bni2 + newl, bni3 + newl, bni1, bni2, bni3]
                result = element.setNodesByIdentifier(eftApex, nodeIdentifiers)
                bodyMeshGroup.addElement(element)
                urinaryBladderMeshGroup.addElement(element)
                elementIdentifier += 1

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            numberInXi1 = refineElementsCountAround
            numberInXi2 = refineElementsCountUp
            numberInXi3 = refineElementsCountThroughWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            element = meshrefinement._sourceElementiterator.next()

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
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("neck of urinary bladder"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("dome of the bladder"))
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

        SerosaOfNeck =findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("serosa of neck of urinary bladder"))
        SerosaOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_serosa)
        lumenOfNeck = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("lumen of neck of urinary bladder"))
        lumenOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_lumen)

