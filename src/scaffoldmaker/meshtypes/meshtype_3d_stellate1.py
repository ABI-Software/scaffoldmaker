"""
Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, \
    findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.stellate_terms import get_stellate_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftLocalNodes
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.matrix import rotateAboutZAxis
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.vector import setMagnitude


class MeshType_3d_stellate1(Scaffold_base):
    """
    Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
    """

    mouseMeanMesh = {
        'meshEdits':
            [[[-242, -142.9, -96.77], [-249, -54.64, -3.137], [77.02, -112.8, 2.094]],
             [[56.72, -109.7, -95.44], [426.9, -47.3, 0.4029], [13.12, 118.7, -1.919]],
             [[463.5, -143.6, -98.9], [410.9, -33.87, -5.144], [13.43, 112, -3.136]],
             [[923.8, -181.1, -110.7], [526, -53.8, -10.86], [11.61, 108.6, -5.169]],
             [[-379.9, 60.5, -101.7], [194.5, -284.8, 7.069], [264.6, 224.9, 2.306]],
             [[70.7, 8.808, -97.41], [425.4, -45.37, 0.214], [15.02, 119.3, -2.076]],
             [[476.7, -31.45, -102.2], [409.7, -37.01, -6.137], [13.2, 112.1, -3.111]],
             [[935.3, -72.22, -115], [526.2, -53.76, -10.76], [9.611, 109.5, -3.776]],
             [[1327, -116.7, -135], [224, -25.64, -10.8], [18.59, 199.8, -8.407]],
             [[-203.5, 220.3, -100.9], [159, -73.66, 3.295], [96.56, 88.12, -0.277]],
             [[87.07, 127.1, -99.29], [421.9, -42.15, -0.837], [16.73, 118.9, -2.162]],
             [[489.5, 80.41, -105], [405.7, -40.24, -6.158], [12.96, 111.6, -2.566]],
             [[943.6, 38.12, -118.9], [516.6, -52.51, -8.843], [7.786, 110.8, -4.435]],
             [[-242.4, -141.4, -23.45], [-248.8, -54.59, -3.241], [77.04, -112.7, 2.425]],
             [[56.34, -108.3, -22.1], [426.7, -47.59, 1.009], [13.15, 118.9, -2.095]],
             [[463.2, -142.3, -25.31], [410.8, -33.73, -5.277], [13.43, 112.2, -3.372]],
             [[923.6, -179.7, -35.87], [526, -53.9, -11.26], [11.57, 108.6, -4.45]],
             [[-380.3, 61.98, -28.76], [194.4, -284.9, 7.396], [264.7, 224.9, 2.089]],
             [[70.32, 10.27, -24.13], [425.3, -45.26, 0.8815], [15, 119.4, -1.745]],
             [[476.4, -30.02, -28.41], [409.8, -37.1, -4.907], [13.19, 112.2, -3.329]],
             [[935.1, -70.86, -40.15], [526.3, -53.94, -9.548], [9.58, 109.5, -3.552]],
             [[1327, -115.3, -58.62], [224, -26.02, -10.31], [18.66, 199.8, -6.379]],
             [[-203.8, 221.8, -28.04], [159.1, -73.78, 3.567], [96.63, 88.24, -0.6312]],
             [[86.66, 128.5, -25.68], [421.7, -42.31, -0.5815], [16.68, 118.9, -1.337]],
             [[489.2, 81.9, -31.44], [405.6, -40.42, -4.792], [12.96, 111.6, -3.211]],
             [[943.3, 39.51, -43.96], [516.6, -52.48, -8.07], [7.717, 110.8, -3.905]],
             [[-408.8, 350.8, -117], [-123.6, 189.7, -39.89], [-129.5, -43.37, -3.048]],
             [[-540.3, 307.4, -120.2], [-108.3, 188.8, -41.46], [-131.6, -41.38, -2.229]],
             [[-596, 406.9, -161.3], [-25.94, 44.84, -26.59], [-239.5, -77.71, -3.878]],
             [[-694.8, 103.5, -108.8], [30.37, 104.1, 0.9179], [-172.7, 23.53, -4.444]],
             [[-676.1, 269.5, -121.6], [-88.31, 186.1, -40.49], [-132.9, -39.08, -1.443]],
             [[-409.2, 352.4, -43.89], [-123.7, 189.9, -39.85], [-129.5, -43.28, -3.529]],
             [[-540.7, 308.9, -47.28], [-108.3, 189, -41.12], [-131.6, -41.41, -2.472]],
             [[-596.5, 408.5, -88.38], [-25.87, 44.65, -26.49], [-239.3, -77.46, -4.69]],
             [[-695.1, 105, -35.83], [30.06, 104.1, 1.195], [-172.7, 23.52, -4.408]],
             [[-676.4, 271.1, -48.62], [-88.47, 186.3, -40.28], [-132.8, -39.08, -1.369]],
             [[-795.9, -39.7, -121.2], [-269.6, -114.8, -32.52], [147.6, -86.96, 6.78]],
             [[-649, -124.3, -115.2], [-263.3, -124.6, -30.02], [146.9, -88.49, 6.997]],
             [[-822.5, -187.8, -141.6], [-91.52, -28.67, -15.87], [269.5, -165.5, 14.64]],
             [[-503.2, -214.7, -108.1], [-263.1, -129.7, -28.88], [145.5, -92.58, 7.724]],
             [[-796.2, -38.2, -48.06], [-269.7, -115, -31.62], [147.5, -86.96, 6.869]],
             [[-649.4, -122.8, -42.02], [-263.2, -124.7, -29.45], [146.9, -88.44, 6.967]],
             [[-822.9, -186.4, -67.74], [-91.59, -28.64, -15.31], [269.5, -165.4, 15.28]],
             [[-503.5, -213.2, -34.94], [-262.9, -129.5, -28.54], [145.5, -92.62, 7.544]]]
                    }
    @staticmethod
    def getName():
        return '3D Stellate 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse cervicothoracic ganglion 1',
            'Mouse cervicothoracic ganglion mean 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        options['Base parameter set'] = parameterSetName

        isMouse = 'Mouse' in parameterSetName

        if isMouse:
            options['Numbers of elements along arms'] = [4,2,2]
        else:
            options['Numbers of elements along arms'] = [4,2,2]
        options['Element width central'] = 0.8
        options['Element length along arm'] = 0.8
        options['Element width across arm'] = 0.3
        options['Element thickness'] = 0.1
        options['Refine'] = False
        options['Refine number of elements along arm'] = 1
        options['Refine number of elements across arm'] = 1
        options['Refine number of elements through thickness'] = 1
        return options


    @staticmethod
    def getOrderedOptionNames():

        return [
            'Numbers of elements along arms',
            'Element width central',
            'Element length along arm',
            'Element width across arm',
            'Element thickness',
            'Refine',
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness'
        ]

    @classmethod
    def checkOptions(cls, options):
        for key in [
            'Refine number of elements along arm',
            'Refine number of elements across arm',
            'Refine number of elements through thickness']:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Element width central',
            'Element length along arm',
            'Element width across arm',
            'Element thickness']:
            if options[key] < 0:
                options[key] = 0.1

        armCountsKey = 'Numbers of elements along arms'
        armCount = 3
        countsLen = len(options[armCountsKey])
        if countsLen < armCount:
            options[armCountsKey] += [2 for i in range(armCount - countsLen)]
        elif countsLen > armCount:
            options[armCountsKey] = options[armCountsKey][:armCount]
        for i, numberOfElements in enumerate(options[armCountsKey]):
            if numberOfElements < 2:
                options[armCountsKey][i] = 2

        parameterSetName = options['Base parameter set']
        isMouse = 'Mouse' in parameterSetName
        isMean = 'mean' in parameterSetName
        if isMouse and isMean:
            options['Numbers of elements along arms'] = [4,2,2]

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isDefault = 'Default' in parameterSetName
        isMouse = 'Mouse' in parameterSetName
        isMean = 'mean' in parameterSetName

        fm = region.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = findOrCreateFieldCoordinates(fm)
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        armCount = 3
        elementLengthCentral = options['Element width central']
        elementLengths = [options['Element length along arm'],
                          options['Element width across arm'],
                          options['Element thickness']]
        elementsCountsAlongArms = options['Numbers of elements along arms']
        elementsCount2 = 2
        elementsCount3 = 1
        useCrossDerivatives = False
        # arm group annotations for user
        armTerms, _ = getAutomaticArmFaceTerms(armCount)
        armGroups = [AnnotationGroup(region, armTerm) for armTerm in armTerms]
        stellateTerm = get_stellate_term("cervicothoracic ganglion") if isMouse else ("stellate", None)
        stellateGroup = AnnotationGroup(region, stellateTerm)
        annotationGroups = [stellateGroup] + armGroups
        armMeshGroups = [a.getMeshGroup(mesh) for a in armGroups]
        stellateMeshGroup = stellateGroup.getMeshGroup(mesh)

        # markers with element number and xi position
        allMarkers = {}
        if isMouse:
            xProportion = {}
            xProportion['ICN'] = 0.9
            xProportion['VA'] = 0.9
            xProportion['DA'] = 0.9
            xProportion['C8'] = 0.9
            xProportion['T1'] = 0.25
            xProportion['T2'] = 0.5
            xProportion['T3'] = 0.75
            xProportion['TST'] = 1
            armNumber = {}
            armNumber['ICN'] = 2
            armNumber['VA'] = 2
            armNumber['DA'] = 3
            armNumber['C8'] = 3
            armNumber['T1'] = 1
            armNumber['T2'] = 1
            armNumber['T3'] = 1
            armNumber['TST'] = 1
            nerveAbbrev = list(xProportion.keys())
            elementIndex = {}
            xi1 = {}
            for nerve in nerveAbbrev:
                elementIndex[nerve] = int(xProportion[nerve] * elementsCountsAlongArms[armNumber[nerve] - 1])
                xi1[nerve] = 1 if xProportion[nerve] == 1 else xProportion[nerve] * elementsCountsAlongArms[
                    armNumber[nerve] - 1] - elementIndex[nerve]
                elementIndex[nerve] += 1 if xProportion[nerve] < 1 else 0

            allMarkers = {"Inferior cardiac nerve": {"elementID": elementIndex['ICN'] + 2 * elementsCountsAlongArms[0],
                                                     "xi": [xi1['ICN'], 0.0, 0.5]},
                          "Ventral ansa subclavia": {
                              "elementID": elementIndex['VA'] + 2 * elementsCountsAlongArms[0] + elementsCountsAlongArms[1],
                              "xi": [xi1['VA'], 1.0, 0.5]},
                          "Dorsal ansa subclavia": {
                              "elementID": elementIndex['DA'] + 2 * (elementsCountsAlongArms[0] + elementsCountsAlongArms[1]),
                              "xi": [xi1['DA'], 0.0, 0.5]},
                          "Cervical spinal nerve 8": {
                              "elementID": elementIndex['C8'] + 2 * (elementsCountsAlongArms[0] + elementsCountsAlongArms[1]) +
                                           elementsCountsAlongArms[2], "xi": [xi1['C8'], 1.0, 0.5]},
                          "Thoracic spinal nerve 1": {"elementID": elementIndex['T1'], "xi": [xi1['T1'], 0.0, 0.5]},
                          "Thoracic spinal nerve 2": {"elementID": elementIndex['T2'], "xi": [xi1['T2'], 0.0, 0.5]},
                          "Thoracic spinal nerve 3": {"elementID": elementIndex['T3'], "xi": [xi1['T3'], 0.0, 0.5]},
                          "Thoracic sympathetic nerve trunk": {"elementID": elementIndex['TST'],
                                                               "xi": [xi1['TST'], 1.0, 0.5]},
                          }
            markerGroup = findOrCreateFieldGroup(fm, "marker")
            markerName = findOrCreateFieldStoredString(fm, name="marker_name")
            markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

            markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
            markerTemplateInternal = nodes.createNodetemplate()
            markerTemplateInternal.defineField(markerName)
            markerTemplateInternal.defineField(markerLocation)

        # Create nodes
        nodeIdentifier = 1
        minArmAngle = 2 * math.pi / armCount
        halfArmArcAngleRadians = minArmAngle / 2
        if not isMean:
            dipMultiplier = 1
            for na in range(armCount):
                elementsCount_i = [elementsCountsAlongArms[na], elementsCount2, elementsCount3]
                x, ds1, ds2, nWheelEdge = createArm(halfArmArcAngleRadians, elementLengths, elementLengthCentral, elementsCount_i, dipMultiplier, armCount, na)
                for ix in range(len(x)):
                    if na == 0 or ix not in nWheelEdge:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[ix])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1[ix])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2[ix])
                        nodeIdentifier += 1
        else:
            x_dx_all = cls.mouseMeanMesh['meshEdits']
            xyz_all = [x[0] for x in x_dx_all]
            dxyz = [[x[1], x[2]] for x in x_dx_all]
            nodeIdentifier = 1
            for i, nx in enumerate(xyz_all):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dxyz[i][0])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dxyz[i][1])
                nodeIdentifier += 1
        nodesCountsPerArm = [0] + [((elementsCount2+1)*e+1)*2 for e in elementsCountsAlongArms]

        # Create elements
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives() #createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementIdentifier = 1

        cumNodesCountsPerArm = [sum(nodesCountsPerArm[:i + 1]) for i in range(len(nodesCountsPerArm))]
        nCentre = [elementsCountsAlongArms[0]+1, int(nodesCountsPerArm[1]/2) + elementsCountsAlongArms[0]+1]
        for na in range(armCount):
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCountsAlongArms[na]):
                        scalefactors = None
                        ### NODES ###
                        no2 = (elementsCountsAlongArms[na] + 1)
                        no3 = (elementsCount2 + 1) * no2 - 2
                        offset = (cumNodesCountsPerArm[na])
                        bni = e3 * no3 + e2 * no2 + e1 + 1 + offset
                        if e2 == 0:
                            if e1 == 0 and na > 0: # and na < armCount -1: # wheelSouth
                                nWh = cumNodesCountsPerArm[na - 1] + (2 * elementsCountsAlongArms[na - 1]) + 2
                                nplUq = int(nodesCountsPerArm[na+1]/2) - elementsCountsAlongArms[na] # unused nodes at centre and shared edge
                                npl = int(nodesCountsPerArm[na+1]/2)  #  nodes at centre and shared edge
                                if na < armCount-1:
                                    cn = cumNodesCountsPerArm[na] + elementsCountsAlongArms[na]-2
                                    no2 = cumNodesCountsPerArm[na]
                                    em = elementsCountsAlongArms[na]
                                    nwPrev = [nWh, nWh + int(nodesCountsPerArm[na]/2)] # previous arm's edge, depends on armCount.
                                    nodeIdentifiers = [nwPrev[0], no2 + 1, nCentre[0], no2 + em,
                                                       nwPrev[1], no2 + em - 1 + nplUq,
                                                       nCentre[1], bni + (4 * em) - 2]
                                else:
                                    nplPrev = int(nodesCountsPerArm[na]/2) - 2
                                    no2 = elementsCountsAlongArms[na]-1
                                    no3 = int(nodesCountsPerArm[na+1]/2) - 3
                                    nwPrev = [cumNodesCountsPerArm[na-1] + 2*(elementsCountsAlongArms[na-1]),
                                              cumNodesCountsPerArm[na-1] + 2*(elementsCountsAlongArms[na-1]) + nplPrev]
                                    start = cumNodesCountsPerArm[na] - 3
                                    nodeIdentifiers = [nwPrev[0], start,
                                                       nCentre[0], start + no2,
                                                       nwPrev[1],  start + no3,
                                                       nCentre[1], start + no2 +no3]
                            elif e1 == elementsCountsAlongArms[na] - 1: # armEnd, south
                                if na == 0:
                                    nodeIdentifiers = [bni, bni + no2 - 1, bni + no2, bni + no3,
                                                       bni + no2 + no3 - 1, bni + no2 + no3]
                                else:
                                    no3 = armCount*elementsCountsAlongArms[na] - 1
                                    no2 = elementsCountsAlongArms[na]
                                    if na > 1:
                                        bni -= 4
                                        no3 -= 1
                                    nodeIdentifiers = [bni-1, bni + no2-2, bni + no2 - 1,
                                                       bni + no3 - 1, bni + no2 -2 + no3, bni + no2 + no3 - 1]
                            elif na > 0 and e1 > 0: #  [na=1+, e1=1+, e2=0] for len=3+
                                bni -= 1 + ((armCount+1)*(na-1))
                                no2 = elementsCountsAlongArms[na]
                                no3 = armCount*no2 - (na-1) - 1
                                nodeIdentifiers = [bni, bni + 1, bni + no2 - 1, bni + no2,
                                                   bni + no3, bni + no3 + 1,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                            else:
                                nodeIdentifiers = [bni, bni + 1, bni + no2 - 1, bni + no2, bni + no3, bni + no3 + 1,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]
                        else:
                            if e1 == 0 and na > 0:  # and na < armCount -1: # wheelNorth
                                if na < armCount - 1:
                                    bni -= armCount
                                    npl = int(nodesCountsPerArm[na+1] / 2) - 2
                                    no2 = elementsCountsAlongArms[na]
                                    nodeIdentifiers = [nCentre[0], bni + 1, bni + no2 + 1, bni + no2 + 2,
                                                       nCentre[1], bni + npl + 1, bni + npl + no2 + 1, bni + npl + no2 + 2]
                                else: # last arm
                                    bni = cumNodesCountsPerArm[na] - 2 - (armCount - elementsCountsAlongArms[na])
                                    nodeIdentifiers = [nCentre[0], bni + 1, 1, bni + no2,
                                                       nCentre[1], bni + no3 - 2,
                                                       int(nodesCountsPerArm[1]/2)+1, bni + no2 + no3 - armCount]
                            elif e1 == elementsCountsAlongArms[na] - 1: # armEnd north
                                if na > 0:
                                    no2 = elementsCountsAlongArms[na]
                                    nplUq = int(nodesCountsPerArm[na + 1] / 2) - 2
                                    if na > 1:
                                        adj = na - 1
                                        bni -= armCount*na + (armCount-elementsCountsAlongArms[na]) + 1
                                        if elementsCountsAlongArms[na] < 3:
                                           bni += 1
                                        if elementsCountsAlongArms[na]>3:
                                            bni -= elementsCountsAlongArms[na] - 3
                                        no2 += 1 - adj
                                        no3 = nplUq - adj
                                        nodeIdentifiers = [bni, bni +1, bni+no2,
                                                           bni+no3, bni +no3+1, bni+no2 +no3]
                                    else:
                                        bni -= armCount
                                        nodeIdentifiers = [bni, bni +1, bni+no2 + 1,
                                                           bni+nplUq, bni +nplUq+1, bni+no2 +nplUq+ 1]
                                else:
                                    nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no3 - 1, bni + no3,
                                                       bni + no2 + no3 - 1]

                            elif na > 0 and e1 > 0: #  [na=1+, e1=1+, e2=1] for len=3+
                                adj = na - 1
                                bni -= armCount *na + adj
                                no2 -= adj
                                k = armCount*elementsCountsAlongArms[na] - na
                                nodeIdentifiers = [bni, bni+1, bni + no2, bni + no2 + 1,
                                                   bni + k, bni + k + 1,
                                                   bni + no2 + k, bni + no2 + k + 1]
                            else:
                                nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no2,
                                                   bni + no3 - 1, bni + no3,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]

                        if e1 == 0:  # wheel
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            if armCount == 3:
                                if e2 == 0:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    scaleEftNodeValueLabels(eft1, [1, 5], [Node.VALUE_LABEL_D_DS1,
                                                             Node.VALUE_LABEL_D_DS2], [1])
                                    ns = [3, 7]
                                else:
                                    ns = [1, 5]
                                if na == 0:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                    elif e2 == 1:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS2, [1])])
                                elif na == 2:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                    elif e2 == 1:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                            elif armCount == 4:
                                if e2 == 0:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    scaleEftNodeValueLabels(eft1, [1, 5], [Node.VALUE_LABEL_D_DS1,
                                                             Node.VALUE_LABEL_D_DS2], [1])
                                    ns = [3, 7]
                                else:
                                    ns = [1, 5]
                                if na == 0:
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 2:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                elif na == 3:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])

                        elif e1 < (elementsCountsAlongArms[na] - 1):
                            eft1 = eft
                            elementtemplate1 = elementtemplate
                        else:
                            # rounded ends of arms. Collapse xi2 at xi1 = 1
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2, [])
                            if e2 == 0:
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                nodeIdentifiers = [nodeIdentifiers[0], nodeIdentifiers[2], nodeIdentifiers[1],
                                                   nodeIdentifiers[3], nodeIdentifiers[5], nodeIdentifiers[4]]
                            else:  # e2 == 1
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                            ln_map = [1, 2, 3, 2, 4, 5, 6, 5]
                            remapEftLocalNodes(eft1, 6, ln_map)

                        if eft1 is not eft:
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        result3 = element.setScaleFactors(eft1, scalefactors) if scalefactors else None

                        # add to meshGroup
                        stellateMeshGroup.addElement(element)
                        armMeshGroups[na].addElement(element)

                        elementIdentifier += 1

        # annotation fiducial points
        if isMouse:
            for key in allMarkers:
                xi = allMarkers[key]["xi"]
                addMarker = {"name": key, "xi": allMarkers[key]["xi"]}

                markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
                nodeIdentifier += 1
                cache.setNode(markerPoint)
                markerName.assignString(cache, addMarker["name"])
                elementID = allMarkers[key]["elementID"]
                element = mesh.findElementByIdentifier(elementID)
                markerLocation.assignMeshLocation(cache, element, addMarker["xi"])

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCount1 = options['Refine number of elements along arm']
        refineElementsCount2 = options['Refine number of elements across arm']
        refineElementsCount3 = options['Refine number of elements through thickness']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCount1, refineElementsCount2, refineElementsCount3)


    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add point annotation groups from the 1D mesh.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New point annotation groups are appended to this list.
        """
        # create  groups
        fm = region.getFieldmodule()
        armCount = len(options['Numbers of elements along arms'])
        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
        is_exterior_face_xi2_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))
        armTerms, faceTerms = getAutomaticArmFaceTerms(armCount)

        armGroups = [getAnnotationGroupForTerm(annotationGroups, armTerm) for armTerm in armTerms]
        isArm =[armGroup.getFieldElementGroup(mesh2d) for armGroup in armGroups]
        for arm in range(armCount):
            is_face = fm.createFieldOr(
            fm.createFieldAnd(isArm[arm - 1], is_exterior_face_xi2_1),
            fm.createFieldAnd(isArm[arm], is_exterior_face_xi2_0))
            faceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, faceTerms[arm - 1])
            faceGroup.getMeshGroup(mesh2d).addElementsConditional(is_face)

def getAutomaticArmFaceTerms(armCount):
    armTerms = [("stellate arm %d"%(i), None) for i in range(1,armCount+1)]
    faceTerms = [("stellate face %d-%d" % (i, (i % armCount) + 1), None) for i in range(1, armCount + 1)]
    return armTerms, faceTerms

def createArm(halfArmArcAngleRadians, elementLengths, elementLengthCentral, elementsCount, dipMultiplier, armCount, armIndex):
    """
    Create single arm unit.
    Base length of element is 1.
    Direction: anticlockwise.
    Minimum arm length is 2 elements.
    :param halfArmArcAngleRadians: angle arising from base node (rad)
    :param elementLengths: [Element length along arm, half Element width across arm, Element thickness]
    :param elementLengthCentral: Radial element length about central node.
    :param elementsCount: list of numbers of elements along arm length, across width and through thickness directions
    :param dipMultiplier: factor that wheel nodes protrude by, relative to unit length
    :param armCount: number of arms in body
    :param armIndex: 0-based
    :return: x: positions of nodes in arm
    :return: nodeIdentifier: number of last node in arm +1
    :return: xnodes_ds1: ds1 derivatives of nodes associated with x
    :return: xnodes_ds2: ds2 derivatives of nodes associated with y
    :return: rmVertexNodes: indices of nodes at arm vertices around central node (including central node)
    """

    elementsCount1, elementsCount2, elementsCount3 = elementsCount
    [elementLength, elementWidth, elementHeight] = elementLengths

    shorterArmEnd = True
    armAngle = 2*halfArmArcAngleRadians*armIndex
    armAngleConst = 2*halfArmArcAngleRadians

    xnodes_ds1 = []
    xnodes_ds2 = []
    dx_ds1 = [elementLength, 0.0, 0.0]
    dx_ds2 = [0.0, elementWidth, 0.0]

    nodes_per_layer = (elementsCount1 + 1) * (elementsCount2 + 1) - 2 # discount 2 armEnd corners
    x = []

    # nCentre and rmVertexNodes are in pythonic indexing
    nCentre = elementsCount1
    nCentre = [nCentre, nCentre + nodes_per_layer]
    if armIndex == 0:
        rmVertexNodes = []
    elif armIndex != armCount - 1:
        rmVertexNodes = nCentre + [0, nodes_per_layer]
    else:
        rmVertexNodes = nCentre + [0, nodes_per_layer,
                          2* (elementsCount1+1) - 1, 2* (elementsCount1+1) - 1 + nodes_per_layer ]
    dcent = [elementLengthCentral * math.cos(armAngleConst / 2), elementLengthCentral * math.sin(armAngleConst / 2), 0.0]
    dvertex = [elementLengthCentral * dipMultiplier * math.cos(halfArmArcAngleRadians),
                elementLengthCentral * dipMultiplier * math.sin(halfArmArcAngleRadians)]

    dipLength = 0.5*(elementLengthCentral + elementWidth)*dipMultiplier
    dvertex = [dipLength * math.cos(halfArmArcAngleRadians),
               dipLength * math.sin(halfArmArcAngleRadians)]
    dipMag = 2*dipLength - elementLengthCentral

    nid = 0
    for e3 in range(elementsCount3 + 1):
        x3 = e3 * elementHeight
        for e2 in range(elementsCount2 + 1):
            for e1 in range(elementsCount1 + 1):
                # ignore if armEnd corner nodes
                if (e1 == elementsCount1) and ((e2 == 0) or (e2 == elementsCount2)):
                    pass
                else:
                    if e1 == 0 and e2 == 0:
                        x1 = dvertex[0]
                        x2 = -dvertex[1]
                    elif e1 == 0 and e2 == elementsCount2:
                        x1 = dvertex[0]
                        x2 = dvertex[1]
                    else:
                        if e1 == elementsCount1 and shorterArmEnd:
                            x1 = 0.5*(elementLength+elementWidth) + elementLength*(e1-1)
                        else:
                            x1 = elementLength*e1
                        x2 = (e2 - 1) * elementWidth
                    x.append([x1, x2, x3])

                    # DERIVATIVES
                    ds1 = dx_ds1.copy()
                    ds2 = dx_ds2.copy()
                    if e2 == 1 and e1 == 0:
                        ds2 = dcent
                        ds1 = [dcent[0], -dcent[1], dcent[2]]
                    elif (e1 == 0) and ((e2 == 0) or (e2 == elementsCount2)):
                        ds2 = [dcent[0], -dcent[1], dcent[2]] if (e2 == 0) else dcent
                        ds2 = setMagnitude(ds2, dipMag)
                        ds1 = rotateAboutZAxis(ds2, -math.pi / 2)
                    elif e1 == elementsCount1 and e2 == elementsCount2-1: # armEnd
                        ds1 = [elementWidth,0,0]
                        ds2 = [0, 1 * (elementLength + elementWidth), 0]
                        if shorterArmEnd:
                            ds2 = [0, 0.5 * (elementLength + elementWidth), 0]
                    xnodes_ds1.append(ds1)
                    xnodes_ds2.append(ds2)
                    nid += 1

    # d1 derivatives at dips around central node
    nidDip = [0, elementsCount1 * 2 + 1,
              nodes_per_layer, elementsCount1 * 2 + 1 + nodes_per_layer]
    for nid in nidDip:
        nx = [x[nid], x[nid+1]]
        nd1 = [xnodes_ds1[nid], xnodes_ds1[nid+1]]
        ds1Smooth = smoothCubicHermiteDerivativesLine(nx, nd1,
                                           fixStartDirection=True, fixEndDerivative=True)
        xnodes_ds1[nid] = ds1Smooth[0]


    # rotate entire arm about origin by armAngle except for centre nodes
    tol = 1e-12
    for j, n in enumerate(x):
        xynew = rotateAboutZAxis(n, armAngle)
        [xnew, ynew] = [r for r in xynew][:2]
        if (abs(xnew) < tol):
            xnew = 0
        if (abs(ynew) < tol):
            ynew = 0
        x[j] = [xnew, ynew, x[j][2]]

        xnodes_ds1[j] = rotateAboutZAxis(xnodes_ds1[j], armAngle)
        xnodes_ds2[j] = rotateAboutZAxis(xnodes_ds2[j], armAngle)

    return (x, xnodes_ds1, xnodes_ds2, rmVertexNodes)
