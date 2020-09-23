"""
Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.stellate_terms import get_stellate_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftLocalNodes
from scaffoldmaker.utils.matrix import rotateAboutZAxis
from scaffoldmaker.utils.vector import magnitude, setMagnitude
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine

class MeshType_3d_stellate1(Scaffold_base):
    """
    Generates a 3-D planar stellate mesh with cross arms radiating from a central node, and variable numbers of elements along each arm.
    """
    @staticmethod
    def getName():
        return '3D Stellate 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1']

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

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        isDefault = 'Default' in options['Base parameter set']
        isMouse = 'Mouse' in options['Base parameter set']

        armCount = 3
        elementLengthCentral = options['Element width central']
        elementLengths = [options['Element length along arm'],
                 options['Element width across arm'],
                 options['Element thickness']]
        elementsCount1 = options['Numbers of elements along arms']
        elementsCount2 = 2
        elementsCount3 = 1
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = findOrCreateFieldCoordinates(fm)

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

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
                elementIndex[nerve] = int(xProportion[nerve] * elementsCount1[armNumber[nerve]-1])
                xi1[nerve] = 1 if xProportion[nerve] == 1 else xProportion[nerve] * elementsCount1[armNumber[nerve]-1] - elementIndex[nerve]
                elementIndex[nerve] += 1 if xProportion[nerve] < 1 else 0
                j = 10

            allMarkers = { "Inferior cardiac nerve" : {"elementID": elementIndex['ICN']+2*elementsCount1[0], "xi": [xi1['ICN'], 0.0, 0.5]},
                            "Ventral ansa subclavia" : {"elementID": elementIndex['VA']+2*elementsCount1[0]+elementsCount1[1], "xi": [xi1['VA'], 1.0, 0.5]},
                            "Dorsal ansa subclavia" : {"elementID": elementIndex['DA']+2*(elementsCount1[0]+elementsCount1[1]), "xi": [xi1['DA'], 0.0, 0.5]},
                            "Cervical spinal nerve 8" : {"elementID": elementIndex['C8']+2*(elementsCount1[0]+elementsCount1[1])+elementsCount1[2], "xi": [xi1['C8'], 1.0, 0.5]},
                            "Thoracic spinal nerve 1" : {"elementID": elementIndex['T1'], "xi": [xi1['T1'], 0.0, 0.5]},
                            "Thoracic spinal nerve 2" : {"elementID": elementIndex['T2'], "xi": [xi1['T2'], 0.0, 0.5]},
                            "Thoracic spinal nerve 3" : {"elementID": elementIndex['T3'], "xi": [xi1['T3'], 0.0, 0.5]},
                            "Thoracic sympathetic nerve trunk" : {"elementID": elementIndex['TST'], "xi": [xi1['TST'], 1.0, 0.5]},
                           }

        # arm group annotations for user
        armNames, faceNames = getUnlabelledStellateArmNames(armCount)
        armGroups = [AnnotationGroup(region, armNames[0]),
                     AnnotationGroup(region, armNames[1]),
                     AnnotationGroup(region, armNames[2])]
        stellateGroup = AnnotationGroup(region, get_stellate_term("cervicothoracic ganglion"))
        annotationGroups = armGroups.copy()
        annotationGroups.append(stellateGroup)

        armMeshGroups = [a.getMeshGroup(mesh) for a in armGroups]
        stellateMeshGroup = stellateGroup.getMeshGroup(mesh)

        # Create nodes
        numNodesPerArm = [0]
        dipMultiplier = 1
        nodeIdentifier = 1
        minArmAngle = 2 * math.pi / armCount
        halfArmArcAngleRadians = minArmAngle/2
        xx = []
        xds1 = []
        xds2 = []
        x_in_nodes = []
        for na in range(armCount):
            elementsCount_i = [elementsCount1[na], elementsCount2, elementsCount3]
            x, ds1, ds2, nWheelEdge = createArm(halfArmArcAngleRadians, elementLengths, elementLengthCentral, elementsCount_i, dipMultiplier, armCount, na)
            for ix in range(len(x)):
                if na == 0 or ix not in nWheelEdge:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[ix])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ds1[ix])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ds2[ix])

                    nodeIdentifier += 1
                    x_in_nodes.append(x[ix])
            numNodesPerArm.append(len(x))
            xx.append(x)
            xds1.append(ds1)
            xds2.append(ds2)

        # Create elements
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = bicubichermitelinear.createEftNoCrossDerivatives() #createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        cumNumNodesPerArm = [sum(numNodesPerArm[:i + 1]) for i in range(len(numNodesPerArm))]
        nCentre = [elementsCount1[0]+1, int(numNodesPerArm[1]/2) + elementsCount1[0]+1]
        for na in range(armCount):
            for e3 in range(elementsCount3):
                for e2 in range(elementsCount2):
                    for e1 in range(elementsCount1[na]):

                        ### NODES ###
                        no2 = (elementsCount1[na] + 1)
                        no3 = (elementsCount2 + 1) * no2 - 2
                        offset = (cumNumNodesPerArm[na])
                        bni = e3 * no3 + e2 * no2 + e1 + 1 + offset
                        if e2 == 0:
                            if e1 == 0 and na > 0: # and na < armCount -1: # wheelSouth
                                nWh = cumNumNodesPerArm[na - 1] + (2 * elementsCount1[na - 1]) + 2
                                nplUq = int(numNodesPerArm[na+1]/2) - elementsCount1[na] # unused nodes at centre and shared edge
                                npl = int(numNodesPerArm[na+1]/2)  #  nodes at centre and shared edge
                                if na < armCount-1:
                                    cn = cumNumNodesPerArm[na] + elementsCount1[na]-2
                                    no2 = cumNumNodesPerArm[na]
                                    em = elementsCount1[na]
                                    nwPrev = [nWh, nWh + int(numNodesPerArm[na]/2)] # previous arm's edge, depends on armCount.
                                    nodeIdentifiers = [nwPrev[0], no2 + 1, nCentre[0], no2 + em,
                                                       nwPrev[1], no2 + em - 1 + nplUq,
                                                       nCentre[1], bni + (4 * em) - 2]
                                else:
                                    nplPrev = int(numNodesPerArm[na]/2) - 2
                                    no2 = elementsCount1[na]-1
                                    no3 = int(numNodesPerArm[na+1]/2) - 3
                                    nwPrev = [cumNumNodesPerArm[na-1] + 2*(elementsCount1[na-1]),
                                              cumNumNodesPerArm[na-1] + 2*(elementsCount1[na-1]) + nplPrev]
                                    start = cumNumNodesPerArm[na] - 3
                                    nodeIdentifiers = [nwPrev[0], start,
                                                       nCentre[0], start + no2,
                                                       nwPrev[1],  start + no3,
                                                       nCentre[1], start + no2 +no3]
                            elif e1 == elementsCount1[na] - 1: # armEnd, south
                                if na == 0:
                                    nodeIdentifiers = [bni, bni + no2 - 1, bni + no2, bni + no3,
                                                       bni + no2 + no3 - 1, bni + no2 + no3]
                                else:
                                    no3 = armCount*elementsCount1[na] - 1
                                    no2 = elementsCount1[na]
                                    if na > 1:
                                        bni -= 4
                                        no3 -= 1
                                    nodeIdentifiers = [bni-1, bni + no2-2, bni + no2 - 1,
                                                       bni + no3 - 1, bni + no2 -2 + no3, bni + no2 + no3 - 1]
                            elif na > 0 and e1 > 0: #  [na=1+, e1=1+, e2=0] for len=3+
                                bni -= 1 + ((armCount+1)*(na-1))
                                no2 = elementsCount1[na]
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
                                    npl = int(numNodesPerArm[na+1] / 2) - 2
                                    no2 = elementsCount1[na]
                                    nodeIdentifiers = [nCentre[0], bni + 1, bni + no2 + 1, bni + no2 + 2,
                                                       nCentre[1], bni + npl + 1, bni + npl + no2 + 1, bni + npl + no2 + 2]
                                else: # last arm
                                    bni = cumNumNodesPerArm[na] - 2 - (armCount - elementsCount1[na])
                                    nodeIdentifiers = [nCentre[0], bni + 1, 1, bni + no2,
                                                       nCentre[1], bni + no3 - 2,
                                                       int(numNodesPerArm[1]/2)+1, bni + no2 + no3 - armCount]
                            elif e1 == elementsCount1[na] - 1: # armEnd north
                                if na > 0:
                                    no2 = elementsCount1[na]
                                    nplUq = int(numNodesPerArm[na + 1] / 2) - 2
                                    if na > 1:
                                        adj = na - 1
                                        bni -= armCount*na + (armCount-elementsCount1[na]) + 1
                                        if elementsCount1[na] < 3:
                                           bni += 1
                                        if elementsCount1[na]>3:
                                            bni -= elementsCount1[na] - 3
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
                                k = armCount*elementsCount1[na] - na
                                nodeIdentifiers = [bni, bni+1, bni + no2, bni + no2 + 1,
                                                   bni + k, bni + k + 1,
                                                   bni + no2 + k, bni + no2 + k + 1]
                            else:
                                nodeIdentifiers = [bni - 1, bni, bni + no2 - 1, bni + no2,
                                                   bni + no3 - 1, bni + no3,
                                                   bni + no2 + no3 - 1, bni + no2 + no3]

                        if e1 == 0:  # wheel
                            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            if armCount == 3:
                                if e2 == 0:
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
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
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
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                elif na == 1:
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
                                    remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    if e2 == 0:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, ns, Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])

                        elif e1 < (elementsCount1[na] - 1):
                            eft1 = eft
                            elementtemplate1 = elementtemplate
                            scalefactors = None
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

        ############################
        # annotation fiducial points
        ############################
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
        stellateGroup = getAnnotationGroupForTerm(annotationGroups, get_stellate_term("cervicothoracic ganglion"))
        armNames, faceNames = getUnlabelledStellateArmNames(armCount)
        arm1Group = getAnnotationGroupForTerm(annotationGroups, armNames[0])
        arm2Group = getAnnotationGroupForTerm(annotationGroups, armNames[1])
        arm3Group = getAnnotationGroupForTerm(annotationGroups, armNames[2])

        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
        is_exterior_face_xi2_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))

        is_arm1 = arm1Group.getFieldElementGroup(mesh2d)
        is_arm2 = arm2Group.getFieldElementGroup(mesh2d)
        is_arm3 = arm3Group.getFieldElementGroup(mesh2d)
        is_left_face = fm.createFieldOr(
            fm.createFieldAnd(is_arm2, is_exterior_face_xi2_1),
            fm.createFieldAnd(is_arm3, is_exterior_face_xi2_0))
        leftStellateGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, faceNames[1])
        leftStellateGroup.getMeshGroup(mesh2d).addElementsConditional(is_left_face)

        is_top_face = fm.createFieldOr(
            fm.createFieldAnd(is_arm1, is_exterior_face_xi2_1),
            fm.createFieldAnd(is_arm2, is_exterior_face_xi2_0))
        topStellateGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, faceNames[0])
        topStellateGroup.getMeshGroup(mesh2d).addElementsConditional(is_top_face)

        is_bottom_face = fm.createFieldOr(
            fm.createFieldAnd(is_arm1, is_exterior_face_xi2_0),
            fm.createFieldAnd(is_arm3, is_exterior_face_xi2_1))
        bottomStellateGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, faceNames[-1])
        bottomStellateGroup.getMeshGroup(mesh2d).addElementsConditional(is_bottom_face)

def getUnlabelledStellateArmNames(armCount):
    armNames = [("stellate arm %d"%(d), None) for d in range(1,armCount+1)]
    faceNames = [("stellate face %d-%d"%(d, d+1), None) for d in range(1,armCount+1)]
    faceNames[-1] = ("stellate face %d-1"%(armCount), None)
    return armNames, faceNames

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

