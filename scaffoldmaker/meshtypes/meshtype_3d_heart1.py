"""
Generates a 3-D heart model including ventricles, base and atria.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase1 import MeshType_3d_heartventriclesbase1
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils import zinc_utils
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heart1(Scaffold_base):
    '''
    Generates a 3-D heart model including ventricles, base and atria.
    '''

    @staticmethod
    def getName():
        return '3D Heart 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Unit Human 1',
            'Unit Mouse 1',
            'Unit Pig 1',
            'Unit Rat 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = MeshType_3d_heartatria1.getDefaultOptions(parameterSetName)
        # ventricles overrides some atria options, so update from it:
        optionsVentricles = MeshType_3d_heartventriclesbase1.getDefaultOptions(parameterSetName)
        options.update(optionsVentricles)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventriclesbase1.getOrderedOptionNames()
        optionNamesAtria = MeshType_3d_heartatria1.getOrderedOptionNames()
        # insert numbers of elements in atria in initial group
        for optionName in [
            'Number of elements around atrial free wall',
            'Number of elements around atrial septum',
            'Number of elements up atria',
            'Number of elements inlet']:
            optionNames.insert(6, optionName)
            optionNamesAtria.remove(optionName)
        # remove dependent or repeated options in atria1
        optionNamesAtria.remove('Aorta outer plus diameter')
        for optionName in optionNames:
            if optionName in optionNamesAtria:
                optionNamesAtria.remove(optionName)
        # add remaining atria options
        optionNames += optionNamesAtria
        # want refinement options last
        for optionName in [
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        '''
        dependentChanges = MeshType_3d_heartventriclesbase1.checkOptions(options) \
            or MeshType_3d_heartatria1.checkOptions(options)
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']
        elementsCountAroundAtrialFreeWall = options['Number of elements around atrial free wall']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundAtria = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
        useCrossDerivatives = False

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = zinc_utils.getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        # generate heartventriclesbase2 model and put atria2 on it
        annotationGroups = MeshType_3d_heartventriclesbase1.generateBaseMesh(region, options)
        annotationGroups += MeshType_3d_heartatria1.generateBaseMesh(region, options)
        lFibrousRingGroup = AnnotationGroup(region, 'left fibrous ring', FMANumber = 77124, lyphID = 'Lyph ID unknown')
        rFibrousRingGroup = AnnotationGroup(region, 'right fibrous ring', FMANumber = 77125, lyphID = 'Lyph ID unknown')
        annotationGroups += [ lFibrousRingGroup, rFibrousRingGroup ]

        # annotation points
        dataCoordinates = zinc_utils.getOrCreateCoordinateField(fm, 'data_coordinates')
        dataLabel = zinc_utils.getOrCreateLabelField(fm, 'data_label')
        dataElementXi = zinc_utils.getOrCreateElementXiField(fm, 'data_element_xi')

        datapoints = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        datapointTemplateInternal = datapoints.createNodetemplate()
        datapointTemplateInternal.defineField(dataCoordinates)
        datapointTemplateInternal.defineField(dataLabel)
        datapointTemplateInternal.defineField(dataElementXi)

        ##############
        # Create nodes
        ##############

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # discover left and right fibrous ring nodes from ventricles and atria
        # because nodes are iterated in identifier order, the lowest and first are on the lv outlet cfb, right and left on lower outer layers
        lavInnerNodeId = [ [], [] ]  # [n2][n1]
        lavOuterNodeId = [ [], [] ]  # [n2][n1]
        iter = lFibrousRingGroup.getNodesetGroup(nodes).createNodeiterator()
        cfbNodeId = iter.next().getIdentifier()
        cfbLeftNodeId = iter.next().getIdentifier()
        for n1 in range(elementsCountAroundAtria):
            lavInnerNodeId[0].append(iter.next().getIdentifier())
        lavOuterNodeId[0].append(cfbNodeId)
        lavOuterNodeId[0].append(cfbLeftNodeId)
        for n1 in range(elementsCountAroundAtrialFreeWall - 1):
            lavOuterNodeId[0].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            lavOuterNodeId[0].append(None)
        for n1 in range(elementsCountAroundAtria):
            lavInnerNodeId[1].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialFreeWall + 1):
            lavOuterNodeId[1].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            lavOuterNodeId[1].append(None)
        ravInnerNodeId = [ [], [] ]  # [n2][n1]
        ravOuterNodeId = [ [], [] ]  # [n2][n1]
        iter = rFibrousRingGroup.getNodesetGroup(nodes).createNodeiterator()
        cfbNodeId = iter.next().getIdentifier()
        cfbRightNodeId = iter.next().getIdentifier()
        for n1 in range(elementsCountAroundAtria):
            ravInnerNodeId[0].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            ravOuterNodeId[0].append(None)
        for n1 in range(elementsCountAroundAtrialFreeWall - 1):
            ravOuterNodeId[0].append(iter.next().getIdentifier())
        ravOuterNodeId[0].append(cfbRightNodeId)
        ravOuterNodeId[0].append(cfbNodeId)
        for n1 in range(elementsCountAroundAtria):
            ravInnerNodeId[1].append(iter.next().getIdentifier())
        for n1 in range(elementsCountAroundAtrialSeptum - 1):
            ravOuterNodeId[1].append(None)
        cfbUpperNodeId = iter.next().getIdentifier()  # cfb from left will be first
        for n1 in range(elementsCountAroundAtrialFreeWall):
            ravOuterNodeId[1].append(iter.next().getIdentifier())
        ravOuterNodeId[1].append(cfbUpperNodeId)

        #for n2 in range(2):
        #    print('n2', n2)
        #    print('lavInnerNodeId', lavInnerNodeId[n2])
        #    print('lavOuterNodeId', lavOuterNodeId[n2])
        #    print('ravInnerNodeId', ravInnerNodeId[n2])
        #    print('ravOuterNodeId', ravOuterNodeId[n2])

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(3)

        lFibrousRingMeshGroup = lFibrousRingGroup.getMeshGroup(mesh)
        rFibrousRingMeshGroup = rFibrousRingGroup.getMeshGroup(mesh)

        elementIdentifier = startElementIdentifier = zinc_utils.getMaximumElementIdentifier(mesh) + 1

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # create fibrous ring elements

        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives, linearAxis = 2, d_ds1 = Node.VALUE_LABEL_D_DS1, d_ds2 = Node.VALUE_LABEL_D_DS3)
        eftFibrousRing = bicubichermitelinear.createEftBasic()

        # left fibrous ring, starting at crux / collapsed posterior interatrial sulcus
        cruxElementId = None
        for e in range(-1, elementsCountAroundAtrialFreeWall):
            eft1 = eftFibrousRing
            n1 = e
            nids = [
                lavInnerNodeId[0][n1], lavInnerNodeId[0][n1 + 1], lavInnerNodeId[1][n1], lavInnerNodeId[1][n1 + 1],
                lavOuterNodeId[0][n1], lavOuterNodeId[0][n1 + 1], lavOuterNodeId[1][n1], lavOuterNodeId[1][n1 + 1]]
            scalefactors = None
            meshGroups = [ lFibrousRingMeshGroup ]

            if e == -1:
                # interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = ravInnerNodeId[0][-1]
                nids[2] = ravInnerNodeId[1][-1]
                nids.pop(6)
                nids.pop(4)
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                # reverse d3 on cfb:
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [0] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ rFibrousRingMeshGroup ]
            elif e == 0:
                # general linear map d3 adjacent to collapsed sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                # reverse d1, d3 on cfb, left cfb:
                scaleEftNodeValueLabels(eft1, [ 6 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
            elif e == 1:
                # reverse d1, d3 on left cfb:
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == (elementsCountAroundAtrialFreeWall - 1):
                # general linear map d3 adjacent to collapsed sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element fibrous ring left', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # right fibrous ring, starting at crux / collapsed posterior interatrial sulcus
        ran1FreeWallStart = elementsCountAroundAtrialSeptum - 1
        for e in range(-1, elementsCountAroundAtrialFreeWall):
            eft1 = eftFibrousRing
            n1 = ran1FreeWallStart + e
            nids = [
                ravInnerNodeId[0][n1], ravInnerNodeId[0][n1 + 1], ravInnerNodeId[1][n1], ravInnerNodeId[1][n1 + 1],
                ravOuterNodeId[0][n1], ravOuterNodeId[0][n1 + 1], ravOuterNodeId[1][n1], ravOuterNodeId[1][n1 + 1]]
            scalefactors = None
            meshGroups = [ rFibrousRingMeshGroup ]

            if e == -1:
                # interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = lavInnerNodeId[0][elementsCountAroundAtrialFreeWall]
                nids[2] = lavInnerNodeId[1][elementsCountAroundAtrialFreeWall]
                nids.pop(6)
                nids.pop(4)
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                meshGroups += [ lFibrousRingMeshGroup ]
                cruxElementId = elementIdentifier
            elif e == 0:
                # general linear map d3 adjacent to collapsed crux/posterior sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
            elif e == (elementsCountAroundAtrialFreeWall - 2):
                # reverse d1, d3 on right cfb:
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1]) ])
            elif e == (elementsCountAroundAtrialFreeWall - 1):
                # general linear map d3 adjacent to collapsed cfb/anterior sulcus
                eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                # reverse d1, d3 on right cfb, cfb:
                scaleEftNodeValueLabels(eft1, [ 5 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1]) ])
                remapEftNodeValueLabel(eft1, [ 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element fibrous ring right', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # fibrous ring septum:
        meshGroups = [ lFibrousRingMeshGroup, rFibrousRingMeshGroup ]
        for e in range(elementsCountAroundAtrialSeptum):
            eft1 = eftFibrousRing
            nlm = elementsCountAroundAtrialFreeWall + e
            nlp = (nlm + 1)%elementsCountAroundAtria
            nrm = elementsCountAroundAtrialSeptum - 1 - e
            nrp = nrm - 1
            nids = [ lavInnerNodeId[0][nlm], lavInnerNodeId[0][nlp], lavInnerNodeId[1][nlm], lavInnerNodeId[1][nlp],
                     ravInnerNodeId[0][nrm], ravInnerNodeId[0][nrp], ravInnerNodeId[1][nrm], ravInnerNodeId[1][nrp] ]

            eft1 = bicubichermitelinear.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            if e == 0:
                # general linear map d3 adjacent to collapsed posterior interventricular sulcus
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 6, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            elif e == (elementsCountAroundAtrialSeptum - 1):
                # general linear map d3 adjacent to cfb
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 5, 7 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            else:
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3 ], [ 1 ])

            result = elementtemplate1.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = 7
            #print('create element fibrous ring septum', elementIdentifier, result, result2, result3, nids)
            elementIdentifier += 1

            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # annotation points
        cruxElement = mesh.findElementByIdentifier(cruxElementId)
        cruxXi = [ 0.0, 0.5, 1.0 ]
        cache.setMeshLocation(cruxElement, cruxXi)
        result, cruxCoordinates = coordinates.evaluateReal(cache, 3)
        datapoint = datapoints.createNode(-1, datapointTemplateInternal)
        cache.setNode(datapoint)
        dataCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cruxCoordinates)
        dataLabel.assignString(cache, 'crux')
        dataElementXi.assignMeshLocation(cache, cruxElement, cruxXi)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        elementsCountAroundAtrialFreeWall = options['Number of elements around atrial free wall']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundAtria = elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        MeshType_3d_heartventriclesbase1.refineMesh(meshrefinement, options)
        MeshType_3d_heartatria1.refineMesh(meshrefinement, options)
        element = meshrefinement._sourceElementiterator.next()
        startFibrousRingElementIdentifier = element.getIdentifier()
        lastFibrousRingElementIdentifier = startFibrousRingElementIdentifier + 2*elementsCountAroundAtrialFreeWall + elementsCountAroundAtrialSeptum + 1
        i = -1
        while element.isValid():
            numberInXi1 = refineElementsCountSurface
            numberInXi2 = 1  # since fibrous ring is thin
            numberInXi3 = refineElementsCountThroughWall
            elementIdentifier = element.getIdentifier()
            if i in [ -1, elementsCountAroundAtrialFreeWall ]:
                numberInXi1 = refineElementsCountThroughWall
            meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            if elementIdentifier == lastFibrousRingElementIdentifier:
                return  # finish on last so can continue elsewhere
            element = meshrefinement._sourceElementiterator.next()
            i += 1

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()
