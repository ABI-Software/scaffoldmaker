"""
Generates a 3-D esophagus mesh along the central line,
with variable numbers of elements around, along and through
wall, with variable radius and thickness along.
"""

import copy
import math
from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldStoredString, \
    findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.esophagus_terms import get_esophagus_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import geometry
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class MeshType_3d_esophagus1(Scaffold_base):
    """
    Generates a 3-D esophagus mesh with variable numbers of elements around, along the central line, and through wall.
    The esophagus is created by a function that generates an elliptical tube segment and uses tubemesh to map the
    segment along a central path profile.
    """

    centralPathDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 4
                },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  7.5, 289.8,  0.0 ], [  2.8,  -33.6,  0.0 ], [  9.9,  0.8, 0.0 ], [  0.1, -1.6, 0.0 ], [ 0.0, 0.0, 10.0 ], [ 0.0, 0.0, 0.0] ] , 
                [ [  8.9, 242.5,  0.0 ], [ -2.6,  -87.5,  0.0 ], [ 10.0, -0.3, 0.0 ], [  0.1, -0.6, 0.0 ], [ 0.0, 0.0, 10.0 ], [ 0.0, 0.0, 0.0] ] , 
                [ [ -2.5, 115.4,  0.0 ], [  3.9, -102.0,  0.0 ], [ 10.0,  0.4, 0.0 ], [ -0.3,  1.9, 0.0 ], [ 0.0, 0.0, 10.0 ], [ 0.0, 0.0, 0.0] ] , 
                [ [ 10.1,  40.3,  0.0 ], [ 17.8,  -57.9,  0.0 ], [  9.6,  2.9, 0.0 ], [ -1.3,  3.2, 0.0 ], [ 0.0, 0.0, 10.0 ], [ 0.0, 0.0, 0.0] ] , 
                [ [ 28.6,  -0.1,  0.0 ], [ 18.5,  -22.1,  0.0 ], [  7.7,  6.5, 0.0 ], [ -2.4,  3.9, 0.0 ], [ 0.0, 0.0, 10.0 ], [ 0.0, 0.0, 0.0] ] ] ),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_esophagus_term('cervical part of esophagus')[0],
                    'ontId': get_esophagus_term('cervical part of esophagus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2-3',
                    'name': get_esophagus_term('thoracic part of esophagus')[0],
                    'ontId': get_esophagus_term('thoracic part of esophagus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4',
                    'name': get_esophagus_term('abdominal part of esophagus')[0],
                    'ontId': get_esophagus_term('abdominal part of esophagus')[1]
                }]
            })
        }

    @staticmethod
    def getName():
        return '3D Esophagus 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around': 8,
            'Number of elements along': 12,
            'Number of elements through wall': 4,
            'Wall thickness': 5.0,  # 3.0mm when relaxed
            'Mucosa relative thickness': 0.35,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.25,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around',
            'Number of elements along',
            'Number of elements through wall',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
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
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
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
        if options['Number of elements through wall'] != (1 or 4):
            options['Number of elements through wall'] = 4
        for key in [
            'Number of elements around',
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Wall thickness'] < 0.0:
            options['Wall thickness'] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        centralPath = options['Central path']
        elementsCountAround = options['Number of elements around']
        elementsCountAlong = options['Number of elements along']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        mucosaRelThickness = options['Mucosa relative thickness']
        submucosaRelThickness = options['Submucosa relative thickness']
        circularRelThickness = options['Circular muscle layer relative thickness']
        longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        esophagusTermsAlong = [None, 'cervical part of esophagus',
                               'thoracic part of esophagus', 'abdominal part of esophagus']
        arcLengthOfGroupsAlong = []
        for i in range(len(esophagusTermsAlong)):
            tmpRegion = region.createRegion()
            centralPath.generate(tmpRegion)
            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                groupName=esophagusTermsAlong[i])
            arcLength = 0.0
            for e in range(len(cxGroup) - 1):
                arcLength += interp.getCubicHermiteArcLength(cxGroup[e], cd1Group[e],
                                                             cxGroup[e + 1], cd1Group[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)

            if i == 0:
                cx = cxGroup
                cd1 = cd1Group
                cd2 = cd2Group
                cd3 = cd3Group
                cd12 = cd12Group
                cd13 = cd13Group

            del tmpRegion

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        centralPathLength = arcLengthOfGroupsAlong[0]
        elementAlongLength = centralPathLength / elementsCountAlong

        elementsCountAlongGroups = []
        groupLength = 0.0
        e = 0
        elementsCount = 1
        length = elementAlongLength
        for i in range(1, len(esophagusTermsAlong)):
            groupLength += arcLengthOfGroupsAlong[i]
            if e == elementsCountAlong - 2:
                elementsCount += 1
                elementsCountAlongGroups.append(elementsCount)
            else:
                while length < groupLength:
                    elementsCount += 1
                    e += 1
                    length += elementAlongLength

                # check which end is grouplength closer to
                distToUpperEnd = abs(length - groupLength)
                distToLowerEnd = abs(groupLength - (length - elementsCountAlong))
                if distToLowerEnd < distToUpperEnd:
                    elementsCount -= 1
                    elementsCountAlongGroups.append(elementsCount)
                    e -= 1
                    length -= elementAlongLength
                else:
                    elementsCountAlongGroups.append(elementsCount)
            elementsCount = 0

        majorRadiusElementList = sd2
        minorRadiusElementList = sd3

        # Create annotation groups along esophagus
        esophagusGroup = AnnotationGroup(region, get_esophagus_term("esophagus"))
        cervicalGroup = AnnotationGroup(region, get_esophagus_term("cervical part of esophagus"))
        thoracicGroup = AnnotationGroup(region, get_esophagus_term("thoracic part of esophagus"))
        abdominalGroup = AnnotationGroup(region, get_esophagus_term("abdominal part of esophagus"))

        annotationGroupAlong = [[esophagusGroup, cervicalGroup],
                                [esophagusGroup, thoracicGroup],
                                [esophagusGroup, abdominalGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        annotationGroupsAround = []
        for i in range(elementsCountAround):
            annotationGroupsAround.append([])

        # Groups through wall
        longitudinalMuscleGroup = AnnotationGroup(region,
                                                  get_esophagus_term("esophagus smooth muscle longitudinal layer"))
        circularMuscleGroup = AnnotationGroup(region, get_esophagus_term("esophagus smooth muscle circular layer"))
        submucosaGroup = AnnotationGroup(region, get_esophagus_term("submucosa of esophagus"))
        mucosaGroup = AnnotationGroup(region, get_esophagus_term("esophagus mucosa"))

        if elementsCountThroughWall == 1:
            relativeThicknessList = [1.0]
            annotationGroupsThroughWall = [[]]
        else:
            relativeThicknessList = [mucosaRelThickness, submucosaRelThickness,
                                     circularRelThickness, longitudinalRelThickness]
            annotationGroupsThroughWall = [[mucosaGroup],
                                           [submucosaGroup],
                                           [circularMuscleGroup],
                                           [longitudinalMuscleGroup]]

        xToSample = []
        d1ToSample = []
        for n2 in range(elementsCountAlong + 1):
            # Create inner points
            cx = [0.0, 0.0, elementAlongLength * n2]
            axis1 = [vector.magnitude(majorRadiusElementList[n2]), 0.0, 0.0]
            axis2 = [0.0, vector.magnitude(minorRadiusElementList[n2]), 0.0]
            xInner, d1Inner = geometry.createEllipsePoints(cx, 2 * math.pi, axis1, axis2,
                                                           elementsCountAround, startRadians=0.0)
            xToSample += xInner
            d1ToSample += d1Inner

        d2ToSample = [[0.0, 0.0, elementAlongLength]] * (elementsCountAround * (elementsCountAlong+1))

        # Sample along length
        xInnerRaw = []
        d2InnerRaw = []
        xToWarp = []
        d1ToWarp = []
        d2ToWarp = []
        flatWidthList = []
        xiList = []

        for n1 in range(elementsCountAround):
            xForSamplingAlong = []
            d2ForSamplingAlong = []
            for n2 in range(elementsCountAlong + 1):
                idx = n2 * elementsCountAround + n1
                xForSamplingAlong.append(xToSample[idx])
                d2ForSamplingAlong.append(d2ToSample[idx])
            xSampled, d2Sampled = interp.sampleCubicHermiteCurves(xForSamplingAlong, d2ForSamplingAlong,
                                                                  elementsCountAlong, arcLengthDerivatives=True)[0:2]
            xInnerRaw.append(xSampled)
            d2InnerRaw.append(d2Sampled)

        # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
        for n2 in range(elementsCountAlong + 1):
            xAround = []
            d2Around = []

            for n1 in range(elementsCountAround):
                x = xInnerRaw[n1][n2]
                d2 = d2InnerRaw[n1][n2]
                xAround.append(x)
                d2Around.append(d2)

            d1Around = []
            for n1 in range(elementsCountAround):
                v1 = xAround[n1]
                v2 = xAround[(n1 + 1) % elementsCountAround]
                d1 = d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c * arcLengthAround for c in vector.normalise(d1)]
                d1Around.append(dx_ds1)
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)

            xToWarp += xAround
            d1ToWarp += d1Smoothed
            d2ToWarp += d2Around

            # Flat width and xi
            flatWidth = 0.0
            xiFace = []
            for n1 in range(elementsCountAround):
                v1 = xAround[n1]
                d1 = d1Smoothed[n1]
                v2 = xAround[(n1 + 1) % elementsCountAround]
                d2 = d1Smoothed[(n1 + 1) % elementsCountAround]
                flatWidth += interp.getCubicHermiteArcLength(v1, d1, v2, d2)
            flatWidthList.append(flatWidth)

            for n1 in range(elementsCountAround + 1):
                xi = 1.0 / elementsCountAround * n1
                xiFace.append(xi)
            xiList.append(xiFace)

        # Project reference point for warping onto central path
        sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
            tubemesh.getPlaneProjectionOnCentralPath(xToWarp, elementsCountAround, elementsCountAlong,
                                                     centralPathLength, sx, sd1, sd2, sd12)

        # Warp points
        segmentAxis = [0.0, 0.0, 1.0]
        closedProximalEnd = False

        innerRadiusAlong = []
        for n2 in range(elementsCountAlong + 1):
            firstNodeAlong = xToWarp[n2 * elementsCountAround]
            midptSegmentAxis = [0.0, 0.0, elementAlongLength * n2]
            radius = vector.magnitude(firstNodeAlong[c] - midptSegmentAxis[c] for c in range(3))
            innerRadiusAlong.append(radius)

        xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
            tubemesh.warpSegmentPoints(xToWarp, d1ToWarp, d2ToWarp, segmentAxis, sxRefList, sd1RefList,
                                       sd2ProjectedListRef, elementsCountAround, elementsCountAlong,
                                       zRefList, innerRadiusAlong, closedProximalEnd)

        # Create coordinates and derivatives
        transitElementList = [0]*elementsCountAround
        xList, d1List, d2List, d3List, curvatureList = \
            tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList,
                                             [wallThickness]*(elementsCountAlong+1), relativeThicknessList,
                                             elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                                             transitElementList)

        # Create flat coordinates
        xFlat, d1Flat, d2Flat = tubemesh.createFlatCoordinates(xiList, flatWidthList, length, wallThickness,
                                                               relativeThicknessList, elementsCountAround,
                                                               elementsCountAlong, elementsCountThroughWall,
                                                               transitElementList)

        # Create nodes and elements
        xOrgan = []
        d1Organ = []
        d2Organ = []
        nodeIdentifier, elementIdentifier, annotationGroups = \
            tubemesh.createNodesAndElements(region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat,
                                            xOrgan, d1Organ, d2Organ, None, elementsCountAround, elementsCountAlong,
                                            elementsCountThroughWall, annotationGroupsAround, annotationGroupsAlong,
                                            annotationGroupsThroughWall, firstNodeIdentifier, firstElementIdentifier,
                                            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd)

        # annotation fiducial points
        fm = region.getFieldmodule()
        fm.beginChange()
        mesh = fm.findMeshByDimension(3)
        cache = fm.createFieldcache()

        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        markerNames = ["proximodorsal midpoint on serosa of upper esophageal sphincter",
                       "proximoventral midpoint on serosa of upper esophageal sphincter",
                       "distal point of lower esophageal sphincter serosa on the greater curvature of stomach",
                       "distal point of lower esophageal sphincter serosa on the lesser curvature of stomach"]

        totalElements = elementIdentifier
        radPerElementAround = math.pi * 2.0 / elementsCountAround
        elementAroundHalfPi = int(0.25 * elementsCountAround)
        xi1HalfPi = (math.pi * 0.5 - radPerElementAround * elementAroundHalfPi)/radPerElementAround
        elementAroundPi = int(0.5 * elementsCountAround)
        xi1Pi = (math.pi - radPerElementAround * elementAroundPi)/radPerElementAround

        markerElementIdentifiers = [elementsCountAround * elementsCountThroughWall - elementAroundHalfPi,
                                    elementAroundHalfPi + 1 + elementsCountAround * (elementsCountThroughWall - 1),
                                    totalElements - elementsCountAround,
                                    totalElements - elementsCountAround + elementAroundPi]

        markerXis = [[1.0 - xi1HalfPi, 0.0, 1.0],
                     [xi1HalfPi, 0.0, 1.0],
                     [0.0, 1.0, 1.0],
                     [xi1Pi, 1.0, 1.0]]

        for n in range(len(markerNames)):
            markerGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_esophagus_term(markerNames[n]))
            markerElement = mesh.findElementByIdentifier(markerElementIdentifiers[n])
            markerXi = markerXis[n]
            cache.setMeshLocation(markerElement, markerXi)
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            nodeIdentifier += 1
            cache.setNode(markerPoint)
            markerName.assignString(cache, markerGroup.getName())
            markerLocation.assignMeshLocation(cache, markerElement, markerXi)
            for group in [esophagusGroup, markerGroup]:
                group.getNodesetGroup(nodes).addNode(markerPoint)

        fm.endChange()

        return annotationGroups

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
        mesh2d = fm.findMeshByDimension(2)

        esophagusGroup = getAnnotationGroupForTerm(annotationGroups, get_esophagus_term("esophagus"))
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_esophagus = esophagusGroup.getFieldElementGroup(mesh2d)
        is_serosa = fm.createFieldAnd(is_esophagus, is_exterior_face_xi3_1)
        serosa = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_esophagus_term("serosa of esophagus"))
        serosa.getMeshGroup(mesh2d).addElementsConditional(is_serosa)
