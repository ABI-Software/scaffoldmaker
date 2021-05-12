"""
Brainstem mesh using a tapered cylinder
"""

from __future__ import division
import copy
from opencmiss.zinc.element import Element
from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.field import Field, findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.brainstem_terms import get_brainstem_annotation_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, Tapered, ConeBaseProgression, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage


class MeshType_3d_brainstem1(Scaffold_base):
    """
    Generates a tapered cylinder for the brainstem based on solid cylinder mesh, with variable numbers of elements in major, minor and length directions. Regions of the brainstem are annotated.
    """

    centralPathDefaultScaffoldPackages = {
        'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues( # dimensional.
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, -1.0, 5.0], [0.0, 0.0, -4.5], [5.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 2.4, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, -1.0, 0.5], [0.0, 0.0, -4.5], [6.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 4, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, -1.0, -4.0], [0.0, 0.0, -4.5], [7.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 4.5, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, -1.0, -8.5], [0.0, 0.0, -4.5], [8.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 5.5, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, -1.0, -13.0], [0.0, 0.0, -4.5], [9.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 0.0]]
            # 'meshEdits': exnodeStringFromNodeValues( # dimensional.
            #     [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            #      Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
            #         [[0.0, -1.0, 5.0], [0.0, 0.0, -4.5], [5.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 2.4, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, -1.0, 0.5], [0.0, 0.0, -4.5], [5.8, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 4, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, -1.0, -4.0], [0.0, 0.0, -4.5], [8.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 4.5, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, -1.0, -8.5], [0.0, 0.0, -4.5], [10.2, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 5.5, 0.0], [0.0, 0.0, 0.0]],
            #         [[0.0, -1.0, -13.0], [0.0, 0.0, -4.5], [9.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Brainstem 1'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 6,
            'Number of elements across minor': 6,
            'Number of elements along': 12,
            'Lower half': False,
            'Taper major increment': 0.2,
            'Taper minor increment': 0.1,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major and minor': 1,
            'Refine number of elements along': 1
        }
        return options

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


    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements along',
            'Lower half',
            'Refine',
            'Refine number of elements across major and minor',
            'Refine number of elements along'
        ]


    @classmethod
    def checkOptions(cls, options):

        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False
        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1

        if options['Number of elements along'] < 2:
            options['Number of elements along'] = 2

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        fm = region.getFieldmodule()
        mesh = fm.findMeshByDimension(3)
        coordinates = findOrCreateFieldCoordinates(fm)

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAlong = options['Number of elements along']

        addNuclearGroupAnnotation = False

        elemPerLayer = ((elementsCountAcrossMajor - 2) * elementsCountAcrossMinor) + (
                    2 * (elementsCountAcrossMinor - 2))
        if addNuclearGroupAnnotation:
            NAGroup = AnnotationGroup(region, get_brainstem_annotation_term('nucleus ambiguus'))
            NTSGroup = AnnotationGroup(region, get_brainstem_annotation_term('nucleus tractus solitarius'))
            rapheGroup = AnnotationGroup(region, get_brainstem_annotation_term('raphe complex'))
            mPBGroup = AnnotationGroup(region, get_brainstem_annotation_term('medial parabrachial nucleus'))
            lPBGroup = AnnotationGroup(region, get_brainstem_annotation_term('lateral parabrachial nucleus'))
            KFNGroup = AnnotationGroup(region, get_brainstem_annotation_term('Kölliker-Fuse nucleus'))
            RTNGroup = AnnotationGroup(region, get_brainstem_annotation_term('retrotrapezoid nucleus'))
            botzGroup = AnnotationGroup(region, get_brainstem_annotation_term('Bötzinger complex'))
            prebotzGroup = AnnotationGroup(region, get_brainstem_annotation_term('pre-Bötzinger complex'))
            rVRGGroup = AnnotationGroup(region, get_brainstem_annotation_term('rostral ventral respiratory group'))
            cVRGGroup = AnnotationGroup(region, get_brainstem_annotation_term('caudal ventral respiratory group'))
            DRGGroup = AnnotationGroup(region, get_brainstem_annotation_term('dorsal respiratory group'))
            xyzProportion = {}
            xyzProportion['NA'] = [[0.3,0.5,0.3], [0.7,0.5,0.3]]
            xyzProportion['NTS'] = [[0.4,0.8,0.3], [0.6,0.8,0.3]]
            xyzProportion['raphe'] = [[0.45,0.1,0.7], [0.55,0.1,0.7]]
            xyzProportion['mPB'] = [[0.2,0.65,0.9], [0.8,0.65,0.9]]
            xyzProportion['lPB'] = [[0.15,0.65,0.9], [0.85,0.65,0.9]]
            xyzProportion['KFN'] = [[0.1,0.65,0.9], [0.9,0.65,0.9]]
            xyzProportion['RTN'] = [[0.1,0.1,0.7], [0.9,0.1,0.7]]
            xyzProportion['botz'] = [[0.1,0.5,0.6], [0.9,0.5,0.6]]
            xyzProportion['prebotz'] = [[0.1,0.5,0.5], [0.9,0.5,0.5]]
            xyzProportion['rVRG'] = [[0.1,0.5,0.4], [0.9,0.5,0.4]]
            xyzProportion['cVRG'] = [[0.1,0.5,0.3], [0.9,0.5,0.3]]
            xyzProportion['DRG'] = [[0.4,0.8,0.3], [0.6,0.8,0.3]] # inside NTS
            nucAbbrev = [key for key in xyzProportion.keys() if isinstance(xyzProportion[key], list)]
            elementIndex = {}
            xi1 = {}
            for nuc in nucAbbrev:
                elementIndex[nuc] = []
                for prow in xyzProportion[nuc]:
                    # element numbering is along minor, descending x, and ends have 2 less elements
                    elemsPerColumn = [elementsCountAcrossMinor-2]+ [elementsCountAcrossMinor]*(elementsCountAcrossMajor-2)+[elementsCountAcrossMinor-2]
                    cumElemsPerColumn = [0] + [sum(elemsPerColumn[:i + 1]) for i in range(len(elemsPerColumn))]
                    for i in range(3):
                        prow[i] -= 0.1 if prow[i] == 1 else 0
                    [nx, ny, nz] = [int(prow[0] * elementsCountAcrossMajor),
                                    int(prow[1] * elementsCountAcrossMinor),
                                    int(prow[2] * elementsCountAlong)]
                    if (nx == 0 or nx >= elementsCountAcrossMajor-1):# and ny > (elementsCountAcrossMinor/2):
                        ny = int(ny/2)
                    nx -= 1 if nx == elementsCountAcrossMajor else 0
                    nz -= 1 if nz == elementsCountAlong else 0
                    elementIndex[nuc].append(ny + cumElemsPerColumn[nx] + (elemPerLayer*nz) + 1)
                    xi1[nuc] = 1 if prow == 1 else prow

            # order the same as nucAbbrev
            annotationGroups = [ NAGroup, NTSGroup, rapheGroup, mPBGroup, lPBGroup, KFNGroup, RTNGroup, botzGroup, prebotzGroup, rVRGGroup, cVRGGroup, DRGGroup ]
        else:
            annotationGroups = []
        brainstemGroup = AnnotationGroup(region, get_brainstem_annotation_term('brainstem'))
        midbrainGroup = AnnotationGroup(region, get_brainstem_annotation_term('midbrain'))
        ponsGroup = AnnotationGroup(region, get_brainstem_annotation_term('pons'))
        medullaGroup = AnnotationGroup(region, get_brainstem_annotation_term('medulla oblongata'))
        brainstemMeshGroup = brainstemGroup.getMeshGroup(mesh)
        midbrainMeshGroup = midbrainGroup.getMeshGroup(mesh)
        ponsMeshGroup = ponsGroup.getMeshGroup(mesh)
        medullaMeshGroup = medullaGroup.getMeshGroup(mesh)
        allMeshGroups = [igroup.getMeshGroup(mesh) for igroup in annotationGroups]
        annotationGroups += [brainstemGroup, midbrainGroup, ponsGroup, medullaGroup]

        #######################
        # CREATE MAIN BODY MESH
        #######################
        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        taperedParams = Tapered(majorRatio=options['Taper major increment'],
                                minorRatio=options['Taper minor increment'],
                                majorProgressionMode=ConeBaseProgression.ARITHMETIC_PROGRESSION,
                                minorProgressionMode=ConeBaseProgression.ARITHMETIC_PROGRESSION)

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor,
                            centre=[0.0, 0.0, 0.0],
                            alongAxis=cylinderCentralPath.alongAxis[0], majorAxis=cylinderCentralPath.majorAxis[0],
                            minorRadius=cylinderCentralPath.minorRadii[0])

        cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
                            cylinderShape=cylinderShape,
                                 tapered = taperedParams,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        if addNuclearGroupAnnotation:
            for jn, nuc in enumerate(nucAbbrev):
                for elementIdentifier in elementIndex[nuc]:
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    allMeshGroups[jn].addElement(element)

        iRegionBoundaries = [int(7*elementsCountAlong/15),int(14*elementsCountAlong/15)]
        for elementIdentifier in range(1, mesh.getSize()+1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            brainstemMeshGroup.addElement(element)
            if elementIdentifier > (iRegionBoundaries[-1]*elemPerLayer):
                midbrainMeshGroup.addElement(element)
            elif elementIdentifier > (iRegionBoundaries[0]*elemPerLayer) and elementIdentifier <= (iRegionBoundaries[-1]*elemPerLayer):
                ponsMeshGroup.addElement(element)
            else:
                medullaMeshGroup.addElement(element)

        ##############################
        # point markers
        ##############################
        eIndexPM = {}
        xiPM = {}
        pointMarkers = {}
        eIndexPM['caudal-dorsal'] = int(elemPerLayer/2)
        eIndexPM['midRostralCaudal-dorsal'] = int(elemPerLayer / 2) + (elemPerLayer * int((elementsCountAlong/2)-1))
        eIndexPM['rostral-dorsal'] = (elemPerLayer*(elementsCountAlong-1)) + int(elemPerLayer/2)
        eIndexPM['caudal-ventral'] = int(elemPerLayer/2) - (elementsCountAcrossMinor-1)
        eIndexPM['midRostralCaudal-ventral'] = eIndexPM['midRostralCaudal-dorsal'] - (elementsCountAcrossMinor-1)
        eIndexPM['rostral-ventral'] = int((elemPerLayer*(elementsCountAlong-1)) + (int(elemPerLayer/2) - elementsCountAcrossMinor + 1))
        xiPM['caudal-ventral'] = [1.0, 0.0, 0.0]
        xiPM['caudal-dorsal'] = [1.0, 0.0, 1.0]
        xiPM['midRostralCaudal-ventral'] = [1.0, 1.0, 0.0]
        xiPM['midRostralCaudal-dorsal'] = [1.0, 1.0, 1.0]
        xiPM['rostral-ventral'] = [1.0, 1.0, 0.0]
        xiPM['rostral-dorsal'] = [1.0, 1.0, 1.0]
        for key in eIndexPM.keys():
            pointMarkers[key] = {"elementID": eIndexPM[key], "xi": xiPM[key]}
        emergentMarkers = createCranialNerveEmergentMarkers(region, mesh, "coordinates")
        pointMarkers.update(emergentMarkers)


        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        cache = fm.createFieldcache()
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        if pointMarkers:
            for key in pointMarkers:
                addMarker = {"name": key, "xi": pointMarkers[key]["xi"]}
                markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
                nodeIdentifier += 1
                cache.setNode(markerPoint)
                markerName.assignString(cache, addMarker["name"])
                elementID = pointMarkers[key]["elementID"]
                element = mesh.findElementByIdentifier(elementID)
                markerLocation.assignMeshLocation(cache, element, addMarker["xi"])

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major and minor']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the 1D mesh.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for ventral-level elements.
        New point annotation groups are appended to this list.
        """
        # create  groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi1 = fm.createFieldOr(
            fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)),
            fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)))
        is_exterior_face_xi3 = fm.createFieldOr(fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)), fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)))

        annoGroup = AnnotationGroup(region, get_brainstem_annotation_term('brainstem'))
        isGroup = annoGroup.getFieldElementGroup(mesh2d)
        is_face1 = fm.createFieldAnd(isGroup, is_exterior_face_xi1)
        is_face3 = fm.createFieldAnd(isGroup, is_exterior_face_xi3)
        is_face_ext = fm.createFieldOr(is_face1, is_face3)
        faceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, ("brainstem_exterior", None))
        faceGroup.getMeshGroup(mesh2d).addElementsConditional(is_face_ext)

        # external regions
        namelist = ['midbrain', 'pons', 'medulla oblongata']
        for subregion in namelist:
            subGroup = AnnotationGroup(region, get_brainstem_annotation_term(subregion))
            issub = subGroup.getFieldElementGroup(mesh2d)
            is_subface = fm.createFieldOr(fm.createFieldAnd(issub, is_exterior_face_xi1), fm.createFieldAnd(issub, is_exterior_face_xi3))
            subFaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, (subregion+'_exterior', None))
            subFaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_subface)

def createCranialNerveEmergentMarkers(region, mesh, coordinatesName):
    # create marker points for locations the cranial nerves emerge from brainstem mesh, based on the USF cat brainstem data.
    # return element xi
    # use findMeshLocation to find the elementxi in an arbitrary mesh of given number of elements.

    # brainstem_coordinates: the left-side nerves
    nerveDict = {'OCULOMOTOR_left':[-0.13912257342955267, -0.5345161733750351, -0.7374762051676923],
                 'TROCHLEAR_left':[-0.13148279719950992, 0.4218745504359067, -0.7375838988856348],
                 'TRIGEMINAL_left': [-0.7605971693047597, -0.4025791045292648, -0.6862730212268676],
                 'ABDUCENS_left': [-0.19517975766630574, -0.6252563181242173, -0.8205128205130072],
                 'FACIAL_left': [-0.5824675040481234, -0.3554448371502354, -0.24509655553058302],
                 'VESTIBULOCOCHLEAR_left': [-0.6147505791411602, -0.32790803815838, -0.24509655403515848],
                 'GLOSSOPHARYNGEAL_left': [-0.7307312460087607, -0.2576952819028721, -0.39215539053073717],
                 'VAGUS_left': [-0.6741855912315219, -0.25981298010131126, -0.24509655277992023],
                 'ACCESSORY_cranialRoot_left':[-0.6741855912315219, -0.25981298010131126, -0.24509655277992023],
                 'HYPOGLOSSAL_left': [-0.044776303107883636, -0.5027870527016534, -0.10510117079651562]
                 }

    rightDict = {}
    for key in nerveDict.keys():
        nerveName = key.split('_')[0]+'_right'
        xyz = [-1*nerveDict[key][0], nerveDict[key][1], nerveDict[key][2]]
        rightDict.update({nerveName:xyz})
    nerveDict.update(rightDict)

    nerveNames = list(nerveDict.keys())

    # add to data coordinates
    markerNameField = 'marker_name'
    fm = region.getFieldmodule()
    cache = fm.createFieldcache()
    datapoints = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    data_coordinates = findOrCreateFieldCoordinates(fm, "data_coordinates")
    markerName = findOrCreateFieldStoredString(fm, name="marker_name")
    dnodetemplate = datapoints.createNodetemplate()
    dnodetemplate.defineField(data_coordinates)
    dnodetemplate.setValueNumberOfVersions(data_coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    dnodetemplate.defineField(markerName)
    dnodeIdentifier = 1
    for nerveName in nerveNames:
        node = datapoints.createNode(dnodeIdentifier, dnodetemplate)
        cache.setNode(node)
        addEnd = nerveDict[nerveName].copy()
        data_coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, addEnd)
        markerName.assignString(cache, nerveName)
        dnodeIdentifier += 1

    # find element-xi for these data_coordinates
    dataNamesField = fm.findFieldByName(markerNameField)
    coordinates = findOrCreateFieldCoordinates(fm, coordinatesName)
    found_mesh_location = fm.createFieldFindMeshLocation(data_coordinates, coordinates, mesh)
    found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)
    xi_projected_data = {}
    nodeIter = datapoints.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        cache.setNode(node)
        element, xi = found_mesh_location.evaluateMeshLocation(cache, 3)
        marker_name = dataNamesField.evaluateString(cache)
        if element.isValid():
            addProjection = {marker_name:{"elementID": element.getIdentifier(), "xi": xi,"nodeID": node.getIdentifier()}}
            xi_projected_data.update(addProjection)
        node = nodeIter.next()

    result = datapoints.destroyAllNodes()

    return xi_projected_data

