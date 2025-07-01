"""
Generates a 3-D uterus mesh from a 1-D network layout, with variable
numbers of elements around, along and through wall.
"""
import copy

from cmlibs.maths.vectorops import add, cross, mult, set_magnitude, sub, normalize, magnitude, \
    axis_angle_to_rotation_matrix
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import setEftScaleFactorIds, remapEftNodeValueLabel
from scaffoldmaker.utils.geometry import sampleEllipsePoints, createEllipsoidPoints
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurvesSmooth, \
    smoothCurveSideCrossDerivatives, smoothCubicHermiteDerivativesLine, \
    getNearestLocationBetweenCurves, interpolateCubicHermite, sampleCubicHermiteCurves
from scaffoldmaker.utils.networkmesh import NetworkMesh, pathValueLabels
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData, \
    PatchTubeNetworkMeshSegment
from scaffoldmaker.utils.zinc_utils import group_add_connected_elements, get_nodeset_path_ordered_field_parameters
from scaffoldmaker.utils.tracksurface import TrackSurface

import math

class UterusTubeNetworkMeshGenerateData(TubeNetworkMeshGenerateData):

    def __init__(self, region, meshDimension, coordinateFieldName="coordinates",
                 startNodeIdentifier=1, startElementIdentifier=1, isLinearThroughShell=False, isShowTrimSurfaces=False):
        """
        :param isLinearThroughWall: Callers should only set if 3-D with no core.
        :param isShowTrimSurfaces: Tells junction generateMesh to make 2-D trim surfaces.
        """
        super(UterusTubeNetworkMeshGenerateData, self).__init__(
            region, meshDimension, coordinateFieldName, startNodeIdentifier, startElementIdentifier,
            isLinearThroughShell, isShowTrimSurfaces)
        self._fundusGroup = self.getOrCreateAnnotationGroup(get_uterus_term("fundus of uterus"))
        self._bodyGroup = self.getOrCreateAnnotationGroup(get_uterus_term("body of uterus"))
        self._upperCervixGroup = self.getOrCreateAnnotationGroup(("upper cervix", ""))
        self._lowerCervixGroup = self.getOrCreateAnnotationGroup(("lower cervix", ""))
        # force these annotation group names in base class
        self._leftGroup = self.getOrCreateAnnotationGroup(("left uterus", ""))
        self._rightGroup = self.getOrCreateAnnotationGroup(("right uterus", ""))
        self._dorsalGroup = self.getOrCreateAnnotationGroup(("dorsal uterus", ""))
        self._ventralGroup = self.getOrCreateAnnotationGroup(("ventral uterus", ""))

    def getFundusMeshGroup(self):
        return self._fundusGroup.getMeshGroup(self._mesh)

    def getBodyMeshGroup(self):
        return self._bodyGroup.getMeshGroup(self._mesh)

    def getUpperCervixMeshGroup(self):
        return self._upperCervixGroup.getMeshGroup(self._mesh)

    def getLowerCervixMeshGroup(self):
        return self._lowerCervixGroup.getMeshGroup(self._mesh)

class UterusTubeNetworkMeshBuilder(TubeNetworkMeshBuilder):
    """
    Adds left, right, dorsal, ventral, fundus, body annotations.
    Future: derive from BodyTubeNetworkMeshBuilder to get left/right/dorsal/ventral.
    """

    def createSegment(self, networkSegment):
        if networkSegment.isPatch():
            pathParametersList = [get_nodeset_path_ordered_field_parameters(
                self._layoutNodes, self._layoutCoordinates, pathValueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())]
            if self._layoutInnerCoordinates:
                pathParametersList.append(get_nodeset_path_ordered_field_parameters(
                    self._layoutNodes, self._layoutInnerCoordinates, pathValueLabels,
                    networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions()))
            elementsCountAround = self._defaultElementsCountAround
            elementsCountCoreBoxMinor = self._defaultElementsCountCoreBoxMinor
            coreBoundaryScalingMode = self._defaultCoreBoundaryScalingMode

            return PatchTubeNetworkMeshSegment(networkSegment, pathParametersList, elementsCountAround,
                                      self._elementsCountThroughShell, self._isCore, elementsCountCoreBoxMinor,
                                      self._elementsCountTransition, coreBoundaryScalingMode)  # args

        return super(UterusTubeNetworkMeshBuilder, self).createSegment(networkSegment)

    def generateMesh(self, generateData):
        super(UterusTubeNetworkMeshBuilder, self).generateMesh(generateData)
        # build temporary left/right dorsal/ventral groups
        mesh = generateData.getMesh()
        fundusMeshGroup = generateData.getFundusMeshGroup()
        bodyMeshGroup = generateData.getBodyMeshGroup()
        upperCervixMeshGroup = generateData.getUpperCervixMeshGroup()
        lowerCervixMeshGroup = generateData.getLowerCervixMeshGroup()
        leftMeshGroup = generateData.getLeftMeshGroup()
        rightMeshGroup = generateData.getRightMeshGroup()
        dorsalMeshGroup = generateData.getDorsalMeshGroup()
        ventralMeshGroup = generateData.getVentralMeshGroup()
        for networkSegment in self._networkMesh.getNetworkSegments():
            segment = self._segments[networkSegment]
            annotationTerms = segment.getAnnotationTerms()
            segment.addSideD3ElementsToMeshGroup(True, ventralMeshGroup)
            segment.addSideD3ElementsToMeshGroup(False, dorsalMeshGroup)
            for annotationTerm in annotationTerms:
                if "left fundus" in annotationTerm[0] or "right fundus" in annotationTerm[0]:
                    elementsCountRim = segment.getElementsCountRim()
                    elementsCountAlong = segment.getSampledElementsCountAlong()
                    elementsCountAround = segment.getElementsCountAround()
                    if "left" in annotationTerm[0]:
                        e1FundusLimitStart = 0
                        e1FundusRange = (elementsCountAround // 2 - (1 if elementsCountAround % 2 == 1 else 0)) // 2
                    else:
                        e1FundusLimitStart = elementsCountAround // 4
                        e1FundusRange = elementsCountAround // 2
                    for e1 in range(elementsCountAround):
                        for e2 in range(elementsCountAlong):
                            for e3 in range(elementsCountRim):
                                elementIdentifier = segment.getRimElementId(e1, e2, e3)
                                if elementIdentifier is not None:
                                    element = mesh.findElementByIdentifier(elementIdentifier)
                                    if (e1 >= e1FundusLimitStart and e1 < e1FundusLimitStart + e1FundusRange) or \
                                            ("left" in annotationTerm[0] and e1 >= elementsCountAround - e1FundusRange):
                                        fundusMeshGroup.addElement(element)
                                    else:
                                        bodyMeshGroup.addElement(element)
                if "body" in annotationTerm[0]:
                    segment.addAllElementsToMeshGroup(bodyMeshGroup)
                if "cervix" in annotationTerm[0]:
                    elementsCountRim = segment.getElementsCountRim()
                    elementsCountAlong = segment.getSampledElementsCountAlong()
                    elementsCountAround = segment.getElementsCountAround()
                    for e1 in range(elementsCountAround):
                        for e2 in range(elementsCountAlong):
                            for e3 in range(elementsCountRim):
                                elementIdentifier = segment.getRimElementId(e1, e2, e3)
                                if elementIdentifier is not None:
                                    element = mesh.findElementByIdentifier(elementIdentifier)
                                    if e2 >= int(elementsCountAlong * 0.5):
                                        lowerCervixMeshGroup.addElement(element)
                                    else:
                                        upperCervixMeshGroup.addElement(element)
                if "left" in annotationTerm[0]:
                    segment.addAllElementsToMeshGroup(leftMeshGroup)
                    break
                if "right" in annotationTerm[0]:
                    segment.addAllElementsToMeshGroup(rightMeshGroup)
                    break
            else:
                # segment on main axis
                segment.addSideD2ElementsToMeshGroup(True, leftMeshGroup)
                segment.addSideD2ElementsToMeshGroup(False, rightMeshGroup)


class MeshType_1d_uterus_network_layout1(MeshType_1d_network_layout1):
    """
    Defines uterus network layout.
    """

    @classmethod
    def getName(cls):
        return "1D Uterus Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default",
                "Human 1",
                "Human Pregnant 1",
                "Mouse 1",
                "Rat 1"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = parameterSetName
        if parameterSetName in ("Mouse 1", "Rat 1"):
            options["Structure"] = (
                "1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-36.1,"
                "17-18-19-20-21-22-23-24-25-26-27-28-29-30-31-32-36.2,"
                "#-33-34-35-36.3,"
                "36.4-37-38-39,"
                "39-40,"
                "40-41-42")
            options["Oviduct diameter"] = 0.5
            options["Oviduct length"] = 4.0
            options["Body length"] = 0.75
            options["Fundus width between oviducts"] = 0.8
            options["Fundus depth between oviducts"] = 0.5
            options["Cervical length"] = 0.25
            options["Cervical width around internal os"] = 0.8
            options["Cervical depth around internal os"] = 0.5
            options["Cervical width around external os"] = 0.8
            options["Cervical depth around external os"] = 0.5
            options["Vagina length"] = 0.5
            options["Vagina width around vagina orifice"] = 0.8
            options["Vagina depth around vagina orifice"] = 0.5
            options["Inner proportion body"] = 0.5
            options["Inner proportion cervix"] = 0.5
            options["Inner proportion vagina"] = 0.5
            options["Angle of anteversion degrees"] = 0.0
        elif parameterSetName == "Human Pregnant 1":
            options["Structure"] = (
                "1-2-3-4-5-6-7-8-9-31.1,"
                "10-11-12-13-14-15-16-17-18-31.2,"                
                "#-19-20-21-22-23-24-25-26-27-28-29-30-31.3,"
                "31.4-32-33-34-35-36-37-38-39-40-41-42-43,"
                "43-44,"
                "44-45-46-47-48")
            options["Oviduct diameter"] = 0.35
            options["Oviduct length"] = 10.0
            options["Body length"] = 30.0
            options["Fundus width between oviducts"] = 28.0
            options["Fundus depth between oviducts"] = 28.0
            options["Cervical length"] = 1.0
            options["Cervical width around internal os"] = 5.5
            options["Cervical depth around internal os"] = 3.8
            options["Cervical width around external os"] = 5.0
            options["Cervical depth around external os"] = 3.0
            options["Vagina length"] = 10.0
            options["Vagina width around vagina orifice"] = 1.25
            options["Vagina depth around vagina orifice"] = 1.25
            options["Inner proportion body"] = 0.8
            options["Inner proportion cervix"] = 0.15
            options["Inner proportion vagina"] = 0.8
            options["Angle of anteversion degrees"] = 70.0
        else:
            options["Structure"] = (
                "1-2-3-4-5-6-7-8-23.1,"
                "9-10-11-12-13-14-15-16-23.2,"
                "#-17-18-19-20-21-22-23.3,"
                "23.4-24-25-26-27-28-29,"
                "29-30-31,"
                "31-32-33-34-35-36-37-38")
            options["Oviduct diameter"] = 0.35
            options["Oviduct length"] = 10.0
            options["Body length"] = 7.0
            options["Fundus width between oviducts"] = 8.0
            options["Fundus depth between oviducts"] = 6.0
            options["Cervical length"] = 1.0
            options["Cervical width around internal os"] = 5.5
            options["Cervical depth around internal os"] = 3.8
            options["Cervical width around external os"] = 5.0
            options["Cervical depth around external os"] = 3.0
            options["Vagina length"] = 10.0
            options["Vagina width around vagina orifice"] = 1.25
            options["Vagina depth around vagina orifice"] = 1.25
            options["Inner proportion body"] = 0.75
            options["Inner proportion cervix"] = 0.15
            options["Inner proportion vagina"] = 0.8
            options["Angle of anteversion degrees"] = 70.0

        options["Define inner coordinates"] = True
        options["Inner proportion oviducts"] = 0.5

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Oviduct diameter",
            "Oviduct length",
            "Body length",
            "Fundus width between oviducts",
            "Fundus depth between oviducts",
            "Cervical length",
            "Cervical width around internal os",
            "Cervical depth around internal os",
            "Cervical width around external os",
            "Cervical depth around external os",
            "Vagina length",
            "Vagina width around vagina orifice",
            "Vagina depth around vagina orifice",
            "Inner proportion oviducts",
            "Inner proportion body",
            "Inner proportion cervix",
            "Inner proportion vagina",
            "Angle of anteversion degrees"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Oviduct diameter",
            "Oviduct length",
            "Body length",
            "Fundus width between oviducts",
            "Fundus depth between oviducts",
            "Cervical length",
            "Cervical width around internal os",
            "Cervical depth around internal os",
            "Cervical width around external os",
            "Cervical depth around external os",
            "Vagina length",
            "Vagina width around vagina orifice",
            "Vagina depth around vagina orifice"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
        for key in [
            "Inner proportion oviducts",
            "Inner proportion body",
            "Inner proportion cervix",
            "Inner proportion vagina",
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 0.9:
                options[key] = 0.9

        for key, angleRange in {
            "Angle of anteversion degrees": (-100.0, 100.0)
        }.items():
            if options[key] < angleRange[0]:
                options[key] = angleRange[0]
            elif options[key] > angleRange[1]:
                options[key] = angleRange[1]
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, NetworkMesh
        """
        parameterSetName = options['Base parameter set']
        structure = options["Structure"]
        oviductLength = options["Oviduct length"]
        oviductRadius = 0.5 * options["Oviduct diameter"]
        bodyLength = options["Body length"]
        halfFundusWidth = 0.5 * options["Fundus width between oviducts"]
        halfFundusDepth = 0.5 * options["Fundus depth between oviducts"]
        cervicalLength = options["Cervical length"]
        halfCervicalWidthInternalOs = 0.5 * options["Cervical width around internal os"]
        halfCervicalDepthInternalOs = 0.5 * options["Cervical depth around internal os"]
        halfCervicalWidthExternalOs = 0.5 * options["Cervical width around external os"]
        halfCervicalDepthExternalOs = 0.5 * options["Cervical depth around external os"]
        vaginaLength = options["Vagina length"]
        halfVaginaOrificeWidth = 0.5 * options["Vagina width around vagina orifice"]
        halfVaginaOrificeDepth = 0.5 * options["Vagina depth around vagina orifice"]
        anteversionAngleRad = math.radians(options["Angle of anteversion degrees"])
        innerProportionOviducts = options["Inner proportion oviducts"]
        innerProportionBody = options["Inner proportion body"]
        innerProportionCervix = options["Inner proportion cervix"]
        innerProportionVagina = options["Inner proportion vagina"]

        isPregnant = parameterSetName in 'Human Pregnant 1'
        isRat = parameterSetName in 'Rat 1'
        isRodent = parameterSetName in ("Mouse 1", "Rat 1")

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)
        # nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # set up element annotations
        uterusGroup = AnnotationGroup(region, get_uterus_term("uterus"))
        if isRodent:
            leftOviductGroup = AnnotationGroup(region, get_uterus_term("left uterine horn"))
            rightOviductGroup = AnnotationGroup(region, get_uterus_term("right uterine horn"))
            leftFundusGroup = AnnotationGroup(region, ("left fundus", ""))
            rightFundusGroup = AnnotationGroup(region, ("right fundus", ""))
        else:
            leftOviductGroup = AnnotationGroup(region, get_uterus_term("left oviduct"))
            rightOviductGroup = AnnotationGroup(region, get_uterus_term("right oviduct"))
            leftFundusGroup = AnnotationGroup(region, ("left fundus", ""))
            rightFundusGroup = AnnotationGroup(region, ("right fundus", ""))
        bodyGroup = AnnotationGroup(region, get_uterus_term("body of uterus"))
        cervixGroup = AnnotationGroup(region, ("cervix", ""))
        vaginaGroup = AnnotationGroup(region, get_uterus_term("vagina"))
        fundusPatchGroup = AnnotationGroup(region, ("fundus patch", ""))
        annotationGroups = [uterusGroup, leftOviductGroup, rightOviductGroup, leftFundusGroup,
                            rightFundusGroup, fundusPatchGroup,
                            bodyGroup, cervixGroup, vaginaGroup]

        uterusMeshGroup = uterusGroup.getMeshGroup(mesh)
        elementIdentifier = 1

        left = 0
        right = 1

        if isRodent:
            oviductElementsCount = 16
            fundusPatchElementsCount = 3
            fundusPostBodyJunctionElementsCount = 3
            cervixElementsCount = 1
            vaginaElementsCount = 2
        elif isPregnant:
            oviductElementsCount = 9
            fundusPatchElementsCount = 12
            fundusPostBodyJunctionElementsCount = 12
            cervixElementsCount = 1
            vaginaElementsCount = 4
        else:
            oviductElementsCount = 8
            fundusPatchElementsCount = 6
            fundusPostBodyJunctionElementsCount = 6
            cervixElementsCount = 2
            vaginaElementsCount = 7

        for side in (left, right):
            sideOviductGroup = leftOviductGroup if (side == left) else rightOviductGroup
            meshGroups = [uterusMeshGroup, sideOviductGroup.getMeshGroup(mesh)]
            for e in range(oviductElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

        meshGroups = [uterusMeshGroup, fundusPatchGroup.getMeshGroup(mesh)]
        for e in range(fundusPatchElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        meshGroups = [uterusMeshGroup, bodyGroup.getMeshGroup(mesh)]
        for e in range(fundusPostBodyJunctionElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        meshGroups = [cervixGroup.getMeshGroup(mesh)]
        for e in range(cervixElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        meshGroups = [vaginaGroup.getMeshGroup(mesh)]
        for e in range(vaginaElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        # need to ensure inner coordinates are at least defined:
        cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=0.75)
        innerCoordinates = findOrCreateFieldCoordinates(fieldmodule, "inner coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = 1

        # Calculate scales
        fundusScalePostBodyJunction = bodyLength / fundusPostBodyJunctionElementsCount
        cervicalScale = cervicalLength / cervixElementsCount
        vaginalScale = vaginaLength / vaginaElementsCount

        zero = [0.0, 0.0, 0.0]
        xBodyJunction = [0.0, 0.0, 0.0]

        if isPregnant:
            aThetaCervicalEnd = math.asin(halfCervicalWidthInternalOs / halfFundusWidth)
            aEllipse = bodyLength / math.cos(aThetaCervicalEnd)
        cThetaCervicalEnd = math.asin(halfCervicalDepthInternalOs / halfFundusDepth)
        cEllipse = bodyLength / math.cos(cThetaCervicalEnd)

        # Oviducts
        d1BodyJunction = []
        d2BodyJunction = []
        d3BodyJunction = []
        d12BodyJunction = []
        d13BodyJunction = []
        id2BodyJunction = []
        id3BodyJunction = []
        id12BodyJunction = []
        id13BodyJunction = []

        if isRodent:
            rC = bodyLength
            thetaLimit = math.radians(55.0)
            for side in (left, right):
                rTheta = rC * thetaLimit
                straightLength = oviductLength - rTheta
                theta = thetaLimit * (-1.0 if side == left else 1.0)
                xCurveEnd = add([rC * math.cos(theta), rC * math.sin(theta), 0.0], [-rC, 0.0, 0.0])
                d1 = [(-1.0 if side == left else 1.0) * rC * math.sin(theta),
                      (1.0 if side == left else -1.0) * rC * math.cos(theta),
                      0.0]
                xStart = sub(xCurveEnd, set_magnitude(d1, straightLength))
                xEnd = xBodyJunction
                nx = [xStart, xCurveEnd, xEnd]
                nd1 = [set_magnitude(d1, straightLength), set_magnitude(d1, rTheta),
                       [0.0, (1.0 if side == left else -1.0) * rTheta, 0.0]]
                xOviduct, d1Oviduct = sampleCubicHermiteCurves(nx, nd1, oviductElementsCount)[0:2]

                for i in range(oviductElementsCount + 1):
                    x = xOviduct[i]
                    d1 = d1Oviduct[i]
                    d3 = [0.0, 0.0, oviductRadius]
                    d2Direction = cross(normalize(d3), normalize([-c for c in d1]))
                    d2 = set_magnitude(d2Direction, -oviductRadius)
                    d12 = zero
                    d13 = zero
                    id2 = mult(d2, innerProportionOviducts)
                    id3 = mult(d3, innerProportionOviducts)

                    if i < oviductElementsCount:
                        node = nodes.findNodeByIdentifier(nodeIdentifier)
                        fieldcache.setNode(node)
                        setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3)
                        setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3)
                        nodeIdentifier += 1
                    else:
                        d1BodyJunction.append(d1)
                        d2BodyJunction.append(d2)
                        d3BodyJunction.append(d3)
                        d12BodyJunction.append(d12)
                        d13BodyJunction.append(d13)
                        id2BodyJunction.append(id2)
                        id3BodyJunction.append(id3)
                        id12BodyJunction.append(d12)
                        id13BodyJunction.append(d13)

        else:
            segment1Length = halfFundusWidth + oviductLength
            segment1LengthScale = segment1Length / oviductElementsCount
            elementsAlongHalfFundusWidth = math.ceil(halfFundusWidth / segment1LengthScale)

            for side in (left, right):
                xStart = [0.0, segment1Length * (-1.0 if side == left else 1.0), 0.0]
                d1Oviduct = [0.0, segment1LengthScale * (1.0 if side == left else -1.0), 0.0]
                for i in range(oviductElementsCount + 1):
                    x = add(xStart, mult(d1Oviduct, i))
                    if segment1LengthScale * i < oviductLength:
                        d2 = [-oviductRadius if side == left else oviductRadius, 0.0, 0.0]
                        d3 = [0.0, 0.0, oviductRadius]
                        d12 = [0.0, 0.0, 0.0]
                        d13 = [0.0, 0.0, 0.0]
                        # id2 = mult(d2, innerProportionOviducts)
                        # id3 = mult(d3, innerProportionOviducts)
                        # id12 = mult(d12, innerProportionOviducts)
                        # id13 = mult(d13, innerProportionOviducts)
                    else: # in ellipse zone
                        theta = math.acos(x[1] / halfFundusWidth)
                        if abs(halfFundusDepth * math.sin(theta)) < oviductRadius:
                            d2 = [-oviductRadius if side == left else oviductRadius, 0.0, 0.0]
                            d3 = [0.0, 0.0, oviductRadius]
                            d12 = [0.0, 0.0, 0.0]
                            d13 = [0.0, 0.0, 0.0]
                        else:
                            d2 = [halfFundusDepth * math.sin(theta) * (-1.0 if side == left else 1.0), 0.0, 0.0] # fundusHeight??
                            d3 = [0.0, 0.0, halfFundusDepth * math.sin(theta)]
                            d12 = [halfFundusDepth * math.cos(theta) * (0.5 * math.pi / elementsAlongHalfFundusWidth),
                                   0.0, 0.0]
                            d13 = [0.0, 0.0,
                                   halfFundusDepth * math.cos(theta) * (0.5 * math.pi / elementsAlongHalfFundusWidth) *
                                   (-1.0 if side == left else 1.0)]

                        # if abs(halfFundusDepth * innerProportionBody * math.sin(theta)) < uterineTubeRadius * innerProportionUterineTubes:
                        #     id2 = [innerProportionUterineTubes * (-uterineTubeRadius if side == left else uterineTubeRadius), 0.0, 0.0]
                        #     id3 = [0.0, 0.0, uterineTubeRadius * innerProportionUterineTubes]
                        #     id12 = [0.0, 0.0, 0.0]
                        #     id13 = [0.0, 0.0, 0.0]
                        # else:
                        #     id2 = [innerProportionBody * halfFundusDepth * math.sin(theta) * (-1.0 if side == left else 1.0), 0.0,
                        #               0.0]  # fundusHeight??
                        #     id3 = [0.0, 0.0, innerProportionBody * halfFundusDepth * math.sin(theta)]
                        #     id12 = [innerProportionBody * halfFundusDepth * math.cos(theta) * (0.5 * math.pi / elementsAlongHalfFundusWidth),
                        #            0.0, 0.0]
                        #     id13 = [0.0, 0.0,
                        #            innerProportionBody * halfFundusDepth * math.cos(theta) * (0.5 * math.pi / elementsAlongHalfFundusWidth) *
                        #            (-1.0 if side == left else 1.0)]

                    id2 = mult(d2, innerProportionOviducts)
                    id3 = mult(d3, innerProportionOviducts)
                    id12 = mult(d12, innerProportionOviducts)
                    id13 = mult(d13, innerProportionOviducts)

                    if i < oviductElementsCount:
                        node = nodes.findNodeByIdentifier(nodeIdentifier)
                        fieldcache.setNode(node)
                        setNodeFieldParameters(coordinates, fieldcache, x, d1Oviduct, d2, d3, d12, d13)
                        setNodeFieldParameters(innerCoordinates, fieldcache, x, d1Oviduct, id2, id3, id12, id13)
                        nodeIdentifier += 1
                    else:
                        d1BodyJunction.append(d1Oviduct)
                        d2BodyJunction.append(d2)
                        d3BodyJunction.append(d3)
                        d12BodyJunction.append(d12)
                        d13BodyJunction.append(d13)
                        id2BodyJunction.append(id2)
                        id3BodyJunction.append(id3)
                        id12BodyJunction.append(id12)
                        id13BodyJunction.append(id13)

        xFundusPatchStart = [-bodyLength, 0.0, 0.0]
        d1FundusPatch = [bodyLength / fundusPatchElementsCount, 0.0, 0.0]
        nxPatch = []
        nd1Patch = []
        nd2Patch = []
        nd3Patch = []

        for i in range(fundusPatchElementsCount + 1):
            x = [xFundusPatchStart[0] + d1FundusPatch[0] * i, 0.0, 0.0]
            nxPatch.append(x)
            nd1Patch.append(d1FundusPatch)
            xi = i / fundusPatchElementsCount
            width = xi * halfFundusWidth + (1.0 - xi) * (halfFundusWidth * 0.01 if isRodent else
                                                         halfCervicalWidthInternalOs)
            if isPregnant:
                thetaA = math.acos(x[0] / aEllipse)
                width = halfFundusWidth * math.sin(thetaA)
            thetaC = math.acos(x[0] / cEllipse)
            depth = halfFundusDepth * math.sin(thetaC)
            nd2Patch.append([0.0, width, 0.0])
            nd3Patch.append([0.0, 0.0, depth])
        nd13Patch = smoothCurveSideCrossDerivatives(nxPatch, nd1Patch, [nd3Patch])[0]
        nd12Patch = smoothCurveSideCrossDerivatives(nxPatch, nd1Patch, [nd2Patch])[0]

        for i in range(fundusPatchElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            id3 = mult(nd3Patch[i], innerProportionBody)
            id13 = mult(nd13Patch[i], innerProportionBody)
            if isPregnant:
                id2 = mult(nd2Patch[i], innerProportionBody)
                id12 = mult(nd12Patch[i], innerProportionBody)
            else:
                xi = i / fundusPatchElementsCount
                width = xi * halfFundusWidth * innerProportionBody + \
                        (1.0 - xi) * (halfFundusWidth * 0.01 * innerProportionBody if isRodent else
                                      halfCervicalWidthInternalOs * innerProportionCervix)

                id2 = [0.0, width, 0.0]
                id12 = [0.0,
                        (halfFundusWidth * innerProportionBody -
                         halfCervicalWidthInternalOs * innerProportionCervix) / fundusPatchElementsCount,
                        0.0]
            setNodeFieldParameters(coordinates, fieldcache, nxPatch[i], nd1Patch[i], nd2Patch[i], nd3Patch[i],
                                   nd12Patch[i], nd13Patch[i])
            setNodeFieldParameters(innerCoordinates, fieldcache, nxPatch[i], nd1Patch[i], id2, id3, id12, id13)
            nodeIdentifier += 1

        # body junction
        node = nodes.findNodeByIdentifier(nodeIdentifier)
        fieldcache.setNode(node)
        for field in (coordinates, innerCoordinates):
            field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, xBodyJunction)
        for side in (left, right):
            version = 1 if side == left else 2
            d1 = d1BodyJunction[version - 1]
            d2 = d2BodyJunction[version - 1]
            d3 = d3BodyJunction[version - 1]
            d12 = d12BodyJunction[version - 1]
            d13 = d13BodyJunction[version - 1]
            id2 = id2BodyJunction[version - 1]
            id3 = id3BodyJunction[version - 1]
            id12 = id12BodyJunction[version - 1]
            id13 = id13BodyJunction[version - 1]
            setNodeFieldVersionDerivatives(coordinates, fieldcache, version, d1, d2, d3, d12, d13)
            setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, d1, id2, id3, id12, id13)

        version = 3
        id3 = mult(nd3Patch[-1], innerProportionBody)
        id13 = mult(nd13Patch[-1], innerProportionBody)
        if isPregnant:
            id2 = mult(nd2Patch[-1], innerProportionBody)
            id12 = mult(nd12Patch[-1], innerProportionBody)
        else:
            id2 = [0.0, halfFundusWidth * innerProportionBody, 0.0]
            id12 = [0.0,
                    (halfCervicalWidthInternalOs * innerProportionCervix -
                     halfFundusWidth * innerProportionBody) / fundusPatchElementsCount,
                    0.0]
        setNodeFieldVersionDerivatives(coordinates, fieldcache, version, nd1Patch[-1], nd2Patch[-1], nd3Patch[-1],
                                       nd12Patch[-1], nd13Patch[-1])
        setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, nd1Patch[-1], id2, id3, id12, id13)

        # Post body junction
        nd12 = []
        nd13 = []
        nxBody = []
        nd1Body = []
        nd2Body = []
        nd3Body = []

        dWidth = (halfCervicalWidthInternalOs - halfFundusWidth) / fundusPostBodyJunctionElementsCount
        dDepth = (halfCervicalDepthInternalOs - halfFundusDepth) / fundusPostBodyJunctionElementsCount

        for i in range(fundusPostBodyJunctionElementsCount + 1):
            x = [fundusScalePostBodyJunction * i, 0.0, 0.0]
            d1 = [fundusScalePostBodyJunction, 0.0, 0.0]
            xi = i / fundusPostBodyJunctionElementsCount
            width = xi * halfCervicalWidthInternalOs + (1.0 - xi) * halfFundusWidth
            if isPregnant:
                thetaA = math.acos(x[0] / aEllipse)
                width = halfFundusWidth * math.sin(thetaA)
            thetaC = math.acos(x[0] / cEllipse)
            depth = halfFundusDepth * math.sin(thetaC)
            d2 = [0.0, width, 0.0]
            d3 = [0.0, 0.0, depth]
            nxBody.append(x)
            nd1Body.append(d1)
            nd2Body.append(d2)
            nd3Body.append(d3)
            nd12.append([0.0, dWidth, 0.0])
            nd13.append([0.0, 0.0, dDepth])

        if dWidth == 0.0:
            nd12Body = nd12
        else:
            nd12Body = smoothCubicHermiteDerivativesLine(nd2Body, nd12)
        if dDepth == 0.0:
            nd13Body = nd13
        else:
            nd13Body = smoothCubicHermiteDerivativesLine(nd3Body, nd13)

        for i in range(fundusPostBodyJunctionElementsCount + 1):
            x = nxBody[i]
            d1 = nd1Body[i]
            d2 = nd2Body[i]
            d3 = nd3Body[i]
            d12 = nd12Body[i]
            d13 = nd13Body[i]
            id3 = mult(d3, innerProportionBody)
            id13 = mult(d13, innerProportionBody)
            if isPregnant:
                id2 = mult(d2, innerProportionBody)
                id12 = mult(d12, innerProportionBody)
            else:
                xi = i / fundusPostBodyJunctionElementsCount
                width = xi * halfCervicalWidthInternalOs * innerProportionCervix + \
                        (1.0 - xi) * halfFundusWidth * innerProportionBody
                id2 = [0.0, width, 0.0]
                id12 = [0.0,
                        (halfCervicalWidthInternalOs * innerProportionCervix -
                         halfFundusWidth * innerProportionBody) / fundusPatchElementsCount,
                        0.0]

            if i:
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            else:
                version = 4
                setNodeFieldVersionDerivatives(coordinates, fieldcache, version, d1, d2, d3, d12, d13)
                setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, d1, id2, id3, id12, id13)
            nodeIdentifier += 1

        # Pre-calculate points in cervix and vagina
        rotMat = axis_angle_to_rotation_matrix([0.0, 1.0, 0.0], anteversionAngleRad)
        vaginalStartX = bodyLength + cervicalLength

        turningPtIdx = 1
        xTurningPt = [vaginalStartX + vaginalScale * turningPtIdx, 0.0, 0.0]
        xTranslateMat = sub(xTurningPt, xBodyJunction)

        li = list(range(2, vaginaElementsCount + 1))

        # one node before vagina
        nx = []
        nd1 = []
        nx.append([vaginalStartX - cervicalScale, 0.0, 0.0])
        nd1.append([cervicalScale, 0.0, 0.0])
        for i in li:
            x = [vaginalStartX + vaginalScale * i, 0.0, 0.0]
            d1 = [0.5 * (vaginalScale + cervicalScale) if i == 0 else vaginalScale, 0.0, 0.0]
            if i > turningPtIdx:
                xTranslate = sub(x, xTranslateMat)
                xRot = [rotMat[j][0] * xTranslate[0] + rotMat[j][1] * xTranslate[1] + rotMat[j][2] * xTranslate[2]
                        for j in range(3)]
                x = add(xRot, xTranslateMat)
                d1 = [rotMat[j][0] * d1[0] + rotMat[j][1] * d1[1] + rotMat[j][2] * d1[2] for j in range(3)]
            nx.append(x)
            nd1.append(d1)

        xSampledCurve, d1SampledCurve = sampleCubicHermiteCurves(nx, nd1, vaginaElementsCount + 1,
                                                                 arcLengthDerivatives=True)[0:2]
        d1SmoothedCurve = smoothCubicHermiteDerivativesLine(xSampledCurve, d1SampledCurve, fixStartDirection=True)

        # assign points to cervix
        nodeIdentifier -= 1
        cervixStartX = bodyLength
        nxCervix = [[cervixStartX, 0.0, 0.0]] + xSampledCurve[0:2]
        nd1Cervix = [[0.5 * (fundusScalePostBodyJunction + cervicalScale), 0.0, 0.0]] + d1SmoothedCurve[0:2]
        nxCervixSampled, nd1CervixSampled = \
            sampleCubicHermiteCurves(nxCervix, nd1Cervix, cervixElementsCount, arcLengthDerivatives=True)[0:2]
        nd1CervixSampled[0] = [0.5 * (fundusScalePostBodyJunction + cervicalScale), 0.0, 0.0]
        nd1CervixSmoothed = smoothCubicHermiteDerivativesLine(nxCervixSampled, nd1CervixSampled,
                                                              fixStartDerivative=True,
                                                              fixEndDerivative=True)
        nd2Cervix = []
        nd3Cervix = []
        for i in range(cervixElementsCount + 1):
            xi = i / cervixElementsCount
            d1 = nd1CervixSmoothed[i]
            halfCervixWidth = xi * halfCervicalWidthExternalOs + (1 - xi) * halfCervicalWidthInternalOs
            halfCervixDepth = xi * halfCervicalDepthExternalOs + (1 - xi) * halfCervicalDepthInternalOs
            d2 = [0.0, halfCervixWidth, 0.0]
            nd2Cervix.append(d2)
            d3Direction = cross(normalize(d1), normalize(d2))
            nd3Cervix.append(set_magnitude(d3Direction, halfCervixDepth))

        nxVagina = []
        nd1Vagina = []
        nd2Vagina = []
        nd3Vagina = []

        del xSampledCurve[0]
        del d1SmoothedCurve[0]

        for i in range(vaginaElementsCount + 1):
            xi = i / vaginaElementsCount
            nxVagina.append(xSampledCurve[i])
            nd1Vagina.append(d1SmoothedCurve[i])
            d2 = \
                interpolateCubicHermite(
                    [0.0, halfCervicalWidthExternalOs, 0.0],
                    [0.0, (halfCervicalWidthExternalOs - halfCervicalWidthInternalOs) / cervixElementsCount, 0.0],
                    [0.0, halfVaginaOrificeWidth, 0.0],
                    [0.0, (halfVaginaOrificeWidth - halfCervicalWidthExternalOs) / vaginaElementsCount, 0.0], xi)

            if i > turningPtIdx:
                d2 = [rotMat[j][0] * d2[0] + rotMat[j][1] * d2[1] + rotMat[j][2] * d2[2] for j in range(3)]
            nd2Vagina.append(d2)

            d3 = cross(normalize(d1SmoothedCurve[i]), normalize(d2))
            d3Interpolated = \
                interpolateCubicHermite(
                    [0.0, 0.0, halfCervicalDepthExternalOs],
                    [0.0, 0.0, (halfCervicalDepthExternalOs - halfCervicalDepthInternalOs) / cervixElementsCount],
                    [0.0, 0.0, halfVaginaOrificeDepth],
                    [0.0, 0.0, (halfVaginaOrificeDepth - halfCervicalDepthExternalOs) / vaginaElementsCount], xi)
            d3 = set_magnitude(d3, magnitude(d3Interpolated))
            nd3Vagina.append(d3)

        nd2 = nd2Cervix + nd2Vagina[1:]
        nd3 = nd3Cervix + nd3Vagina[1:]
        nd12 = [sub(nd2[c + 1], nd2[c]) for c in range(len(nd2) - 1)]
        nd12.append(nd12[-1])
        nd13 = [sub(nd3[c + 1], nd3[c]) for c in range(len(nd3) - 1)]
        nd13.append(nd13[-1])

        if halfCervicalWidthExternalOs - halfCervicalWidthInternalOs:
            nd12 = smoothCubicHermiteDerivativesLine(nd2, nd12)
        if halfCervicalDepthExternalOs - halfCervicalDepthInternalOs:
            nd13 = smoothCubicHermiteDerivativesLine(nd3, nd13)

        nd12Cervix = nd12[0:cervixElementsCount + 1]
        nd12Vagina = nd12[cervixElementsCount:]
        nd13Cervix = nd13[0:cervixElementsCount + 1]
        nd13Vagina = nd13[cervixElementsCount:]

        for i in range(cervixElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = nxCervixSampled[i]
            d1 = nd1CervixSmoothed[i]
            d2 = nd2Cervix[i]
            d3 = nd3Cervix[i]
            if i == 0:
                d12 = mult(add(nd12Body[-1], nd12Cervix[i]), 0.5)
                d13 = mult(add(nd13Body[-1], nd13Cervix[i]), 0.5)
            else:
                d12 = nd12Cervix[i]
                d13 = nd13Cervix[i]

            id2 = mult(d2, innerProportionCervix)
            id3 = mult(d3, innerProportionCervix)
            id12 = mult(d12, innerProportionCervix)
            id13 = mult(d13, innerProportionCervix)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1

        for i in range(vaginaElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = nxVagina[i]
            d1 = nd1Vagina[i]
            d2 = nd2Vagina[i]
            d3 = nd3Vagina[i]
            d12 = nd12Vagina[i]
            d13 = nd13Vagina[i]
            if i == 0:
                innerProportion = innerProportionCervix
            else:
                innerProportion = innerProportionVagina
            id2 = mult(d2, innerProportion)
            id3 = mult(d3, innerProportion)
            id12 = mult(d12, innerProportion)
            id13 = mult(d13, innerProportion)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1

        return annotationGroups, networkMesh


    @classmethod
    def getInteractiveFunctions(cls):
        """
        Edit base class list to include only valid functions.
        """
        interactiveFunctions = super(MeshType_1d_uterus_network_layout1, cls).getInteractiveFunctions()
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Edit structure...":
                interactiveFunctions.remove(interactiveFunction)
                break
        return interactiveFunctions


class MeshType_3d_uterus2(Scaffold_base):
    """
    Generates a 3-D uterus mesh from a 1-D network layout with variable numbers of elements around, along and through
    wall.
    Magnitude of D2 and D3 are the radii of the uterus in the respective directions.
    """

    @classmethod
    def getName(cls):
        return '3D Uterus 2'

    @classmethod
    def getParameterSetNames(cls):
        return [
            'Default',
            'Human 1',
            'Human Pregnant 1',
            'Mouse 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):

        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options = {
            'Base parameter set': useParameterSetName,
            'Network layout': ScaffoldPackage(MeshType_1d_uterus_network_layout1,
                                              defaultParameterSetName=useParameterSetName),
            'Number of elements around': 24,
            'Number of elements around oviduct': 8,
            'Number of elements through wall': 1,
            'Number of elements along oviduct': 6,
            'Number of elements along body': 4,
            'Number of elements along cervix': 2,
            'Number of elements along vagina': 6,
            'Target element density along longest segment': 5.5,
            'Use linear through wall': False, # True needs work
            'Show trim surfaces': False,
            'Show pseudo-segment': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Mouse' in parameterSetName or 'Rat' in parameterSetName:
            options['Number of elements around'] = 12
            options['Number of elements around oviduct'] = 8
            options['Number of elements along oviduct'] = 12
            options['Number of elements along body'] = 2
            options['Number of elements along cervix'] = 1
            options['Number of elements along vagina'] = 2

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            'Network layout',
            'Number of elements around',
            'Number of elements around oviduct',
            'Number of elements through wall',
            'Number of elements along oviduct',
            'Number of elements along body',
            'Number of elements along cervix',
            'Number of elements along vagina',
            'Use linear through wall',
            'Show trim surfaces',
            'Show pseudo-segment',
            'Refine',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [MeshType_1d_uterus_network_layout1]
        return []

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
        if optionName == 'Network layout':
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_uterus_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if (options["Network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Network layout")):
            options["Network layout"] = ScaffoldPackage(MeshType_1d_uterus_network_layout1)
        for key in [
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Number of elements around',
            'Number of elements around oviduct']:
            if options[key] < 4:
                options[key] = 4
            elif (options[key] % 4) > 0:
                options[key] += options[key] % 4
        if options["Number of elements through wall"] < 1:
            options["Number of elements through wall"] = 1

        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic hermite or bicubic hermite-linear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ("Default", "Human 1")
        isPregnant = parameterSetName in ("Human Pregnant 1")

        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        # elementsCountThroughWall = options["Number of elements through wall"]
        # elementsCountAroundUterineTube = options["Number of elements around oviduct"]
        # elementsCountAroundBody = options["Number of elements around"]
        # useLinerThroughWall = options["Use linear through wall"]

        annotationElementsCountsAlong = []
        annotationElementsCountsAround = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            elementsCountAlong = 0
            elementsCountAround = 0
            name = layoutAnnotationGroup.getName()
            if "oviduct" in name or "uterine horn" in name:
                elementsCountAlong = options['Number of elements along oviduct']
            elif "body" in name or "fundus patch" in name:
                elementsCountAlong = options['Number of elements along body']
            elif "cervix" in name:
                elementsCountAlong = options['Number of elements along cervix']
            elif "vagina" in name:
                elementsCountAlong = options['Number of elements along vagina']
            annotationElementsCountsAlong.append(elementsCountAlong)

            if "oviduct" in name or "uterine horn" in name or "left fundus" in name or "right fundus" in name:
                elementsCountAround = options['Number of elements around oviduct']
            elif "body" in name or "fundus patch" in name:
                elementsCountAround = options['Number of elements around']
            annotationElementsCountsAround.append(elementsCountAround)

        uterusTubeNetworkMeshBuilder = UterusTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=2.0, # not used
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAlong=annotationElementsCountsAlong,
            defaultElementsCountAround=options['Number of elements around'],
            annotationElementsCountsAround=annotationElementsCountsAround,
            elementsCountThroughShell=options["Number of elements through wall"],
            useOuterTrimSurfaces=False)
        uterusTubeNetworkMeshBuilder.build()

        generateData = UterusTubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughShell=options["Use linear through wall"],
            isShowTrimSurfaces=options["Show trim surfaces"])

        uterusTubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()
        nodeIdentifier, elementIdentifier = generateData.getNodeElementIdentifiers()

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # fieldcache = fieldmodule.createFieldcache()
        # nodetemplate = nodes.createNodetemplate()
        # nodetemplate.defineField(coordinates)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fieldmodule.findMeshByDimension(3)
        # useCrossDerivatives = False
        # if useLinerThroughWall:
        #     eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        # else:
        #     eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        #
        # eft = eftfactory.createEftBasic()
        #
        # elementtemplate = mesh.createElementtemplate()
        # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        # elementtemplate.defineField(coordinates, -1, eft)
        # elementtemplateX = mesh.createElementtemplate()
        # elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        #
        # mesh2d = fieldmodule.findMeshByDimension(2)
        #
        # # Identify boundary nodes
        # fieldmodule.defineAllFaces()
        # fundusPatchGroup = getAnnotationGroupForTerm(annotationGroups, ("fundus patch", "None"))
        # leftFundusGroup = getAnnotationGroupForTerm(annotationGroups, ("left fundus", "None"))
        # rightFundusGroup = getAnnotationGroupForTerm(annotationGroups, ("right fundus", "None"))
        # bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("body of uterus"))
        #
        # for annotationGroup in [fundusPatchGroup, leftFundusGroup, rightFundusGroup, bodyGroup]:
        #     annotationGroup.addSubelements()
        #
        # is_fundusPatch = fundusPatchGroup.getGroup()
        # is_leftFundus = leftFundusGroup.getGroup()
        # is_rightFundus = rightFundusGroup.getGroup()
        # is_body = bodyGroup.getGroup()
        #
        # is_leftSegmentBoundary = fieldmodule.createFieldAnd(is_leftFundus, is_fundusPatch)
        # is_rightSegmentBoundary = fieldmodule.createFieldAnd(is_rightFundus, is_fundusPatch)
        # is_bodySegmentBoundary = fieldmodule.createFieldAnd(is_body, is_fundusPatch)
        #
        # leftSegmentBoundary = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
        #                                                          ["left segment boundary", "None"])
        # leftSegmentBoundary.getMeshGroup(mesh2d).addElementsConditional(is_leftSegmentBoundary)
        # leftSegmentBoundaryNodeset = leftSegmentBoundary.getNodesetGroup(nodes)
        #
        # rightSegmentBoundary = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
        #                                                          ["right segment boundary", "None"])
        # rightSegmentBoundary.getMeshGroup(mesh2d).addElementsConditional(is_rightSegmentBoundary)
        # rightSegmentBoundaryNodeset = rightSegmentBoundary.getNodesetGroup(nodes)
        #
        # bodySegmentBoundary = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
        #                                                          ["body segment boundary", "None"])
        # bodySegmentBoundary.getMeshGroup(mesh2d).addElementsConditional(is_bodySegmentBoundary)
        # bodySegmentBoundaryNodeset = bodySegmentBoundary.getNodesetGroup(nodes)
        #
        # for nodeset in [leftSegmentBoundaryNodeset, rightSegmentBoundaryNodeset, bodySegmentBoundaryNodeset]:
        #     nodeIDs = []
        #     xBoundary = []
        #     d2Boundary = []
        #     nodeIter = nodeset.createNodeiterator()
        #     node = nodeIter.next()
        #     fieldcache.setNode(node)
        #     while node.isValid():
        #         nodeIDs.append(node.getIdentifier())
        #         xBoundary.append(coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)[1])
        #         d2Boundary.append(coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)[1])
        #         node = nodeIter.next()
        #         fieldcache.setNode(node)
        #
        #     if nodeset == leftSegmentBoundaryNodeset:
        #         nodeIdsTmp = copy.deepcopy(nodeIDs)
        #         xTmp = copy.deepcopy(xBoundary)
        #         d2Tmp = copy.deepcopy(d2Boundary)
        #         # reorder
        #         nodeIdsLeftSegmentBoundary = reorderNodesForMapping(nodeIdsTmp, elementsCountThroughWall)
        #         xLeftSegmentBoundary = reorderNodesForMapping(xTmp, elementsCountThroughWall)
        #         d2LeftSegmentBoundary = reorderNodesForMapping(d2Tmp, elementsCountThroughWall)
        #
        #     elif nodeset == rightSegmentBoundaryNodeset:
        #         nodeIdsRightSegmentBoundary = copy.deepcopy(nodeIDs)
        #         xRightSegmentBoundary = copy.deepcopy(xBoundary)
        #         d2RightSegmentBoundary = copy.deepcopy(d2Boundary)
        #
        #     elif nodeset == bodySegmentBoundaryNodeset:
        #         # print(nodeIDs)
        #         nodeIdsBodySegmentBoundary = copy.deepcopy(nodeIDs[4 * (elementsCountThroughWall + 1):]) # remove the 2 corners on both sides
        #         xBodySegmentBoundary = copy.deepcopy(xBoundary[4 * (elementsCountThroughWall + 1):])
        #         d2BodySegmentBoundary = copy.deepcopy(d2Boundary[4 * (elementsCountThroughWall + 1):])
        #
        # # print("left", nodeIdsLeftSegmentBoundary)
        # # print("right", nodeIdsRightSegmentBoundary)
        # # print("body", nodeIdsBodySegmentBoundary)
        #
        # # Remove fundus patch
        # if not options['Show pseudo-segment']:
        #     destroyGroup = fundusPatchGroup.getGroup()
        #     mesh.destroyElementsConditional(destroyGroup)
        #     nodes.destroyNodesConditional(destroyGroup)
        #     del destroyGroup
        #
        uterusGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("uterus"))
        # fundusGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("fundus of uterus"))
        #
        # add human specific annotations
        allMarkers = {}
        if isHuman:
            allMarkers = {"junction of left round ligament with uterus": {"x": [-0.495368, -4.54665, 0.0]},
                          "junction of right round ligament with uterus": {"x": [-0.495368, 4.54665, 0.0]},
                          "junction of left uterosacral ligament with uterus": {"x": [7.0, -2.75, 0.0]},
                          "junction of right uterosacral ligament with uterus": {"x": [7.0, 2.75, 0.0]}}
        elif isPregnant:
            allMarkers = {"junction of left round ligament with uterus": {"x": [-2.81269, -16.7248, 0.0]},
                          "junction of right round ligament with uterus": {"x": [-2.81269, 16.7248, 0.0]},
                          "junction of left uterosacral ligament with uterus": {"x": [30, -2.75, 0.0]},
                          "junction of right uterosacral ligament with uterus": {"x": [30, 2.75, 0]}}

        for key in allMarkers:
            x = allMarkers[key]["x"]
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                        get_uterus_term(key), isMarker=True)
            markerNode = group.createMarkerNode(nodeIdentifier, coordinates, x)
            nodeIdentifier = markerNode.getIdentifier() + 1
            for group in annotationGroups:
                group.getNodesetGroup(nodes).addNode(markerNode)

        # # Make patch
        # networkLayoutOptions = options['Network layout']
        # networkSettings = networkLayoutOptions.getScaffoldSettings()
        # fundusHeight = networkSettings["Fundus height scale"] * networkSettings["Body length"]
        # halfFundusDepth = networkSettings["Fundus depth between oviducts"] * 0.5
        # halfFundusWidth = networkSettings["Fundus width between oviducts"] * 0.5
        # innerProportionBody = networkSettings["Inner proportion body"]
        #
        # networkSegments = networkMesh.getNetworkSegments()
        # layoutFieldmodule = layoutRegion.getFieldmodule()
        # layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # layoutCoordinates = layoutFieldmodule.findFieldByName('coordinates')
        # bodyJunctionGroup = layoutFieldmodule.findFieldByName("body junction").castGroup()
        # bodyJunctionNodeset = bodyJunctionGroup.getNodesetGroup(layoutNodes)
        #
        # bodyJunctionNodeID = []
        # nodeiterator = bodyJunctionNodeset.createNodeiterator()
        # node = nodeiterator.next()
        # while node.isValid():
        #     bodyJunctionNodeID.append(node.getIdentifier())
        #     node = nodeiterator.next()
        #
        # nodes_IDs = []
        # versions = []
        # for segment in networkSegments:
        #     segmentNodesID = segment.getNodeIdentifiers()
        #     for i in range(len(segmentNodesID)):
        #         if segmentNodesID[i] in bodyJunctionNodeID:
        #             nodes_IDs.append(segmentNodesID[i])
        #             versions.append(segment.getNodeVersions()[i])
        #
        # cx, cd1, cd2, cd3, cd12, cd13 = get_nodeset_path_ordered_field_parameters(
        #     bodyJunctionNodeset, layoutCoordinates,
        #     [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
        #      Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
        #      Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3], nodes_IDs, versions)
        #
        # xBodyJunction = cx[0]
        # d2BodyJunctionSegment1 = set_magnitude(cd2[1], fundusHeight)
        # d3BodyJunctionSegment1 = set_magnitude(cd3[1], halfFundusDepth)
        # zero = [0.0, 0.0, 0.0]
        #
        # # dorsal = 0
        # # ventral = 1
        # # xAll = []
        # # d1All = []
        # # d2All = []
        # #
        # # for n3 in range(elementsCountThroughWall + 1):
        # #     sf = (1.0 - innerProportionBody) / elementsCountThroughWall * n3 + innerProportionBody
        # #     d2 = mult(d2BodyJunctionSegment1, sf)
        # #     d3 = mult(d3BodyJunctionSegment1, sf)
        # #     xMid, d1Mid = sampleEllipsePoints(xBodyJunction, d2, d3, 0.5 * math.pi, 1.5 * math.pi,
        # #                                       int(elementsCountAroundUterineTube * 0.5))
        # #
        # #     elementsCountAlongPatch = int((elementsCountAroundBody - elementsCountAroundUterineTube) * 0.5)
        # #     newNodesCountAlongPatch = elementsCountAlongPatch - 1
        # #
        # #     # print('nodesRightSegmentBoundary = ', nodeIdsRightSegmentBoundary)
        # #     # print('len(xLeftSegmentBoundary) = ', len(xLeftSegmentBoundary))
        # #     # for n in range(1, len(xMid) - 1):
        # #     n = int(0.5 * len(xMid))
        # #     idxAtSegment = n + int(((elementsCountAroundUterineTube * 0.5) + 1) * n3)
        # #     nStart = xRightSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else xLeftSegmentBoundary[
        # #         idxAtSegment]
        # #     nEnd = xLeftSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else xRightSegmentBoundary[
        # #         idxAtSegment]
        # #     d2Start = d2RightSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else d2LeftSegmentBoundary[
        # #         idxAtSegment]
        # #     d2End = [-1 * c for c in d2LeftSegmentBoundary[idxAtSegment]] \
        # #         if n <= int(0.5 * len(xMid)) else [-1 * c for c in d2RightSegmentBoundary[idxAtSegment]]
        # #
        # #     sf = 2.0 if isPregnant else 1.5 # make top round
        # #     nx = [nStart, xMid[n], nEnd]
        # #     nd1 = [d2Start,
        # #            [0.0, sf * (-halfFundusWidth if n <= int(0.5 * len(xMid)) else halfFundusWidth), 0.0],
        # #            d2End]
        # #
        # #     px, pd1 = sampleCubicHermiteCurvesSmooth(nx, nd1, elementsCountAlongPatch)[:2]
        # #     pd1 = smoothCubicHermiteDerivativesLine(px, pd1)
        # #
        # #     xAlongApex = px[1:-1]
        # #     # d1AlongApex = pd1[1:-1]
        # #     #
        # #     # xRowList = [xBodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        # #     #                                  newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]]
        # #     # d1RowList = [[zero for c in range(newNodesCountAlongPatch)]]
        # #     # d2RowList = [d2BodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        # #     #                                    newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]]
        # #     # frontBodySegmentBoundaryNodes = nodeIdsBodySegmentBoundary[
        # #     #                                 newNodesCountAlongPatch * 2 * n3:
        # #     #                                 newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]
        # #     xDorsalBodySegmentBoundary = xBodySegmentBoundary[
        # #                                 newNodesCountAlongPatch * 2 * n3:
        # #                                 newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]
        # #     d2DorsalBodySegmentBoundary = d2BodySegmentBoundary[
        # #                                  newNodesCountAlongPatch * 2 * n3:
        # #                                  newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]
        # #
        # #     xVentralBodySegmentBoundary = \
        # #        xBodySegmentBoundary[
        # #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        # #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2]
        # #     d2VentralBodySegmentBoundary = \
        # #         d2BodySegmentBoundary[
        # #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        # #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2]
        # #     # nodeIdsVentralBodySegmentBoundary = \
        # #     #     list(reversed(nodeIdsBodySegmentBoundary[
        # #     #                   newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        # #     #                   newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2]))
        # #     # print(nodeIdsVentralBodySegmentBoundary)
        # #
        # #     xLayer = []
        # #     d1Layer = []
        # #     d2Layer = []
        # #     for side in (dorsal, ventral):
        # #         if side == ventral:
        # #             xAlongApex.reverse()
        # #             xSideOfBodySegmentBoundary = xVentralBodySegmentBoundary
        # #             d2SideOfBodySegmentBoundary = d2VentralBodySegmentBoundary
        # #         else:
        # #             xSideOfBodySegmentBoundary = xDorsalBodySegmentBoundary
        # #             d2SideOfBodySegmentBoundary = d2DorsalBodySegmentBoundary
        # #         pxAlongAll = []
        # #         pd2AlongAll = []
        # #         for n1 in range(len(xSideOfBodySegmentBoundary)):
        # #             sf = 1.0 #2.0 if isPregnant else 1.5
        # #             nx = [xAlongApex[n1], xSideOfBodySegmentBoundary[n1]]
        # #             nd2 = [[0.0, 0.0, sf * halfFundusDepth * (1.0 if side == dorsal else -1.0)],
        # #                    d2SideOfBodySegmentBoundary[n1]]
        # #             px, pd2 = sampleCubicHermiteCurves(nx, nd2, int(0.25 * elementsCountAroundUterineTube))[:2]
        # #             pd2 = smoothCubicHermiteDerivativesLine(px, pd2)
        # #             pxAlongAll.append(px)
        # #             pd2AlongAll.append(pd2)
        # #
        # #         for n2 in range(int(0.25 * elementsCountAroundUterineTube)):
        # #             nx = []
        # #             nd1 = []
        # #             nd2 = []
        # #             for n1 in range(len(xDorsalBodySegmentBoundary)):
        # #                 nx.append(pxAlongAll[n1][n2])
        # #                 nd2.append(pd2AlongAll[n1][n2])
        # #                 if n1 == len(xDorsalBodySegmentBoundary) - 1:
        # #                     nd1.append(sub(pxAlongAll[n1][n2], pxAlongAll[n1 - 1][n2]))
        # #                 else:
        # #                     nd1.append(sub(pxAlongAll[n1 + 1][n2], pxAlongAll[n1][n2]))
        # #             nd1 = smoothCubicHermiteDerivativesLine(nx, nd1)
        # #             if side == ventral and n2 == 0:
        # #                 pass
        # #             else:
        # #                 xLayer.append(nx)
        # #                 d1Layer.append(nd1)
        # #                 d2Layer.append(nd2)
        # #
        # #     xAll.append(xLayer)
        # #     d1All.append(d1Layer)
        # #     d2All.append(d2Layer)
        # #
        # # # Re-arrange rows to get into correct layers and order from front to back
        # # xNewNodes = []
        # # d1NewNodes = []
        # # d2NewNodes = []
        # # for n3 in range(elementsCountThroughWall + 1):
        # #     nxLayer = []
        # #     nd1Layer = []
        # #     nd2Layer = []
        # #     for n2 in range(int(0.5 * elementsCountAroundUterineTube - 1)):
        # #         if n2 < int(0.25 * elementsCountAroundUterineTube):
        # #             n2Idx = int(0.25 * elementsCountAroundUterineTube - (n2 + 1))
        # #             # print('if,', n2Idx)
        # #         else:
        # #             n2Idx = n2
        # #             # print('else,', n2Idx)
        # #         nxLayer.append(xAll[n3][n2Idx])
        # #         nd1Layer.append(d1All[n3][n2Idx])
        # #         nd2Layer.append(d2All[n3][n2Idx])
        # #     xNewNodes.append(nxLayer)
        # #     d1NewNodes.append(nd1Layer)
        # #     d2NewNodes.append(nd2Layer)
        # #
        # # # for p2 in range(len(xNewNodes[0])):
        # # #     for p3 in range(elementsCountThroughWall + 1):
        # # #         for p1 in range(len(xNewNodes[0][p2])):
        # # #             node = nodes.createNode(nodeIdentifier, nodetemplate)
        # # #             fieldcache.setNode(node)
        # # #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, xNewNodes[p3][p2][p1])
        # # #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1NewNodes[p3][p2][p1])
        # # #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2NewNodes[p3][p2][p1])
        # # #             nodeIdentifier += 1
        # #
        # #     # Make nodes
        # #     nodeIDsNewNodes = []
        # #     frontBodyNodeIDs = []
        # #     backBodyNodeIDs = []
        # #     for layer in range(elementsCountThroughWall + 1):
        # #         frontBodyNodeIDs.append([nodeIdsRightSegmentBoundary[layer * (int(elementsCountAroundUterineTube * 0.5 + 1))]] +
        # #                             nodeIdsBodySegmentBoundary[newNodesCountAlongPatch * 2 * layer : newNodesCountAlongPatch * 2 * layer + newNodesCountAlongPatch] +
        # #                             [nodeIdsLeftSegmentBoundary[layer * int(elementsCountAroundUterineTube * 0.5 + 1)]])
        # #
        # #         backBodyNodeIDs.append(
        # #             [nodeIdsLeftSegmentBoundary[(layer + 1) * int(elementsCountAroundUterineTube * 0.5) + layer]] +
        # #              nodeIdsBodySegmentBoundary[newNodesCountAlongPatch * 2 * layer + newNodesCountAlongPatch : newNodesCountAlongPatch * 2 * layer + 2 * newNodesCountAlongPatch] +
        # #             [nodeIdsRightSegmentBoundary[(layer + 1) * int(elementsCountAroundUterineTube * 0.5) + layer]])
        # #
        # #     nodeIDsNewNodes.append(frontBodyNodeIDs)
        # #     for row in range(len(xNewNodes[0])):
        # #         nodeIDsLayer = []
        # #         for layer in range(len(xNewNodes)):
        # #             if row <= int(0.5 * (len(xNewNodes[0]) - 1)):
        # #                 nodeIDsAround = [nodeIdsRightSegmentBoundary[1 + row + layer * int(0.5 * elementsCountAroundUterineTube + 1)]]
        # #                 print(nodeIDsAround)
        # #             else:
        # #                 nodeIDsAround = [nodeIdsLeftSegmentBoundary[row + 1 + layer * int(0.5 * elementsCountAroundUterineTube + 1)]]
        # #             for col in range(len(xNewNodes[0][row])):
        # #                 colIdx = col #if row <= int(0.5 * (len(xNewNodes[0]) - 1)) else -(col + 1)
        # #                 node = nodes.createNode(nodeIdentifier, nodetemplate)
        # #                 fieldcache.setNode(node)
        # #                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, xNewNodes[layer][row][colIdx])
        # #                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1NewNodes[layer][row][colIdx])
        # #                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2NewNodes[layer][row][colIdx])
        # #                 if not useLinerThroughWall:
        # #                     if layer < len(xNewNodes) - 1:
        # #                         d3 = sub(xNewNodes[layer + 1][row][colIdx], xNewNodes[layer][row][colIdx])
        # #                     else:
        # #                         d3 = sub(xNewNodes[layer][row][colIdx], xNewNodes[layer - 1][row][colIdx])
        # #                     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
        # #                 nodeIDsAround.append(nodeIdentifier)
        # #                 nodeIdentifier += 1
        # #             if row <= int(0.5 * (len(xNewNodes[0]) - 1)):
        # #                 nodeIDsAround.append(nodeIdsLeftSegmentBoundary[1 + row + layer * int(0.5 * elementsCountAroundUterineTube + 1)])
        # #             else:
        # #                 nodeIDsAround.append(
        # #                     nodeIdsRightSegmentBoundary[1 + row + layer * int(0.5 * elementsCountAroundUterineTube + 1)])
        # #             nodeIDsLayer.append(nodeIDsAround)
        # #         nodeIDsNewNodes.append(nodeIDsLayer)
        # #     nodeIDsNewNodes.append(backBodyNodeIDs)
        #
        # # original
        # xNewNodes = []
        # d1NewNodes = []
        # d2NewNodes = []
        # for n3 in range(elementsCountThroughWall + 1):
        #     sf = (1.0 - innerProportionBody) / elementsCountThroughWall * n3 + innerProportionBody
        #     d2 = mult(d2BodyJunctionSegment1, sf)
        #     d3 = mult(d3BodyJunctionSegment1, sf)
        #     xMid, d1Mid = sampleEllipsePoints(xBodyJunction, d2, d3, 0.5 * math.pi, 1.5*math.pi,
        #                                       int(elementsCountAroundUterineTube * 0.5))
        #
        #     elementsCountAlongPatch = int((elementsCountAroundBody - elementsCountAroundUterineTube) * 0.5)
        #     newNodesCountAlongPatch = elementsCountAlongPatch - 1
        #     xRowList = [xBodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        #                                      newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]]
        #     d1RowList = [[zero for c in range(newNodesCountAlongPatch)]]
        #     d2RowList = [d2BodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        #                                        newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]]
        #     frontBodySegmentBoundaryNodes = nodeIdsBodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        #                                      newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]
        #     d2FrontBodySegmentBoundary = d2BodySegmentBoundary[newNodesCountAlongPatch * 2 * n3:
        #                                        newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch]
        #
        #     # print('nodesRightSegmentBoundary = ', nodeIdsRightSegmentBoundary)
        #     # print('len(xLeftSegmentBoundary) = ', len(xLeftSegmentBoundary))
        #     for n in range(1, len(xMid) - 1):
        #         idxAtSegment = n + int(((elementsCountAroundUterineTube * 0.5) + 1) * n3)
        #         nStart = xRightSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else xLeftSegmentBoundary[idxAtSegment]
        #         nEnd = xLeftSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else xRightSegmentBoundary[idxAtSegment]
        #         d2Start = d2RightSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else d2LeftSegmentBoundary[idxAtSegment]
        #         d2End = [-1 * c for c in d2LeftSegmentBoundary[idxAtSegment]] if n <= int(0.5 * len(xMid)) else [-1 * c for c in d2RightSegmentBoundary[idxAtSegment]]
        #         nodeIDStart = nodeIdsRightSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else nodeIdsLeftSegmentBoundary[
        #             idxAtSegment]
        #         nodeIDEnd = nodeIdsLeftSegmentBoundary[idxAtSegment] if n <= int(0.5 * len(xMid)) else nodeIdsRightSegmentBoundary[
        #             idxAtSegment]
        #         sf = 2.0 if isPregnant else 1.5
        #         nx = [nStart, xMid[n], nEnd]
        #         nd1 = [d2Start,
        #                [0.0, sf * (-halfFundusWidth if n <= int(0.5 * len(xMid)) else halfFundusWidth), 0.0],
        #                d2End]
        #
        #         px, pd1 = sampleCubicHermiteCurvesSmooth(nx, nd1, elementsCountAlongPatch)[:2]
        #         pd1 = smoothCubicHermiteDerivativesLine(px, pd1)
        #
        #         # if isPregnant: # CHECK!
        #         #     node = nodes.findNodeByIdentifier(nodeIDStart)
        #         #     fieldcache.setNode(node)
        #         #     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, pd1[0])
        #         #
        #         #     node = nodes.findNodeByIdentifier(nodeIDEnd)
        #         #     fieldcache.setNode(node)
        #         #     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, mult(pd1[-1], -1))
        #
        #         pxNew = px[1:-1]
        #         pd1New = pd1[1:-1]
        #
        #         if n > int(0.5*len(xMid)):
        #             pxNew.reverse()
        #             pd1New.reverse()
        #
        #         xRowList.append(pxNew)
        #         d1RowList.append(pd1New)
        #         d2RowList.append([[] for c in range(len(pxNew))])
        #
        #     backBodySegmentBoundaryNodes = \
        #         list(reversed(nodeIdsBodySegmentBoundary[
        #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2]))
        #     d2BackBodySegmentBoundary = \
        #         list(reversed(d2BodySegmentBoundary[
        #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        #                       newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2]))
        #     xRowList.append(
        #         list(reversed(xBodySegmentBoundary[newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        #                                            newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2])))
        #     d1RowList.append([zero for c in range(newNodesCountAlongPatch)])
        #     d2RowList.append(
        #         list(reversed(d2BodySegmentBoundary[newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch:
        #                                             newNodesCountAlongPatch * 2 * n3 + newNodesCountAlongPatch * 2])))
        #
        #     for n1 in range(newNodesCountAlongPatch):
        #         # nodes in front
        #         d2AlongAll = []
        #         xAlong = [xRowList[0][n1]]
        #         d2Along = [d2RowList[0][n1]]
        #         for n2 in range(1, int(0.5 * len(xRowList) + 1)):
        #             xAlong += [xRowList[n2][n1]]
        #             d2Along += [sub(xRowList[n2 - 1][n1], xRowList[n2][n1])]
        #         d2Along[-1] = set_magnitude([0.0, 0.0, 1.0], magnitude(d2Along[-1]))
        #         d2Along = smoothCubicHermiteDerivativesLine(list(reversed(xAlong)),
        #                                                     list(reversed(d2Along)), fixStartDirection=True)
        #         # Replace with end derivative
        #         node = nodes.findNodeByIdentifier(frontBodySegmentBoundaryNodes[n1])
        #         fieldcache.setNode(node)
        #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
        #                                       mult((add(d2Along[-1], d2FrontBodySegmentBoundary[n1])), 0.5))
        #
        #         d2AlongAll += list(reversed(d2Along))
        #         # nodes at the back
        #         xAlong = [xRowList[n2][n1]]
        #         d2AlongBack = [[-1 * c for c in d2Along[-1]]]
        #         for n2 in range(int(0.5 * len(xRowList) + 1), len(xRowList) - 1):
        #             xAlong += [xRowList[n2][n1]]
        #             d2AlongBack += [sub(xRowList[n2 + 1][n1], xRowList[n2][n1])]
        #         xAlong += [xRowList[-1][n1]]
        #         d2AlongBack += [d2RowList[-1][n1]]
        #         d2AlongBack = smoothCubicHermiteDerivativesLine(xAlong, d2AlongBack, fixStartDirection=True)
        #
        #         # Replace with end derivative
        #         node = nodes.findNodeByIdentifier(backBodySegmentBoundaryNodes[n1])
        #         fieldcache.setNode(node)
        #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1,
        #                                       mult((add(d2AlongBack[-1], d2BackBodySegmentBoundary[n1])), 0.5))
        #         d2AlongAll += d2AlongBack[1:]
        #
        #         for n2 in range(1, len(d2AlongAll) - 1):
        #             d2RowList[n2][n1] = d2AlongAll[n2]
        #
        #     xNewNodes.append(xRowList)
        #     d1NewNodes.append(d1RowList)
        #     d2NewNodes.append(d2RowList)
        #
        # # print('nodeIdsLeftSegmentBoundary', nodeIdsLeftSegmentBoundary)
        # # Make nodes
        # nodeIDsNewNodes = []
        # frontBodyNodeIDs = []
        # backBodyNodeIDs = []
        # for layer in range(elementsCountThroughWall + 1):
        #     frontBodyNodeIDs.append([nodeIdsRightSegmentBoundary[layer * (int(elementsCountAroundUterineTube * 0.5 + 1))]] +
        #                         nodeIdsBodySegmentBoundary[newNodesCountAlongPatch * 2 * layer : newNodesCountAlongPatch * 2 * layer + newNodesCountAlongPatch] +
        #                         [nodeIdsLeftSegmentBoundary[layer * int(elementsCountAroundUterineTube * 0.5 + 1)]])
        #
        #     backBodyNodeIDs.append(
        #         [nodeIdsLeftSegmentBoundary[(layer + 1) * int(elementsCountAroundUterineTube * 0.5) + layer]] +
        #          nodeIdsBodySegmentBoundary[newNodesCountAlongPatch * 2 * layer + newNodesCountAlongPatch : newNodesCountAlongPatch * 2 * layer + 2 * newNodesCountAlongPatch] +
        #         [nodeIdsRightSegmentBoundary[(layer + 1) * int(elementsCountAroundUterineTube * 0.5) + layer]])
        #
        # nodeIDsNewNodes.append(frontBodyNodeIDs)
        # for row in range(1, len(xNewNodes[0]) - 1):
        #     nodeIDsLayer = []
        #     for layer in range(len(xNewNodes)):
        #         if row <= int(0.5 * (len(xNewNodes[0]) - 1)):
        #             nodeIDsAround = [nodeIdsRightSegmentBoundary[row + layer * int(0.5 * elementsCountAroundUterineTube + 1)]]
        #         else:
        #             nodeIDsAround = [nodeIdsLeftSegmentBoundary[row + layer * int(0.5 * elementsCountAroundUterineTube + 1)]]
        #         for col in range(len(xNewNodes[0][row])):
        #             colIdx = col if row <= int(0.5 * (len(xNewNodes[0]) - 1)) else -(col + 1)
        #             node = nodes.createNode(nodeIdentifier, nodetemplate)
        #             fieldcache.setNode(node)
        #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, xNewNodes[layer][row][colIdx])
        #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1NewNodes[layer][row][colIdx])
        #             coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2NewNodes[layer][row][colIdx])
        #             if not useLinerThroughWall:
        #                 if layer < len(xNewNodes) - 1:
        #                     d3 = sub(xNewNodes[layer + 1][row][colIdx], xNewNodes[layer][row][colIdx])
        #                 else:
        #                     d3 = sub(xNewNodes[layer][row][colIdx], xNewNodes[layer - 1][row][colIdx])
        #                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
        #             nodeIDsAround.append(nodeIdentifier)
        #             nodeIdentifier += 1
        #         if row <= int(0.5 * (len(xNewNodes[0]) - 1)):
        #             nodeIDsAround.append(nodeIdsLeftSegmentBoundary[row + layer * int(0.5 * elementsCountAroundUterineTube + 1)])
        #         else:
        #             nodeIDsAround.append(
        #                 nodeIdsRightSegmentBoundary[row + layer * int(0.5 * elementsCountAroundUterineTube + 1)])
        #         nodeIDsLayer.append(nodeIDsAround)
        #     nodeIDsNewNodes.append(nodeIDsLayer)
        # nodeIDsNewNodes.append(backBodyNodeIDs)
        #
        # # print(nodeIDsNewNodes)
        # # for n2 in range(len(nodeIDsNewNodes)):
        # #     for n3 in range(len(nodeIDsNewNodes[n2])):
        # #         print(nodeIDsNewNodes[n2][n3])
        #
        # # print(len(nodeIDsNewNodes), len(nodeIDsNewNodes[0]), len(nodeIDsNewNodes[0][0]))
        #
        # leftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("left uterus"))
        # rightGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("right uterus"))
        # dorsalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("dorsal uterus"))
        # ventralGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("ventral uterus"))
        #
        # # Make elements
        # for e2 in range(len(nodeIDsNewNodes) - 1):
        #     for e3 in range(elementsCountThroughWall):
        #         for e1 in range(newNodesCountAlongPatch + 1):
        #             e1Idx = -(e1 + 1) if e2 == int(0.5 * len(nodeIDsNewNodes)) else e1
        #             e1IdxNext = -(e1 + 2) if e2 == int(0.5 * len(nodeIDsNewNodes)) else e1 + 1
        #             bni1 = nodeIDsNewNodes[(e2 + 1) if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2][e3][e1Idx]
        #             bni2 = nodeIDsNewNodes[(e2 + 1) if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2][e3][e1IdxNext]
        #             bni3 = nodeIDsNewNodes[e2 if e2 < 0.5*len(nodeIDsNewNodes) -1  else e2 + 1][e3][e1]
        #             bni4 = nodeIDsNewNodes[e2 if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2 + 1][e3][e1 + 1]
        #             bni5 = nodeIDsNewNodes[(e2 + 1) if e2 < 0.5*len(nodeIDsNewNodes) - 1  else e2][e3 + 1][e1Idx]
        #             bni6 = nodeIDsNewNodes[(e2 + 1) if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2][e3 + 1][e1IdxNext]
        #             bni7 = nodeIDsNewNodes[e2 if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2 + 1][e3 + 1][e1]
        #             bni8 = nodeIDsNewNodes[e2 if e2 < 0.5*len(nodeIDsNewNodes) - 1 else e2 + 1][e3 + 1][e1 + 1]
        #             nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
        #             elementtemplate1 = elementtemplate
        #             eft1 = eft
        #             scaleFactors = []
        #
        #             if e2 == 0:
        #                 if e1 == 0:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #                     # print('A', elementIdentifier) # , "-", nodeIdentifiers)
        #                 elif e1 == newNodesCountAlongPatch:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #                     # print('B', elementIdentifier)
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #                 # print('check', elementIdentifier, "-", nodeIdentifiers)
        #
        #             elif len(nodeIDsNewNodes) * 0.5 < 2:
        #                 scaleFactors = [-1.0]
        #                 eft1 = eftfactory.createEftNoCrossDerivatives()
        #                 setEftScaleFactorIds(eft1, [1], [])
        #                 if e1 == 0: # Right
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     # print('C', elementIdentifier)
        #                 elif e1 == newNodesCountAlongPatch: # Left
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #                     # print('D', elementIdentifier)
        #                 else:
        #                     remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     # print('E', elementIdentifier)
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #
        #             elif 0 < e2 <= int(len(nodeIDsNewNodes) * 0.5) - 1: # the one between body boundary and touching row in the middle
        #                 if e1 == 0: # Right
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     # print('F', elementIdentifier)
        #                 elif e1 == newNodesCountAlongPatch: # Left
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     # print('G', elementIdentifier)
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #             elif e2 == int(len(nodeIDsNewNodes) * 0.5):
        #                 # print(elementIdentifier, '-', nodeIdentifiers)
        #                 scaleFactors = [-1.0]
        #                 eft1 = eftfactory.createEftNoCrossDerivatives()
        #                 setEftScaleFactorIds(eft1, [1], [])
        #                 if e1 == 0:
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     # print('H', elementIdentifier, '-', nodeIdentifiers)
        #                 elif e1 == newNodesCountAlongPatch:
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     # print('I', elementIdentifier, nodeIdentifiers)
        #                 else:
        #                     remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     # print('J', elementIdentifier)
        #
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #
        #             elif len(nodeIDsNewNodes) * 0.5 > 2 and int(0.5 * len(nodeIDsNewNodes)) < e2 < len(nodeIDsNewNodes)-2:
        #                 if e1 == 0:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     # print('K', elementIdentifier)
        #
        #                 elif e1 == newNodesCountAlongPatch:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     # print('L', elementIdentifier)
        #
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #
        #             elif e2 == len(nodeIDsNewNodes) - 2:
        #                 if e1 == 0:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [1])])
        #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #
        #                 elif e1 == newNodesCountAlongPatch:
        #                     scaleFactors = [-1.0]
        #                     eft1 = eftfactory.createEftNoCrossDerivatives()
        #                     setEftScaleFactorIds(eft1, [1], [])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, [])])
        #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
        #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
        #                                            [(Node.VALUE_LABEL_D_DS2, [1])])
        #
        #                 elementtemplateX.defineField(coordinates, -1, eft1)
        #                 elementtemplate1 = elementtemplateX
        #
        #             element = mesh.createElement(elementIdentifier, elementtemplate1)
        #             result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        #             if scaleFactors:
        #                 element.setScaleFactors(eft1, scaleFactors)
        #             groups = [uterusGroup, fundusGroup]
        #             if e2 < int(0.5 * len(nodeIDsNewNodes)):
        #                 groups.append(dorsalGroup)
        #                 if e1 < newNodesCountAlongPatch * 0.5:
        #                     groups.append(rightGroup)
        #                 else:
        #                     groups.append(leftGroup)
        #             else:
        #                 groups.append(ventralGroup)
        #                 if e1 < newNodesCountAlongPatch * 0.5:
        #                     groups.append(leftGroup)
        #                 else:
        #                     groups.append(rightGroup)
        #
        #             for annotationGroup in groups:
        #                 meshGroup = annotationGroup.getMeshGroup(mesh)
        #                 meshGroup.addElement(element)
        #             elementIdentifier = elementIdentifier + 1
        #
        # # annotationGroups.remove(leftSegmentBoundary)
        # # annotationGroups.remove(rightSegmentBoundary)
        # # annotationGroups.remove(bodySegmentBoundary)
        # # annotationGroups.remove(fundusPatchGroup)

        myometriumGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_uterus_term("myometrium"))
        myometriumGroup.getMeshGroup(mesh).addElementsConditional(uterusGroup.getGroup())

        return annotationGroups, None

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
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ("Default", "Human 1", "Human Pregnant 1")
        isRodent = parameterSetName in ("Mouse 1", "Rat 1")

        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        uterusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterus"))
        leftOviductGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left oviduct"))
        rightOviductGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right oviduct"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("body of uterus"))
        fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("fundus of uterus"))
        cervixGroup = getAnnotationGroupForTerm(annotationGroups, ("cervix", ""))
        upperCervixGroup = getAnnotationGroupForTerm(annotationGroups, ("upper cervix", ""))
        lowerCervixGroup = getAnnotationGroupForTerm(annotationGroups, ("lower cervix", ""))
        vaginaGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("vagina"))

        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
        is_exterior_face_xi2_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))

        is_uterus = uterusGroup.getGroup()
        is_uterus_outer = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_1)
        is_uterus_inner = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_0)

        is_leftOviduct = leftOviductGroup.getGroup()
        is_rightOviduct = rightOviductGroup.getGroup()

        is_fundus = fundusGroup.getGroup()
        is_fundus_outer = fm.createFieldAnd(is_fundus, is_exterior_face_xi3_1)
        is_fundus_inner = fm.createFieldAnd(is_fundus, is_exterior_face_xi3_0)

        is_body = bodyGroup.getGroup()
        is_body_outer = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_inner = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_cervix = cervixGroup.getGroup()
        is_cervix_outer = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_1)
        is_cervix_inner = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_0)

        is_upperCervix = upperCervixGroup.getGroup()
        is_upperCervix_outer = fm.createFieldAnd(is_upperCervix, is_exterior_face_xi3_1)

        is_lowerCervix = lowerCervixGroup.getGroup()
        is_lowerCervix_outer = fm.createFieldAnd(is_lowerCervix, is_exterior_face_xi3_1)

        is_vagina = vaginaGroup.getGroup()
        is_vagina_outer = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_1)
        is_vagina_inner = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_0)
        is_vagina_xi2_0 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_0)
        is_vagina_xi2_1 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_1)
        is_vagina_xi2_01 = fm.createFieldXor(is_vagina_xi2_0, is_vagina_xi2_1)

        serosaOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterus"))
        serosaOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_outer)

        uterineCavity = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("uterine cavity"))
        uterineCavity.getMeshGroup(mesh2d).addElementsConditional(is_uterus_inner)
        uterineCavity.getMeshGroup(mesh2d).removeElementsConditional(is_leftOviduct)
        uterineCavity.getMeshGroup(mesh2d).removeElementsConditional(is_rightOviduct)

        serosaOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                          get_uterus_term("serosa of body of uterus"))
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_outer)
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_upperCervix_outer)

        lumenOfFundus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_uterus_term("lumen of fundus of uterus"))
        lumenOfFundus.getMeshGroup(mesh2d).addElementsConditional(is_fundus_inner)

        serosaOfFundus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                          get_uterus_term("serosa of fundus of uterus"))
        serosaOfFundus.getMeshGroup(mesh2d).addElementsConditional(is_fundus_outer)

        lumenOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_uterus_term("lumen of body of uterus"))
        lumenOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_inner)

        lumenOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("uterine cervix"))
        lumenOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_inner)

        serosaOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of vagina"))
        serosaOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_outer)
        serosaOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_lowerCervix_outer)

        lumenOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("vaginal canal"))
        lumenOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_inner)

        leftGroup = getAnnotationGroupForTerm(annotationGroups, ("left uterus", ""))
        rightGroup = getAnnotationGroupForTerm(annotationGroups, ("right uterus", ""))
        dorsalGroup = getAnnotationGroupForTerm(annotationGroups, ("dorsal uterus", ""))
        ventralGroup = getAnnotationGroupForTerm(annotationGroups, ("ventral uterus", ""))

        # leftFundusGroup = getAnnotationGroupForTerm(annotationGroups, ("left fundus", "None"))
        # rightFundusGroup = getAnnotationGroupForTerm(annotationGroups, ("right fundus", "None"))

        if isHuman:
            leftOviductGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left oviduct"))
            is_leftOviduct = leftOviductGroup.getGroup()
            is_leftOviduct_outer = fm.createFieldAnd(is_leftOviduct, is_exterior_face_xi3_1)
            is_leftOviduct_inner = fm.createFieldAnd(is_leftOviduct, is_exterior_face_xi3_0)

            serosaOfLeftOviduct = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of left oviduct"))
            serosaOfLeftOviduct.getMeshGroup(mesh2d).addElementsConditional(is_leftOviduct_outer)

            rightOviductGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right oviduct"))
            is_rightOviduct = rightOviductGroup.getGroup()
            is_rightOviduct_outer = fm.createFieldAnd(is_rightOviduct, is_exterior_face_xi3_1)
            is_rightOviduct_inner = fm.createFieldAnd(is_rightOviduct, is_exterior_face_xi3_0)

            serosaOfRightOviduct = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                         get_uterus_term("serosa of right oviduct"))
            serosaOfRightOviduct.getMeshGroup(mesh2d).addElementsConditional(is_rightOviduct_outer)

            lumenOfLeftOviduct = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of left oviduct"))
            lumenOfLeftOviduct.getMeshGroup(mesh2d).addElementsConditional(is_leftOviduct_inner)

            lumenOfRightOviduct = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                        get_uterus_term("lumen of right oviduct"))
            lumenOfRightOviduct.getMeshGroup(mesh2d).addElementsConditional(is_rightOviduct_inner)

            is_pubocervical = fm.createFieldAnd(is_body_outer, is_cervix_outer)
            pubocervical = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("pubocervical ligament"))
            pubocervical.getMeshGroup(mesh1d).addElementsConditional(is_pubocervical)

            is_internal_os = fm.createFieldAnd(is_body_inner, is_cervix_inner)
            internalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_uterus_term("internal cervical os"))
            internalOs.getMeshGroup(mesh1d).addElementsConditional(is_internal_os)

            is_external_os = fm.createFieldAnd(is_vagina_inner, is_cervix_inner)
            externalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("external cervical os"))
            externalOs.getMeshGroup(mesh1d).addElementsConditional(is_external_os)

            is_vagina_orifice = fm.createFieldAnd(is_vagina_xi2_01, is_exterior_face_xi3_0)
            vaginaOrifice = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("vagina orifice"))
            vaginaOrifice.getMeshGroup(mesh1d).addElementsConditional(is_vagina_orifice)

            # ligaments
            is_dorsalVentral = fm.createFieldAnd(dorsalGroup.getGroup(), ventralGroup.getGroup())
            is_dorsalVentralSerosa = fm.createFieldAnd(is_dorsalVentral, is_exterior_face_xi3_1)
            is_leftDorsalVentralSerosa = fm.createFieldAnd(leftGroup.getGroup(), is_dorsalVentralSerosa)
            is_rightDorsalVentralSerosa = fm.createFieldAnd(rightGroup.getGroup(), is_dorsalVentralSerosa)
            fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("fundus of uterus"))
            is_bodyNotFundus = fm.createFieldAnd(bodyGroup.getGroup(), fm.createFieldNot(fundusGroup.getGroup()))

            # Broad ligament of uterus
            is_leftBroadLigament = fm.createFieldAnd(is_bodyNotFundus, is_leftDorsalVentralSerosa)
            leftBroadLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("left broad ligament of uterus"))
            leftBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_leftBroadLigament)
            # add connected edges from left oviduct, avoiding adding dorsal-ventral edges on the superior edge
            leftBroadLigament.addSubelements()  # need current nodes in ligament for group_add_connected_elements
            tmpGroup = fm.createFieldGroup()
            tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
            tmpMeshGroup.addElementsConditional(fm.createFieldAnd(is_leftOviduct, is_leftDorsalVentralSerosa))
            group_add_connected_elements(leftBroadLigament.getGroup(), tmpMeshGroup)
            del tmpMeshGroup
            del tmpGroup

            is_rightBroadLigament = fm.createFieldAnd(is_bodyNotFundus, is_rightDorsalVentralSerosa)
            rightBroadLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("right broad ligament of uterus"))
            rightBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_rightBroadLigament)
            # add connected edges from right oviduct, avoiding adding dorsal-ventral edges on the superior edge
            rightBroadLigament.addSubelements()  # need current nodes in ligament for group_add_connected_elements
            tmpGroup = fm.createFieldGroup()
            tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
            tmpMeshGroup.addElementsConditional(fm.createFieldAnd(is_rightOviduct, is_rightDorsalVentralSerosa))
            group_add_connected_elements(rightBroadLigament.getGroup(), tmpMeshGroup)
            del tmpMeshGroup
            del tmpGroup

            # Transverse cervical ligament
            is_leftTransverseCervicalLigament = fm.createFieldAnd(cervixGroup.getGroup(), is_leftDorsalVentralSerosa)
            leftTransverseCervicalLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("left transverse cervical ligament"))
            leftTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(
                is_leftTransverseCervicalLigament)

            is_rightTransverseCervicalLigament = fm.createFieldAnd(cervixGroup.getGroup(), is_rightDorsalVentralSerosa)
            rightTransverseCervicalLigament = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_uterus_term("right transverse cervical ligament"))
            rightTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(
                is_rightTransverseCervicalLigament)

        if isRodent:
            rightHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right uterine horn"))
            leftHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left uterine horn"))

            is_rightHorn = rightHornGroup.getGroup()
            is_rightHorn_outer = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_1)
            is_rightHorn_inner = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_0)

            is_leftHorn = leftHornGroup.getGroup()
            is_leftHorn_outer = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_1)
            is_leftHorn_inner = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_0)

            serosaOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of right horn"))
            serosaOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_outer)

            lumenOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of right horn"))
            lumenOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_inner)

            serosaOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("serosa of left horn"))
            serosaOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_outer)

            lumenOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_uterus_term("lumen of left horn"))
            lumenOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_inner)

        # annotationGroups.remove(leftFundusGroup)
        # annotationGroups.remove(rightFundusGroup)
        annotationGroups.remove(cervixGroup)
        annotationGroups.remove(upperCervixGroup)
        annotationGroups.remove(lowerCervixGroup)

def setNodeFieldParameters(field, fieldcache, x, d1, d2, d3, d12=None, d13=None):
    """
    Assign node field parameters x, d1, d2, d3 of field.
    :param field: Field parameters to assign.
    :param fieldcache: Fieldcache with node set.
    :param x: Parameters to set for Node.VALUE_LABEL_VALUE.
    :param d1: Parameters to set for Node.VALUE_LABEL_D_DS1.
    :param d2: Parameters to set for Node.VALUE_LABEL_D_DS2.
    :param d3: Parameters to set for Node.VALUE_LABEL_D_DS3.
    :param d12: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS2.
    :param d13: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS3.
    :return:
    """
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
    if d12:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d12)
    if d13:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, d13)


def setNodeFieldVersionDerivatives(field, fieldcache, version, d1, d2, d3, d12=None, d13=None):
    """
    Assign node field parameters d1, d2, d3 of field.
    :param field: Field to assign parameters of.
    :param fieldcache: Fieldcache with node set.
    :param version: Version of d1, d2, d3 >= 1.
    :param d1: Parameters to set for Node.VALUE_LABEL_D_DS1.
    :param d2: Parameters to set for Node.VALUE_LABEL_D_DS2.
    :param d3: Parameters to set for Node.VALUE_LABEL_D_DS3.
    :param d12: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS2.
    :param d13: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS3.
    :return:
    """
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, d1)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, d2)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, d3)
    if d12:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, d12)
    if d13:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, d13)

def reorderNodesForMapping(xList, elementsCountThroughWall):
    """

    """
    elementsCountPerLayer = len(xList) // (elementsCountThroughWall + 1)
    xReordered = []
    for n3 in range(elementsCountThroughWall + 1):
        xLayer = xList[elementsCountPerLayer * n3:elementsCountPerLayer * (n3 + 1)]
        xLayerReversed = list(reversed(xLayer))
        xLayerOrdered = xLayerReversed[len(xLayerReversed) // 2:] + \
                        xLayerReversed[0: len(xLayerReversed) // 2]
        xReordered += xLayerOrdered

    return xReordered


