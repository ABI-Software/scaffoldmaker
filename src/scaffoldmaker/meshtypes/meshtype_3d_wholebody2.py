"""
Generates a 3D body coordinates using tube network mesh.
"""
from cmlibs.maths.vectorops import cross, mult, set_magnitude
from cmlibs.utils.zinc.field import Field, find_or_create_field_coordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.utils.interpolation import (
    computeCubicHermiteEndDerivative, interpolateLagrangeHermiteDerivative, sampleCubicHermiteCurvesSmooth)
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
import math

class MeshType_1d_human_body_network_layout1(MeshType_1d_network_layout1):
    """
    Defines body network layout.
    """

    @classmethod
    def getName(cls):
        return "1D Human Body Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = parameterSetName
        options["Structure"] = (
            "1-2-3-4,"
            "4-5-6.1," 
            "6.2-14-15-16-17-18-19,19-20,"
            "6.3-21-22-23-24-25-26,26-27,"
            "6.1-7-8-9,"
            "9-10-11-12-13.1,"
            "13.2-28-29-30,30-31-32-33,33-34,"
            "13.3-35-36-37,37-38-39-40,40-41")
        options["Define inner coordinates"] = True
        options["Head depth"] = 2.0
        options["Head length"] = 2.5
        options["Head width"] = 2.0
        options["Neck length"] = 1.5
        options["Shoulder drop"] = 0.7
        options["Shoulder width"] = 5.0
        options["Arm lateral angle degrees"] = 10.0
        options["Arm length"] = 7.0
        options["Arm top diameter"] = 1.0
        options["Wrist thickness"] = 0.5
        options["Wrist width"] = 0.7
        options["Hand length"] = 1.5
        options["Hand thickness"] = 0.3
        options["Hand width"] = 1.0
        options["Torso depth"] = 2.5
        options["Torso length"] = 6.0
        options["Torso width"] = 3.5
        options["Pelvis drop"] = 1.0
        options["Pelvis width"] = 2.0
        options["Leg lateral angle degrees"] = 10.0
        options["Leg length"] = 9.0
        options["Leg top diameter"] = 1.75
        options["Leg bottom diameter"] = 0.75
        options["Foot length"] = 2.0
        options["Foot thickness"] = 0.4
        options["Foot width"] = 1.2
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Head depth",
            "Head length",
            "Head width",
            "Neck length",
            "Shoulder drop",
            "Shoulder width",
            "Arm lateral angle degrees",
            "Arm length",
            "Arm top diameter",
            "Wrist thickness",
            "Wrist width",
            "Hand length",
            "Hand thickness",
            "Hand width",
            "Torso depth",
            "Torso length",
            "Torso width",
            "Pelvis drop",
            "Pelvis width",
            "Leg lateral angle degrees",
            "Leg length",
            "Leg top diameter",
            "Leg bottom diameter",
            "Foot length",
            "Foot thickness",
            "Foot width"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Head depth",
            "Head length",
            "Head width",
            "Neck length",
            "Shoulder drop",
            "Shoulder width",
            "Arm length",
            "Arm top diameter",
            "Wrist thickness",
            "Wrist width",
            "Hand length",
            "Hand thickness",
            "Hand width",
            "Pelvis drop",
            "Pelvis width",
            "Torso depth",
            "Torso width",
            "Leg length",
            "Leg top diameter",
            "Leg bottom diameter",
            "Foot length",
            "Foot thickness",
            "Foot width"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
        for key in [
            "Arm lateral angle degrees"
        ]:
            if options[key] < -60.0:
                options[key] = -60.0
            elif options[key] > 200.0:
                options[key] = 200.0
        for key in [
            "Leg lateral angle degrees"
        ]:
            if options[key] < -20.0:
                options[key] = -20.0
            elif options[key] > 60.0:
                options[key] = 60.0
        return dependentChanges


    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, NetworkMesh
        """
        # parameterSetName = options['Base parameter set']
        structure = options["Structure"]
        headDepth = options["Head depth"]
        headLength = options["Head length"]
        headWidth = options["Head width"]
        neckLength = options["Neck length"]
        shoulderDrop = options["Shoulder drop"]
        halfShoulderWidth = 0.5 * options["Shoulder width"]
        armAngleRadians = math.radians(options["Arm lateral angle degrees"])
        armAngleDrop = shoulderDrop * math.cos(armAngleRadians)
        armLength = options["Arm length"]
        armTopRadius = 0.5 * options["Arm top diameter"]
        halfWristThickness = 0.5 * options["Wrist thickness"]
        halfWristWidth = 0.5 * options["Wrist width"]
        handLength = options["Hand length"]
        halfHandThickness = 0.5 * options["Hand thickness"]
        halfHandWidth = 0.5 * options["Hand width"]
        halfTorsoDepth = 0.5 * options["Torso depth"]
        torsoLength = options["Torso length"]
        halfTorsoWidth = 0.5 * options["Torso width"]
        pelvisDrop = options["Pelvis drop"]
        halfPelvisWidth = 0.5 * options["Pelvis width"]
        legAngleRadians = math.radians(options["Leg lateral angle degrees"])
        legLength = options["Leg length"]
        legTopRadius = 0.5 * options["Leg top diameter"]
        legBottomRadius = 0.5 * options["Leg bottom diameter"]

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)

        # set up element annotations
        bodyGroup = AnnotationGroup(region, get_body_term("body"))
        headGroup = AnnotationGroup(region, get_body_term("head"))
        neckGroup = AnnotationGroup(region, get_body_term("neck"))
        armGroup = AnnotationGroup(region, get_body_term("arm"))
        leftArmGroup = AnnotationGroup(region, get_body_term("left arm"))
        rightArmGroup = AnnotationGroup(region, get_body_term("right arm"))
        handGroup = AnnotationGroup(region, get_body_term("hand"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        legGroup = AnnotationGroup(region, get_body_term("leg"))
        leftLegGroup = AnnotationGroup(region, get_body_term("left leg"))
        rightLegGroup = AnnotationGroup(region, get_body_term("right leg"))
        annotationGroups = [bodyGroup, headGroup, neckGroup, armGroup, leftArmGroup, rightArmGroup, handGroup,
                            thoraxGroup, abdomenGroup, legGroup, leftLegGroup, rightLegGroup]
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        elementIdentifier = 1
        headElementsCount = 3
        meshGroups = [bodyMeshGroup, headGroup.getMeshGroup(mesh)]
        for e in range(headElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        neckElementsCount = 2
        meshGroups = [bodyMeshGroup, neckGroup.getMeshGroup(mesh)]
        for e in range(neckElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        left = 0
        right = 1
        armElementsCount = 7
        armMeshGroup = armGroup.getMeshGroup(mesh)
        handMeshGroup = handGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideArmGroup = leftArmGroup if (side == left) else rightArmGroup
            meshGroups = [bodyMeshGroup, armMeshGroup, sideArmGroup.getMeshGroup(mesh)]
            for e in range(armElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                if e == (armElementsCount - 1):
                    handMeshGroup.addElement(element)
                elementIdentifier += 1
        thoraxElementsCount = 3
        abdomenElementsCount = 4
        torsoElementsCount = thoraxElementsCount + abdomenElementsCount
        meshGroups = [bodyMeshGroup, thoraxGroup.getMeshGroup(mesh)]
        for e in range(thoraxElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        meshGroups = [bodyMeshGroup, abdomenGroup.getMeshGroup(mesh)]
        for e in range(abdomenElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        legElementsCount = 7
        legMeshGroup = legGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideLegGroup = leftLegGroup if (side == left) else rightLegGroup
            meshGroups = [bodyMeshGroup, legMeshGroup, sideLegGroup.getMeshGroup(mesh)]
            for e in range(legElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = find_or_create_field_coordinates(fieldmodule).castFiniteElement()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        headScale = headLength / headElementsCount
        nodeIdentifier = 1
        d1 = [headScale, 0.0, 0.0]
        d2 = [0.0, 0.5 * headWidth, 0.0]
        d3 = [0.0, 0.0, 0.5 * headDepth]
        for i in range(headElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, [headScale * i, 0.0, 0.0])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
            nodeIdentifier += 1
        neckScale = neckLength / neckElementsCount
        d1 = [neckScale, 0.0, 0.0]
        d2 = [0.0, 0.5 * headWidth, 0.0]
        d3 = [0.0, 0.0, 0.5 * headWidth]
        for i in range(1, neckElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [headLength + neckScale * i, 0.0, 0.0]
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
            nodeIdentifier += 1
        armJunctionNodeIdentifier = nodeIdentifier
        torsoScale = torsoLength / torsoElementsCount
        d1 = [torsoScale, 0.0, 0.0]
        d2 = [0.0, halfTorsoWidth, 0.0]
        d3 = [0.0, 0.0, halfTorsoDepth]
        torsoStartX = headLength + neckLength
        sx = [torsoStartX, 0.0, 0.0]
        for i in range(torsoElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [torsoStartX + torsoScale * i, 0.0, 0.0]
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
            nodeIdentifier += 1
        legJunctionNodeIdentifier = nodeIdentifier - 1
        px = x

        # arms
        armStartX = torsoStartX + armAngleDrop
        nonHandArmLength = armLength - handLength
        armScale = nonHandArmLength / (armElementsCount - 3)
        sd3 = [0.0, 0.0, armTopRadius]
        hd3 = [0.0, 0.0, halfHandWidth]
        for side in (left, right):
            armAngle = armAngleRadians if (side == left) else -armAngleRadians
            cosArmAngle = math.cos(armAngle)
            sinArmAngle = math.sin(armAngle)
            armStartY = halfShoulderWidth if (side == left) else -halfShoulderWidth
            x = [armStartX, armStartY, 0.0]
            d1 = [armScale * cosArmAngle, armScale * sinArmAngle, 0.0]
            # set leg versions 2 (left) and 3 (right) on leg junction node, and intermediate shoulder node
            sd1 = interpolateLagrangeHermiteDerivative(sx, x, d1, 0.0)
            nx, nd1 = sampleCubicHermiteCurvesSmooth([sx, x], [sd1, d1], 2, derivativeMagnitudeEnd=armScale)[0:2]
            for n in range(2):
                node = nodes.findNodeByIdentifier(nodeIdentifier if (n > 0) else armJunctionNodeIdentifier)
                fieldcache.setNode(node)
                version = 1 if (n > 0) else 2 if (side == left) else 3
                sd1 = nd1[n]
                sd2 = set_magnitude(cross(sd3, sd1), armTopRadius)
                if n > 0:
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, version, nx[n])
                    nodeIdentifier += 1
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, sd1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, sd2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, sd3)
            # main part of arm to wrist
            for i in range(armElementsCount - 2):
                xi = i / (armElementsCount - 3)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = [armStartX + d1[0] * i, armStartY + d1[1] * i, d1[2] * i]
                halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
                halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
                d2 = [-halfThickness * sinArmAngle, halfThickness * cosArmAngle, 0.0]
                d3 = [0.0, 0.0, halfWidth]
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                nodeIdentifier += 1
            # hand
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            h = (armElementsCount - 3) + handLength / armScale
            hx = [armStartX + armLength * cosArmAngle, armStartY + armLength * sinArmAngle, 0.0]
            hd1 = computeCubicHermiteEndDerivative(x, d1, hx, d1)
            hd2 = set_magnitude(d2, halfHandThickness)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, hx)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, hd1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, hd2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, hd3)
            nodeIdentifier += 1

        # legs
        legStartX = torsoStartX + torsoLength + pelvisDrop
        legScale = legLength / (legElementsCount - 1)
        pd3 = [0.0, 0.0, halfTorsoDepth]
        for side in (left, right):
            legAngle = legAngleRadians if (side == left) else -legAngleRadians
            coslegAngle = math.cos(legAngle)
            sinlegAngle = math.sin(legAngle)
            legStartY = halfPelvisWidth if (side == left) else -halfPelvisWidth
            x = [legStartX, legStartY, 0.0]
            d1 = [legScale * coslegAngle, legScale * sinlegAngle, 0.0]

            # set leg versions 2 (left) and 3 (right) on leg junction node
            node = nodes.findNodeByIdentifier(legJunctionNodeIdentifier)
            fieldcache.setNode(node)
            pd1 = interpolateLagrangeHermiteDerivative(px, x, d1, 0.0)
            pd2 = set_magnitude(cross(pd3, pd1), legTopRadius)
            version = 2 if (side == left) else 3
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, pd1)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, pd2)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, pd3)

            for i in range(legElementsCount):
                xi = i / (legElementsCount - 1)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = [legStartX + d1[0] * i, legStartY + d1[1] * i, d1[2] * i]
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = [-radius * sinlegAngle, radius * coslegAngle, 0.0]
                d3 = [0.0, 0.0, radius]
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                nodeIdentifier += 1

        smoothOptions = {
            "Field": {"coordinates": True, "inner coordinates": False},
            "Smooth D12": True,
            "Smooth D13": True}
        cls.smoothSideCrossDerivatives(region, options, networkMesh, smoothOptions, None)

        cls.defineInnerCoordinates(region, coordinates, options, networkMesh)

        return annotationGroups, networkMesh

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Edit base class list to include only valid functions.
        """
        interactiveFunctions = super(MeshType_1d_human_body_network_layout1, cls).getInteractiveFunctions()
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Edit structure...":
                interactiveFunctions.remove(interactiveFunction)
                break
        return interactiveFunctions

class MeshType_3d_wholebody2(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network with core representing the human body.
    """

    @classmethod
    def getName(cls):
        return "3D Whole Body 2"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1 Coarse",
            "Human 1 Medium",
            "Human 1 Fine"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        useParameterSetName = "Human 1 Coarse" if (parameterSetName == "Default") else parameterSetName
        options["Base parameter set"] = useParameterSetName
        options["Body network layout"] = ScaffoldPackage(MeshType_1d_human_body_network_layout1)
        options["Number of elements around head"] = 12
        options["Number of elements around torso"] = 12
        options["Number of elements around arm"] = 8
        options["Number of elements around leg"] = 8
        options["Number of elements through shell"] = 1
        options["Target element density along longest segment"] = 5.0
        options["Show trim surfaces"] = False
        options["Use Core"] = True
        options["Number of elements across core box minor"] = 2
        options["Number of elements across core transition"] = 1
        if "Medium" in useParameterSetName:
            options["Number of elements around head"] = 16
            options["Number of elements around torso"] = 16
            options["Number of elements around leg"] = 12
            options["Target element density along longest segment"] = 8.0
        elif "Fine" in useParameterSetName:
            options["Number of elements around head"] = 24
            options["Number of elements around torso"] = 24
            options["Number of elements around arm"] = 12
            options["Number of elements around leg"] = 16
            options["Number of elements through shell"] = 1
            options["Target element density along longest segment"] = 10.0
            options["Number of elements across core box minor"] = 4

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Body network layout",
            "Number of elements around head",
            "Number of elements around torso",
            "Number of elements around arm",
            "Number of elements around leg",
            "Number of elements through shell",
            "Target element density along longest segment",
            "Show trim surfaces",
            "Use Core",
            "Number of elements across core box minor",
            "Number of elements across core transition"]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Body network layout":
            return [MeshType_1d_human_body_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                "Invalid parameter set " + str(parameterSetName) + " for scaffold " + str(scaffoldType.getName()) + \
                " in option " + str(optionName) + " of scaffold " + cls.getName()
        if optionName == "Body network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_human_body_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if not options["Body network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Body network layout"):
            options["Body network layout"] = ScaffoldPackage(MeshType_1d_human_body_network_layout1)
        minElementsCountAround = None
        for key in [
            "Number of elements around head",
            "Number of elements around torso",
            "Number of elements around arm",
            "Number of elements around leg"
        ]:
            if options[key] < 8:
                options[key] = 8
            elif options[key] % 4:
                options[key] += 4 - (options[key] % 4)
            if (minElementsCountAround is None) or (options[key] < minElementsCountAround):
                minElementsCountAround = options[key]

        if options["Number of elements through shell"] < 0:
            options["Number of elements through shell"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0

        if options["Number of elements across core transition"] < 1:
            options["Number of elements across core transition"] = 1

        maxElementsCountCoreBoxMinor = minElementsCountAround // 2 - 2
        for key in [
            "Number of elements across core box minor"
        ]:
            if options[key] < 2:
                options[key] = 2
            elif options[key] > maxElementsCountCoreBoxMinor:
                options[key] = maxElementsCountCoreBoxMinor
                dependentChanges = True
            elif options[key] % 2:
                options[key] += options[key] % 2

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        # parameterSetName = options['Base parameter set']
        layoutRegion = region.createRegion()
        networkLayout = options["Body network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        coreBoxMinorCount = options["Number of elements across core box minor"]
        coreTransitionCount = options['Number of elements across core transition']
        annotationAroundCounts = []
        # implementation currently uses major count including transition
        annotationCoreMajorCounts = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            aroundCount = 0
            coreMajorCount = 0
            name = layoutAnnotationGroup.getName()
            if ("head" in name) or ("neck" in name):
                aroundCount = options["Number of elements around head"]
            elif ("abdomen" in name) or ("thorax" in name) or ("torso" in name):
                aroundCount = options["Number of elements around torso"]
            elif "arm" in name:
                aroundCount = options["Number of elements around arm"]
            elif "leg" in name:
                aroundCount = options["Number of elements around leg"]
            if aroundCount:
                coreMajorCount = aroundCount // 2 - coreBoxMinorCount + 2 * coreTransitionCount
            annotationAroundCounts.append(aroundCount)
            annotationCoreMajorCounts.append(coreMajorCount)
        isCore = options["Use Core"]

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Number of elements around head"],
            elementsCountThroughWall=options["Number of elements through shell"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=annotationAroundCounts,
            defaultElementsCountAcrossMajor=annotationCoreMajorCounts[-1],
            elementsCountTransition=coreTransitionCount,
            annotationElementsCountsAcrossMajor=annotationCoreMajorCounts,
            isCore=isCore)

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None

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
        mesh2d = fm.findMeshByDimension(2)

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))

        is_exterior = fm.createFieldIsExterior()
        is_on_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_skin = fm.createFieldAnd(is_exterior, is_on_face_xi3_1)

        skinMeshGroup = skinGroup.getMeshGroup(mesh2d)
        skinMeshGroup.addElementsConditional(is_skin)
