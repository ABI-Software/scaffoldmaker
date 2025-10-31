"""
Generates a 3D body coordinates using tube network mesh.
"""
from cmlibs.maths.vectorops import add, cross, mult, set_magnitude, sub, magnitude, axis_angle_to_rotation_matrix, matrix_vector_mult, dot
from cmlibs.utils.zinc.field import Field, find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_maximum_node_identifier
from cmlibs.zinc.element import Element
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import (
    AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm, evaluateAnnotationMarkerNearestMeshLocation)
from scaffoldmaker.annotation.body_terms import get_body_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.interpolation import (
    computeCubicHermiteEndDerivative, getCubicHermiteArcLength, interpolateLagrangeHermiteDerivative,
    sampleCubicHermiteCurvesSmooth, smoothCubicHermiteDerivativesLine)
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import BodyTubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.human_network_layout import constructNetworkLayoutStructure, humanElementCounts
# from scaffoldmaker.utils.human_network_layout import 
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
        options["Structure"] = constructNetworkLayoutStructure(humanElementCounts)
        options["Define inner coordinates"] = True
        options["Head depth"] = 2.0
        options["Head length"] = 2.2
        options["Head width"] = 2.0
        options["Neck length"] = 1.3
        options["Shoulder drop"] = 1.0
        options["Shoulder width"] = 4.5
        options["Left shoulder flexion degrees"] = 0.0
        options["Right shoulder flexion degrees"] = 0.0
        options["Left shoulder abduction degrees"] = 10.0
        options["Right shoulder abduction degrees"] = 10.0
        options["Left elbow flexion degrees"] = 110.0
        options["Right elbow flexion degrees"] = 90.0
        options["Arm length"] = 7.5
        options["Arm top diameter"] = 1.0
        options["Arm twist angle degrees"] = 0.0
        options["Wrist thickness"] = 0.5
        options["Wrist width"] = 0.7
        options["Hand length"] = 2.0
        options["Hand thickness"] = 0.2
        options["Hand width"] = 1.0
        options["Thorax length"] = 2.5 
        options["Abdomen length"] = 3.0
        options["Torso depth"] = 2.5
        options["Torso width"] = 3.2
        options["Pelvis drop"] = 1.5
        options["Pelvis width"] = 2.0
        options["Left leg abduction degrees"] = 10.0
        options["Right leg abduction degrees"] = 10.0
        options["Leg length"] = 10.0
        options["Leg top diameter"] = 2.0
        options["Leg bottom diameter"] = 0.7
        options["Left knee flexion degrees"] = 0.0
        options["Right knee flexion degrees"] = 0.0
        options["Left ankle flexion degrees"] = 90.0
        options["Right ankle flexion degrees"] = 90.0
        options["Foot height"] = 1.25
        options["Foot length"] = 2.5
        options["Foot thickness"] = 0.3
        options["Foot width"] = 1.0
        options["Inner proportion default"] = 0.7
        options["Inner proportion head"] = 0.35
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
            "Left shoulder abduction degrees",
            "Right shoulder abduction degrees",
            "Left shoulder flexion degrees",
            "Right shoulder flexion degrees",
            "Left elbow flexion degrees",
            "Right elbow flexion degrees",
            "Arm length",
            "Arm top diameter",
            "Arm twist angle degrees",
            "Wrist thickness",
            "Wrist width",
            "Hand length",
            "Hand thickness",
            "Hand width",
            "Thorax length",
            "Abdomen length",
            "Torso depth",
            "Torso width",
            "Pelvis drop",
            "Pelvis width",
            "Left leg abduction degrees",
            "Right leg abduction degrees",
            "Leg length",
            "Leg top diameter",
            "Leg bottom diameter",
            "Left knee flexion degrees",
            "Right knee flexion degrees",
            "Foot height",
            "Foot length",
            "Foot thickness",
            "Foot width",
            "Left ankle flexion degrees",
            "Right ankle flexion degrees",
            "Inner proportion default",
            "Inner proportion head"
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
            "Thorax length",
            "Abdomen length",
            "Torso depth",
            "Torso width",
            "Leg length",
            "Leg top diameter",
            "Leg bottom diameter",
            "Foot height",
            "Foot length",
            "Foot thickness",
            "Foot width"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
        for key in [
            "Inner proportion default",
            "Inner proportion head"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 0.9:
                options[key] = 0.9
        for key, angleRange in {
            "Left shoulder abduction degrees": (-60.0, 180.0),
            "Right shoulder abduction degrees": (-60.0, 180.0),
            "Left shoulder flexion degrees": (-60.0, 200.0),
            "Right shoulder flexion degrees": (-60.0, 200.0),
            "Left elbow flexion degrees": (0.0, 150.0),
            "Right elbow flexion degrees": (0.0, 150.0),
            "Left knee flexion degrees": (0.0, 140.0),
            "Right knee flexion degrees": (0.0, 140.0),
            "Left ankle flexion degrees": (60.0, 140.0),
            "Right ankle flexion degrees": (60.0, 140.0),
            "Arm twist angle degrees": (-90.0, 90.0),
            "Left leg abduction degrees": (-20.0, 60.0),
            "Right leg abduction degrees": (-20.0, 60.0)
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
        # parameterSetName = options['Base parameter set']
        structure = options["Structure"]
        halfHeadDepth = 0.5 * options["Head depth"]
        headLength = options["Head length"]
        halfHeadWidth = 0.5 * options["Head width"]
        neckLength = options["Neck length"]
        shoulderDrop = options["Shoulder drop"]
        halfShoulderWidth = 0.5 * options["Shoulder width"]
        shoulderLeftFlexionRadians = math.radians(options["Left shoulder flexion degrees"])
        shoulderRightFlexionRadians = math.radians(options["Right shoulder flexion degrees"])
        armLeftAngleRadians = math.radians(options["Left shoulder abduction degrees"])
        armRightAngleRadians = math.radians(options["Right shoulder abduction degrees"])
        elbowLeftFlexionRadians = math.radians(options["Left elbow flexion degrees"])
        elbowRightFlexionRadians = math.radians(options["Right elbow flexion degrees"])
        armLength = options["Arm length"]
        armTopRadius = 0.5 * options["Arm top diameter"]
        armTwistAngleRadians = math.radians(options["Arm twist angle degrees"])
        halfWristThickness = 0.5 * options["Wrist thickness"]
        halfWristWidth = 0.5 * options["Wrist width"]
        handLength = options["Hand length"]
        halfHandThickness = 0.5 * options["Hand thickness"]
        halfHandWidth = 0.5 * options["Hand width"]
        halfTorsoDepth = 0.5 * options["Torso depth"]
        thoraxLength = options["Thorax length"]
        abdomenLength = options["Abdomen length"]
        halfTorsoWidth = 0.5 * options["Torso width"]
        pelvisDrop = options["Pelvis drop"]
        halfPelvisWidth = 0.5 * options["Pelvis width"]
        leftLegAbductionRadians = math.radians(options["Left leg abduction degrees"])
        rightLegAbductionRadians = math.radians(options["Right leg abduction degrees"])
        legLength = options["Leg length"]
        legTopRadius = 0.5 * options["Leg top diameter"]
        legBottomRadius = 0.5 * options["Leg bottom diameter"]
        kneeLeftFlexionRadians = math.radians(options["Left knee flexion degrees"])
        kneeRightFlexionRadians = math.radians(options["Right knee flexion degrees"])
        ankleLeftFlexionRadians = math.radians( options["Left ankle flexion degrees"])
        ankleRightFlexionRadians = math.radians(options["Right ankle flexion degrees"])
        footHeight = options["Foot height"]
        footLength = options["Foot length"]
        halfFootThickness = 0.5 * options["Foot thickness"]
        halfFootWidth = 0.5 * options["Foot width"]
        innerProportionDefault = options["Inner proportion default"]
        innerProportionHead = options["Inner proportion head"]
        # Store coordinates for kinematic tree markers
        options['Kinematic tree'] = {}
        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)
        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)
        # set up element annotations
        bodyGroup = AnnotationGroup(region, get_body_term("body"))
        headGroup = AnnotationGroup(region, get_body_term("head"))
        neckGroup = AnnotationGroup(region, get_body_term("neck"))
        armGroup = AnnotationGroup(region, get_body_term("upper limb"))
        # armToHandGroup = AnnotationGroup(region, ("arm to hand", ""))
        leftArmGroup = AnnotationGroup(region, get_body_term("left upper limb"))
        leftBrachiumGroup = AnnotationGroup(region, get_body_term("left brachium"))
        leftAntebrachiumGroup = AnnotationGroup(region, get_body_term("left antebrachium"))
        # leftElbowGroup = AnnotationGroup(region, get_body_term("left elbow"))
        leftHandGroup = AnnotationGroup(region, get_body_term("left hand"))
        rightArmGroup = AnnotationGroup(region, get_body_term("right upper limb"))
        rightBrachiumGroup = AnnotationGroup(region, get_body_term("right brachium"))
        rightAntebrachiumGroup = AnnotationGroup(region, get_body_term("right antebrachium"))
        # rightElbowGroup = AnnotationGroup(region, get_body_term("right elbow"))
        rightHandGroup = AnnotationGroup(region, get_body_term("right hand"))
        handGroup = AnnotationGroup(region, get_body_term("hand"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        legGroup = AnnotationGroup(region, get_body_term("lower limb"))
        legToFootGroup = AnnotationGroup(region, ("leg to foot", ""))
        leftLegGroup = AnnotationGroup(region, get_body_term("left lower limb"))
        leftUpperLegGroup = AnnotationGroup(region, get_body_term("left upper leg"))
        leftLowerLegGroup = AnnotationGroup(region, get_body_term("left lower leg"))
        leftFootGroup = AnnotationGroup(region, get_body_term("left foot"))
        rightLegGroup = AnnotationGroup(region, get_body_term("right lower limb"))
        rightUpperLegGroup = AnnotationGroup(region, get_body_term("right upper leg"))
        rightLowerLegGroup = AnnotationGroup(region, get_body_term("right lower leg"))
        rightFootGroup = AnnotationGroup(region, get_body_term("right foot"))
        footGroup = AnnotationGroup(region, get_body_term("foot"))
        annotationGroups = [bodyGroup, headGroup, neckGroup,
                            thoraxGroup, abdomenGroup,
                            leftBrachiumGroup, leftAntebrachiumGroup, leftHandGroup, 
                            rightBrachiumGroup, rightAntebrachiumGroup, rightHandGroup,
                            leftLegGroup, leftUpperLegGroup, leftLowerLegGroup, leftFootGroup,
                            rightLegGroup, rightUpperLegGroup, rightLowerLegGroup, rightFootGroup,
                            armGroup, leftArmGroup, rightArmGroup, handGroup,
                            legGroup, footGroup]
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        elementIdentifier = 1
        headElementsCount = humanElementCounts['headElementsCount']
        meshGroups = [bodyMeshGroup, headGroup.getMeshGroup(mesh)]
        for e in range(headElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        neckElementsCount = humanElementCounts['neckElementsCount']
        meshGroups = [bodyMeshGroup, neckGroup.getMeshGroup(mesh)]
        for e in range(neckElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        left = 0
        right = 1
        brachiumElementsCount = humanElementCounts['brachiumElementsCount']
        antebrachiumElementsCount = humanElementCounts['antebrachiumElementsCount']
        handElementsCount = humanElementCounts['handElementsCount']
        armToHandElementsCount = brachiumElementsCount + antebrachiumElementsCount #all elbow nodes count as 1
        armMeshGroup = armGroup.getMeshGroup(mesh)
        # armToHandMeshGroup = armToHandGroup.getMeshGroup(mesh)
        handMeshGroup = handGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideArmGroup = leftArmGroup if (side == left) else rightArmGroup
            sideBrachiumGroup = leftBrachiumGroup if (side == left) else rightBrachiumGroup
            sideAntebrachiumGroup = leftAntebrachiumGroup if (side == left) else rightAntebrachiumGroup
            # sideElbowGroup = leftElbowGroup if (side == left) else rightElbowGroup
            sideHandGroup = leftHandGroup if (side == left) else rightHandGroup
            # Setup brachium elements
            meshGroups = [bodyMeshGroup, 
                          armMeshGroup, 
                          sideArmGroup.getMeshGroup(mesh), sideBrachiumGroup.getMeshGroup(mesh)]
            for e in range(brachiumElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
            # Setup antebrachium elements
            meshGroups = [bodyMeshGroup, 
                          armMeshGroup,
                           sideArmGroup.getMeshGroup(mesh), sideAntebrachiumGroup.getMeshGroup(mesh)]
            for e in range(antebrachiumElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
            # Setup hand elements
            meshGroups = [bodyMeshGroup, 
                          armMeshGroup, 
                          handMeshGroup, sideHandGroup.getMeshGroup(mesh)]
            for e in range(handElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
        thoraxElementsCount = humanElementCounts['thoraxElementsCount']
        abdomenElementsCount = humanElementCounts['abdomenElementsCount']
        # Setup thorax elements
        meshGroups = [bodyMeshGroup, thoraxGroup.getMeshGroup(mesh)]
        for e in range(thoraxElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        # Setup abdomen elements 
        meshGroups = [bodyMeshGroup, abdomenGroup.getMeshGroup(mesh)]
        for e in range(abdomenElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        upperLegElementsCount = humanElementCounts['upperLegElementsCount']
        lowerLegElementsCount = humanElementCounts['lowerLegElementsCount']
        footElementsCount = humanElementCounts['footElementsCount']
        legToFootElementsCount = upperLegElementsCount + lowerLegElementsCount
        legMeshGroup = legGroup.getMeshGroup(mesh)
        legToFootMeshGroup = legToFootGroup.getMeshGroup(mesh)
        footMeshGroup = footGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideLegGroup = leftLegGroup if (side == left) else rightLegGroup
            sideUpperLegGroup = leftUpperLegGroup if (side == left) else rightUpperLegGroup
            sideLowerLegGroup = leftLowerLegGroup if (side == left) else rightLowerLegGroup
            sideFootGroup = leftFootGroup if (side == left) else rightFootGroup
            # Upper leg
            meshGroups = [bodyMeshGroup, legMeshGroup, 
                          sideLegGroup.getMeshGroup(mesh), sideUpperLegGroup.getMeshGroup(mesh)]
            for e in range(upperLegElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
            # Lower leg
            meshGroups = [bodyMeshGroup, legMeshGroup, 
                          sideLegGroup.getMeshGroup(mesh), sideLowerLegGroup.getMeshGroup(mesh)]
            for e in range(lowerLegElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
            # Foot
            meshGroups = [bodyMeshGroup, legMeshGroup, footMeshGroup, sideFootGroup.getMeshGroup(mesh)]
            for e in range(footElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        # need to ensure inner coordinates are at least defined:
        cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=0.75)
        innerCoordinates = find_or_create_field_coordinates(fieldmodule, "inner coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        headScale = headLength / headElementsCount
        nodeIdentifier = 1
        d1 = [headScale, 0.0, 0.0]
        d2 = [0.0, halfHeadWidth, 0.0]
        d3 = [0.0, 0.0, halfHeadDepth]
        id2 = mult(d2, innerProportionHead)
        id3 = mult(d3, innerProportionHead)
        for i in range(headElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [headScale * i, 0.0, 0.0]
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3)
            nodeIdentifier += 1

        neckScale = neckLength / neckElementsCount
        d2 = [0.0, halfHeadWidth, 0.0]
        d3 = [0.0, 0.0, halfHeadWidth]
        id2 = mult(d2, innerProportionHead)
        id3 = mult(d3, innerProportionHead)
        for i in range(neckElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [headLength + neckScale * i, 0.0, 0.0]
            d1 = [0.5 * (headScale + neckScale) if (i == 0) else neckScale, 0.0, 0.0]
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3)
            nodeIdentifier += 1
        armJunctionNodeIdentifier = nodeIdentifier
        options['Kinematic tree']['thorax_top'] = x

        thoraxScale = thoraxLength / thoraxElementsCount
        thoraxStartX = headLength + neckLength
        sx = [thoraxStartX, 0.0, 0.0]
        for i in range(thoraxElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [thoraxStartX + thoraxScale * i, 0.0, 0.0]
            if i == 0:
                d1 = [0.5 * (neckScale + thoraxScale), 0.0, 0.0]
                d2 = [0.0, 0.5 * (halfTorsoWidth + halfHeadWidth), 0.0]
                d12 = [0.0, halfTorsoWidth - halfHeadWidth, 0.0]
                d3 = [0.0, 0.0, 0.5 * (halfHeadWidth + halfTorsoDepth)]
                id2 = [0.0, 0.5 * (innerProportionHead * halfHeadWidth + innerProportionDefault * halfTorsoWidth), 0.0]
                id12 = [0.0, innerProportionDefault * halfTorsoWidth - innerProportionHead * halfHeadWidth, 0.0]
                id3 = mult(d3, 0.5 * (innerProportionHead + innerProportionDefault))
            else:
                d1 = [thoraxScale, 0.0, 0.0]
                d2 = [0.0, halfTorsoWidth, 0.0]
                d12 = None
                d3 = [0.0, 0.0, halfTorsoDepth]
                id2 = mult(d2, innerProportionDefault)
                id12 = None
                id3 = mult(d3, innerProportionDefault)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
            nodeIdentifier += 1
        
        abdomenScale = abdomenLength / abdomenElementsCount
        d2 = [0.0, halfTorsoWidth, 0.0]
        d3 = [0.0, 0.0, halfTorsoDepth]
        id2 = mult(d2, innerProportionDefault)
        id3 = mult(d3, innerProportionDefault)
        abdomenStartX = thoraxStartX + thoraxLength
        for i in range(abdomenElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [abdomenStartX + abdomenScale * i, 0.0, 0.0]
            d1 = [0.5 * (thoraxScale + abdomenScale) if (i == 0) else abdomenScale, 0.0, 0.0]
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3)
            nodeIdentifier += 1
        legJunctionNodeIdentifier = nodeIdentifier - 1
        px = [abdomenStartX + abdomenLength, 0.0, 0.0]
        options['Kinematic tree']['lumbar_body'] = px
        # arms
        for side in (left, right):
            # Shoulder rotation
            # rotate shoulder with arm, pivoting about shoulder drop below arm junction on network
            # this has the realistic effect of shoulders becoming narrower with higher angles
            # initial shoulder rotation with arm is negligible, hence:
            armAngleRadians = armLeftAngleRadians if (side == left) else armRightAngleRadians
            shoulderRotationFactor = 1.0 - math.cos(0.5 * armAngleRadians)
            # assume shoulder drop is half shrug distance to get limiting shoulder angle for 180 degree arm rotation
            shoulderLimitAngleRadians = math.asin(1.5 * shoulderDrop / halfShoulderWidth)
            shoulderAngleRadians = shoulderRotationFactor * shoulderLimitAngleRadians
            nonHandArmLength = armLength - handLength
            shoulderElementCount = 2
            armScale = nonHandArmLength / (armToHandElementsCount - shoulderElementCount)
            d12_mag = (halfWristThickness - armTopRadius) / (armToHandElementsCount - shoulderElementCount)
            d13_mag = (halfWristWidth - armTopRadius) / (armToHandElementsCount - shoulderElementCount)
            armAngle = armAngleRadians if (side == left) else -armAngleRadians
            cosArmAngle = math.cos(armAngle)
            sinArmAngle = math.sin(armAngle)
            armStartX = thoraxStartX + shoulderDrop - halfShoulderWidth * math.sin(shoulderAngleRadians)
            armStartY = (halfShoulderWidth if (side == left) else -halfShoulderWidth) * math.cos(shoulderAngleRadians)
            armStart = [armStartX, armStartY, 0.0]
            x = armStart
            armDirn = [cosArmAngle, sinArmAngle, 0.0]
            armSide = [-sinArmAngle, cosArmAngle, 0.0]
            armFront = cross(armDirn, armSide)
            d1 = mult(armDirn, armScale)
            # set arm versions 2 (left) and 3 (right) on arm junction node, and intermediate shoulder node
            sd1 = interpolateLagrangeHermiteDerivative(sx, x, d1, 0.0)
            nx, nd1 = sampleCubicHermiteCurvesSmooth([sx, x], [sd1, d1], 2, derivativeMagnitudeEnd=armScale)[0:2]
            arcLengths = [getCubicHermiteArcLength(nx[i], nd1[i], nx[i + 1], nd1[i + 1]) for i in range(2)]
            sd2_list = []
            sd3_list = []
            sNodeIdentifiers = []
            side_label = 'l' if (side == left) else 'r'
            options['Kinematic tree']['humerus_' + side_label] = nx[1]
            for i in range(2):
                sNodeIdentifiers.append(nodeIdentifier if (i > 0) else armJunctionNodeIdentifier)
                node = nodes.findNodeByIdentifier(sNodeIdentifiers[-1])
                fieldcache.setNode(node)
                version = 1 if (i > 0) else 2 if (side == left) else 3
                sd1 = nd1[i]
                sDistance = sum(arcLengths[i:])
                sHalfHeight = armTopRadius + sDistance * -d12_mag
                sHalfDepth = armTopRadius + sDistance * -d13_mag
                sd3 = [0.0, 0.0, sHalfDepth]
                sid3 = mult(sd3, innerProportionDefault)
                sd2 = set_magnitude(cross(sd3, sd1), sHalfHeight)
                sid2 = mult(sd2, innerProportionDefault)
                sd2_list.append(sd2)
                sd3_list.append(sd3)
                if i > 0:
                    for field in (coordinates, innerCoordinates):
                        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, nx[i])
                    nodeIdentifier += 1
                setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, sd3)
                setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, sid3)
            sd2_list.append([-armTopRadius * sinArmAngle, armTopRadius * cosArmAngle, 0.0])
            sd3_list.append([0.0, 0.0, armTopRadius])
            for i in range(2):
                node = nodes.findNodeByIdentifier(sNodeIdentifiers[i])
                fieldcache.setNode(node)
                version = 1 if (i > 0) else 2 if (side == left) else 3
                sd12 = sub(sd2_list[i + 1], sd2_list[i])
                sd13 = sub(sd3_list[i + 1], sd3_list[i])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sd12)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sd13)
                sid12 = mult(sd12, innerProportionDefault)
                sid13 = mult(sd13, innerProportionDefault)
                innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sid12)
                innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sid13)
            # Arm twist
            elementTwistAngle = ((armTwistAngleRadians if (side == left) else -armTwistAngleRadians) /
                                 (armToHandElementsCount - 3))
            # Shoulder flexion
            d2 = mult(armSide, armScale)
            d3 = mult(armFront, armScale)
            # Updating frame of reference wrt flexion angle (using d2 as rotation axis)
            shoulderFlexionRadians = shoulderLeftFlexionRadians if (side == left) else shoulderRightFlexionRadians
            shoulderJointAngleRadians = math.pi - shoulderFlexionRadians
            shoulderRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), shoulderFlexionRadians)
            shoulderHalfRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), shoulderFlexionRadians/2)
            armDirn = matrix_vector_mult(shoulderRotationMatrix, armDirn)
            armSide = armSide
            armFront = cross(armDirn, armSide)
            # The d3 direction in the shoulder node is rotated by half this angle 
            # To ensure a better transition at this node. 
            shoulderDirn = armDirn
            shoulderSide = armSide
            shoulderFront = matrix_vector_mult(shoulderRotationMatrix, d3)
            # This rotation factor is used to adjust the position of the knee node relative 
            # to the angle of flexion, and ensures a proper transition between the upper and lower leg
            rotationFactor = 1.0*math.sin(shoulderFlexionRadians)*(math.sqrt(2)-1)      
            i = 0 
            xi = i / (armToHandElementsCount - 2)
            halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
            halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
            # halfWidth = halfWidth/math.sin(shoulderJointAngleRadians/2)
            shoulderPosition = add(armStart, set_magnitude(d1, 0))
            x = shoulderPosition
            d1 = set_magnitude(shoulderDirn, armScale)
            d2 = set_magnitude(shoulderSide, halfThickness)
            d3 = set_magnitude(shoulderFront, halfWidth)
            d12 = set_magnitude(shoulderSide, d12_mag)
            d13 = set_magnitude(shoulderFront, d13_mag)
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1
            # armStart = add(shoulderPosition,d1)
            # d1 = mult(armDirn, armScale)
            # Setting brachium coordinates
            for i in range(1, brachiumElementsCount - shoulderElementCount):
                xi = i / (armToHandElementsCount - 2)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(armStart, mult(d1, i))
                halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
                halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
                if i == 0:
                    twistAngle = 0.0
                else:
                    twistAngle = -0.5 * elementTwistAngle + elementTwistAngle * i
                if twistAngle == 0.0:
                    d2 = mult(armSide, halfThickness)
                    d3 = mult(armFront, halfWidth)
                    d12 = mult(armSide, d12_mag)
                    d13 = mult(armFront, d13_mag)
                else:
                    cosTwistAngle = math.cos(twistAngle)
                    sinTwistAngle = math.sin(twistAngle)
                    d2 = sub(mult(armSide, halfThickness * cosTwistAngle),
                             mult(armFront, halfThickness * sinTwistAngle))
                    d3 = add(mult(armFront, halfWidth * cosTwistAngle),
                             mult(armSide, halfWidth * sinTwistAngle))
                    d12 = set_magnitude(d2, d12_mag)
                    d13 = set_magnitude(d3, d13_mag)
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id12 = mult(d12, innerProportionDefault)
                id13 = mult(d13, innerProportionDefault)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
            # Diagram to calculate elbow flexion
            # 1 
            # |
            # 2 -- 3 
            # 1 is brachium, 2 is the elbow, 4 is antebrachium
            # we fix the position of the nodes 1 and 3, and calculate 
            # the position of 2 using the sampleCubicHermiteCurvesSmooth function. 
            # This process also gives us the correct d1 and d3 directions for nodes 1 to 5.

            # Calculating initial d2 and d3 before rotation
            # Necessary in case there is a non-zero twist angle
            if twistAngle == 0.0:
                d2 = armSide 
                d3 = armFront
            else:
                cosTwistAngle = math.cos(twistAngle)
                sinTwistAngle = math.sin(twistAngle)
                d2 = sub(mult(armSide, cosTwistAngle),
                            mult(armFront, sinTwistAngle))
                d3 = add(mult(armFront,  cosTwistAngle),
                            mult(armSide, sinTwistAngle))
            # Updating frame of reference wrt flexion angle (using d2 as rotation axis)
            elbowFlexionRadians = elbowLeftFlexionRadians if (side == left) else elbowRightFlexionRadians
            elbowJointAngleRadians = math.pi - elbowFlexionRadians
            elbowRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), elbowFlexionRadians)
            elbowHalfRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), elbowFlexionRadians/2)
            # Calculating initial estimation for directions at the elbow node
            elbowDirn = matrix_vector_mult(elbowHalfRotationMatrix, d1)
            elbowFront = matrix_vector_mult(elbowHalfRotationMatrix, d3)
            elbowSide = d2
            # Calculating direction for the antebrachium
            antebrachiumDirn = matrix_vector_mult(elbowRotationMatrix, armDirn)
            antebrachiumSide = armSide
            antebrachiumFront = cross(antebrachiumDirn, antebrachiumSide)
            # This rotation factor is used to adjust the position of the joint node relative 
            # to the angle of flexion, and ensures a proper transition between the two parts
            rotationFactor = math.sin(elbowFlexionRadians)*(math.sqrt(2)-1)     
            jointPositions = []
            # jointPositions.append(sub(x, d1)) #1 
            jointPositions.append(x)
            exi = (i+1) / (armToHandElementsCount - 2)
            ehalfdWidth = exi * halfWristWidth + (1.0 - exi) * armTopRadius
            ehalfdWidth = ehalfdWidth * rotationFactor
            eDir = add(set_magnitude(armDirn, armScale), set_magnitude(elbowFront, ehalfdWidth))
            eDir = set_magnitude(eDir, armScale)
            jointPositions.append(add(x, eDir)) # 2
            eDir =  add(set_magnitude(elbowFront, -ehalfdWidth), set_magnitude(antebrachiumDirn, armScale))
            eDir = set_magnitude(eDir, armScale)
            jointPositions.append(add(jointPositions[-1], eDir)) #3
            jointDir = [armDirn, elbowDirn, antebrachiumDirn] 
            jointPositions, jointDirn = sampleCubicHermiteCurvesSmooth(
                jointPositions, jointDir, 2, derivativeMagnitudeStart=armScale, derivativeMagnitudeEnd=armScale)[0:2]  
            j = 1
            i += 1
            xi = i / (armToHandElementsCount - 2)
            halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
            # if (j == 1): 
            halfWidth = (halfWidth)*(1/math.sin(elbowJointAngleRadians/2) - 0.5*rotationFactor)
            halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
            x = jointPositions[j]
            d1 = jointDirn[j]
            elbowFront = cross(d1, d2)
            d2 = set_magnitude(elbowSide, halfThickness)
            d3 = set_magnitude(elbowFront, halfWidth)
            d12 = set_magnitude(elbowSide, d12_mag)
            d13 = set_magnitude(elbowFront, d13_mag)
            d13 = add(d13, set_magnitude(d1, -2*halfWidth*math.sin(elbowJointAngleRadians/2)))
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1
            options['Kinematic tree']['ulna_' + side_label] = x
            # Antebrachium nodes starts after the elbow node
            antebrachiumStart = jointPositions[-1]
            d1 = mult(antebrachiumDirn, armScale)
            # Change d1 to the antebrachium direction
            for i in range(brachiumElementsCount - shoulderElementCount + 1, armToHandElementsCount - 1):
                xi = (i) / (armToHandElementsCount - 2)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(antebrachiumStart, mult(d1, i - (brachiumElementsCount - shoulderElementCount + 1))) 
                halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
                halfWidth =  xi * halfWristWidth + (1.0 - xi) * armTopRadius
                if i == 0:
                    twistAngle = armTwistAngleRadians if (side == left) else -armTwistAngleRadians
                else:
                    twistAngle = -0.5 * elementTwistAngle + elementTwistAngle * (i - 1)
                if twistAngle == 0.0:
                    d2 = mult(antebrachiumSide, halfThickness)
                    d3 = mult(antebrachiumFront, halfWidth)
                    d12 = mult(antebrachiumSide, d12_mag)
                    d13 = mult(antebrachiumFront, d13_mag)
                else:
                    cosTwistAngle = math.cos(twistAngle)
                    sinTwistAngle = math.sin(twistAngle)
                    d2 = sub(mult(antebrachiumSide, halfThickness * cosTwistAngle),
                             mult(antebrachiumFront, halfThickness * sinTwistAngle))
                    d3 = add(mult(antebrachiumFront, halfWidth * cosTwistAngle),
                             mult(antebrachiumSide, halfWidth * sinTwistAngle))
                    d12 = set_magnitude(d2, d12_mag)
                    d13 = set_magnitude(d3, d13_mag)
                    if i < (antebrachiumElementsCount - 1):
                        d12 = add(d12, set_magnitude(d3, -halfThickness * elementTwistAngle))
                        d13 = add(d13, set_magnitude(d2, halfWidth * elementTwistAngle))
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id12 = mult(d12, innerProportionDefault)
                id13 = mult(d13, innerProportionDefault)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
            # Hand twist
            options['Kinematic tree']['wrist_' + side_label] = x
            assert handElementsCount == 1
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            hx = add(x, mult(antebrachiumDirn, handLength))
            hd1 = computeCubicHermiteEndDerivative(x, d1, hx, d1)
            twistAngle = armTwistAngleRadians if (side == left) else -armTwistAngleRadians
            if twistAngle == 0.0:
                hd2 = set_magnitude(d2, halfHandThickness)
                hd3 = set_magnitude(d3, halfHandWidth)
            else:
                cosTwistAngle = math.cos(twistAngle)
                sinTwistAngle = math.sin(twistAngle)
                hd2 = sub(mult(antebrachiumSide, halfHandThickness * cosTwistAngle),
                          mult(antebrachiumFront, halfHandThickness * sinTwistAngle))
                hd3 = add(mult(antebrachiumFront, halfHandWidth * cosTwistAngle),
                          mult(antebrachiumSide, halfHandWidth * sinTwistAngle))
            hid2 = mult(hd2, innerProportionDefault)
            hid3 = mult(hd3, innerProportionDefault)
            setNodeFieldParameters(coordinates, fieldcache, hx, hd1, hd2, hd3)
            setNodeFieldParameters(innerCoordinates, fieldcache, hx, hd1, hid2, hid3)
            nodeIdentifier += 1
        # legs
        legStartX = abdomenStartX + abdomenLength + pelvisDrop
        nonFootLegLength = legLength - footHeight
        legScale = nonFootLegLength / (legToFootElementsCount - 1)
        d12_mag = (legBottomRadius - legTopRadius) / (armToHandElementsCount - 2)
        d13_mag = (legBottomRadius - legTopRadius) / (armToHandElementsCount - 2)
        pd3 = [0.0, 0.0, 0.5 * legTopRadius + 0.5 * halfTorsoDepth]
        pid3 = mult(pd3, innerProportionDefault)
        for side in (left, right):
            legAngle = leftLegAbductionRadians if (side == left) else -rightLegAbductionRadians
            cosLegAngle = math.cos(legAngle)
            sinLegAngle = math.sin(legAngle)
            legStartY = halfPelvisWidth if (side == left) else -halfPelvisWidth
            x = legStart = [legStartX, legStartY, 0.0]
            legDirn = [cosLegAngle, sinLegAngle, 0.0]
            legSide = [-sinLegAngle, cosLegAngle, 0.0]
            legFront = cross(legDirn, legSide)
            d1 = mult(legDirn, legScale)
            # set leg versions 2 (left) and 3 (right) on leg junction node
            node = nodes.findNodeByIdentifier(legJunctionNodeIdentifier)
            fieldcache.setNode(node)
            pd1 = interpolateLagrangeHermiteDerivative(px, x, d1, 0.0)
            pd2 = set_magnitude(cross(pd3, pd1), 0.5 * legTopRadius + 0.5 * halfTorsoWidth)
            pid2 = mult(pd2, innerProportionDefault)
            pd12 = sub(mult(legSide, legTopRadius), pd2)
            pd13 = sub([0.0, 0.0, legTopRadius], pd3)
            pid12 = mult(pd12, innerProportionDefault)
            pid13 = mult(pd13, innerProportionDefault)
            version = 2 if (side == left) else 3
            setNodeFieldVersionDerivatives(coordinates, fieldcache, version, pd1, pd2, pd3, pd12, pd13)
            setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, pd1, pid2, pid3, pid12, pid13)
            d12 = [-d12_mag * sinLegAngle, d12_mag * cosLegAngle, 0.0]
            id12 = mult(d12, innerProportionDefault)
            d13 = [0.0, 0.0, d13_mag]
            id13 = mult(d13, innerProportionDefault)
            # Upper leg
            for i in range(upperLegElementsCount):
                xi = i / legToFootElementsCount
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(legStart, mult(d1, i))
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = mult(legSide, radius)
                d3 = [0.0, 0.0, radius]
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
            # knee
            # Updating frame of reference wrt rotation angle (using d2 as rotation axis)
            kneeFlexionRadians = kneeLeftFlexionRadians if (side == left) else kneeRightFlexionRadians
            # Angle between upper and lower leg
            kneeJointAngleRadians = math.pi - kneeFlexionRadians
            kneeRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, 1), kneeFlexionRadians)
            kneeHalfrotationMatrix = axis_angle_to_rotation_matrix(mult(d2, 1), kneeFlexionRadians/2)
            lowerLegDirn = matrix_vector_mult(kneeRotationMatrix, d1)
            lowerLegSide = legSide
            lowerLegFront = cross(lowerLegDirn, lowerLegSide)
            # The d3 direction in the knee node is rotated by half the flexion angle
            # To ensure a better transition at this node. 
            kneeDirn = lowerLegDirn
            kneeSide = d2 
            kneeFront = matrix_vector_mult(kneeHalfrotationMatrix, d3)
            # This rotation factor is used to adjust the position of the knee node relative 
            # to the angle of flexion, and ensures a proper transition between the upper and lower leg
            rotationFactor = 2.0*math.sin(kneeFlexionRadians)*(math.sqrt(2)-1)
            i += 1
            xi = i / legToFootElementsCount
            radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
            # d3 direction is fattened at the joint to ensure a proper transition in the tube network
            kneeRadius = radius/math.sin(kneeJointAngleRadians/2)
            kneePosition = add(x, set_magnitude(d1, legScale - rotationFactor*kneeRadius))
            x = kneePosition
            d1 = set_magnitude(kneeDirn, legScale + rotationFactor*kneeRadius)
            d2 = set_magnitude(kneeSide, radius)
            d3 = set_magnitude(kneeFront, kneeRadius)
            d12 = set_magnitude(kneeSide, d12_mag)
            d13 = set_magnitude(kneeFront, d13_mag)
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1
            # Lower leg
            lowerLegStart = add(kneePosition, d1)
            d1 = set_magnitude(lowerLegDirn, legScale)
            for i in range(upperLegElementsCount + 1, legToFootElementsCount):
                xi = i / legToFootElementsCount
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(lowerLegStart, mult(d1, i - (upperLegElementsCount + 1)))
                # if (i == legToFootElementsCount): 
                #     # x = add(x, set_magnitude(d1, 1.5 * radius))
                #     d1 = set_magnitude(d1, legScale + 1.5 * radius)
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = set_magnitude(lowerLegSide, radius)
                d3 = set_magnitude(lowerLegFront, radius)
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
            # foot
            # Updating frame of reference wrt rotation angle (using d2 as rotation axis)
            ankleFlexionRadians = ankleLeftFlexionRadians if (side == left) else ankleRightFlexionRadians
            # Angle between lower leg and foot
            ankleJointAngleRadians = math.pi - ankleFlexionRadians
            ankleRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), (ankleFlexionRadians))
            ankleHalfRotationMatrix = axis_angle_to_rotation_matrix(mult(d2, -1), (ankleFlexionRadians/2))
            # The d3 direction in the ankle node is rotated by half this angle 
            # To ensure a better transition at this node.
            ankleDirn = matrix_vector_mult(ankleRotationMatrix, d1)
            ankleSide = d2
            ankleFront = cross(ankleDirn, ankleSide)
            footDirn = ankleDirn
            footSide = d2 
            footFront = matrix_vector_mult(ankleHalfRotationMatrix, d3)
            footd1 = footDirn
            footd2 = set_magnitude(footSide, halfFootWidth)
            footd3 = footFront
            # This rotation factor is used to adjust 'fatten' the d3 direction at the ankle node
            # As well as adjusting the position of the node depending on the 
            # angle of flexion, and ensures a proper transition between lower leg and the foot
            rotationFactor = 2.0*math.sin(ankleFlexionRadians)*(math.sqrt(2)- 1)
            ankleThickness = halfFootThickness/math.sin(ankleJointAngleRadians/2)
            # This positioning of the food nodes bends edge connecting the leg and the foot 
            # Which allows the scaffold to better capture the shape of the calcaneus
            anklePosition = add(x, set_magnitude(d1, legScale - 1*rotationFactor*footHeight))
            fx = [
                x, 
                add(anklePosition, set_magnitude(footd1, 0)), 
                add(anklePosition, set_magnitude(footd1, footLength - legBottomRadius)), 

            ]
            fd1 = [d1, set_magnitude(footd1, 0.5*footLength), set_magnitude(footd1, 0.5*footLength)]
            fd1 = smoothCubicHermiteDerivativesLine(
                fx, fd1, fixAllDirections=True, fixStartDerivative=True
                )
            fd2 = [d2, footd2, footd2]
            fd3 = [d3,
                    set_magnitude(footd3, ankleThickness + legBottomRadius),
                   set_magnitude(cross(fd1[2], fd2[2]), halfFootThickness)
            ]
            fd12 = sub(fd2[2], fd2[1])
            fd13 = sub(fd3[2], fd3[1]) 
            fid12 = mult(fd12, innerProportionDefault)
            fid13 = mult(fd13, innerProportionDefault)
            for i in range(1, 3):
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                setNodeFieldParameters(coordinates, fieldcache, fx[i], fd1[i], fd2[i], fd3[i], fd12, fd13)
                fid2 = mult(fd2[i], innerProportionDefault)
                fid3 = mult(fd3[i], innerProportionDefault)
                setNodeFieldParameters(innerCoordinates, fieldcache, fx[i], fd1[i], fid2, fid3, fid12, fid13)
                nodeIdentifier += 1
            side_label = 'l' if (side == left) else 'r'
            options['Kinematic tree']['toes_' + side_label] = fx[i]

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
        options["Number of elements along head"] = 2
        options["Number of elements along neck"] = 1
        options["Number of elements along thorax"] = 2
        options["Number of elements along abdomen"] = 2
        options["Number of elements along brachium"] = 5
        options["Number of elements along antebrachium"] = 3
        options["Number of elements along hand"] = 1
        options["Number of elements along upper leg"] = 3
        options["Number of elements along lower leg"] = 2
        options["Number of elements along foot"] = 2
        options["Number of elements around head"] = 12
        options["Number of elements around torso"] = 12
        options["Number of elements around arm"] = 8
        options["Number of elements around leg"] = 8
        options["Number of elements through shell"] = 1
        options["Show trim surfaces"] = False
        options["Use Core"] = True
        options["Number of elements across core box minor"] = 2
        options["Number of elements across core transition"] = 1
        if "Medium" in useParameterSetName:
            options["Number of elements along head"] = 3
            options["Number of elements along neck"] = 2
            options["Number of elements along thorax"] = 3
            options["Number of elements along abdomen"] = 3
            options["Number of elements along brachium"] = 4
            options["Number of elements along antebrachium"] = 3
            options["Number of elements along hand"] = 1
            options["Number of elements along upper leg"] = 2
            options["Number of elements along lower leg"] = 2
            options["Number of elements along foot"] = 2
            options["Number of elements around head"] = 16
            options["Number of elements around torso"] = 16
            options["Number of elements around leg"] = 12
        elif "Fine" in useParameterSetName:
            options["Number of elements along head"] = 4
            options["Number of elements along neck"] = 2
            options["Number of elements along thorax"] = 4
            options["Number of elements along abdomen"] = 4
            options["Number of elements along brachium"] = 5
            options["Number of elements along antebrachium"] = 4
            options["Number of elements along hand"] = 2
            options["Number of elements along upper leg"] = 3
            options["Number of elements along lower leg"] = 2
            options["Number of elements along foot"] = 3
            options["Number of elements around head"] = 20
            options["Number of elements around torso"] = 20
            options["Number of elements around arm"] = 12
            options["Number of elements around leg"] = 16
            options["Number of elements through shell"] = 2
            options["Number of elements across core box minor"] = 4

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Body network layout",
            "Number of elements along head",
            "Number of elements along neck",
            "Number of elements along thorax",
            "Number of elements along abdomen",
            "Number of elements along brachium",
            "Number of elements along antebrachium",
            "Number of elements along hand",
            "Number of elements along upper leg",
            "Number of elements along lower leg",
            "Number of elements along foot",
            "Number of elements around head",
            "Number of elements around torso",
            "Number of elements around arm",
            "Number of elements around leg",
            "Number of elements through shell",
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
        if (options["Body network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Body network layout")):
            options["Body network layout"] = ScaffoldPackage(MeshType_1d_human_body_network_layout1)
        for key in [
            "Number of elements along head",
            "Number of elements along neck",
            "Number of elements along thorax",
            "Number of elements along abdomen",
            "Number of elements along brachium",
            "Number of elements along antebrachium",
            "Number of elements along hand",
            "Number of elements along upper leg",
            "Number of elements along lower leg",
            "Number of elements along foot"
        ]:
            if options[key] < 1:
                options[key] = 1
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
        networkLayout = options["Body network layout"]
        elementsCountAlongHead = options["Number of elements along head"]
        elementsCountAlongNeck = options["Number of elements along neck"]
        elementsCountAlongThorax = options["Number of elements along thorax"]
        elementsCountAlongAbdomen = options["Number of elements along abdomen"]
        elementsCountAlongBrachium = options["Number of elements along brachium"]
        elementsCountAlongAntebrachium = options["Number of elements along antebrachium"]
        elementsCountAlongHand = options["Number of elements along hand"]
        elementsCountAlongUpperLeg = options["Number of elements along upper leg"]
        elementsCountAlongLowerLeg = options["Number of elements along lower leg"]
        elementsCountAlongFoot = options["Number of elements along foot"]
        elementsCountAroundHead = options["Number of elements around head"]
        elementsCountAroundTorso = options["Number of elements around torso"]
        elementsCountAroundArm = options["Number of elements around arm"]
        elementsCountAroundLeg = options["Number of elements around leg"]
        isCore = options["Use Core"]

        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        annotationAlongCounts = []
        annotationAroundCounts = []
        defaultCoreBoundaryScalingMode = 1
        annotationCoreBoundaryScalingMode = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            alongCount = 0
            aroundCount = 0
            coreBoundaryScalingMode = 0
            name = layoutAnnotationGroup.getName()
            if "head" in name:
                alongCount = elementsCountAlongHead
                aroundCount = elementsCountAroundHead
            elif "neck" in name:
                alongCount = elementsCountAlongNeck
                aroundCount = elementsCountAroundHead
            elif "thorax" in name:
                alongCount = elementsCountAlongThorax
                aroundCount = elementsCountAroundTorso
                coreBoundaryScalingMode = 2
            elif "abdomen" in name:
                alongCount = elementsCountAlongAbdomen
                aroundCount = elementsCountAroundTorso
                coreBoundaryScalingMode = 2
            elif " brachium" in name:
                alongCount = elementsCountAlongBrachium
                aroundCount = elementsCountAroundArm
            elif " antebrachium" in name:
                alongCount = elementsCountAlongAntebrachium
                aroundCount = elementsCountAroundArm
            elif "hand" in name:
                alongCount = elementsCountAlongHand
                aroundCount = elementsCountAroundArm
            elif "upper leg" in name:
                alongCount = elementsCountAlongUpperLeg
                aroundCount = elementsCountAroundLeg
            elif "lower leg" in name:
                alongCount = elementsCountAlongLowerLeg
                aroundCount = elementsCountAroundLeg
            elif "foot" in name:
                alongCount = elementsCountAlongFoot
                aroundCount = elementsCountAroundLeg
            annotationAlongCounts.append(alongCount)
            annotationAroundCounts.append(aroundCount)
            annotationCoreBoundaryScalingMode.append(coreBoundaryScalingMode)

        tubeNetworkMeshBuilder = BodyTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=2.0,  # not used for body
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAlong=annotationAlongCounts,
            defaultElementsCountAround=options["Number of elements around head"],
            annotationElementsCountsAround=annotationAroundCounts,
            elementsCountThroughShell=options["Number of elements through shell"],
            isCore=isCore,
            elementsCountTransition=options['Number of elements across core transition'],
            defaultElementsCountCoreBoxMinor=options["Number of elements across core box minor"],
            annotationElementsCountsCoreBoxMinor=[],
            defaultCoreBoundaryScalingMode=defaultCoreBoundaryScalingMode,
            annotationCoreBoundaryScalingMode=annotationCoreBoundaryScalingMode,
            useOuterTrimSurfaces=True)

        meshDimension = 3
        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, meshDimension,
            isLinearThroughShell=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        if isCore:
            fieldmodule = region.getFieldmodule()
            mesh = fieldmodule.findMeshByDimension(meshDimension)
            thoraxGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thorax"))
            abdomenGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdomen"))
            coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))

            thoracicCavityGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_body_term("thoracic cavity"))
            is_thoracic_cavity = fieldmodule.createFieldAnd(thoraxGroup.getGroup(), coreGroup.getGroup())
            thoracicCavityGroup.getMeshGroup(mesh).addElementsConditional(is_thoracic_cavity)

            abdominalCavityGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_body_term("abdominal cavity"))
            is_abdominal_cavity = fieldmodule.createFieldAnd(abdomenGroup.getGroup(), coreGroup.getGroup())
            abdominalCavityGroup.getMeshGroup(mesh).addElementsConditional(is_abdominal_cavity)

        # Kinematic tree markers 
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        node_identifier = max(1, get_maximum_node_identifier(nodes) + 1)
        coordinates = find_or_create_field_coordinates(fieldmodule)
        stickman_markers = networkLayout._scaffoldSettings['Kinematic tree']
        for marker_name, marker_position in stickman_markers.items():
            marker_group = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, (marker_name, ""), isMarker=True
                )
            marker_group.createMarkerNode(
                node_identifier, coordinates, marker_position
                )
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
        isCore = options["Use Core"]

        # create 2-D surface mesh groups, 1-D spinal cord
        fieldmodule = region.getFieldmodule()
        mesh2d = fieldmodule.findMeshByDimension(2)
        mesh1d = fieldmodule.findMeshByDimension(1)

        is_exterior = fieldmodule.createFieldIsExterior()
        is_face_xi3_0 = fieldmodule.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)

        skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis outer surface"))
        is_skin = is_exterior if isCore else fieldmodule.createFieldAnd(
            is_exterior, fieldmodule.createFieldNot(is_face_xi3_0))
        skinGroup.getMeshGroup(mesh2d).addElementsConditional(is_skin)

        leftArmGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left upper limb"))
        leftArmSkinGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_body_term("left upper limb skin epidermis outer surface"))
        leftArmSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
            fieldmodule.createFieldAnd(leftArmGroup.getGroup(), is_exterior))
        rightArmGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right upper limb"))
        rightArmSkinGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_body_term("right upper limb skin epidermis outer surface"))
        rightArmSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
            fieldmodule.createFieldAnd(rightArmGroup.getGroup(), is_exterior))
        leftLegGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left lower limb"))
        leftLegSkinGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_body_term("left lower limb skin epidermis outer surface"))
        leftLegSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
            fieldmodule.createFieldAnd(leftLegGroup.getGroup(), is_exterior))
        rightLegGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right lower limb"))
        rightLegSkinGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_body_term("right lower limb skin epidermis outer surface"))
        rightLegSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
            fieldmodule.createFieldAnd(rightLegGroup.getGroup(), is_exterior))

        if isCore:
            coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))
            shellGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("shell"))
            leftGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left"))
            rightGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right"))
            dorsalGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("dorsal"))

            is_core_shell = fieldmodule.createFieldAnd(coreGroup.getGroup(), shellGroup.getGroup())
            is_left_right = fieldmodule.createFieldAnd(leftGroup.getGroup(), rightGroup.getGroup())
            is_left_right_dorsal = fieldmodule.createFieldAnd(is_left_right, dorsalGroup.getGroup())

            neckGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("neck"))
            thoracicCavityGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thoracic cavity"))
            abdominalCavityGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdominal cavity"))
            armGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("upper limb"))
            legGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("lower limb"))

            thoracicCavityBoundaryGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_body_term("thoracic cavity boundary surface"))
            is_thoracic_cavity_boundary = fieldmodule.createFieldAnd(
                thoracicCavityGroup.getGroup(),
                fieldmodule.createFieldOr(
                    fieldmodule.createFieldOr(neckGroup.getGroup(), armGroup.getGroup()),
                    fieldmodule.createFieldOr(shellGroup.getGroup(), abdominalCavityGroup.getGroup())))
            thoracicCavityBoundaryGroup.getMeshGroup(mesh2d).addElementsConditional(is_thoracic_cavity_boundary)

            abdominalCavityBoundaryGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_body_term("abdominal cavity boundary surface"))
            is_abdominal_cavity_boundary = fieldmodule.createFieldAnd(
                abdominalCavityGroup.getGroup(),
                fieldmodule.createFieldOr(
                    thoracicCavityGroup.getGroup(),
                    fieldmodule.createFieldOr(shellGroup.getGroup(), legGroup.getGroup())))
            abdominalCavityBoundaryGroup.getMeshGroup(mesh2d).addElementsConditional(is_abdominal_cavity_boundary)

            diaphragmGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("diaphragm"))
            is_diaphragm = fieldmodule.createFieldAnd(thoracicCavityGroup.getGroup(), abdominalCavityGroup.getGroup())
            diaphragmGroup.getMeshGroup(mesh2d).addElementsConditional(is_diaphragm)

            spinalCordGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("spinal cord"))
            is_spinal_cord = fieldmodule.createFieldAnd(is_core_shell, is_left_right_dorsal)
            spinalCordGroup.getMeshGroup(mesh1d).addElementsConditional(is_spinal_cord)


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
