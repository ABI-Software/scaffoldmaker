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
        options["Left elbow flexion degrees"] = 90.0
        options["Right elbow flexion degrees"] = 110.0
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
        options["Pelvis drop"] = 0.8
        options["Pelvis width"] = 2.0
        options["Left leg abduction degrees"] = 10.0
        options["Right leg abduction degrees"] = 10.0
        options["Left hip flexion degrees"] = 90.0
        options["Right hip flexion degrees"] = 90.0
        options["Leg length"] = 11.0
        options["Leg top diameter"] = 2.0
        options["Leg bottom diameter"] = 0.7
        options["Left knee flexion degrees"] = 90.0
        options["Right knee flexion degrees"] = 90.0
        options["Left ankle flexion degrees"] = 90.0
        options["Right ankle flexion degrees"] = 90.0
        options["Foot height"] = 1.25
        options["Foot length"] = 2.0
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
            "Left hip flexion degrees",
            "Right hip flexion degrees",
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
            "Left hip flexion degrees": (0.0, 150.0),
            "Right hip flexion degrees": (0.0, 150.0),
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
        hipLeftFlexionRadians = math.radians(options["Left hip flexion degrees"])
        hipRightFlexionRadians = math.radians(options["Right hip flexion degrees"])
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
        leftShoulderGroup = AnnotationGroup(region, get_body_term("left shoulder"))
        leftBrachiumGroup = AnnotationGroup(region, get_body_term("left brachium"))
        leftAntebrachiumGroup = AnnotationGroup(region, get_body_term("left antebrachium"))
        # leftElbowGroup = AnnotationGroup(region, get_body_term("left elbow"))
        leftHandGroup = AnnotationGroup(region, get_body_term("left hand"))
        rightArmGroup = AnnotationGroup(region, get_body_term("right upper limb"))
        rightShoulderGroup = AnnotationGroup(region, get_body_term("right shoulder"))
        rightBrachiumGroup = AnnotationGroup(region, get_body_term("right brachium"))
        rightAntebrachiumGroup = AnnotationGroup(region, get_body_term("right antebrachium"))
        # rightElbowGroup = AnnotationGroup(region, get_body_term("right elbow"))
        rightHandGroup = AnnotationGroup(region, get_body_term("right hand"))
        handGroup = AnnotationGroup(region, get_body_term("hand"))
        thoraxGroup = AnnotationGroup(region, get_body_term("thorax"))
        abdomenGroup = AnnotationGroup(region, get_body_term("abdomen"))
        hipGroup = AnnotationGroup(region, get_body_term("hip"))
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
                            thoraxGroup, abdomenGroup, hipGroup,
                            leftShoulderGroup, leftBrachiumGroup, leftAntebrachiumGroup, leftHandGroup, 
                            rightShoulderGroup, rightBrachiumGroup, rightAntebrachiumGroup, rightHandGroup,
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
        shoulderElementsCount = humanElementCounts['shoulderElementsCount']
        brachiumElementsCount = humanElementCounts['brachiumElementsCount']
        antebrachiumElementsCount = humanElementCounts['antebrachiumElementsCount']
        handElementsCount = humanElementCounts['handElementsCount']
        armToHandElementsCount = shoulderElementsCount + brachiumElementsCount + antebrachiumElementsCount 
        armMeshGroup = armGroup.getMeshGroup(mesh)
        # armToHandMeshGroup = armToHandGroup.getMeshGroup(mesh)
        handMeshGroup = handGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideArmGroup = leftArmGroup if (side == left) else rightArmGroup
            sideShoulderGroup = leftShoulderGroup if (side == left) else rightShoulderGroup
            sideBrachiumGroup = leftBrachiumGroup if (side == left) else rightBrachiumGroup
            sideAntebrachiumGroup = leftAntebrachiumGroup if (side == left) else rightAntebrachiumGroup
            # sideElbowGroup = leftElbowGroup if (side == left) else rightElbowGroup
            sideHandGroup = leftHandGroup if (side == left) else rightHandGroup
            # Setup shoulder elements
            meshGroups = [bodyMeshGroup, 
                          armMeshGroup, 
                          sideArmGroup.getMeshGroup(mesh), sideShoulderGroup.getMeshGroup(mesh)]
            for e in range(shoulderElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
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
        hipElementsCount = humanElementCounts['hipElementsCount']
        upperLegElementsCount = humanElementCounts['upperLegElementsCount']
        lowerLegElementsCount = humanElementCounts['lowerLegElementsCount']
        footElementsCount = humanElementCounts['footElementsCount']
        legToFootElementsCount = hipElementsCount + upperLegElementsCount + lowerLegElementsCount
        legMeshGroup = legGroup.getMeshGroup(mesh)
        hipMeshGroup = hipGroup.getMeshGroup(mesh)
        legToFootMeshGroup = legToFootGroup.getMeshGroup(mesh)
        footMeshGroup = footGroup.getMeshGroup(mesh)
        for side in (left, right):
            sideLegGroup = leftLegGroup if (side == left) else rightLegGroup
            sideUpperLegGroup = leftUpperLegGroup if (side == left) else rightUpperLegGroup
            sideLowerLegGroup = leftLowerLegGroup if (side == left) else rightLowerLegGroup
            sideFootGroup = leftFootGroup if (side == left) else rightFootGroup
            # Hip
            meshGroups = [bodyMeshGroup, legMeshGroup, hipMeshGroup,
                          sideLegGroup.getMeshGroup(mesh), sideUpperLegGroup.getMeshGroup(mesh)]
            for e in range(hipElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
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
        

        thoraxScale = thoraxLength / thoraxElementsCount
        thoraxStartX = headLength + neckLength
        sx = [thoraxStartX, 0.0, 0.0]
        options['Kinematic tree']['thorax_top'] = sx
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
        # options['Kinematic tree']['lumbar_body'] = px
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
            # shoulderElementsCount = 2
            armScale = nonHandArmLength / (armToHandElementsCount - shoulderElementsCount)
            d12_mag = (halfWristThickness - armTopRadius) / (armToHandElementsCount - shoulderElementsCount)
            d13_mag = (halfWristWidth - armTopRadius) / (armToHandElementsCount - shoulderElementsCount)
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
            flexionRotFactor = 1.0*math.sin(shoulderFlexionRadians)*(math.sqrt(2)-1)      
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
            for i in range(1, brachiumElementsCount):
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
            # Elbow
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
                armSide = d2 
                armFront = d3
            # Updating frame of reference wrt flexion angle (using d2 as rotation axis)
            elbowFlexionRadians = elbowLeftFlexionRadians if (side == left) else elbowRightFlexionRadians
            armDir = [armDirn, armSide, armFront]
            i += 1
            xi = i / (armToHandElementsCount - 2)
            halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
            halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
            rotationCoeff = 0.25
            [x, elbowDir, antebrachiumStart, antebrachiumDir, elbowd3_mag, elbowd13_mag] =\
                  getJointAndDistalFlexionFrames(elbowFlexionRadians, x, armDir,halfWidth, armScale, rotationCoeff,ventralFlexion=True)
            # Set coordiantes for joint node
            # x = jointPositions[1]
            elbowDirn, elbowSide, elbowFront = elbowDir
            d1 = mult(elbowDirn, armScale)
            # elbowFront = cross(d1, d2)
            d2 = mult(elbowSide, halfThickness)
            d3 = mult(elbowFront, elbowd3_mag)
            d12 = mult(elbowSide, d12_mag)
            d13 = add(
                mult(elbowFront, d13_mag), 
                mult(elbowDirn, elbowd13_mag)
                )
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
            # antebrachiumStart = jointPositions[-1]
            antebrachiumDirn, antebrachiumSide, antebrachiumFront = antebrachiumDir
            d1 = set_magnitude(antebrachiumDirn, armScale)
            # Change d1 to the antebrachium direction
            j=0
            for i in range(brachiumElementsCount + 1, armToHandElementsCount - 1):
                xi = (i) / (armToHandElementsCount - 2)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(antebrachiumStart, mult(d1, j)) 
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
                j += 1
            # Hand twist
            options['Kinematic tree']['hand_' + side_label] = x
            assert handElementsCount == 1
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            hx = add(x, mult(antebrachiumDirn, handLength))
            hd1 = computeCubicHermiteEndDerivative(x, d1, hx, d1)
            twistAngle = armTwistAngleRadians if (side == left) else -armTwistAngleRadians
            if twistAngle >= 0.0:
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
        d12_mag = (legBottomRadius - legTopRadius) / (legToFootElementsCount)
        d13_mag = (legBottomRadius - legTopRadius) / (legToFootElementsCount)
        pd3 = [0.0, 0.0, 0.5 * legTopRadius + 0.5 * halfTorsoDepth]
        pid3 = mult(pd3, innerProportionDefault)
        for side in (left, right):
            side_label = 'l' if (side == left) else 'r'
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
            # options['Kinematic tree']['femur_' + side_label] = x
            # Upper leg
            for i in range(hipElementsCount-1):
                xi = i / legToFootElementsCount
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(legStart, mult(d1, i))
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = mult(legSide, radius)
                d3 = mult(legFront, radius)
                d13 = mult(legFront, d13_mag)
                d12 = mult(legSide, d12_mag)
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id13 = mult(d13, innerProportionDefault)
                id12 = mult(d12, innerProportionDefault)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
            # Frontal hip flexion
            # Updating frame of reference wrt flexion angle (using d2 as rotation axis)
            hipFlexionRadians = hipLeftFlexionRadians if (side == left) else hipRightFlexionRadians
            i += 1
            xi = i / legToFootElementsCount
            radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
            legDir = [legDirn, legSide, legFront]
            rotationCoeff = 0.18
            [x, hipDir, upperLegStart, upperLegDir, hipd3_mag, hipd13_mag] =\
                  getJointAndDistalFlexionFrames(hipFlexionRadians, x, legDir,radius, legScale, rotationCoeff,ventralFlexion=True)
            hipDirn, hipSide, hipFront = hipDir
            d1 = mult(hipDirn, legScale)
            d2 = mult(hipSide, radius)
            d3 = mult(hipFront, 0.9*hipd3_mag)
            d12 = mult(hipSide, d12_mag)
            d13 = add(
                mult(hipFront, d13_mag), 
                mult(hipDirn, 0.7*hipd13_mag)
                )
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1
            options['Kinematic tree']['femur_' + side_label] = x
            upperLegDirn, upperLegSide, upperLegFront = upperLegDir
            d1 = set_magnitude(upperLegDirn, armScale)
            # Rest of upper leg
            j = 0
            for i in range(hipElementsCount, hipElementsCount+upperLegElementsCount-1):
                xi = i / legToFootElementsCount
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(upperLegStart, mult(d1, j))
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = mult(upperLegSide, radius)
                d3 = mult(upperLegFront, radius)
                d12 = mult(upperLegSide, d12_mag)
                d13 = mult(upperLegFront, d13_mag)
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id12 = set_magnitude(d12, d12_mag)
                id13 = set_magnitude(d13, d13_mag)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
                j += 1
            # knee
            kneeFlexionRadians = kneeLeftFlexionRadians if (side == left) else kneeRightFlexionRadians
            # # Set coordiantes for joint node
            i += 1
            xi = i / legToFootElementsCount
            radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
            upperLegDir = [upperLegDirn, upperLegSide, upperLegFront]
            rotationCoeff = 0.25
            [x, kneeDir, lowerLegStart, lowerLegDir, kneed3_mag, kneed13_mag] =\
                  getJointAndDistalFlexionFrames(kneeFlexionRadians, x, upperLegDir,radius, legScale, rotationCoeff,ventralFlexion=False)
            kneeDirn, kneeSide, kneeFront = kneeDir
            d1 = mult(kneeDirn, legScale)
            d2 = mult(kneeSide, radius)
            d3 = mult(kneeFront, kneed3_mag)
            d12 = mult(kneeSide, d12_mag)
            d13 = add(
                set_magnitude(d3, d13_mag),
                set_magnitude(d1, kneed13_mag)
            )
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1
            options['Kinematic tree']['tibia_' + side_label] = x
            # Lower leg
            lowerLegDirn, lowerLegSide, lowerLegFront = lowerLegDir
            d1 = set_magnitude(lowerLegDirn, legScale)
            j = 0 
            for i in range(hipElementsCount+upperLegElementsCount, legToFootElementsCount-1):
                xi = i / legToFootElementsCount
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(lowerLegStart, mult(d1, j))
                radius = xi * legBottomRadius + (1.0 - xi) * legTopRadius
                d2 = set_magnitude(lowerLegSide, radius)
                d3 = set_magnitude(lowerLegFront, radius)
                d12 = set_magnitude(d2, d12_mag)
                d13 = set_magnitude(d3, d13_mag)
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id12 = set_magnitude(d12, d12_mag)
                id13 = set_magnitude(d13, d13_mag)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
                j+=1
            # foot
            ankleFlexionRadians = ankleLeftFlexionRadians if (side == left) else ankleRightFlexionRadians
            i += 1
            radius = (halfFootThickness+legBottomRadius)/2
            rotationCoeff = 0.15
            [x, ankleDir, footStart, footDir, ankled3_mag, ankled13_mag] =\
                  getJointAndDistalFlexionFrames(
                      ankleFlexionRadians, x, lowerLegDir,radius, 
                      legScale, rotationCoeff,ventralFlexion=True, distalnScale=footLength)
            ankleDirn, ankleSide, ankleFront = ankleDir
            footDirn, footSide, footFront = footDir
            d1 = mult(footDirn, footLength)
            d2 = set_magnitude(ankleSide, halfFootWidth)
            d3 = set_magnitude(ankleFront, ankled3_mag)
            d12 = set_magnitude(ankleSide, d12_mag)
            d13 = add(
                set_magnitude(ankleFront, d13_mag), 
                set_magnitude(footDirn, ankled13_mag)
            )
            id2 = mult(d2, innerProportionDefault)
            id3 = mult(d3, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            id13 = mult(d13, innerProportionDefault)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            nodeIdentifier += 1                 
            # Foot end nodes
            
            d1 = mult(footDirn, footLength)
            j = 0 
            for i in range(footElementsCount):
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                x = add(footStart, mult(d1, j))
                d2 = set_magnitude(footSide, halfFootWidth)
                d3 = set_magnitude(footFront, halfFootThickness)
                d12 = set_magnitude(d2, d12_mag)
                d13 = sub(d3, set_magnitude(ankleFront, ankled3_mag))
                id2 = mult(d2, innerProportionDefault)
                id3 = mult(d3, innerProportionDefault)
                id12 = set_magnitude(d12, d12_mag)
                id13 = set_magnitude(d13, d13_mag)
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
                nodeIdentifier += 1
                j+=1
            # x = footStart
            # footDirn, footSide, footFront = footDir
            # d1 = mult(footDirn, footLength)
            # # d1 = computeCubicHermiteEndDerivative(jointPositions[1], jointDirn[1], jointPositions[2], jointDirn[2])
            # d2 = set_magnitude(footSide, halfFootWidth)
            # d3 = set_magnitude(footFront, halfFootThickness)
            # d12 = set_magnitude(d2, d12_mag)
            # # d13 = set_magnitude(d3, d13_mag)
            # d13 = sub(d3, set_magnitude(ankleFront, ankled3_mag))
            # id2 = mult(d2, innerProportionDefault)
            # id3 = mult(d3, innerProportionDefault)
            # id12 = set_magnitude(d12, d12_mag)
            # id13 = set_magnitude(d13, d13_mag)
            # node = nodes.findNodeByIdentifier(nodeIdentifier)
            # fieldcache.setNode(node)
            # setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
            # setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
            # nodeIdentifier += 1
            options['Kinematic tree']['toes_' + side_label] = x
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
        options["Number of elements along shoulder"] = 2
        options["Number of elements along brachium"] = 3
        options["Number of elements along antebrachium"] = 3
        options["Number of elements along hand"] = 1
        options["Number of elements along hip"] = 2
        options["Number of elements along upper leg"] = 3
        options["Number of elements along lower leg"] = 3
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
            options["Number of elements along shoulder"] = 2
            options["Number of elements along brachium"] = 3
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
            options["Number of elements along shoulder"] = 2
            options["Number of elements along brachium"] = 3
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
            "Number of elements along shoulder",
            "Number of elements along brachium",
            "Number of elements along antebrachium",
            "Number of elements along hand",
            "Number of elements along hip",
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
        elementsCountAlongShoulder = options["Number of elements along shoulder"]
        elementsCountAlongBrachium = options["Number of elements along brachium"]
        elementsCountAlongAntebrachium = options["Number of elements along antebrachium"]
        elementsCountAlongHand = options["Number of elements along hand"]
        elementsCountAlongHip = options["Number of elements along hip"]
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
            elif "shoulder" in name:
                alongCount = elementsCountAlongShoulder
                aroundCount = elementsCountAroundArm
            elif " brachium" in name:
                alongCount = elementsCountAlongBrachium
                aroundCount = elementsCountAroundArm
            elif " antebrachium" in name:
                alongCount = elementsCountAlongAntebrachium
                aroundCount = elementsCountAroundArm
            elif "hand" in name:
                alongCount = elementsCountAlongHand
                aroundCount = elementsCountAroundArm
            elif "hip" in name:
                alongCount = elementsCountAlongHip
                aroundCount = elementsCountAroundLeg
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

def getJointAndDistalFlexionFrames(jointFlexionRadians, proximalNodePosition, proximalDir, 
                             frontScale, dirnScale, rotationCoeff, distalnScale = False, ventralFlexion=True):
    jointAngleRadians = math.pi - jointFlexionRadians
    proximalDirn, proximalSide, proximalFront = proximalDir
    ventral = 1 if ventralFlexion else -1 
    jointRotationMatrix = axis_angle_to_rotation_matrix(mult(proximalSide, -1*ventral), jointFlexionRadians)
    jointHalfRotationMatrix = axis_angle_to_rotation_matrix(mult(proximalSide, -1*ventral), jointFlexionRadians/2)
    # Joint directions (frame is rotated by half the flexion angle)
    jointDirn = matrix_vector_mult(jointHalfRotationMatrix, proximalDirn)
    jointSide = proximalSide
    jointFront = matrix_vector_mult(jointHalfRotationMatrix, proximalFront)
    # Proximal directions
    distalDirn = matrix_vector_mult(jointRotationMatrix, proximalDirn)
    distalSide = proximalSide
    distalFront = matrix_vector_mult(jointRotationMatrix, proximalFront)
    # These rotation factors are used to adjust the position of the joint node relative 
    # to the angle of flexion, and ensures a proper transition between the two parts
    flexionRotFactor = 1*math.sin(jointAngleRadians)     
    jointRotFactor = 1/math.sin(jointAngleRadians/2)
    d13RotFactor = math.sqrt(2)*math.tan(jointFlexionRadians/2)
    # Diagram to calculate position of joint node with flexion
    # 1 
    # |
    # 2 -- 3 
    # 1 is the proximal node, 2 is the joint node, 3 is distal node
    # we fix the position of the nodes 1 and 3, and calculate 
    # the position of 2 using the sampleCubicHermiteCurvesSmooth function.
    nodePositions = []
    nodePositions.append(proximalNodePosition) #1
    # The joint node is not set directly at the corner
    # Instead it is nudged forward towards the center of the tube
    # To get as smooth of a line as possible between node #1 and #3
    rotDisplacementFactor = rotationCoeff*frontScale*flexionRotFactor
    proximalScale = dirnScale
    distalScale = distalnScale if (distalnScale) else dirnScale
    jointAdjustDir = add(
                    set_magnitude(jointFront, ventral*rotDisplacementFactor), 
                    set_magnitude(distalDirn, rotDisplacementFactor), 
                )
    jointAdjustDir = add(
                set_magnitude(proximalDirn, proximalScale), 
                jointAdjustDir
            )
    jointAdjustDir = set_magnitude(jointAdjustDir, proximalScale)
    nodePositions.append(add(nodePositions[-1], jointAdjustDir)) #2
    jointAdjustDir = add(
            set_magnitude(jointDirn, rotDisplacementFactor), 
            set_magnitude(proximalDirn, rotDisplacementFactor), 
        )
    jointAdjustDir =  add(
        set_magnitude(distalDirn, distalScale), 
        jointAdjustDir
        )
    jointAdjustDir = set_magnitude(jointAdjustDir, distalScale)
    nodePositions.append(add(nodePositions[-1], jointAdjustDir)) #3
    # nodeDir = [proximalDirn, jointDirn, distalDirn]
    # nodePositions, nodeDirn = sampleCubicHermiteCurvesSmooth(
    # nodePositions, jo, 2, 
    # derivativeMagnitudeStart=dirnScale, derivativeMagnitudeEnd=dirnScale
    # )[0:2]  
    jointNodePosition = nodePositions[1]
    jointDir = [jointDirn, jointSide, jointFront]
    distalNodePosition = nodePositions[2]
    distalDir = [distalDirn, distalSide, distalFront]
    # Magnitudes
    jointd3_mag = frontScale*(jointRotFactor-(rotationCoeff*flexionRotFactor))
    jointd13_mag = -1*ventral*frontScale*d13RotFactor
    return [jointNodePosition, jointDir, distalNodePosition, distalDir, jointd3_mag, jointd13_mag]