"""
Generates a 3D renal pelvis using tube network mesh.
"""
import math

from cmlibs.maths.vectorops import mult, cross, add, sub, set_magnitude, rotate_about_z_axis, \
    rotate_vector_around_vector
from cmlibs.utils.zinc.field import find_or_create_field_coordinates, findOrCreateFieldCoordinates
from cmlibs.zinc.field import Field

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.annotation.ureter_terms import get_ureter_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshGenerateData, RenalPelvisTubeNetworkMeshBuilder
from cmlibs.zinc.node import Node


class MeshType_1d_renal_pelvis_network_layout1(MeshType_1d_network_layout1):
    """
    Defines renal pelvis network layout.
    """

    isHuman = True
    isRat = False

    @classmethod
    def getName(cls):
        return "1D Renal Pelvis Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1",
            "Rat 1"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = "Human 1" if (parameterSetName == "Default") else parameterSetName

        options["Define inner coordinates"] = True
        options["Inner proportion default"] = 0.8
        options["Inner proportion ureter"] = 0.7

        if parameterSetName in ["Default", "Human 1"]:
            cls.isHuman = True
            cls.isRat = False

            options["Top major calyx"] = True
            options["Middle major calyx"] = True
            options["Bottom major calyx"] = True
            options["Upper minor calyx"] = True
            options["Lower minor calyx"] = True
            options["Rotate upper, middle and lower minor calyxes"] = True
            options["Number of calyxes at top minor calyx"] = 2
            options["Number of calyxes at middle major calyx"] = 2
            options["Number of calyxes at bottom minor calyx"] = 2
            options["Number of calyxes at upper minor calyx"] = 2
            options["Number of calyxes at lower minor calyx"] = 2
            options["Ureter length"] = 4.0
            options["Ureter radius"] = 0.1
            options["Ureteropelvic junction width"] = 0.3
            options["Ureter bend angle degrees"] = 45
            options["Major calyx length"] = 0.7
            options["Major calyx radius"] = 0.15
            options["Major calyx angle degrees"] = 200
            options["Middle major calyx length"] = 0.4
            options["Major to bottom/top minor calyx length"] = 0.4
            options["Major to lower/upper minor calyx length"] = 0.3
            options["Bottom/top minor calyx length"] = 0.2
            options["Lower/upper minor calyx length"] = 0.2
            options["Minor calyx radius"] = 0.1
            options["Bottom/top minor calyx bifurcation angle degrees"] = 90
            options["Bottom/top minor calyx rotate angle degrees"] = 30
            options["Lower/upper minor calyx bifurcation angle degrees"] = 90
            options["Lower/upper minor calyx bend angle degrees"] = 10
            options["Renal pyramid length"] = 0.6
            options["Renal pyramid width"] = 0.6
            options["Renal pyramid thickness"] = 0.2

            options["Structure"] = cls.getPelvisLayoutStructure(options)

        elif "Rat 1" in parameterSetName:
            options["Structure"] = (
                "1-2, 2-3.1,3.2-4,4-5,5-6-7-8"
            )
            cls.isHuman = False
            cls.isRat = True

            options["Define inner coordinates"] = True
            options["Top major calyx"] = False
            options["Middle major calyx"] = True
            options["Bottom major calyx"] = False
            options["Upper minor calyx"] = False
            options["Lower minor calyx"] = False
            options["Rotate upper, middle and lower minor calyxes"] = True
            options["Number of calyxes at top minor calyx"] = 1
            options["Number of calyxes at middle major calyx"] = 1
            options["Number of calyxes at bottom minor calyx"] = 1
            options["Number of calyxes at upper minor calyx"] = 1
            options["Number of calyxes at lower minor calyx"] = 1
            options["Ureter length"] = 3.0
            options["Ureter radius"] = 0.1
            options["Ureteropelvic junction width"] = options["Ureter radius"]
            options["Ureter bend angle degrees"] = 45
            options["Major calyx length"] = 0.6
            options["Major calyx radius"] = 0.1
            options["Major calyx angle degrees"] = 170
            options["Middle major calyx length"] = 0.4
            options["Major to bottom/top minor calyx length"] = 0.3
            options["Major to lower/upper minor calyx length"] = 0.3
            options["Bottom/top minor calyx length"] = 0.2
            options["Lower/upper minor calyx length"] = 0.2
            options["Minor calyx radius"] = 0.1
            options["Bottom/top minor calyx bifurcation angle degrees"] = 90
            options["Bottom/top minor calyx rotate angle degrees"] = 0
            options["Lower/upper minor calyx bifurcation angle degrees"] = 90
            options["Lower/upper minor calyx bend angle degrees"] = 10
            options["Renal pyramid length"] = 0.6
            options["Renal pyramid width"] = 0.6
            options["Renal pyramid thickness"] = 0.2

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        if cls.isHuman:
            return [
                "Top major calyx",
                "Middle major calyx",
                "Bottom major calyx",
                "Upper minor calyx",
                "Lower minor calyx",
                "Rotate upper, middle and lower minor calyxes",
                "Number of calyxes at top minor calyx",
                "Number of calyxes at upper minor calyx",
                "Number of calyxes at middle major calyx",
                "Number of calyxes at lower minor calyx",
                "Number of calyxes at bottom minor calyx",
                "Ureter length",
                "Ureter radius",
                "Ureteropelvic junction width",
                "Ureter bend angle degrees",
                "Major calyx length",
                "Major calyx radius",
                "Major calyx angle degrees",
                "Middle major calyx length",
                "Major to bottom/top minor calyx length",
                "Major to lower/upper minor calyx length",
                "Bottom/top minor calyx length",
                "Lower/upper minor calyx length",
                "Minor calyx radius",
                "Bottom/top minor calyx bifurcation angle degrees",
                "Bottom/top minor calyx rotate angle degrees",
                "Lower/upper minor calyx bifurcation angle degrees",
                "Lower/upper minor calyx bend angle degrees",
                "Renal pyramid length",
                "Renal pyramid width",
                "Renal pyramid thickness",
                "Inner proportion default",
                "Inner proportion ureter"
            ]
        elif cls.isRat:
            return [
                "Ureter length",
                "Ureter radius",
                "Ureter bend angle degrees",
                "Middle major calyx length",
                "Minor calyx radius",
                "Renal pyramid length",
                "Renal pyramid width",
                "Renal pyramid thickness",
                "Inner proportion default",
                "Inner proportion ureter"
            ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Ureter length",
            "Ureter radius",
            "Ureteropelvic junction width",
            "Ureter bend angle degrees",
            "Major calyx length",
            "Major calyx radius",
            "Major calyx angle degrees",
            "Middle major calyx length",
            "Major to bottom/top minor calyx length",
            "Major to lower/upper minor calyx length",
            "Bottom/top minor calyx length",
            "Lower/upper minor calyx length",
            "Minor calyx radius",
            "Bottom/top minor calyx bifurcation angle degrees",
            "Lower/upper minor calyx bifurcation angle degrees",
            "Lower/upper minor calyx bend angle degrees",
            "Renal pyramid length",
            "Renal pyramid width",
            "Renal pyramid thickness"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1

        for key in [
            "Number of calyxes at top minor calyx",
            "Number of calyxes at upper minor calyx",
            "Number of calyxes at middle major calyx",
            "Number of calyxes at lower minor calyx",
            "Number of calyxes at bottom minor calyx"
        ]:
            if options[key] < 1:
                options[key] = 1
            elif options[key] > 2:
                options[key] = 2

        for key in [
            "Inner proportion default",
            "Inner proportion ureter"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 0.9:
                options[key] = 0.9

        for key, angleRange in {
            "Ureter bend angle degrees": (0.0, 45.0),
            "Major calyx angle degrees": (130.0, 200.0),
            "Bottom/top minor calyx bifurcation angle degrees": (60.0, 120.0),
            "Lower/upper minor calyx bifurcation angle degrees": (60.0, 120.0),
            "Lower/upper minor calyx bend angle degrees": (0.0, 10.0)
        }.items():
            if options[key] < angleRange[0]:
                options[key] = angleRange[0]
            elif options[key] > angleRange[1]:
                options[key] = angleRange[1]

        options["Structure"] = cls.getPelvisLayoutStructure(options)

        return dependentChanges

    @classmethod
    def getPelvisLayoutStructure(cls, options):
        """
        Returns the 1D layout structure of the renal pelvis depending on user inputs.
        """
        isTopMC = options["Top major calyx"]
        isMidMC = options["Middle major calyx"]
        isBottomMC = options["Bottom major calyx"]
        isUpperMC = options["Upper minor calyx"]
        isLowerMC = options["Lower minor calyx"]

        nTopMinorCalyxes = options["Number of calyxes at top minor calyx"]
        nMidMajorCalyxes = options["Number of calyxes at middle major calyx"]
        nBottomMinorCalyxes = options["Number of calyxes at bottom minor calyx"]
        nUpperMinorCalyxes = options["Number of calyxes at upper minor calyx"]
        nLowerMinorCalyxes = options["Number of calyxes at lower minor calyx"]

        nTopMajorCalyxes = 2 if isUpperMC else 1
        nBottomMajorCalyxes = 2 if isLowerMC else 1

        structure = "1-2, 2-3.1"
        nodeIdentifier = 4
        nodeDerivative = 2

        # Build major calyx connections
        majorCalyxNodeIdentifiers = []
        majorCalyxFlags = [isBottomMC, isMidMC, isTopMC]

        for i, hasMajorCalyx in enumerate(majorCalyxFlags):
            if hasMajorCalyx:
                majorCalyxNodeIdentifiers.append(nodeIdentifier)
                layout = f"3.{nodeDerivative}-{nodeIdentifier}"
                structure += f",{layout}"
                nodeDerivative += 1
                nodeIdentifier += 1

        # Build bottom major calyx to minor calyx connections
        bmcNodeIdentifiers = []
        if isBottomMC:
            startNode = majorCalyxNodeIdentifiers[0]
            for n in range(nBottomMajorCalyxes):
                bmcNodeIdentifiers.append(nodeIdentifier)
                layout = f"{startNode}.{n + 2}-{nodeIdentifier}" if isLowerMC else f"{startNode}-{nodeIdentifier}"
                structure += f",{layout}"
                nodeIdentifier += 1

        # Build top major calyx to minor calyx connections
        tmcNodeIdentifiers = []
        if isTopMC:
            startNode = majorCalyxNodeIdentifiers[-1]
            for n in range(nTopMajorCalyxes):
                tmcNodeIdentifiers.append(nodeIdentifier)
                layout = f"{startNode}.{n + 2}-{nodeIdentifier}" if isUpperMC else f"{startNode}-{nodeIdentifier}"
                structure += f",{layout}"
                nodeIdentifier += 1

        # Build minor calyx connections
        minorCalyxNodeIdentifiers = []

        # Bottom minor calyxes
        if isBottomMC:
            calyxCounts = [nBottomMinorCalyxes, nLowerMinorCalyxes]
            for i in range(nBottomMajorCalyxes):
                nCalyxes = calyxCounts[i]
                startNode = bmcNodeIdentifiers[i]
                for j in range(nCalyxes):
                    minorCalyxNodeIdentifiers.append(nodeIdentifier)
                    layout = f"{startNode}.{j + 2}-{nodeIdentifier}" if nCalyxes > 1 else f"{startNode}-{nodeIdentifier}"
                    structure += f",{layout}"
                    nodeIdentifier += 1

        # Middle minor calyxes
        if isMidMC:
            # Determine correct index for middle major calyx
            if len(majorCalyxNodeIdentifiers) == 3:  # All three major calyxes present
                midIndex = 1
            elif len(majorCalyxNodeIdentifiers) == 2 and isTopMC:  # Mid and Top present
                midIndex = 0
            else:  # Mid and Bottom present, or only Mid present
                midIndex = -1 if len(majorCalyxNodeIdentifiers) == 2 else 0

            startNode = majorCalyxNodeIdentifiers[midIndex]
            for n in range(nMidMajorCalyxes):
                minorCalyxNodeIdentifiers.append(nodeIdentifier)
                layout = f"{startNode}.{n + 2}-{nodeIdentifier}" if nMidMajorCalyxes > 1 else f"{startNode}-{nodeIdentifier}"
                structure += f",{layout}"
                nodeIdentifier += 1

        # Top minor calyxes
        if tmcNodeIdentifiers:
            calyxCounts = [nTopMinorCalyxes, nUpperMinorCalyxes]
            for i in range(nTopMajorCalyxes):
                nCalyxes = calyxCounts[i]
                startNode = tmcNodeIdentifiers[i]
                for j in range(nCalyxes):
                    minorCalyxNodeIdentifiers.append(nodeIdentifier)
                    layout = f"{startNode}.{j + 2}-{nodeIdentifier}" if nCalyxes > 1 else f"{startNode}-{nodeIdentifier}"
                    structure += f",{layout}"
                    nodeIdentifier += 1

        # Add final connections for each minor calyx
        for nid in minorCalyxNodeIdentifiers:
            layout = f"{nid}-{nodeIdentifier}-{nodeIdentifier + 1}-{nodeIdentifier + 2}"
            structure += f",{layout}"
            nodeIdentifier += 3

        return structure

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return [] empty list of AnnotationGroup, NetworkMesh
        """
        # parameters
        structure = options["Structure"]

        isTopMC = options["Top major calyx"]
        isMidMC = options["Middle major calyx"]
        isBottomMC = options["Bottom major calyx"]
        isUpperMC = options["Upper minor calyx"]
        isLowerMC = options["Lower minor calyx"]
        isRotateMinorCalyx = options["Rotate upper, middle and lower minor calyxes"]

        nTopMinorCalyxes = options["Number of calyxes at top minor calyx"]
        nUpperMinorCalyxes = options["Number of calyxes at upper minor calyx"]
        nMidMinorCalyxes = options["Number of calyxes at middle major calyx"]
        nLowerMinorCalyxes = options["Number of calyxes at lower minor calyx"]
        nBottomMinorCalyxes = options["Number of calyxes at bottom minor calyx"]
        nMinorCalyxesList = [nBottomMinorCalyxes, nLowerMinorCalyxes, nMidMinorCalyxes, nTopMinorCalyxes, nUpperMinorCalyxes]

        ureterLength = options["Ureter length"]
        ureterRadius = options["Ureter radius"]
        upjWidth = options["Ureteropelvic junction width"]
        ureterBendAngle = options["Ureter bend angle degrees"]
        majorCalyxLength = options["Major calyx length"]
        majorCalyxRadius = options["Major calyx radius"]
        majorCalyxAngle = options["Major calyx angle degrees"]
        majorToBottomMinorCalyxLength = options["Major to bottom/top minor calyx length"]
        majorToLowerMinorCalyxLength = options["Major to lower/upper minor calyx length"]
        middleMajorLength = options["Middle major calyx length"]
        bottomMinorCalyxLength = options["Bottom/top minor calyx length"]
        minorCalyxRadius = options["Minor calyx radius"]
        lowerMinorCalyxLength = options["Lower/upper minor calyx length"]
        bottomMinorCalyxAngle = options["Bottom/top minor calyx bifurcation angle degrees"]
        minorCalyxRotateAngle = options["Bottom/top minor calyx rotate angle degrees"]
        lowerMinorCalyxAngle = options["Lower/upper minor calyx bifurcation angle degrees"]
        minorCalyxBendAngle = options["Lower/upper minor calyx bend angle degrees"]
        pyramidLength = options["Renal pyramid length"]
        pyramidWidth = options["Renal pyramid width"]
        pyramidThickness = options["Renal pyramid thickness"]
        innerProportionDefault = options["Inner proportion default"]
        innerProportionUreter = options["Inner proportion ureter"]

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)

        # set up element annotations
        renalPelvisGroup = AnnotationGroup(region, get_kidney_term("renal pelvis"))
        ureterGroup = AnnotationGroup(region, get_ureter_term("ureter"))
        majorCalyxGroup = AnnotationGroup(region, get_kidney_term("major calyx"))
        minorCalyxGroup = AnnotationGroup(region, get_kidney_term("minor calyx"))
        renalPyramidGroup = AnnotationGroup(region, get_kidney_term("renal pyramid"))
        annotationGroups = [renalPelvisGroup, ureterGroup, majorCalyxGroup, minorCalyxGroup, renalPyramidGroup]

        # Ureter elements
        renalPelvisMeshGroup = renalPelvisGroup.getMeshGroup(mesh)
        elementIdentifier = 1
        ureterElementsCount = 2
        meshGroups = [renalPelvisMeshGroup, ureterGroup.getMeshGroup(mesh)]
        for e in range(ureterElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        # Major calyx elements
        majorCalyxElementsCount = sum([isBottomMC, isMidMC, isTopMC])
        meshGroups = [renalPelvisMeshGroup, majorCalyxGroup.getMeshGroup(mesh)]
        majorCalyxFlags = [isBottomMC, isMidMC, isTopMC]

        totalMajorCalyxes = sum(majorCalyxFlags)
        startIndex = 1 if (totalMajorCalyxes == 2 and isTopMC) or (totalMajorCalyxes == 1) else 0

        middleMajorCalyxIdentifier = None
        for e in range(startIndex, majorCalyxElementsCount + startIndex):
            element = mesh.findElementByIdentifier(elementIdentifier)
            if isMidMC and e == 1:  # middle major calyx index
                middleMajorCalyxIdentifier = elementIdentifier
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        # Minor calyx elements
        meshGroups = [renalPelvisMeshGroup, minorCalyxGroup.getMeshGroup(mesh)]

        bottomMinor, lowerMinor, middleMajor, topMinor, upperMinor = 0, 1, 2, 3, 4
        calyxes = [bottomMinor, lowerMinor, middleMajor, upperMinor, topMinor]
        flags = [isBottomMC, isLowerMC, isMidMC, isUpperMC, isTopMC]
        minorCalyxList = [c for c, f in zip(calyxes, flags) if f]

        for calyx in minorCalyxList:
            cElementIdentifier = middleMajorCalyxIdentifier if (calyx == middleMajor and isMidMC) else elementIdentifier
            element = mesh.findElementByIdentifier(cElementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            if not (calyx == middleMajor and isMidMC):
                elementIdentifier += 1

        # Minor calyx to pyramid connections
        renalPyramidMeshGroup = renalPyramidGroup.getMeshGroup(mesh)
        minorCalyxMeshGroup = minorCalyxGroup.getMeshGroup(mesh)
        meshGroups = [renalPelvisMeshGroup, minorCalyxMeshGroup]
        minorCalyxElementsCount = 1

        for calyx in minorCalyxList:
            for side in range(nMinorCalyxesList[calyx]):
                for e in range(minorCalyxElementsCount):
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)
                    elementIdentifier += 1

        # Pyramid elements
        pyramidElementsCount = 3
        for calyx in minorCalyxList:
            for side in range(nMinorCalyxesList[calyx]):
                for e in range(pyramidElementsCount):
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    renalPyramidMeshGroup.addElement(element)
                    elementIdentifier += 1

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=innerProportionDefault)
        innerCoordinates = find_or_create_field_coordinates(fieldmodule, "inner coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # Generate ureter coordinates
        nodeIdentifier = 1
        ureterScale = ureterLength / ureterElementsCount
        ureterBendAngleRadians = math.radians(ureterBendAngle)
        sinUreterAngle = math.sin(ureterBendAngleRadians)
        cosUreterAngle = math.cos(ureterBendAngleRadians)

        endX = [0.0, 0.0, 0.0]
        tx = endX[0] - ureterLength * cosUreterAngle
        ty = endX[1] - ureterLength * sinUreterAngle
        startX = [tx, ty, 0.0]

        d1 = [ureterScale, 0.0, 0.0]
        td1 = rotate_about_z_axis(d1, 2 * ureterBendAngleRadians)

        sx, sd1 = sampleCubicHermiteCurves([startX, endX], [td1, d1], ureterElementsCount,  arcLengthDerivatives=True)[0:2]
        sd1 = smoothCubicHermiteDerivativesLine(sx, sd1, fixEndDirection=True)

        for e in range(ureterElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            if e < ureterElementsCount:
                d3 = [0.0, 0.0, ureterRadius]
                sd2 = set_magnitude(cross(d3, sd1[e]), ureterRadius)
            else:
                d3 = [0.0, 0.0, majorCalyxRadius]
                sd2 = set_magnitude(cross(d3, sd1[e]), upjWidth)
            id3 = mult(d3, innerProportionUreter)
            sid2 = mult(sd2, innerProportionUreter)

            for field, derivatives in ((coordinates, (sd1[e], sd2, d3)), (innerCoordinates, (sd1[e], sid2, id3))):
                setNodeFieldParameters(field, fieldcache, sx[e], *derivatives)
            nodeIdentifier += 1

        majorCalyxJunctionNodeIdentifier = nodeIdentifier - 1
        majorCalyxStartX = endX

        # Generate major calyx coordinates
        majorCalyxAngleRadians = math.radians(majorCalyxAngle / 2)
        sx = majorCalyxStartX
        bottomMajor, middleMajor, topMajor = 0, 1, 2
        version = 2
        majorCalyxNodeIdentifiers, majorCalyxXList, majorCalyxD1List = [], [], []

        for calyx in (bottomMajor, middleMajor, topMajor):
            if majorCalyxFlags[calyx]:
                majorCalyxNodeIdentifiers.append(nodeIdentifier)
                calyxLength = middleMajorLength if calyx == middleMajor else majorCalyxLength
                majorCalyxAngle = 0 if calyx == middleMajor else majorCalyxAngleRadians
                cosMajorCalyxAngle = math.cos(majorCalyxAngle)
                sinMajorCalyxAngle = math.sin(-majorCalyxAngle if calyx == bottomMajor else majorCalyxAngle)

                majorCalyxEndX = sx[0] + calyxLength * cosMajorCalyxAngle
                majorCalyxEndY = sx[1] + calyxLength * sinMajorCalyxAngle
                x = [majorCalyxEndX, majorCalyxEndY, 0.0]
                sd1 = sub(x, sx)

                majorCalyxXList.append(x)
                majorCalyxD1List.append(sd1)

                for i in range(2):
                    isJunctionNode = i == 0
                    nodeId = majorCalyxJunctionNodeIdentifier if isJunctionNode else nodeIdentifier
                    node = nodes.findNodeByIdentifier(nodeId)
                    fieldcache.setNode(node)

                    sd3 = [0.0, 0.0, majorCalyxRadius]
                    sid3 = mult(sd3, innerProportionDefault)
                    sd2 = set_magnitude(cross(sd3, sd1), majorCalyxRadius)
                    sid2 = mult(sd2, innerProportionDefault)

                    if not isJunctionNode:
                        for field, derivatives in (
                                (coordinates, [x, sd1, sd2, sd3]),
                                (innerCoordinates, [x, sd1, sid2, sid3])
                        ):
                            for valueLabel, value in zip(
                                    (Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                     Node.VALUE_LABEL_D_DS3),
                                    derivatives
                            ):
                                field.setNodeParameters(fieldcache, -1, valueLabel, 1, value)
                        nodeIdentifier += 1

                    setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, sd3)
                    setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, sd3)
                version += 1
            else:
                majorCalyxNodeIdentifiers.append(None)
                majorCalyxXList.append(None)
                majorCalyxD1List.append(None)

        # Generate top and bottom major calyx to minor calyx connections
        lowerMinor, bottomMinor = 1, 0
        topMinor, upperMinor = 0, 1
        minorCalyxNodeIdentifiers, minorCalyxXList, minorCalyxD1List = [], [], []

        for calyx in [bottomMajor, middleMajor, topMajor]:
            if calyx == middleMajor:
                minorCalyxNodeIdentifiers.append(majorCalyxNodeIdentifiers[1])
                minorCalyxXList.append(majorCalyxXList[1])
                minorCalyxD1List.append(majorCalyxD1List[1])
                continue

            sides = [bottomMinor, lowerMinor] if calyx == bottomMajor else [topMinor, upperMinor]

            if majorCalyxFlags[calyx]:
                cNodeIdentifier = majorCalyxNodeIdentifiers[calyx]
                signValue = -1.0 if calyx == bottomMajor else 1.0
                sx = majorCalyxXList[calyx]

                for side in sides:
                    # Skip if minor calyx is not present
                    if side in [lowerMinor, upperMinor]:
                        if (calyx == bottomMajor and not isLowerMC) or (calyx == topMajor and not isUpperMC):
                            minorCalyxNodeIdentifiers.append(None)
                            minorCalyxXList.append(None)
                            minorCalyxD1List.append(None)
                            continue

                    calyxLength = majorToBottomMinorCalyxLength if side in (
                    bottomMinor, topMinor) else majorToLowerMinorCalyxLength

                    if side in (bottomMinor, topMinor):
                        x = [sx[0], sx[1] + calyxLength * signValue, sx[2]]
                    else:
                        x = [sx[0] + calyxLength, sx[1], sx[2]]

                    # Apply rotation if needed
                    theta = math.radians(-minorCalyxRotateAngle if calyx == bottomMajor else minorCalyxRotateAngle)
                    rx = sx[0] + (x[0] - sx[0]) * math.cos(theta) - (x[1] - sx[1]) * math.sin(theta)
                    ry = sx[1] + (x[0] - sx[0]) * math.sin(theta) + (x[1] - sx[1]) * math.cos(theta)
                    x = [rx, ry, 0.0]

                    # Apply bend angle for lower/upper minor calyxes
                    if side in (lowerMinor, upperMinor):
                        ec = nMinorCalyxesList[1] if side == lowerMinor else nMinorCalyxesList[4]
                        bendAngle = 0 if ec < 2 else minorCalyxBendAngle
                        theta = math.radians(-bendAngle if calyx == bottomMajor else bendAngle)
                        rx = sx[0] + (x[0] - sx[0]) * math.cos(theta) - (x[1] - sx[1]) * math.sin(theta)
                        ry = sx[1] + (x[0] - sx[0]) * math.sin(theta) + (x[1] - sx[1]) * math.cos(theta)
                        x = [rx, ry, 0.0]

                    sd1 = sub(x, sx)
                    minorCalyxNodeIdentifiers.append(nodeIdentifier)
                    minorCalyxXList.append(x)
                    minorCalyxD1List.append(sd1)

                    for i in range(2):
                        isJunctionNode = i == 0
                        nodeId = cNodeIdentifier if isJunctionNode else nodeIdentifier
                        node = nodes.findNodeByIdentifier(nodeId)
                        fieldcache.setNode(node)

                        version = 2 if side == 0 else 3
                        sd3 = [0.0, 0.0, minorCalyxRadius]
                        if calyx == bottomMajor and version == 2:
                            sd2 = [minorCalyxRadius, 0.0, 0.0]
                        else:
                            sd2 = set_magnitude(cross(sd3, sd1), minorCalyxRadius)

                        sid2 = mult(sd2, innerProportionDefault)
                        sid3 = mult(sd3, innerProportionDefault)

                        if not isJunctionNode:
                            for field, derivatives in (
                                    (coordinates, [x, sd1, sd2, sd3]),
                                    (innerCoordinates, [x, sd1, sid2, sid3])
                            ):
                                for valueLabel, value in zip(
                                        (Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                         Node.VALUE_LABEL_D_DS3),
                                        derivatives
                                ):
                                    field.setNodeParameters(fieldcache, -1, valueLabel, 1, value)
                            nodeIdentifier += 1

                        setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, sd3)
                        setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, sid3)
            else:
                for side in sides:
                    minorCalyxNodeIdentifiers.append(None)
                    minorCalyxXList.append(None)
                    minorCalyxD1List.append(None)

        # Generate minor calyx to renal pyramid connections
        anterior, posterior = 0, 1
        bottomMinor, lowerMinor, middleMajor, topMinor, upperMinor = 0, 1, 2, 3, 4
        renalPyramidStartX, renalPyramidNodeIdentifiers = [], []
        connection_sd2_list, connection_sd3_list = [], []

        for calyx in (bottomMinor, lowerMinor, middleMajor, topMinor, upperMinor):
            if minorCalyxNodeIdentifiers[calyx] is None:
                renalPyramidNodeIdentifiers.append(None)
                renalPyramidStartX.append(None)
                connection_sd2_list.append(None)
                connection_sd3_list.append(None)
                continue
            else:
                renalPyramidNodeIdentifiers.append([])
                renalPyramidStartX.append([])
                connection_sd2_list.append([])
                connection_sd3_list.append([])

            index = 0 if calyx <= lowerMinor else (1 if calyx == middleMajor else 2)
            if [isBottomMC, isMidMC, isTopMC][index]:
                cNodeIdentifier = minorCalyxNodeIdentifiers[calyx]
                sx = minorCalyxXList[calyx]
                minorCalyxLength = bottomMinorCalyxLength if calyx in (bottomMinor, topMinor) else lowerMinorCalyxLength

                if nMinorCalyxesList[calyx] > 1:
                    minorCalyxAngle = bottomMinorCalyxAngle if calyx in (bottomMinor, topMinor) else lowerMinorCalyxAngle
                else:
                    minorCalyxAngle = 0

                minorCalyxHalfAngle = 0.5 * minorCalyxAngle
                minorCalyxHalfAngleRadians = math.radians(minorCalyxHalfAngle)
                sinMinorCalyxAngle = math.sin(minorCalyxHalfAngleRadians)
                cosMinorCalyxAngle = math.cos(minorCalyxHalfAngleRadians)

                for side in range(nMinorCalyxesList[calyx]):
                    if calyx in [bottomMinor, topMinor]:
                        nx = sx[0] + minorCalyxLength * (-sinMinorCalyxAngle if side == anterior else sinMinorCalyxAngle)
                        ny = sx[1] + minorCalyxLength * (-cosMinorCalyxAngle if calyx == bottomMinor else cosMinorCalyxAngle)
                        x = [nx, ny, 0.0]
                    elif calyx == lowerMinor and isLowerMC or calyx == upperMinor and isUpperMC:
                        nx = sx[0] + minorCalyxLength * cosMinorCalyxAngle
                        nz = sx[2] + minorCalyxLength * (sinMinorCalyxAngle if side == anterior else -sinMinorCalyxAngle)
                        x = [nx, sx[1], nz]
                        if isRotateMinorCalyx:
                            rotateAxis = [1, 0, 0]
                            tx = rotate_vector_around_vector((sub(x, sx)), rotateAxis, math.radians(90))
                            x = add(tx, sx)
                    else:
                        nx = sx[0] + minorCalyxLength * cosMinorCalyxAngle
                        nz = sx[2] + minorCalyxLength * (sinMinorCalyxAngle if side == anterior else -sinMinorCalyxAngle)
                        x = [nx, sx[1], nz]
                        if isRotateMinorCalyx:
                            rotateAxis = [1, 0, 0]
                            tx = rotate_vector_around_vector((sub(x, sx)), rotateAxis, math.radians(90))
                            x = add(tx, sx)

                    if calyx != middleMajor:
                        theta = math.radians(-minorCalyxRotateAngle) if calyx in [lowerMinor, bottomMinor] else math.radians(minorCalyxRotateAngle)
                        rx = sx[0] + (x[0] - sx[0]) * math.cos(theta) - (x[1] - sx[1]) * math.sin(theta)
                        ry = sx[1] + (x[0] - sx[0]) * math.sin(theta) + (x[1] - sx[1]) * math.cos(theta)
                        x = [rx, ry, x[2]]

                    renalPyramidNodeIdentifiers[-1].append(nodeIdentifier)
                    renalPyramidStartX[-1].append(x)
                    sd1 = sub(x, sx)
                    if calyx in [bottomMinor, topMinor]:
                        sd3 = [0.0, 0.0, minorCalyxRadius]
                        sd2 = set_magnitude(cross(sd3, sd1), minorCalyxRadius)
                    else:
                        sd2 = set_magnitude(cross([0.0, 0.0, 1.0], sd1), minorCalyxRadius)
                        sd3 = set_magnitude(cross(sd1, sd2), minorCalyxRadius)
                    connection_sd2_list[-1].append(sd2)
                    connection_sd3_list[-1].append(sd3)

                    sid2 = mult(sd2, innerProportionDefault)
                    sid3 = mult(sd3, innerProportionDefault)
                    sNodeIdentifiers = []
                    version = 2 if side == anterior else 3
                    for i in range(2):
                        isJunctionNode = i == 0
                        nodeId = cNodeIdentifier if isJunctionNode else nodeIdentifier
                        sNodeIdentifiers.append(nodeId)
                        node = nodes.findNodeByIdentifier(nodeId)
                        fieldcache.setNode(node)
                        if not isJunctionNode:
                            for field, derivatives in (
                                    (coordinates, [x, sd1, sd2, sd3]),
                                    (innerCoordinates, [x, sd1, sid2, sid3])
                            ):
                                for valueLabel, value in zip(
                                        (Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                         Node.VALUE_LABEL_D_DS3),
                                        derivatives
                                ):
                                    field.setNodeParameters(fieldcache, -1, valueLabel, 1, value)
                            nodeIdentifier += 1
                        setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, sd3)
                        setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, sid3)
            else:
                renalPyramidNodeIdentifiers.append([])
                renalPyramidStartX.append([])
                connection_sd2_list.append([])
                connection_sd3_list.append([])
                for side in range(nMinorCalyxesList[calyx]):
                    renalPyramidNodeIdentifiers[-1].append(None)
                    renalPyramidStartX[-1].append(None)
                    connection_sd2_list[-1].append(None)
                    connection_sd3_list[-1].append(None)

        # Generate minor calyx to renal pyramids
        pyramidHalfWidth = 0.5 * pyramidWidth

        for calyx in [bottomMinor, lowerMinor, middleMajor, topMinor, upperMinor]:
            if renalPyramidNodeIdentifiers[calyx] is None:
                continue

            index = 0 if calyx <= lowerMinor else (1 if calyx == middleMajor else 2)
            if majorCalyxFlags[index]:
                for side in range(nMinorCalyxesList[calyx]):
                    sx = renalPyramidStartX[calyx][side]

                    # Calculate pyramid direction
                    if calyx in [bottomMinor, topMinor]:
                        tx = [0.0, -1.0, 0.0] if calyx == bottomMinor else [0.0, 1.0, 0.0]
                        angle = 25 if nMinorCalyxesList[calyx + 1] > 1 else bottomMinorCalyxAngle * 0.5
                        rotateAngle = 0 if nMinorCalyxesList[calyx] <= 1 else (-angle if side == anterior else angle)

                        rotateAngle = (-rotateAngle + minorCalyxRotateAngle if calyx == topMinor else rotateAngle - minorCalyxRotateAngle)
                        theta = math.radians(rotateAngle)
                        rx = tx[0] * math.cos(theta) - tx[1] * math.sin(theta)
                        ry = tx[0] * math.sin(theta) + tx[1] * math.cos(theta)
                        pyramidDirn = [rx, ry, 0.0]
                    else:
                        tx = [1.0, 0.0, 0.0]
                        theta = 0
                        if nMinorCalyxesList[calyx] > 1:
                            if calyx == lowerMinor:
                                theta = math.radians(-minorCalyxBendAngle)
                            elif calyx == upperMinor:
                                theta = math.radians(minorCalyxBendAngle)

                        rx = tx[0] * math.cos(theta) - tx[1] * math.sin(theta)
                        ry = tx[0] * math.sin(theta) + tx[1] * math.cos(theta)
                        pyramidDirn = [rx, ry, 0.0]

                        if isRotateMinorCalyx:
                            rotateAngle = 0 if nMinorCalyxesList[calyx] <= 1 else (-25 if side == anterior else 25)
                            if calyx in [lowerMinor, upperMinor]:
                                rotateAngle += (minorCalyxRotateAngle if calyx == upperMinor else -minorCalyxRotateAngle)
                            rotateAngleRad = math.radians(rotateAngle)
                            pyramidDirn = rotate_about_z_axis(pyramidDirn, rotateAngleRad)

                    # Generate pyramid node positions
                    xList = []
                    for e in range(pyramidElementsCount):
                        pyramidLengthScale = [0.3, 0.7, 1.0][e] * pyramidLength
                        x = add(sx, mult(pyramidDirn, pyramidLengthScale))
                        xList.append(x)

                    # Generate pyramid nodes
                    pyramid_sd2_list = []
                    pyramid_sd3_list = []
                    sNodeIdentifiers = []

                    for e in range(pyramidElementsCount):
                        sNodeIdentifiers.append(nodeIdentifier)
                        node = nodes.findNodeByIdentifier(nodeIdentifier)
                        fieldcache.setNode(node)

                        sd1 = sub(xList[1], xList[0])
                        pyramidWidthScale = [minorCalyxRadius, pyramidHalfWidth, 0.9 * pyramidHalfWidth][e]
                        thickness = minorCalyxRadius if e == 0 else (0.9 * pyramidThickness if e == pyramidElementsCount - 1 else pyramidThickness)

                        sd3 = [0.0, 0.0, thickness]
                        sd2 = set_magnitude(cross(sd3, sd1), pyramidWidthScale)

                        pyramid_sd2_list.append(sd2)
                        pyramid_sd3_list.append(sd3)
                        sid2 = mult(sd2, innerProportionDefault)
                        sid3 = mult(sd3, innerProportionDefault)

                        for field, derivatives in (
                                (coordinates, [xList[e], sd1, sd2, sd3]),
                                (innerCoordinates, [xList[e], sd1, sid2, sid3])
                        ):
                            for valueLabel, value in zip(
                                    (Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                     Node.VALUE_LABEL_D_DS3),
                                    derivatives
                            ):
                                field.setNodeParameters(fieldcache, -1, valueLabel, 1, value)
                        nodeIdentifier += 1

                    # Set cross derivatives for pyramid nodes
                    pyramid_sd2_list.append(connection_sd2_list[calyx][side])
                    pyramid_sd3_list.append(connection_sd3_list[calyx][side])

                    for e in range(pyramidElementsCount):
                        node = nodes.findNodeByIdentifier(sNodeIdentifiers[e])
                        fieldcache.setNode(node)
                        sd12 = sub(pyramid_sd2_list[e + 1], pyramid_sd2_list[e])
                        sd13 = sub(pyramid_sd3_list[e + 1], pyramid_sd3_list[e])
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13)
                        sid12 = mult(sd12, innerProportionDefault)
                        sid13 = mult(sd13, innerProportionDefault)
                        innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sid12)
                        innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sid13)

        return annotationGroups, networkMesh

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Edit base class list to include only valid functions.
        """
        interactiveFunctions = super(MeshType_1d_renal_pelvis_network_layout1, cls).getInteractiveFunctions()
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Edit structure...":
                interactiveFunctions.remove(interactiveFunction)
                break
        return interactiveFunctions


class MeshType_3d_renal_pelvis1(Scaffold_base):
    """
    Generates a 3-D renal pelvis.
    """

    @classmethod
    def getName(cls):
        return "3D Renal Pelvis 1"

    @classmethod
    def getParameterSetNames(cls):
        return MeshType_1d_renal_pelvis_network_layout1.getParameterSetNames()

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Base parameter set"] = useParameterSetName
        options["Network layout"] = ScaffoldPackage(MeshType_1d_renal_pelvis_network_layout1,
                                                    defaultParameterSetName=useParameterSetName)
        options["Left"] = True
        options["Elements count around"] = 8
        options["Elements count through shell"] = 1
        options["Annotation elements counts around"] = [0]
        options["Target element density along longest segment"] = 4.0
        options["Use linear through shell"] = False
        options["Use outer trim surfaces"] = True
        options["Show trim surfaces"] = False
        options["Refine"] = False
        options["Refine number of elements"] = 4
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Network layout",
            "Left",
            "Elements count around",
            "Elements count through shell",
            "Annotation elements counts around",
            "Target element density along longest segment",
            "Use linear through shell",
            "Use outer trim surfaces",
            "Show trim surfaces",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Network layout":
            return [MeshType_1d_renal_pelvis_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + ".getOptionScaffoldTypeParameterSetNames.  " + \
            "Invalid option \"" + optionName + "\" scaffold type " + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()  # use the defaults from the network layout scaffold

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
        if optionName == "Network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_renal_pelvis_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if (options["Network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Network layout")):
            options["Body network layout"] = ScaffoldPackage(MeshType_1d_renal_pelvis_network_layout1)

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        elementsCountAlongPyramid = 4

        annotationAlongCounts = []
        defaultCoreBoundaryScalingMode = 1
        annotationCoreBoundaryScalingMode = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            alongCount = 0
            coreBoundaryScalingMode = 0
            name = layoutAnnotationGroup.getName()
            if "renal pyramid" in name:
                alongCount = elementsCountAlongPyramid
            annotationAlongCounts.append(alongCount)
            annotationCoreBoundaryScalingMode.append(coreBoundaryScalingMode)

        tubeNetworkMeshBuilder = RenalPelvisTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAlong=annotationAlongCounts,
            defaultElementsCountAround=options["Elements count around"],
            elementsCountThroughShell=options["Elements count through shell"],
            annotationElementsCountsAround=options["Annotation elements counts around"],
            defaultCoreBoundaryScalingMode=defaultCoreBoundaryScalingMode,
            annotationCoreBoundaryScalingMode=annotationCoreBoundaryScalingMode)

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughShell=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        mesh = generateData.getMesh()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        coreGroup = getAnnotationGroupForTerm(annotationGroups, ("core", "")).getGroup()
        shellGroup = getAnnotationGroupForTerm(annotationGroups, ("shell", "")).getGroup()

        tempGroup = fm.createFieldAdd(shellGroup, coreGroup)
        allGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,("all", ""))
        allGroup.getMeshGroup(mesh).addElementsConditional(tempGroup)
        allNodeset = allGroup.getNodesetGroup(nodes)

        isLeft = options["Left"]
        rotateMeshAboutAxis(90, fm, coordinates, allNodeset, axis=3)
        if not isLeft:
            rotateMeshAboutAxis(180, fm, coordinates, allNodeset, axis=1)

        return annotationGroups, None

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCount = options['Refine number of elements']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)


def rotateMeshAboutAxis(rotateAngle, fm, coordinates, nodeset, axis):
    """
    Rotates the mesh coordinates about a specified axis using the right-hand rule.
    :param rotateAngle: Angle of rotation in degrees.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param nodeset: Zinc NodesetGroup containing nodes to transform.
    :param axis: Axis of rotation.
    :return: None
    """
    rotateAngle = -math.radians(rotateAngle)  # negative value due to right handed rule

    if axis == 1:
        # Rotation about x-axis
        rotateMatrix = fm.createFieldConstant([1.0, 0.0, 0.0,
                                               0.0, math.cos(rotateAngle), math.sin(rotateAngle),
                                               0.0, -math.sin(rotateAngle), math.cos(rotateAngle)])
    elif axis == 2:
        # Rotation about y-axis
        rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), 0.0, -math.sin(rotateAngle),
                                               0.0, 1.0, 0.0,
                                               math.sin(rotateAngle), 0.0, math.cos(rotateAngle)])
    elif axis == 3:
        # Rotation about z-axis
        rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), math.sin(rotateAngle), 0.0,
                                               -math.sin(rotateAngle), math.cos(rotateAngle), 0.0,
                                               0.0, 0.0, 1.0])

    rotated_coordinates = fm.createFieldMatrixMultiply(3, rotateMatrix, coordinates)

    fieldassignment = coordinates.createFieldassignment(rotated_coordinates)
    fieldassignment.setNodeset(nodeset)
    fieldassignment.assign()


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
