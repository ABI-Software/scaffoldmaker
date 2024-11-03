"""
Generates a 3D body coordinates using tube network mesh.
"""
from cmlibs.maths.vectorops import add, cross, mult, set_magnitude, sub
from cmlibs.utils.zinc.field import Field, find_or_create_field_coordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import (
    AnnotationGroup, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm)
from scaffoldmaker.annotation.trigeminal_nerve_terms import get_trigeminal_nerve_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.interpolation import (
    computeCubicHermiteEndDerivative, getCubicHermiteArcLength, interpolateLagrangeHermiteDerivative,
    sampleCubicHermiteCurvesSmooth, smoothCubicHermiteDerivativesLine)
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import BodyTubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
import math


class MeshType_1d_human_trigeminal_nerve_network_layout1(MeshType_1d_network_layout1):
    """
    Defines trigeminal nerve network layout.
    """

    @classmethod
    def getName(cls):
        return "1D Human Trigeminal Nerve Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = parameterSetName
        options["Structure"] = (
            "1-2-3,"
            "3-4-5.1," 
            "5.2-6-7-8-9,"
            "5.3-10-11-12-13,"
            "5.4-14-15-16-17")
        options["Define inner coordinates"] = True
        options["Trigeminal nerve root depth"] = 0.5
        options["Trigeminal nerve root length"] = 3.0
        options["Trigeminal nerve root width"] = 1.0
        options["Trigeminal ganglion depth"] = 0.5
        options["Trigeminal ganglion length"] = 2.0
        options["Trigeminal ganglion width"] = 3.5
        options["Mandibular nerve depth"] = 0.5
        options["Mandibular nerve lateral angle degrees"] = -45
        options["Mandibular nerve length"] = 3.0
        options["Mandibular nerve start width"] = 2.0
        options["Mandibular nerve end width"] = 1.0
        options["Maxillary nerve depth"] = 0.5
        options["Maxillary nerve lateral angle degrees"] = 0.0
        options["Maxillary nerve length"] = 3.0
        options["Maxillary nerve start width"] = 2.0
        options["Maxillary nerve end width"] = 1.0
        options["Ophthalmic nerve depth"] = 0.5
        options["Ophthalmic nerve lateral angle degrees"] = 45
        options["Ophthalmic nerve length"] = 2.0
        options["Ophthalmic nerve start width"] = 2.0
        options["Ophthalmic nerve end width"] = 1.0
        options["Inner proportion default"] = 0.75

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Trigeminal nerve root depth",
            "Trigeminal nerve root length",
            "Trigeminal nerve root width",
            "Trigeminal ganglion depth",
            "Trigeminal ganglion length",
            "Trigeminal ganglion width",
            "Mandibular nerve depth",
            "Mandibular nerve lateral angle degrees",
            "Mandibular nerve length",
            "Mandibular nerve start width",
            "Mandibular nerve end width",
            "Maxillary nerve depth",
            "Maxillary nerve lateral angle degrees",
            "Maxillary nerve length",
            "Maxillary nerve start width",
            "Maxillary nerve end width",
            "Ophthalmic nerve depth",
            "Ophthalmic nerve lateral angle degrees",
            "Ophthalmic nerve length",
            "Ophthalmic nerve start width",
            "Ophthalmic nerve end width",
            "Inner proportion default"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Trigeminal nerve root depth",
            "Trigeminal nerve root length",
            "Trigeminal nerve root width",
            "Trigeminal ganglion depth",
            "Trigeminal ganglion length",
            "Trigeminal ganglion width",
            "Mandibular nerve depth",
            "Mandibular nerve length",
            "Mandibular nerve start width",
            "Mandibular nerve end width",
            "Maxillary nerve depth",
            "Maxillary nerve length",
            "Maxillary nerve start width",
            "Maxillary nerve end width",
            "Ophthalmic nerve depth",
            "Ophthalmic nerve length",
            "Ophthalmic nerve start width",
            "Ophthalmic nerve end width"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
        for key in [
            "Inner proportion default",
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 0.9:
                options[key] = 0.9
        # for key in [
        #     "Mandibular nerve lateral angle degrees",
        #     "Maxillary nerve lateral angle degrees",
        #     "Ophthalmic nerve lateral angle degrees"
        # ]:
        #     if options[key] < -60.0:
        #         options[key] = -60.0
        #     elif options[key] > 200.0:
        #         options[key] = 200.0
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
        halfNerveRootDepth = 0.5 * options["Trigeminal nerve root depth"]
        nerveRootLength = options["Trigeminal nerve root length"]
        halfNerveRootWidth = 0.5 * options["Trigeminal nerve root width"]

        halfTrigeminalGanglionDepth = 0.5 * options["Trigeminal ganglion depth"]
        trigeminalGanglionLength = options["Trigeminal ganglion length"]
        halfTrigeminalGanglionWidth = 0.5 * options["Trigeminal ganglion width"]

        halfMandibularDepth = 0.5 * options["Mandibular nerve depth"]
        mandibularLength = options["Mandibular nerve length"]
        halfMandibularStartWidth = 0.5 * options["Mandibular nerve start width"]
        halfMandibularEndWidth = 0.5 * options["Mandibular nerve end width"]
        mandibularAngleRadians = math.radians(options["Mandibular nerve lateral angle degrees"])

        maxillaryLength = options["Maxillary nerve length"]
        halfMaxillaryStartWidth = 0.5 * options["Maxillary nerve start width"]
        halfMaxillaryEndWidth = 0.5 * options["Maxillary nerve end width"]
        halfMaxillaryDepth = 0.5 * options["Maxillary nerve depth"]
        maxillaryAngleRadians = math.radians(options["Maxillary nerve lateral angle degrees"])

        ophthalmicLength = options["Ophthalmic nerve length"]
        halfOphthalmicStartWidth = 0.5 * options["Ophthalmic nerve start width"]
        halfOphthalmicEndWidth = 0.5 * options["Ophthalmic nerve end width"]
        halfOphthalmicDepth = 0.5 * options["Ophthalmic nerve depth"]
        ophthalmicAngleRadians = math.radians(options["Ophthalmic nerve lateral angle degrees"])

        innerProportionDefault = options["Inner proportion default"]

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)

        # set up element annotations
        trigeminalNerveGroup = AnnotationGroup(region, get_trigeminal_nerve_term("trigeminal nerve"))
        nerveRootGroup = AnnotationGroup(region, get_trigeminal_nerve_term("trigeminal nerve root"))
        trigeminalGanglionGroup = AnnotationGroup(region, get_trigeminal_nerve_term("trigeminal ganglion"))
        mandibularNerveGroup = AnnotationGroup(region, get_trigeminal_nerve_term("mandibular nerve"))
        maxillaryNerveGroup = AnnotationGroup(region, get_trigeminal_nerve_term("maxillary nerve"))
        ophthalmicNerveGroup = AnnotationGroup(region, get_trigeminal_nerve_term("ophthalmic nerve"))
        annotationGroups = [trigeminalNerveGroup, nerveRootGroup, trigeminalGanglionGroup,
                            mandibularNerveGroup, maxillaryNerveGroup, ophthalmicNerveGroup]
        trigeminalNerveMeshGroup = trigeminalNerveGroup.getMeshGroup(mesh)
        elementIdentifier = 1
        nerveRootElementsCount = 2
        meshGroups = [trigeminalNerveMeshGroup, nerveRootGroup.getMeshGroup(mesh)]
        for e in range(nerveRootElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1
        trigeminalGanglionElementsCount = 2
        meshGroups = [trigeminalNerveMeshGroup, trigeminalGanglionGroup.getMeshGroup(mesh)]
        for e in range(trigeminalGanglionElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        nerveBranchElementsCount = 4
        for nerveBranch in ('mandibular','maxillary','ophthalmic'):
            if nerveBranch == 'mandibular':
                nerveGroup = mandibularNerveGroup
            elif nerveBranch == 'maxillary':
                nerveGroup = maxillaryNerveGroup
            elif nerveBranch == 'ophthalmic':
                nerveGroup = ophthalmicNerveGroup
            meshGroups = [trigeminalNerveMeshGroup, nerveGroup.getMeshGroup(mesh)]
            for e in range(nerveBranchElementsCount):
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

        nerveRootScale = nerveRootLength / nerveRootElementsCount
        nodeIdentifier = 1
        d1 = [nerveRootScale, 0.0, 0.0]
        d2 = [0.0, halfNerveRootWidth, 0.0]
        d3 = [0.0, 0.0, halfNerveRootDepth]
        id2 = mult(d2, innerProportionDefault)
        id3 = mult(d3, innerProportionDefault)
        for i in range(nerveRootElementsCount):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [nerveRootScale * i, 0.0, 0.0]
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3)
            print('root nodes', nodeIdentifier)
            nodeIdentifier += 1

        trigeminalGanglionScale = trigeminalGanglionLength / trigeminalGanglionElementsCount
        # halfTrigeminalGanglionWidth = 0.5 * (halfMaxillaryStartWidth + halfMandibularStartWidth)
        d12_mag = (halfTrigeminalGanglionWidth - halfNerveRootWidth) / trigeminalGanglionElementsCount
        d3 = [0.0, 0.0, halfTrigeminalGanglionDepth]
        id3 = mult(d3, innerProportionDefault)
        for i in range(trigeminalGanglionElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [nerveRootLength + trigeminalGanglionScale * i, 0.0, 0.0]
            d1 = [0.5 * (nerveRootScale + trigeminalGanglionScale) if (i == 0) else trigeminalGanglionScale, 0.0, 0.0]
            sHalfWidth = halfNerveRootWidth + i * trigeminalGanglionScale * d12_mag
            d2 = [0.0, sHalfWidth, 0.0]
            d12 = [0.0, d12_mag, 0.0]
            id2 = mult(d2, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
            print('TG nodes', nodeIdentifier)
            nodeIdentifier += 1

        nerveJunctionNodeIdentifier = nodeIdentifier - 1
        xJunction = x
        # print('Junction', nerveJunctionNodeIdentifier)

        m1ShoulderDrop = 0.7
        m1ShoulderRotationFactor = 1.0 - math.cos(0.5 * mandibularAngleRadians)
        # assume shoulder drop is half shrug distance to get limiting shoulder angle for 180 degree arm rotation
        m1ShoulderLimitAngleRadians = math.asin(2.0 * m1ShoulderDrop / halfTrigeminalGanglionWidth)
        m1ShoulderAngleRadians = m1ShoulderRotationFactor * m1ShoulderLimitAngleRadians
        mandibularStartX = nerveRootLength + trigeminalGanglionLength + m1ShoulderDrop - \
                           halfTrigeminalGanglionWidth * math.sin(m1ShoulderAngleRadians)
        mandibularStartY = -halfTrigeminalGanglionWidth * math.cos(m1ShoulderAngleRadians)
        x = [mandibularStartX, mandibularStartY, 0.0]

        mandibularScale = mandibularLength / (nerveBranchElementsCount - 2)
        cosMandibularAngle = math.cos(mandibularAngleRadians)
        sinMandibularAngle = math.sin(mandibularAngleRadians)
        d1 = [mandibularScale * cosMandibularAngle, mandibularScale * sinMandibularAngle, 0.0]

        sd1 = interpolateLagrangeHermiteDerivative(xJunction, x, d1, 0.0)
        nx, nd1 = sampleCubicHermiteCurvesSmooth([xJunction, x], [sd1, d1], 2, derivativeMagnitudeEnd=mandibularScale)[
                  0:2]
        arcLengths = [getCubicHermiteArcLength(nx[i], nd1[i], nx[i + 1], nd1[i + 1]) for i in range(2)]

        d12_mag = (halfMandibularEndWidth - halfMandibularStartWidth) / (nerveBranchElementsCount - 2)
        d3 = [0.0, 0.0, halfMandibularDepth]
        id3 = mult(d3, innerProportionDefault)

        node = nodes.findNodeByIdentifier(nodeIdentifier)
        fieldcache.setNode(node)
        sd2_list = []

        sNodeIdentifiers = []
        for i in range(2):
            sNodeIdentifiers.append(nodeIdentifier if (i > 0) else nerveJunctionNodeIdentifier)
            node = nodes.findNodeByIdentifier(sNodeIdentifiers[-1])
            fieldcache.setNode(node)
            version = 1 if (i > 0) else 2
            sd1 = nd1[i]
            sDistance = sum(arcLengths[i:])
            sHalfWidth = halfMandibularStartWidth + sDistance * -d12_mag
            sd2 = set_magnitude(cross(d3, sd1), sHalfWidth)
            sid2 = mult(sd2, innerProportionDefault)
            sd2_list.append(sd2)
            if i > 0:
                for field in (coordinates, innerCoordinates):
                    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, nx[i])
                print('mandibular nodes A', nodeIdentifier)
                nodeIdentifier += 1
            setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, d3)
            setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, id3)
        sd2_list.append([-halfOphthalmicStartWidth * sinMandibularAngle,
                         halfOphthalmicStartWidth * cosMandibularAngle, 0.0])

        for i in range(2):
            node = nodes.findNodeByIdentifier(sNodeIdentifiers[i])
            fieldcache.setNode(node)
            version = 1 if (i > 0) else 2
            sd12 = sub(sd2_list[i + 1], sd2_list[i])
            # sd13 = sub(sd3_list[i + 1], sd3_list[i])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sd12)
            # coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sd13)
            sid12 = mult(sd12, innerProportionDefault)
            # sid13 = mult(sd13, innerProportionDefault)
            innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sid12)
            # innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sid13)

        # remaining part of nerve
        for i in range(nerveBranchElementsCount - 1):
            # xi = i / (nerveBranchElementsCount - 2)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [mandibularStartX + d1[0] * i, mandibularStartY + d1[1] * i, d1[2] * i]
            sHalfWidth = halfMandibularStartWidth + i * mandibularScale * d12_mag
            d2 = set_magnitude(cross(d3, d1), sHalfWidth)
            d12 = [-d12_mag * sinMandibularAngle, d12_mag * cosMandibularAngle, 0.0]
            id2 = mult(d2, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
            print('mandibular nodes', nodeIdentifier)
            nodeIdentifier += 1

        # mandibularScale = mandibularLength / nerveBranchElementsCount
        # mandibularStartX = nerveRootLength + trigeminalGanglionLength
        # sx = [mandibularStartX, 0.0, 0.0]
        # cosMandibularAngle = math.cos(mandibularAngleRadians)
        # sinMandibularAngle = math.sin(mandibularAngleRadians)
        # d12_mag = (halfMandibularEndWidth - halfMandibularStartWidth) / nerveBranchElementsCount
        # for i in range(nerveBranchElementsCount + 1):
        #     node = nodes.findNodeByIdentifier(nodeIdentifier)
        #     fieldcache.setNode(node)
        #     x = [sx[0] + mandibularScale * i * cosMandibularAngle,
        #          sx[1] + mandibularScale * i * sinMandibularAngle,
        #          0.0]
        #     d1 = [mandibularScale*cosMandibularAngle, mandibularScale*sinMandibularAngle, 0.0]
        #     d3 = [0.0, 0.0, halfMandibularDepth]
        #     id3 = mult(d3, innerProportionDefault)
        #     sHalfWidth = halfMandibularStartWidth + i * mandibularScale * d12_mag
        #     d2 = set_magnitude(cross(d3, d1), sHalfWidth)
        #     d12 = [-d12_mag * sinMandibularAngle, d12_mag * cosMandibularAngle, 0.0]
        #     id2 = mult(d2, innerProportionDefault)
        #     id12 = mult(d12, innerProportionDefault)
        #
        #     if i == 0:
        #         node = nodes.findNodeByIdentifier(nerveJunctionNodeIdentifier)
        #         fieldcache.setNode(node)
        #         version = 2
        #         setNodeFieldVersionDerivatives(coordinates, fieldcache, version, d1, d2, d3)
        #         setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, d1, id2, id3)
        #     else:
        #         setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
        #         setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
        #         nodeIdentifier += 1
        #     print('mandibular nodes', nodeIdentifier)

        maxillaryScale = maxillaryLength / nerveBranchElementsCount
        maxillaryStartX = nerveRootLength + trigeminalGanglionLength
        sx = [maxillaryStartX, 0.0, 0.0]
        cosMaxillaryAngle = math.cos(maxillaryAngleRadians)
        sinMaxillaryAngle = math.sin(maxillaryAngleRadians)
        d12_mag = (halfMaxillaryEndWidth - halfMaxillaryStartWidth) / nerveBranchElementsCount
        for i in range(nerveBranchElementsCount + 1):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [sx[0] + maxillaryScale * i * cosMaxillaryAngle,
                 sx[1] + maxillaryScale * i * sinMaxillaryAngle,
                 0.0]
            d1 = [maxillaryScale * cosMaxillaryAngle, maxillaryScale * sinMaxillaryAngle, 0.0]
            d3 = [0.0, 0.0, halfMaxillaryDepth]
            id3 = mult(d3, innerProportionDefault)
            sHalfWidth = halfMaxillaryStartWidth + i * maxillaryScale * d12_mag
            d2 = set_magnitude(cross(d3, d1), sHalfWidth)
            d12 = [-d12_mag * sinMaxillaryAngle, d12_mag * cosMaxillaryAngle, 0.0]
            id2 = mult(d2, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            if i == 0:
                node = nodes.findNodeByIdentifier(nerveJunctionNodeIdentifier)
                fieldcache.setNode(node)
                version = 3
                setNodeFieldVersionDerivatives(coordinates, fieldcache, version, d1, d2, d3)
                setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, d1, id2, id3)
            else:
                setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
                setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
                nodeIdentifier += 1
            print('maxillary nodes', nodeIdentifier)

        oShoulderDrop = 0.7
        oShoulderRotationFactor = 1.0 - math.cos(0.5 * ophthalmicAngleRadians)
        # assume shoulder drop is half shrug distance to get limiting shoulder angle for 180 degree arm rotation
        oShoulderLimitAngleRadians = math.asin(2.0 * oShoulderDrop / halfTrigeminalGanglionWidth)
        oShoulderAngleRadians = oShoulderRotationFactor * oShoulderLimitAngleRadians
        ophthalmicStartX = nerveRootLength + trigeminalGanglionLength + oShoulderDrop - \
                           halfTrigeminalGanglionWidth * math.sin(oShoulderAngleRadians)
        ophthalmicStartY = halfTrigeminalGanglionWidth * math.cos(oShoulderAngleRadians)
        x = [ophthalmicStartX, ophthalmicStartY, 0.0]

        ophthalmicScale = ophthalmicLength / (nerveBranchElementsCount - 2)
        cosOphthalmicAngle = math.cos(ophthalmicAngleRadians)
        sinOphthalmicAngle = math.sin(ophthalmicAngleRadians)
        d1 = [ophthalmicScale * cosOphthalmicAngle, ophthalmicScale * sinOphthalmicAngle, 0.0]

        sd1 = interpolateLagrangeHermiteDerivative(xJunction, x, d1, 0.0)
        nx, nd1 = sampleCubicHermiteCurvesSmooth([xJunction, x], [sd1, d1], 2, derivativeMagnitudeEnd=ophthalmicScale)[
                  0:2]
        arcLengths = [getCubicHermiteArcLength(nx[i], nd1[i], nx[i + 1], nd1[i + 1]) for i in range(2)]

        d12_mag = (halfOphthalmicEndWidth - halfOphthalmicStartWidth) / (nerveBranchElementsCount - 2)
        d3 = [0.0, 0.0, halfOphthalmicDepth]
        id3 = mult(d3, innerProportionDefault)

        node = nodes.findNodeByIdentifier(nodeIdentifier)
        fieldcache.setNode(node)
        sd2_list = []

        sNodeIdentifiers = []
        for i in range(2):
            sNodeIdentifiers.append(nodeIdentifier if (i > 0) else nerveJunctionNodeIdentifier)
            node = nodes.findNodeByIdentifier(sNodeIdentifiers[-1])
            fieldcache.setNode(node)
            version = 1 if (i > 0) else 4
            sd1 = nd1[i]
            sDistance = sum(arcLengths[i:])
            sHalfWidth = halfOphthalmicStartWidth + sDistance * -d12_mag
            sd2 = set_magnitude(cross(d3, sd1), sHalfWidth)
            sid2 = mult(sd2, innerProportionDefault)
            sd2_list.append(sd2)
            if i > 0:
                for field in (coordinates, innerCoordinates):
                    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, nx[i])
                print('ophthalmic nodes A', nodeIdentifier)
                nodeIdentifier += 1
            setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, d3)
            setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, id3)
        sd2_list.append([-halfOphthalmicStartWidth * sinOphthalmicAngle,
                         halfOphthalmicStartWidth * cosOphthalmicAngle, 0.0])

        for i in range(2):
            node = nodes.findNodeByIdentifier(sNodeIdentifiers[i])
            fieldcache.setNode(node)
            version = 1 if (i > 0) else 4
            sd12 = sub(sd2_list[i + 1], sd2_list[i])
            # sd13 = sub(sd3_list[i + 1], sd3_list[i])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sd12)
            # coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sd13)
            sid12 = mult(sd12, innerProportionDefault)
            # sid13 = mult(sd13, innerProportionDefault)
            innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sid12)
            # innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sid13)

        # remaining part of nerve
        for i in range(nerveBranchElementsCount - 1):
            # xi = i / (nerveBranchElementsCount - 2)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            x = [ophthalmicStartX + d1[0] * i, ophthalmicStartY + d1[1] * i, d1[2] * i]
            sHalfWidth = halfOphthalmicStartWidth + i * ophthalmicScale * d12_mag
            d2 = set_magnitude(cross(d3, d1), sHalfWidth)
            d12 = [-d12_mag * sinOphthalmicAngle, d12_mag * cosOphthalmicAngle, 0.0]
            id2 = mult(d2, innerProportionDefault)
            id12 = mult(d12, innerProportionDefault)
            setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12)
            setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12)
            print('ophthalmic nodes', nodeIdentifier)
            nodeIdentifier += 1


        # # nerve branches
        # mandibularStartX = thoraxStartX + shoulderDrop - halfShoulderWidth * math.sin(shoulderAngleRadians)
        # nonHandArmLength = armLength - handLength
        # armScale = nonHandArmLength / (armToHandElementsCount - 2)  # 2 == shoulder elements count
        # d12_mag = (halfWristThickness - armTopRadius) / (armToHandElementsCount - 2)
        # d13_mag = (halfWristWidth - armTopRadius) / (armToHandElementsCount - 2)
        # hd3 = [0.0, 0.0, halfHandWidth]
        # hid3 = mult(hd3, innerProportionDefault)
        # for side in (left, right):
        #     armAngle = armAngleRadians if (side == left) else -armAngleRadians
        #     cosArmAngle = math.cos(armAngle)
        #     sinArmAngle = math.sin(armAngle)
        #     armStartY = (halfShoulderWidth if (side == left) else -halfShoulderWidth) * math.cos(shoulderAngleRadians)
        #     x = [armStartX, armStartY, 0.0]
        #     d1 = [armScale * cosArmAngle, armScale * sinArmAngle, 0.0]
        #     # set leg versions 2 (left) and 3 (right) on leg junction node, and intermediate shoulder node
        #     sd1 = interpolateLagrangeHermiteDerivative(sx, x, d1, 0.0)
        #     nx, nd1 = sampleCubicHermiteCurvesSmooth([sx, x], [sd1, d1], 2, derivativeMagnitudeEnd=armScale)[0:2]
        #     arcLengths = [getCubicHermiteArcLength(nx[i], nd1[i], nx[i + 1], nd1[i + 1]) for i in range(2)]
        #     sd2_list = []
        #     sd3_list = []
        #     sNodeIdentifiers = []
        #     for i in range(2):
        #         sNodeIdentifiers.append(nodeIdentifier if (i > 0) else armJunctionNodeIdentifier)
        #         node = nodes.findNodeByIdentifier(sNodeIdentifiers[-1])
        #         fieldcache.setNode(node)
        #         version = 1 if (i > 0) else 2 if (side == left) else 3
        #         sd1 = nd1[i]
        #         sDistance = sum(arcLengths[i:])
        #         sHalfHeight = armTopRadius + sDistance * -d12_mag
        #         sHalfDepth = armTopRadius + sDistance * -d13_mag
        #         sd3 = [0.0, 0.0, sHalfDepth]
        #         sid3 = mult(sd3, innerProportionDefault)
        #         sd2 = set_magnitude(cross(sd3, sd1), sHalfHeight)
        #         sid2 = mult(sd2, innerProportionDefault)
        #         sd2_list.append(sd2)
        #         sd3_list.append(sd3)
        #         if i > 0:
        #             for field in (coordinates, innerCoordinates):
        #                 field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, nx[i])
        #             nodeIdentifier += 1
        #         setNodeFieldVersionDerivatives(coordinates, fieldcache, version, sd1, sd2, sd3)
        #         setNodeFieldVersionDerivatives(innerCoordinates, fieldcache, version, sd1, sid2, sid3)
        #     sd2_list.append([-armTopRadius * sinArmAngle, armTopRadius * cosArmAngle, 0.0])
        #     sd3_list.append([0.0, 0.0, armTopRadius])
        #     for i in range(2):
        #         node = nodes.findNodeByIdentifier(sNodeIdentifiers[i])
        #         fieldcache.setNode(node)
        #         version = 1 if (i > 0) else 2 if (side == left) else 3
        #         sd12 = sub(sd2_list[i + 1], sd2_list[i])
        #         sd13 = sub(sd3_list[i + 1], sd3_list[i])
        #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sd12)
        #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sd13)
        #         sid12 = mult(sd12, innerProportionDefault)
        #         sid13 = mult(sd13, innerProportionDefault)
        #         innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sid12)
        #         innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sid13)
        #     d12 = [-d12_mag * sinArmAngle, d12_mag * cosArmAngle, 0.0]
        #     id12 = mult(d12, innerProportionDefault)
        #     d13 = [0.0, 0.0, d13_mag]
        #     id13 = mult(d13, innerProportionDefault)
        #     # main part of arm to wrist
        #     for i in range(armToHandElementsCount - 1):
        #         xi = i / (armToHandElementsCount - 2)
        #         node = nodes.findNodeByIdentifier(nodeIdentifier)
        #         fieldcache.setNode(node)
        #         x = [armStartX + d1[0] * i, armStartY + d1[1] * i, d1[2] * i]
        #         halfThickness = xi * halfWristThickness + (1.0 - xi) * armTopRadius
        #         halfWidth = xi * halfWristWidth + (1.0 - xi) * armTopRadius
        #         d2 = [-halfThickness * sinArmAngle, halfThickness * cosArmAngle, 0.0]
        #         d3 = [0.0, 0.0, halfWidth]
        #         id2 = mult(d2, innerProportionDefault)
        #         id3 = mult(d3, innerProportionDefault)
        #         setNodeFieldParameters(coordinates, fieldcache, x, d1, d2, d3, d12, d13)
        #         setNodeFieldParameters(innerCoordinates, fieldcache, x, d1, id2, id3, id12, id13)
        #         nodeIdentifier += 1
        #     # hand
        #     assert handElementsCount == 1
        #     node = nodes.findNodeByIdentifier(nodeIdentifier)
        #     fieldcache.setNode(node)
        #     hx = [armStartX + armLength * cosArmAngle, armStartY + armLength * sinArmAngle, 0.0]
        #     hd1 = computeCubicHermiteEndDerivative(x, d1, hx, d1)
        #     hd2 = set_magnitude(d2, halfHandThickness)
        #     hid2 = mult(hd2, innerProportionDefault)
        #     setNodeFieldParameters(coordinates, fieldcache, hx, hd1, hd2, hd3)
        #     setNodeFieldParameters(innerCoordinates, fieldcache, hx, hd1, hid2, hid3)
        #     nodeIdentifier += 1

        return annotationGroups, networkMesh

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Edit base class list to include only valid functions.
        """
        interactiveFunctions = super(MeshType_1d_human_trigeminal_nerve_network_layout1, cls).getInteractiveFunctions()
        for interactiveFunction in interactiveFunctions:
            if interactiveFunction[0] == "Edit structure...":
                interactiveFunctions.remove(interactiveFunction)
                break
        return interactiveFunctions


class MeshType_3d_trigeminalnerve2(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network with core representing the human trigeminal nerve.
    """

    @classmethod
    def getName(cls):
        return "3D Trigeminal Nerve 2"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Base parameter set"] = useParameterSetName
        options["Trigeminal nerve network layout"] = ScaffoldPackage(MeshType_1d_human_trigeminal_nerve_network_layout1)
        options["Number of elements along nerve root"] = 4
        options["Number of elements along trigeminal ganglion"] = 2
        options["Number of elements along mandibular nerve"] = 4
        options["Number of elements along maxillary nerve"] = 4
        options["Number of elements along ophthalmic nerve"] = 4
        options["Number of elements around nerve root"] = 8
        # options["Number of elements around trigeminal ganglion"] = 12
        options["Number of elements around maxillary nerve"] = 8
        options["Number of elements around mandibular nerve"] = 8
        options["Number of elements around ophthalmic nerve"] = 8
        options["Number of elements through shell"] = 1
        options["Show trim surfaces"] = False
        options["Use Core"] = True
        options["Number of elements across core box minor"] = 2
        options["Number of elements across core transition"] = 1

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Trigeminal nerve network layout",
            "Number of elements along nerve root",
            "Number of elements along trigeminal ganglion",
            "Number of elements along maxillary nerve",
            "Number of elements along mandibular nerve",
            "Number of elements along ophthalmic nerve",
            "Number of elements around nerve root",
            # "Number of elements around trigeminal ganglion",
            "Number of elements around maxillary nerve",
            "Number of elements around mandibular nerve",
            "Number of elements around ophthalmic nerve",
            "Number of elements through shell",
            "Show trim surfaces",
            "Use Core",
            "Number of elements across core box minor",
            "Number of elements across core transition"]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Trigeminal nerve network layout":
            return [MeshType_1d_human_trigeminal_nerve_network_layout1]
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
        if optionName == "Trigeminal nerve network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_human_trigeminal_nerve_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if (options["Trigeminal nerve network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Trigeminal nerve network layout")):
            options["Trigeminal nerve network layout"] = ScaffoldPackage(MeshType_1d_human_trigeminal_nerve_network_layout1)
        for key in [
            "Number of elements along nerve root",
            "Number of elements along trigeminal ganglion",
            "Number of elements along maxillary nerve",
            "Number of elements along mandibular nerve",
            "Number of elements along ophthalmic nerve"
        ]:
            if options[key] < 1:
                options[key] = 1
        minElementsCountAround = None
        for key in [
            "Number of elements around nerve root",
            # "Number of elements around trigeminal ganglion",
            "Number of elements around maxillary nerve",
            "Number of elements around mandibular nerve",
            "Number of elements around ophthalmic nerve"
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
        networkLayout = options["Trigeminal nerve network layout"]
        elementsCountAlongNerveRoot = options["Number of elements along nerve root"]
        elementsCountAlongTrigeminalGanglion = options["Number of elements along trigeminal ganglion"]
        elementsCountAlongMaxillaryNerve = options["Number of elements along maxillary nerve"]
        elementsCountAlongMandibularNerve = options["Number of elements along mandibular nerve"]
        elementsCountAlongOphthalmicNerve = options["Number of elements along ophthalmic nerve"]
        elementsCountAroundNerveRoot = options["Number of elements around nerve root"]
        # elementsCountAroundTrigeminalGanglion = options["Number of elements around trigeminal ganglion"]
        elementsCountAroundMaxillaryNerve = options["Number of elements around maxillary nerve"]
        elementsCountAroundMandibularNerve = options["Number of elements around mandibular nerve"]
        elementsCountAroundOphthalmicNerve = options["Number of elements around ophthalmic nerve"]
        isCore = options["Use Core"]

        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        annotationAlongCounts = []
        annotationAroundCounts = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            alongCount = 0
            aroundCount = 0
            name = layoutAnnotationGroup.getName()
            if "nerve root" in name:
                alongCount = elementsCountAlongNerveRoot
                aroundCount = elementsCountAroundNerveRoot
            elif "ganglion" in name:
                alongCount = elementsCountAlongTrigeminalGanglion
                aroundCount = elementsCountAroundNerveRoot
            elif "maxillary" in name:
                alongCount = elementsCountAlongMaxillaryNerve
                aroundCount = elementsCountAroundMaxillaryNerve
            elif "mandibular" in name:
                alongCount = elementsCountAlongMandibularNerve
                aroundCount = elementsCountAroundMandibularNerve
            elif "ophthalmic" in name:
                alongCount = elementsCountAlongOphthalmicNerve
                aroundCount = elementsCountAroundOphthalmicNerve
            annotationAlongCounts.append(alongCount)
            annotationAroundCounts.append(aroundCount)

        tubeNetworkMeshBuilder = BodyTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=2.0,  # not used for body
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAlong=annotationAlongCounts,
            defaultElementsCountAround=options["Number of elements around nerve root"],
            annotationElementsCountsAround=annotationAroundCounts,
            elementsCountThroughShell=options["Number of elements through shell"],
            isCore=isCore,
            elementsCountTransition=options['Number of elements across core transition'],
            defaultElementsCountCoreBoxMinor=options["Number of elements across core box minor"],
            annotationElementsCountsCoreBoxMinor=[],
            useOuterTrimSurfaces=True)

        meshDimension = 3
        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, meshDimension,
            isLinearThroughShell=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None

    # @classmethod
    # def defineFaceAnnotations(cls, region, options, annotationGroups):
    #     """
    #     Add face annotation groups from the highest dimension mesh.
    #     Must have defined faces and added subelements for highest dimension groups.
    #     :param region: Zinc region containing model.
    #     :param options: Dict containing options. See getDefaultOptions().
    #     :param annotationGroups: List of annotation groups for top-level elements.
    #     New face annotation groups are appended to this list.
    #     """
    #     isCore = options["Use Core"]
    #
    #     # create 2-D surface mesh groups, 1-D spinal cord
    #     fieldmodule = region.getFieldmodule()
    #     mesh2d = fieldmodule.findMeshByDimension(2)
    #     mesh1d = fieldmodule.findMeshByDimension(1)
    #
    #     is_exterior = fieldmodule.createFieldIsExterior()
    #     is_face_xi3_0 = fieldmodule.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)
    #
    #     skinGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("skin epidermis"))
    #     is_skin = is_exterior if isCore else fieldmodule.createFieldAnd(
    #         is_exterior, fieldmodule.createFieldNot(is_face_xi3_0))
    #     skinGroup.getMeshGroup(mesh2d).addElementsConditional(is_skin)
    #
    #     leftArmGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left arm"))
    #     leftArmSkinGroup = findOrCreateAnnotationGroupForTerm(
    #         annotationGroups, region, get_body_term("left arm skin epidermis"))
    #     leftArmSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
    #         fieldmodule.createFieldAnd(leftArmGroup.getGroup(), is_exterior))
    #     rightArmGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right arm"))
    #     rightArmSkinGroup = findOrCreateAnnotationGroupForTerm(
    #         annotationGroups, region, get_body_term("right arm skin epidermis"))
    #     rightArmSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
    #         fieldmodule.createFieldAnd(rightArmGroup.getGroup(), is_exterior))
    #     leftLegGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left leg"))
    #     leftLegSkinGroup = findOrCreateAnnotationGroupForTerm(
    #         annotationGroups, region, get_body_term("left leg skin epidermis"))
    #     leftLegSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
    #         fieldmodule.createFieldAnd(leftLegGroup.getGroup(), is_exterior))
    #     rightLegGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right leg"))
    #     rightLegSkinGroup = findOrCreateAnnotationGroupForTerm(
    #         annotationGroups, region, get_body_term("right leg skin epidermis"))
    #     rightLegSkinGroup.getMeshGroup(mesh2d).addElementsConditional(
    #         fieldmodule.createFieldAnd(rightLegGroup.getGroup(), is_exterior))
    #
    #     if isCore:
    #         coreGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("core"))
    #         shellGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("shell"))
    #         leftGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("left"))
    #         rightGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("right"))
    #         dorsalGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("dorsal"))
    #
    #         is_core_shell = fieldmodule.createFieldAnd(coreGroup.getGroup(), shellGroup.getGroup())
    #         is_left_right = fieldmodule.createFieldAnd(leftGroup.getGroup(), rightGroup.getGroup())
    #         is_left_right_dorsal = fieldmodule.createFieldAnd(is_left_right, dorsalGroup.getGroup())
    #
    #         neckGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("neck"))
    #         thoracicCavityGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("thoracic cavity"))
    #         abdominalCavityGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("abdominal cavity"))
    #         armGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("arm"))
    #         legGroup = getAnnotationGroupForTerm(annotationGroups, get_body_term("leg"))
    #
    #         thoracicCavityBoundaryGroup = findOrCreateAnnotationGroupForTerm(
    #             annotationGroups, region, get_body_term("thoracic cavity boundary"))
    #         is_thoracic_cavity_boundary = fieldmodule.createFieldAnd(
    #             thoracicCavityGroup.getGroup(),
    #             fieldmodule.createFieldOr(
    #                 fieldmodule.createFieldOr(neckGroup.getGroup(), armGroup.getGroup()),
    #                 fieldmodule.createFieldOr(shellGroup.getGroup(), abdominalCavityGroup.getGroup())))
    #         thoracicCavityBoundaryGroup.getMeshGroup(mesh2d).addElementsConditional(is_thoracic_cavity_boundary)
    #
    #         abdominalCavityBoundaryGroup = findOrCreateAnnotationGroupForTerm(
    #             annotationGroups, region, get_body_term("abdominal cavity boundary"))
    #         is_abdominal_cavity_boundary = fieldmodule.createFieldAnd(
    #             abdominalCavityGroup.getGroup(),
    #             fieldmodule.createFieldOr(
    #                 thoracicCavityGroup.getGroup(),
    #                 fieldmodule.createFieldOr(shellGroup.getGroup(), legGroup.getGroup())))
    #         abdominalCavityBoundaryGroup.getMeshGroup(mesh2d).addElementsConditional(is_abdominal_cavity_boundary)
    #
    #         diaphragmGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("diaphragm"))
    #         is_diaphragm = fieldmodule.createFieldAnd(thoracicCavityGroup.getGroup(), abdominalCavityGroup.getGroup())
    #         diaphragmGroup.getMeshGroup(mesh2d).addElementsConditional(is_diaphragm)
    #
    #         spinalCordGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_body_term("spinal cord"))
    #         is_spinal_cord = fieldmodule.createFieldAnd(is_core_shell, is_left_right_dorsal)
    #         spinalCordGroup.getMeshGroup(mesh1d).addElementsConditional(is_spinal_cord)


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
