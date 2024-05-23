"""
Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches
"""

import copy
import math

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, get_group_list, find_or_create_field_group
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    findAnnotationGroupByName, getAnnotationGroupForTerm
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

from cmlibs.maths.vectorops import add, sub, mult, div, dot, matrix_det, matrix_mult, matrix_inv, normalize
from scaffoldmaker.utils.vector import magnitude_squared, magnitude, crossproduct3, vectorRejection, setMagnitude
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, remapEftNodeValueLabelVersion, \
    remapEftNodeValueLabelWithNodes, setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, getCubicHermiteBasisDerivatives, \
    interpolateCubicHermite, interpolateHermiteLagrange, smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters, set_nodeset_field_parameters, \
    make_nodeset_derivatives_orthogonal

from scaffoldfitter.fitter import Fitter
from scaffoldfitter.fitterstepalign import FitterStepAlign
from scaffoldfitter.fitterstepfit import FitterStepFit


class MeshType_3d_vagus_box2(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @staticmethod
    def getName():
        return "3D Vagus Box 2"

    @staticmethod
    def getParameterSetNames():
        return [
            'Human Left Trunk 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk': 45,
            'Apply fitting': True,
            'Add branches': True,
            'Add 3D box': True,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along the trunk',
            'Apply fitting',
            'Add branches',
            'Add 3D box'
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options['Number of elements along the trunk'] < 10:
            options['Number of elements along the trunk'] = 10
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        applyFitting = options['Apply fitting']
        addBranches = options['Add branches']
        add3D = options['Add 3D box']
        xmlData = False

        # load data from file
        if xmlData:
            marker_data, trunk_group_name, trunk_data, trunk_radius, branch_data, branch_radius_data = load_data_xml(region)
        else:
            marker_data, trunk_group_name, trunk_data, trunk_radius, branch_data, branch_radius_data = load_data(region)
        assert len(marker_data) >= 2, f"At least two landmarks are expected in the data. Incomplete data."

        # get fitted vagus region
        #get_fitted_vagus_nerve(region, marker_data, trunk_group_name, trunk_data, branch_data,
        #                       elementsAlongTrunk

        # setup
        annotationGroups = []

        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                       Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                       Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        # node - geometric coordinates
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for valueLabel in valueLabels[1:]:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 1)

        # branch special node
        nodeTemplateNoValue = nodes.createNodetemplate()
        nodeTemplateNoValue.defineField(coordinates)
        for valueLabel in valueLabels:
            if valueLabel == Node.VALUE_LABEL_VALUE:
                nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, valueLabel, 0)
            else:
                nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, valueLabel, 1)

        # vagus centroid
        mesh1d = fieldmodule.findMeshByDimension(1)
        hermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft1d = mesh1d.createElementfieldtemplate(hermiteBasis)
        linetemplate = mesh1d.createElementtemplate()
        linetemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplate.defineField(coordinates, -1, eft1d)

        # vagus centroid - branch root
        cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft1dNV = mesh1d.createElementfieldtemplate(hermiteBasis)
        eft1dNV.setNumberOfLocalNodes(4)
        eft1dNV.setNumberOfLocalScaleFactors(4)
        for i in range(4):
            eft1dNV.setScaleFactorType(i + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
        remapEftNodeValueLabelWithNodes(eft1dNV, 1, Node.VALUE_LABEL_VALUE,
                                        [(3, Node.VALUE_LABEL_VALUE, [1]),
                                         (3, Node.VALUE_LABEL_D_DS1, [2]),
                                         (4, Node.VALUE_LABEL_VALUE, [3]),
                                         (4, Node.VALUE_LABEL_D_DS1, [4])])
        linetemplateBranchRoot = mesh1d.createElementtemplate()
        linetemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplateBranchRoot.defineField(coordinates, -1, eft1dNV)

        # vagus box
        mesh3d = fieldmodule.findMeshByDimension(3)
        hermiteBilinearBasis = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        hermiteBilinearBasis.setFunctionType(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft3d = mesh3d.createElementfieldtemplate(hermiteBilinearBasis)
        setEftScaleFactorIds(eft3d, [1], [])
        ln = 1
        for n3 in range(2):
            s3 = [1] if (n3 == 0) else []
            for n2 in range(2):
                s2 = [1] if (n2 == 0) else []
                for n1 in range(2):
                    remapEftNodeValueLabel(eft3d, [ln], Node.VALUE_LABEL_VALUE,
                                           [(Node.VALUE_LABEL_VALUE, []),
                                            (Node.VALUE_LABEL_D_DS2, s2),
                                            (Node.VALUE_LABEL_D_DS3, s3)])
                    remapEftNodeValueLabel(eft3d, [ln], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []),
                                            (Node.VALUE_LABEL_D2_DS1DS2, s2),
                                            (Node.VALUE_LABEL_D2_DS1DS3, s3)])
                    ln += 1
        remapEftLocalNodes(eft3d, 2, [1, 2, 1, 2, 1, 2, 1, 2])
        elementtemplate = mesh3d.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft3d)

        # vagus box - branch root
        # 1 & 3 - trunk nodes, 2 - branch 2nd element
        cubicHermiteBilinearBasis = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        cubicHermiteBilinearBasis.setFunctionType(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft3dNV = mesh3d.createElementfieldtemplate(hermiteBilinearBasis)
        # setEftScaleFactorIds(eft3dNV, [1], [], 72)
        # ln = 1
        # for n3 in range(2):
        #     s3 = [1] if (n3 == 0) else []
        #     for n2 in range(2):
        #         s2 = [1] if (n2 == 0) else []
        #         for n1 in range(2):
        #             if ln % 2 == 0:
        #                 # local nodes 2, 4, 6, 8
        #                 remapEftNodeValueLabel(eft3d, [ln], Node.VALUE_LABEL_VALUE,
        #                                        [(Node.VALUE_LABEL_VALUE, []),
        #                                         (Node.VALUE_LABEL_D_DS2, s2),
        #                                         (Node.VALUE_LABEL_D_DS3, s3)])
        #                 remapEftNodeValueLabel(eft3d, [ln], Node.VALUE_LABEL_D_DS1,
        #                                        [(Node.VALUE_LABEL_D_DS1, []),
        #                                         (Node.VALUE_LABEL_D2_DS1DS2, s2),
        #                                         (Node.VALUE_LABEL_D2_DS1DS3, s3)])
        #                 ln += 1
        #             else:
        #                 # local nodes 1, 3, 5, 7
        #                 remapEftNodeValueLabelWithNodes(eft3dNV, ln, Node.VALUE_LABEL_VALUE,
        #                                                 [(1, Node.VALUE_LABEL_VALUE, si + 1),
        #                                                  (1, Node.VALUE_LABEL_D_DS1, si + 2),
        #                                                  (1, Node.VALUE_LABEL_D_DS2, si + 3),
        #                                                  (1, Node.VALUE_LABEL_D2_DS1DS2, si + 4),
        #                                                  (1, Node.VALUE_LABEL_D_DS3, si + 5),
        #                                                  (1, Node.VALUE_LABEL_D2_DS1DS3, si + 6),
        #                                                  (3, Node.VALUE_LABEL_VALUE, si + 7),
        #                                                  (3, Node.VALUE_LABEL_D_DS1, si + 8),
        #                                                  (3, Node.VALUE_LABEL_D_DS2, si + 9),
        #                                                  (3, Node.VALUE_LABEL_D2_DS1DS2, si + 10),
        #                                                  (3, Node.VALUE_LABEL_D_DS3, si + 11),
        #                                                  (3, Node.VALUE_LABEL_D2_DS1DS3, si + 12)])
        #                 remapEftNodeValueLabelWithNodes(eft3dNV, ln, Node.VALUE_LABEL_VALUE,
        #                                                 [(1, Node.VALUE_LABEL_VALUE, si + 1),
        #                                                  (1, Node.VALUE_LABEL_D_DS1, si + 2),
        #                                                  (1, Node.VALUE_LABEL_D_DS2, si + 3),
        #                                                  (1, Node.VALUE_LABEL_D2_DS1DS2, si + 4),
        #                                                  (1, Node.VALUE_LABEL_D_DS3, si + 5),
        #                                                  (1, Node.VALUE_LABEL_D2_DS1DS3, si + 6),
        #                                                  (3, Node.VALUE_LABEL_VALUE, si + 7),
        #                                                  (3, Node.VALUE_LABEL_D_DS1, si + 8),
        #                                                  (3, Node.VALUE_LABEL_D_DS2, si + 9),
        #                                                  (3, Node.VALUE_LABEL_D2_DS1DS2, si + 10),
        #                                                  (3, Node.VALUE_LABEL_D_DS3, si + 11),
        #                                                  (3, Node.VALUE_LABEL_D2_DS1DS3, si + 12)])

        setEftScaleFactorIds(eft3dNV, [1], [], 96)
        si = 1
        for ln in [1, 3, 5, 7]:
            for oldValueLabel in [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]:
                remapEftNodeValueLabelWithNodes(eft3dNV, ln, oldValueLabel,
                                                [(1, Node.VALUE_LABEL_VALUE, si + 1),
                                                 (1, Node.VALUE_LABEL_D_DS1, si + 2),
                                                 (1, Node.VALUE_LABEL_D_DS2, si + 3),
                                                 (1, Node.VALUE_LABEL_D2_DS1DS2, si + 4),
                                                 (1, Node.VALUE_LABEL_D_DS3, si + 5),
                                                 (1, Node.VALUE_LABEL_D2_DS1DS3, si + 6),
                                                 (3, Node.VALUE_LABEL_VALUE, si + 7),
                                                 (3, Node.VALUE_LABEL_D_DS1, si + 8),
                                                 (3, Node.VALUE_LABEL_D_DS2, si + 9),
                                                 (3, Node.VALUE_LABEL_D2_DS1DS2, si + 10),
                                                 (3, Node.VALUE_LABEL_D_DS3, si + 11),
                                                 (3, Node.VALUE_LABEL_D2_DS1DS3, si + 12)])
                si += 12
        ln = 2
        for n3 in range(2):
            s3 = [1] if (n3 == 0) else []
            for n2 in range(2):
                s2 = [1] if (n2 == 0) else []
                remapEftNodeValueLabel(eft3dNV, [ln], Node.VALUE_LABEL_VALUE,
                                       [(Node.VALUE_LABEL_VALUE, []),
                                        (Node.VALUE_LABEL_D_DS2, s2),
                                        (Node.VALUE_LABEL_D_DS3, s3)])
                remapEftNodeValueLabel(eft3dNV, [ln], Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, []),
                                        (Node.VALUE_LABEL_D2_DS1DS2, s2),
                                        (Node.VALUE_LABEL_D2_DS1DS3, s3)])
                ln += 2
        remapEftLocalNodes(eft3dNV, 3, [1, 2, 3, 2, 2, 2, 2, 2, 2])
        elementtemplateBranchRoot = mesh3d.createElementtemplate()
        elementtemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplateBranchRoot.defineField(coordinates, -1, eft3dNV)

        # vagus epineurium
        mesh2d = fieldmodule.findMeshByDimension(2)
        bicubichermiteSerendipityBasis = (
            fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
        scalefactors2d = [-1.0, 0.5, 0.25 * math.pi]
        # 4 elements around circle
        facetemplate_and_eft_list = [None] * 4
        for e in range(4):
            eft2d = mesh2d.createElementfieldtemplate(bicubichermiteSerendipityBasis)
            setEftScaleFactorIds(eft2d, [1, 2, 3], [])
            ln = 1
            for n2 in range(2):
                for n1 in range(2):
                    valueExpression = [(Node.VALUE_LABEL_VALUE, [])]
                    d_ds1Expression = [(Node.VALUE_LABEL_D_DS1, [])]
                    d_ds2Expression = []
                    pole = (e + n2) % 4
                    if pole == 0:
                        valueExpression.append((Node.VALUE_LABEL_D_DS2, [2]))
                        d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [2]))
                        d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [3]))
                    elif pole == 1:
                        valueExpression.append((Node.VALUE_LABEL_D_DS3, [2]))
                        d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [2]))
                        d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [1, 3]))
                    elif pole == 2:
                        valueExpression.append((Node.VALUE_LABEL_D_DS2, [1, 2]))
                        d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [1, 2]))
                        d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [1, 3]))
                    elif pole == 3:
                        valueExpression.append((Node.VALUE_LABEL_D_DS3, [1, 2]))
                        d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [1, 2]))
                        d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [3]))

                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_VALUE, valueExpression)
                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_D_DS1, d_ds1Expression)
                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_D_DS2, d_ds2Expression)
                    ln += 1

            remapEftLocalNodes(eft2d, 2, [1, 2, 1, 2])
            facetemplate = mesh2d.createElementtemplate()
            facetemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplate.defineField(coordinates, -1, eft2d)
            facetemplate_and_eft_list[e] = (facetemplate, eft2d)

        # vagus epineureum - branch root
        # 1 & 3 - trunk nodes, 2 - branch 2nd element
        bicubichermiteSerendipityBasisNV = (
            fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
        # 4 elements around circle
        facetemplate_and_eft_list_NV = [None] * 4
        for e in range(4):
            eft2dNV = mesh2d.createElementfieldtemplate(bicubichermiteSerendipityBasisNV)
            setEftScaleFactorIds(eft2dNV, [1, 2, 3], [], 72)
            # ln = 1
            # si = 4
            # for n2 in range(2):
            #     for n1 in range(2):
            #         pole = (e + n2) % 4
            #         if ln % 2 == 0:
            #             valueExpression = [(Node.VALUE_LABEL_VALUE, [])]
            #             d_ds1Expression = [(Node.VALUE_LABEL_D_DS1, [])]
            #             d_ds2Expression = []
            #             if pole == 0:
            #                 valueExpression.append((Node.VALUE_LABEL_D_DS2, [2]))
            #                 d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [2]))
            #                 d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [3]))
            #             elif pole == 1:
            #                 valueExpression.append((Node.VALUE_LABEL_D_DS3, [2]))
            #                 d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [2]))
            #                 d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [1, 3]))
            #             elif pole == 2:
            #                 valueExpression.append((Node.VALUE_LABEL_D_DS2, [1, 2]))
            #                 d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [1, 2]))
            #                 d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [1, 3]))
            #             elif pole == 3:
            #                 valueExpression.append((Node.VALUE_LABEL_D_DS3, [1, 2]))
            #                 d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [1, 2]))
            #                 d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [3]))
            #
            #             remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_VALUE, valueExpression)
            #             remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_D_DS1, d_ds1Expression)
            #             remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_D_DS2, d_ds2Expression)
            #         else:
            #             sf_count = 12
            #             for pn in range(2):
            #                 for si in range(12):
            #         ln += 1

            si = 3
            for ln in [1, 3]:
                for oldValueLabel in [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2]:
                    remapEftNodeValueLabelWithNodes(eft2dNV, ln, oldValueLabel,
                                                    [(1, Node.VALUE_LABEL_VALUE, si + 1),
                                                     (1, Node.VALUE_LABEL_D_DS1, si + 2),
                                                     (1, Node.VALUE_LABEL_D_DS2, si + 3),
                                                     (1, Node.VALUE_LABEL_D2_DS1DS2, si + 4),
                                                     (1, Node.VALUE_LABEL_D_DS3, si + 5),
                                                     (1, Node.VALUE_LABEL_D2_DS1DS3, si + 6),
                                                     (3, Node.VALUE_LABEL_VALUE, si + 7),
                                                     (3, Node.VALUE_LABEL_D_DS1, si + 8),
                                                     (3, Node.VALUE_LABEL_D_DS2, si + 9),
                                                     (3, Node.VALUE_LABEL_D2_DS1DS2, si + 10),
                                                     (3, Node.VALUE_LABEL_D_DS3, si + 11),
                                                     (3, Node.VALUE_LABEL_D2_DS1DS3, si + 12)])
                    si += 12
            ln = 2
            for n2 in range(2):
                valueExpression = [(Node.VALUE_LABEL_VALUE, [])]
                d_ds1Expression = [(Node.VALUE_LABEL_D_DS1, [])]
                d_ds2Expression = []
                pole = (e + n2) % 4
                if pole == 0:
                    valueExpression.append((Node.VALUE_LABEL_D_DS2, [2]))
                    d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [2]))
                    d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [3]))
                elif pole == 1:
                    valueExpression.append((Node.VALUE_LABEL_D_DS3, [2]))
                    d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [2]))
                    d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [1, 3]))
                elif pole == 2:
                    valueExpression.append((Node.VALUE_LABEL_D_DS2, [1, 2]))
                    d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS2, [1, 2]))
                    d_ds2Expression.append((Node.VALUE_LABEL_D_DS3, [1, 3]))
                elif pole == 3:
                    valueExpression.append((Node.VALUE_LABEL_D_DS3, [1, 2]))
                    d_ds1Expression.append((Node.VALUE_LABEL_D2_DS1DS3, [1, 2]))
                    d_ds2Expression.append((Node.VALUE_LABEL_D_DS2, [3]))

                remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_VALUE, valueExpression)
                remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_D_DS1, d_ds1Expression)
                remapEftNodeValueLabel(eft2dNV, [ln], Node.VALUE_LABEL_D_DS2, d_ds2Expression)
                ln += 2

            remapEftLocalNodes(eft2dNV, 3, [1, 2, 3, 2])
            facetemplateNV = mesh2d.createElementtemplate()
            facetemplateNV.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplateNV.defineField(coordinates, -1, eft2dNV)
            facetemplate_and_eft_list_NV[e] = (facetemplateNV, eft2dNV)


        termNameVagusLengthList = {
            # "level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
            "level of superior border of jugular foramen on the vagus nerve": 8.6342,
            "level of inferior border of jugular foramen on the vagus nerve": 16.7227,
            "level of C1 transverse process on the vagus nerve": 32.1129,
            "level of angle of mandible on the vagus nerve": 42.2450,
            "level of greater horn of hyoid on the vagus nerve": 45.6122,
            "level of carotid bifurcation on the vagus nerve": 48.3581,
            "level of laryngeal prominence on the vagus nerve": 68.8431,
            "level of superior border of the clavicle on the vagus nerve": 117.5627,
            "level of jugular notch on the vagus nerve": 124.6407,
            "level of sternal angle on the vagus nerve": 151.2352,
            "1 cm superior to start of esophageal plexus on the vagus nerve": 165.5876,
            "level of esophageal hiatus on the vagus nerve": 254.32879,
            "level of aortic hiatus on the vagus nerve": 291.3695,
            # "level of end of trunk on Japanese dataset": 312.5 # note this term is also not on the list of annotations
        }

        lengthToDiameterRatio = 312.5  # calculated from total length of nerve/average diameter of nerve
        rescaledlengthToDiameterRatio = 100.0
        rescaledTermNameVagusLengthList = {}
        for term in termNameVagusLengthList:
            rescaledTermNameVagusLengthList[term] = termNameVagusLengthList[term] / lengthToDiameterRatio \
                                                    * rescaledlengthToDiameterRatio

        # choose markers for building initial scaffold
        use_marker_names = [list(marker_data.keys())[0],
                            list(marker_data.keys())[-1]]
        assert [name.lower() in marker_data for name in use_marker_names] and \
               [name.lower() in rescaledTermNameVagusLengthList for name in use_marker_names]

        # build 1d trunk centroid line
        elementsAlongTrunk = options['Number of elements along the trunk']
        elementLength = rescaledlengthToDiameterRatio / (elementsAlongTrunk - 1)
        tx, td1, vtx, vtd1 = estimate_trunk_coordinates(elementsAlongTrunk, elementLength,
                                                        marker_data, rescaledTermNameVagusLengthList)

        trunkCentroidGroup = AnnotationGroup(region, (trunk_group_name, ""))
        annotationGroups.append(trunkCentroidGroup)
        trunkCentroidMeshGroup = trunkCentroidGroup.getMeshGroup(mesh1d)

        nodeIdentifier = 1
        lineIdentifier = 1
        for n in range(elementsAlongTrunk):
            sx = tx[n]
            sd1 = td1

            node = nodes.createNode(nodeIdentifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

            if n > 0:
                nids = [nodeIdentifier - 1, nodeIdentifier]
                line = mesh1d.createElement(lineIdentifier, linetemplate)
                line.setNodesByIdentifier(eft1d, nids)
                trunkCentroidMeshGroup.addElement(line)
                lineIdentifier += 1
            nodeIdentifier += 1

        # set markers
        for marker_name, marker_coordinate in marker_data.items():
            annotationGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
            annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
            nodeIdentifier += 1

        # create temporary model file
        sir = region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf")
        region.write(sir)

        fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
        fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf"
        fitter = fit_trunk_model(fitter_model_file, fitter_data_file)
        set_fitted_group_nodes(region, fitter)

        trunkNodesetGroup = trunkCentroidGroup.getNodesetGroup(nodes)
        trunk_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        trunk_parameters = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name, valueLabels)

        # branch building
        if addBranches:
            branch_root_params = {}
            for branch_name in branch_data.keys():
                print(branch_name)

                # determine branch approximate start and closest trunk node index
                branch_coordinates = [branch_node[0] for branch_node in branch_data[branch_name]]
                bx, bd1, trunk_node_id, branch_root_xi, branch_root_x, elementsAlongBranch = \
                    estimate_branch_coordinates(branch_coordinates, trunk_parameters, elementLength)
                branch_root_params[branch_name] = [trunk_node_id, trunk_node_id + 1, branch_root_xi, branch_root_x]

                branchCentroidGroup = AnnotationGroup(region, (branch_name, 'None'))
                annotationGroups.append(branchCentroidGroup)
                branchCentroidMeshGroup = branchCentroidGroup.getMeshGroup(mesh1d)

                for n in range(elementsAlongBranch):
                    sx = bx[n]
                    sd1 = bd1

                    if n == 0:
                        # create branch special node
                        node = nodes.createNode(nodeIdentifier, nodeTemplateNoValue)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

                        if n == 1:
                            # create branch special element
                            nids = [nodeIdentifier - 1, nodeIdentifier,
                                    trunk_node_id, trunk_node_id + 1]
                            line = mesh1d.createElement(lineIdentifier, linetemplateBranchRoot)
                            line.setNodesByIdentifier(eft1dNV, nids)
                            scalefactorsNV = getCubicHermiteBasis(branch_root_xi)
                            line.setScaleFactors(eft1dNV, list(scalefactorsNV))
                            branchCentroidMeshGroup.addElement(line)
                            lineIdentifier += 1
                        else:
                            nids = [nodeIdentifier - 1, nodeIdentifier]
                            line = mesh1d.createElement(lineIdentifier, linetemplate)
                            line.setNodesByIdentifier(eft1d, nids)
                            branchCentroidMeshGroup.addElement(line)
                            lineIdentifier += 1
                    nodeIdentifier += 1

                # remove trunk nodes from branch group
                branchNodesetGroup = branchCentroidGroup.getNodesetGroup(nodes)
                if branchNodesetGroup.isValid():
                    branchNodesetGroup.removeNodesConditional(trunk_group)

            sir = region.createStreaminformationRegion()
            srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_model.exf")
            region.write(sir)

            if applyFitting:
                # geometry fitting - branches
                fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
                fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_model.exf"

                fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
                set_fitted_group_nodes(region, fitter, branch_name)

            # extract first branch node - x & d1 value
            branch_group = find_or_create_field_group(fieldmodule, branch_name)
            branch_nodes = branch_group.getNodesetGroup(nodes)
            node_iter = branch_nodes.createNodeiterator()
            node = node_iter.next()
            fieldcache.setNode(node)
            _, sd1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            branch_root_params[branch_name].append(sd1)

        # add box and epineureum - trunk and branches
        vagusBoxTrunkGroup = AnnotationGroup(region, ("left vagus box", 'None'))
        annotationGroups.append(vagusBoxTrunkGroup)
        vagusBoxMeshGroup = vagusBoxTrunkGroup.getMeshGroup(mesh3d)

        vagusEpineuriumAnnotationGroup = AnnotationGroup(region, ("left vagus epineurium", ""))
        annotationGroups.append(vagusEpineuriumAnnotationGroup)
        vagusEpineuriumMeshGroup = vagusEpineuriumAnnotationGroup.getMeshGroup(mesh2d)

        # set side and cross derivatives - d2, d3, d12, d13
        setSideCrossDerivatives(trunkNodesetGroup, coordinates, False)

        elementIdentifier = 1
        faceIdentifier = 1
        for n in range(1, elementsAlongTrunk):
            node = nodes.findNodeByIdentifier(n + 1)
            node_id = node.getIdentifier()
            fieldcache.setNode(node)

            if n > 0:
                # line = mesh1d.findElementByIdentifier(elementIdentifier)
                # print(line.getNode(eft1d, 1).getIdentifier(), line.getNode(eft1d, 2).getIdentifier())
                nids = [node_id - 1, node_id]

                element = mesh3d.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft3d, nids)
                element.setScaleFactors(eft3d, [-1.0])
                vagusBoxMeshGroup.addElement(element)
                elementIdentifier += 1

                for e in range(4):
                    facetemplate, eft2d = facetemplate_and_eft_list[e]
                    face = mesh2d.createElement(faceIdentifier, facetemplate)
                    face.setNodesByIdentifier(eft2d, nids)
                    face.setScaleFactors(eft2d, scalefactors2d)
                    vagusEpineuriumMeshGroup.addElement(face)
                    faceIdentifier += 1

        if addBranches:
            for branch_name in branch_data.keys():
                branchBoxTrunkGroup = AnnotationGroup(region, (branch_name + ' box', 'None'))
                annotationGroups.append(branchBoxTrunkGroup)
                branchBoxMeshGroup = branchBoxTrunkGroup.getMeshGroup(mesh3d)

                branchEpineuriumAnnotationGroup = AnnotationGroup(region, (branch_name + " epineurium", ""))
                annotationGroups.append(branchEpineuriumAnnotationGroup)
                branchEpineuriumMeshGroup = branchEpineuriumAnnotationGroup.getMeshGroup(mesh2d)

                branch_group = find_or_create_field_group(fieldmodule, branch_name)
                branch_nodes = branch_group.getNodesetGroup(nodes)

                # set derivatives d2, d3, d12, d13 for nodes
                setSideCrossDerivatives(branch_nodes, coordinates, True)

                node_iter = branch_nodes.createNodeiterator()
                for n in range(0, branch_nodes.getSize()):
                    node = node_iter.next()
                    node_id = node.getIdentifier()
                    fieldcache.setNode(node)

                    if n == 1:
                        print(branch_name)
                        print(branch_root_params[branch_name])

                        trunk_segment_start_id = branch_root_params[branch_name][0]
                        trunk_segment_end_id = branch_root_params[branch_name][1]
                        branch_root_xi = branch_root_params[branch_name][2]

                        # branch start data
                        # assume d12, d13 are both zero
                        bd1 = branch_root_params[branch_name][4]
                        bd2, bd3 = set_group_nodes_derivatives_orthogonal([bd1])

                        # trunk interpolation
                        fns = list(getCubicHermiteBasis(branch_root_xi))  # for x, d2, d3
                        dfns = list(getCubicHermiteBasisDerivatives(branch_root_xi))  # for d1, d12, d13

                        # access to two trunk segments to get interpolated values?
                        node = nodes.findNodeByIdentifier(trunk_segment_start_id)
                        fieldcache.setNode(node)
                        _, tx_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                        _, td1_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                        _, td2_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                        _, td12_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
                        _, td3_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                        _, td13_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)
                        node = nodes.findNodeByIdentifier(trunk_segment_end_id)
                        fieldcache.setNode(node)
                        _, tx_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                        _, td1_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                        _, td2_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                        _, td12_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
                        _, td3_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                        _, td13_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)

                        td1 = [dot(dfns, [tx_1[0], td1_1[0], tx_2[0], td1_2[0]]),
                               dot(dfns, [tx_1[1], td1_1[1], tx_2[1], td1_2[1]]),
                               dot(dfns, [tx_1[2], td1_1[2], tx_2[2], td1_2[2]])]
                        td2 = [dot(fns, [td2_1[0], td12_1[0], td2_2[0], td12_2[0]]),
                               dot(fns, [td2_1[1], td12_1[1], td2_2[1], td12_2[1]]),
                               dot(fns, [td2_1[2], td12_1[2], td2_2[2], td12_2[2]])]
                        td3 = [dot(fns, [td3_1[0], td13_1[0], td3_2[0], td13_2[0]]),
                               dot(fns, [td3_1[1], td13_1[1], td3_2[1], td13_2[1]]),
                               dot(fns, [td3_1[2], td13_1[2], td3_2[2], td13_2[2]])]

                        basis_from = [td1, td2, td3]
                        basis_to = [bd1, bd2[0], bd3[0]]
                        coefs = matrix_mult(basis_to, matrix_inv(basis_from))

                        # value, ds1, ds2, ds12, ds3, ds13
                        scalefactorsX = [fns[0], fns[1], 0, 0, 0, 0,
                                         fns[2], fns[3], 0, 0, 0, 0]
                        scalefactorsD1 = [coefs[0][0] * dfns[0], coefs[0][0] * dfns[1],  # value, ds1
                                          coefs[0][1] * fns[0], coefs[0][1] * fns[1],
                                          coefs[0][2] * fns[0], coefs[0][2] * fns[1],
                                          coefs[0][0] * dfns[2], coefs[0][0] * dfns[3],  # value, ds1
                                          coefs[0][1] * fns[2], coefs[0][1] * fns[3],
                                          coefs[0][2] * fns[2], coefs[0][2] * fns[3]]
                        scalefactorsD2 = [coefs[1][0] * dfns[0], coefs[1][0] * dfns[1],  # value, ds1
                                          coefs[1][1] * fns[0], coefs[1][1] * fns[1],
                                          coefs[1][2] * fns[0], coefs[1][2] * fns[1],
                                          coefs[1][0] * dfns[2], coefs[1][0] * dfns[3],  # value, ds1
                                          coefs[1][1] * fns[2], coefs[1][1] * fns[3],
                                          coefs[1][2] * fns[2], coefs[1][2] * fns[3]]
                        scalefactorsD3 = [coefs[2][0] * dfns[0], coefs[2][0] * dfns[1],  # value, ds1
                                          coefs[2][1] * fns[0], coefs[2][1] * fns[1],
                                          coefs[2][2] * fns[0], coefs[2][2] * fns[1],
                                          coefs[2][0] * dfns[2], coefs[2][0] * dfns[3],  # value, ds1
                                          coefs[2][1] * fns[2], coefs[2][1] * fns[3],
                                          coefs[2][2] * fns[2], coefs[2][2] * fns[3]]
                        scalefactorsD12 = [0, 0, dfns[0], dfns[1], 0, 0,
                                           0, 0, dfns[2], dfns[3], 0, 0]
                        scalefactorsD13 = [0, 0, 0, 0, dfns[0], dfns[1],
                                           0, 0, 0, 0, dfns[2], dfns[3]]

                        scalefactors = [-1]
                        scalefactors.extend(sub(sub(scalefactorsX, scalefactorsD2), scalefactorsD3))
                        scalefactors.extend(sub(sub(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                        scalefactors.extend(sub(add(scalefactorsX, scalefactorsD2), scalefactorsD3))
                        scalefactors.extend(sub(add(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                        scalefactors.extend(add(sub(scalefactorsX, scalefactorsD2), scalefactorsD3))
                        scalefactors.extend(add(sub(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                        scalefactors.extend(add(add(scalefactorsX, scalefactorsD2), scalefactorsD3))
                        scalefactors.extend(add(add(scalefactorsD1, scalefactorsD12), scalefactorsD13))

                        nids = [trunk_segment_start_id, node_id, trunk_segment_end_id]

                        element = mesh3d.createElement(elementIdentifier, elementtemplateBranchRoot)
                        element.setNodesByIdentifier(eft3dNV, nids)
                        element.setScaleFactors(eft3dNV, scalefactors)
                        branchBoxMeshGroup.addElement(element)
                        elementIdentifier += 1

                        for e in range(4):
                            scalefactors = scalefactors2d[:]
                            if e == 0:
                                # poles 0 & 1
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD2, 0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD12, 0.5)))
                                scalefactors.extend(mult(scalefactorsD3, 0.25 * math.pi))
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD3, 0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD13, 0.5)))
                                scalefactors.extend(mult(scalefactorsD2, -0.25 * math.pi))
                            if e == 1:
                                # poles 1 & 2
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD3, 0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD13, 0.5)))
                                scalefactors.extend(mult(scalefactorsD2, -0.25 * math.pi))
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD2, -0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD12, -0.5)))
                                scalefactors.extend(mult(scalefactorsD3, -0.25 * math.pi))
                            if e == 2:
                                # poles 2 & 3
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD2, -0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD12, -0.5)))
                                scalefactors.extend(mult(scalefactorsD3, -0.25 * math.pi))
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD3, -0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD13, -0.5)))
                                scalefactors.extend(mult(scalefactorsD2, 0.25 * math.pi))
                            if e == 3:
                                # poles 3 & 0
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD3, -0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD13, -0.5)))
                                scalefactors.extend(mult(scalefactorsD2, 0.25 * math.pi))
                                scalefactors.extend(add(scalefactorsX, mult(scalefactorsD2, 0.5)))
                                scalefactors.extend(add(scalefactorsD1, mult(scalefactorsD12, 0.5)))
                                scalefactors.extend(mult(scalefactorsD3, 0.25 * math.pi))

                            facetemplateNV, eft2dNV = facetemplate_and_eft_list_NV[e]
                            face = mesh2d.createElement(faceIdentifier, facetemplateNV)
                            face.setNodesByIdentifier(eft2dNV, nids)
                            face.setScaleFactors(eft2dNV, scalefactors)
                            branchEpineuriumMeshGroup.addElement(face)
                            faceIdentifier += 1

                    if n > 1:
                        nids = [node_id - 1, node_id]

                        element = mesh3d.createElement(elementIdentifier, elementtemplate)
                        element.setNodesByIdentifier(eft3d, nids)
                        element.setScaleFactors(eft3d, [-1.0])
                        branchBoxMeshGroup.addElement(element)
                        elementIdentifier += 1

                        for e in range(4):
                            facetemplate, eft2d = facetemplate_and_eft_list[e]
                            face = mesh2d.createElement(faceIdentifier, facetemplate)
                            face.setNodesByIdentifier(eft2d, nids)
                            face.setScaleFactors(eft2d, scalefactors2d)
                            branchEpineuriumMeshGroup.addElement(face)
                            faceIdentifier += 1

        sir = region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_final_model.exf")
        region.write(sir)

        return annotationGroups, None

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Override in classes with face annotation groups.
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """

        # Create 2d surface mesh groups
        fieldmodule = region.getFieldmodule()
        mesh2d = fieldmodule.findMeshByDimension(2)
        mesh1d = fieldmodule.findMeshByDimension(1)

        vagusEpineuriumAnnotationGroup = getAnnotationGroupForTerm(annotationGroups, ("left vagus epineurium", ""))
        vagusEpineuriumMeshGroup = vagusEpineuriumAnnotationGroup.getMeshGroup(mesh2d)
        vagusAnteriorLineAnnotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups,
                                                                              region, ("left vagus anterior line", ""))
        vagusAnteriorLineMeshGroup = vagusAnteriorLineAnnotationGroup.getMeshGroup(mesh1d)

        faceIterator = vagusEpineuriumMeshGroup.createElementiterator()
        quadrant = 0
        face = faceIterator.next()
        while face.isValid():
            if quadrant == 0:
                line = face.getFaceElement(4)
                vagusAnteriorLineMeshGroup.addElement(line)
            quadrant = (quadrant + 1) % 4
            face = faceIterator.next()


# supplementary functions
def get_nodeset_fieldgroup_parameters(nodeset, field, group_name, valueLabels):
    """

    """

    fieldmodule = nodeset.getFieldmodule()
    finite_element_field = field.castFiniteElement()
    assert finite_element_field.isValid(), "get_nodeset_fieldgroup_parameters:  Field is not finite element type"

    components_count = field.getNumberOfComponents()
    fieldcache = fieldmodule.createFieldcache()

    group = fieldmodule.findFieldByName(group_name).castGroup()

    node_fieldgroup_parameters = []
    if group.isValid():
        group_nodes = group.getNodesetGroup(nodeset)
        node_iterator = group_nodes.createNodeiterator()
        node = node_iterator.next()
        while node.isValid():
            fieldcache.setNode(node)
            node_parameters = []
            field_defined_at_node = False
            for valueLabel in valueLabels:
                result, parameters = finite_element_field.getNodeParameters(fieldcache, -1, valueLabel, 1,
                                                                            components_count)
                field_defined_at_node = True
                node_parameters.append(parameters)
            node_fieldgroup_parameters.append(node_parameters)
            node = node_iterator.next()

    return node_fieldgroup_parameters


def setSideCrossDerivatives(nodeset, field, isBranch=False):
    """

    """

    getValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
    _, node_field_parameters = get_nodeset_field_parameters(nodeset, field, getValueLabels)

    setValueLabels = [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                      Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

    if isBranch:
        node_field_parameters.pop(0)

    x = [nodeParameter[1][0][0] for nodeParameter in node_field_parameters]
    d1 = [nodeParameter[1][1][0] for nodeParameter in node_field_parameters]
    d2, d3 = set_group_nodes_derivatives_orthogonal(d1)
    d12, d13 = smoothCurveSideCrossDerivatives(x, d1, [d2, d3])

    index = 0
    new_node_field_parameters = []
    for node_identifier, node_parameters in node_field_parameters:
        new_node_parameters = [[d2[index]], [d12[index]], [d3[index]], [d13[index]]]
        new_node_field_parameters.append((node_identifier, new_node_parameters))
        index += 1
    set_nodeset_field_parameters(nodeset, field, setValueLabels, new_node_field_parameters)


def set_group_nodes_derivatives_orthogonal(d1):
    """

    """

    yx = [0.0, 1.0, 0.0]
    zx = [0.0, 0.0, 1.0]

    d2 = []
    d3 = []
    for i in range(len(d1)):
        td2 = vectorRejection(yx, d1[i])
        td2 = setMagnitude(td2, magnitude(yx))  # change magnitude later with radius
        d2.append(td2)

        td3 = crossproduct3(d1[i], td2)
        td3 = setMagnitude(td3, magnitude(zx))  # change magnitude later with radius
        d3.append(td3)

    return d2, d3


def estimate_trunk_coordinates(elementsAlongTrunk, step, marker_data, termNameVagusLengthList):
    """

    """

    # choose markers for building initial scaffold
    # for now just use 1st and last markers from the list,
    # later should be using markers corresponding to either left or right vagus
    use_marker_names = [list(marker_data.keys())[0],
                        list(marker_data.keys())[-1]]
    assert [name.lower() in termNameVagusLengthList for name in use_marker_names] or \
           ['left ' + name.lower() in termNameVagusLengthList for name in use_marker_names] or \
           ['right ' + name.lower() in termNameVagusLengthList for name in use_marker_names]

    x1, y1, z1 = marker_data[use_marker_names[0]]
    x2, y2, z2 = marker_data[use_marker_names[1]]
    t1 = termNameVagusLengthList[use_marker_names[0].replace('left ', '', 1).replace('right ', '', 1)]
    t2 = termNameVagusLengthList[use_marker_names[1].replace('left ', '', 1).replace('right ', '', 1)]

    dx, dy, dz = [(x2 - x1) / (t2 - t1),
                  (y2 - y1) / (t2 - t1),
                  (z2 - z1) / (t2 - t1)]
    trunk_ld1 = [dx * step, dy * step, dz * step]

    trunk_nodes = []
    vagus_trunk_nodes = []
    vagus_trunk_ld1 = [0, 0, step]
    for i in range(elementsAlongTrunk):
        trunk_nodes.append([x1 + dx * (i * step - t1),
                            y1 + dy * (i * step - t1),
                            z1 + dz * (i * step - t1)])
        vagus_trunk_nodes.append([0, 0, i * step])

    return trunk_nodes, trunk_ld1, vagus_trunk_nodes, vagus_trunk_ld1


def estimate_branch_coordinates(region, branch_coordinates, step):
    """

    """

    fm = region.getFieldmodule()
    fieldcache = fm.createFieldcache()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    branch_start_x, branch_end_x, parent_s_nid, parent_f_nid, parent_group_name = find_branch_start_parent_segment(
        region, branch_coordinates)

    # determine parent hermite curve parameters
    node = nodes.findNodeByIdentifier(parent_s_nid)
    fieldcache.setNode(node)
    _, px_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, pd1_1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

    node = nodes.findNodeByIdentifier(parent_f_nid)
    fieldcache.setNode(node)
    _, px_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
    _, pd1_2 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

    # find xi closest to branch_start on a cubic Hermite curve by bisection
    xi_a = 0
    xi_b = 1
    eps = 0.005
    while (xi_b - xi_a) > eps:
        dsq_a = magnitude_squared(sub(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_a)))
        dsq_b = magnitude_squared(sub(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_b)))
        if dsq_a >= dsq_b:
            xi_a = (xi_a + xi_b) / 2
        else:
            xi_b = (xi_a + xi_b) / 2
    branch_root_xi = (xi_a + xi_b) / 2
    #print('\t xi = %s' % branch_root_xi)

    # recalculate branch start parameters
    branch_start_x = interpolateHermiteLagrange(px_1, pd1_1, px_2, branch_root_xi)
    #branch_start_x = interpolateLagrangeHermite(px_1, px_2, pd1_2, branch_root_xi)

    # determine branch approximate end
    # max_distance_squared = 0
    # for branch_node in branch_coordinates:
    #     distance_squared = magnitude_squared(sub(branch_node, branch_start_x))
    #     if distance_squared >= max_distance_squared:
    #         max_distance_squared = distance_squared
    #         branch_end_x = branch_node

    branch_length = magnitude(sub(branch_start_x, branch_end_x))
    elementsAlongBranch = int(branch_length / step)
    if elementsAlongBranch < 3:
        elementsAlongBranch = 3
    dx, dy, dz = div(sub(branch_end_x, branch_start_x), (elementsAlongBranch - 1))

    branch_coordinates = []
    for i in range(elementsAlongBranch):
        branch_coordinates.append([branch_start_x[0] + dx * i,
                                   branch_start_x[1] + dy * i,
                                   branch_start_x[2] + dz * i])

    return branch_coordinates, [dx, dy, dz], parent_s_nid, parent_f_nid, parent_group_name, branch_root_xi, \
           elementsAlongBranch


def find_branch_start_parent_segment(region, branch_coordinates):
    """

    """

    def find_branch_ends(branch_coordinates):
        """
        Given list of node coordinates, find two furthest from each other nodes.
        Returns a list of two node coordinates.
        """

        max_distance_squared = 0
        for branch_node_x in branch_coordinates:
            for branch_node_y in branch_coordinates:
                distance_squared = magnitude_squared(sub(branch_node_x, branch_node_y))
                if distance_squared >= max_distance_squared:
                    max_distance_squared = distance_squared
                    branch_ends_points = [branch_node_x, branch_node_y]
        return branch_ends_points


    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    branch_ends_points = find_branch_ends(branch_coordinates)

    group_list = get_group_list(fm)
    dsq_to_branch = []
    tmp = []
    for group in group_list:
        if group.getName() == 'marker':
            group_list.remove(group)
        else:
            group_nodeset = group.getNodesetGroup(nodes)
            _, group_parameters = get_nodeset_field_parameters(group_nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
            group_ids = [parameter[0] for parameter in group_parameters]
            group_x = [parameter[1][0][0] for parameter in group_parameters]

            # find closest to branch distance, branch start, group node closest to the branch
            min_dsq = float('inf')
            for i in range(len(group_x)):
                node_x = group_x[i]
                if node_x is None:
                    continue
                for branch_point in branch_ends_points:
                    distance_squared = magnitude_squared(sub(node_x, branch_point))
                    if distance_squared <= min_dsq:
                        min_dsq = distance_squared
                        branch_start = branch_point
                        closest_index = i

            # determine segment closest to branch  (previous or next to the node)
            if closest_index == len(group_x) - 1:
                closest_index -= 1
            if closest_index > 0 and (closest_index < len(group_x) - 1):
                # TODO: should think about the case when closest_index = 0:
                # actual starting point might not have a node,
                # therefore first segment is not always taken into consideration
                dsq_prev_node = magnitude_squared(sub(branch_start, group_x[closest_index - 1]))
                dsq_next_node = magnitude_squared(sub(branch_start, group_x[closest_index + 1]))
                if dsq_prev_node < dsq_next_node:
                    closest_index -= 1

            closest_s_node_id = group_ids[closest_index]
            closest_f_node_id = group_ids[closest_index + 1]

            dsq_to_branch.append(min_dsq)
            tmp.append([closest_s_node_id, closest_f_node_id, branch_start])

    parent_index = dsq_to_branch.index(min(dsq_to_branch))
    parent_s_node_id = tmp[parent_index][0]
    parent_f_node_id = tmp[parent_index][1]
    branch_start = tmp[parent_index][2]
    branch_end = branch_ends_points[0] if branch_ends_points[1] == branch_start else branch_ends_points[1]
    parent_group_name = group_list[parent_index].getName()

    return branch_start, branch_end, parent_s_node_id, parent_f_node_id, parent_group_name


def set_fitted_group_nodes(region, fitter_region, group_name=None):
    """
    Updates node values for a particular field group or all region after fitting.
    """

    fieldmodule = region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()
    coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    fitter_fieldmodule = fitter_region.getFieldmodule()
    fitter_fieldcache = fitter_fieldmodule.createFieldcache()
    fitter_coordinates = fitter_region.getModelCoordinatesField().castFiniteElement()
    fitter_nodes = fitter_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    if group_name:
        group = find_or_create_field_group(fitter_fieldmodule, group_name)
        if group.isValid():
            fitter_nodes = group.getNodesetGroup(fitter_nodes)

    # reset trunk nodes with the fitted nodes
    fitter_node_iter = fitter_nodes.createNodeiterator()
    fitter_node = fitter_node_iter.next()
    while fitter_node.isValid():
        fitter_fieldcache.setNode(fitter_node)
        _, lx = fitter_coordinates.getNodeParameters(fitter_fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        _, ld1 = fitter_coordinates.getNodeParameters(fitter_fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

        node = nodes.findNodeByIdentifier(fitter_node.getIdentifier())
        fieldcache.setNode(node)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, lx)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
        fitter_node = fitter_node_iter.next()


def fit_trunk_model(modelfile, datafile, trunk_group_name=None):
    """
    Initialises scaffold fitter and runs through its steps for vagus nerve trunk.
    """
    fitter = Fitter(modelfile, datafile)
    fitter.load()

    # initial configuration
    fitter_fieldmodule = fitter.getFieldmodule()
    fitter.setModelCoordinatesFieldByName('coordinates')
    fitter.setDataCoordinatesFieldByName('coordinates')
    if trunk_group_name:
        fitter.setModelFitGroupByName(trunk_group_name)
    fitter.setFibreField(fitter_fieldmodule.findFieldByName("zero fibres"))
    fitter.setMarkerGroupByName('marker')  # not necessary, it's marker by default
    fitter.setDiagnosticLevel(0)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupDataWeight('marker', 200.0)
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight(None, 5.0)
    fit1.setNumberOfIterations(10)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    fit2 = FitterStepFit()
    fit2.setGroupDataWeight('marker', 400.0)
    fit1.setGroupDataWeight(None, None)
    fit2.setNumberOfIterations(5)
    fit2.setUpdateReferenceState(True)
    fitter.addFitterStep(fit2)

    # align step
    # align = FitterStepAlign()
    # align.setAlignMarkers(True)
    # align.setAlignGroups(True)
    # align.setScaleProportion(1.0)
    # fitter.addFitterStep(align)
    # align.run()

    # fit step 1
    # fit1 = FitterStepFit()
    # fit1.setGroupStrainPenalty(None, [15.0])
    # fit1.setGroupCurvaturePenalty(None, [50.0])
    # fit1.setGroupDataWeight('marker', 200.0)
    # fit1.setNumberOfIterations(5)
    # fit1.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit1)
    #
    # # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionError()
    rmsTrunkError, maxTrunkError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('left vagus X nerve trunk')
    rmsMarkerError, maxMarkerError = fitter.getDataRMSAndMaximumProjectionErrorForGroup('marker')

    # print('(all) RMS error: ' + str(rmsError))
    # print('(all) Max error: ' + str(maxError))
    # print('(trunk) RMS error: ' + str(rmsTrunkError))
    # print('(trunk) Max error: ' + str(maxTrunkError))
    # print('(marker) RMS error: ' + str(rmsMarkerError))
    # print('(marker) Max error: ' + str(maxMarkerError))

    return fitter


def fit_branches_model(modelfile, datafile, branch_name=None):
    """
    Initialises scaffold fitter and runs through its steps for a given vagus nerve branch.
    """
    fitter = Fitter(modelfile, datafile)
    fitter.load()

    # initial configuration
    fitter_fieldmodule = fitter.getFieldmodule()
    fitter.setModelCoordinatesFieldByName('coordinates')
    if branch_name:
        fitter.setModelFitGroupByName(branch_name)
    fitter.setFibreField(fitter_fieldmodule.findFieldByName("zero fibres"))
    fitter.setDataCoordinatesFieldByName('coordinates')
    fitter.setMarkerGroupByName('marker')  # not necessary, it's marker by default
    fitter.setDiagnosticLevel(0)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight('marker', 200.0)
    fit1.setNumberOfIterations(5)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionErrorForGroup(branch_name)
    # print('(%s) RMS error: %f' % (branch_name, rmsError))
    # print('(%s) Max error: %f' % (branch_name, maxError))

    return fitter


def get_fitted_vagus_nerve(region, marker_data, trunk_group_name, trunk_data_coordinates, branch_data, options):
    """

    """

    applyFitting = options['Apply fitting']
    addBranches = options['Add branches']
    add3D = options['Add 3D box']
    xmlData = False

    termNameVagusLengthList = {
        #"level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
        "level of superior border of the jugular foramen on the vagus nerve": 8.6342,
        "level of inferior border of the jugular foramen on the vagus nerve": 16.7227,
        "level of C1 transverse process on the vagus nerve": 32.1129,
        "level of angle of mandible on the vagus nerve": 42.2450,
        "level of greater horn of hyoid on the vagus nerve": 45.6122,
        "level of carotid bifurcation on the vagus nerve": 48.3581,
        "level of laryngeal prominence on the vagus nerve": 68.8431,
        "level of superior border of the clavicle on the vagus nerve": 117.5627,
        "level of jugular notch on the vagus nerve": 124.6407,
        "level of sternal angle on the vagus nerve": 151.2352,
        "1 cm superior to start of esophageal plexus on the vagus nerve": 165.5876,
        "level of esophageal hiatus on the vagus nerve": 254.32879,
        "level of aortic hiatus on the vagus nerve": 291.3695,
        #"level of end of trunk on Japanese dataset": 312.5  # note this term is also not on the list of annotations
    }

    lengthToDiameterRatio = 312.5  # calculated from total length of nerve/average diameter of nerve
    rescaledlengthToDiameterRatio = 100.0
    rescaledTermNameVagusLengthList = {}
    for term in termNameVagusLengthList:
        rescaledTermNameVagusLengthList[term] = termNameVagusLengthList[ term] \
                                                / lengthToDiameterRatio * rescaledlengthToDiameterRatio


    annotationGroups = []
    fieldmodule = region.getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]

    # node - geometric coordinates
    coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    for valueLabel in valueLabels[1:]:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 1)

    # branch special node
    nodeTemplateNoValue = nodes.createNodetemplate()
    nodeTemplateNoValue.defineField(coordinates)
    for valueLabel in valueLabels:
        if valueLabel == Node.VALUE_LABEL_VALUE:
            nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, valueLabel, 0)
        else:
            nodeTemplateNoValue.setValueNumberOfVersions(coordinates, -1, valueLabel, 1)

            # vagus centroid
            mesh = fieldmodule.findMeshByDimension(1)
            hermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft1d = mesh.createElementfieldtemplate(hermiteBasis)
            linetemplate = mesh.createElementtemplate()
            linetemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
            linetemplate.defineField(coordinates, -1, eft1d)

            # vagus centroid - branch root
            cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft1dNV = mesh.createElementfieldtemplate(hermiteBasis)
            eft1dNV.setNumberOfLocalNodes(4)
            eft1dNV.setNumberOfLocalScaleFactors(4)
            for i in range(4):
                eft1dNV.setScaleFactorType(i + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
            remapEftNodeValueLabelWithNodes(eft1dNV, 1, Node.VALUE_LABEL_VALUE,
                                            [(3, Node.VALUE_LABEL_VALUE, [1]),
                                             (3, Node.VALUE_LABEL_D_DS1, [2]),
                                             (4, Node.VALUE_LABEL_VALUE, [3]),
                                             (4, Node.VALUE_LABEL_D_DS1, [4])])
            linetemplateBranchRoot = mesh.createElementtemplate()
            linetemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
            linetemplateBranchRoot.defineField(coordinates, -1, eft1dNV)

    # build 1d trunk centroid line
    elementsAlongTrunk = options['Number of elements along the trunk']
    elementLength = rescaledlengthToDiameterRatio / (elementsAlongTrunk - 1)
    tx, td1, vtx, vtd1 = estimate_trunk_coordinates(elementsAlongTrunk, elementLength,
                                                    marker_data, rescaledTermNameVagusLengthList)

    trunkCentroidGroup = AnnotationGroup(region, (trunk_group_name, ""))
    annotationGroups.append(trunkCentroidGroup)
    trunkCentroidMeshGroup = trunkCentroidGroup.getMeshGroup(mesh)

    nodeIdentifier = 1
    lineIdentifier = 1
    for n in range(elementsAlongTrunk):
        sx = tx[n]
        sd1 = td1

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        fieldcache.setNode(node)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

        if n > 0:
            nids = [nodeIdentifier - 1, nodeIdentifier]
            line = mesh.createElement(lineIdentifier, linetemplate)
            line.setNodesByIdentifier(eft1d, nids)
            trunkCentroidMeshGroup.addElement(line)
            lineIdentifier += 1
        nodeIdentifier += 1

    # set markers
    for marker_name, marker_coordinate in marker_data.items():
        annotationGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
        annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
        nodeIdentifier += 1

        # create temporary model file
        sir = region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf")
        region.write(sir)

        fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
        fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf"
        fitter = fit_trunk_model(fitter_model_file, fitter_data_file)
        set_fitted_group_nodes(region, fitter)

        # branch building
        branch_root_params = {}
        for branch_name in branch_data.keys():
            print(branch_name)
            branch_coordinates = [branch_node[0] for branch_node in branch_data[branch_name]]

            # determine branch approximate start and closest trunk node index
            bx, bd1, parent_s_nid, parent_f_nid, parent_group_name, branch_root_xi, elementsAlongBranch = \
                estimate_branch_coordinates(region, branch_coordinates, elementLength)
            branch_root_params[branch_name] = [parent_s_nid, parent_f_nid, branch_root_xi, bx[0]]
            print(' parent: ', parent_group_name)

            branchCentroidGroup = AnnotationGroup(region, (branch_name, 'None'))
            annotationGroups.append(branchCentroidGroup)
            branchCentroidMeshGroup = branchCentroidGroup.getMeshGroup(mesh)

            for n in range(elementsAlongBranch):
                sx = bx[n]
                sd1 = bd1

                if n == 0:
                    # create branch special node
                    node = nodes.createNode(nodeIdentifier, nodeTemplateNoValue)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

                    if n == 1:
                        # create branch special element
                        nids = [nodeIdentifier - 1, nodeIdentifier,
                                parent_s_nid, parent_f_nid]
                        line = mesh.createElement(lineIdentifier, linetemplateBranchRoot)
                        line.setNodesByIdentifier(eft1dNV, nids)
                        scalefactorsNV = getCubicHermiteBasis(branch_root_xi)
                        line.setScaleFactors(eft1dNV, list(scalefactorsNV))
                        branchCentroidMeshGroup.addElement(line)
                        lineIdentifier += 1
                    else:
                        nids = [nodeIdentifier - 1, nodeIdentifier]
                        line = mesh.createElement(lineIdentifier, linetemplate)
                        line.setNodesByIdentifier(eft1d, nids)
                        branchCentroidMeshGroup.addElement(line)
                        lineIdentifier += 1
                nodeIdentifier += 1

                # remove trunk nodes from branch group
                parent_group = find_or_create_field_group(fieldmodule, parent_group_name)
                branchNodesetGroup = branchCentroidGroup.getNodesetGroup(nodes)
                if branchNodesetGroup.isValid():
                    branchNodesetGroup.removeNodesConditional(parent_group)

                sir = region.createStreaminformationRegion()
                srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_model.exf")
                region.write(sir)

                # geometry fitting - branches
                fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
                fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_model.exf"

                fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
                set_fitted_group_nodes(region, fitter, branch_name)

                # extract first branch node - x & d1 value
                branch_group = find_or_create_field_group(fieldmodule, branch_name)
                branch_nodes = branch_group.getNodesetGroup(nodes)
                node_iter = branch_nodes.createNodeiterator()
                node = node_iter.next()
                fieldcache.setNode(node)
                _, sd1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                branch_root_params[branch_name].append(sd1)






def load_data(region):
    """
    Extract data from supplied exf datafile, separate out data related to
    vagus trunk, vagus branches, fascicles, markers (anatomical landmarks)
    """

    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid()

    fm = data_region.getFieldmodule()
    fc = fm.createFieldcache()

    group_list = get_group_list(fm)
    group_map = {}
    trunk_group_name = None
    branch_group_names = []
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group
        if 'trunk' in group_name:
            trunk_group_name = group_name
        if 'branch' and ('left A') in group_name:
            branch_group_names.append(group_name)

    # extract marker data
    markers = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    marker_coordinates = fm.findFieldByName("marker_data_coordinates").castFiniteElement()
    marker_name_field = fm.findFieldByName("marker_data_name")
    assert marker_coordinates.isValid() and (marker_coordinates.getNumberOfComponents() == 3)

    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fc.setNode(marker_node)
            result, x = marker_coordinates.evaluateReal(fc, 3)
            marker_name = marker_name_field.evaluateString(fc)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    # extract other vagus data
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)
    radius = fm.findFieldByName("radius").castFiniteElement()

    trunk_data_coordinates = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name,
                                                               [Node.VALUE_LABEL_VALUE])
    if radius.isValid():
        trunk_radius = get_nodeset_fieldgroup_parameters(nodes, radius, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    else:
        trunk_radius = [2 for i in range(1, len(trunk_data_coordinates))]

    branch_data = {}
    branch_radius_data = {}
    for branch_name in branch_group_names:
        branch_parameters = get_nodeset_fieldgroup_parameters(nodes, coordinates, branch_name, [Node.VALUE_LABEL_VALUE])
        branch_data[branch_name] = branch_parameters
        if radius.isValid():
            branch_radius = get_nodeset_fieldgroup_parameters(nodes, radius, branch_name, [Node.VALUE_LABEL_VALUE])
        else:
            branch_radius = [1 for i in range(1, len(branch_parameters))]
        branch_radius_data[branch_name] = branch_radius

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    return marker_data, trunk_group_name, trunk_data_coordinates, trunk_radius, branch_data, branch_radius_data


def load_data_xml(region):
    """
    Extract data from supplied datafile (converted from xml file),
    separate out data related to vagus trunk, vagus branches, fascicles, markers (anatomical landmarks)
    """

    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid()

    fm = data_region.getFieldmodule()
    fc = fm.createFieldcache()

    group_list = get_group_list(fm)
    group_map = {}
    trunk_group_name = None
    branch_group_names = []
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group
        if 'trunk' in group_name:
            trunk_group_name = group_name
        if 'branch' and ('left A' or 'right A') in group_name:
            branch_group_names.append(group_name)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)
    radius = fm.findFieldByName("radius").castFiniteElement()

    markers = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    marker_names = fm.findFieldByName("marker_name")

    # extract data related to vagus markers, trunk and A-branches
    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fc.setNode(marker_node)
            _, x = coordinates.evaluateReal(fc, 3)
            marker_name = marker_names.evaluateString(fc)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    # extract other vagus data
    trunk_data_coordinates = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name,
                                                               [Node.VALUE_LABEL_VALUE])
    if radius.isValid():
        trunk_radius = get_nodeset_fieldgroup_parameters(nodes, radius, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    else:
        trunk_radius = [2 for i in range(1, len(trunk_data_coordinates))]

    branch_data = {}
    branch_radius_data = {}
    for branch_name in branch_group_names:
        branch_parameters = get_nodeset_fieldgroup_parameters(nodes, coordinates, branch_name, [Node.VALUE_LABEL_VALUE])
        branch_data[branch_name] = branch_parameters
        if radius.isValid():
            branch_radius = get_nodeset_fieldgroup_parameters(nodes, radius, branch_name, [Node.VALUE_LABEL_VALUE])
        else:
            branch_radius = [1 for i in range(1, len(branch_parameters))]
        branch_radius_data[branch_name] = branch_radius

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    return trunk_group_name, trunk_data_coordinates, trunk_radius, marker_data, branch_data

