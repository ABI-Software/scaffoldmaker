"""
Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches
"""

import math

from cmlibs.utils.zinc.field import find_or_create_field_group, findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field, FieldGroup
from cmlibs.zinc.node import Node

from cmlibs.maths.vectorops import add, cross, div, dot, magnitude, matrix_mult, matrix_inv, mult, rejection, \
    set_magnitude, sub

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.vagus_terms import get_vagus_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, \
    remapEftNodeValueLabelWithNodes, setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, getCubicHermiteBasisDerivatives, \
    interpolateCubicHermite, interpolateHermiteLagrange, smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters, set_nodeset_field_parameters, \
    make_nodeset_derivatives_orthogonal
from scaffoldmaker.utils.read_vagus_data import load_exf_data, load_exf_data_contours_japanese_dataset

from scaffoldfitter.fitter import Fitter
from scaffoldfitter.fitterstepalign import FitterStepAlign
from scaffoldfitter.fitterstepfit import FitterStepFit

class MeshType_3d_vagus_box1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @staticmethod
    def getName():
        return "3D Vagus Box 1"

    @staticmethod
    def getParameterSetNames():
        return [
            'Human Trunk 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk': 30,
            'Iterations (fit trunk)': 1,
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along the trunk',
            'Iterations (fit trunk)',
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options['Number of elements along the trunk'] < 10:
            options['Number of elements along the trunk'] = 10
        if options['Iterations (fit trunk)'] < 1:
            options['Iterations (fit trunk)'] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the 3d mesh for a vagus with branches, incorporating 1d central line,
        2d epineureum and 3d box based on vagus radius.

        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        # Zinc setup for vagus scaffold
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()

        # node - geometric coordinates
        valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                       Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                       Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for valueLabel in valueLabels[1:]:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, 1)

        # vagus box (3d)
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

        # vagus box (3d) branch root
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
        #                                         (Node.VALUE_LABEL_D2_DS1DS3ё s3)])
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

        # vagus centroid (1d)
        mesh1d = fieldmodule.findMeshByDimension(1)
        hermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft1d = mesh1d.createElementfieldtemplate(hermiteBasis)
        linetemplate = mesh1d.createElementtemplate()
        linetemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplate.defineField(coordinates, -1, eft1d)

        # vagus centroid (1d) branch root
        #cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft1dNV = mesh1d.createElementfieldtemplate(hermiteBasis)
        eft1dNV.setNumberOfLocalNodes(3)
        eft1dNV.setNumberOfLocalScaleFactors(24)
        for si in range(24):
            eft1dNV.setScaleFactorType(si + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
        si = 0
        for oldValueLabel in [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]:
            remapEftNodeValueLabelWithNodes(eft1dNV, 1, oldValueLabel,
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
        linetemplateBranchRoot = mesh1d.createElementtemplate()
        linetemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplateBranchRoot.defineField(coordinates, -1, eft1dNV)

        # vagus (2d) epineurium
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

        # vagus epineureum (2d) branch root
        # 1 & 3 - trunk nodes, 2 - branch 2nd element
        bicubichermiteSerendipityBasisNV = (
            fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
        # 4 elements around circle
        facetemplate_and_eft_list_BranchRoot = [None] * 4
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
            facetemplateBranchRoot = mesh2d.createElementtemplate()
            facetemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplateBranchRoot.defineField(coordinates, -1, eft2dNV)
            facetemplate_and_eft_list_BranchRoot[e] = (facetemplateBranchRoot, eft2dNV)


        elementsAlongTrunk = options['Number of elements along the trunk']

        # load data from file
        print('Extracting data...')
        data_region = region.getParent().findChildByName('data')
        if data_region.isValid():
            marker_data, trunk_group_name, trunk_data, trunk_radius, branch_data, branch_parents, branch_radius_data = \
                load_exf_data(data_region)
        assert len(marker_data) >= 2, f"At least two landmarks are expected in the data. Incomplete data."

        # evaluate & fit centroid lines for trunk and branches
        print('Building centerlines for scaffold...')
        fit_region, branches_order, branch_root_parameters = evaluate_vagus_1d_coordinates(
            region, marker_data, trunk_group_name, trunk_data, branch_data, branch_parents, options)
        fit_fieldmodule = fit_region.getFieldmodule()
        fit_fieldcache = fit_fieldmodule.createFieldcache()
        fit_coordinates = fit_fieldmodule.findFieldByName("coordinates").castFiniteElement()
        fit_nodes = fit_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        annotationGroups = []
        # vagus annotation groups
        vagusCentroidGroup = AnnotationGroup(region, ("Vagus centroid", ""))
        annotationGroups.append(vagusCentroidGroup)
        vagusCentroidMeshGroup = vagusCentroidGroup.getMeshGroup(mesh1d)

        vagusEpineuriumAnnotationGroup = AnnotationGroup(region, ("Vagus epineurium", ""))
        annotationGroups.append(vagusEpineuriumAnnotationGroup)
        vagusEpineuriumMeshGroup = vagusEpineuriumAnnotationGroup.getMeshGroup(mesh2d)

        node_map = {}
        print('Building trunk...')

        # read nodes
        fit_trunk_group = find_or_create_field_group(fit_fieldmodule, trunk_group_name)
        fit_trunk_nodes = fit_trunk_group.getNodesetGroup(fit_nodes)
        sn = []
        sx = []
        sd1 = []
        fit_node_iter = fit_trunk_nodes.createNodeiterator()
        fit_node = fit_node_iter.next()
        while fit_node.isValid():
            fit_fieldcache.setNode(fit_node)
            fit_node_id = fit_node.getIdentifier()
            _, tx = fit_coordinates.getNodeParameters(fit_fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            _, td1 = fit_coordinates.getNodeParameters(fit_fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            sn.append(fit_node_id)
            sx.append(tx)
            sd1.append(td1)
            fit_node = fit_node_iter.next()

        # calculate side and cross derivatives - d2, d3, d12, d13
        sd2, sd3 = set_group_nodes_derivatives_orthogonal(sd1)
        sd12, sd13 = smoothCurveSideCrossDerivatives(sx, sd1, [sd2, sd3])

        #  trunk annotation groups
        # TODO: use get_vagus_term for name when terms are finalised
        #trunkBoxGroup = AnnotationGroup(region, get_vagus_term(trunk_group_name.lower()))
        trunkBoxGroup = AnnotationGroup(region, (trunk_group_name, ''))
        annotationGroups.append(trunkBoxGroup)
        trunkBoxMeshGroup = trunkBoxGroup.getMeshGroup(mesh3d)

        nodeIdentifier = 1
        lineIdentifier = 1
        elementIdentifier = 1
        faceIdentifier = 1

        # create nodes
        for n in range(elementsAlongTrunk):
            node_map[sn[n]] = nodeIdentifier
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, sd2[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[n])
            #trunkCentroidNodesetGroup.addNode(node)
            nodeIdentifier += 1

        # create elements
        for n in range(1, elementsAlongTrunk):
            node_id = n + 1
            nids = [node_id - 1, node_id]

            line = mesh1d.createElement(lineIdentifier, linetemplate)
            line.setNodesByIdentifier(eft1d, nids)
            #trunkCentroidMeshGroup.addElement(line)
            vagusCentroidMeshGroup.addElement(line)
            lineIdentifier += 1

            element = mesh3d.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft3d, nids)
            element.setScaleFactors(eft3d, [-1.0])
            trunkBoxMeshGroup.addElement(element)
            elementIdentifier += 1

            for e in range(4):
                facetemplate, eft2d = facetemplate_and_eft_list[e]
                face = mesh2d.createElement(faceIdentifier, facetemplate)
                face.setNodesByIdentifier(eft2d, nids)
                face.setScaleFactors(eft2d, scalefactors2d)
                vagusEpineuriumMeshGroup.addElement(face)
                faceIdentifier += 1

        print('Building branches...')
        for branch_name in branches_order:
            #print(' ', branch_name)

            # read nodes
            fit_branch_group = find_or_create_field_group(fit_fieldmodule, branch_name)
            fit_branch_nodes = fit_branch_group.getNodesetGroup(fit_nodes)
            elementsAlongBranch = fit_branch_nodes.getSize()
            sn = []
            sx = []
            sd1 = []
            fit_node_iter = fit_branch_nodes.createNodeiterator()
            fit_node = fit_node_iter.next()  # ignore branch root node
            fit_node = fit_node_iter.next()
            while fit_node.isValid():
                fit_fieldcache.setNode(fit_node)
                fit_node_id = fit_node.getIdentifier()
                _, tx = fit_coordinates.getNodeParameters(fit_fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, td1 = fit_coordinates.getNodeParameters(fit_fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                sn.append(fit_node_id)
                sx.append(tx)
                sd1.append(td1)
                fit_node = fit_node_iter.next()

            # calculate side and cross derivatives - d2, d3, d12, d13
            sd2, sd3 = set_group_nodes_derivatives_orthogonal(sd1)
            sd12, sd13 = smoothCurveSideCrossDerivatives(sx, sd1, [sd2, sd3])

            # branch root parameters (assume d12, d13 are both zero for now)
            trunk_segment_start_id = node_map[branch_root_parameters[branch_name][0]]
            trunk_segment_end_id = node_map[branch_root_parameters[branch_name][1]]
            branch_root_xi = branch_root_parameters[branch_name][2]
            bx = branch_root_parameters[branch_name][3]
            bd1 = branch_root_parameters[branch_name][4]
            bd2, bd3 = set_group_nodes_derivatives_orthogonal([bd1])

            # trunk interpolation
            fns = list(getCubicHermiteBasis(branch_root_xi))  # for x, d2, d3
            dfns = list(getCubicHermiteBasisDerivatives(branch_root_xi))  # for d1, d12, d13

            # access two trunk segments to get interpolated values
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

            # branch annotation groups
            # TODO: use get_vagus_term for name when terms are finalised
            # branchBoxGroup = AnnotationGroup(region, get_vagus_term(branch_name.lower()))
            branchBoxGroup = AnnotationGroup(region, (branch_name, ""))
            annotationGroups.append(branchBoxGroup)
            branchBoxMeshGroup = branchBoxGroup.getMeshGroup(mesh3d)

            # create nodes (excluding branch root)
            for n in range(elementsAlongBranch - 1):
                node_map[sn[n]] = nodeIdentifier
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, sd2[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[n])

                if n == 0:
                    # branch root element
                    nids = [trunk_segment_start_id, nodeIdentifier, trunk_segment_end_id]

                    scalefactors = []
                    scalefactors.extend(scalefactorsX)
                    scalefactors.extend(scalefactorsD1)

                    line = mesh1d.createElement(lineIdentifier, linetemplateBranchRoot)
                    line.setNodesByIdentifier(eft1dNV, nids)
                    line.setScaleFactors(eft1dNV, list(scalefactors))
                    vagusCentroidMeshGroup.addElement(line)
                    lineIdentifier += 1

                else:
                    # other elements
                    nids = [nodeIdentifier - 1, nodeIdentifier]

                    line = mesh1d.createElement(lineIdentifier, linetemplate)
                    line.setNodesByIdentifier(eft1d, nids)
                    vagusCentroidMeshGroup.addElement(line)
                    lineIdentifier += 1

                nodeIdentifier += 1


            # create elements
            for n in range(0, elementsAlongBranch - 1):
                node_id = node_map[sn[n]]
                if n == 0:
                    # branch root element
                    nids = [trunk_segment_start_id, node_id, trunk_segment_end_id]

                    scalefactors = [-1]
                    scalefactors.extend(sub(sub(scalefactorsX, scalefactorsD2), scalefactorsD3))
                    scalefactors.extend(sub(sub(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                    scalefactors.extend(sub(add(scalefactorsX, scalefactorsD2), scalefactorsD3))
                    scalefactors.extend(sub(add(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                    scalefactors.extend(add(sub(scalefactorsX, scalefactorsD2), scalefactorsD3))
                    scalefactors.extend(add(sub(scalefactorsD1, scalefactorsD12), scalefactorsD13))
                    scalefactors.extend(add(add(scalefactorsX, scalefactorsD2), scalefactorsD3))
                    scalefactors.extend(add(add(scalefactorsD1, scalefactorsD12), scalefactorsD13))

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

                        facetemplateBranchRoot, eft2dNV = facetemplate_and_eft_list_BranchRoot[e]
                        face = mesh2d.createElement(faceIdentifier, facetemplateBranchRoot)
                        face.setNodesByIdentifier(eft2dNV, nids)
                        face.setScaleFactors(eft2dNV, scalefactors)
                        vagusEpineuriumMeshGroup.addElement(face)
                        faceIdentifier += 1

                else:
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
                        vagusEpineuriumMeshGroup.addElement(face)
                        faceIdentifier += 1

            # remove trunk nodes from branch group
            parentFieldGroup = find_or_create_field_group(fieldmodule, branch_parents[branch_name])
            branchNodesetGroup = branchBoxGroup.getNodesetGroup(nodes)
            if branchNodesetGroup.isValid():
                branchNodesetGroup.removeNodesConditional(parentFieldGroup)

        # set markers
        print('Adding anatomical landmarks...')
        for marker_name, marker_coordinate in marker_data.items():
            # TODO: use get_vagus_term for name when terms are finalised
            # annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_vagus_term(marker_name), isMarker=True)
            annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, (marker_name,''), isMarker=True)

            annotationGroup.createMarkerNode(nodeIdentifier, coordinates, marker_coordinate)
            nodeIdentifier += 1

        print('Adding material coordinates...')

        # add material coordinates
        coordinates.setName("vagus coordinates")  # temporarily rename
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        region.write(sir)
        result, buffer = srm.getBuffer()
        coordinates.setName("coordinates")  # restore name before reading vagus coordinates back in
        sir = region.createStreaminformationRegion()
        sir.createStreamresourceMemoryBuffer(buffer)
        region.read(sir)  # read and merge with region, thus having coordinates and vagus coordinates in the region together
        vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()

        rescaledVagusTrunkLength = 1.0
        vagusAspectRatio = 0.005  # vagus approx diameter (5mm) / vagus length (85mm)
        vagusRadius = vagusAspectRatio * rescaledVagusTrunkLength / 2

        # calculate derivatives
        derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
        derivative_xi2 = mesh3d.getChartDifferentialoperator(1, 2)
        derivative_xi3 = mesh3d.getChartDifferentialoperator(1, 3)

        zero = [0.0, 0.0, 0.0]

        # vagus trunk goes along z axis with origin as top of the vagus and rescaledVagusTrunkLength as bottom of the trunk
        elem_iter = mesh3d.createElementiterator()
        element = elem_iter.next()
        while element.isValid():
            # trunk elements first,
            # followed by branch elements (with first element having more than 2 local nodes)

            element_id = element.getIdentifier()
            eft = element.getElementfieldtemplate(coordinates, -1)
            local_nodes_count = eft.getNumberOfLocalNodes()
            if local_nodes_count == 2:
                if element_id == 1:
                    # first trunk element
                    x = [0.0, 0.0, 0.0]
                    d1 = [0, 0, rescaledVagusTrunkLength / (elementsAlongTrunk - 1)]
                    [d2], [d3] = set_group_nodes_derivatives_orthogonal([d1])
                    d2 = set_magnitude(d2, vagusRadius)
                    d3 = set_magnitude(d3, vagusRadius)

                    node = element.getNode(eft, 1)
                    fieldcache.setNode(node)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)

                x = add(x, d1)

                node = element.getNode(eft, 2)
                fieldcache.setNode(node)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)

            elif local_nodes_count > 2:
                # first branch element
                fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5])  # assuming xi1 is along the branch
                _, x = vagus_coordinates.evaluateReal(fieldcache, 3)
                _, d1 = vagus_coordinates.evaluateDerivative(derivative_xi1, fieldcache, 3)
                _, d2 = vagus_coordinates.evaluateDerivative(derivative_xi2, fieldcache, 3)
                _, d3 = vagus_coordinates.evaluateDerivative(derivative_xi3, fieldcache, 3)

                #[d2], [d3] = set_group_nodes_derivatives_orthogonal([d1])
                d2 = set_magnitude(d2, vagusRadius)
                d3 = set_magnitude(d3, vagusRadius)
                x = add(x, d1)

                node = element.getNode(eft, 2)
                fieldcache.setNode(node)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                vagus_coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)

            element = elem_iter.next()

        print('Done')
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

        vagusEpineuriumAnnotationGroup = getAnnotationGroupForTerm(annotationGroups, ("Vagus epineurium", ""))
        vagusEpineuriumMeshGroup = vagusEpineuriumAnnotationGroup.getMeshGroup(mesh2d)
        vagusAnteriorLineAnnotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups,
                                                                              region, ("Vagus anterior line", ""))
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
def magnitude_squared(v):
    '''
    return: squared scalar magnitude of vector v
    '''

    # TODO: proposed function to cmlibs.maths
    return sum(c * c for c in v)

def find_dataset_endpoints(coordinate_dataset):
    """
    Given list of node coordinates, find two furthest from each other nodes.
    Returns a list of two node coordinates.
    """

    max_distance_squared = 0
    ends_points = []
    for node_x in coordinate_dataset:
        for node_y in coordinate_dataset:
            distance_squared = magnitude_squared(sub(node_x, node_y))
            if distance_squared >= max_distance_squared:
                max_distance_squared = distance_squared
                ends_points = [node_y, node_x]
    return ends_points


def setSideCrossDerivatives(nodeset, field, isBranch = False):
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
        td2 = rejection(yx, d1[i])
        td2 = set_magnitude(td2, magnitude(yx)) # change magnitude later with radius
        d2.append(td2)

        td3 = cross(d1[i], td2)
        td3 = set_magnitude(td3, magnitude(zx)) # change magnitude later with radius
        d3.append(td3)

    return d2, d3


def evaluate_vagus_1d_coordinates(region, marker_data, trunk_group_name, trunk_data, branch_data, branch_parents,
                                  options):
    """
    Generates an intermediate hermite x bilinear 1-D central line mesh for a vagus nerve with branches, markers.
    First, it evaluates trunk coordinates based on top and bottom of the supplied markers, then fits the coordinates to
    supplied trunk data.
    Second, it creates branches based on the
    """

    fit_region = region.createRegion()

    elementsAlongTrunk = options['Number of elements along the trunk']
    iterationsNumber = options['Iterations (fit trunk)']

    fit_fm = fit_region.getFieldmodule()
    fit_fc = fit_fm.createFieldcache()
    fit_nodes = fit_fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    fit_coordinates = findOrCreateFieldCoordinates(fit_fm).castFiniteElement()

    # 1d centroid line
    fit_nodetemplate = fit_nodes.createNodetemplate()
    fit_nodetemplate.defineField(fit_coordinates)
    value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
    for value_label in value_labels:
        fit_nodetemplate.setValueNumberOfVersions(fit_coordinates, -1, value_label, 1)

    fit_mesh1d = fit_fm.findMeshByDimension(1)
    fit_hermiteBasis = fit_fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    fit_eft1d = fit_mesh1d.createElementfieldtemplate(fit_hermiteBasis)
    fit_eltemplate = fit_mesh1d.createElementtemplate()
    fit_eltemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
    fit_eltemplate.defineField(fit_coordinates, -1, fit_eft1d)

    # trunk group
    trunkCentroidGroup = find_or_create_field_group(fit_fm, trunk_group_name)
    trunkCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
    trunkCentroidMeshGroup = trunkCentroidGroup.getOrCreateMeshGroup(fit_mesh1d)

    # used for fitting only
    trunkFitCentroidGroup = find_or_create_field_group(fit_fm, trunk_group_name + '-fit')
    trunkFitCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
    trunkFitCentroidMeshGroup = trunkFitCentroidGroup.getOrCreateMeshGroup(fit_mesh1d)

    for ii in range(iterationsNumber):
        annotationGroups = []

        if ii == 0:
            tx, td1, elementLength = estimate_trunk_coordinates(elementsAlongTrunk, marker_data)
        else:
            # read tx from fit_coordinates
            _, node_field_parameters = get_nodeset_field_parameters(fit_nodes, fit_coordinates, value_labels)
            tx = [nodeParameter[1][0][0] for nodeParameter in node_field_parameters]
            td1 = [nodeParameter[1][1][0] for nodeParameter in node_field_parameters]

        trunk_nodes_data_bounds = estimate_trunk_data_boundaries(tx, trunk_data, elementsAlongTrunk)

        nodeIdentifier = 1
        lineIdentifier = 1
        nodes_before = []
        nodes_after = []
        for n in range(elementsAlongTrunk):
            sx = tx[n]
            sd1 = td1[n]

            if ii == 0:
                node = fit_nodes.createNode(nodeIdentifier, fit_nodetemplate)
            else:
                node = fit_nodes.findNodeByIdentifier(nodeIdentifier)
            fit_fc.setNode(node)
            fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, sx)
            fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

            # add to trunk group used for data fitting
            if trunk_nodes_data_bounds[0] <= nodeIdentifier <= trunk_nodes_data_bounds[-1]:
                pass
            elif nodeIdentifier < trunk_nodes_data_bounds[0]:
                nodes_before.append(nodeIdentifier)
            else:
                nodes_after.append(nodeIdentifier)

            if n > 0:
                nids = [nodeIdentifier - 1, nodeIdentifier]
                if ii == 0:
                    line = fit_mesh1d.createElement(lineIdentifier, fit_eltemplate)
                else:
                    line = fit_mesh1d.findElementByIdentifier(lineIdentifier)
                line.setNodesByIdentifier(fit_eft1d, nids)
                trunkCentroidMeshGroup.addElement(line)
                # add element to trunk group used for data fitting
                if nodeIdentifier - 1 >= trunk_nodes_data_bounds[0] and nodeIdentifier <= trunk_nodes_data_bounds[-1]:
                    trunkFitCentroidMeshGroup.addElement(line)
                lineIdentifier += 1
            nodeIdentifier += 1

        if ii == 0:
            # set markers
            for marker_name, marker_coordinate in marker_data.items():
                # TODO: use get_vagus_term for name when terms are finalised
                #annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, fit_region, get_vagus_term(marker_name), isMarker=True)
                annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, fit_region, (marker_name, ''), isMarker=True)

                annotationGroup.createMarkerNode(nodeIdentifier, fit_coordinates, marker_coordinate)
                nodeIdentifier += 1
        else:
            nodeIdentifier += len(marker_data)

        # create temporary model file
        sir = fit_region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf")
        fit_region.write(sir)

        print('... Fitting trunk, iteration', str(ii + 1))
        fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
        fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_trunk_model.exf"
        fitter = fit_trunk_model(fitter_model_file, fitter_data_file, trunk_group_name + '-fit')
        set_fitted_group_nodes(fit_region, fitter, trunk_group_name + '-fit')

        if len(nodes_before) > 0:
            # recalculate unfitted nodes by the first fitted node
            node = fit_nodes.findNodeByIdentifier(nodes_before[-1] + 1)
            fit_fc.setNode(node)
            _, lx = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            _, ld1 = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

            node_count = 1
            for i in range(len(nodes_before) - 1, -1, -1):
                node_id = nodes_before[i]
                x = [lx[j] - node_count * ld1[j] for j in range(3)]

                node = fit_nodes.findNodeByIdentifier(node_id)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, x)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                node_count += 1

        if len(nodes_after) > 0:
            # recalculate unfitted nodes by the last fitted node
            node = fit_nodes.findNodeByIdentifier(nodes_after[0] - 1)
            fit_fc.setNode(node)
            _, lx = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            _, ld1 = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, 3)

            node_count = 1
            for i in range(len(nodes_after)):
                node_id = nodes_after[i]
                x = [lx[j] + node_count * ld1[j] for j in range(3)]

                node = fit_nodes.findNodeByIdentifier(node_id)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, x)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, ld1)
                node_count += 1

    visited_branches = []
    branch_root_parameters = {}
    # 1d centroid line - branch root node
    fit_nodetemplateBranchRoot = fit_nodes.createNodetemplate()
    fit_nodetemplateBranchRoot.defineField(fit_coordinates)
    for value_label in value_labels:
        if value_label == Node.VALUE_LABEL_VALUE:
            fit_nodetemplateBranchRoot.setValueNumberOfVersions(fit_coordinates, -1, value_label, 0)
        else:
            fit_nodetemplateBranchRoot.setValueNumberOfVersions(fit_coordinates, -1, value_label, 1)

    # branch root element
    fit_cubicHermiteBasis = fit_fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    fit_eft1dBranchRoot = fit_mesh1d.createElementfieldtemplate(fit_cubicHermiteBasis)
    fit_eft1dBranchRoot.setNumberOfLocalNodes(4)
    fit_eft1dBranchRoot.setNumberOfLocalScaleFactors(4)
    for i in range(4):
        fit_eft1dBranchRoot.setScaleFactorType(i + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_ELEMENT_GENERAL)
    remapEftNodeValueLabelWithNodes(fit_eft1dBranchRoot, 1, Node.VALUE_LABEL_VALUE,
                                    [(3, Node.VALUE_LABEL_VALUE, [1]),
                                     (3, Node.VALUE_LABEL_D_DS1, [2]),
                                     (4, Node.VALUE_LABEL_VALUE, [3]),
                                     (4, Node.VALUE_LABEL_D_DS1, [4])])
    fit_eltemplateBranchRoot = fit_mesh1d.createElementtemplate()
    fit_eltemplateBranchRoot.setElementShapeType(Element.SHAPE_TYPE_LINE)
    fit_eltemplateBranchRoot.defineField(fit_coordinates, -1, fit_eft1dBranchRoot)

    print('... Adding branches')
    queue = [branch for branch in branch_parents.keys() if branch_parents[branch] == trunk_group_name]
    while queue:
        branch_name = queue.pop(0)
        print(branch_name)

        if branch_name in visited_branches:
            continue
        visited_branches.append(branch_name)
        queue.extend([branch for branch in branch_parents.keys() if branch_parents[branch] == branch_name])

        branch_coordinates = [branch_node[0] for branch_node in branch_data[branch_name]]
        branch_parent_name = branch_parents[branch_name]
        #print('  parent: ', branch_parent_name)

        # determine branch approximate start and closest trunk node index
        bx, bd1, parent_s_nid, parent_f_nid, branch_root_xi, elementsAlongBranch = \
            estimate_branch_coordinates(fit_region, branch_coordinates, elementLength, branch_parent_name)
        branch_root_parameters[branch_name] = [parent_s_nid, parent_f_nid, branch_root_xi, bx[0]]
        # print('  branch between nodes: ', parent_s_nid, parent_f_nid)

        # branchCentroidGroup = AnnotationGroup(fit_region, (branch_name, 'None'))
        # annotationGroups.append(branchCentroidGroup)
        # branchCentroidMeshGroup = branchCentroidGroup.getMeshGroup(fit_mesh1d)

        branchCentroidGroup = find_or_create_field_group(fit_fm, branch_name)
        branchCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
        branchCentroidMeshGroup = branchCentroidGroup.getOrCreateMeshGroup(fit_mesh1d)

        for n in range(elementsAlongBranch):
            sx = bx[n]
            sd1 = bd1

            if n == 0:
                # create branch special node
                node = fit_nodes.createNode(nodeIdentifier, fit_nodetemplateBranchRoot)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
            else:
                node = fit_nodes.createNode(nodeIdentifier, fit_nodetemplate)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

                if n == 1:
                    # create branch root element
                    nids = [nodeIdentifier - 1, nodeIdentifier,
                            parent_s_nid, parent_f_nid]
                    line = fit_mesh1d.createElement(lineIdentifier, fit_eltemplateBranchRoot)
                    line.setNodesByIdentifier(fit_eft1dBranchRoot, nids)
                    scalefactorsNV = getCubicHermiteBasis(branch_root_xi)
                    line.setScaleFactors(fit_eft1dBranchRoot, list(scalefactorsNV))
                    branchCentroidMeshGroup.addElement(line)
                    lineIdentifier += 1
                else:
                    nids = [nodeIdentifier - 1, nodeIdentifier]
                    line = fit_mesh1d.createElement(lineIdentifier, fit_eltemplate)
                    line.setNodesByIdentifier(fit_eft1d, nids)
                    branchCentroidMeshGroup.addElement(line)
                    lineIdentifier += 1
            nodeIdentifier += 1

        # remove trunk nodes from branch group
        parentGroup = find_or_create_field_group(fit_fm, branch_parent_name)
        branchNodesetGroup = branchCentroidGroup.getNodesetGroup(fit_nodes)
        if branchNodesetGroup.isValid():
            branchNodesetGroup.removeNodesConditional(parentGroup)

        sir = fit_region.createStreaminformationRegion()
        srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_model.exf")
        fit_region.write(sir)

        #print('  ... branch fitting')
        fitter_data_file = "C:/MAP/output/vagus_scaffold_temp/vagus_data.exf"
        fitter_model_file = "C:/MAP/output/vagus_scaffold_temp/vagus_model.exf"
        fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
        set_fitted_group_nodes(fit_region, fitter, branch_name)

        # extract first branch node - x & d1 fitted value
        branch_group = find_or_create_field_group(fit_fm, branch_name)
        branch_nodes = branch_group.getNodesetGroup(fit_nodes)
        node_iter = branch_nodes.createNodeiterator()
        node = node_iter.next()
        fit_fc.setNode(node)
        _, sd1 = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        branch_root_parameters[branch_name].append(sd1)

    return fit_region, visited_branches, branch_root_parameters


def estimate_trunk_coordinates(elementsAlongTrunk, marker_data):
    """

    """

    # choose markers for building initial scaffold
    # at the moment uses the first and the last markers in the data
    termNameVagusLengthList = {
        "level of exiting brainstem": 0.0,  # note this term is not on the list of annotations
        "level of superior border of jugular foramen": 8.6342,
        "level of inferior border of jugular foramen": 16.7227,
        "level of C1 transverse process": 32.1129,
        "level of angle of mandible": 42.2450,
        "level of greater horn of hyoid": 45.6122,
        "level of carotid bifurcation": 48.3581,
        "level of laryngeal prominence": 68.8431,
        "level of superior border of the clavicle": 117.5627,
        "level of jugular notch": 124.6407,
        "level of carina": 149.5929, # not on the list of annotations yet!
        "level of sternal angle": 151.2352,
        "level of 1 cm superior to start of esophageal plexus": 165.5876,
        "level of esophageal hiatus": 254.32879,
        "level of aortic hiatus": 291.3695,
        "level of end of trunk": 312.5  # note this term is also not on the list of annotations
    }

    totalVagusLength = 312.5  # calculated from total length of nerve/average diameter of nerve

    use_markers = [list(marker_data.keys())[0],
                   list(marker_data.keys())[-1]]

    pts = []
    params = []

    for marker in use_markers:
        use_marker_name = marker.replace('left ', '', 1).replace('right ', '', 1).replace(' on the vagus nerve', '', 1)
        assert use_marker_name in termNameVagusLengthList

        pts.append(marker_data[marker])
        params.append(termNameVagusLengthList[use_marker_name])

    step = totalVagusLength / (elementsAlongTrunk - 1)
    dx, dy, dz = [(pts[1][dim] - pts[0][dim]) / (params[1] - params[0]) for dim in range(3)]

    trunk_coords = []
    trunk_d1 = []
    for i in range(elementsAlongTrunk):
        trunk_coords.append([pts[0][0] + dx * (i * step - params[0]),
                             pts[0][1] + dy * (i * step - params[0]),
                             pts[0][2] + dz * (i * step - params[0])])
        trunk_d1.append([dx * step, dy * step, dz * step])

    return trunk_coords, trunk_d1, step


def estimate_trunk_data_boundaries(trunk_nodes, trunk_data, elementsAlongTrunk):
    """
    """

    # finding trunk nodes at the start and end of trunk data - by projections
    # TODO: change brute force to something (takes too long for CASE data),
    # at the moment assumes trunk_data is numbered top to bottom

    #trunk_data_endpoints = find_dataset_endpoints([trunk_pt[0] for trunk_pt in trunk_data])
    trunk_data_endpoints = [trunk_data[0][0], trunk_data[-1][0]]

    trunk_nodes_data_bounds = []
    for ind, endpoint in enumerate(trunk_data_endpoints):
        ap = sub(endpoint, trunk_nodes[0])
        ab = sub(trunk_nodes[-1], trunk_nodes[0])
        param = dot(ap, ab) / dot(ab, ab)  # between 0 and 1
        nearby_node = param * (elementsAlongTrunk - 1) + 1

        if ind == 0:
            # trunk node near the start of data
            trunk_nodes_data_bounds.append(math.floor(nearby_node))
        else:
            # trunk node near the end of data
            if nearby_node >= elementsAlongTrunk:
                trunk_nodes_data_bounds.append(elementsAlongTrunk)
            else:
                trunk_nodes_data_bounds.append(math.ceil(nearby_node))

    return trunk_nodes_data_bounds


def estimate_branch_coordinates(region, branch_coordinates, elementLength, branch_parent_name):
    """

    """

    fm = region.getFieldmodule()
    fieldcache = fm.createFieldcache()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    # assumes parent_group_name is known
    branch_start_x, branch_end_x, parent_s_nid, parent_f_nid = find_branch_start_segment(region,
                                                                                         branch_coordinates,
                                                                                         branch_parent_name)

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

    branch_length = magnitude(sub(branch_end_x, branch_start_x))
    elementsAlongBranch = int(branch_length / elementLength - 1)
    if elementsAlongBranch < 3:
        # need to have at least 3 nodes (including fictionary), otherwise branch
        elementsAlongBranch = 3

    branch_coordinates = []
    dx, dy, dz = div(sub(branch_end_x, branch_start_x), (elementsAlongBranch - 1))
    for i in range(elementsAlongBranch):
        branch_coordinates.append([branch_start_x[0] + dx * i,
                                   branch_start_x[1] + dy * i,
                                   branch_start_x[2] + dz * i])

    return branch_coordinates, [dx, dy, dz], parent_s_nid, parent_f_nid, branch_root_xi, elementsAlongBranch


def find_branch_start_segment(region, branch_coordinates, parent_group_name):
    """

    """
    branch_ends_points = find_dataset_endpoints(branch_coordinates)

    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    parent_group = find_or_create_field_group(fm, parent_group_name)
    parent_nodeset = parent_group.getNodesetGroup(nodes)

    _, group_parameters = get_nodeset_field_parameters(parent_nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
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

    parent_s_node_id = group_ids[closest_index]
    parent_f_node_id = group_ids[closest_index + 1]
    branch_end = branch_ends_points[0] if branch_ends_points[1] == branch_start else branch_ends_points[1]

    return branch_start, branch_end, parent_s_node_id, parent_f_node_id


# fitter functions
def fit_trunk_model(modelfile, datafile, trunk_group_name = None):
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

    # align step
    # align = FitterStepAlign()
    # align.setAlignMarkers(True)
    # align.setAlignGroups(True)
    # align.setScaleProportion(1.0)
    # fitter.addFitterStep(align)
    # align.run()

    # fit step 1
    # fit1 = FitterStepFit()
    # fit1.setGroupDataWeight('marker', 200.0)
    # fit1.setGroupStrainPenalty(None, [15.0])
    # fit1.setGroupCurvaturePenalty(None, [50.0])
    # fit1.setGroupDataWeight(None, 5.0)
    # fit1.setNumberOfIterations(10)
    # fit1.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit1)
    #
    # # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupDataWeight('marker', 400.0)
    # fit1.setGroupDataWeight(None, None)
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    # fit step 1
    fit1 = FitterStepFit()
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setGroupDataWeight(None, 10.0)
    fit1.setNumberOfIterations(10)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    fit2 = FitterStepFit()
    fit1.setGroupDataWeight(None, [10.0])
    fit1.setGroupStrainPenalty(None, [5.0])
    fit1.setGroupCurvaturePenalty(None, [100.0])
    fit2.setNumberOfIterations(5)
    fit2.setUpdateReferenceState(True)
    fitter.addFitterStep(fit2)

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


def fit_branches_model(modelfile, datafile, branch_name = None):
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
    fit1.setGroupDataWeight(None, 5.0)
    fit1.setNumberOfIterations(5)
    fit1.setUpdateReferenceState(True)
    fitter.addFitterStep(fit1)

    # fit step 2
    # fit2 = FitterStepFit()
    # fit2.setGroupStrainPenalty(None, [50.0])
    # fit2.setGroupCurvaturePenalty(None, [10.0])
    # fit2.setNumberOfIterations(5)
    # fit2.setUpdateReferenceState(True)
    # fitter.addFitterStep(fit2)

    fitter.run()

    rmsError, maxError = fitter.getDataRMSAndMaximumProjectionErrorForGroup(branch_name)
    #print('(%s) RMS error: %f' % (branch_name, rmsError))
    #print('(%s) Max error: %f' % (branch_name, maxError))

    return fitter


def set_fitted_group_nodes(region, fitter_region, group_name = None):
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

