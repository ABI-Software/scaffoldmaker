import os
import math
import logging
import tempfile

from cmlibs.maths.vectorops import (
    add, cross, distance_squared, div, dot, magnitude, matrix_mult, matrix_inv, mult, normalize, rejection,
    set_magnitude, sub)
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.field import (
    find_or_create_field_stored_mesh_location, find_or_create_field_group, find_or_create_field_coordinates,
    orphan_field_by_name)
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field, FieldFindMeshLocation, FieldGroup
from cmlibs.zinc.node import Node
from scaffoldfitter.fitter import Fitter as GeometryFitter
from scaffoldfitter.fitterstepfit import FitterStepFit
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    findAnnotationGroupByName
from scaffoldmaker.annotation.vagus_terms import get_vagus_branch_term, get_vagus_marker_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, remapEftNodeValueLabelWithNodes, \
    setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import (
    computeLagrangeHermiteDerivative, evaluateCoordinatesOnCurve, evaluateScalarOnCurve, getCubicHermiteBasis,
    getCubicHermiteBasisDerivatives, getCubicHermiteArcLength, getCubicHermiteCurvesLength,
    getCubicHermiteTrimmedCurvesLengths, getNearestLocationOnCurve, get_curve_from_points, interpolateCubicHermite,
    interpolateCubicHermiteDerivative, interpolateHermiteLagrange, interpolateLagrangeHermiteDerivative,
    sampleCubicHermiteCurves, sampleCubicHermiteCurvesSmooth, smoothCurveSideCrossDerivatives,
    track_curve_side_direction)
from scaffoldmaker.utils.read_vagus_data import load_vagus_data
from scaffoldmaker.utils.zinc_utils import (
    define_and_fit_field, find_or_create_field_zero_fibres, fit_hermite_curve, generate_curve_mesh, generate_datapoints,\
    generate_mesh_marker_points, get_nodeset_field_parameters)


logger = logging.getLogger(__name__)


def get_left_vagus_marker_locations_list():
    # vagus markers location in material coordinates between 0 to 1
    left_termNameVagusLengthList = {
        # cervical region
        "left level of superior border of jugular foramen on the vagus nerve": 0.02737296,
        "left level of inferior border of jugular foramen on the vagus nerve": 0.04434952,
        # "left level of inferior border of cranium on the vagus nerve": 0.0588,
        # "left level of C1 transverse process on the vagus nerve": 0.10276128,
        "left level of angle of the mandible on the vagus nerve": 0.12533074,
        # "left level of greater horn of hyoid on the vagus nerve": 0.14595904,
        "left level of carotid bifurcation on the vagus nerve": 0.15738364,
        "left level of laryngeal prominence on the vagus nerve": 0.20541934,
        # thoracic region
        "left level of superior border of the clavicle on the vagus nerve": 0.33847976,
        "left level of jugular notch on the vagus nerve": 0.38062311,
        "left level of sternal angle on the vagus nerve": 0.48395264,
        # "left level of 1 cm superior to start of esophageal plexus on the vagus nerve": 0.52988032,
        # abdominal region
        # "left level of esophageal hiatus on the vagus nerve": 0.813852428,
        # "left level of aortic hiatus on the vagus nerve": 0.9323824,
        # "left level of end of trunk": 1.0  # note this term is also not on the list of annotations
    }
    return left_termNameVagusLengthList


def get_right_vagus_marker_locations_list():
    # vagus markers location in material coordinates between 0 to 1
    right_termNameVagusLengthList = {
        # cervical region
        "right level of superior border of jugular foramen on the vagus nerve": 0.02762944,
        "right level of inferior border of jugular foramen on the vagus nerve": 0.04434952,
        # "right level of inferior border of cranium on the vagus nerve": 0.0588,
        # "right level of C1 transverse process on the vagus nerve": 0.10276128,
        "right level of angle of the mandible on the vagus nerve": 0.12648368,
        # "right level of greater horn of hyoid on the vagus nerve": 0.14595904,
        "right level of carotid bifurcation on the vagus nerve": 0.17798550,
        "right level of laryngeal prominence on the vagus nerve": 0.23144827,
        # thoracic region
        "right level of superior border of the clavicle on the vagus nerve": 0.33948916,
        "right level of jugular notch on the vagus nerve": 0.38937585,
        "right level of sternal angle on the vagus nerve": 0.48764507,
        # "right level of 1 cm superior to start of esophageal plexus on the vagus nerve": 0.52988032,
        # abdominal region
        # "right level of esophageal hiatus on the vagus nerve": 0.813852428,
        # "right level of aortic hiatus on the vagus nerve": 0.9323824,
    }
    return right_termNameVagusLengthList


def is_bony_landmark(marker_name):
    """
    Checks if supplied marker_name is a bony landmark.
    Not used currently.
    """
    bony_landmarks_names = [
        "level of superior border of jugular foramen on the vagus nerve",
        "level of inferior border of jugular foramen on the vagus nerve",
        "level of inferior border of cranium on the vagus nerve",
        "level of C1 transverse process on the vagus nerve",
        "level of angle of the mandible on the vagus nerve",
        "level of greater horn of hyoid on the vagus nerve",
        "level of superior border of the clavicle on the vagus nerve",
        "level of jugular notch on the vagus nerve",
        "level of sternal angle on the vagus nerve"
    ]
    if marker_name in bony_landmarks_names:
        return True
    return False


class MeshType_3d_nerve1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh based on data supplied by an input file.
    """

    @classmethod
    def getName(cls):
        return "3D Nerve 1"

    @classmethod
    def getParameterSetNames(cls):
        return [
            'Default',
            'Human Left Vagus 1',
            'Human Right Vagus 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            'Number of elements along the trunk pre-fit': 20,
            'Number of elements along the trunk': 40,
            'Trunk proportion': 0.6,
            'Trunk fit number of iterations': 5,
            'Default trunk diameter mm': 3.0,
            'Branch diameter trunk proportion': 0.5
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            'Number of elements along the trunk pre-fit',
            'Number of elements along the trunk',
            'Trunk proportion',
            'Trunk fit number of iterations',
            'Default trunk diameter mm',
            'Branch diameter trunk proportion'
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            'Number of elements along the trunk',
            'Number of elements along the trunk pre-fit'
        ]:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Trunk proportion',
            'Branch diameter trunk proportion'
        ]:
            if options[key] < 0.1:
                options[key] = 0.1
            elif options[key] > 1.0:
                options[key] = 1.0
        if options['Trunk fit number of iterations'] < 0:
            options['Trunk fit number of iterations'] = 0
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the 3d mesh for a vagus with branches, incorporating 1d central line,
        2d epineureum and 3d box based on constant vagus radius.

        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        return: list of AnnotationGroup, None
        """
        trunk_elements_count_prefit = options['Number of elements along the trunk pre-fit']
        trunk_elements_count = options['Number of elements along the trunk']
        trunk_proportion = options['Trunk proportion']
        trunk_fit_iterations = options['Trunk fit number of iterations']
        default_trunk_radius_mm = 0.5 * options['Default trunk diameter mm']
        branch_radius_proportion = options['Branch diameter trunk proportion']

        # Zinc setup for vagus scaffold
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()

        # node - geometric coordinates
        value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                        Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for value_label in value_labels[1:]:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, value_label, 1)

        # 3D box
        mesh3d = fieldmodule.findMeshByDimension(3)
        hermite_bilinear_basis = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        hermite_bilinear_basis.setFunctionType(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft3d = mesh3d.createElementfieldtemplate(hermite_bilinear_basis)
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

        # 3D box branch root
        # scale factors:
        # + 1 global scale factor (-1.0)
        # + 2 offsets of centroid in parent xi2, xi3
        # + 4 element scale factors for parent Hermite basis functions
        # + 4 element scale factors for parent Hermite basis function derivatives
        # + 9 coefficients of dB/dT, derivatives of branch derivatives w.r.t. trunk derivatives at root point
        eft3dBR = mesh3d.createElementfieldtemplate(hermite_bilinear_basis)
        setEftScaleFactorIds(eft3dBR, [1], [], 19)
        fns = [2, 3, 4, 5]
        dfns = [6, 7, 8, 9]
        # these are values in [-1.0, 1.0] which scale d2 and d3 to add to centroid x for non-centroid branch start
        cxd2 = 10
        cxd3 = 11
        # si_db_dt = [[10, 11, 12], [13, 14, 15], [16, 17, 18]]
        c11 = 12
        c12 = 13
        c13 = 14
        c21 = 15
        c22 = 16
        c23 = 17
        c31 = 18
        c32 = 19
        c33 = 20
        pd2 = [[None, None, fns[0], fns[1], None, None],
               [None, None, fns[2], fns[3], None, None]]
        pd3 = [[None, None, None, None, fns[0], fns[1]],
               [None, None, None, None, fns[2], fns[3]]]
        bx = [[[fns[0]], [fns[1]], None, None, None, None],
              [[fns[2]], [fns[3]], None, None, None, None]]
        bd1 = [[[c11, dfns[0]], [c11, dfns[1]], [c12, fns[0]], [c12, fns[1]], [c13, fns[0]], [c13, fns[1]]],
               [[c11, dfns[2]], [c11, dfns[3]], [c12, fns[2]], [c12, fns[3]], [c13, fns[2]], [c13, fns[3]]]]
        bd2 = [[[c21, dfns[0]], [c21, dfns[1]], [c22, fns[0]], [c22, fns[1]], [c23, fns[0]], [c23, fns[1]]],
               [[c21, dfns[2]], [c21, dfns[3]], [c22, fns[2]], [c22, fns[3]], [c23, fns[2]], [c23, fns[3]]]]
        bd3 = [[[c31, dfns[0]], [c31, dfns[1]], [c32, fns[0]], [c32, fns[1]], [c33, fns[0]], [c33, fns[1]]],
               [[c31, dfns[2]], [c31, dfns[3]], [c32, fns[2]], [c32, fns[3]], [c33, fns[2]], [c33, fns[3]]]]

        local_nodes = [1, 3]
        value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                       Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                       Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]
        ln = 1
        for n3 in range(2):
            s3 = [1] if (n3 == 0) else []
            for n2 in range(2):
                s2 = [1] if (n2 == 0) else []
                for fv in range(2):
                    old_value_label = value_labels[fv]
                    expression_terms = []
                    ubx, ubd2, ubd3 = (bx, bd2, bd3) if (fv == 0) else (bd1, None, None)
                    for tn in range(2):
                        tln = local_nodes[tn]
                        for pv in range(6):
                            value_label = value_labels[pv]
                            if ubx[tn][pv] is not None:
                                expression_terms.append((tln, value_label, ubx[tn][pv]))
                            if ubx is bx:
                                # offset from parent centroid
                                if pd2[tn][pv]:
                                    expression_terms.append((tln, value_label, [cxd2, pd2[tn][pv]]))
                                if pd3[tn][pv]:
                                    expression_terms.append((tln, value_label, [cxd3, pd3[tn][pv]]))
                            if ubd2:
                                expression_terms.append((tln, value_label, s2 + ubd2[tn][pv]))
                            if ubd3:
                                expression_terms.append((tln, value_label, s3 + ubd3[tn][pv]))
                    remapEftNodeValueLabelWithNodes(eft3dBR, ln, old_value_label, expression_terms)
                ln += 2
        ln = 2
        for n3 in range(2):
            s3 = [1] if (n3 == 0) else []
            for n2 in range(2):
                s2 = [1] if (n2 == 0) else []
                remapEftNodeValueLabel(eft3dBR, [ln], Node.VALUE_LABEL_VALUE,
                                       [(Node.VALUE_LABEL_VALUE, []),
                                        (Node.VALUE_LABEL_D_DS2, s2),
                                        (Node.VALUE_LABEL_D_DS3, s3)])
                remapEftNodeValueLabel(eft3dBR, [ln], Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, []),
                                        (Node.VALUE_LABEL_D2_DS1DS2, s2),
                                        (Node.VALUE_LABEL_D2_DS1DS3, s3)])
                ln += 2
        # remap nodes so that parent nodes are local nodes 1 and 2, last regular node is local node 3
        remapEftLocalNodes(eft3dBR, 3, [1, 3, 2, 3, 3, 3, 3, 3])
        elementtemplate_branch_root = mesh3d.createElementtemplate()
        elementtemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate_branch_root.defineField(coordinates, -1, eft3dBR)

        # 1D centroid
        mesh1d = fieldmodule.findMeshByDimension(1)
        hermite_basis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft1d = mesh1d.createElementfieldtemplate(hermite_basis)
        linetemplate = mesh1d.createElementtemplate()
        linetemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplate.defineField(coordinates, -1, eft1d)

        # 1D centroid branch root
        # scale factors:
        # + 1 global scale factor (-1.0)
        # + 2 offsets of centroid in parent xi2, xi3
        # + 4 element scale factors for parent Hermite basis functions
        # + 4 element scale factors for parent Hermite basis function derivatives
        # + 3 coefficients of dB/dT for derivatives of branch derivatives w.r.t. trunk derivatives at root, d1 only
        eft1dBR = mesh1d.createElementfieldtemplate(hermite_basis)
        eft1dBR.setNumberOfLocalNodes(3)
        setEftScaleFactorIds(eft1dBR, [], [], 13)
        ln = 1
        for fv in range(2):
            old_value_label = value_labels[fv]
            expression_terms = []
            ubx = bx if (fv == 0) else bd1
            for tn in range(2):
                tln = local_nodes[tn]
                for pv in range(6):
                    value_label = value_labels[pv]
                    if ubx[tn][pv] is not None:
                        expression_terms.append((tln, value_label, [c - 1 for c in ubx[tn][pv]]))
                    if ubx is bx:
                        # offset from parent centroid
                        if pd2[tn][pv]:
                            expression_terms.append((tln, value_label, [cxd2 - 1, pd2[tn][pv] - 1]))
                        if pd3[tn][pv]:
                            expression_terms.append((tln, value_label, [cxd3 - 1, pd3[tn][pv] - 1]))
            remapEftNodeValueLabelWithNodes(eft1dBR, ln, old_value_label, expression_terms)
        ln += 2
        # remap nodes so that parent nodes are local nodes 1 and 2, last regular node is local node 3
        remapEftLocalNodes(eft1dBR, 3, [1, 3, 2])
        linetemplate_branch_root = mesh1d.createElementtemplate()
        linetemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplate_branch_root.defineField(coordinates, -1, eft1dBR)

        # 2D epineurium - 4 elements around circle
        # scale factors:
        # + 3 global scale factor (-1.0, 0.5, pi/4)
        scalefactors2d = [-1.0, 0.5, 0.25 * math.pi]
        c_minus1 = 1
        c_half = 2
        c_pi__4 = 3
        mesh2d = fieldmodule.findMeshByDimension(2)
        bicubic_hermite_serendipity_basis = (
            fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
        facetemplate_and_eft_list = []
        for f in range(4):
            eft2d = mesh2d.createElementfieldtemplate(bicubic_hermite_serendipity_basis)
            setEftScaleFactorIds(eft2d, [1, 2, 3], [])
            ln = 1
            for n2 in range(2):
                for n1 in range(2):
                    value_expression_terms = [(Node.VALUE_LABEL_VALUE, [])]
                    d_ds1_expression_terms = [(Node.VALUE_LABEL_D_DS1, [])]
                    d_ds2_expression_terms = []
                    pole = (f + n2) % 4
                    if pole == 0:
                        value_expression_terms.append((Node.VALUE_LABEL_D_DS2, [c_half]))
                        d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS2, [c_half]))
                        d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS3, [c_pi__4]))
                    elif pole == 1:
                        value_expression_terms.append((Node.VALUE_LABEL_D_DS3, [c_half]))
                        d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS3, [c_half]))
                        d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS2, [c_minus1, c_pi__4]))
                    elif pole == 2:
                        value_expression_terms.append((Node.VALUE_LABEL_D_DS2, [c_minus1, c_half]))
                        d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS2, [c_minus1, c_half]))
                        d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS3, [c_minus1, c_pi__4]))
                    else:  # pole == 3:
                        value_expression_terms.append((Node.VALUE_LABEL_D_DS3, [c_minus1, c_half]))
                        d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS3, [c_minus1, c_half]))
                        d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS2, [c_pi__4]))
                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_VALUE, value_expression_terms)
                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_D_DS1, d_ds1_expression_terms)
                    remapEftNodeValueLabel(eft2d, [ln], Node.VALUE_LABEL_D_DS2, d_ds2_expression_terms)
                    ln += 1
            remapEftLocalNodes(eft2d, 2, [1, 2, 1, 2])
            facetemplate = mesh2d.createElementtemplate()
            facetemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplate.defineField(coordinates, -1, eft2d)
            facetemplate_and_eft_list.append((facetemplate, eft2d))

        # 2D epineureum branch root - 4 elements around circle
        # scale factors:
        # + 3 global scale factor (-1.0, 0.5, pi/4)
        # + 2 offsets of centroid in parent xi2, xi3
        # + 4 element scale factors for parent Hermite basis functions
        # + 4 element scale factors for parent Hermite basis function derivatives
        # + 9 coefficients of dB/dT, derivatives of branch derivatives w.r.t. trunk derivatives at root point
        # need to offset the following by 2 to fit 2 more global scale factors
        fns = [4, 5, 6, 7]
        dfns = [8, 9, 10, 11]
        # these are values in [-1.0, 1.0] which scale d2 and d3 to add to centroid x for non-centroid branch start
        cxd2 = 12
        cxd3 = 13
        # si_db_dt = [[10, 11, 12], [13, 14, 15], [16, 17, 18]]
        c11 = 14
        c12 = 15
        c13 = 16
        c21 = 17
        c22 = 18
        c23 = 19
        c31 = 20
        c32 = 21
        c33 = 22
        pd2 = [[None, None, fns[0], fns[1], None, None],
               [None, None, fns[2], fns[3], None, None]]
        pd3 = [[None, None, None, None, fns[0], fns[1]],
               [None, None, None, None, fns[2], fns[3]]]
        bx = [[[fns[0]], [fns[1]], None, None, None, None],
              [[fns[2]], [fns[3]], None, None, None, None]]
        bd1 = [[[c11, dfns[0]], [c11, dfns[1]], [c12, fns[0]], [c12, fns[1]], [c13, fns[0]], [c13, fns[1]]],
               [[c11, dfns[2]], [c11, dfns[3]], [c12, fns[2]], [c12, fns[3]], [c13, fns[2]], [c13, fns[3]]]]
        bd2 = [[[c21, dfns[0]], [c21, dfns[1]], [c22, fns[0]], [c22, fns[1]], [c23, fns[0]], [c23, fns[1]]],
               [[c21, dfns[2]], [c21, dfns[3]], [c22, fns[2]], [c22, fns[3]], [c23, fns[2]], [c23, fns[3]]]]
        bd3 = [[[c31, dfns[0]], [c31, dfns[1]], [c32, fns[0]], [c32, fns[1]], [c33, fns[0]], [c33, fns[1]]],
               [[c31, dfns[2]], [c31, dfns[3]], [c32, fns[2]], [c32, fns[3]], [c33, fns[2]], [c33, fns[3]]]]
        bicubic_value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2]
        facetemplate_and_eft_list_branch_root = []
        for f in range(4):
            eft2dBR = mesh2d.createElementfieldtemplate(bicubic_hermite_serendipity_basis)
            setEftScaleFactorIds(eft2dBR, [1, 2, 3], [], 19)
            ln = 1
            for n2 in range(2):
                ubx = bx  # common centroid for all branch root nodes
                ubd1 = bd1  # all d1 derivatives are parallel at branch root
                pole = (f + n2) % 4
                if pole == 0:
                    ubxo, ubxos = bd2, [c_half]
                    ubd2, ubd2s = bd3, [c_pi__4]
                elif pole == 1:
                    ubxo, ubxos = bd3, [c_half]
                    ubd2, ubd2s = bd2, [c_minus1, c_pi__4]
                elif pole == 2:
                    ubxo, ubxos = bd2, [c_minus1, c_half]
                    ubd2, ubd2s = bd3, [c_minus1, c_pi__4]
                else:  # pole == 3:
                    ubxo, ubxos = bd3, [c_minus1, c_half]
                    ubd2, ubd2s = bd2, [c_pi__4]
                for fv in range(3):
                    old_value_label = bicubic_value_labels[fv]
                    expression_terms = []
                    for tn in range(2):
                        tln = local_nodes[tn]
                        for pv in range(6):
                            value_label = value_labels[pv]
                            if old_value_label == Node.VALUE_LABEL_VALUE:
                                # start at parent centroid
                                if ubx[tn][pv] is not None:
                                    expression_terms.append((tln, value_label, ubx[tn][pv]))
                                # offset from parent centroid
                                if pd2[tn][pv]:
                                    expression_terms.append((tln, value_label, [cxd2, pd2[tn][pv]]))
                                if pd3[tn][pv]:
                                    expression_terms.append((tln, value_label, [cxd3, pd3[tn][pv]]))
                                # offset from branch centroid
                                expression_terms.append((tln, value_label, ubxos + ubxo[tn][pv]))
                            elif old_value_label == Node.VALUE_LABEL_D_DS1:
                                expression_terms.append((tln, value_label, ubd1[tn][pv]))
                            else:  # old_value_label == Node.VALUE_LABEL_D_DS2:
                                expression_terms.append((tln, value_label, ubd2s + ubd2[tn][pv]))
                    remapEftNodeValueLabelWithNodes(eft2dBR, ln, old_value_label, expression_terms)
                ln += 2
            ln = 2
            for n2 in range(2):
                value_expression_terms = [(Node.VALUE_LABEL_VALUE, [])]
                d_ds1_expression_terms = [(Node.VALUE_LABEL_D_DS1, [])]
                d_ds2_expression_terms = []
                pole = (f + n2) % 4
                if pole == 0:
                    value_expression_terms.append((Node.VALUE_LABEL_D_DS2, [2]))
                    d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS2, [2]))
                    d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS3, [3]))
                elif pole == 1:
                    value_expression_terms.append((Node.VALUE_LABEL_D_DS3, [2]))
                    d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS3, [2]))
                    d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS2, [1, 3]))
                elif pole == 2:
                    value_expression_terms.append((Node.VALUE_LABEL_D_DS2, [1, 2]))
                    d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS2, [1, 2]))
                    d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS3, [1, 3]))
                else:  # pole == 3:
                    value_expression_terms.append((Node.VALUE_LABEL_D_DS3, [1, 2]))
                    d_ds1_expression_terms.append((Node.VALUE_LABEL_D2_DS1DS3, [1, 2]))
                    d_ds2_expression_terms.append((Node.VALUE_LABEL_D_DS2, [3]))
                remapEftNodeValueLabel(eft2dBR, [ln], Node.VALUE_LABEL_VALUE, value_expression_terms)
                remapEftNodeValueLabel(eft2dBR, [ln], Node.VALUE_LABEL_D_DS1, d_ds1_expression_terms)
                remapEftNodeValueLabel(eft2dBR, [ln], Node.VALUE_LABEL_D_DS2, d_ds2_expression_terms)
                ln += 2

            # remap nodes so that parent nodes are local nodes 1 and 2, last regular node is local node 3
            remapEftLocalNodes(eft2dBR, 3, [1, 3, 2, 3])
            facetemplate_branch_root = mesh2d.createElementtemplate()
            facetemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplate_branch_root.defineField(coordinates, -1, eft2dBR)
            facetemplate_and_eft_list_branch_root.append((facetemplate_branch_root, eft2dBR))

        # =========
        # Load Data
        # =========

        vagus_data = load_vagus_data(region)
        invalid_data = not vagus_data
        if not invalid_data:
            marker_data = vagus_data.get_level_markers()
            trunk_data = vagus_data.get_trunk_coordinates()
            if len(marker_data) < 2:
                logger.error("Missing or incomplete data. At least two landmarks are expected in the data.")
                invalid_data = True
            if len(trunk_data) < 2:
                logger.error("Missing or incomplete data. At least two trunk points are required.")
                invalid_data = True
            if vagus_data.get_side_label() not in ['left', 'right']:
                logger.error("Missing left or right marker side indication.")
                invalid_data = True
        if invalid_data:
            return [], None

        # ===========
        # Build Trunk
        # ===========

        # fit trunk with radius and orientation
        only_1d_trunk = False
        region1d = region if only_1d_trunk else region.createRegion()
        tx, td1, td2, td12, td3, td13, default_trunk_radius = generate_trunk_1d(
            vagus_data, trunk_proportion, trunk_elements_count_prefit, trunk_elements_count,
            trunk_fit_iterations, default_trunk_radius_mm, region1d)
        trunk_length = getCubicHermiteCurvesLength(tx, td1)
        trunk_mean_element_length = trunk_length / trunk_elements_count

        default_branch_radius = branch_radius_proportion * default_trunk_radius
        default_branch_size = 2.0 * default_branch_radius  # half span of box is double radius

        annotation_groups = []
        annotation_term_map = vagus_data.get_annotation_term_map()
        trunk_group_name = vagus_data.get_trunk_group_name()
        trunk_group = AnnotationGroup(region, (trunk_group_name, annotation_term_map[trunk_group_name]))
        annotation_groups.append(trunk_group)

        if only_1d_trunk:
            return annotation_groups, None

        # vagus annotation groups
        centroid_annotation_group = AnnotationGroup(region, get_vagus_marker_term("vagus centroid"))
        annotation_groups.append(centroid_annotation_group)
        centroid_mesh_group = centroid_annotation_group.getMeshGroup(mesh1d)

        epineurium_annotation_group = AnnotationGroup(region, get_vagus_marker_term("vagus epineureum"))
        annotation_groups.append(epineurium_annotation_group)
        epineurium_mesh_group = epineurium_annotation_group.getMeshGroup(mesh2d)

        # trunk annotation groups
        trunk_mesh_group = trunk_group.getMeshGroup(mesh3d)
        trunk_face_mesh_group = trunk_group.getMeshGroup(mesh2d)
        trunk_line_mesh_group = trunk_group.getMeshGroup(mesh1d)

        node_identifier = 1
        element_identifier = 1
        face_identifier = 1
        line_identifier = 1

        # create trunk nodes
        trunk_nodes_count = trunk_elements_count + 1
        tnid = []
        for n in range(trunk_elements_count + 1):
            node = nodes.createNode(node_identifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, td1[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, td2[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, td12[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, td3[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, td13[n])
            tnid.append(node_identifier)
            node_identifier += 1

        # create trunk elements
        for e in range(trunk_elements_count):
            nids = [e + 1, e + 2]

            element = mesh3d.createElement(element_identifier, elementtemplate)
            element.setNodesByIdentifier(eft3d, nids)
            element.setScaleFactors(eft3d, [-1.0])
            trunk_mesh_group.addElement(element)
            element_identifier += 1

            line = mesh1d.createElement(line_identifier, linetemplate)
            line.setNodesByIdentifier(eft1d, nids)
            centroid_mesh_group.addElement(line)
            trunk_line_mesh_group.addElement(line)
            line_identifier += 1

            for f in range(4):
                facetemplate, eft2d = facetemplate_and_eft_list[f]
                face = mesh2d.createElement(face_identifier, facetemplate)
                face.setNodesByIdentifier(eft2d, nids)
                face.setScaleFactors(eft2d, scalefactors2d)
                epineurium_mesh_group.addElement(face)
                trunk_face_mesh_group.addElement(face)
                face_identifier += 1

        # create trunk markers
        is_left = vagus_data.get_side_label() == 'left'
        vagus_level_terms = get_left_vagus_marker_locations_list() if is_left \
            else get_right_vagus_marker_locations_list()
        ordered_marker_data = []  # list from top to bottom of nerve of (name, material_coordinate)
        for marker_term_name, material_coordinate in vagus_level_terms.items():
            for idx, data in enumerate(ordered_marker_data):
                if material_coordinate < data[1]:
                    break
            else:
                idx = len(ordered_marker_data)
            ordered_marker_data.insert(idx, (marker_term_name, material_coordinate))
        for marker_name, material_coordinate in ordered_marker_data:
            if material_coordinate > trunk_proportion:
                continue
            annotationGroup = findOrCreateAnnotationGroupForTerm(
                annotation_groups, region, get_vagus_marker_term(marker_name), isMarker=True)
            element_index_real = (material_coordinate / trunk_proportion) * trunk_elements_count
            if element_index_real >= trunk_elements_count:
                element_index, xi = trunk_elements_count - 1, 1.0
            else:
                element_index = math.floor(element_index_real)
                xi = element_index_real - element_index
            annotationGroup.createMarkerNode(
                node_identifier, element=mesh3d.findElementByIdentifier(element_index + 1), xi=[xi, 0.5, 0.5])
            node_identifier += 1

        # ==============
        # Build Branches
        # ==============

        parent_parameters = {}
        parent_parameters[trunk_group_name] = (tx, td1, td2, td12, td3, td13, tnid)

        # following is assigned to branch coordinates before finding trunk location
        branch_start_coordinates = fieldmodule.createFieldConstant([0.0, 0.0, 0.0])
        # note zinc will not use the fast find cache due to field modifications with change caching on
        find_trunk_location = fieldmodule.createFieldFindMeshLocation(
            branch_start_coordinates, coordinates, trunk_mesh_group)
        find_trunk_location.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_NEAREST)

        # iterate over branches off trunk, and branches of branches
        visited_branches_order = []
        branch_root_parameters = {}
        branch_data = vagus_data.get_branch_data()
        branch_parent_map = vagus_data.get_branch_parent_map()
        queue = [branch for branch in branch_parent_map.keys() if branch_parent_map[branch] == trunk_group_name]
        while queue:
            branch_name = queue.pop(0)
            if branch_name in visited_branches_order:
                logger.warning("already processed branch " + branch_name)
                continue
            visited_branches_order.append(branch_name)

            branch_px = [branch_node[0] for branch_node in branch_data[branch_name]]
            branch_parent_name = branch_parent_map[branch_name]
            trunk_is_parent = branch_parent_name == trunk_group_name
            # print(branch_name, '<--', branch_parent_name)

            tx, td1, td2, td12, td3, td13, tnid = parent_parameters[branch_parent_name]

            # get point in trunk volume closest to first point in branch data
            # parent_group = trunk_group
            parent_mesh_group = trunk_mesh_group
            find_parent_location = find_trunk_location
            if not trunk_is_parent:
                parent_group = fieldmodule.findFieldByName(branch_parent_name).castGroup()
                parent_mesh_group = parent_group.getMeshGroup(mesh3d)
                find_parent_location = fieldmodule.createFieldFindMeshLocation(
                    branch_start_coordinates, coordinates, parent_mesh_group)
                find_parent_location.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_NEAREST)
                del parent_group
            branch_start_x = branch_px[1 if (len(branch_px) > 2) else 0]
            fieldcache.clearLocation()
            branch_start_coordinates.assignReal(fieldcache, branch_start_x)
            parent_element, parent_xi = find_parent_location.evaluateMeshLocation(fieldcache, 3)
            if not parent_element.isValid():
                logger.error("Nerve: branch " + branch_name + " start point could not be found in parent nerve")
                continue
            cxd2 = 2.0 * (parent_xi[1] - 0.5)
            cxd3 = 2.0 * (parent_xi[2] - 0.5)
            parent_first_element = parent_mesh_group.createElementiterator().next()
            parent_location = (parent_element.getIdentifier() - parent_first_element.getIdentifier(), parent_xi[0])
            if (not trunk_is_parent) and (parent_location[0] == 0):
                # can't have branch from the root element of a branch
                parent_location = (1, 0.0)
            branch_max_element_length = trunk_mean_element_length
            # unless only 2 points, don't want first point which only marks connection on trunk centroid
            px = branch_px[1:] if (len(branch_px) > 2) else branch_px
            ax, ad1 = get_curve_from_points(px, maximum_element_length=branch_max_element_length)
            bx, bd1 = fit_hermite_curve(ax, ad1, px)
            branch_length = getCubicHermiteCurvesLength(bx, bd1)
            branch_elements_count = max(2, math.ceil(branch_length / branch_max_element_length))
            cx, cd1 = sampleCubicHermiteCurves(bx, bd1, branch_elements_count, arcLengthDerivatives=True)[0:2]

            # parent interpolation
            pn1 = parent_location[0]
            pn2 = pn1 + 1
            pxi = parent_location[1]
            fns = list(getCubicHermiteBasis(pxi))  # for x, d2, d3
            dfns = list(getCubicHermiteBasisDerivatives(pxi))  # for d1
            # first derivatives interpolated on parent:
            pd1 = [dot(dfns, [tx[pn1][c], td1[pn1][c], tx[pn2][c], td1[pn2][c]]) for c in range(3)]
            pd2 = [dot(fns, [td2[pn1][c], td12[pn1][c], td2[pn2][c], td12[pn2][c]]) for c in range(3)]
            pd3 = [dot(fns, [td3[pn1][c], td13[pn1][c], td3[pn2][c], td13[pn2][c]]) for c in range(3)]
            # first derivatives required on branch:
            bd1 = cd1[0]
            branch_root_size = branch_radius_proportion * magnitude(pd2) if trunk_is_parent else default_branch_size
            bd3 = set_magnitude(cross(pd1, bd1), branch_root_size)
            bd2 = set_magnitude(cross(bd3, bd1), branch_root_size)

            basis_from = [pd1, pd2, pd3]
            basis_to = [bd1, bd2, bd3]
            coefs = matrix_mult(basis_to, matrix_inv(basis_from))

            # branch annotation groups
            branch_box_group = AnnotationGroup(region, (branch_name, annotation_term_map[branch_name]))
            annotation_groups.append(branch_box_group)
            branch_box_mesh_group = branch_box_group.getMeshGroup(mesh3d)
            branch_box_face_mesh_group = branch_box_group.getMeshGroup(mesh2d)
            branch_box_line_mesh_group = branch_box_group.getMeshGroup(mesh1d)

            # get side derivatives, minimising rotation from trunk
            # dir2 = normalize(bd2)
            dir3 = normalize(bd3)
            cd2 = [bd2]
            cd3 = [bd3]
            for e in range(len(cx) - 1):
                dir1, dir2, dir3 = track_curve_side_direction(cx, cd1, dir3, (e, 0.0), (e, 1.0))
                cd2.append(set_magnitude(dir2, default_branch_size))
                cd3.append(set_magnitude(dir3, default_branch_size))
            cd12, cd13 = smoothCurveSideCrossDerivatives(cx, cd1, [cd2, cd3])

            # create branch elements and nodes past root
            cnid = [None]  # no first node on branch
            for e in range(len(cx) - 1):
                n = e + 1
                node = nodes.createNode(node_identifier, nodetemplate)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, cd12[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, cd13[n])

                if e == 0:
                    # branch root 3D element
                    nids = [tnid[pn1], tnid[pn2], node_identifier]
                    scalefactors = [-1] + fns + dfns + [cxd2, cxd3] + coefs[0] + coefs[1] + coefs[2]
                    element = mesh3d.createElement(element_identifier, elementtemplate_branch_root)
                    element.setNodesByIdentifier(eft3dBR, nids)
                    element.setScaleFactors(eft3dBR, scalefactors)
                    # branch root 1D line
                    scalefactors = fns + dfns + [cxd2, cxd3] + coefs[0]
                    line = mesh1d.createElement(line_identifier, linetemplate_branch_root)
                    line.setNodesByIdentifier(eft1dBR, nids)
                    line.setScaleFactors(eft1dBR, scalefactors)
                else:
                    # branch regular 3D element
                    nids = [node_identifier - 1, node_identifier]
                    element = mesh3d.createElement(element_identifier, elementtemplate)
                    element.setNodesByIdentifier(eft3d, nids)
                    element.setScaleFactors(eft3d, [-1.0])
                    # branch regular 1D line
                    line = mesh1d.createElement(line_identifier, linetemplate)
                    line.setNodesByIdentifier(eft1d, nids)
                    line.setScaleFactors(eft1d, [-1.0])
                branch_box_mesh_group.addElement(element)
                centroid_mesh_group.addElement(line)
                branch_box_line_mesh_group.addElement(line)
                cnid.append(node_identifier)
                element_identifier += 1
                line_identifier += 1

                # 2D epineurium
                for f in range(4):
                    if e == 0:
                        # branch root 2D face
                        facetemplate_branch_root, eft2dBR = facetemplate_and_eft_list_branch_root[f]
                        nids = [tnid[pn1], tnid[pn2], node_identifier]
                        scalefactors = scalefactors2d + fns + dfns + [cxd2, cxd3] + coefs[0] + coefs[1] + coefs[2]
                        face = mesh3d.createElement(face_identifier, facetemplate_branch_root)
                        face.setNodesByIdentifier(eft2dBR, nids)
                        face.setScaleFactors(eft2dBR, scalefactors)
                    else:
                        # branch regular 2D face
                        facetemplate, eft2d = facetemplate_and_eft_list[f]
                        nids = [node_identifier - 1, node_identifier]
                        face = mesh2d.createElement(face_identifier, facetemplate)
                        face.setNodesByIdentifier(eft2d, nids)
                        face.setScaleFactors(eft2d, scalefactors2d)
                    epineurium_mesh_group.addElement(face)
                    branch_box_face_mesh_group.addElement(face)
                    face_identifier += 1

                node_identifier += 1

            # add branches of branches, storing parameters for embedding sub-branch root
            child_branches = [branch for branch in branch_parent_map.keys() if branch_parent_map[branch] == branch_name]
            if child_branches:
                queue = child_branches + queue
                parent_parameters[branch_name] = (cx, cd1, cd2, cd12, cd3, cd13, cnid)

        return annotation_groups, None

        # add material coordinates
        # print('Adding material coordinates...')
        # vagus trunk goes along z axis with origin as top of the vagus
        coordinates.setName("vagus coordinates")  # temporarily rename
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        region.write(sir)
        result, buffer = srm.getBuffer()
        coordinates.setName("coordinates")  # restore name before reading vagus coordinates back in
        sir = region.createStreaminformationRegion()
        sir.createStreamresourceMemoryBuffer(buffer)
        # read and merge with region, thus having coordinates and vagus coordinates in the region together
        region.read(sir)
        vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()

        # calculate derivatives
        derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
        derivative_xi2 = mesh3d.getChartDifferentialoperator(1, 2)
        derivative_xi3 = mesh3d.getChartDifferentialoperator(1, 3)

        zero = [0.0, 0.0, 0.0]
        rescaled_step = rescaled_vagus_trunk_length / (trunk_nodes_count - 1)

        elem_iter = mesh3d.createElementiterator()
        element = elem_iter.next()
        x = None
        while element.isValid():
            # trunk elements first, followed by branch elements (with first element having more than 2 local nodes)
            element_id = element.getIdentifier()
            eft = element.getElementfieldtemplate(vagus_coordinates, -1)
            local_nodes_count = eft.getNumberOfLocalNodes()
            if local_nodes_count == 2:
                if element_id == 1:
                    # first trunk element
                    x = [0.0, 0.0, 0.0]
                    d1 = [0, 0, rescaled_step]
                    [d2], [d3] = set_group_nodes_derivatives_orthogonal([d1], vagus_radius)

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

                d1 = set_magnitude(d1, rescaled_step)
                d2 = set_magnitude(d2, vagus_branch_radius)
                d3 = set_magnitude(d3, vagus_branch_radius)
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

        # calculate & add markers vagus coordinates
        # marker_names = fieldmodule.findFieldByName("marker_name")
        marker_location = fieldmodule.findFieldByName("marker_location")
        marker_coordinates = fieldmodule.findFieldByName("marker coordinates")
        host_vagus_coordinates = fieldmodule.createFieldEmbedded(vagus_coordinates, marker_location)

        marker_nodetemplate = nodes.createNodetemplate()
        marker_nodetemplate.defineField(vagus_coordinates)
        marker_nodetemplate.undefineField(marker_coordinates)

        marker_group = fieldmodule.findFieldByName("marker").castGroup()
        marker_nodeset_group = marker_group.getNodesetGroup(nodes)
        nodeiterator = marker_nodeset_group.createNodeiterator()
        node = nodeiterator.next()
        while node.isValid():
            node.merge(marker_nodetemplate)
            node = nodeiterator.next()
        fieldassignment = vagus_coordinates.createFieldassignment(host_vagus_coordinates)
        fieldassignment.setNodeset(marker_nodeset_group)
        fieldassignment.assign()

        marker_coordinates.setManaged(False)
        del marker_coordinates

        # print('Adding extra visualisation groups...')
        branch_common_groups = vagus_data.get_branch_common_group_map()
        for branch_common_name, branch_names in branch_common_groups.items():
            branch_common_group = AnnotationGroup(region, (branch_common_name, ""))
            branch_common_mesh_group = branch_common_group.getMeshGroup(mesh3d)
            annotation_groups.append(branch_common_group)

            for branch_name in branch_names:
                branch_group = findAnnotationGroupByName(annotation_groups, branch_name)
                branch_mesh_group = branch_group.getMeshGroup(mesh3d)

                el_iter = branch_mesh_group.createElementiterator()
                element = el_iter.next()
                while element.isValid():
                    branch_common_mesh_group.addElement(element)
                    element = el_iter.next()

        # add cervical and thoracic trunk groups
        side_label = vagus_data.get_side_label()
        if side_label == 'left':
            cervical_trunk_group = AnnotationGroup(region, get_vagus_branch_term('left cervical vagus nerve'))
            thoracic_trunk_group = AnnotationGroup(region, get_vagus_branch_term('left thoracic vagus nerve'))
            boundary_marker = findAnnotationGroupByName(
                annotation_groups, "left level of superior border of the clavicle on the vagus nerve")

        if side_label == 'right':
            cervical_trunk_group = AnnotationGroup(region, get_vagus_branch_term('right cervical vagus nerve'))
            thoracic_trunk_group = AnnotationGroup(region, get_vagus_branch_term('right thoracic vagus nerve'))
            boundary_marker = findAnnotationGroupByName(
                annotation_groups, "right level of superior border of the clavicle on the vagus nerve")

        element, _ = boundary_marker.getMarkerLocation()
        boundary_element_id = element.getIdentifier()

        cervical_trunk_mesh_group = cervical_trunk_group.getMeshGroup(mesh3d)
        thoracic_trunk_mesh_group = thoracic_trunk_group.getMeshGroup(mesh3d)
        annotation_groups.append(cervical_trunk_group)
        annotation_groups.append(thoracic_trunk_group)

        trunk_group = findAnnotationGroupByName(annotation_groups, trunk_group_name)
        trunk_mesh_group = trunk_group.getMeshGroup(mesh3d)

        el_iter = trunk_mesh_group.createElementiterator()
        element = el_iter.next()
        while element.isValid():
            element_id = element.getIdentifier()
            if element_id < boundary_element_id:
                cervical_trunk_mesh_group.addElement(element)
            else:
                thoracic_trunk_mesh_group.addElement(element)
            element = el_iter.next()

        # print('Done\n')

        return annotation_groups, None

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

        epineurium_annotation_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                            get_vagus_marker_term("vagus epineureum"))
        epineurium_mesh_group = epineurium_annotation_group.getMeshGroup(mesh2d)
        vagusAnteriorLineAnnotationGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_vagus_marker_term("vagus anterior line"))
        vagusAnteriorLineMeshGroup = vagusAnteriorLineAnnotationGroup.getMeshGroup(mesh1d)

        faceIterator = epineurium_mesh_group.createElementiterator()
        quadrant = 0
        face = faceIterator.next()
        while face.isValid():
            if quadrant == 0:
                line = face.getFaceElement(4)
                vagusAnteriorLineMeshGroup.addElement(line)
            quadrant = (quadrant + 1) % 4
            face = faceIterator.next()


def generate_trunk_1d(vagus_data, trunk_proportion, trunk_elements_count_prefit, trunk_elements_count,
                      trunk_fit_iterations, default_trunk_radius_mm, region):
    """
    Build and fit a 1-D trunk curve to trunk data, calibrated to marker point positions.
    :param vagus_data: Vagus data extracted from input data region.
    :param trunk_proportion: Proportion up to 1.0 of whole vagus to build.
    :param trunk_elements_count_prefit: Number of elements in pre-fit mesh to trunk data.
    :param trunk_elements_count: Number of elements in final 1-D mesh.
    :param trunk_fit_iterations: Number of iterations in main trunk fit >= 1.
    :param default_trunk_radius_mm: Radius to use if no radius parameters, in mm. This is scaled in magnitude by
    factors of 10 until it is within the expected aspect ratio.
    :param region: Region to put the fitted 1-D geometry including marker points in.
    :return: tx, td1, td2, td12, td3, td13 (parameters for 1-D fitted trunk geometry, left and anterior side
    directions and rates of change w.r.t. d1), default_trunk_radius (in same units as data)
    """

    # 1. pre-fit to range of trunk data

    trunk_data_coordinates = vagus_data.get_trunk_coordinates()
    is_left = vagus_data.get_side_label() == 'left'
    raw_marker_data = vagus_data.get_level_markers()
    px = [e[0] for e in trunk_data_coordinates]
    bx, bd1 = get_curve_from_points(px, number_of_elements=trunk_elements_count_prefit)
    length = getCubicHermiteCurvesLength(bx, bd1)
    # outlier_length = 0.025 * length
    # # needs to be bigger if fewer elements:
    # if trunk_elements_count_prefit < 20:
    #     outlier_length += 0.075 * (20 - trunk_elements_count_prefit) / 19.0
    cx, cd1 = fit_hermite_curve(bx, bd1, px)  # , outlier_length=outlier_length)
    # resample to even size
    dx, dd1 = sampleCubicHermiteCurvesSmooth(cx, cd1, trunk_elements_count_prefit)[0:2]
    # generate_curve_mesh(region, dx, dd1, group_name=vagus_data.get_trunk_group_name())[0]
    # return dx, dd1, [], []

    # 2. make initial full trunk extending pre-fit based on furthest marker material coordinates

    vagus_level_terms = get_left_vagus_marker_locations_list() if is_left \
        else get_right_vagus_marker_locations_list()
    marker_data = []  # list from top to bottom of nerve of (name, material_coordinate, data_coordinates)
    for marker_term_name, material_coordinate in vagus_level_terms.items():
        if marker_term_name in raw_marker_data.keys():
            data_coordinates =raw_marker_data[marker_term_name]
            for idx, data in enumerate(marker_data):
                if material_coordinate < data[1]:
                    break
            else:
                idx = len(marker_data)
            marker_data.insert(idx, (marker_term_name, material_coordinate, data_coordinates))

    start_marker_material_coordinate = marker_data[0][1]
    start_marker_curve_location = getNearestLocationOnCurve(dx, dd1, marker_data[0][2])[0]
    start_marker_extra_length = 0.0
    if (start_marker_curve_location[0] == 0) and (start_marker_curve_location[1] < 1.0E-4):
        # case when start marker point is before start of curve
        start_marker_extra_length = max(0.0, dot(sub(dx[0], marker_data[0][2]), normalize(dd1[0])))
    end_marker_material_coordinate = marker_data[-1][1]
    end_marker_curve_location = getNearestLocationOnCurve(dx, dd1, marker_data[-1][2])[0]
    end_marker_extra_length = 0.0
    if (end_marker_curve_location[0] == (trunk_elements_count_prefit - 1)) and (
            end_marker_curve_location[1] > 0.9999):
        # case when end marker point is after end of curve
        end_marker_extra_length = max(0.0, dot(sub(marker_data[-1][2], dx[-1]), normalize(dd1[-1])))
    start_length, length, end_length = getCubicHermiteTrimmedCurvesLengths(
        dx, dd1, start_marker_curve_location, end_marker_curve_location)[0:3]

    # determine material coordinates of start and end of curve, from range of marker material coordinates
    marker_delta_material_length = end_marker_material_coordinate - start_marker_material_coordinate
    marker_delta_length = length + start_marker_extra_length + end_marker_extra_length
    material_length_per_length = marker_delta_material_length / marker_delta_length
    start_curve_material_coordinate = (start_marker_material_coordinate -
                                       (start_length - start_marker_extra_length) * material_length_per_length)
    end_curve_material_coordinate = (end_marker_material_coordinate +
                                     (end_length - end_marker_extra_length) * material_length_per_length)
    # extend curves at each end, by moving end node if short extension, or adding node if large
    start_extra_length = start_curve_material_coordinate / material_length_per_length
    start_direction = normalize(dd1[0])
    mag_d1 = magnitude(dd1[0])
    start_dx = sub(dx[0], mult(start_direction, start_extra_length))
    if start_extra_length < (0.51 * mag_d1):
        dx[0] = start_dx
        dd1[0] = set_magnitude(dd1[0], mag_d1 + 2.0 * start_extra_length)
    else:
        dx.insert(0, start_dx)
        dd1.insert(0, mult(start_direction, 2.0 * start_extra_length - mag_d1))
    end_extra_length = (trunk_proportion - end_curve_material_coordinate) / material_length_per_length
    end_direction = normalize(dd1[-1])
    mag_d1 = magnitude(dd1[-1])
    end_dx = add(dx[-1], mult(end_direction, end_extra_length))
    if start_extra_length < (0.51 * mag_d1):
        dx[-1] = end_dx
        dd1[-1] = set_magnitude(dd1[-1], mag_d1 + 2.0 * end_extra_length)
    else:
        dx.append(end_dx)
        dd1.append(mult(end_direction, 2.0 * end_extra_length - mag_d1))

    # resample to even size and final elements count - this is the initial state for full trunk fit with markers
    ex, ed1 = sampleCubicHermiteCurvesSmooth(dx, dd1, trunk_elements_count)[0:2]

    # 3. full trunk centroid fit

    fit_region = region.createRegion()
    fieldmodule = fit_region.getFieldmodule()
    coordinates = find_or_create_field_coordinates(fieldmodule)
    components_count = coordinates.getNumberOfComponents()
    mesh1d = fieldmodule.findMeshByDimension(1)
    trunk_group_name = vagus_data.get_trunk_group_name()
    with ChangeManager(fieldmodule):
        node_identifier = generate_curve_mesh(fit_region, ex, ed1, group_name=trunk_group_name)[0]
        trunk_group = find_or_create_field_group(fieldmodule, trunk_group_name)
        pr = vagus_data.get_trunk_radius()
        field_names_and_values = [("radius", pr)] if pr else []
        data_identifier = 1
        data_identifier = generate_datapoints(
            fit_region, px, data_identifier, field_names_and_values=field_names_and_values, group_name=trunk_group_name)
        radius = fieldmodule.findFieldByName("radius").castFiniteElement()

        # add marker points in order along trunk
        ordered_marker_data = []  # list of (name, curve_location)
        ordered_material_coordinate = []  # list of material coordinate
        for marker_term_name, material_coordinate in vagus_level_terms.items():
            for idx, prior_material_coordinate in enumerate(ordered_material_coordinate):
                if material_coordinate < prior_material_coordinate:
                    break
            else:
                idx = len(ordered_marker_data)
            real_element_xi = (material_coordinate / trunk_proportion) * trunk_elements_count
            element_index = min(int(real_element_xi), trunk_elements_count - 1)
            element_xi = (element_index + 1, real_element_xi - element_index)
            ordered_marker_data.insert(idx, (marker_term_name, element_xi))
            ordered_material_coordinate.insert(idx, material_coordinate)
        node_identifier = generate_mesh_marker_points(mesh1d, ordered_marker_data, node_identifier)

        # add data marker points:
        data_identifier = generate_datapoints(
            fit_region, list(raw_marker_data.values()), data_identifier,
            field_names_and_values=[("marker_name", list(raw_marker_data.keys()))], group_name="marker")

        zero_fibres = find_or_create_field_zero_fibres(fieldmodule)

    # note that fitting is very slow if done within ChangeManager as find mesh location is slow
    # this includes working with the user-supplied region which is called with ChangeManager on.
    fitter = GeometryFitter(region=fit_region)
    length = getCubicHermiteCurvesLength(ex, ed1)
    outlier_length = 0.025 * length
    fitter.getInitialFitterStepConfig().setGroupOutlierLength(None, outlierLength=outlier_length)
    # fitter.setDiagnosticLevel(1)
    fitter.setModelCoordinatesField(coordinates)
    fitter.setFibreField(zero_fibres)
    del zero_fibres
    fitter.defineCommonMeshFields()
    fitter.setDataCoordinatesField(coordinates)
    fitter.defineDataProjectionFields()
    fitter.setMarkerGroupByName("marker")
    fitter.initializeFit()

    # calibration_points_count = 23954
    points_count_calibration_factor = len(px) / 25000
    # calibration_length = 27840.0
    length_calibration_factor = length / 25000.0
    strain_penalty = 1000.0 * points_count_calibration_factor * length_calibration_factor
    curvature_penalty = 1.0E+8 * points_count_calibration_factor * (length_calibration_factor ** 3)
    marker_weight = 10.0 * points_count_calibration_factor
    sliding_factor = 0.0001

    if trunk_fit_iterations > 0:
        fit1 = FitterStepFit()
        fitter.addFitterStep(fit1)
        fit1.setGroupDataWeight("marker", marker_weight)
        fit1.setGroupStrainPenalty(None, [strain_penalty])
        fit1.setGroupCurvaturePenalty(None, [curvature_penalty])
        fit1.setGroupDataSlidingFactor(None, sliding_factor)
        fit1.run()
        del fit1

    if trunk_fit_iterations > 1:
        fit2 = FitterStepFit()
        fitter.addFitterStep(fit2)
        fit2.setGroupStrainPenalty(None, [0.1 * strain_penalty])
        fit2.setGroupCurvaturePenalty(None, [0.1 * curvature_penalty])
        fit2.setGroupDataSlidingFactor(None, 0.1 * sliding_factor)
        fit2.setNumberOfIterations(trunk_fit_iterations - 1)
        fit2.run()
        del fit2

    rms_error, max_error = fitter.getDataRMSAndMaximumProjectionError()

    fitter.cleanup()
    del fitter

    # fit radius
    if pr:
        gradient1_penalty = 1000.0 * points_count_calibration_factor * length_calibration_factor
        gradient2_penalty = 1.0E+8 * points_count_calibration_factor * (length_calibration_factor ** 3)
        define_and_fit_field(fit_region, "coordinates", "coordinates", "radius",
                             gradient1_penalty, gradient2_penalty, group_name=trunk_group_name)

    # extract fitted trunk parameters from trunk nodes (not marker points)
    length = getCubicHermiteCurvesLength(ex, ed1)
    default_trunk_radius = default_trunk_radius_mm
    scale_step = 10.0  # because data is frequently in 1/100 mm.
    max_radius = 0.02 * length
    while (default_trunk_radius * scale_step) < max_radius:
        default_trunk_radius *= scale_step
    tx = []
    td1 = []
    rx = []
    rd1 = []
    fieldcache = fieldmodule.createFieldcache()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodeiterator = nodes.createNodeiterator()
    node = nodeiterator.next()
    for n in range(trunk_elements_count + 1):
        fieldcache.setNode(node)
        result, x = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, components_count)
        result, d1 = coordinates.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, components_count)
        tx.append(x)
        td1.append(d1)
        if pr:
            result, x = radius.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, d1 = radius.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 1)
        else:
            x = default_trunk_radius
            d1 = 0.0
        rx.append(x)
        rd1.append(d1)
        node = nodeiterator.next()

    # remove all datapoints from previous fits.
    datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    datapoints.destroyAllNodes()

    # fit orientation
    debug_print = False
    orientation_dct = vagus_data.get_orientation_data()
    orientation_x = []
    orientation_twist_angles = []
    if orientation_dct:
        # convert to list in order down the trunk
        length = getCubicHermiteCurvesLength(tx, td1)
        orientation_directions = []  # list of anterior directions
        orientation_names = []  # list of orientation names
        orientation_locations = []  # list of curve location (element index, xi) for orientation points
        one_sqrt2 = 1.0 / math.sqrt(2.0)
        # weights of d2, d3 to give anterior direction for a given orientation direction name
        orientation_anterior_weights = {
            "left": (-1.0, 0.0),
            "left anterior": (-one_sqrt2, one_sqrt2),
            "anterior": (0.0, 1.0),
            "right anterior": (one_sqrt2, one_sqrt2),
            "right": (1.0, 0.0),
            "right posterior": (one_sqrt2, -one_sqrt2),
            "posterior": (0.0, -1.0),
            "left posterior": (-one_sqrt2, -one_sqrt2)
        }
        for name, x_list in orientation_dct.items():
            direction_name = name.split("orientation ", 1)[1]
            weights = orientation_anterior_weights.get(direction_name)
            if not weights:
                logger.warning("Nerve: Ignoring unrecognized orientation points with name '" + name + "'")
                continue
            wt2, wt3 = weights
            for data_x in x_list:
                curve_location, x = getNearestLocationOnCurve(tx, td1, data_x)
                e1 = curve_location[0]
                e2 = curve_location[0] + 1
                xi = curve_location[1]
                d1 = interpolateCubicHermiteDerivative(tx[e1], td1[e1], tx[e2], td1[e2], xi)
                dir1 = normalize(d1)
                projection = sub(data_x, x)
                normal_projection = rejection(projection, dir1)
                projection_error = magnitude(normal_projection)
                if projection_error < rms_error:
                    # skip orientation points within rmsError of trunk centroid as inaccurate
                    logger.warning("Nerve: Ignoring orientation point '" + name + "' at location", curve_location,
                                   "as projection error", projection_error, "is less than trunk RMS error", rms_error)
                    continue
                dirp = normalize(projection)
                if math.fabs(dot(dir1, dirp)) > 0.5:
                    # skip orientation points with severely non-normal projection e.g. past end of curve
                    logger.warning("Nerve: Ignoring orientation point '" + name + "' at location", curve_location,
                                    "as projection is oblique or co-linear with centroid curve")
                    continue
                dir2 = normalize(cross(dirp, dir1))
                dir3 = cross(dir1, dir2)
                anterior_direction = add(mult(dir2, wt2), mult(dir3, wt3))
                for idx, orientation_location in enumerate(orientation_locations):
                    if curve_location < orientation_location:
                        break
                else:
                    idx = len(orientation_locations)
                orientation_x.insert(idx, x)
                orientation_directions.insert(idx, anterior_direction)
                orientation_names.insert(idx, name)
                orientation_locations.insert(idx, curve_location)
        if not orientation_directions:
            logger.warning("Nerve: All orientation points ignored, using default orientation")
        else:
            twist_angle = 0.0
            orientation_twist_angles.append(twist_angle)  # top orientation point is at 0 radians
            direction = orientation_directions[0]
            for i in range(1, len(orientation_locations)):
                dir1, dir2, dir3 = track_curve_side_direction(
                    tx, td1, direction, orientation_locations[i - 1], orientation_locations[i])
                direction = orientation_directions[i]
                y = dot(direction, dir2)
                x = dot(direction, dir3)
                delta_twist_angle = -math.atan2(y, x)
                twist_angle += delta_twist_angle
                orientation_twist_angles.append(twist_angle)
            if debug_print:
                print("Orientation data:")
                for location, direction, twist_angle, name, in zip(
                        orientation_locations, orientation_directions, orientation_twist_angles, orientation_names):
                    print(location, direction, twist_angle, name)

    if orientation_twist_angles:
        data_identifier = 1
        twist_angle_field_name = "twist angle"
        twist_points_count_calibration_factor = len(orientation_twist_angles) / 27.0
        data_identifier = generate_datapoints(
            fit_region, orientation_x, data_identifier,
            field_names_and_values=[(twist_angle_field_name, orientation_twist_angles)])

        gradient1_penalty = 100.0 * twist_points_count_calibration_factor * length_calibration_factor
        gradient2_penalty = 1.0E+8 * twist_points_count_calibration_factor * (length_calibration_factor ** 3)
        define_and_fit_field(fit_region, "coordinates", "coordinates", twist_angle_field_name,
                             gradient1_penalty, gradient2_penalty, group_name=trunk_group_name)
        twist_angle = fieldmodule.findFieldByName(twist_angle_field_name).castFiniteElement()
        # extract fitted twist angle parameters from nodes
        ax = []
        ad1 = []
        for n in range(trunk_elements_count + 1):
            node = nodes.findNodeByIdentifier(n + 1)
            fieldcache.setNode(node)
            result, x = twist_angle.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 1)
            result, d1 = twist_angle.getNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 1)
            ax.append(x)
            ad1.append(d1)
        if debug_print:
            print("Orientation fit:")
            for i in range(len(orientation_locations)):
                fit_twist_angle = evaluateScalarOnCurve(ax, ad1, orientation_locations[i])
                print(orientation_locations[i], orientation_twist_angles[i], "-->", fit_twist_angle)
            print("node, ax, ad1:")
            for i in range(len(ax)):
                print(i + 1, ax[i], ad1[i])

        # get best-fit anterior direction at first orientation point
        first_point = 0
        first_location = orientation_locations[first_point]
        dir1, dir2, dir3 = track_curve_side_direction(
            tx, td1, orientation_directions[first_point], first_location, first_location)
        first_twist_data = orientation_twist_angles[first_point]
        first_twist_radians = evaluateScalarOnCurve(ax, ad1, first_location)
        delta_twist_radians = first_twist_radians - first_twist_data
        first_direction = normalize(add(mult(dir3, math.cos(delta_twist_radians)),
                                        mult(dir2, -math.sin(delta_twist_radians))))
    else:
        ad1 = ax = [0.0] * (trunk_elements_count + 1)
        first_location = (0, 0.0)
        first_twist_radians = 0.0
        # default anterior direction
        first_direction = [0.0, 1.0, 0.0]

    # compute side/anterior directions and their rates of change
    td2 = []
    td12 = []
    td3 = []
    td13 = []
    # track from first orientation point along nodes back to start, then forward from first orientation point to end
    prev_twist_radians = first_twist_radians
    prev_direction = first_direction
    prev_location = first_location
    next_location = (prev_location[0], 0.0)
    forward = False
    node_indexes = []
    tmp_ta = []
    tmp_dta = []
    for n in range(trunk_elements_count + 1):
        dir1, dir2, dir3 = track_curve_side_direction(
            tx, td1, prev_direction, prev_location, next_location, forward)
        node_index = (next_location[0] + 1) if forward else next_location[0]
        x = tx[node_index]
        d1 = td1[node_index]
        rv = 2.0 * rx[node_index]  # span of box is double radius
        rd = 2.0 * rd1[node_index]  # as above
        next_twist_radians = ax[node_index]
        ad = ad1[node_index]
        if forward:
            delta_twist_radians = next_twist_radians - prev_twist_radians
        else:
            delta_twist_radians = prev_twist_radians - next_twist_radians
        ante = normalize(add(mult(dir3, math.cos(delta_twist_radians)),
                             mult(dir2, -math.sin(delta_twist_radians))))
        left = normalize(cross(ante, dir1))

        d2 = set_magnitude(left, rv)
        d12 = add(mult(left, rd), mult(ante, rv * ad))
        d3 = set_magnitude(ante, rv)
        d13 = add(mult(ante, rd), mult(left, -rv * ad))
        ix = n if forward else 0
        node_indexes.insert(ix, node_index)
        td2.insert(ix, d2)
        td12.insert(ix, d12)
        td3.insert(ix, d3)
        td13.insert(ix, d13)
        tmp_ta.insert(ix, next_twist_radians)
        tmp_dta.insert(ix, delta_twist_radians)

        prev_twist_radians = next_twist_radians
        prev_direction = ante
        prev_location = next_location
        if forward:
            next_location = (next_location[0] + 1, 1.0)
        elif next_location[0] > 0:
            next_location = (next_location[0] - 1, 0.0)
        else:
            prev_twist_radians = first_twist_radians
            prev_direction = first_direction
            prev_location = first_location
            next_location = (prev_location[0], 1.0)
            forward = True
    if debug_print:
        print("node, ta, dta, d3:")
        for i in range(len(tmp_ta)):
            print(i + 1, tmp_ta[i], tmp_dta[i], td3[i])

    # copy model to user-supplied region
    sir = fit_region.createStreaminformationRegion()
    srm = sir.createStreamresourceMemory()
    sir.setResourceDomainTypes(srm, Field.DOMAIN_TYPE_NODES | Field.DOMAIN_TYPE_MESH1D)
    fit_region.write(sir)
    result, buffer = srm.getBuffer()
    sir = region.createStreaminformationRegion()
    srm = sir.createStreamresourceMemoryBuffer(buffer)
    region.read(sir)

    return tx, td1, td2, td12, td3, td13, default_trunk_radius


def set_group_nodes_derivatives_orthogonal(d1, radius):
    """
    Create derivatives d2 and d3, given first derivative d1 and a constant radius r.
    All three vectors should be an orthogonal set.
    :param d1: d1 derivative.
    :param radius: Constant radius used as magnitude of vectors d2 and d3.
    :return: Derivatives d2 and d3.
    """
    def set_first_v(v, v1):
        return rejection(v, v1)

    def set_second_v(v1, v2):
        return cross(v1, v2)

    yz = [1.0, 0.0, 0.0]
    yx = [0.0, 1.0, 0.0]
    zx = [0.0, 0.0, 1.0]

    # initial guess
    # td2 = set_d2(d1[0], yz)
    # td3 = set_d3(d1[0], td2)
    # print(d1[0], td2, td3)
    # td2 = set_d2(d1[0], yx)
    # td3 = set_d3(d1[0], td2)
    # print(d1[0], td2, td3)
    # td2 = set_d2(d1[0], zx)
    # td3 = set_d3(d1[0], td2)
    # print(d1[0], td2, td3)

    d2 = []
    d3 = []
    for td1 in d1:
        td3 = set_first_v(yx, td1)
        td3 = set_magnitude(td3, radius)
        d3.append(td3)

        td2 = set_second_v(td1, td3)
        td2 = set_magnitude(td2, radius)
        d2.append(td2)

    return d2, d3
