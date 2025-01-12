"""
Generates a hermite x bilinear 1-D central line mesh for a vagus nerve with branches
"""
import os
import math
import tempfile

from cmlibs.maths.vectorops import add, cross, div, dot, magnitude, matrix_mult, matrix_inv, mult, rejection, \
    set_magnitude, sub
from cmlibs.utils.zinc.field import find_or_create_field_group, findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis, Elementfieldtemplate
from cmlibs.zinc.field import Field, FieldGroup
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, remapEftNodeValueLabelWithNodes, \
    setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, getCubicHermiteBasisDerivatives, \
    interpolateCubicHermite, interpolateHermiteLagrange, smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters

from scaffoldfitter.fitter import Fitter
from scaffoldfitter.fitterstepalign import FitterStepAlign
from scaffoldfitter.fitterstepfit import FitterStepFit

from scaffoldmaker.utils.read_vagus_data import load_vagus_data
from scaffoldmaker.annotation.vagus_terms import get_vagus_branch_term, get_vagus_marker_term


def get_vagus_marker_locations_list():
    # vagus markers location in material coordinates between 0 to 1
    termNameVagusLengthList = {
        # cervical region
        # "level of exiting brainstem on the vagus nerve": 0.0,  # note this term is not on the list of annotations
        "level of superior border of jugular foramen on the vagus nerve": 0.02762944,
        # # "level of inferior border of jugular foramen on the vagus nerve": 0.05351264,
        "level of inferior border of jugular foramen on the vagus nerve": 0.047304398,
        # "level of inferior border of cranium on the vagus nerve": 0.0588,
        # "level of C1 transverse process on the vagus nerve": 0.10276128,
        # "level of angle of the mandible on the vagus nerve": 0.135184,
        # # "level of greater horn of hyoid on the vagus nerve": 0.14595904,
        "level of greater horn of hyoid on the vagus nerve": 0.171654481,
        # "level of carotid bifurcation on the vagus nerve": 0.15474592,
        # "level of laryngeal prominence on the vagus nerve": 0.22029792,
        # thoracic region
        # # "level of superior border of the clavicle on the vagus nerve": 0.37620064,
        "level of superior border of the clavicle on the vagus nerve": 0.348953222,
        # "level of jugular notch on the vagus nerve": 0.39885024,
        # "level of carina": 0.47869728,  # not on the list of annotations yet!
        # "level of sternal angle on the vagus nerve": 0.48395264,
        # "level of 1 cm superior to start of esophageal plexus on the vagus nerve": 0.52988032,
        # abdominal region
        # "level of esophageal hiatus on the vagus nerve": 0.813852428,
        # "level of aortic hiatus on the vagus nerve": 0.9323824,
        # "level of end of trunk": 1.0  # note this term is also not on the list of annotations
    }
    return termNameVagusLengthList


def is_bony_landmark(name):
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
    if name in bony_landmarks_names:
        return True
    return False


class MeshType_3d_vagus_box1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @staticmethod
    def getName():
        return "3D Vagus Box 1"

    @staticmethod
    def getParameterSetNames():
        return ['Human Trunk 1']

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
        value_labels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                        Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = findOrCreateFieldCoordinates(fieldmodule).castFiniteElement()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for value_label in value_labels[1:]:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, value_label, 1)

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
        elementtemplate_branch_root = mesh3d.createElementtemplate()
        elementtemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate_branch_root.defineField(coordinates, -1, eft3dNV)

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
        linetemplate_branch_root = mesh1d.createElementtemplate()
        linetemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linetemplate_branch_root.defineField(coordinates, -1, eft1dNV)

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
        facetemplate_and_eft_list_branch_root = [None] * 4
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
            facetemplate_branch_root = mesh2d.createElementtemplate()
            facetemplate_branch_root.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            facetemplate_branch_root.defineField(coordinates, -1, eft2dNV)
            facetemplate_and_eft_list_branch_root[e] = (facetemplate_branch_root, eft2dNV)

        # load data from file
        print('Extracting data...')
        vagus_data = load_vagus_data(region)
        marker_data = vagus_data.get_level_markers()
        trunk_data = vagus_data.get_trunk_coordinates()
        assert len(marker_data) >= 2, f"At least two landmarks are expected in the data. Incomplete data."
        assert len(trunk_data) >= 2, f"At least two trunk points are required."

        if any(['level of esophageal hiatus' in marker_name or
                'level of aortic hiatus' in marker_name for marker_name in marker_data.keys()]):
            is_full_vagus = True
            print('Estimate_trunk: Vagus top to bottom')
        else:
            is_full_vagus = False
            print('Estimate_trunk: Vagus top to esophageal plexus')

        # settings for radius (assume equal sides) and material coordinates
        vagus_aspect_ratio = 0.005  # assuming vagus approx diameter (5mm) / vagus length (85mm)
        branch_to_trunk_ratio = 0.5
        vagus_trunk_length = estimate_real_trunk_length(trunk_data)
        radius = math.ceil(vagus_aspect_ratio * vagus_trunk_length / 2)
        branch_radius = branch_to_trunk_ratio * radius

        rescaled_vagus_trunk_length = estimate_parametrised_vagus_length(is_full_vagus)
        vagus_radius = vagus_aspect_ratio * rescaled_vagus_trunk_length / 2
        vagus_branch_radius = branch_to_trunk_ratio * vagus_radius

        # evaluate & fit centroid lines for trunk and branches
        print('Building centerlines for scaffold...')
        fit_region, marker_fit_groups, branches_order, \
            branch_root_parameters = evaluate_vagus_1d_coordinates(region, vagus_data, is_full_vagus, options)
        fit_fieldmodule = fit_region.getFieldmodule()
        fit_fieldcache = fit_fieldmodule.createFieldcache()
        fit_coordinates = fit_fieldmodule.findFieldByName("coordinates").castFiniteElement()
        fit_nodes = fit_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        annotation_groups = []
        annotation_term_map = vagus_data.get_annotation_term_map()
        # vagus annotation groups
        vagusCentroidGroup = AnnotationGroup(region, ("Vagus centroid", ""))
        annotation_groups.append(vagusCentroidGroup)
        vagusCentroidMeshGroup = vagusCentroidGroup.getMeshGroup(mesh1d)

        vagusEpineuriumAnnotationGroup = AnnotationGroup(region, ("Vagus epineurium", ""))
        annotation_groups.append(vagusEpineuriumAnnotationGroup)
        vagusEpineuriumMeshGroup = vagusEpineuriumAnnotationGroup.getMeshGroup(mesh2d)

        node_map = {}
        print('Building trunk...')

        # read trunk nodes
        trunk_group_name = vagus_data.get_trunk_group_name()
        fit_trunk_group = find_or_create_field_group(fit_fieldmodule, trunk_group_name)
        fit_trunk_nodes = fit_trunk_group.getNodesetGroup(fit_nodes)
        elements_along_trunk = fit_trunk_nodes.getSize()
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
        sd2, sd3 = set_group_nodes_derivatives_orthogonal(sd1, radius)
        sd12, sd13 = smoothCurveSideCrossDerivatives(sx, sd1, [sd2, sd3])
        trunk_d3 = sd3[0] # use for branch orientation

        # trunk annotation groups
        trunk_box_group = AnnotationGroup(region, (trunk_group_name, annotation_term_map[trunk_group_name]))
        annotation_groups.append(trunk_box_group)
        trunk_box_mesh_group = trunk_box_group.getMeshGroup(mesh3d)

        node_identifier = 1
        line_identifier = 1
        element_identifier = 1
        face_identifier = 1

        # create nodes
        for n in range(elements_along_trunk):
            node_map[sn[n]] = node_identifier
            node = nodes.createNode(node_identifier, nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, sd2[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3[n])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[n])
            node_identifier += 1

        # create elements
        for n in range(1, elements_along_trunk):
            node_id = n + 1
            nids = [node_id - 1, node_id]

            line = mesh1d.createElement(line_identifier, linetemplate)
            line.setNodesByIdentifier(eft1d, nids)
            vagusCentroidMeshGroup.addElement(line)
            line_identifier += 1

            element = mesh3d.createElement(element_identifier, elementtemplate)
            element.setNodesByIdentifier(eft3d, nids)
            element.setScaleFactors(eft3d, [-1.0])
            trunk_box_mesh_group.addElement(element)
            element_identifier += 1

            for e in range(4):
                facetemplate, eft2d = facetemplate_and_eft_list[e]
                face = mesh2d.createElement(face_identifier, facetemplate)
                face.setNodesByIdentifier(eft2d, nids)
                face.setScaleFactors(eft2d, scalefactors2d)
                vagusEpineuriumMeshGroup.addElement(face)
                face_identifier += 1

        print('Building branches...')
        branch_parent_map = vagus_data.get_branch_parent_map()
        for branch_name in branches_order:
            # read nodes
            fit_branch_group = find_or_create_field_group(fit_fieldmodule, branch_name)
            fit_branch_nodes = fit_branch_group.getNodesetGroup(fit_nodes)
            elements_along_branch = fit_branch_nodes.getSize()
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

            # branch root parameters (assume d12, d13 are both zero for now)
            trunk_segment_start_id = node_map[branch_root_parameters[branch_name][0]]
            trunk_segment_end_id = node_map[branch_root_parameters[branch_name][1]]
            branch_root_xi = branch_root_parameters[branch_name][2]
            bx = branch_root_parameters[branch_name][3]
            bd1 = branch_root_parameters[branch_name][4]

            # temporarily add branch root parameters to node parameters
            sx.insert(0, bx)
            sd1.insert(0, bd1)

            # calculate side and cross derivatives - d2, d3, d12, d13
            sd2, sd3 = set_group_nodes_derivatives_orthogonal_branches(sd1, trunk_d3, branch_radius)
            sd12, sd13 = smoothCurveSideCrossDerivatives(sx, sd1, [sd2, sd3])

            # get & remove branch root parameters from node parameters
            bd2 = sd2[0]
            bd3 = sd3[0]
            bd12 = sd12[0]
            bd13 = sd13[0]

            sx.pop(0)
            sd1.pop(0)
            sd2.pop(0)
            sd3.pop(0)
            sd12.pop(0)
            sd13.pop(0)

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
            basis_to = [bd1, bd2, bd3]
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
            branch_box_group = AnnotationGroup(region, (branch_name, annotation_term_map[branch_name]))
            annotation_groups.append(branch_box_group)
            branch_box_mesh_group = branch_box_group.getMeshGroup(mesh3d)

            # create nodes (excluding branch root)
            for n in range(elements_along_branch - 1):
                node_map[sn[n]] = node_identifier
                node = nodes.createNode(node_identifier, nodetemplate)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, sd2[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[n])

                if n == 0:
                    # branch root element
                    nids = [trunk_segment_start_id, node_identifier, trunk_segment_end_id]

                    scalefactors = []
                    scalefactors.extend(scalefactorsX)
                    scalefactors.extend(scalefactorsD1)

                    line = mesh1d.createElement(line_identifier, linetemplate_branch_root)
                    line.setNodesByIdentifier(eft1dNV, nids)
                    line.setScaleFactors(eft1dNV, list(scalefactors))
                    vagusCentroidMeshGroup.addElement(line)
                    line_identifier += 1

                else:
                    # other elements
                    nids = [node_identifier - 1, node_identifier]

                    line = mesh1d.createElement(line_identifier, linetemplate)
                    line.setNodesByIdentifier(eft1d, nids)
                    vagusCentroidMeshGroup.addElement(line)
                    line_identifier += 1

                node_identifier += 1

            # create elements
            for n in range(0, elements_along_branch - 1):
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

                    element = mesh3d.createElement(element_identifier, elementtemplate_branch_root)
                    element.setNodesByIdentifier(eft3dNV, nids)
                    element.setScaleFactors(eft3dNV, scalefactors)
                    branch_box_mesh_group.addElement(element)
                    element_identifier += 1

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

                        facetemplate_branch_root, eft2dNV = facetemplate_and_eft_list_branch_root[e]
                        face = mesh2d.createElement(face_identifier, facetemplate_branch_root)
                        face.setNodesByIdentifier(eft2dNV, nids)
                        face.setScaleFactors(eft2dNV, scalefactors)
                        vagusEpineuriumMeshGroup.addElement(face)
                        face_identifier += 1

                else:
                    nids = [node_id - 1, node_id]

                    element = mesh3d.createElement(element_identifier, elementtemplate)
                    element.setNodesByIdentifier(eft3d, nids)
                    element.setScaleFactors(eft3d, [-1.0])
                    branch_box_mesh_group.addElement(element)
                    element_identifier += 1

                    for e in range(4):
                        facetemplate, eft2d = facetemplate_and_eft_list[e]
                        face = mesh2d.createElement(face_identifier, facetemplate)
                        face.setNodesByIdentifier(eft2d, nids)
                        face.setScaleFactors(eft2d, scalefactors2d)
                        vagusEpineuriumMeshGroup.addElement(face)
                        face_identifier += 1

            # remove trunk nodes from branch group
            parent_field_group = find_or_create_field_group(fieldmodule, branch_parent_map[branch_name])
            branch_nodeset_group = branch_box_group.getNodesetGroup(nodes)
            if branch_nodeset_group.isValid():
                branch_nodeset_group.removeNodesConditional(parent_field_group)

        # set markers
        print('Adding anatomical landmarks...')
        for marker_group in marker_fit_groups:
            marker_name = marker_group.getName()
            annotationGroup = findOrCreateAnnotationGroupForTerm(annotation_groups, region,
                                                                 get_vagus_marker_term(marker_name),
                                                                 isMarker=True)
            element, xi = marker_group.getMarkerLocation()
            annotationGroup.createMarkerNode(node_identifier,
                                             element=mesh3d.findElementByIdentifier(element.getIdentifier()),
                                             xi=[xi[0], 0.5, 0.5])
            node_identifier += 1

        # add material coordinates
        print('Adding material coordinates...')
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
        rescaled_step = rescaled_vagus_trunk_length / (elements_along_trunk - 1)

        elem_iter = mesh3d.createElementiterator()
        element = elem_iter.next()
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
        marker_names = fieldmodule.findFieldByName("marker_name")
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

        print('Done\n')

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
    """
    return: squared scalar magnitude of vector v
    """

    # TODO: proposed function to cmlibs.maths
    return sum(c * c for c in v)


def distance_squared(u, v):
    """
    return: squared distance to avoid the cost of a square root
    """

    # TODO: proposed function to cmlibs.maths
    return sum((u_i - v_i) ** 2 for u_i, v_i in zip(u, v))


def find_dataset_endpoints(points):
    """
    Given list of XYZ coordinates, find two furthest apart from each other.
    return: a list of two coordinates.
    """

    # Step 1: Deterministically select points for the sample (every step_size-th point)
    step_size = 1 if len(points) < 100 else 250
    sampled_indices = list(range(0, len(points), step_size))
    if (len(points) - 1) not in sampled_indices:
        sampled_indices.append(len(points) - 1)

    # Step 2: Find the furthest points from within the sampled subset
    max_distance = 0
    furthest_pair = (sampled_indices[0], sampled_indices[1])

    for i in range(len(sampled_indices)):
        for j in range(i + 1, len(sampled_indices)):
            dist = distance_squared(points[sampled_indices[i]], points[sampled_indices[j]])
            if dist > max_distance:
                max_distance = dist
                furthest_pair = (sampled_indices[i], sampled_indices[j])

    # Step 3: Compare sampled endpoints against all points, maintaining original order
    for i in range(len(points)):
        for endpoint_index in furthest_pair:
            if i < endpoint_index:
                dist = distance_squared(points[i], points[endpoint_index])
                if dist > max_distance:
                    max_distance = dist
                    furthest_pair = (i, endpoint_index)

    endpoints = [points[furthest_pair[0]], points[furthest_pair[1]]]
    return endpoints


def find_point_projection_relative_to_segment(point, segment_start, segment_end):
    """
    return: the position of the projection relative to the line segment, defined by segment_start and segment_end
    """

    ap = sub(point, segment_start)
    ab = sub(segment_end, segment_start)
    projection_scalar = dot(ap, ab) / dot(ab, ab) if dot(ab, ab) != 0 else 0
    return projection_scalar


def set_group_nodes_derivatives_orthogonal(d1, radius):
    """

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


def set_group_nodes_derivatives_orthogonal_branches(d1, trunk_d3, radius):
    """

    """
    def set_first_v(v, v1):
        return rejection(v, v1)

    def set_second_v(v1, v2):
        return cross(v1, v2)

    d2 = []
    d3 = []
    for td1 in d1:
        td3 = set_first_v(trunk_d3, td1)
        td3 = set_magnitude(td3, radius)
        d3.append(td3)

        td2 = set_second_v(td3, td1)
        td2 = set_magnitude(td2, radius)
        d2.append(td2)

    return d2, d3


def evaluate_vagus_1d_coordinates(region, vagus_data, is_full_vagus, options):
    """
    Generates an intermediate hermite x bilinear 1-D central line mesh for a vagus nerve with branches, markers.
    First, it evaluates trunk coordinates based on top and bottom of the supplied markers, then fits the coordinates to
    supplied trunk data.
    Second, it creates branches based on the
    """

    elements_along_trunk = options['Number of elements along the trunk']
    iterations_number = options['Iterations (fit trunk)']

    fit_region = region.createRegion()
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
    trunk_group_name = vagus_data.get_trunk_group_name()
    trunkCentroidGroup = find_or_create_field_group(fit_fm, trunk_group_name)
    trunkCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
    trunkCentroidMeshGroup = trunkCentroidGroup.getOrCreateMeshGroup(fit_mesh1d)

    # used for fitting only
    trunkFitCentroidGroup = find_or_create_field_group(fit_fm, trunk_group_name + '-fit')
    trunkFitCentroidGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
    trunkFitCentroidMeshGroup = trunkFitCentroidGroup.getOrCreateMeshGroup(fit_mesh1d)

    trunk_data = vagus_data.get_trunk_coordinates()
    marker_data = vagus_data.get_level_markers()

    annotation_groups = []

    for ii in range(iterations_number):

        if ii == 0:
            trunk_data_endpoints = find_dataset_endpoints([trunk_pt[0] for trunk_pt in trunk_data])
            marker_locations, tx, td1, \
                step, element_length = estimate_trunk_coordinates(elements_along_trunk, marker_data,
                                                                  trunk_data_endpoints, is_full_vagus)
        else:
            # read tx from fit_coordinates
            _, node_field_parameters = get_nodeset_field_parameters(fit_nodes, fit_coordinates, value_labels)
            tx = [nodeParameter[1][0][0] for nodeParameter in node_field_parameters]
            td1 = [nodeParameter[1][1][0] for nodeParameter in node_field_parameters]

        trunk_nodes_data_bounds = estimate_trunk_data_boundaries(tx, elements_along_trunk, trunk_data_endpoints)
        print(trunk_nodes_data_bounds)

        node_identifier = 1
        line_identifier = 1

        nodes_before = []
        nodes_after = []
        for n in range(elements_along_trunk):
            sx = tx[n]
            sd1 = td1[n]

            if ii == 0:
                node = fit_nodes.createNode(node_identifier, fit_nodetemplate)
            else:
                node = fit_nodes.findNodeByIdentifier(node_identifier)
            fit_fc.setNode(node)
            fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, sx)
            fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

            # add to trunk group used for data fitting
            if trunk_nodes_data_bounds[0] <= node_identifier <= trunk_nodes_data_bounds[-1]:
                pass
            elif node_identifier < trunk_nodes_data_bounds[0]:
                nodes_before.append(node_identifier)
            else:
                nodes_after.append(node_identifier)

            if n > 0:
                nids = [node_identifier - 1, node_identifier]
                if ii == 0:
                    line = fit_mesh1d.createElement(line_identifier, fit_eltemplate)
                else:
                    line = fit_mesh1d.findElementByIdentifier(line_identifier)
                line.setNodesByIdentifier(fit_eft1d, nids)
                trunkCentroidMeshGroup.addElement(line)
                # add element to trunk group used for data fitting
                if node_identifier - 1 >= trunk_nodes_data_bounds[0] and node_identifier <= trunk_nodes_data_bounds[-1]:
                    trunkFitCentroidMeshGroup.addElement(line)
                line_identifier += 1
            node_identifier += 1

        if ii == 0:
            # set markers
            for marker_name, marker_location in marker_locations.items():
                print(marker_name)
                annotation_group = findOrCreateAnnotationGroupForTerm(annotation_groups, fit_region,
                                                                      get_vagus_marker_term(marker_name),
                                                                      isMarker=True)
                annotation_group.createMarkerNode(node_identifier,
                                                  element=fit_mesh1d.findElementByIdentifier(marker_location[0]),
                                                  xi=[marker_location[1]])
                node_identifier += 1
        else:
            node_identifier += len(marker_data)

        # create temporary model file
        sir = fit_region.createStreaminformationRegion()
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            fitter_model_file = temp_file.name
            srf = sir.createStreamresourceFile(fitter_model_file)
            fit_region.write(sir)

        print('... Fitting trunk, iteration', str(ii + 1))
        fitter_data_file = vagus_data.get_datafile_path()
        if is_full_vagus:
            # non-REVA data (current scaffold based on Japanese dataset)
            fitter = fit_full_trunk_model(fitter_model_file, fitter_data_file, trunk_group_name + '-fit')
        else:
            # REVA data
            fitter = fit_trunk_model(fitter_model_file, fitter_data_file, trunk_group_name + '-fit')
        set_fitted_group_nodes(fit_region, fitter, trunk_group_name + '-fit')

        # remove temporary model file
        os.remove(fitter_model_file)

        if len(nodes_before) > 0 or len(nodes_after) > 0:
            # calculate average derivative d1 along the fitted vagus trunk
            trunk_fit_nodes = trunkFitCentroidGroup.getNodesetGroup(fit_nodes)
            trunk_fit_size = trunk_fit_nodes.getSize()

            node_iter = fit_nodes.createNodeiterator()
            node = node_iter.next() # ignore the first node
            avg_d1 = [0, 0, 0]
            nid = 2
            while node.isValid() and nid < trunk_fit_size:
                fit_fc.setNode(node)
                _, ld1 = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                avg_d1 = add(avg_d1, ld1)
                node = node_iter.next()
                nid += 1
            avg_d1 = [dim / trunk_fit_size for dim in avg_d1]
            avg_d1 = set_magnitude(avg_d1, element_length)

            if len(nodes_before) > 0:
                # start unfitted nodes from the first fitted node coordinate
                node_id = trunk_nodes_data_bounds[0]
                node = fit_nodes.findNodeByIdentifier(node_id)
                fit_fc.setNode(node)
                _, lx = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, avg_d1)

                node_count = 1
                for i in range(len(nodes_before) - 1, -1, -1):
                    node_id = nodes_before[i]
                    x = [lx[j] - node_count * avg_d1[j] for j in range(3)]

                    node = fit_nodes.findNodeByIdentifier(node_id)
                    fit_fc.setNode(node)
                    fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, avg_d1)
                    node_count += 1

            if len(nodes_after) > 0:
                # start unfitted nodes from the last fitted node coordinate
                node_id = trunk_nodes_data_bounds[-1]
                node = fit_nodes.findNodeByIdentifier(node_id)
                fit_fc.setNode(node)
                _, lx = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, avg_d1)

                node_count = 1
                for i in range(len(nodes_after)):
                    node_id = nodes_after[i]
                    x = [lx[j] + node_count * avg_d1[j] for j in range(3)]

                    node = fit_nodes.findNodeByIdentifier(node_id)
                    fit_fc.setNode(node)
                    fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, avg_d1)
                    node_count += 1

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

    visited_branches_order = []
    branch_root_parameters = {}

    print('... Adding branches')
    branch_data = vagus_data.get_branch_data()
    branch_parent_map = vagus_data.get_branch_parent_map()
    queue = [branch for branch in branch_parent_map.keys() if branch_parent_map[branch] == trunk_group_name]
    while queue:
        branch_name = queue.pop(0)
        print(branch_name)

        if branch_name in visited_branches_order:
            continue
        visited_branches_order.append(branch_name)
        queue.extend([branch for branch in branch_parent_map.keys() if branch_parent_map[branch] == branch_name])

        branch_coordinates = [branch_node[0] for branch_node in branch_data[branch_name]]
        branch_parent_name = branch_parent_map[branch_name]
        #print(branch_name, ' -> ', branch_parent_name)

        # determine branch approximate start and closest trunk node index
        bx, bd1, parent_s_nid, parent_f_nid, branch_root_xi, elements_along_branch = \
            estimate_branch_coordinates(fit_region, branch_coordinates, element_length, branch_parent_name)
        branch_root_parameters[branch_name] = [parent_s_nid, parent_f_nid, branch_root_xi, bx[0]]
        # print('  branch between nodes:', parent_s_nid, parent_f_nid, 'at loc =', branch_root_xi)

        branch_centroid_group = find_or_create_field_group(fit_fm, branch_name)
        branch_centroid_group.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
        branch_centroid_mesh_group = branch_centroid_group.getOrCreateMeshGroup(fit_mesh1d)

        for n in range(elements_along_branch):
            sx = bx[n]
            sd1 = bd1

            if n == 0:
                # create branch special node
                node = fit_nodes.createNode(node_identifier, fit_nodetemplateBranchRoot)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)
            else:
                node = fit_nodes.createNode(node_identifier, fit_nodetemplate)
                fit_fc.setNode(node)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_VALUE, 1, sx)
                fit_coordinates.setNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, sd1)

                if n == 1:
                    # create branch root element
                    nids = [node_identifier - 1, node_identifier,
                            parent_s_nid, parent_f_nid]
                    line = fit_mesh1d.createElement(line_identifier, fit_eltemplateBranchRoot)
                    line.setNodesByIdentifier(fit_eft1dBranchRoot, nids)
                    scalefactorsNV = getCubicHermiteBasis(branch_root_xi)
                    line.setScaleFactors(fit_eft1dBranchRoot, list(scalefactorsNV))
                    branch_centroid_mesh_group.addElement(line)
                    line_identifier += 1
                else:
                    nids = [node_identifier - 1, node_identifier]
                    line = fit_mesh1d.createElement(line_identifier, fit_eltemplate)
                    line.setNodesByIdentifier(fit_eft1d, nids)
                    branch_centroid_mesh_group.addElement(line)
                    line_identifier += 1
            node_identifier += 1

        # remove trunk nodes from branch group
        parent_group = find_or_create_field_group(fit_fm, branch_parent_name)
        branch_nodeset_group = branch_centroid_group.getNodesetGroup(fit_nodes)
        if branch_nodeset_group.isValid():
            branch_nodeset_group.removeNodesConditional(parent_group)

        # create temporary model file
        sir = fit_region.createStreaminformationRegion()
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            fitter_model_file = temp_file.name
            srf = sir.createStreamresourceFile(fitter_model_file)
            fit_region.write(sir)

        # print('fitting %s' % branch_name)
        fitter = fit_branches_model(fitter_model_file, fitter_data_file, branch_name)
        set_fitted_group_nodes(fit_region, fitter, branch_name)

        # remove temporary model file
        os.remove(fitter_model_file)

        # extract first branch node - d1 fitted value
        branch_group = find_or_create_field_group(fit_fm, branch_name)
        branch_nodes = branch_group.getNodesetGroup(fit_nodes)
        node_iter = branch_nodes.createNodeiterator()
        node = node_iter.next()
        fit_fc.setNode(node)
        _, sd1 = fit_coordinates.getNodeParameters(fit_fc, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        branch_root_parameters[branch_name].append(sd1)

    # evaluate marker geometric coordinates
    # print('----')
    # for annotationGroup in annotationGroups:
    #     xyz = annotationGroup.evaluateMarkerMaterialCoordinatesFromElementXi(fit_coordinates)
    #     element, xi = annotationGroup.getMarkerLocation()
    #     print(annotationGroup.getName(), element.getIdentifier(), xi, xyz)
    # print('----')

    # remove temporary data file
    fitter_data_file = vagus_data.get_datafile_path()
    vagus_data.reset_datafile_path()
    os.remove(fitter_data_file)

    return fit_region, annotation_groups, visited_branches_order, branch_root_parameters


def estimate_trunk_coordinates(elements_along_trunk, raw_marker_data, trunk_data_endpoints, is_full_vagus):
    """

    """

    # choose markers for building initial scaffold
    # at the moment uses the first and the last markers in the data
    term_name_vagus_length_list = get_vagus_marker_locations_list()
    total_vagus_length = estimate_parametrised_vagus_length(is_full_vagus)
    step = total_vagus_length / (elements_along_trunk - 1)
    element_length = magnitude(sub(trunk_data_endpoints[0], trunk_data_endpoints[-1])) / (elements_along_trunk - 1)

    # filter out markers not chosen for fitting
    marker_data = {k: v for k, v in raw_marker_data.items()
                   if k.replace('left ', '', 1).replace('right ', '', 1) in term_name_vagus_length_list.keys()}

    use_markers = [list(marker_data.keys())[0],
                   list(marker_data.keys())[-1]]

    pts = []
    params = []

    # calculate markers xi locations
    marker_locations = {}
    for marker_name in marker_data.keys():
        use_marker_name = marker_name.replace('left ', '', 1).replace('right ', '', 1)
        if use_marker_name in term_name_vagus_length_list:
            marker_material_coord = term_name_vagus_length_list[use_marker_name]
            tmp = math.modf(marker_material_coord / step)
            marker_locations[marker_name] = (int(tmp[1] + 1), tmp[0])  # (element, xi)
            # print(marker_name, marker_material_coord, marker_locations[marker_name])

            if marker_name in use_markers:
                pts.append(marker_data[marker_name])
                params.append(term_name_vagus_length_list[use_marker_name])

    # calculate trunk initial coordinates
    trunk_coordinates = []
    trunk_d1 = []
    dx, dy, dz = [(pts[1][dim] - pts[0][dim]) / (params[1] - params[0]) for dim in range(3)]

    for i in range(elements_along_trunk):
        trunk_coordinates.append([pts[0][0] + dx * (i * step - params[0]),
                                  pts[0][1] + dy * (i * step - params[0]),
                                  pts[0][2] + dz * (i * step - params[0])])
        trunk_d1.append([dx * step, dy * step, dz * step])

    return marker_locations, trunk_coordinates, trunk_d1, step, element_length


def estimate_parametrised_vagus_length(is_full_vagus):
    if is_full_vagus:
        total_vagus_length = 1.0
    else:
        # approximate parametrised length top to esophageal plexus
        total_vagus_length = 0.6
    return total_vagus_length


def estimate_real_trunk_length(trunk_data):
    trunk_data_endpoints = find_dataset_endpoints([trunk_pt[0] for trunk_pt in trunk_data])
    return magnitude(sub(trunk_data_endpoints[0],trunk_data_endpoints[-1]))


def estimate_trunk_data_boundaries(trunk_nodes, elements_along_trunk, trunk_data_endpoints):
    """
    """

    trunk_nodes_data_bounds = []
    for ind, endpoint in enumerate(trunk_data_endpoints):
        param = find_point_projection_relative_to_segment(endpoint, trunk_nodes[0], trunk_nodes[-1])
        if param < 0:
            nearby_node = 1
        elif param > 1:
            nearby_node = elements_along_trunk
        else:
            nearby_node = param * (elements_along_trunk - 1) + 1

        if ind == 0:
            # trunk node near the start of data
            trunk_nodes_data_bounds.append(math.floor(nearby_node))
        else:
            # trunk node near the end of data
            trunk_nodes_data_bounds.append(math.ceil(nearby_node))

    return trunk_nodes_data_bounds


def estimate_branch_coordinates(region, branch_coordinates, element_length, branch_parent_name):
    """

    """

    fm = region.getFieldmodule()
    fieldcache = fm.createFieldcache()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    branch_start_x, branch_end_x, parent_s_nid, parent_f_nid = find_branch_start_segment(region, branch_coordinates,
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
        dsq_a = distance_squared(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_a))
        dsq_b = distance_squared(branch_start_x, interpolateCubicHermite(px_1, pd1_1, px_2, pd1_2, xi_b))
        if dsq_a >= dsq_b:
            xi_a = (xi_a + xi_b) / 2
        else:
            xi_b = (xi_a + xi_b) / 2
    branch_root_xi = (xi_a + xi_b) / 2

    # recalculate branch start parameters
    branch_start_x = interpolateHermiteLagrange(px_1, pd1_1, px_2, branch_root_xi)
    branch_length = magnitude(sub(branch_end_x, branch_start_x))
    elements_along_branch = math.floor(branch_length / element_length) + 1
    if elements_along_branch < 3:
        # need to have at least 3 nodes (including fictionary node) to allow child branches
        elements_along_branch = 3
    if elements_along_branch > 10:
        elements_along_branch = 10

    branch_coordinates = []
    dx, dy, dz = div(sub(branch_end_x, branch_start_x), (elements_along_branch - 1))
    for i in range(elements_along_branch):
        branch_coordinates.append([branch_start_x[0] + dx * i,
                                   branch_start_x[1] + dy * i,
                                   branch_start_x[2] + dz * i])

    return branch_coordinates, [dx, dy, dz], parent_s_nid, parent_f_nid, branch_root_xi, elements_along_branch


def find_branch_start_segment(region, branch_coordinates, parent_group_name):
    """

    """

    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

    parent_group = find_or_create_field_group(fm, parent_group_name)
    parent_nodes = parent_group.getNodesetGroup(nodes)
    _, group_parameters = get_nodeset_field_parameters(parent_nodes, coordinates, [Node.VALUE_LABEL_VALUE])
    parent_nodes_ids = [parameter[0] for parameter in group_parameters]
    parent_nodes_x = [parameter[1][0][0] for parameter in group_parameters]

    # find branch ends in data
    branch_ends_points = find_dataset_endpoints(branch_coordinates)

    # find closest to branch distance, branch start, group node closest to the branch
    min_dsq = float('inf')
    for i in range(len(parent_nodes_x)):
        node_x = parent_nodes_x[i]
        if node_x is None:
            continue

        for branch_point in branch_ends_points:
            dist = distance_squared(node_x, branch_point)
            if dist <= min_dsq:
                min_dsq = dist
                branch_start = branch_point
                closest_index = i

    # determine segment closest to branch (previous or next to the node)
    if closest_index == 0:
        parent_start_index = closest_index
    elif closest_index == len(parent_nodes_x) - 1:
        parent_start_index = closest_index - 1
    else:
        proj_before = find_point_projection_relative_to_segment(branch_start,
                                                                parent_nodes_x[closest_index - 1],
                                                                parent_nodes_x[closest_index])
        proj_after = find_point_projection_relative_to_segment(branch_start,
                                                               parent_nodes_x[closest_index],
                                                               parent_nodes_x[closest_index + 1])
        if 0 <= proj_before <= 1:
            parent_start_index = closest_index - 1
        elif 0 <= proj_after <= 1:
            parent_start_index = closest_index
        elif abs(proj_before) < abs(proj_after):
            parent_start_index = closest_index - 1
        else:
            parent_start_index = closest_index

    parent_s_node_id = parent_nodes_ids[parent_start_index]
    parent_f_node_id = parent_nodes_ids[parent_start_index + 1]
    branch_end = branch_ends_points[0] if branch_ends_points[1] == branch_start else branch_ends_points[1]

    return branch_start, branch_end, parent_s_node_id, parent_f_node_id


# fitter functions
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

    # updated temporary trunk fitting workflow
    # fit step 1
    fit1 = FitterStepFit()
    fitter.addFitterStep(fit1)
    fit1.setGroupDataWeight('marker', [80.0])
    fit1.setGroupDataSlidingFactor('marker', 0.01)
    fit1.setGroupStrainPenalty(None, [5000.0])
    fit1.setGroupCurvaturePenalty(None, [0.0])
    fit1.setNumberOfIterations(1)

    # fit step 2
    fit2 = FitterStepFit()
    fitter.addFitterStep(fit2)
    fit2.setGroupDataWeight('marker', [70.0])
    fit2.setGroupDataSlidingFactor('marker', None)
    fit2.setGroupStrainPenalty(None, [2500.0])
    fit2.setGroupCurvaturePenalty(None, [5.0])
    fit2.setNumberOfIterations(3)

    # fit step 3
    fit3 = FitterStepFit()
    fitter.addFitterStep(fit3)
    fit3.setGroupDataWeight('marker', None)
    fit3.setGroupStrainPenalty(None, [1500.0])
    fit3.setGroupCurvaturePenalty(None, [15.0])
    fit3.setNumberOfIterations(5)

    fitter.run()

    #rms_trunk_error, max_trunk_error = fitter.getDataRMSAndMaximumProjectionErrorForGroup('left vagus X nerve trunk')
    #rms_marker_error, max_marker_error = fitter.getDataRMSAndMaximumProjectionErrorForGroup('marker')

    return fitter


def fit_full_trunk_model(modelfile, datafile, trunk_group_name=None):
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

    # trunk fitting workflow used for japanese dataset
    # fit step 1
    fit1 = FitterStepFit()
    fitter.addFitterStep(fit1)
    fit1.setGroupDataWeight('marker', [100.0])
    fit1.setGroupDataSlidingFactor('marker', 0.01)
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setNumberOfIterations(10)
    fit1.setUpdateReferenceState(True)

    # fit step 2
    fit2 = FitterStepFit()
    fitter.addFitterStep(fit2)
    fit2.setGroupDataWeight('marker', [100.0])
    fit2.setGroupDataSlidingFactor('marker', 0.01)
    fit2.setGroupStrainPenalty(None, [5.0])
    fit2.setGroupCurvaturePenalty(None, [100.0])
    fit2.setNumberOfIterations(5)
    fit2.setUpdateReferenceState(True)

    fitter.run()

    return fitter


def fit_branches_model(modelfile, datafile, branch_name=None):
    """
    Initialises scaffold fitter and runs through its steps for a given vagus nerve branch.
    """

    # initial configuration
    fitter = Fitter(modelfile, datafile)
    fitter.load()
    fitter.setModelCoordinatesFieldByName('coordinates')
    if branch_name:
        fitter.setModelFitGroupByName(branch_name)
    fitter.setFibreField(fitter.getFieldmodule().findFieldByName("zero fibres"))
    fitter.setDataCoordinatesFieldByName('coordinates')
    fitter.setMarkerGroupByName('marker')  # not necessary, it's marker by default
    fitter.setDiagnosticLevel(0)

    # fit step 1
    fit1 = FitterStepFit()
    fitter.addFitterStep(fit1)
    fit1.setGroupStrainPenalty(None, [15.0])
    fit1.setGroupCurvaturePenalty(None, [50.0])
    fit1.setNumberOfIterations(5)
    fit1.setUpdateReferenceState(True)

    fitter.run()

    rms_error, _ = fitter.getDataRMSAndMaximumProjectionErrorForGroup(branch_name)
    # print(branch_name, rms_error)

    return fitter


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
