"""
Generates a 3-D vertebra model including body, spinous process, transverse processes,
articular processes, vertebral arch and pedicle.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase1 import MeshType_3d_heartventriclesbase1
from scaffoldmaker.utils.eft_utils import *
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from opencmiss.zinc.element import Element


class MeshType_3d_vertebra1(object):
    """
    Generates a 3-D vertebral model.
    """

    @staticmethod
    def getName():
        return '3D Vertebra 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Type of vertebra' : "Thoracic",  # Cervical, Thoracic, Lumbar
            'Number of elements for vertebral body' : 8,
            'Unit scale': 1.0,
            'Vertebral body semi-axis a length factor' : 1,
            'Vertebral body semi-axis b length factor' : 1,
            'Vertebral body semi-axis c length factor' : 1,
            'Vertebral body posterior curvature factor' : 1,
            'Refine': False,
            'Refine number of elements surface': 4,
            'Refine number of elements through body': 1,
            'Use cross derivatives': False,
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Type of vertebra',
            'Number of elements around atrial free wall',
            'Unit scale',
            'Vertebral body semi-axis a length factor',
            'Vertebral body semi-axis b length factor',
            'Vertebral body semi-axis c length factor',
            'Vertebral body posterior curvature factor',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through body',
            #,'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        """
        :return:  True if dependent options changed, otherwise False. This
        happens where two or more options must change together to be valid.
        """
        return None

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """

        vertebralType = options['Type of vertebra']
        elementsCountForVertebralBody = options['Number of elements for vertebral body']
        unitScale = options['Unit scale']
        bodyAxisALengthFactor = options['Vertebral body semi-axis a length factor']
        bodyAxisBLengthFactor = options['Vertebral body semi-axis b length factor']
        bodyAxisCLengthFactor = options['Vertebral body semi-axis c length factor']
        bodyPosteriorCurvatureFactor = options['Vertebral body posterior curvature factor']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        ##############
        # Create nodes
        ##############

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        # LA/RA inlet elements are linear through the wall, hence their nodes do not have D_DS3 parameters
        nodetemplateLinearS3 = nodes.createNodetemplate()
        nodetemplateLinearS3.defineField(coordinates)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateLinearS3.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)





