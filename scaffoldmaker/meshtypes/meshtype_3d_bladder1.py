"""
Generates a 3-D bladder mesh.
"""

from __future__ import division
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_3d_bladder1(Scaffold_base):
    '''
    3-D bladder scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Bladder 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        return {
        }

    def getOrderedOptionNames():
        return []

    @classmethod
    def checkOptions(cls, options):
        pass

    @classmethod
    def generateMesh(cls, region, options):
        return None
