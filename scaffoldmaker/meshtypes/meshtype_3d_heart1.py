"""
Generates a 3-D heart model including ventricles, base and atria.
"""

from __future__ import division
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.meshtype_3d_heartatria1 import MeshType_3d_heartatria1
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase1 import MeshType_3d_heartventriclesbase1
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.meshrefinement import MeshRefinement

class MeshType_3d_heart1(object):
    '''
    Generates a 3-D heart model including ventricles, base and atria.
    '''

    @staticmethod
    def getName():
        return '3D Heart 1'

    @staticmethod
    def getDefaultOptions():
        options = MeshType_3d_heartatria1.getDefaultOptions()
        # ventricles overrides some atria options, so update from it:
        optionsVentricles = MeshType_3d_heartventriclesbase1.getDefaultOptions()
        options.update(optionsVentricles)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventriclesbase1.getOrderedOptionNames()
        optionNamesAtria = MeshType_3d_heartatria1.getOrderedOptionNames()
        # insert numbers of elements in atria in initial group
        for optionName in [
            'Number of elements around atrial free wall',
            'Number of elements around atrial septum',
            'Number of elements up atria',
            'Number of elements inlet']:
            optionNames.insert(6, optionName)
            optionNamesAtria.remove(optionName)
        # remove dependent or repeated options in atria1
        optionNamesAtria.remove('Aorta outer plus diameter')
        for optionName in optionNames:
            if optionName in optionNamesAtria:
                optionNamesAtria.remove(optionName)
        # add remaining atria options
        optionNames += optionNamesAtria
        # want refinement options last
        for optionName in [
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through LV wall',
            'Refine number of elements through wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        MeshType_3d_heartventriclesbase1.checkOptions(options)
        MeshType_3d_heartatria1.checkOptions(options)
        # only works with particular numbers of elements around
        #options['Number of elements around atrial septum'] = 2
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        # set dependent outer diameter used in atria2
        options['Aorta outer plus diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        # generate heartventriclesbase2 model and put atria2 on it
        annotationGroups = MeshType_3d_heartventriclesbase1.generateBaseMesh(region, options)
        annotationGroups += MeshType_3d_heartatria1.generateBaseMesh(region, options)

        fm.endChange()
        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        MeshType_3d_heartventriclesbase1.refineMesh(meshrefinement, options)
        MeshType_3d_heartatria1.refineMesh(meshrefinement, options)

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()