"""
Generates a 3-D heart model including ventricles, base and atria.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.annotation.annotationgroup import mergeAnnotationGroups
from scaffoldmaker.meshtypes.meshtype_3d_heartatria2 import MeshType_3d_heartatria2
from scaffoldmaker.meshtypes.meshtype_3d_heartventriclesbase2 import MeshType_3d_heartventriclesbase2
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_heart2(Scaffold_base):
    '''
    Generates a 3-D heart model including ventricles, base and atria.
    '''

    @staticmethod
    def getName():
        return '3D Heart 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = MeshType_3d_heartventriclesbase2.getDefaultOptions(parameterSetName)
        optionsAtria = MeshType_3d_heartatria2.getDefaultOptions(parameterSetName)
        options.update(optionsAtria)
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = MeshType_3d_heartventriclesbase2.getOrderedOptionNames()
        optionNamesAtria = MeshType_3d_heartatria2.getOrderedOptionNames()
        # insert numbers of elements in atria in initial group
        for optionName in [
            'Number of elements up atria',
            'Number of elements around atrial septum']:
            optionNames.insert(5, optionName)
            optionNamesAtria.remove(optionName)
        # remove dependent or repeated options in atria2
        optionNamesAtria.remove('LV outlet outer diameter')
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
            'Refine number of elements through RV wall',
            'Refine number of elements through atrial wall']:
            optionNames.remove(optionName)
            optionNames.append(optionName)
        return optionNames

    @staticmethod
    def checkOptions(options):
        dependentChanges = MeshType_3d_heartventriclesbase2.checkOptions(options) \
            or MeshType_3d_heartatria2.checkOptions(options)
        # only works with particular numbers of elements around
        options['Number of elements around atrial septum'] = 2
        # set dependent outer diameter used in atria2
        options['LV outlet outer diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        # set dependent outer diameter used in atria2
        options['LV outlet outer diameter'] = options['LV outlet inner diameter'] + 2.0*options['LV outlet wall thickness']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        # generate heartventriclesbase2 model and put atria2 on it
        ventriclesAnnotationGroups = MeshType_3d_heartventriclesbase2.generateBaseMesh(region, options)
        atriaAnnotationGroups = MeshType_3d_heartatria2.generateBaseMesh(region, options)
        annotationGroups = mergeAnnotationGroups(ventriclesAnnotationGroups, atriaAnnotationGroups)

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
        MeshType_3d_heartventriclesbase2.refineMesh(meshrefinement, options)
        MeshType_3d_heartatria2.refineMesh(meshrefinement, options)
