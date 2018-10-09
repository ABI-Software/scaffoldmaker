"""
Generates a 3-D elliptical tube mesh from a central line
with variable numbers of elements around, along and
through wall, with variable major and minor axis length and
wall thickness.
"""

from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh import *

class MeshType_3d_ellipticaltube1(object):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Elliptical Tube 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements around' : 8,
            'Number of elements along' : 1,
            'Number of elements through wall' : 1,
            'Major axis length of inner ellipse': 2.0,
            'Minor axis length of inner ellipse': 1.0,
            'Wall thickness': 0.25,
            'Use cross derivatives' : False,
            'Use linear through wall' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements along',
            'Number of elements through wall',
            'Major axis length of inner ellipse',
            'Minor axis length of inner ellipse',
            'Wall thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements around'] < 2) :
            options['Number of elements around'] = 2
        for key in [
            'Major axis length of inner ellipse',
            'Minor axis length of inner ellipse',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountAlong = options['Number of elements along']
        elementsCountThroughWall = options['Number of elements through wall']
        a = options['Major axis length of inner ellipse']
        b = options['Minor axis length of inner ellipse']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])

        nextNodeIdentifier, nextElementIdentifier = generatetubemesh(region, elementsCountAlong, elementsCountAround, elementsCountThroughWall, 
            a, b, wallThickness, useCrossDerivatives, useCubicHermiteThroughWall)

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
