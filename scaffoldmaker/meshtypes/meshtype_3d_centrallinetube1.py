"""
Generates a 3-D tubular mesh from a central line
with variable numbers of elements around, along and
through wall, with variable major and minor axis length and
wall thickness.
"""

from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh import *

class MeshType_3d_centrallinetube1(Scaffold_base):
    '''
    Generates a 3-D tubular mesh with variable numbers
    of elements around, along the central line, and through wall.
    The ellipsoidal tube is created from a central line and
    lateral axes data
    '''
    @staticmethod
    def getName():
        return '3D Central Line Tube 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 8,
            'Number of elements along' : 6,
            'Number of elements through wall' : 1,
            'Tube type' : 1,
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
            'Tube type',
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
        if (options['Tube type'] < 1 or options['Tube type'] > 4 ) :
            options['Tube type'] = 1

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
        tubeType = options['Tube type']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])

        if tubeType == 1:
            # Straight tube
            cx = [[1.0, 3.0, 0.0], [ 2.0, 0.0, 4.0 ] ]
            cd1 = [[ 1.0, -3.0, 4.0 ], [ 1.0, -3.0, 4.0 ]]
            cd2 = [ [ 0.0, 0.2, 0.0 ], [ 0.0, 0.2, 0.0 ]]
            cd3 = [ [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.5 ]]

            # thickness in cd2 and cd3 directions and derivatives (rate of change)
            t2 = [ 0.1, 0.5 ]
            t2d = [ 0.0, 0.0 ]
            t3 = [ 0.5, 0.5 ]
            t3d = [ 0.0, 0.0 ]

        elif tubeType == 2:
            # Curved tube 1
            cx = [ [ 0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [ 2.0, 0.0, 0.0 ] ]
            cd1 = [ [ 1.0, 1.0, 0.0 ], [ 1.0, 0.0, 0.0 ], [1.0, -1.0, 0.0] ]
            cd2 = [ [ 0.0, 0.1, 0.0 ], [ 0.0, 0.1, 0.0 ], [ 0.0, 0.1, 0.0 ] ]
            cd3 = [ [ 0.0, 0.0, 0.2 ], [ 0.0, 0.0, 0.2 ], [ 0.0, 0.0, 0.2] ]

            # thickness in cd2 and cd3 directions and derivatives (rate of change)
            t2 = [ 0.1, 0.1, 0.1 ]
            t2d = [ 0.0, 0.0, 0.0 ]
            t3 = [ 0.1, 0.1, 0.1 ]
            t3d = [ 0.0, 0.0, 0.0 ]

        elif tubeType == 3:
            # Curved tube 2
            cx = [ [ 0.0, 0.0, 1.0], [1.5, 1.0, 0.0], [ 3.0, -1.0, 0.0 ], [ 5.0, 1.5, 1.0]]
            cd1 = [ [ 4.0, 0.0, 0.0 ], [ 2.0, 0.0, 0.0 ], [3.0, 0.0, 0.0], [ 3.0, 0.0, 0.0 ]]
            cd2 = [ [ 0.0, 0.2, 0.0 ], [ 0.0, 0.2, 0.0 ], [ 0.0, 0.2, 0.0 ], [ 0.0, 0.2, 0.0] ]
            cd3 = [ [ 0.0, 0.0, 0.2 ], [ 0.0, 0.0, 0.2 ], [ 0.0, 0.0, 0.2], [ 0.0, 0.0, 0.2 ]]

            # thickness in cd2 and cd3 directions and derivatives (rate of change)
            t2 = [ 0.1, 0.1, 0.1, 0.1 ]
            t2d = [ 0.0, 0.0, 0.0, 0.0]
            t3 = [ 0.1, 0.1, 0.1, 0.1]
            t3d = [ 0.0, 0.0, 0.0, 0.0]
        
        elif tubeType == 4:
            # Colon
            cx = [ [ 0.0, 0.0, 0.0], [0.0, 10.0, 0.0], [5.0, 9.0, 0.0], [ 10.0, 10.0, 0.0 ], [ 10.0, -2.0, 0.0], [ 7.0, -4.0, 0.0] ]
            cd1 = [ [ 0.0, 10.0, 0.0 ], [ 5.0, 5.0, 0.0 ], [5.0, 0.0, 0.0], [ 5.0, -5.0, 0.0 ], [ -3.0, -5.0, 0.0 ], [ -3.0, 0.0, 0.0 ]]
            cd2 = [ [ 0.0, 0.5, 0.0 ], [ 0.0, 0.5, 0.0 ], [ 0.0, 0.5, 0.0 ], [ 0.0, 0.5, 0.0], [ 0.0, 0.5, 0.0 ], [ 0.0, 0.5, 0.0 ]]
            cd3 = [ [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.5], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.5], [ 0.0, 0.0, 0.5]]

            # thickness in cd2 and cd3 directions and derivatives (rate of change)
            t2 = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ]
            t2d = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
            t3 = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
            t3d = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        nextNodeIdentifier, nextElementIdentifier = generatetubemesh(region, elementsCountAlong, elementsCountAround, elementsCountThroughWall, 
            cx, cd1, cd2, cd3, t2, t2d, t3, t3d, useCrossDerivatives, useCubicHermiteThroughWall)

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
