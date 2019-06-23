"""
Generates a 3-D heart atria model, suitable for attachment to the
3-D Heart Ventricles with Base 1.
"""

from __future__ import division
import copy
import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getCircleProjectionAxes, getEllipseAngleFromVector, getEllipseArcLength, getEllipseRadiansToX, updateEllipseAngleByArcLength, createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import zinc_utils
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartatria1(Scaffold_base):
    '''
    3-D heart atria model, suitable for attachment to the 3-D Heart Ventricles with Base 1.
    '''

    lpvOstiumDefaultScaffoldPackages = {
        'LPV Human 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 2,
                'Number of elements across common' : 2,
                'Number of elements around vessel' : 8,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.2,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.16,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.11,
                'Vessel wall thickness' : 0.009,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 30.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'LPV Pig 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 1,
                'Number of elements across common' : 2,
                'Number of elements around vessel' : 8,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.23,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.16,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.16,
                'Vessel wall thickness' : 0.011,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 30.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } )
        }

    rpvOstiumDefaultScaffoldPackages = {
        'RPV Human 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 2,
                'Number of elements across common' : 2,
                'Number of elements around vessel' : 8,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.2,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.2,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.12,
                'Vessel wall thickness' : 0.009,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 30.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'RPV Pig 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 1,
                'Number of elements across common' : 2,
                'Number of elements around vessel' : 8,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.24,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.16,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.17,
                'Vessel wall thickness' : 0.011,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 30.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } )
        }

    @staticmethod
    def getName():
        return '3D Heart Atria 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Unit Human 1',
            'Unit Mouse 1',
            'Unit Pig 1',
            'Unit Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Pig' in parameterSetName:
            lpvOstium = cls.lpvOstiumDefaultScaffoldPackages['LPV Pig 1']
            rpvOstium = cls.rpvOstiumDefaultScaffoldPackages['RPV Pig 1']
        else:
            lpvOstium = cls.lpvOstiumDefaultScaffoldPackages['LPV Human 1']
            rpvOstium = cls.rpvOstiumDefaultScaffoldPackages['RPV Human 1']
        options = {}
        options['Number of elements around atrial septum'] = 3
        options['Number of elements around left atrium free wall'] = 8
        options['Number of elements around right atrium free wall'] = 6
        options['Number of elements over atria'] = 8
        options['Unit scale'] = 1.0
        options['Atria base inner major axis length'] = 0.55
        options['Atria base inner minor axis length'] = 0.45
        options['Atria major axis rotation degrees'] = 40.0
        options['Atria outer height'] = 0.45
        options['Atrial septum height'] = 0.3
        options['Atrial septum length'] = 0.3
        options['Atrial septum thickness'] = 0.075
        options['Coronary sinus height'] = 0.07
        options['Fossa ovalis height'] = 0.1
        options['Fossa ovalis length'] = 0.15
        options['Fossa ovalis thickness'] = 0.025
        options['Fossa ovalis midpoint height'] = 0.17
        options['Left atrium venous free wall thickness'] = 0.025
        options['Right atrium venous free wall thickness'] = 0.015
        options['Crista terminalis thickness'] = 0.035
        options['Atrial base wall thickness'] = 0.07
        options['Atrial base slope degrees'] = 30.0
        options['Aorta outer plus diameter'] = 0.35
        options['Atrial base front incline degrees'] = 15.0
        options['Atrial base back incline degrees'] = 20.0
        options['Atrial base side incline degrees'] = 20.0
        options['Left atrial appendage left'] = 0.9
        options['Left atrium venous anterior over'] = 0.7
        options['Left atrium venous midpoint posterior left'] = 0.5
        options['Right atrium venous midpoint over'] = 0.41
        options['Right atrium venous right'] = 0.4
        options['Left pulmonary vein ostium'] = copy.deepcopy(lpvOstium)
        options['Left pulmonary vein ostium angle degrees'] = 20.0
        options['Left pulmonary vein ostium position left'] = 0.62
        options['Left pulmonary vein ostium position over'] = 0.5
        options['Right pulmonary vein ostium'] = copy.deepcopy(rpvOstium)
        options['Right pulmonary vein ostium angle degrees'] = 80.0
        options['Right pulmonary vein ostium position left'] = 0.16
        options['Right pulmonary vein ostium position over'] = 0.4
        options['Inferior vena cava inlet position over'] = 0.18
        options['Inferior vena cava inlet position right'] = 0.2
        options['Inferior vena cava inlet angle left degrees'] = -15.0
        options['Inferior vena cava inlet angle over degrees'] = -30.0
        options['Inferior vena cava inlet derivative factor'] = 1.0
        options['Inferior vena cava inlet length'] = 0.1
        options['Inferior vena cava inlet inner diameter'] = 0.22
        options['Inferior vena cava inlet wall thickness'] = 0.015
        options['Superior vena cava inlet position over'] = 0.65
        options['Superior vena cava inlet position right'] = 0.2
        options['Superior vena cava inlet angle left degrees'] = -20.0
        options['Superior vena cava inlet angle over degrees'] = 10.0
        options['Superior vena cava inlet derivative factor'] = 1.0
        options['Superior vena cava inlet length'] = 0.1
        options['Superior vena cava inlet inner diameter'] = 0.18
        options['Superior vena cava inlet wall thickness'] = 0.015
        options['Refine'] = False
        options['Refine number of elements surface'] = 4
        options['Refine number of elements through wall'] = 1
        options['Use cross derivatives'] = False

        if 'Human' in parameterSetName:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 80.0
        elif 'Mouse' in parameterSetName:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 5.0
            #options['Left pulmonary vein angle up degrees'] = 30.0
            #options['Left pulmonary vein inner diameter'] = 0.16
            #options['Left pulmonary vein wall thickness'] = 0.011
            #options['Right pulmonary vein angle up degrees'] = 10.0
            #options['Right pulmonary vein inner diameter'] = 0.15
            #options['Right pulmonary vein wall thickness'] = 0.011
            options['Superior vena cava inlet inner diameter'] = 0.17
            options['Superior vena cava inlet wall thickness'] = 0.012
        elif 'Pig' in parameterSetName:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 80.0
            options['Inferior vena cava inlet position over'] = 0.18
            options['Inferior vena cava inlet position right'] = 0.18
            options['Inferior vena cava inlet angle left degrees'] = 30.0
            options['Inferior vena cava inlet angle over degrees'] = -30.0
            options['Inferior vena cava inlet derivative factor'] = 0.5
            options['Superior vena cava inlet position over'] = 0.6
            options['Superior vena cava inlet angle left degrees'] = -15.0
            options['Superior vena cava inlet angle over degrees'] = -10.0
        elif 'Rat' in parameterSetName:
            if 'Unit' not in parameterSetName:
                options['Unit scale'] = 12.0
            #options['Number of left pulmonary veins'] = 1
            #options['Left pulmonary vein angle up degrees'] = 30.0
            #options['Left pulmonary vein inner diameter'] = 0.16
            #options['Left pulmonary vein wall thickness'] = 0.011
            #options['Right pulmonary vein angle up degrees'] = 10.0
            #options['Right pulmonary vein inner diameter'] = 0.15
            #options['Right pulmonary vein wall thickness'] = 0.011
            options['Superior vena cava inlet inner diameter'] = 0.17
            options['Superior vena cava inlet wall thickness'] = 0.012
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around atrial septum',
            'Number of elements around left atrium free wall',
            'Number of elements around right atrium free wall',
            'Number of elements over atria',
            'Unit scale',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atria outer height',
            'Atrial septum height',
            'Atrial septum length',
            'Atrial septum thickness',
            'Coronary sinus height',
            'Fossa ovalis midpoint height',
            'Fossa ovalis height',
            'Fossa ovalis length',
            'Fossa ovalis thickness',
            'Left atrium venous free wall thickness',
            'Right atrium venous free wall thickness',
            'Crista terminalis thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Aorta outer plus diameter',
            'Atrial base front incline degrees',
            'Atrial base back incline degrees',
            'Atrial base side incline degrees',
            'Left atrial appendage left',
            'Left atrium venous anterior over',
            'Left atrium venous midpoint posterior left',
            'Right atrium venous midpoint over',
            'Right atrium venous right',
            'Left pulmonary vein ostium',
            'Left pulmonary vein ostium angle degrees',
            'Left pulmonary vein ostium position left',
            'Left pulmonary vein ostium position over',
            'Right pulmonary vein ostium',
            'Right pulmonary vein ostium angle degrees',
            'Right pulmonary vein ostium position left',
            'Right pulmonary vein ostium position over',
            'Inferior vena cava inlet position over',
            'Inferior vena cava inlet position right',
            'Inferior vena cava inlet angle left degrees',
            'Inferior vena cava inlet angle over degrees',
            'Inferior vena cava inlet derivative factor',
            'Inferior vena cava inlet length',
            'Inferior vena cava inlet inner diameter',
            'Inferior vena cava inlet wall thickness',
            'Superior vena cava inlet position over',
            'Superior vena cava inlet position right',
            'Superior vena cava inlet angle left degrees',
            'Superior vena cava inlet angle over degrees',
            'Superior vena cava inlet derivative factor',
            'Superior vena cava inlet length',
            'Superior vena cava inlet inner diameter',
            'Superior vena cava inlet wall thickness',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall',
            #,'Use cross derivatives'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName in [ 'Left pulmonary vein ostium', 'Right pulmonary vein ostium' ]:
            return [ MeshType_3d_ostium1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Left pulmonary vein ostium':
            return list(cls.lpvOstiumDefaultScaffoldPackages.keys())
        if optionName == 'Right pulmonary vein ostium':
            return list(cls.rpvOstiumDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Left pulmonary vein ostium':
            if not parameterSetName:
                parameterSetName = list(cls.lpvOstiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.lpvOstiumDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Right pulmonary vein ostium':
            if not parameterSetName:
                parameterSetName = list(cls.rpvOstiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.rpvOstiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @staticmethod
    def checkOptions(options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        if options['Number of elements around atrial septum'] < 2:
            options['Number of elements around atrial septum'] = 2
        for key in [
            'Number of elements around left atrium free wall',
            'Number of elements around right atrium free wall']:
            if options[key] < 6:
                options[key] = 6
            elif options[key] > 10:
                options[key] = 10
        for key in ['Number of elements over atria']:
            if options[key] < 6:
                options[key] = 6
            elif options[key] > 6:
                options[key] = 8
        for key in [
            'Unit scale',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria outer height',
            'Atrial septum height',
            'Atrial septum length',
            'Atrial septum thickness',
            'Coronary sinus height',
            'Fossa ovalis height',
            'Fossa ovalis length',
            'Fossa ovalis thickness',
            'Fossa ovalis midpoint height',
            'Left atrium venous free wall thickness',
            'Right atrium venous free wall thickness',
            'Crista terminalis thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Inferior vena cava inlet derivative factor',
            'Inferior vena cava inlet length',
            'Inferior vena cava inlet inner diameter',
            'Inferior vena cava inlet wall thickness',
            'Superior vena cava inlet derivative factor',
            'Superior vena cava inlet length',
            'Superior vena cava inlet inner diameter',
            'Superior vena cava inlet wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'Left atrial appendage left',
            'Left atrium venous anterior over',
            'Left atrium venous midpoint posterior left',
            'Right atrium venous midpoint over',
            'Right atrium venous right',
            'Left pulmonary vein ostium position left',
            'Left pulmonary vein ostium position over',
            'Right pulmonary vein ostium position left',
            'Right pulmonary vein ostium position over',
            'Inferior vena cava inlet position over',
            'Inferior vena cava inlet position right',
            'Superior vena cava inlet position over',
            'Superior vena cava inlet position right']:
            if options[key] < 0.001:
                options[key] = 0.001
            elif options[key] > 0.999:
                options[key] = 0.999
        if options['Aorta outer plus diameter'] < options['Atrial septum thickness']:
            options['Aorta outer plus diameter'] = options['Atrial septum thickness']
        for key in [
            'Atria major axis rotation degrees']:
            if options[key] < -75.0:
                options[key] = -75.0
            elif options[key] > 75.0:
                options[key] = 75.0
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountAroundLeftAtrium = elementsCountAroundLeftAtriumFreeWall + elementsCountAroundAtrialSeptum
        elementsCountAroundRightAtriumFreeWall = options['Number of elements around right atrium free wall']
        elementsCountOverAtria = options['Number of elements over atria']
        unitScale = options['Unit scale']

        aBaseInnerMajorMag = unitScale*0.5*options['Atria base inner major axis length']
        aBaseInnerMinorMag = unitScale*0.5*options['Atria base inner minor axis length']
        aMajorAxisRadians = math.radians(options['Atria major axis rotation degrees'])
        aOuterHeight = unitScale*options['Atria outer height']
        aortaOuterPlusRadius = unitScale*0.5*options['Aorta outer plus diameter']
        aBaseFrontInclineRadians = math.radians(options['Atrial base front incline degrees'])
        aBaseSideInclineRadians = math.radians(options['Atrial base side incline degrees'])
        aBaseBackInclineRadians = math.radians(options['Atrial base back incline degrees'])
        aSeptumHeight = unitScale*options['Atrial septum height']
        aSeptumLength = unitScale*options['Atrial septum length']
        aSeptumThickness = unitScale*options['Atrial septum thickness']
        coronarySinusHeight = unitScale*options['Coronary sinus height']
        foMidpointZ = unitScale*options['Fossa ovalis midpoint height']
        foMagZ = unitScale*0.5*options['Fossa ovalis height']
        foMagY = unitScale*0.5*options['Fossa ovalis length']
        foThickness = unitScale*options['Fossa ovalis thickness']
        laVenousFreeWallThickness = unitScale*options['Left atrium venous free wall thickness']
        raVenousFreeWallThickness = unitScale*options['Right atrium venous free wall thickness']
        cristaTerminalisThickness = unitScale*options['Crista terminalis thickness']
        aBaseWallThickness = unitScale*options['Atrial base wall thickness']
        aBaseSlopeRadians = math.radians(options['Atrial base slope degrees'])
        laaLeft = options['Left atrial appendage left']
        laVenousLimitAnterior = options['Left atrium venous anterior over']
        laVenousMidpointPosteriorLeft = options['Left atrium venous midpoint posterior left']
        raVenousRight = options['Right atrium venous right']
        raVenousMidpointOver = options['Right atrium venous midpoint over']
        lpvOstium = options['Left pulmonary vein ostium']
        lpvOstiumAngleRadians = math.radians(options['Left pulmonary vein ostium angle degrees'])
        lpvOstiumPositionLeft = options['Left pulmonary vein ostium position left']
        lpvOstiumPositionOver = options['Left pulmonary vein ostium position over']
        rpvOstium = options['Right pulmonary vein ostium']
        rpvOstiumAngleRadians = math.radians(options['Right pulmonary vein ostium angle degrees'])
        rpvOstiumPositionLeft = options['Right pulmonary vein ostium position left']
        rpvOstiumPositionOver = options['Right pulmonary vein ostium position over']
        ivcPositionOver = options['Inferior vena cava inlet position over']
        ivcPositionRight = options['Inferior vena cava inlet position right']
        ivcAngleLeftRadians = math.radians(options['Inferior vena cava inlet angle left degrees'])
        ivcAngleOverRadians = math.radians(options['Inferior vena cava inlet angle over degrees'])
        ivcDerivativeFactor = options['Inferior vena cava inlet derivative factor']
        ivcLength = unitScale*options['Inferior vena cava inlet length']
        ivcInnerRadius = unitScale*0.5*options['Inferior vena cava inlet inner diameter']
        ivcWallThickness = unitScale*options['Inferior vena cava inlet wall thickness']
        svcPositionOver = options['Superior vena cava inlet position over']
        svcPositionRight = options['Superior vena cava inlet position right']
        svcAngleLeftRadians = math.radians(options['Superior vena cava inlet angle left degrees'])
        svcAngleOverRadians = math.radians(options['Superior vena cava inlet angle over degrees'])
        svcDerivativeFactor = options['Superior vena cava inlet derivative factor']
        svcLength = unitScale*options['Superior vena cava inlet length']
        svcInnerRadius = unitScale*0.5*options['Superior vena cava inlet inner diameter']
        svcWallThickness = unitScale*options['Superior vena cava inlet wall thickness']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = zinc_utils.getOrCreateCoordinateField(fm)
        cache = fm.createFieldcache()

        laGroup = AnnotationGroup(region, 'left atrium', FMANumber = 7097, lyphID = 'Lyph ID unknown')
        raGroup = AnnotationGroup(region, 'right atrium', FMANumber = 7096, lyphID = 'Lyph ID unknown')
        aSeptumGroup = AnnotationGroup(region, 'interatrial septum', FMANumber = 7108, lyphID = 'Lyph ID unknown')
        fossaGroup = AnnotationGroup(region, 'fossa ovalis', FMANumber = 9246, lyphID = 'Lyph ID unknown')
        lipvGroup = AnnotationGroup(region, 'left inferior pulmonary vein', FMANumber = 49913, lyphID = 'Lyph ID unknown')
        lspvGroup = AnnotationGroup(region, 'left superior pulmonary vein', FMANumber = 49916, lyphID = 'Lyph ID unknown')
        ripvGroup = AnnotationGroup(region, 'right inferior pulmonary vein', FMANumber = 49911, lyphID = 'Lyph ID unknown')
        rspvGroup = AnnotationGroup(region, 'right superior pulmonary vein', FMANumber = 49914, lyphID = 'Lyph ID unknown')
        ivcInletGroup = AnnotationGroup(region, 'inferior vena cava inlet', FMANumber = 10951, lyphID = 'Lyph ID unknown')
        svcInletGroup = AnnotationGroup(region, 'superior vena cava inlet', FMANumber = 4720, lyphID = 'Lyph ID unknown')
        annotationGroups = [ laGroup, raGroup, aSeptumGroup, fossaGroup, lipvGroup, lspvGroup, ripvGroup, rspvGroup, ivcInletGroup, svcInletGroup ]
        # av boundary nodes are put in left and right fibrous ring groups only so they can be found by heart1
        lFibrousRingGroup = AnnotationGroup(region, 'left fibrous ring', FMANumber = 77124, lyphID = 'Lyph ID unknown')
        rFibrousRingGroup = AnnotationGroup(region, 'right fibrous ring', FMANumber = 77125, lyphID = 'Lyph ID unknown')

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

        nodeIdentifier = max(1, zinc_utils.getMaximumNodeIdentifier(nodes) + 1)

        mesh = fm.findMeshByDimension(3)

        elementIdentifier = max(1, zinc_utils.getMaximumElementIdentifier(mesh) + 1)

        laMeshGroup = laGroup.getMeshGroup(mesh)
        raMeshGroup = raGroup.getMeshGroup(mesh)
        aSeptumMeshGroup = aSeptumGroup.getMeshGroup(mesh)
        fossaMeshGroup = fossaGroup.getMeshGroup(mesh)
        lipvMeshGroup = lipvGroup.getMeshGroup(mesh)
        lspvMeshGroup = lspvGroup.getMeshGroup(mesh)
        ripvMeshGroup = ripvGroup.getMeshGroup(mesh)
        rspvMeshGroup = rspvGroup.getMeshGroup(mesh)
        ivcInletMeshGroup = ivcInletGroup.getMeshGroup(mesh)
        svcInletMeshGroup = svcInletGroup.getMeshGroup(mesh)

        # get elements count over atria, around left and right free wall base, and around ostia
        # note elementsCountOverAtriaCoronarySinus is assumed to be 1
        elementsCountOverAtriaCoronarySinus, \
        elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior, \
        elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior \
            = getOverAtriaElementsCounts(elementsCountOverAtria)
        elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftArialAppendageBase, elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV \
            = getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall)
        elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightArialAppendageBase, elementsCountAroundRightAtriumAorta \
            = getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall)
        elementsCountAroundLpvOstium = 10 if (elementsCountOverAtria == 8) else 6
        #elementsCountAroundLpvOstium = elementsCountOverLeftAtriumVenous + elementsCountOverLeftAtriumNonVenousPosterior - elementsCountOverAtriaCoronarySinus \
        #    + 2*elementsCountAroundLeftAtriumLPV + 1
        #if (elementsCountAroundLpvOstium % 2) == 1:
        #    elementsCountAroundLpvOstium += 1
        elementsCountOverSideLeftAtriumLPV = elementsCountAroundLpvOstium - elementsCountAroundLeftAtriumLPV - elementsCountOverLeftAtriumVenous
        elementsCountAroundRpvOstium = 2*(elementsCountOverLeftAtriumVenous + elementsCountAroundLeftAtriumRPV)

        # GRC fudge factors:
        aOuterSeptumHeight = 1.2*aSeptumHeight
        iaGrooveDerivative = 0.25*aSeptumThickness

        aBaseSlopeHeight = aBaseWallThickness*math.sin(aBaseSlopeRadians)
        aBaseSlopeLength = aBaseWallThickness*math.cos(aBaseSlopeRadians)
        aBaseOuterMajorMag = aBaseInnerMajorMag + aBaseSlopeLength
        aBaseOuterMinorMag = aBaseInnerMinorMag + aBaseSlopeLength

        elementsCountAroundTrackSurface = 20  # must be even, twice number of elements along
        elementsCountAcrossTrackSurface = 10
        labx, labd1, labd2, labd3, rabx, rabd1, rabd2, rabd3, ltBaseOuterx, ltBaseOuterd1, ltBaseOuterd2, aSeptumBaseCentre, laCentre, laSeptumRadians, = \
            getAtriumBasePoints(elementsCountAroundAtrialSeptum, elementsCountAroundLeftAtriumFreeWall, elementsCountAroundRightAtriumFreeWall,
                aBaseInnerMajorMag, aBaseInnerMinorMag, aMajorAxisRadians,
                aBaseWallThickness, aBaseSlopeHeight, aBaseSlopeLength, aSeptumLength, aSeptumThickness,
                aortaOuterPlusRadius, aBaseFrontInclineRadians, aBaseSideInclineRadians, aBaseBackInclineRadians,
                laaLeft, laVenousMidpointPosteriorLeft, raVenousRight, elementsCountAroundTrackSurface)

        laTrackSurface = getAtriumTrackSurface(elementsCountAroundTrackSurface, elementsCountAcrossTrackSurface,
            ltBaseOuterx, ltBaseOuterd1, ltBaseOuterd2, aSeptumBaseCentre, aOuterHeight, aOuterSeptumHeight, iaGrooveDerivative)
        raTrackSurface = laTrackSurface.createMirrorX()

        # need to create pulmonary vein ostia early because other derivatives are smoothed to fit them

        # create left pulmonary vein ostium
        lpvOstiumSettings = copy.deepcopy(lpvOstium.getScaffoldSettings())
        lpvOstiumSettings['Unit scale'] *= unitScale
        lpvOstiumSettings['Ostium wall thickness'] = laVenousFreeWallThickness
        lpvOstiumSettings['Outlet'] = False
        lpvOstiumSettings['Use linear through ostium wall'] = False
        lpvCount = lpvOstiumSettings['Number of vessels']
        if lpvCount == 1:
            lpvOstiumSettings['Number of elements around vessel'] = elementsCountAroundLpvOstium
        else:  # 2 vessels
            elementsCountAcrossLpvOstium = 2
            lpvOstiumSettings['Number of elements across common'] = elementsCountAcrossLpvOstium
            lpvOstiumSettings['Number of elements around vessel'] = elementsCountAroundLpvOstium//2 + elementsCountAcrossLpvOstium
        lpvOstiumPosition = laTrackSurface.createPositionProportion(lpvOstiumPositionOver, lpvOstiumPositionLeft)
        # get absolute direction on surface corresponding to chosen angle
        cx, cd1, cd2 = laTrackSurface.evaluateCoordinates(lpvOstiumPosition, derivatives = True)
        td1, td2, td3 = calculate_surface_axes(cd1, cd2, [ 0.0, 1.0, 0.0 ])
        zAngleRadians = math.atan2(td1[0], -td2[0])
        #print('zAngleRadians',zAngleRadians)
        cosAngle = math.cos(zAngleRadians + lpvOstiumAngleRadians)
        sinAngle = math.sin(zAngleRadians + lpvOstiumAngleRadians)
        lpvOstiumDirection = [ (cosAngle*-td2[c] + sinAngle*td1[c]) for c in range(3) ]
        vesselMeshGroups = [ [ laMeshGroup, lipvMeshGroup ], [ laMeshGroup, lspvMeshGroup ] ] if (lpvCount == 2) else [ [ laMeshGroup, lipvMeshGroup, lspvMeshGroup ] ]
        nodeIdentifier, elementIdentifier, (lpvox, lpvod1, lpvod2, lpvod3, lpvoNodeId, lpvoPositions) = \
            generateOstiumMesh(region, lpvOstiumSettings, laTrackSurface, lpvOstiumPosition, lpvOstiumDirection, nodeIdentifier, elementIdentifier, vesselMeshGroups)

        # create right pulmonary vein ostium
        rpvOstiumSettings = copy.deepcopy(rpvOstium.getScaffoldSettings())
        rpvOstiumSettings['Unit scale'] *= unitScale
        rpvOstiumSettings['Ostium wall thickness'] = laVenousFreeWallThickness
        rpvOstiumSettings['Outlet'] = False
        rpvOstiumSettings['Use linear through ostium wall'] = False
        rpvCount = rpvOstiumSettings['Number of vessels']
        if rpvCount == 1:
            rpvOstiumSettings['Number of elements around vessel'] = elementsCountAroundRpvOstium
        else:  # 2 vessels
            elementsCountAcrossRpvOstium = 2
            rpvOstiumSettings['Number of elements across common'] = elementsCountAcrossRpvOstium
            rpvOstiumSettings['Number of elements around vessel'] = elementsCountAroundRpvOstium//2 + elementsCountAcrossRpvOstium
        rpvOstiumPosition = laTrackSurface.createPositionProportion(rpvOstiumPositionOver, rpvOstiumPositionLeft)
        # get absolute direction on surface corresponding to chosen angle
        cx, cd1, cd2 = laTrackSurface.evaluateCoordinates(rpvOstiumPosition, derivatives = True)
        td1, td2, td3 = calculate_surface_axes(cd1, cd2, [ 0.0, 1.0, 0.0 ])
        zAngleRadians = math.atan2(td1[0], -td2[0])
        #print('zAngleRadians',zAngleRadians)
        cosAngle = math.cos(zAngleRadians + rpvOstiumAngleRadians)
        sinAngle = math.sin(zAngleRadians + rpvOstiumAngleRadians)
        rpvOstiumDirection = [ (cosAngle*-td2[c] + sinAngle*td1[c]) for c in range(3) ]
        vesselMeshGroups = [ [ laMeshGroup, ripvMeshGroup ], [ laMeshGroup, rspvMeshGroup ] ] if (lpvCount == 2) else [ [ laMeshGroup, ripvMeshGroup, rspvMeshGroup ] ]
        nodeIdentifier, elementIdentifier, (rpvox, rpvod1, rpvod2, rpvod3, rpvoNodeId, rpvoPositions) = \
            generateOstiumMesh(region, rpvOstiumSettings, laTrackSurface, rpvOstiumPosition, rpvOstiumDirection, nodeIdentifier, elementIdentifier, vesselMeshGroups)

        # get points over interatrial septum on exterior groove
        agn1Mid = elementsCountOverRightAtriumNonVenousAnterior + elementsCountOverRightAtriumVenous//2
        # 1. go up coronary sinus height on posterior
        nx = laTrackSurface.nx [:laTrackSurface.elementsCount1 + 1]
        nd = laTrackSurface.nd1[:laTrackSurface.elementsCount1 + 1]
        csx, csd2, e, xi = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd, coronarySinusHeight)
        lagcsProportion = (e + xi)/laTrackSurface.elementsCount1
        agLength = sum(interp.getCubicHermiteArcLength(nx[e], nd[e], nx[e + 1], nd[e + 1]) for e in range(laTrackSurface.elementsCount1))
        # arbitrarily set midpoint derivative to reduce element spacing to fit nearby RPV ostium
        agMidpointDerivative = 0.75/elementsCountOverAtria  # GRC fudge factor tweak was 1.0
        agx  = [ labx [1][0] ]
        agd1 = [ labd1[1][0] ]
        agd2 = [ labd2[1][0] ]  # rescale later
        # add lengths over groove from anterior to posterior, intersecting laVenousLimitAnterior, raVenousMidpointOver, lagcsProportion
        ragProportionLengths1 = interp.sampleCubicElementLengths(1.0 - laVenousLimitAnterior, elementsCountOverLeftAtriumNonVenousAnterior)
        ragProportionLengths2 = interp.sampleCubicElementLengths(laVenousLimitAnterior - raVenousMidpointOver, elementsCountOverLeftAtriumVenous//2, \
            startDerivative = ragProportionLengths1[-1], endDerivative = agMidpointDerivative)
        ragProportionLengths3 = interp.sampleCubicElementLengths(raVenousMidpointOver - lagcsProportion, \
            elementsCountOverLeftAtriumVenous//2 + elementsCountOverLeftAtriumNonVenousPosterior - elementsCountOverAtriaCoronarySinus, \
            startDerivative = ragProportionLengths2[-1], endDerivative = lagcsProportion)
        ragProportionLengths = ragProportionLengths1 + ragProportionLengths2 + ragProportionLengths3 + [ lagcsProportion ]
        # get d1 magnitudes over crest at posterior, middle and anterior/aorta
        d1a = vector.magnitude(labd1[1][0])
        d1p = vector.magnitude(labd1[1][elementsCountAroundLeftAtriumFreeWall])
        d1m = 0.25*d1p  # GRC fudge factor - not necessarily used to get trackSurface!
        ragProportion = 0.0
        ragProportions = [ 0.0 ]
        for e in range(elementsCountOverAtria - 1):
            ragProportion += ragProportionLengths[e]
            ragProportions.append(ragProportion)
            trackPosition = raTrackSurface.createPositionProportion(ragProportion, 0.0)
            x, d2 = raTrackSurface.evaluateCoordinates(trackPosition, derivatives = True)[0:2]
            agx.append(x)
            if ragProportion < 0.5:
                d1s = d1a
                d1f = d1m
                xid1 = 2.0*ragProportion
            else:
                d1s = d1m
                d1f = d1p
                xid1 = 2.0*(ragProportion - 0.5)
            f1, _, f3, _ = interp.getCubicHermiteBasis(xid1)
            agd1.append([ -(f1*d1s + f3*d1f), 0.0, 0.0 ])
            agd2.append(vector.setMagnitude(d2, agLength*0.5*(ragProportionLengths[e] + ragProportionLengths[e + 1])))
        ragProportions.append(1.0)
        laVenousLimitPosterior = ragProportionLengths[-1] + ragProportionLengths[-2]
        # Get heights of elements on aorta up interatrial groove, for building venous limit curves
        aoHeight1 = ragProportionLengths[0]*agLength
        # fix coronary sinus d2 magnitude
        agd2[-1] = vector.setMagnitude(agd2[-1], coronarySinusHeight)
        # smooth d2 up to coronary sinus
        agd2 = interp.smoothCubicHermiteDerivativesLine(agx, agd2, fixAllDirections = True, fixEndDerivative = True)  # , magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        # GRC , magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        # reverse derivative on posterior coronary sinus
        agd1[-1] = [ -d for d in agd1[-1] ]
        agd2[-1] = [ -d for d in agd2[-1] ]
        # add posterior crux point
        agx .append(labx [1][elementsCountAroundLeftAtriumFreeWall])
        agd1.append(labd1[1][elementsCountAroundLeftAtriumFreeWall])
        agd2.append(vector.setMagnitude(labd2[1][elementsCountAroundLeftAtriumFreeWall], coronarySinusHeight))
        agd3 = [ None ]*(elementsCountOverAtria + 1)  # set later, using adjacent points
        # first and last d3 are known:
        agd3[0] = labd3[1][0]
        agd3[-1] = labd3[1][elementsCountAroundLeftAtriumFreeWall]
        # copy derivatives to labd2[1], rabd2[1]
        rabd2[1][elementsCountAroundRightAtriumFreeWall] = labd2[1][0] = agd2[0]
        rabd2[1][0] = labd2[1][elementsCountAroundLeftAtriumFreeWall] = agd2[-1]

        # start getting points on interatrial septum next to coronary sinus, then septum "arch"
        # need these to get fossa angles
        halffoThickness = 0.5*foThickness
        halfaSeptumThickness = 0.5*aSeptumThickness
        xia = 0.35  # GRC fudge factor
        coronarySinusHeightAnterior = (1.0 - xia)*coronarySinusHeight + xia*foMidpointZ
        aSeptumPosteriorY = aSeptumBaseCentre[1] - 0.5*aSeptumLength
        aSeptumAnteriorY = aSeptumBaseCentre[1] + 0.5*aSeptumLength
        x1 = [ 0.0, aSeptumPosteriorY, coronarySinusHeight ]
        d1 = [ 0.0, aSeptumLength, 0.0 ]
        x2 = [ 0.0, aSeptumAnteriorY, coronarySinusHeightAnterior ]
        d2 = interp.interpolateHermiteLagrangeDerivative(x1, d1, x2, 1.0)  # GRC was d1
        asx, asd1 = interp.sampleCubicHermiteCurves([ x1, x2 ], [ d1, d2 ], elementsCountAroundAtrialSeptum, arcLengthDerivatives = True)[0:2]

        # get fossa ovalis points at centre and around
        elementsCountAroundFossa = elementsCountOverAtria + elementsCountAroundAtrialSeptum - 2
        fossaPerimeterLength = getApproximateEllipsePerimeter(foMagY, foMagZ)
        estElementSizeAroundFossa = fossaPerimeterLength/elementsCountAroundFossa
        fossaInnerDerivativeRatio = 1.0  # 4.0/3.0  # GRC fudge factor
        fossaOuterDerivativeRatio = 2.0 - fossaInnerDerivativeRatio
        foMidpointY = aSeptumBaseCentre[1]
        fossaRadiansAround = []
        fox = []
        fod1 = []
        fod2 = []
        for nf in range(elementsCountAroundFossa):
            if nf <= elementsCountAroundAtrialSeptum:
                y = asx[nf][1]
                z = asx[nf][2]
            else:
                na = nf - elementsCountAroundAtrialSeptum + 1
                y = agx[na][1]
                z = agx[na][2]
            radiansAround = getEllipseAngleFromVector(foMagY, foMagZ, y - foMidpointY, z - foMidpointZ)
            #print('fossa y', y, 'z', z, 'angle', math.degrees(radiansAround))
            fossaRadiansAround.append(radiansAround)
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            fox .append([ 0.0, foMidpointY + foMagY*cosRadiansAround, foMidpointZ + foMagZ*sinRadiansAround ])
            fod1.append(vector.setMagnitude([ 0.0, -foMagY*sinRadiansAround, foMagZ*cosRadiansAround ], estElementSizeAroundFossa))
            fod2.append([ 0.0, -fossaOuterDerivativeRatio*foMagY*cosRadiansAround, -fossaOuterDerivativeRatio*foMagZ*sinRadiansAround ])
        fod1 = interp.smoothCubicHermiteDerivativesLoop(fox, fod1, fixAllDirections = True)  # GRC was , magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        fossad3const = [ foThickness, 0.0, 0.0 ]
        foCentrex = []
        foCentred1 = []
        foCentred2 = []
        foCentred3 = []
        fox  = [ fox , copy.deepcopy(fox ) ]
        fod1 = [ fod1, copy.deepcopy(fod1) ]
        fod2 = [ fod2, copy.deepcopy(fod2) ]
        fod3 = [ [ fossad3const ]*elementsCountAroundFossa ]*2
        for n3 in range(2):
            fossaX = (-0.5 if (n3 == 0) else + 0.5)*foThickness
            foCentrex.append([ fossaX, foMidpointY, foMidpointZ ])
            foCentred1.append([ 0.0, fossaInnerDerivativeRatio*foMagY, 0.0 ])
            foCentred2.append([ 0.0, 0.0, fossaInnerDerivativeRatio*foMagZ ])
            foCentred3.append(fossad3const)
            for nf in range(elementsCountAroundFossa):
                fox[n3][nf][0] = fossaX

        # complete getting points on interatrial septum next to coronary sinus, then septum "arch"
        archMagY = 0.5*aSeptumLength
        archMagZ = aSeptumHeight - foMidpointZ
        halfArchEllipsePerimeter = 0.5*getApproximateEllipsePerimeter(archMagY, archMagZ)
        archLength = (2.0*foMidpointZ - coronarySinusHeight - coronarySinusHeightAnterior) + halfArchEllipsePerimeter
        elementsCountOverArch = elementsCountOverAtria - 2
        estArchElementLength = archLength/elementsCountOverArch
        asd2 = [ [ 0.0, 0.0, -estArchElementLength ] ]
        for ns in range(1, elementsCountAroundAtrialSeptum):
            nf = ns
            # make d2 in normal d3, d1
            d2 = interp.interpolateLagrangeHermiteDerivative([ asx[ns][0] - halffoThickness, asx[ns][1], asx[ns][2] ], fox[0][nf], fod2[0][nf], 0.0)
            d2 = vector.setMagnitude([ 0.0, -asd1[ns][2], asd1[ns][1] ], vector.magnitude(d2))
            d2 = interp.smoothCubicHermiteDerivativesLine([ [ asx[ns][0] - halffoThickness, asx[ns][1], asx[ns][2] ], fox[0][nf] ], [ d2, fod2[0][nf] ],
                fixStartDirection = True, fixEndDerivative = True)[0]
            asd2.append(d2)
        asd2.append([ 0.0, 0.0, estArchElementLength ])
        # add points over arch:
        archMagY = 0.5*aSeptumLength
        archMagZ = aSeptumHeight - foMidpointZ
        for na in range(1, elementsCountOverArch):
            nf = elementsCountAroundAtrialSeptum + na
            ng = na + 1
            #print('na', na, 'ng', ng, 'agx', agx[ng])
            if agx[ng][2] > foMidpointZ:
                radiansAround = getEllipseAngleFromVector(archMagY, archMagZ, agx[ng][1] - foMidpointY, agx[ng][2] - foMidpointZ)
                #print('arch y', agx[ng][1], 'z', agx[ng][2], 'angle', math.degrees(radiansAround))
                if radiansAround > 0.0:
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [ 0.0, foMidpointY + archMagY*cosRadiansAround, foMidpointZ + archMagZ*sinRadiansAround ]
                    d2 = vector.setMagnitude([ 0.0, -archMagY*sinRadiansAround, archMagZ*cosRadiansAround ], estArchElementLength)
            elif agx[ng][1] > foMidpointY:
                xi = (aSeptumAnteriorY - foMidpointY)/(agx[ng][1] - foMidpointY)
                x = [ 0.0, aSeptumAnteriorY, (1.0 - xi)*foMidpointZ + xi*agx[ng][2] ]
                d2 = [ 0.0, 0.0, estArchElementLength ]
            else:  # agx[ng][1] > foMidpointY
                xi = (aSeptumPosteriorY - foMidpointY)/(agx[ng][1] - foMidpointY)
                x = [ 0.0, aSeptumPosteriorY, (1.0 - xi)*foMidpointZ + xi*agx[ng][2] ]
                d2 = [ 0.0, 0.0, -estArchElementLength ]
            # make d1 normal to d2, d3
            d1 = interp.interpolateHermiteLagrangeDerivative([ 0.0, fox[0][nf][1], fox[0][nf][2] ], [ -d for d in fod2[0][nf] ], x, 1.0)
            d1 = vector.setMagnitude([ 0.0, d2[2], -d2[1] ], vector.magnitude(d1))
            d1 = interp.smoothCubicHermiteDerivativesLine([ [ 0.0, fox[0][nf][1], fox[0][nf][2] ], x ], [ [ -d for d in fod2[0][nf] ], d1 ],
                fixStartDerivative = True, fixEndDirection = True)[1]
            asx .append(x )
            asd1.append(d1)
            asd2.append(d2)
        asd2 = interp.smoothCubicHermiteDerivativesLoop(asx, asd2, fixAllDirections = True, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        # reverse first d2:
        asd2[0] = [ -d for d in asd2[0] ]
        # copy and displace for la, ra
        asx  = [ asx  , copy.deepcopy(asx ) ]
        asd1 = [ asd1 , copy.deepcopy(asd1) ]
        asd2 = [ asd2 , copy.deepcopy(asd2) ]
        septumd3const = [ aSeptumThickness, 0.0, 0.0 ]
        septumd3minus = [ -aSeptumThickness, 0.0, 0.0 ]
        asd3 = [ [ septumd3const ]*elementsCountAroundFossa, [ septumd3minus ]*elementsCountAroundFossa ]
        for ns in range(elementsCountAroundFossa):
            asx [0][ns][0] -= halfaSeptumThickness
            asx [1][ns][0] += halfaSeptumThickness
            asd1[1][ns] = [ -d for d in asd1[1][ns] ]
            asd2[1][ns][0] = -asd2[1][ns][0]
        # fix up derivative 2 on atrial septum base
        for ns in range(elementsCountAroundAtrialSeptum + 1):
            nl = -elementsCountAroundAtrialSeptum + ns
            labd2[0][nl] = interp.interpolateLagrangeHermiteDerivative(labx[0][nl], [ asx[0][ns][0], asx[0][ns][1], asx[0][ns][2] ], asd2[0][ns], 0.0)
            nr = -ns
            rabd2[0][nr] = [ -labd2[0][nl][0], labd2[0][nl][1], labd2[0][nl][2] ]
        # fix derivative 3 on interatrial groove
        for na in range(elementsCountOverArch + 1):
            ng = 1 + na
            # this is a kludge, but looks alright:
            agd3[ng] = vector.setMagnitude(vector.crossproduct3(agd1[ng], agd2[ng]), raVenousFreeWallThickness)

        # get points on external coronary sinus, around entire free wall of both atria
        # not all of these will become nodes, but they are always used to set base derivatives
        # left atrium
        lacsx  = [ agx [1] ]
        lacsd1 = [ agd1[1] ]
        lacsd2 = [ agd2[1] ]
        lacsd3 = [ agd3[1] ]
        # also get points on left atrium venous limit posterior row above cs
        lavpHeightAboveCs = interp.getCubicHermiteArcLength(agx[-3], agd2[-3], agx[-2], [ -d for d in agd2[-2] ])
        lavpx  = []
        lavpd1 = []
        lavpProportions = []
        lan1Mid = elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftArialAppendageBase
        lan1MidVenous = lan1Mid + elementsCountAroundLeftAtriumLPV
        lacsProportions = [ [ 1.0, 0.0 ] ]
        laApexx = laTrackSurface.nx[-1]
        for n1 in range(1, elementsCountAroundLeftAtriumFreeWall):
            # find position on laTrackSurface corresponding to base outer node
            # avoid point at end of laTrackSurface where derivative 1 is zero
            onEnd = True
            for c in range(3):
                if math.fabs(labx[1][n1][c] - laApexx[c]) > 0.0001:
                    onEnd = False
                    break
            if onEnd:
                startPosition = laTrackSurface.createPositionProportion(0.5, 0.9999)
            else:
                if n1 <= lan1Mid:
                    startProportion1 = 1.0
                    startProportion2 = laaLeft*n1/lan1Mid
                else:
                    startProportion1 = 0.0
                    startProportion2 = (elementsCountAroundLeftAtriumFreeWall - n1)/(elementsCountAroundLeftAtriumFreeWall - lan1Mid)
                startPosition = laTrackSurface.findNearestPosition(labx[1][n1], laTrackSurface.createPositionProportion(startProportion1, startProportion2))
            onAorta = n1 == 1
            direction = [ 0.0, 0.0, 1.0 ] if onAorta else vector.normalise(labd2[1][n1])
            trackDistance1 = aoHeight1 if onAorta else coronarySinusHeight
            position = laTrackSurface.trackVector(startPosition, direction, trackDistance1)
            lacsProportions.append([ (position.e1 + position.xi1)/laTrackSurface.elementsCount1, (position.e2 + position.xi2)/laTrackSurface.elementsCount2 ])
            x, d1, d2 = laTrackSurface.evaluateCoordinates(position, derivatives = True)
            ax1, ax2, ax3 = calculate_surface_axes(d1, d2, direction)
            lacsx .append(x)
            lacsd1.append(vector.setMagnitude(ax2, -vector.magnitude(labd1[1][n1])))
            lacsd2.append(vector.setMagnitude(ax1, trackDistance1))
            lacsd3.append(vector.setMagnitude(ax3, laVenousFreeWallThickness))
            # fix d2 on outer base
            labd2[1][n1] = vector.setMagnitude(labd2[1][n1], trackDistance1)
            if n1 >= lan1Mid:
                position = laTrackSurface.trackVector(position, direction, lavpHeightAboveCs)
                lavpProportions.append([ (position.e1 + position.xi1)/laTrackSurface.elementsCount1, (position.e2 + position.xi2)/laTrackSurface.elementsCount2 ])
                x, d1, d2 = laTrackSurface.evaluateCoordinates(position, derivatives = True)
                ax1, ax2, ax3 = calculate_surface_axes(d1, d2, direction)
                lavpx .append(x)
                lavpd1.append(vector.setMagnitude(ax2, vector.magnitude(labd1[1][n1])))  # approximate; smoothed later
        # add end points and smooth d1
        lacsx .append(agx [-2])
        lacsd1.append(agd1[-2])
        lacsd2.append(agd2[-2])
        lacsd3.append(agd3[-2])
        lacsProportions.append([ 0.0, 0.0 ])
        lacsd1 = interp.smoothCubicHermiteDerivativesLine(lacsx, lacsd1, fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True)
        # get inner points
        lacsx  = [ [agx [0]], lacsx  ]
        lacsd1 = [ [agd1[0]], lacsd1 ]
        lacsd2 = [ [agd2[0]], lacsd2 ]
        lacsd3 = [ [agd3[0]], lacsd3 ]
        for n1 in range(1, elementsCountAroundLeftAtriumFreeWall):
            x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(lacsx[1], lacsd1[1], lacsd2[1], n1, -laVenousFreeWallThickness)
            # do same upwards to get proper value of d2
            nx  = [ labx [1][n1], lacsx [1][n1] ]
            nd1 = [ labd1[1][n1], lacsd1[1][n1] ]
            nd2 = [ labd2[1][n1], lacsd2[1][n1] ]
            _, d2, _, _ = interp.projectHermiteCurvesThroughWall(nx, nd2, [ [ -d for d in d1 ] for d1 in nd1 ], 1, -laVenousFreeWallThickness)
            lacsx [0].append(x)
            lacsd1[0].append(d1)
            lacsd2[0].append(d2)
            lacsd3[0].append(d3)
            # fix d2 on inner base
            labd2[0][n1] = interp.interpolateLagrangeHermiteDerivative(labx[0][n1], x, d2, 0.0)
        lacsx [0].append(agx [-1])
        lacsd1[0].append(agd1[-1])
        lacsd2[0].append(agd2[-1])
        lacsd3[0].append(agd3[-1])
        # right atrium
        racsx  = [ agx [-2] ]
        racsd1 = [ agd1[-2] ]
        racsd2 = [ agd2[-2] ]
        racsd3 = [ agd3[-2] ]
        racsProportions = [ [ 1.0, 0.0 ] ]
        ran1Ctp = elementsCountAroundRightAtriumPosteriorVenous
        #ran1raa1 = elementsCountAroundRightAtriumPosteriorVenous + 1
        ran1Aorta = elementsCountAroundRightAtriumFreeWall - elementsCountAroundRightAtriumAorta
        raApexx = raTrackSurface.nx[-1]
        for n1 in range(1, elementsCountAroundRightAtriumFreeWall):
            # find position on raTrackSurface corresponding to base outer node
            # avoid point at end of raTrackSurface where derivative 1 is zero
            onEnd = True
            for c in range(3):
                if math.fabs(rabx[1][n1][c] - raApexx[c]) > 0.0001:
                    onEnd = False
                    break
            if onEnd:
                startPosition = raTrackSurface.createPositionProportion(0.5, 0.9999)
            else:
                startProportion1 = (elementsCountAroundRightAtriumFreeWall - n1)/elementsCountAroundRightAtriumFreeWall
                startPosition = raTrackSurface.findNearestPosition(rabx[1][n1], raTrackSurface.createPositionProportion(startProportion1, 0.5))
            onAorta = n1 == (elementsCountAroundRightAtriumFreeWall - 1)
            direction = [ 0.0, 0.0, 1.0 ] if onAorta else vector.normalise(rabd2[1][n1])
            trackDistance1 = aoHeight1 if onAorta else coronarySinusHeight
            position = raTrackSurface.trackVector(startPosition, direction, trackDistance1)
            racsProportions.append([ (position.e1 + position.xi1)/raTrackSurface.elementsCount1, (position.e2 + position.xi2)/raTrackSurface.elementsCount2 ])
            x, d1, d2 = raTrackSurface.evaluateCoordinates(position, derivatives = True)
            ax1, ax2, ax3 = calculate_surface_axes(d1, d2, direction)
            racsx .append(x)
            racsd1.append(vector.setMagnitude(ax2, -vector.magnitude(rabd1[1][n1])))
            racsd2.append(vector.setMagnitude(ax1, trackDistance1))
            racsd3.append(vector.setMagnitude(ax3, raVenousFreeWallThickness))
            # fix d2 on outer base
            rabd2[1][n1] = vector.setMagnitude(rabd2[1][n1], trackDistance1)
        # add end points and smooth d1
        racsx .append(agx [1])
        racsd1.append(agd1[1])
        racsd2.append(agd2[1])
        racsd3.append(agd3[1])
        lacsProportions.append([ 0.0, 0.0 ])
        racsd1 = interp.smoothCubicHermiteDerivativesLine(racsx, racsd1, fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True)
        # get inner points
        racsx  = [ [agx [-1]], racsx  ]
        racsd1 = [ [agd1[-1]], racsd1 ]
        racsd2 = [ [agd2[-1]], racsd2 ]
        racsd3 = [ [agd3[-1]], racsd3 ]
        for n1 in range(1, elementsCountAroundRightAtriumFreeWall):
            x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(racsx[1], racsd1[1], racsd2[1], n1, -raVenousFreeWallThickness)
            # do same upwards to get proper value of d2
            nx  = [ rabx [1][n1], racsx [1][n1] ]
            nd1 = [ rabd1[1][n1], racsd1[1][n1] ]
            nd2 = [ rabd2[1][n1], racsd2[1][n1] ]
            _, d2, _, _ = interp.projectHermiteCurvesThroughWall(nx, nd2, [ [ -d for d in d1 ] for d1 in nd1 ], 1, -raVenousFreeWallThickness)
            racsx [0].append(x)
            racsd1[0].append(d1)
            racsd2[0].append(d2)
            racsd3[0].append(d3)
            # fix d2 on inner base
            rabd2[0][n1] = interp.interpolateLagrangeHermiteDerivative(rabx[0][n1], x, d2, 0.0)
        racsx [0].append(agx [0])
        racsd1[0].append(agd1[0])
        racsd2[0].append(agd2[0])
        racsd3[0].append(agd3[0])

        # get points on left atrium over appendage, from aorta to laa end on coronary sinus
        laoax, laoad1, laoad2, laoad3, laoaProportions = laTrackSurface.createHermiteCurvePoints(
            lacsProportions[1][0], lacsProportions[1][1], lacsProportions[lan1Mid][0], lacsProportions[lan1Mid][1],
            elementsCount = elementsCountAroundLeftAtriumRPV - 1 + elementsCountOverSideLeftAtriumLPV,
            derivativeStart = [ (1.0*lacsd1[1][1][c] + lacsd2[1][1][c]) for c in range(3) ],
            derivativeEnd = [ -1.0*d for d in lacsd2[1][lan1Mid] ])
        # get inner points
        laoax  = [ [ None ], laoax  ]
        laoad1 = [ [ None ], laoad1 ]
        laoad2 = [ [ None ], laoad2 ]
        laoad3 = [ [ None ], laoad3 ]
        for n in range(1, len(laoax[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(laoax[1], laoad1[1], laoad2[1], n, -laVenousFreeWallThickness)
            laoax [0].append(x)
            laoad1[0].append(d1)
            laoad2[0].append(d2)
            laoad3[0].append(d3)
            laoad3[1][n] = d3
        # substitute known start and end coordinates
        for n3 in range(2):
            laoax [n3][ 0] = lacsx [n3][1]
            laoad1[n3][ 0] = lacsd1[n3][1]
            laoad2[n3][ 0] = lacsd2[n3][1]
            laoad3[n3][ 0] = lacsd3[n3][1]
            laoax [n3][-1] = lacsx [n3][lan1Mid]
            laoad1[n3][-1] = lacsd1[n3][lan1Mid]
            laoad2[n3][-1] = lacsd2[n3][lan1Mid]
            laoad3[n3][-1] = lacsd3[n3][lan1Mid]
        # smooth d2 to fit adjacent LPV derivative 2
        for n1 in range(2, elementsCountOverSideLeftAtriumLPV + 1):
            n1lpv = -elementsCountOverLeftAtriumVenous//2 - n1
            for n3 in range(2):
                nx  = [ laoax [n3][n1], lpvox [n3][n1lpv] ]
                nd1 = [ laoad2[n3][n1], [ -d for d in lpvod2[n3][n1lpv] ] ]
                laoad2[n3][n1] = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixAllDirections = True, fixEndDerivative = True)[0]

        # get points from left atrium venous anterior on interatrial septum to nearest point on LPV ostium
        # need index of nearest point around LPV ostium
        n1lpv = -1 - elementsCountOverLeftAtriumVenous//2
        position = lpvoPositions[n1lpv]
        lpvoProportion1 = (position.e1 + position.xi1)/laTrackSurface.elementsCount1
        lpvoProportion2 = (position.e2 + position.xi2)/laTrackSurface.elementsCount2
        agn1va = elementsCountOverLeftAtriumNonVenousAnterior
        lavbx, lavbd1, lavbd2, lavbd3, lavbProportions = laTrackSurface.createHermiteCurvePoints(laVenousLimitAnterior, 0.0, lpvoProportion1, lpvoProportion2, \
            elementsCount = elementsCountAroundLeftAtriumRPV + 2,
            derivativeStart = [ d for d in agd1[agn1va] ],
            derivativeEnd = [ -0.5*d for d in lpvod2[1][n1lpv] ])  # fudge factor to control closeness to LPV
        # smooth to fit actual end derivatives
        lavbd1[0] = agd1[agn1va]
        lavbd1[-1] = [ -d for d in lpvod2[1][n1lpv] ]
        lavbd1 = interp.smoothCubicHermiteDerivativesLine(lavbx, lavbd1, fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True)
        # get inner points
        lavbx  = [ [ None ], lavbx  ]
        lavbd1 = [ [ None ], lavbd1 ]
        lavbd2 = [ [ None ], lavbd2 ]
        lavbd3 = [ [ None ], lavbd3 ]
        for n1 in range(1, len(lavbx[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(lavbx[1], lavbd1[1], lavbd2[1], n1, -laVenousFreeWallThickness)
            lavbx [0].append(x)
            lavbd1[0].append(d1)
            lavbd2[0].append(d2)
            lavbd3[0].append(d3)
            lavbd3[1][n1] = d3
        # substitute known start and end coordinates
        asn1va = elementsCountAroundAtrialSeptum - 1 + elementsCountOverLeftAtriumNonVenousAnterior
        lavbx [0][0] = asx [0][asn1va]
        lavbd1[0][0] = asd1[0][asn1va]
        lavbd2[0][0] = asd2[0][asn1va]
        lavbd3[0][0] = asd3[0][asn1va]
        lavbx [1][0] = agx [0][agn1va]
        lavbd1[1][0] = agd1[0][agn1va]
        lavbd2[1][0] = agd2[0][agn1va]
        lavbd3[1][0] = agd3[0][agn1va]
        for n3 in range(2):
            lavbx [n3][-1] = lpvox [n3][n1lpv]
            lavbd1[n3][-1] = lpvod1[n3][n1lpv]
            lavbd2[n3][-1] = lpvod2[n3][n1lpv]
            lavbd3[n3][-1] = lpvod3[n3][n1lpv]


        # get points on left atrium non-venous posterior row above coronary sinus
        # start with lavp* created in cs loop
        # start from ag and reverse
        lavpx  = [ agx [-3] ] + list(reversed(lavpx ))
        lavpd1 = [ agd1[-3] ] + list(reversed(lavpd1))
        lavpProportions = [ [ laVenousLimitPosterior, 0.0 ] ] + list(reversed(lavpProportions))
        # smooth d1 left of ag:
        lavpd1 = interp.smoothCubicHermiteDerivativesLine(lavpx, lavpd1, fixAllDirections = True, fixStartDerivative = True, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)

        # get points on left atrium venous posterior, from venous posterior on interatrial septum to nearest point on LPV ostium
        # need index of nearest point around LPV ostium
        n1lpv = 1
        position = lpvoPositions[n1lpv]
        lpvoProportion1 = (position.e1 + position.xi1)/laTrackSurface.elementsCount1
        lpvoProportion2 = (position.e2 + position.xi2)/laTrackSurface.elementsCount2
        agn1vp = elementsCountOverLeftAtriumNonVenousAnterior + elementsCountOverLeftAtriumVenous
        lavqx, lavqd1, lavqd2, lavqd3, lavqProportions = laTrackSurface.createHermiteCurvePoints(lavpProportions[1][0], lavpProportions[1][1], lpvoProportion1, lpvoProportion2, \
            elementsCount = elementsCountAroundLeftAtriumRPV + 1,
            derivativeStart = lavpd1[1],
            derivativeEnd = [ -1.0*d for d in lpvod2[1][n1lpv] ])  # fudge factor to control closeness to LPV
        lavqx .insert(0, agx [agn1vp])
        lavqd1.insert(0, agd1[agn1vp])
        lavqd2.insert(0, agd2[agn1vp])
        lavqd3.insert(0, agd3[agn1vp])
        # smooth to fit actual end derivatives
        #lavqd1[0] = agd1[agn1vp]
        #lavqd1[-1] = [ -d for d in lpvod2[1][n1lpv] ]
        #lavqd1 = interp.smoothCubicHermiteDerivativesLine(lavqx, lavqd1, fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True)
        # get inner points
        lavqx  = [ [ None ], lavqx  ]
        lavqd1 = [ [ None ], lavqd1 ]
        lavqd2 = [ [ None ], lavqd2 ]
        lavqd3 = [ [ None ], lavqd3 ]
        for n1 in range(1, len(lavqx[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(lavqx[1], lavqd1[1], lavqd2[1], n1, -laVenousFreeWallThickness)
            lavqx [0].append(x)
            lavqd1[0].append(d1)
            lavqd2[0].append(d2)
            lavqd3[0].append(d3)
            lavqd3[1][n1] = d3
        # substitute known start and end coordinates
        asn1vp = -1
        lavqx [0][0] = asx [0][asn1vp]
        lavqd1[0][0] = asd1[0][asn1vp]
        lavqd2[0][0] = asd2[0][asn1vp]
        lavqd3[0][0] = asd3[0][asn1vp]
        lavqx [1][0] = asx [0][agn1vp]
        lavqd1[1][0] = asd1[0][agn1vp]
        lavqd2[1][0] = asd2[0][agn1vp]
        lavqd3[1][0] = asd3[0][agn1vp]
        for n3 in range(2):
            lavqx [n3][-1] = lpvox [n3][n1lpv]
            lavqd1[n3][-1] = lpvod1[n3][n1lpv]
            lavqd2[n3][-1] = lpvod2[n3][n1lpv]
            lavqd3[n3][-1] = lpvod3[n3][n1lpv]

        # get left atrium venous midpoint line from lavb to lavq
        # use points between middle/nearest points on left and right PVs
        n1 = elementsCountOverLeftAtriumVenous//2
        n1lpv = n1 - elementsCountOverLeftAtriumVenous//2 - 1
        n1rpv = -elementsCountOverLeftAtriumVenous//2 - elementsCountAroundLeftAtriumRPV - n1
        positionr = rpvoPositions[n1rpv]
        rpvoProportion1 = (positionr.e1 + positionr.xi1)/laTrackSurface.elementsCount1
        rpvoProportion2 = (positionr.e2 + positionr.xi2)/laTrackSurface.elementsCount2
        positionl = lpvoPositions[n1lpv]
        lpvoProportion1 = (positionl.e1 + positionl.xi1)/laTrackSurface.elementsCount1
        lpvoProportion2 = (positionl.e2 + positionl.xi2)/laTrackSurface.elementsCount2
        lamRowCount = 1
        mpx, mpd1, mpd2, _, nProportions = laTrackSurface.createHermiteCurvePoints(rpvoProportion1, rpvoProportion2, lpvoProportion1, lpvoProportion2,
            elementsCount = lamRowCount + 1, derivativeStart = rpvod2[1][n1rpv], derivativeEnd = [ -d for d in lpvod2[1][n1lpv] ])
        # scale mid derivative 2 to be mean of d1 in LPV, RPV
        d2mag = 0.75*vector.magnitude(lpvod1[1][n1lpv]) + 0.75*vector.magnitude(rpvod1[1][n1rpv])  # GRC fudge factor to enlarge midpoint derivative 1
        for n1 in range(1, lamRowCount + 1):
            mpd2[n1] = vector.setMagnitude(mpd2[n1], d2mag)
        # la midline left
        lamlx, lamld2, lamld1, lamld3 = laTrackSurface.createHermiteCurvePoints(lavbProportions[-2][0], lavbProportions[-2][1], nProportions[-2][0], nProportions[-2][1], \
            elementsCount = elementsCountOverLeftAtriumVenous//2, derivativeStart = lavbd2[1][-2], derivativeEnd = mpd2[-2])[0:4]
        _lamlx, _lamld2, _lamld1, _lamld3 = laTrackSurface.createHermiteCurvePoints(nProportions[-2][0], nProportions[-2][1], lavqProportions[-2][0], lavqProportions[-2][1],
            elementsCount = elementsCountOverLeftAtriumVenous//2, derivativeStart = mpd2[-2], derivativeEnd = lavqd2[1][-2])[0:4]
        # get inner points
        lamlx  = [ [ None ], lamlx  + _lamlx [1:] ]
        lamld1 = [ [ None ], [ [ -d for d in d1 ] for d1 in (lamld1 + _lamld1[1:]) ] ]
        lamld2 = [ [ None ], lamld2  + _lamld2 [1:] ]
        lamld3 = [ [ None ], lamld3  + _lamld3 [1:] ]
        for n2 in range(1, len(lamlx[1])):
            x, d2, d1, d3 = interp.projectHermiteCurvesThroughWall(lamlx[1], lamld2[1], [ [ -d for d in d1 ] for d1 in lamld1[1] ], n2, -laVenousFreeWallThickness)
            lamlx [0].append(x)
            lamld1[0].append([ -d for d in d1 ])
            lamld2[0].append(d2)
            lamld3[0].append(d3)
            lamld3[1][n2] = d3
        # substitute known start and end coordinates
        for n3 in range(2):
            lamlx [n3][ 0] = lavbx [n3][-2]
            lamld1[n3][ 0] = lavbd1[n3][-2]
            lamld2[n3][ 0] = lavbd2[n3][-2]
            lamld3[n3][ 0] = lavbd3[n3][-2]
            lamlx [n3][-1] = lavqx [n3][-2]
            lamld1[n3][-1] = lavqd1[n3][-2]
            lamld2[n3][-1] = lavqd2[n3][-2]
            lamld3[n3][-1] = lavqd3[n3][-2]
        # la midline right
        if lamRowCount == 1:
            lamrx  = [ copy.copy(lamlx [0]), copy.copy(lamlx [1]) ]
            lamrd1 = [ copy.copy(lamld1[0]), copy.copy(lamld1[1]) ]
            lamrd2 = [ copy.copy(lamld2[0]), copy.copy(lamld2[1]) ]
            lamrd3 = [ copy.copy(lamld3[0]), copy.copy(lamld3[1]) ]
        else:
            lamrx, lamrd2, lamrd1, lamrd3 = laTrackSurface.createHermiteCurvePoints(lavbProportions[-3][0], lavbProportions[-3][1], nProportions[1][0], nProportions[1][1], 
                elementsCount = elementsCountOverLeftAtriumVenous//2, derivativeStart = [ (lavbd1[1][-2][c] + lavbd2[1][-2][c]) for c in range(3) ], derivativeEnd = mpd2[1])[0:4]
            _lamrx, _lamrd2, _lamrd1, _lamrd3 = laTrackSurface.createHermiteCurvePoints(nProportions[1][0], nProportions[1][1], lavqProportions[-3][0], lavqProportions[-3][1],
                elementsCount = elementsCountOverLeftAtriumVenous//2, derivativeStart = mpd2[1], derivativeEnd = [ (-lavqd1[1][-3][c] + lavqd2[1][-3][c]) for c in range(3) ])[0:4]
            # get inner points
            lamrx  = [ [ None ], lamrx  + _lamrx [1:] ]
            lamrd1 = [ [ None ], [ [ -d for d in d1 ] for d1 in (lamrd1 + _lamrd1[1:]) ] ]
            lamrd2 = [ [ None ], lamrd2  + _lamrd2 [1:] ]
            lamrd3 = [ [ None ], lamrd3  + _lamrd3 [1:] ]
            for n2 in range(1, len(lamrx[1])):
                x, d2, d1, d3 = interp.projectHermiteCurvesThroughWall(lamrx[1], lamrd2[1], [ [ -d for d in d1 ] for d1 in lamrd1[1] ], n2, -laVenousFreeWallThickness)
                lamrx [0].append(x)
                lamrd1[0].append([ -d for d in d1 ])
                lamrd2[0].append(d2)
                lamrd3[0].append(d3)
                lamrd3[1][n2] = d3
        # substitute known start and end coordinates
        for n3 in range(2):
            lamrx [n3][ 0] = lavbx [n3][-3]
            lamrd1[n3][ 0] = lavbd1[n3][-3]
            lamrd2[n3][ 0] = lavbd2[n3][-3]
            lamrd3[n3][ 0] = lavbd3[n3][-3]
            lamrx [n3][-1] = lavqx [n3][-3]
            lamrd1[n3][-1] = lavqd1[n3][-3]
            lamrd2[n3][-1] = lavqd2[n3][-3]
            lamrd3[n3][-1] = lavqd3[n3][-3]
        # smooth d1 to fit RPV, LPV
        for n2 in range(1, elementsCountOverLeftAtriumVenous):
            n1lpv = n1 - elementsCountOverLeftAtriumVenous//2 - 1
            n1rpv = n1 - elementsCountOverLeftAtriumVenous//2 - elementsCountAroundLeftAtriumRPV
            for n3 in range(2):
                if lamRowCount == 1:
                    nx  = [ rpvox [n3][n1rpv], lamlx [n3][n2], lpvox [n3][n1lpv] ]
                    nd1 = [ rpvod1[n3][n1rpv], lamld1[n3][n2], [ -d for d in lpvod1[n3][n1lpv] ] ]
                    d1 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixAllDirections = True,
                        fixStartDerivative = True, fixEndDerivative = True, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)[1]
                    magd2 = vector.magnitude(lamld1[n3][n2])
                    if vector.magnitude(d1) > magd2:
                        d1 = vector.setMagnitude(d1, magd2)
                    lamrd1[n3][n2] = lamld1[n3][n2] = d1
                else:
                    nx  = [ rpvox [n3][n1rpv], lamrx [n3][n2], lamlx [n3][n2], lpvox [n3][n1lpv] ]
                    nd1 = [ rpvod1[n3][n1rpv], lamrd1[n3][2], lamld1[n3][n2], [ -d for d in lpvod1[n3][n1lpv] ] ]
                    lamrd1[n3][n2], lamld1[n3][n2] = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixAllDirections = True,
                        fixStartDerivative = True, fixEndDerivative = True, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)[1:3]

        # get points on right atrium along crista terminalis from aorta to posterior venous limit
        xi = (1.0 - raVenousMidpointOver - racsProportions[ran1Aorta][0])/(racsProportions[ran1Ctp][0] - racsProportions[ran1Aorta][0])
        rctmpProportion1 = 1.0 - raVenousMidpointOver
        rctmpProportion2 = raVenousRight  # GRC was (failing) (1.0 - xi)*racsProportions[ran1Aorta][1] + xi*racsProportions[ran1Ctp][1]
        elementsCountOverCristaTerminalisAnterior = elementsCountOverRightAtriumVenous//2 + 1
        elementsCountOverCristaTerminalisPosterior = elementsCountOverRightAtriumVenous//2
        _ractx, _ractd2, _ractd1, _ractd3, _ractProportions = raTrackSurface.createHermiteCurvePoints(
            rctmpProportion1, rctmpProportion2, racsProportions[ran1Ctp][0], racsProportions[ran1Ctp][1],
            elementsCount = elementsCountOverCristaTerminalisPosterior,
            derivativeStart = None,
            derivativeEnd = [ -d for d in racsd2[1][ran1Ctp] ])
        ractx, ractd2, ractd1, ractd3, ractProportions = raTrackSurface.createHermiteCurvePoints(
            ragProportions[1], 0.0,
            rctmpProportion1, rctmpProportion2,
            elementsCount = elementsCountOverCristaTerminalisAnterior,
            derivativeStart = [ -1.5*d for d in agd1[1] ],  # GRC fudge factor
            derivativeEnd = _ractd2[0])
        ractx  += _ractx [1:]
        ractd1 += _ractd1[1:]
        ractd2 += _ractd2[1:]
        ractd3 += _ractd3[1:]
        ractProportions += _ractProportions[1:]
        # get distance around posterior coronary sinus to crista terminalis, to scale set ct derivative 1 below
        raVenousWidth = sum(interp.getCubicHermiteArcLength(racsx[1][e], racsd1[1][e], racsx[1][e + 1], racsd1[1][e + 1]) for e in range(elementsCountAroundRightAtriumPosteriorVenous))
        d1mag = -0.2*raVenousWidth  # GRC fudge factor
        ractx  = [ [], ractx  ]
        ractd1 = [ [], [ vector.setMagnitude(d1, d1mag) for d1 in ractd1 ] ]
        ractd2 = [ [], ractd2 ]
        ractd3 = [ [], ractd3 ]
        for n in range(len(ractx[1])):
            x, d2, d1, d3 = interp.projectHermiteCurvesThroughWall(ractx[1], ractd2[1], [ [ -d for d in d1 ] for d1 in ractd1[1] ], n, -cristaTerminalisThickness)
            ractx [0].append(x)
            ractd1[0].append([ -d for d in d1 ])
            ractd2[0].append(d2)
            ractd3[0].append(d3)
            ractd3[1][n] = d3
        # substitute known end coordinates
        for n3 in range(2):
            # overwrite venous right x, d1 on coronary sinus
            racsx [n3][ran1Ctp] = ractx [n3][-1]
            d1mag = min(vector.magnitude(racsd1[n3][ran1Ctp]), 1.0*vector.magnitude(ractd1[n3][-1]))  # GRC fudge factor
            racsd1[n3][ran1Ctp] = vector.setMagnitude(racsd1[n3][ran1Ctp], d1mag)
            #ractd1[n3][-1] = racsd1[n3][ran1Ctp]
            ractd2[n3][-1] = racsd2[n3][ran1Ctp]
            ractd3[n3][-1] = racsd3[n3][ran1Ctp]

        # get points on right atrium ridge midway between inferior and superior vena cavae from crista terminalis to interatrial groove
        # minimum of 2 points over top of venous component
        elementsCountOverSideRightAtriumVC = max(elementsCountAroundRightAtriumPosteriorVenous, 2)
        ravmx, ravmd1, ravmd2, ravmd3 = raTrackSurface.createHermiteCurvePoints(
            ractProportions[elementsCountOverCristaTerminalisAnterior][0], ractProportions[elementsCountOverCristaTerminalisAnterior][1],
            1.0 - raVenousMidpointOver, 0.0,
            elementsCount = elementsCountOverSideRightAtriumVC,
            derivativeStart = ractd1[1][elementsCountOverCristaTerminalisAnterior],
            derivativeEnd = agd1[agn1Mid])[0:4]
        # get inner points
        ravmx  = [ [], ravmx  ]
        ravmd1 = [ [], ravmd1 ]
        ravmd2 = [ [], ravmd2 ]
        ravmd3 = [ [], ravmd3 ]
        # blend d2 between ends:
        magc = vector.magnitude(ractd2[1][elementsCountOverCristaTerminalisAnterior])
        maga = vector.magnitude(agd2[agn1Mid])
        for n in range(1, elementsCountOverSideRightAtriumVC):
            xi = n/elementsCountOverSideRightAtriumVC
            ravmd2[1][n] = vector.setMagnitude(ravmd2[1][n], (1.0 - xi)*magc + xi*maga)
        for n in range(len(ravmx[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(ravmx[1], ravmd1[1], ravmd2[1], n, -raVenousFreeWallThickness)
            ravmx [0].append(x)
            ravmd1[0].append(d1)
            ravmd2[0].append(d2)
            ravmd3[0].append(d3)
            ravmd3[1][n] = d3
        # substitute known end coordinates
        for n3 in range(2):
            ravmd2[n3][0] = ractd2[n3][elementsCountOverCristaTerminalisAnterior]
        asn1Mid = elementsCountAroundAtrialSeptum + agn1Mid - 1
        ravmx [0][-1]  = asx [1][asn1Mid]
        ravmd1[0][-1]  = asd1[1][asn1Mid]
        ravmd2[0][-1]  = asd2[1][asn1Mid]
        ravmd3[0][-1]  = asd3[1][asn1Mid]
        ravmx [1][-1]  = agx [agn1Mid]
        ravmd1[1][-1]  = agd1[agn1Mid]
        ravmd2[1][-1]  = agd2[agn1Mid]
        ravmd3[1][-1]  = agd3[agn1Mid]

        # create la base nodes, adding to lFibrousRingGroup
        labNodeId = [ [], [] ]
        lFibrousRingNodesetGroup = lFibrousRingGroup.getNodesetGroup(nodes)
        # create ra base nodes, adding to rFibrousRingGroup
        rabNodeId = [ [], [] ]
        nodesetGroup = rFibrousRingGroup.getNodesetGroup(nodes)
        for n3 in range(2):
            for n1 in range(len(labx[n3])):
                if not labx[n3][n1]:
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, labx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, labd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, labd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, labd3[n3][n1])
                lFibrousRingNodesetGroup.addNode(node)
                labNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

            for n1 in range(len(rabx[n3])):
                if not rabx[n3][n1]:
                    continue
                if n3 == 1:
                    # find common nodes on left atrium base
                    nodeId = None
                    if n1 == 0:  # crux
                        nodeId = labNodeId[1][elementsCountAroundLeftAtriumFreeWall]
                    elif n1 == elementsCountAroundRightAtriumFreeWall:  # cfb
                        nodeId = labNodeId[1][0]
                    if nodeId:
                        rabNodeId[n3].append(nodeId)
                        nodesetGroup.addNode(nodes.findNodeByIdentifier(nodeId))
                        continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rabx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rabd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rabd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rabd3[n3][n1])
                nodesetGroup.addNode(node)
                rabNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

        # create interatrial groove nodes:
        # start and end with common nodes at cfb and crux
        agNodeId = [ labNodeId[1][0] ]  # cfb
        for n1 in range(1, len(agx) - 1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, agx [n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, agd1[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, agd2[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, agd3[n1])
            agNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
        agNodeId.append(rabNodeId[1][0])  # crux

        # create septum nodes, along coronary sinus and over arch
        asNodeId = [ [], [] ]
        for n3 in range(2):
            for ns in range(len(asx[n3])):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, asx [n3][ns])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, asd1[n3][ns])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, asd2[n3][ns])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, asd3[n3][ns])
                asNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

        # create fossa ovalis nodes
        foCentreNodeId = []
        foNodeId = [ [], [] ]
        for n3 in range(2):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, foCentrex [n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, foCentred1[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, foCentred2[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, foCentred3[n3])
            foCentreNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
            for nf in range(elementsCountAroundFossa):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, fox [n3][nf])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, fod1[n3][nf])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, fod2[n3][nf])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, fod3[n3][nf])
                foNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1

        # create left atrium coronary sinus nodes
        # start and end with common nodes on interatrial groove or septum arch
        lacsNodeId = [ [ asNodeId[0][elementsCountAroundAtrialSeptum] ], [ agNodeId[1] ] ]
        for n3 in range(2):
            for n1 in range(1, len(lacsx[n3]) - 1):
                if elementsCountAroundLeftAtriumAorta < n1 < lan1Mid:
                    # left atrial appendage
                    lacsNodeId[n3].append(None)
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                lacsNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lacsx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lacsd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lacsd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lacsd3[n3][n1])
                nodeIdentifier += 1
        lacsNodeId[0].append(asNodeId[0][0])
        lacsNodeId[1].append(agNodeId[-2])

        # create right atrium coronary sinus nodes
        racsNodeId = [ [ asNodeId[1][0] ], [ agNodeId[-2] ] ]
        for n3 in range(2):
            for n1 in range(1, len(racsx[n3]) - 1):
                if ran1Ctp < n1 <= ran1Aorta:
                    # right atrial appendage
                    racsNodeId[n3].append(None)
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                racsNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, racsx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, racsd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, racsd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, racsd3[n3][n1])
                nodeIdentifier += 1
        racsNodeId[0].append(asNodeId[1][elementsCountAroundAtrialSeptum])
        racsNodeId[1].append(agNodeId[1])

        # create nodes on left atrium over appendage
        laoaNodeId = [ [], [] ]
        for n3 in range(2):
            laoaNodeId[n3].append(lacsNodeId[n3][1])
            for n1 in range(1, len(laoax[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, laoax [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, laoad1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, laoad2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, laoad3[n3][n1])
                laoaNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            laoaNodeId[n3].append(lacsNodeId[n3][lan1Mid])

        # create nodes on left atrium venous anterior to LPV ostium
        lavbNodeId = [ [ asNodeId[0][asn1va] ], [ agNodeId[elementsCountOverLeftAtriumNonVenousAnterior] ] ]
        n1lpv = -1 - elementsCountOverLeftAtriumVenous//2
        for n3 in range(2):
            for n1 in range(1, len(lavbx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lavbx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lavbd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lavbd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lavbd3[n3][n1])
                lavbNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            lavbNodeId[n3].append(lpvoNodeId[n3][n1lpv])

        # create nodes on left atrium venous posterior to LPV ostium
        lavqNodeId = [ [ asNodeId[0][-1] ], [ agNodeId[elementsCountOverLeftAtriumNonVenousAnterior + elementsCountOverLeftAtriumVenous ] ] ]
        n1lpv = 1
        for n3 in range(2):
            for n1 in range(1, len(lavqx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lavqx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lavqd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lavqd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lavqd3[n3][n1])
                lavqNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            lavqNodeId[n3].append(lpvoNodeId[n3][n1lpv])

        # create left atrium venous midpoint nodes, left and right
        lamlNodeId = [ [], [] ]
        for n3 in range(2):
            lamlNodeId[n3].append(lavbNodeId[n3][-2])
            for n2 in range(1, len(lamlx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lamlx [n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lamld1[n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lamld2[n3][n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lamld3[n3][n2])
                lamlNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            lamlNodeId[n3].append(lavqNodeId[n3][-2])
        lamrNodeId = [ [], [] ]
        for n3 in range(2):
            lamrNodeId[n3].append(lavbNodeId[n3][-3])
            if lamRowCount == 1:
                lamrNodeId[n3] += lamlNodeId[n3][1:len(lamrx[n3]) - 1]
            else:
                for n2 in range(1, len(lamrx[n3]) - 1):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lamrx [n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lamrd1[n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lamrd2[n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lamrd3[n3][n2])
                    lamrNodeId[n3].append(nodeIdentifier)
                    nodeIdentifier += 1
            lamrNodeId[n3].append(lavqNodeId[n3][-3])

        # create right atrium crista terminalis nodes
        ractNodeId = [ [], [] ]
        for n3 in range(2):
            ractNodeId[n3].append(agNodeId[1])
            for n1 in range(1, len(ractx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ractx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ractd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ractd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, ractd3[n3][n1])
                ractNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            ractNodeId[n3].append(racsNodeId[n3][ran1Ctp])
            # use second ct node as aorta coronary sinus node
            racsNodeId[n3][ran1Aorta] = ractNodeId[n3][1]

        # create right atrium venous midline nodes
        ravmNodeId = [ [], [] ]
        for n3 in range(2):
            ravmNodeId[n3].append(ractNodeId[n3][elementsCountOverCristaTerminalisAnterior])
            for n1 in range(1, len(ravmx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ravmx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ravmd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ravmd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, ravmd3[n3][n1])
                ravmNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
        ravmNodeId[0].append(asNodeId[1][asn1Mid])
        ravmNodeId[1].append(agNodeId[agn1Mid])

        if False:
            # create lt nodes:
            for n1 in range(len(ltBaseOuterx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ltBaseOuterx [n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ltBaseOuterd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ltBaseOuterd2[n1])
                nodeIdentifier += 1

        drawLaTrackSurface = False
        if drawLaTrackSurface:
            # create track surface nodes:
            laTrackSurfaceFirstNodeIdentifier = nodeIdentifier
            for n in range((laTrackSurface.elementsCount2 + 1)*(laTrackSurface.elementsCount1 + 1)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, laTrackSurface.nx[n] )
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, laTrackSurface.nd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, laTrackSurface.nd2[n])
                nodeIdentifier += 1
        drawRaTrackSurface = False
        if drawRaTrackSurface:
            # create track surface nodes:
            raTrackSurfaceFirstNodeIdentifier = nodeIdentifier
            for n in range((raTrackSurface.elementsCount2 + 1)*(raTrackSurface.elementsCount1 + 1)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, raTrackSurface.nx[n] )
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, raTrackSurface.nd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, raTrackSurface.nd2[n])
                nodeIdentifier += 1

        #################
        # Create elements
        #################

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # left atrium free wall elements to coronary sinus, starting at cfb / anterior interatrial sulcus
        eftBaseSulcus = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftBaseSulcus, [1], [])
        remapEftNodeValueLabel(eftBaseSulcus, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
        remapEftNodeValueLabel(eftBaseSulcus, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftBaseSulcus, [ 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        ln_map = [ 1, 2, 3, 4, 5, 5, 6, 6 ]
        remapEftLocalNodes(eftBaseSulcus, 6, ln_map)
        # general linear map d3 adjacent to collapsed cfb/crux
        eftBaseSulcusNext = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftBaseSulcusNext, [1], [])
        remapEftNodeValueLabel(eftBaseSulcusNext, [ 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftBaseSulcusNext, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        # general linear map d3 adjacent to collapsed cfb/crux
        eftBaseSulcusPrev = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftBaseSulcusPrev, [1], [])
        remapEftNodeValueLabel(eftBaseSulcusPrev, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        remapEftNodeValueLabel(eftBaseSulcusPrev, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
        for e1 in range(-1, elementsCountAroundLeftAtriumFreeWall):
            eft1 = eft
            elementtemplate1 = elementtemplate
            nids = [
                labNodeId[0][e1], labNodeId[0][e1 + 1], lacsNodeId[0][e1], lacsNodeId[0][e1 + 1],
                labNodeId[1][e1], labNodeId[1][e1 + 1], lacsNodeId[1][e1], lacsNodeId[1][e1 + 1]]
            if None in nids:
                continue  # left atrial appendage
            scalefactors = None
            meshGroups = [ laMeshGroup ]
            if e1 == -1:
                # cfb/anterior interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = rabNodeId[0][-elementsCountAroundAtrialSeptum]
                nids[2] = racsNodeId[0][-1]
                nids.pop(6)
                nids.pop(4)
                meshGroups += [ raMeshGroup ]
                eft1 = eftBaseSulcus
                scalefactors = [ -1.0 ]
            elif e1 == 0:
                eft1 = eftBaseSulcusNext
                scalefactors = [ -1.0 ]
            elif e1 == (elementsCountAroundLeftAtriumFreeWall - 1):
                eft1 = eftBaseSulcusPrev
                scalefactors = [ -1.0 ]
            if eft1 is not eft:
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateX
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            #print('create element la', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # right atrium free wall elements to coronary sinus, starting at crux / posterior interatrial sulcus
        for e1 in range(-1, elementsCountAroundRightAtriumFreeWall):
            eft1 = eft
            elementtemplate1 = elementtemplate
            nids = [
                rabNodeId[0][e1], rabNodeId[0][e1 + 1], racsNodeId[0][e1], racsNodeId[0][e1 + 1],
                rabNodeId[1][e1], rabNodeId[1][e1 + 1], racsNodeId[1][e1], racsNodeId[1][e1 + 1]]
            if None in nids:
                continue  # right atrial appendage
            scalefactors = None
            meshGroups = [ raMeshGroup ]
            if e1 == -1:
                # crux/posterior interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = labNodeId[0][elementsCountAroundLeftAtriumFreeWall]
                nids[2] = lacsNodeId[0][elementsCountAroundLeftAtriumFreeWall]
                nids.pop(6)
                nids.pop(4)
                meshGroups += [ laMeshGroup ]
                eft1 = eftBaseSulcus
                scalefactors = [ -1.0 ]
            elif e1 == 0:
                eft1 = eftBaseSulcusNext
                scalefactors = [ -1.0 ]
            elif e1 == (elementsCountAroundRightAtriumFreeWall - 1):
                # similar to eftBaseSulcusPrev, but general linear map node from crista terminalis
                # general linear map d3 adjacent to collapsed cfb/crux
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])  # GRC , ( Node.VALUE_LABEL_D_DS2, []) ])
                remapEftNodeValueLabel(eft1, [ 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                scalefactors = [ -1.0 ]
            if eft1 is not eft:
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateX
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            #print('create element ra', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # left atrium non-venous anterior elements (first row above coronary sinus)
        meshGroups = [ laMeshGroup ]
        for e1 in range(elementsCountAroundLeftAtriumRPV + 1):
            eft1 = eft
            elementtemplate1 = elementtemplate
            nids = [ laoaNodeId[0][e1 - 1], laoaNodeId[0][e1], lavbNodeId[0][e1], lavbNodeId[0][e1 + 1],
                     laoaNodeId[1][e1 - 1], laoaNodeId[1][e1], lavbNodeId[1][e1], lavbNodeId[1][e1 + 1] ]
            scalefactors = None
            if e1 == 0:
                nids[0:2] = [ lacsNodeId[0][e1], lacsNodeId[0][e1 + 1] ]
                nids[4:6] = [ lacsNodeId[1][e1], lacsNodeId[1][e1 + 1] ]
                eft1 = eftBaseSulcusNext
                scalefactors = [ -1.0 ]
            else:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                if e1 == elementsCountAroundLeftAtriumRPV:
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                if e1 == 1:
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
            elementtemplateX.defineField(coordinates, -1, eft1)
            elementtemplate1 = elementtemplateX
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            #print('create element laoa', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # left atrium non-venous posterior elements (first row above coronary sinus)
        meshGroups = [ laMeshGroup ]
        scalefactors = [ -1.0 ]
        for e1 in range(elementsCountAroundLeftAtriumRPV + 1):
            nc = elementsCountAroundLeftAtriumFreeWall - e1
            nids = [ lavqNodeId[0][e1], lavqNodeId[0][e1 + 1], lacsNodeId[0][nc], lacsNodeId[0][nc - 1],
                     lavqNodeId[1][e1], lavqNodeId[1][e1 + 1], lacsNodeId[1][nc], lacsNodeId[1][nc - 1] ]
            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
            if e1 < elementsCountAroundLeftAtriumRPV:
                scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
            if e1 == 0:
                # general linear map d3 adjacent to collapsed inter-atrial groove
                remapEftNodeValueLabel(eft1, [ 1, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                remapEftNodeValueLabel(eft1, [ 3, 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
            elif e1 == elementsCountAroundLeftAtriumRPV:
                scaleEftNodeValueLabels(eft1, [ 3, 7 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [1]) ])
            elementtemplateX.defineField(coordinates, -1, eft1)
            elementtemplate1 = elementtemplateX
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, scalefactors)
            #print('create element lavq', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # left atrium mid line elements, depending on lamRowCount
        lamrDiagonalSide = (lamRowCount == 1)
        meshGroups = [ laMeshGroup ]
        for e1 in range(elementsCountOverLeftAtriumVenous):
            eft1 = eft
            elementtemplate1 = elementtemplate
            nids = [ lamrNodeId[0][e1], lamlNodeId[0][e1], lamrNodeId[0][e1 + 1], lamlNodeId[0][e1 + 1],
                     lamrNodeId[1][e1], lamlNodeId[1][e1], lamrNodeId[1][e1 + 1], lamlNodeId[1][e1 + 1] ]
            scalefactors = None
            if e1 == 0:
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                # d2 diagonal
                remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
            elif e1 == (elementsCountOverLeftAtriumVenous - 1):
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scalefactors = [ -1.0 ]
                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
            if lamRowCount == 1:
                if e1 == 0:
                    remapEftNodeValueLabel(eft1, [ 3, 4, 7, 8 ], Node.VALUE_LABEL_D_DS1, [])
                    if lamrDiagonalSide:
                        # d2 diagonal
                        remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                elif e1 == (elementsCountOverLeftAtriumVenous - 1):
                    remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [])
                    if lamrDiagonalSide:
                        # d2 diagonal
                        remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                else:
                    continue
            elementtemplateX.defineField(coordinates, -1, eft1)
            elementtemplate1 = elementtemplateX
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            #print('create element lam', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # create atrial septum base row of elements
        meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup ]
        for e1 in range(elementsCountAroundAtrialSeptum):
            n1l = elementsCountAroundLeftAtriumFreeWall + e1 - elementsCountAroundLeftAtrium
            n1r = -e1
            nids = [ labNodeId[0][n1l], labNodeId[0][n1l + 1], asNodeId[0][e1], asNodeId[0][e1 + 1], \
                     rabNodeId[0][n1r], rabNodeId[0][n1r - 1], asNodeId[1][e1], asNodeId[1][e1 + 1] ]
            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
            if e1 == 0:
                scaleEftNodeValueLabels(eft1, [ 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            elif e1 == (elementsCountAroundAtrialSeptum - 1):
                scaleEftNodeValueLabels(eft1, [ 5, 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
            else:
                scaleEftNodeValueLabels(eft1, [ 5, 6, 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
            elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            #print('create element as base', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # create fossa ovalis elements
        meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup, fossaMeshGroup ]
        radiansAround0 = fossaRadiansAround[-1]
        if radiansAround0 > fossaRadiansAround[0]:
            radiansAround0 -= 2.0*math.pi
        radiansAround1 = fossaRadiansAround[0]
        radiansAround2 = fossaRadiansAround[1]
        if radiansAround2 < radiansAround1:
            radiansAround2 += 2.0*math.pi
        for e1 in range(elementsCountAroundFossa):
            va = e1
            vb = (e1 + 1)%elementsCountAroundFossa
            eft1 = tricubichermite.createEftShellPoleTop(va*100, vb*100)
            elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            nids = [ foNodeId[0][va], foNodeId[0][vb], foCentreNodeId[0], foNodeId[1][va], foNodeId[1][vb], foCentreNodeId[1] ]
            result2 = element.setNodesByIdentifier(eft1, nids)
            radiansAround3 = fossaRadiansAround[va + 2 - elementsCountAroundFossa]
            if radiansAround3 < radiansAround2:
                radiansAround3 += 2.0*math.pi
            dRadiansAround1 = 0.5*(radiansAround2 - radiansAround0)
            dRadiansAround2 = 0.5*(radiansAround3 - radiansAround1)
            scalefactors = [
                -1.0,
                -math.cos(radiansAround1), -math.sin(radiansAround1), dRadiansAround1,
                -math.cos(radiansAround2), -math.sin(radiansAround2), dRadiansAround2,
                -math.cos(radiansAround1), -math.sin(radiansAround1), dRadiansAround1,
                -math.cos(radiansAround2), -math.sin(radiansAround2), dRadiansAround2
            ]
            result3 = element.setScaleFactors(eft1, scalefactors)
            #print('create element fossa', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier = elementIdentifier + 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            radiansAround0, radiansAround1, radiansAround2 = radiansAround1, radiansAround2, radiansAround3

        # create atrial septum elements around fossa
        meshGroups = [ laMeshGroup, raMeshGroup, aSeptumMeshGroup ]
        for e1 in range(elementsCountAroundFossa):
            e1p = (e1 + 1)%elementsCountAroundFossa
            nids = [ asNodeId[0][e1], asNodeId[0][e1p], foNodeId[0][e1], foNodeId[0][e1p], \
                     asNodeId[1][e1], asNodeId[1][e1p], foNodeId[1][e1], foNodeId[1][e1p] ]
            eft1 = tricubichermite.createEftNoCrossDerivatives()
            setEftScaleFactorIds(eft1, [1], [])
            scalefactors = [ -1.0 ]
            if e1 < elementsCountAroundAtrialSeptum:
                scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 5, 6 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                if e1 == 0:
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                elif e1 == (elementsCountAroundAtrialSeptum - 1):
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
            else:
                remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # use temporary to swap later
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D2_DS1DS2, [] ) ])  # use temporary to swap later
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                if e1 == elementsCountAroundAtrialSeptum:
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                elif e1 == (elementsCountAroundFossa - 1):
                    remapEftNodeValueLabel(eft1, [ 1 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 5 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])  # set again to negate
                else:
                    remapEftNodeValueLabel(eft1, [ 1, 2 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 5, 6 ], Node.VALUE_LABEL_D2_DS1DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
            elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            result2 = element.setNodesByIdentifier(eft1, nids)
            if eft1.getNumberOfLocalScaleFactors() == 1:
                result3 = element.setScaleFactors(eft1, [ -1.0 ])
            #print('create element as', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # create elements around septum arch collapsed on interatrial groove

        eftAGroove = tricubichermite.createEftNoCrossDerivatives()
        setEftScaleFactorIds(eftAGroove, [1], [])
        scaleEftNodeValueLabels(eftAGroove, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
        remapEftNodeValueLabel(eftAGroove, [ 1, 2, 7, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
        remapEftNodeValueLabel(eftAGroove, [ 3, 4, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
        remapEftNodeValueLabel(eftAGroove, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [])
        remapEftNodeValueLabel(eftAGroove, [ 1, 2, 3, 4, 5, 6, 7, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
        ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
        remapEftLocalNodes(eftAGroove, 6, ln_map)
        meshGroups = [ laMeshGroup, raMeshGroup ]
        elementsCountOverArch = elementsCountOverAtria - 2
        elementtemplateX.defineField(coordinates, -1, eftAGroove)
        scalefactors = [ -1.0 ]
        for e1 in range(elementsCountOverArch):
            ns = e1 + elementsCountAroundAtrialSeptum - elementsCountAroundFossa
            nids = [ agNodeId[e1 + 1], agNodeId[e1 + 2], asNodeId[0][ns], asNodeId[0][ns + 1], asNodeId[1][ns], asNodeId[1][ns + 1] ]
            eft1 = eftAGroove
            if e1 == (elementsCountOverArch - 1):
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scaleEftNodeValueLabels(eft1, [ 7, 8 ], [ Node.VALUE_LABEL_D_DS3 ], [ 1 ])
                remapEftNodeValueLabel(eft1, [ 1, 4, 6, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 3, 5, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS3, [])
                remapEftNodeValueLabel(eft1, [ 1, 3, 5, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                ln_map = [ 1, 2, 3, 4, 1, 2, 5, 6 ]
                remapEftLocalNodes(eft1, 6, ln_map)
                elementtemplateX.defineField(coordinates, -1, eft1)
            element = mesh.createElement(elementIdentifier, elementtemplateX)
            result2 = element.setNodesByIdentifier(eft1, nids)
            result3 = element.setScaleFactors(eft1, scalefactors)
            #print('create element ag', element.isValid(), elementIdentifier, result2, result3, nids)
            elementIdentifier += 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

        # create left pulmonary vein annulus
        lpvax  = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        lpvad1 = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        lpvad2 = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        lpvad3 = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        lpvaNodeId = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        lpvaDerivativesMap = [ [ None ]*elementsCountAroundLpvOstium, [ None ]*elementsCountAroundLpvOstium ]
        # set points clockwise from venous midpoint at anterior venous limit
        # insert at indexes such that 0 is one past the midpoint on venous midline
        ix = -(elementsCountOverLeftAtriumVenous//2 + 1)
        # left atrium venous midpoint line left
        for n1 in range(elementsCountOverLeftAtriumVenous + 1):
            for n3 in range(2):
                lpvax [n3][ix] = lamlx [n3][n1]
                lpvad1[n3][ix] = lamld1[n3][n1]
                lpvad2[n3][ix] = lamld2[n3][n1]
                lpvad3[n3][ix] = lamld3[n3][n1]
                lpvaNodeId[n3][ix] = lamlNodeId[n3][n1]
                if n1 == 0:
                    lpvaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None, (0, 1, 0 ) )
                #elif n1 == elementsCountOverLeftAtriumVenous:
                #    lpvaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 1, 0), None, (1, 0, 0 ) )
                else:
                    lpvaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
            ix += 1
        # left around left cs on LPV venous zone, including 2 corners,
        for n1 in range(elementsCountAroundLeftAtriumLPV):
            nc = lan1Mid + elementsCountAroundLeftAtriumLPV - n1 - 1
            for n3 in range(2):
                lpvax [n3][ix] = lacsx [n3][nc]
                lpvad1[n3][ix] = lacsd1[n3][nc]
                lpvad2[n3][ix] = lacsd2[n3][nc]
                lpvad3[n3][ix] = lacsd3[n3][nc]
                lpvaNodeId[n3][ix] = lacsNodeId[n3][nc]
                if n1 == 0:
                    lpvaDerivativesMap[n3][ix] = ( (-1, -1, 0), (0, -1, 0), None, (-1, 0, 0 ) )
                elif n1 == (elementsCountAroundLeftAtriumLPV - 1):
                    lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0 ) )
                else:
                    lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
            ix += 1
        # up left atrium venous side line
        for n1 in range(1, elementsCountOverSideLeftAtriumLPV):
            no = -1 - n1
            for n3 in range(2):
                lpvax [n3][ix] = laoax [n3][no]
                lpvad1[n3][ix] = laoad1[n3][no]
                lpvad2[n3][ix] = laoad2[n3][no]
                lpvad3[n3][ix] = laoad3[n3][no]
                lpvaNodeId[n3][ix] = laoaNodeId[n3][no]
                if n1 == (elementsCountOverSideLeftAtriumLPV - 1):
                    lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None, (-1, 1, 0) )
                else:
                    lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
            ix += 1
        #print('lpvaNodeId[0]',lpvaNodeId[0])
        #print('lpvaNodeId[1]',lpvaNodeId[1])
        #print('lpvaDerivativesMap[0]',lpvaDerivativesMap[0])
        #print('lpvaDerivativesMap[1]',lpvaDerivativesMap[1])
        nodeIdentifier, elementIdentifier = tricubichermite.createAnnulusMesh3d(
            lpvox, lpvod1, lpvod2, lpvod3, lpvoNodeId, None,
            lpvax, lpvad1, lpvad2, lpvad3, lpvaNodeId, lpvaDerivativesMap,
            nodetemplate, nodetemplateLinearS3, nodeIdentifier, elementIdentifier,
            elementsCountRadial = 1, meshGroups = [ laMeshGroup ])

        # create right pulmonary vein annulus
        rpvax  = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        rpvad1 = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        rpvad2 = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        rpvad3 = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        rpvaNodeId = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        rpvaDerivativesMap = [ [ None ]*elementsCountAroundRpvOstium, [ None ]*elementsCountAroundRpvOstium ]
        # set points clockwise from interatrial groove at anterior venous limit
        # insert at indexes such that 0 is the midpoint on interatrial groove
        ix = -elementsCountOverLeftAtriumVenous//2
        # down interatrial groove from anterior venous limit, including both corners
        for n1 in range(elementsCountOverLeftAtriumVenous + 1):
            ns = n1 - elementsCountOverLeftAtriumVenous - 1
            ng = elementsCountOverLeftAtriumNonVenousAnterior + n1
            rpvax [0][ix] = asx [0][ns]
            rpvad1[0][ix] = asd1[0][ns]
            rpvad2[0][ix] = asd2[0][ns]
            rpvad3[0][ix] = asd3[0][ns]
            rpvaNodeId[0][ix] = asNodeId[0][ns]
            rpvax [1][ix] = agx [ng]
            rpvad1[1][ix] = agd1[ng]
            rpvad2[1][ix] = agd2[ng]
            rpvad3[1][ix] = agd3[ng]
            rpvaNodeId[1][ix] = agNodeId[ng]
            if n1 == 0:
                rpvaDerivativesMap[0][ix] = ( (-1, 0, 0), (-1, -1, 0), (1, 0, 1), (0, 1, 0 ) )
                rpvaDerivativesMap[1][ix] = ( (-1, 0, 0), (-1, -1, 0), (-1, 0, 1), (0, 1, 0 ) )
            elif n1 == elementsCountOverLeftAtriumVenous:
                rpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 1, 0), (1, 0, 1), (1, 0, 0 ) )
                rpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 1, 0), (-1, 0, 1), (1, 0, 0 ) )
            else:
                rpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 0, 0), (1, 0, 1) )
                rpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 0, 0), (-1, 0, 1) )
            ix += 1
        # left over posterior venous limit
        for n1 in range(1, elementsCountAroundLeftAtriumRPV):
            for n3 in range(2):
                rpvax [n3][ix] = lavqx [n3][n1]
                rpvad1[n3][ix] = lavqd1[n3][n1]
                rpvad2[n3][ix] = lavqd2[n3][n1]
                rpvad3[n3][ix] = lavqd3[n3][n1]
                rpvaNodeId[n3][ix] = lavqNodeId[n3][n1]
                rpvaDerivativesMap[n3][ix] = ( None, None, None )
            ix += 1
        # up left atrium venous midline right
        for n1 in range(0, elementsCountOverLeftAtriumVenous + 1):
            nm = elementsCountOverLeftAtriumVenous - n1
            for n3 in range(2):
                rpvax [n3][ix] = lamrx [n3][nm]
                rpvad1[n3][ix] = lamrd1[n3][nm]
                rpvad2[n3][ix] = lamrd2[n3][nm]
                rpvad3[n3][ix] = lamrd3[n3][nm]
                rpvaNodeId[n3][ix] = lamrNodeId[n3][nm]
                if n1 == 0:
                    # d1 diagonal out
                    rpvaDerivativesMap[n3][ix] = ( (1, 0, 0), (0, 1, 0), None, (1, -1, 0) )  # GRC ( (1, 0, 0), (1, 1, 0), None, (1, -1, 0) )
                elif lamrDiagonalSide and (n1 == 1):
                    if elementsCountOverLeftAtriumVenous == 2:
                        # diagonal in and out
                        rpvaDerivativesMap[n3][ix] = ( (1, -1, 0), (1, 0, 0), None, (-1, -1, 0) )
                    else:
                        # d1 diagonal in
                        rpvaDerivativesMap[n3][ix] = ( (1, -1, 0), (1, 0, 0), None, (0, -1, 0) )
                elif lamrDiagonalSide and (n1 == (elementsCountOverLeftAtriumVenous - 1)):
                    # d1 diagonal out
                    rpvaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 0, 0), None, (-1, -1, 0) )
                elif n1 == elementsCountOverLeftAtriumVenous:
                    # d1 diagonal in
                    rpvaDerivativesMap[n3][ix] = ( (-1, -1, 0), (0, -1, 0), None, (-1, 0, 0) )  # GRC ( (-1, -1, 0), (1, -1, 0), None, (-1, 0, 0) )
                else:
                    rpvaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 0, 0), None )
            ix += 1
        # right over anterior venous limit
        for n1 in range(1, elementsCountAroundLeftAtriumRPV):
            na = elementsCountAroundLeftAtriumRPV - n1
            for n3 in range(2):
                rpvax [n3][ix] = lavbx [n3][na]
                rpvad1[n3][ix] = lavbd1[n3][na]
                rpvad2[n3][ix] = lavbd2[n3][na]
                rpvad3[n3][ix] = lavbd3[n3][na]
                rpvaNodeId[n3][ix] = lavbNodeId[n3][na]
                rpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
            ix += 1
        #print('rpvaNodeId[0]',rpvaNodeId[0])
        #print('rpvaNodeId[1]',rpvaNodeId[1])
        #print('rpvaDerivativesMap[0]',rpvaDerivativesMap[0])
        #print('rpvaDerivativesMap[1]',rpvaDerivativesMap[1])
        nodeIdentifier, elementIdentifier = tricubichermite.createAnnulusMesh3d(
            rpvox, rpvod1, rpvod2, rpvod3, rpvoNodeId, None,
            rpvax, rpvad1, rpvad2, rpvad3, rpvaNodeId, rpvaDerivativesMap,
            nodetemplate, nodetemplateLinearS3, nodeIdentifier, elementIdentifier,
            elementsCountRadial = 1, meshGroups = [ laMeshGroup ])

        # create inferior and superior vena cavae inlets
        elementsCountAlongVCInlet = 2  # GRC make into a setting?
        for v in range(2):
            if v == 0:
                proportion1 = 1.0 - ivcPositionOver
                proportion2 = ivcPositionRight
                vcAngle1Radians = -ivcAngleOverRadians
                vcAngle2Radians = -ivcAngleLeftRadians
                vcDerivativeFactor = ivcDerivativeFactor
                vcLength = ivcLength
                vcInnerRadius = ivcInnerRadius
                vcWallThickness = ivcWallThickness
                elementsCountAroundVC = elementsCountAroundRightAtriumPosteriorVenous + elementsCountOverSideRightAtriumVC + elementsCountOverRightAtriumVenous
                startRadians = math.pi*elementsCountAroundRightAtriumPosteriorVenous/elementsCountAroundVC
            else:
                proportion1 = 1.0 - svcPositionOver
                proportion2 = svcPositionRight
                vcAngle1Radians = -svcAngleOverRadians
                vcAngle2Radians = -svcAngleLeftRadians
                vcDerivativeFactor = svcDerivativeFactor
                vcLength = svcLength
                vcInnerRadius = svcInnerRadius
                vcWallThickness = svcWallThickness
                elementsCountAroundVC = elementsCountOverRightAtriumVenous//2 + elementsCountOverCristaTerminalisAnterior + elementsCountOverSideRightAtriumVC
                startRadians = math.pi*elementsCountOverSideRightAtriumVC/elementsCountAroundVC
            # vessel
            vcOuterRadius = vcInnerRadius + vcWallThickness
            vcEndDerivative = vcDerivativeFactor*vcLength/elementsCountAlongVCInlet
            ocxPosition = raTrackSurface.createPositionProportion(proportion1, proportion2)
            ocx, d1, d2 = raTrackSurface.evaluateCoordinates(ocxPosition, derivatives = True)
            ocd1, ocd2, ocd3 = calculate_surface_axes(d1, d2, vector.normalise(d1))
            vcx, vd1, vd2, vd3 = getCircleProjectionAxes(ocx, ocd1, ocd2, ocd3, vcAngle1Radians, vcAngle2Radians, vcLength)
            vcd1 = vd1
            vcd2 = [ -d for d in vd2 ]
            vcd3 = [ -vcEndDerivative*d for d in vd3 ]
            vcvx  = [ None, None ]
            vcvd1 = [ None, None ]
            vcvd2 = [ None, None ]
            for n3 in range(2):
                radius = vcInnerRadius if (n3 == 0) else vcOuterRadius
                px, pd1 = createCirclePoints(vcx, vector.setMagnitude(vcd1, radius), vector.setMagnitude(vcd2, radius), elementsCountAroundVC, startRadians)
                vcvx [n3] = px
                vcvd1[n3] = pd1
                vcvd2[n3] = [ vcd3 ]*elementsCountAroundVC
            # annulus/aperture
            vcax  = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad1 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad2 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad3 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcaNodeId = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcaDerivativesMap = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            if v == 0:  # ivc
                # set points clockwise from interatrial groove at venous midpoint
                ix = 0
                # over interatrial groove from cristvenous midpoint to anterior, including both corners
                for n1 in range(elementsCountOverRightAtriumVenous//2, -1, -1):
                    ns = (elementsCountAroundAtrialSeptum - 1 + elementsCountOverRightAtriumNonVenousAnterior + elementsCountOverRightAtriumVenous//2 + n1) % elementsCountAroundFossa
                    ng = elementsCountOverRightAtriumNonVenousAnterior + elementsCountOverRightAtriumVenous//2 + n1
                    #print('v',v,'ix', ix, 'n1', n1, 'ns', ns, 'ng', ng)
                    vcax [0][ix] = asx [1][ns]
                    vcad1[0][ix] = asd1[1][ns]
                    vcad2[0][ix] = asd2[1][ns]
                    vcad3[0][ix] = asd3[1][ns]
                    vcaNodeId[0][ix] = asNodeId[1][ns]
                    vcax [1][ix] = agx [ng]
                    vcad1[1][ix] = agd1[ng]
                    vcad2[1][ix] = agd2[ng]
                    vcad3[1][ix] = agd3[ng]
                    vcaNodeId[1][ix] = agNodeId[ng]
                    if n1 == (elementsCountOverRightAtriumVenous//2):
                        # on coronary sinus, d1 and d2 are reversed
                        vcaDerivativesMap[0][ix] = ( (-1, 0, 0), (-1, -1, 0), (1, 0, 1), (0, 1, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (-1, 0, 0), (-1, -1, 0), (-1, 0, 1), (0, 1, 0 ) )
                    elif n1 == 0:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, -1, 0), (-1, 0, 1), (-1, 0, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 0 ) )
                    else:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, 0, 0), (-1, 0, 1) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, 0, 0), (1, 0, 1) )
                    ix += 1
                # right over venous midline to crista terminalis, excluding corners
                for n1 in range(1, elementsCountOverSideRightAtriumVC):
                    nm = elementsCountOverSideRightAtriumVC - n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nm', nm)
                    for n3 in range(2):
                        vcax [n3][ix] = ravmx [n3][nm]
                        vcad1[n3][ix] = ravmd1[n3][nm]
                        vcad2[n3][ix] = ravmd2[n3][nm]
                        vcad3[n3][ix] = ravmd3[n3][nm]
                        vcaNodeId[n3][ix] = ravmNodeId[n3][nm]
                        vcaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                    ix += 1
                # back over crista terminalis to coronary sinus including corners
                for n1 in range(elementsCountOverCristaTerminalisPosterior + 1):
                    nc = elementsCountOverCristaTerminalisAnterior + n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nc', nc)
                    for n3 in range(2):
                        vcax [n3][ix] = ractx [n3][nc]
                        vcad1[n3][ix] = ractd1[n3][nc]
                        vcad2[n3][ix] = ractd2[n3][nc]
                        vcad3[n3][ix] = ractd3[n3][nc]
                        vcaNodeId[n3][ix] = ractNodeId[n3][nc]
                        if n1 == 0:
                            vcaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0) )
                        elif n1 == elementsCountOverCristaTerminalisPosterior:
                            # on coronary sinus, d1 and d2 are reversed
                            vcaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, -1, 0), None, (-1, 0, 0) )
                        else:
                            vcaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
                    ix += 1
                # left around coronary sinus on venous posterior
                for n1 in range(1, elementsCountAroundRightAtriumPosteriorVenous):
                    nc = elementsCountAroundRightAtriumPosteriorVenous - n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nc', nc)
                    for n3 in range(2):
                        vcax [n3][ix] = racsx [n3][nc]
                        vcad1[n3][ix] = racsd1[n3][nc]
                        vcad2[n3][ix] = racsd2[n3][nc]
                        vcad3[n3][ix] = racsd3[n3][nc]
                        vcaNodeId[n3][ix] = racsNodeId[n3][nc]
                        vcaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                    ix += 1
            else:  # svc
                # set points clockwise from interatrial groove at venous limit posterior
                ix = 0
                # over interatrial groove from venous midpoint to anterior, including both corners
                for n1 in range(elementsCountOverRightAtriumVenous//2, -1, -1):
                    ns = elementsCountAroundAtrialSeptum - 1 + elementsCountOverRightAtriumNonVenousAnterior + n1
                    ng = elementsCountOverRightAtriumNonVenousAnterior + n1
                    #print('v',v,'ix', ix, 'n1', n1, 'ns', ns, 'ng', ng)
                    vcax [0][ix] = asx [1][ns]
                    vcad1[0][ix] = asd1[1][ns]
                    vcad2[0][ix] = asd2[1][ns]
                    vcad3[0][ix] = asd3[1][ns]
                    vcaNodeId[0][ix] = asNodeId[1][ns]
                    vcax [1][ix] = agx [ng]
                    vcad1[1][ix] = agd1[ng]
                    vcad2[1][ix] = agd2[ng]
                    vcad3[1][ix] = agd3[ng]
                    vcaNodeId[1][ix] = agNodeId[ng]
                    if n1 == (elementsCountOverRightAtriumVenous//2):
                        vcaDerivativesMap[0][ix] = ( (1, 0, 0), (1, 1, 0), (-1, 0, 1), (0, -1, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (1, 0, 0), (1, 1, 0), (1, 0, 1), (0, -1, 0 ) )
                    elif n1 == 0:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, -1, 0), (-1, 0, 1), (-1, 0, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 0 ) )
                    else:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, 0, 0), (-1, 0, 1) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, 0, 0), (1, 0, 1) )
                    ix += 1
                # back over crista terminalis to coronary sinus including corners
                for n1 in range(1, elementsCountOverCristaTerminalisAnterior + 1):
                    nc = n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nc', nc)
                    for n3 in range(2):
                        vcax [n3][ix] = ractx [n3][nc]
                        vcad1[n3][ix] = ractd1[n3][nc]
                        vcad2[n3][ix] = ractd2[n3][nc]
                        vcad3[n3][ix] = ractd3[n3][nc]
                        vcaNodeId[n3][ix] = ractNodeId[n3][nc]
                        if n1 == 0:
                            vcaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0) )
                        elif n1 == elementsCountOverCristaTerminalisAnterior:
                            vcaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 1, 0), None, (1, 0, 0) )
                        else:
                            vcaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
                    ix += 1
                # left over venous midline to interatrial groove, excluding corners
                for n1 in range(1, elementsCountOverSideRightAtriumVC):
                    nm = n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nm', nm)
                    for n3 in range(2):
                        vcax [n3][ix] = ravmx [n3][nm]
                        vcad1[n3][ix] = ravmd1[n3][nm]
                        vcad2[n3][ix] = ravmd2[n3][nm]
                        vcad3[n3][ix] = ravmd3[n3][nm]
                        vcaNodeId[n3][ix] = ravmNodeId[n3][nm]
                        vcaDerivativesMap[n3][ix] = ( None, None, None )
                    ix += 1
            #print('vcaNodeId[0]',vcaNodeId[0])
            #print('vcaNodeId[1]',vcaNodeId[1])
            #print('vcaDerivativesMap[0]',vcaDerivativesMap[0])
            #print('vcaDerivativesMap[1]',vcaDerivativesMap[1])
            nodeIdentifier, elementIdentifier = tricubichermite.createAnnulusMesh3d(
                vcvx, vcvd1, vcvd2, None, None, None,
                vcax, vcad1, vcad2, vcad3, vcaNodeId, vcaDerivativesMap,
                nodetemplate, nodetemplateLinearS3, nodeIdentifier, elementIdentifier,
                elementsCountRadial = elementsCountAlongVCInlet, maxEndThickness = 1.5*raVenousFreeWallThickness,
                meshGroups = [ raMeshGroup, ivcInletMeshGroup if (v == 0) else svcInletMeshGroup])

        if drawLaTrackSurface:
            mesh2d = fm.findMeshByDimension(2)
            bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft2d = mesh2d.createElementfieldtemplate(bicubicHermiteBasis)
            # remove cross derivative 12
            for n in range(4):
                r = eft2d.setFunctionNumberOfTerms(n*4 + 4, 0)
            elementtemplate2d = mesh2d.createElementtemplate()
            elementtemplate2d.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            elementtemplate2d.defineField(coordinates, -1, eft2d)
            nodesCount1 = laTrackSurface.elementsCount1 + 1
            for e2 in range(laTrackSurface.elementsCount2):
                for e1 in range(laTrackSurface.elementsCount1):
                    element = mesh2d.createElement(-1, elementtemplate2d)  # since on 2-D mesh
                    nid1 = laTrackSurfaceFirstNodeIdentifier + e2*nodesCount1 + e1
                    element.setNodesByIdentifier(eft2d, [ nid1, nid1 + 1, nid1 + nodesCount1, nid1 + nodesCount1 + 1 ])
        if drawRaTrackSurface:
            mesh2d = fm.findMeshByDimension(2)
            bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft2d = mesh2d.createElementfieldtemplate(bicubicHermiteBasis)
            # remove cross derivative 12
            for n in range(4):
                r = eft2d.setFunctionNumberOfTerms(n*4 + 4, 0)
            elementtemplate2d = mesh2d.createElementtemplate()
            elementtemplate2d.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            elementtemplate2d.defineField(coordinates, -1, eft2d)
            nodesCount1 = raTrackSurface.elementsCount1 + 1
            for e2 in range(raTrackSurface.elementsCount2):
                for e1 in range(raTrackSurface.elementsCount1):
                    element = mesh2d.createElement(-1, elementtemplate2d)  # since on 2-D mesh
                    nid1 = raTrackSurfaceFirstNodeIdentifier + e2*nodesCount1 + e1
                    element.setNodesByIdentifier(eft2d, [ nid1, nid1 + 1, nid1 + nodesCount1, nid1 + nodesCount1 + 1 ])

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
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountAroundRightAtriumFreeWall = options['Number of elements around right atrium free wall']
        refineElementsCountSurface = options['Refine number of elements surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        sourceFm = meshrefinement._sourceFm
        annotationGroups = meshrefinement._sourceAnnotationGroups
        laGroup = findAnnotationGroupByName(annotationGroups, 'left atrium')
        laElementGroupField = laGroup.getFieldElementGroup(meshrefinement._sourceMesh)
        raGroup = findAnnotationGroupByName(annotationGroups, 'right atrium')
        raElementGroupField = raGroup.getFieldElementGroup(meshrefinement._sourceMesh)
        aSeptumGroup = findAnnotationGroupByName(annotationGroups, 'interatrial septum')
        aSeptumElementGroupField = aSeptumGroup.getFieldElementGroup(meshrefinement._sourceMesh)
        isSeptumEdgeWedge = sourceFm.createFieldXor(sourceFm.createFieldAnd(laElementGroupField, raElementGroupField), aSeptumElementGroupField)

        # last atria element is last element in following group:
        lastGroup = findAnnotationGroupByName(annotationGroups, 'superior vena cava inlet')
        lastMeshGroup = lastGroup.getMeshGroup(meshrefinement._sourceMesh)
        lastElementIdentifier = -1
        elementIter = lastMeshGroup.createElementiterator()
        element = elementIter.next()
        while element.isValid():
            lastElementIdentifier = element.getIdentifier()
            element = elementIter.next()

        cache = sourceFm.createFieldcache()
        refineElements3 = refineElementsCountThroughWall
        element = meshrefinement._sourceElementiterator.next()
        wedgeElementCount = 0
        while element.isValid():
            elementIdentifier = element.getIdentifier()
            refineElements1 = refineElementsCountSurface
            refineElements2 = refineElementsCountSurface
            cache.setElement(element)
            result, isWedge = isSeptumEdgeWedge.evaluateReal(cache, 1)
            if isWedge:
                wedgeElementCount += 1
                # the first two around the base are collapsed on 1-3, remainder on 2-3 
                if wedgeElementCount <= 2:
                    refineElements1 = refineElementsCountThroughWall
                else:
                    refineElements2 = refineElementsCountThroughWall
            meshrefinement.refineElementCubeStandard3d(element, refineElements1, refineElements2, refineElements3)
            if elementIdentifier == lastElementIdentifier:
                return  # finish on last so can continue elsewhere
            element = meshrefinement._sourceElementiterator.next()

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


def getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall):
    '''
    Get the number of elements in each section of the left atrium free wall.
    :param elementsCountAroundLeftAtriumFreeWall: Valid range 6-10.
    :return: elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftArialAppendageBase,
        elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV
    '''
    assert 6 <= elementsCountAroundLeftAtriumFreeWall <= 10, \
        'getLeftAtriumBaseFreewallElementsCounts: elements count out of range: ' + str(elementsCountAroundLeftAtriumFreeWall)
    elementsCountAroundLeftAtriumAorta = 1
    elementsCountAroundLeftAtriumVP = (elementsCountAroundLeftAtriumFreeWall + 1)//2
    elementsCountAroundLeftArialAppendageBase = elementsCountAroundLeftAtriumFreeWall - elementsCountAroundLeftAtriumVP - elementsCountAroundLeftAtriumAorta
    elementsCountAroundLeftAtriumRPV = elementsCountAroundLeftAtriumVP//2
    elementsCountAroundLeftAtriumLPV = elementsCountAroundLeftAtriumVP - elementsCountAroundLeftAtriumRPV
    return elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftArialAppendageBase, \
        elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV


def getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall):
    '''
    Get the number of elements in each section of the right atrium free wall.
    :param elementsCountAroundRightAtriumFreeWall: Valid range 6-10.
    :return: elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightArialAppendageBase,
        elementsCountAroundRightAtriumAorta
    '''
    assert 6 <= elementsCountAroundRightAtriumFreeWall <= 10, \
        'getRightAtriumBaseFreewallElementsCounts: elements count out of range: ' + str(elementsCountAroundRightAtriumFreeWall)
    elementsCountAroundRightAtriumAorta = 1
    elementsCountAroundRightAtriumPosteriorVenous = elementsCountAroundRightAtriumFreeWall//4
    elementsCountAroundRightArialAppendageBase = elementsCountAroundRightAtriumFreeWall \
        - elementsCountAroundRightAtriumPosteriorVenous - elementsCountAroundRightAtriumAorta
    return elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightArialAppendageBase, \
        elementsCountAroundRightAtriumAorta


def getOverAtriaElementsCounts(elementsCountOverAtria):
    '''
    Get the number of elements in each section over the atria, at the outer
    interatrial groove.
    :param elementsCountOverAtria: Valid values 6 or 8.
    :return: elementsCountOverAtriaCoronarySinus,
        elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior,
        elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior
    '''
    assert elementsCountOverAtria in [ 6, 8 ], \
        'getOverAtriaElementsCounts: elements count not 6 or 8: ' + str(elementsCountOverAtria)
    elementsCountOverAtriaCoronarySinus = 1
    elementsCountOverLeftAtriumNonVenousAnterior = 2
    elementsCountOverLeftAtriumNonVenousPosterior = 2
    elementsCountOverLeftAtriumVenous = elementsCountOverAtria - elementsCountOverLeftAtriumNonVenousPosterior - elementsCountOverLeftAtriumNonVenousAnterior
    elementsCountOverRightAtriumNonVenousAnterior = 1
    elementsCountOverRightAtriumNonVenousPosterior = 1
    elementsCountOverRightAtriumVenous = elementsCountOverAtria - elementsCountOverRightAtriumNonVenousPosterior - elementsCountOverRightAtriumNonVenousAnterior
    return elementsCountOverAtriaCoronarySinus, \
        elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior, \
        elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior

def getAtriumBasePoints(elementsCountAroundAtrialSeptum, elementsCountAroundLeftAtriumFreeWall, elementsCountAroundRightAtriumFreeWall,
    aBaseInnerMajorMag, aBaseInnerMinorMag, aMajorAxisRadians,
    aBaseWallThickness, aBaseSlopeHeight, aBaseSlopeLength, aSeptumLength, aSeptumThickness,
    aortaOuterPlusRadius, aBaseFrontInclineRadians, aBaseSideInclineRadians, aBaseBackInclineRadians,
    laaLeft, laVenousMidpointPosteriorLeft, raVenousRight, elementsCountAroundTrackSurface):
    """
    Get points around left and right atria based on an ellipse.
    Left atria points start from central fibrous body and wind anticlockwise.
    Right atria points start from crux and wind anticlockwise.
    around LA. Both the cfb and crux are collapsed at the septum.
    Also return la free wall sample outer points used in construction and for defining track surface.
    :return: laBasex[n3], laBased1[n3], laBased2[n3], laBased3[n3],
             raBasex[n3], raBased1[n3], raBased2[n3], raBased3[n3],
             ltBaseOuterx, ltBaseOuterd1, ltBaseOuterd2,
             aSeptumBaseCentre, laCentre, laSeptumRadians
    """
    lvOutletFrontInclineRadians = aBaseFrontInclineRadians  # for now

    aBaseOuterMajorMag = aBaseInnerMajorMag + aBaseSlopeLength
    aBaseOuterMinorMag = aBaseInnerMinorMag + aBaseSlopeLength

    # following are angles in radians around LA ellipse from major axis
    axInner = aBaseInnerMajorMag*math.cos(aMajorAxisRadians)
    ayInner = -aBaseInnerMajorMag*math.sin(aMajorAxisRadians)
    bxInner = aBaseInnerMinorMag*math.sin(aMajorAxisRadians)
    byInner = aBaseInnerMinorMag*math.cos(aMajorAxisRadians)
    laSeptumRadians = math.atan2(bxInner, axInner)
    laCentreX = -0.5*aSeptumThickness - axInner*math.cos(laSeptumRadians) - bxInner*math.sin(laSeptumRadians)
    axOuter = aBaseOuterMajorMag*math.cos(aMajorAxisRadians)
    ayOuter = -aBaseOuterMajorMag*math.sin(aMajorAxisRadians)
    bxOuter = aBaseOuterMinorMag*math.sin(aMajorAxisRadians)
    byOuter = aBaseOuterMinorMag*math.cos(aMajorAxisRadians)

    # get points on central fibrous body centre and cfbLeft (60 degrees clockwise around aorta)
    # rotates about centre of aorta by lvOutletFrontInclineRadians
    cosFrontInclineRadians = math.cos(lvOutletFrontInclineRadians)
    sinFrontInclineRadians = math.sin(lvOutletFrontInclineRadians)
    pi_3 = math.pi/3.0
    cosPi_3 = math.cos(pi_3)
    sinPi_3 = math.sin(pi_3)
    cfbSideOffset = aortaOuterPlusRadius*sinPi_3
    rLeft = aortaOuterPlusRadius*cosPi_3
    cfbLeftX = -cfbSideOffset
    cfbLeftY = -rLeft*cosFrontInclineRadians
    cfbLeftZ = rLeft*sinFrontInclineRadians
    cfbX = 0.0
    cfbY = -aortaOuterPlusRadius*cosFrontInclineRadians
    cfbZ = aortaOuterPlusRadius*sinFrontInclineRadians
    lvOutletDerivativeAround = aortaOuterPlusRadius*pi_3
    laCfbLeftRadians = getEllipseRadiansToX(axOuter, bxOuter, cfbLeftX - laCentreX, 0.5*math.pi)
    aBaseOuterMinorMagPlus = aBaseOuterMinorMag + 0.5*aBaseSlopeLength  # GRC fudge factor
    laCentreY = cfbLeftY - math.cos(laCfbLeftRadians)*aBaseOuterMajorMag*math.sin(-aMajorAxisRadians) \
                         - math.sin(laCfbLeftRadians)*aBaseOuterMinorMagPlus*math.cos(-aMajorAxisRadians)
    laCentreZ = 0.0

    # Convert aSeptumLength to arc angles each side of septum centre
    aSeptumBaseCentreX = -0.5*aSeptumThickness
    aSeptumBaseCentreY = laCentreY + ayInner*math.cos(laSeptumRadians) + byInner*math.sin(laSeptumRadians)
    #print('aSeptumBaseCentreY',aSeptumBaseCentreY,'radians',laSeptumRadians)
    aSeptumBaseCentreZ = -aBaseSlopeHeight
    # GRC ensure aSeptumLength < min(innerminormag, innermajormag)
    aSeptumAnteriorY = aSeptumBaseCentreY + 0.5*aSeptumLength
    #print('ayInner',ayInner)
    #print('byInner',byInner)
    #print('0.5*aSeptumLength', 0.5*aSeptumLength)
    #print('aSeptumBaseCentreY', aSeptumBaseCentreY)
    #print('aSeptumAnteriorY', aSeptumAnteriorY)
    #print('aSeptumAnteriorY - laCentreY', aSeptumAnteriorY - laCentreY)
    laSeptumAnteriorRadians = getEllipseRadiansToX(ayInner, byInner, aSeptumAnteriorY - laCentreY, 0.5*(laSeptumRadians + laCfbLeftRadians))
    #print('aSeptumAnteriorY', aSeptumAnteriorY, 'laSeptumAnteriorRadians', laSeptumAnteriorRadians, 'y',
    #      laCentreY + ayInner*math.cos(laSeptumAnteriorRadians) + byInner*math.sin(laSeptumAnteriorRadians))
    aSeptumPosteriorY = aSeptumBaseCentreY - 0.5*aSeptumLength
    laSeptumPosteriorRadians = getEllipseRadiansToX(ayInner, byInner, aSeptumPosteriorY - laCentreY, 1.5*laSeptumRadians - 0.5*laCfbLeftRadians)
    #print('aSeptumPosteriorY', aSeptumPosteriorY, 'laSeptumPosteriorRadians', laSeptumPosteriorRadians, 'y',
    #      laCentreY + ayInner*math.cos(laSeptumPosteriorRadians) + byInner*math.sin(laSeptumPosteriorRadians))

    # compute common lengths around outer
    atrialPerimeterLength = getApproximateEllipsePerimeter(aBaseOuterMajorMag, aBaseOuterMinorMag)
    atrialSeptumInnerElementLength = getEllipseArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, laSeptumPosteriorRadians, laSeptumAnteriorRadians)/elementsCountAroundAtrialSeptum
    atrialSeptumOuterElementLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumPosteriorRadians, laSeptumAnteriorRadians)/elementsCountAroundAtrialSeptum
    aSeptumCfbLeftElementLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumAnteriorRadians, laCfbLeftRadians)
    aCfbLeftDerivativeLength = 2.0*aSeptumCfbLeftElementLength - atrialSeptumOuterElementLength
    aRemainingLength = atrialPerimeterLength - aSeptumCfbLeftElementLength - elementsCountAroundAtrialSeptum*atrialSeptumOuterElementLength

    #testradians = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, laSeptumPosteriorRadians, atrialSeptumInnerElementLength)
    #print('test radians 1', testradians,'vs',laSeptumRadians,'y',laCentreY + ayInner*math.cos(testradians) + byInner*math.sin(testradians), 'vs',aSeptumBaseCentreY)
    #testradians = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, testradians, atrialSeptumInnerElementLength)
    #print('test radians 2', testradians,'vs',laSeptumAnteriorRadians,'| y',laCentreY + ayInner*math.cos(testradians) + byInner*math.sin(testradians), 'vs',aSeptumAnteriorY)

    baseDerivative2Scale = aortaOuterPlusRadius
    # GRC use these earlier
    sinMajorAxisRadians = math.sin(-aMajorAxisRadians)
    cosMajorAxisRadians = math.cos(-aMajorAxisRadians)

    # generate curves with fixed number of elements around free wall, from which track surface will be sampled
    # this is needed now to sample key points defined relative to the track surface, used to set element size and spacing around

    elementsCountAroundLeftAtriumFreeWallFixed = 8
    ltFreeWallElementLength = (aRemainingLength - 0.5*(atrialSeptumOuterElementLength + aCfbLeftDerivativeLength))/(elementsCountAroundLeftAtriumFreeWallFixed - 2)
    # first two points are around aorta
    ltBaseOuterx  = [ [ cfbX, cfbY, cfbZ ], [ cfbLeftX, cfbLeftY, cfbLeftZ ] ]
    ltBaseOuterd1 = [ [ -lvOutletDerivativeAround, 0.0, 0.0 ], [ -lvOutletDerivativeAround*cosPi_3, lvOutletDerivativeAround*sinPi_3*cosFrontInclineRadians, -lvOutletDerivativeAround*sinPi_3*sinFrontInclineRadians ]]
    ltBaseOuterd2 = [ [ 0.0, baseDerivative2Scale*sinFrontInclineRadians, baseDerivative2Scale*cosFrontInclineRadians ], [ 0.0, baseDerivative2Scale*sinFrontInclineRadians, baseDerivative2Scale*cosFrontInclineRadians ] ]
    aMajorX =  aBaseOuterMajorMag*cosMajorAxisRadians
    aMajorY =  aBaseOuterMajorMag*sinMajorAxisRadians
    aMinorX = -aBaseOuterMinorMag*sinMajorAxisRadians
    aMinorY =  aBaseOuterMinorMag*cosMajorAxisRadians
    radiansAround = laCfbLeftRadians
    sideRadians = laSeptumRadians + math.pi
    backRadians = sideRadians + 0.5*math.pi
    up = [ 0.0, 0.0, 1.0 ]
    for n1 in range(elementsCountAroundLeftAtriumFreeWallFixed - 2):
        elementLength = 0.5*(aCfbLeftDerivativeLength + ltFreeWallElementLength) if (n1 == 0) else ltFreeWallElementLength
        radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        x = [ laCentreX + cosRadiansAround*aMajorX + sinRadiansAround*aMinorX,
              laCentreY + cosRadiansAround*aMajorY + sinRadiansAround*aMinorY,
              0.0 ]
        d1 = vector.setMagnitude([ -sinRadiansAround*aMajorX + cosRadiansAround*aMinorX,
                                   -sinRadiansAround*aMajorY + cosRadiansAround*aMinorY,
                                   0.0 ], ltFreeWallElementLength)
        if radiansAround < sideRadians:
            xi = (radiansAround - laCfbLeftRadians)/(sideRadians - laCfbLeftRadians)
            baseInclineRadians = (1.0 - xi)*aBaseFrontInclineRadians + xi*aBaseSideInclineRadians
        elif radiansAround < backRadians:
            xi = (radiansAround - sideRadians)/(0.5*math.pi)
            baseInclineRadians = (1.0 - xi)*aBaseSideInclineRadians + xi*aBaseBackInclineRadians
        else:
            baseInclineRadians = aBaseBackInclineRadians
        side = vector.normalise([ d1[1], -d1[0], 0.0 ])
        d2 = vector.setMagnitude([ (up[c]*math.cos(baseInclineRadians) + side[c]*math.sin(baseInclineRadians)) for c in range(3) ], baseDerivative2Scale)
        ltBaseOuterx .append(x)
        ltBaseOuterd1.append(d1)
        ltBaseOuterd2.append(d2)
    # generate point on posterior base septum
    xi = 0.9  # GRC fudge factor
    baseSeptumPosteriorx  = [ 0.0, (1.0 - xi)*ltBaseOuterx[0][1] + xi*ltBaseOuterx[-1][1], 0.0 ]
    nx  = [ ltBaseOuterx[-1], baseSeptumPosteriorx ]
    nd1 = interp.smoothCubicHermiteDerivativesLine(nx, [ ltBaseOuterd1[-1], [ vector.magnitude(ltBaseOuterd1[-1]), 0.0, 0.0 ] ],
        fixStartDerivative = True, fixEndDirection = True )
    ltBaseOuterx .append(nx[1])
    # GRC fudge factor: derivative must be lower to fit inlets:
    ltBaseOuterd1.append([ 0.35*d for d in nd1[1] ])
    # derivative 2 slopes directly back = no x component
    ltBaseOuterd2.append([ 0.0, -baseDerivative2Scale*math.sin(aBaseBackInclineRadians), baseDerivative2Scale*math.cos(aBaseBackInclineRadians) ])

    # get key points around la, ra to put element boundaries on (converted to arc angle below)
    vx, vd1 = interp.sampleCubicHermiteCurves(ltBaseOuterx, ltBaseOuterd1, elementsCountAroundTrackSurface)[0:2]
    # laa end
    er = 0.5*laaLeft*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    laaEndY = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[1]
    laaEndRadians = getEllipseRadiansToX(ayOuter, byOuter, laaEndY - laCentreY, sideRadians)
    #print('laaEnd y', laaEndY, 'radians', laaEndRadians)
    # venous midpoint
    er = 0.5*(2.0 - laVenousMidpointPosteriorLeft)*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    laVenousMidpointX = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[0]
    laVenousMidpointRadians = getEllipseRadiansToX(axOuter, bxOuter, laVenousMidpointX - laCentreX, backRadians)
    #print('laVenousMidpointPosteriorLeft', laVenousMidpointPosteriorLeft, 'e', e, 'xi', xi, 'x', laVenousMidpointX, 'radians', laVenousMidpointRadians)
    er = 0.5*(2.0 - raVenousRight)*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    # note these are computed on the left atrium and mirrored at the end
    raVenousPosteriorRightX = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[0]
    raVenousPosteriorRightRadians = getEllipseRadiansToX(axOuter, bxOuter, raVenousPosteriorRightX - laCentreX, backRadians)
    #print('raVenousRight', raVenousRight, 'e', e, 'xi', xi, 'x', raVenousPosteriorRightX, 'radians', raVenousPosteriorRightRadians)

    # get numbers of elements and lengths of sections of left atrium (outer)
    elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftArialAppendageBase, elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV = \
        getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall)
    laaLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laCfbLeftRadians, laaEndRadians)
    laVenousLeftLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laaEndRadians, laVenousMidpointRadians)
    laVenousRightLength = aRemainingLength - laVenousLeftLength - laaLength
    # = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laVenousMidpointRadians, laSeptumPosteriorRadians + 2.0*math.pi)

    # get element lengths/derivatives at edges of each section and transition element sizes between
    laaEndDerivative = (laaLength - 0.5*aCfbLeftDerivativeLength)/(elementsCountAroundLeftArialAppendageBase - 0.5)
    laVenousMidpointDerivative = (laVenousRightLength - 0.5*atrialSeptumOuterElementLength)/(elementsCountAroundLeftAtriumRPV - 0.5)
    laaElementLengths = interp.sampleCubicElementLengths(laaLength, elementsCountAroundLeftArialAppendageBase, startDerivative = aCfbLeftDerivativeLength, endDerivative = laaEndDerivative)
    lvlElementLengths = interp.sampleCubicElementLengths(laVenousLeftLength, elementsCountAroundLeftAtriumLPV, startDerivative = laaEndDerivative, endDerivative = laVenousMidpointDerivative)
    lvrElementLengths = interp.sampleCubicElementLengths(laVenousRightLength, elementsCountAroundLeftAtriumRPV, startDerivative = laVenousMidpointDerivative, endDerivative = atrialSeptumOuterElementLength)

    # get radians of nodes around left atrium, starting at cfb
    elementsCountAroundLeftAtrium = elementsCountAroundLeftAtriumFreeWall + elementsCountAroundAtrialSeptum
    laRadians = []
    radiansAround = laSeptumAnteriorRadians
    for n1 in range(elementsCountAroundLeftAtrium):
        laRadians.append(radiansAround)
        if n1 == 0:
            elementLength = aSeptumCfbLeftElementLength
        elif n1 < (elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftArialAppendageBase):
            elementLength = laaElementLengths[n1 - elementsCountAroundLeftAtriumAorta]
        elif n1 < (elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftArialAppendageBase + elementsCountAroundLeftAtriumLPV):
            elementLength = lvlElementLengths[n1 - elementsCountAroundLeftAtriumAorta - elementsCountAroundLeftArialAppendageBase]
        elif n1 == (elementsCountAroundLeftAtriumFreeWall - 1):
            radiansAround = laSeptumPosteriorRadians + 2.0*math.pi
            continue
        elif n1 < elementsCountAroundLeftAtriumFreeWall:
            elementLength = lvrElementLengths[n1 - elementsCountAroundLeftAtriumAorta - elementsCountAroundLeftArialAppendageBase - elementsCountAroundLeftAtriumLPV]
        else:
            radiansAround = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, radiansAround, atrialSeptumInnerElementLength)
            continue
        radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)

    #print('la radiansAround start', laSeptumAnteriorRadians, 'end', radiansAround - 2.0*math.pi)
    #print('laRadians', laRadians)

    # get numbers of elements and lengths of sections of right atrium (outer)

    elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightArialAppendageBase, elementsCountAroundRightAtriumAorta = \
        getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall)
    raaLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laCfbLeftRadians, raVenousPosteriorRightRadians)
    raPosteriorVenousLength = aRemainingLength - raaLength
    # = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, raVenousPosteriorRightRadians, laSeptumPosteriorRadians + 2.0*math.pi)

    # get element lengths/derivatives at edges of each section and transition element sizes between
    raPosteriorVenousLimitDerivative = (raPosteriorVenousLength - 0.5*atrialSeptumOuterElementLength)/(elementsCountAroundRightAtriumPosteriorVenous - 0.5)
    raaElementLengths = interp.sampleCubicElementLengths(raaLength, elementsCountAroundRightArialAppendageBase, startDerivative = aCfbLeftDerivativeLength, endDerivative = raPosteriorVenousLimitDerivative)
    ravElementLengths = interp.sampleCubicElementLengths(raPosteriorVenousLength, elementsCountAroundRightAtriumPosteriorVenous, startDerivative = raPosteriorVenousLimitDerivative, endDerivative = atrialSeptumOuterElementLength)

    # get radians of nodes around right atrium (computed on left and mirrored at the end), starting at cfb
    elementsCountAroundRightAtrium = elementsCountAroundRightAtriumFreeWall + elementsCountAroundAtrialSeptum
    raRadians = []
    radiansAround = laSeptumAnteriorRadians
    for n1 in range(elementsCountAroundRightAtrium):
        raRadians.append(radiansAround)
        if n1 == 0:
            elementLength = aSeptumCfbLeftElementLength
        elif n1 < (elementsCountAroundRightAtriumAorta + elementsCountAroundRightArialAppendageBase):
            elementLength = raaElementLengths[n1 - elementsCountAroundRightAtriumAorta]
        elif n1 == (elementsCountAroundRightAtriumFreeWall - 1):
            radiansAround = laSeptumPosteriorRadians + 2.0*math.pi
            continue
        elif n1 < elementsCountAroundRightAtriumFreeWall:
            elementLength = ravElementLengths[n1 - elementsCountAroundRightAtriumAorta - elementsCountAroundRightArialAppendageBase]
        else:
            radiansAround = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, radiansAround, atrialSeptumInnerElementLength)
            continue
        radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)

    #print('ra radiansAround start', laSeptumAnteriorRadians, 'end', radiansAround - 2.0*math.pi)
    #print('raRadians', raRadians)

    # get base points on inside and outside of left and right atria
    for a in [ 'la', 'ra' ]:
        if a == 'la':
            aRadians = laRadians
            elementsCountAroundAtrium = elementsCountAroundLeftAtrium
            elementsCountAroundAtriumFreeWall = elementsCountAroundLeftAtriumFreeWall
        else:  # a == 'ra':
            aRadians = raRadians
            elementsCountAroundAtrium = elementsCountAroundRightAtrium
            elementsCountAroundAtriumFreeWall = elementsCountAroundRightAtriumFreeWall
        aBaseInnerx = []
        aBaseInnerd1 = []
        aBaseInnerd2 = []
        aBaseOuterx  = copy.deepcopy(ltBaseOuterx [0:2])
        aBaseOuterd1 = copy.deepcopy(ltBaseOuterd1[0:2])
        aBaseOuterd2 = copy.deepcopy(ltBaseOuterd2[0:2])
        for n3 in range(2):
            if n3 == 0:
                aMajorMag = aBaseInnerMajorMag
                aMinorMag = aBaseInnerMinorMag
                z = -aBaseSlopeHeight
            else:
                aMajorMag = aBaseOuterMajorMag
                aMinorMag = aBaseOuterMinorMag
                z = 0.0

            aMajorX =  aMajorMag*cosMajorAxisRadians
            aMajorY =  aMajorMag*sinMajorAxisRadians
            aMinorX = -aMinorMag*sinMajorAxisRadians
            aMinorY =  aMinorMag*cosMajorAxisRadians

            finalArcLength = prevArcLength = getEllipseArcLength(aMajorMag, aMinorMag, aRadians[-1] - 2.0*math.pi, aRadians[0])
            n1Start = 0 if (n3 == 0) else 2
            n1Limit = elementsCountAroundAtrium if (n3 == 0) else elementsCountAroundAtriumFreeWall
            for n1 in range(n1Start, n1Limit):
                radiansAround = aRadians[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)

                # get derivative around
                if n1 == (elementsCountAroundAtrium - 1):
                    nextArcLength = finalArcLength
                else:
                    nextArcLength = getEllipseArcLength(aMajorMag, aMinorMag, aRadians[n1], aRadians[n1 + 1])
                # GRC review; can use exact values on outer. Compare?
                derivativeLength = 0.5*(prevArcLength + nextArcLength)
                prevArcLength = nextArcLength

                x = [ laCentreX + cosRadiansAround*aMajorX + sinRadiansAround*aMinorX,
                      laCentreY + cosRadiansAround*aMajorY + sinRadiansAround*aMinorY,
                      z ]
                d1 = vector.setMagnitude([ -sinRadiansAround*aMajorX + cosRadiansAround*aMinorX,
                                           -sinRadiansAround*aMajorY + cosRadiansAround*aMinorY,
                                           0.0 ], derivativeLength)
                if (n1 < 1) or (n1 >= elementsCountAroundAtriumFreeWall):
                    d2 = [ 0.0, 0.0, baseDerivative2Scale ]  # calculated later
                else:
                    if radiansAround < sideRadians:
                        xi = (radiansAround - laCfbLeftRadians)/(sideRadians - laCfbLeftRadians)
                        baseInclineRadians = (1.0 - xi)*aBaseFrontInclineRadians + xi*aBaseSideInclineRadians
                    elif radiansAround < backRadians:
                        xi = (radiansAround - sideRadians)/(0.5*math.pi)
                        baseInclineRadians = (1.0 - xi)*aBaseSideInclineRadians + xi*aBaseBackInclineRadians
                    else:
                        baseInclineRadians = aBaseBackInclineRadians
                    side = vector.normalise([ d1[1], -d1[0], 0.0 ])
                    d2 = vector.setMagnitude([ (up[c]*math.cos(baseInclineRadians) + side[c]*math.sin(baseInclineRadians)) for c in range(3) ], baseDerivative2Scale)
                if n3 == 0:
                    aBaseInnerx.append(x)
                    aBaseInnerd1.append(d1)
                    aBaseInnerd2.append(d2)
                else:
                    aBaseOuterx.append(x)
                    aBaseOuterd1.append(d1)
                    aBaseOuterd2.append(d2)

        aBaseOuterx .append(ltBaseOuterx [-1])
        aBaseOuterd1.append(ltBaseOuterd1[-1])
        aBaseOuterd2.append(ltBaseOuterd2[-1])
        for n1 in range(elementsCountAroundAtriumFreeWall + 1, elementsCountAroundAtrium):
            aBaseOuterx .append(None)
            aBaseOuterd1.append(None)
            aBaseOuterd2.append(None)

        # calculate d3 from difference across wall
        aBaseInnerd3 = []
        aBaseOuterd3 = []
        for n1 in range(elementsCountAroundAtrium):
            if aBaseOuterx[n1]:
                d3 = [ (aBaseOuterx[n1][c] - aBaseInnerx[n1][c]) for c in range(3) ]
            else:
                d3 = [ -2.0*aBaseInnerx[n1][0], 0.0, 0.0 ]
            aBaseInnerd3.append(d3)
            aBaseOuterd3.append(copy.deepcopy(d3))
        # fix outer d3 on cfb and crux
        aBaseOuterd3[0][0] = 0.0
        aBaseOuterd3[elementsCountAroundAtriumFreeWall][0] = 0.0

        if a == 'la':
            laBasex  = [ aBaseInnerx , aBaseOuterx  ]
            laBased1 = [ aBaseInnerd1, aBaseOuterd1 ]
            laBased2 = [ aBaseInnerd2, aBaseOuterd2 ]
            laBased3 = [ aBaseInnerd3, aBaseOuterd3 ]
        else:  # a == 'ra':
            raBasex  = [ aBaseInnerx , aBaseOuterx  ]
            raBased1 = [ aBaseInnerd1, aBaseOuterd1 ]
            raBased2 = [ aBaseInnerd2, aBaseOuterd2 ]
            raBased3 = [ aBaseInnerd3, aBaseOuterd3 ]
            # reverse and mirror about x == 0
            # reverse all x components, but only y, z components of d1 as winds in opposite direction
            for li in (raBasex + raBased2 + raBased3):
                for n1 in range(elementsCountAroundAtrium):
                    if li[n1]:
                        li[n1][0] = -li[n1][0]
            for li in raBased1:
                for n1 in range(elementsCountAroundAtrium):
                    if li[n1]:
                        li[n1][1] = -li[n1][1]
                        li[n1][2] = -li[n1][2]
            for li in (raBasex + raBased1 + raBased2 + raBased3):
                li.reverse()
                for n1 in range(elementsCountAroundAtrialSeptum - 1):
                    li.append(li.pop(0))

    return laBasex, laBased1, laBased2, laBased3, raBasex, raBased1, raBased2, raBased3, \
           ltBaseOuterx, ltBaseOuterd1, ltBaseOuterd2, \
           [ aSeptumBaseCentreX, aSeptumBaseCentreY, aSeptumBaseCentreZ ], [ laCentreX, laCentreY, laCentreZ ], laSeptumRadians


def getAtriumTrackSurface(elementsCountAroundTrackSurface, elementsCountAcrossTrackSurface,
        laBaseOuterx, laBaseOuterd1, laBaseOuterd2, aSeptumBaseCentre, aOuterHeight, aOuterSeptumHeight, iaGrooveDerivative):
    '''
    Create a TrackSurface covering the outer surface of the left atrium on which
    inlets will be placed. Elements vary fastest from posterior to anterior, then from septum to outer left.
    :param laBaseOuterx, laBaseOuterd1, laBaseOuterd2: coordinates, derivatives and transverse derivatives
    around left atrium outside wall from cfb to crux.
    :return: TrackSurface
    '''
    # resample to get equal spaced elements around
    elementsCountAlongTrackSurface = elementsCountAroundTrackSurface//2
    vx, vd1, ve, vxi = interp.sampleCubicHermiteCurves(laBaseOuterx, laBaseOuterd1, elementsCountAroundTrackSurface)[0:4]
    # GRC these are not necessarily orthogonal to d1 any more: fix?
    vd2 = interp.interpolateSampleLinear(laBaseOuterd2, ve, vxi)

    # get la ridge points from cubic functions from ax = septum groove centre through cx on peak to dx on mid outer LV base
    ax = [ 0.0, aSeptumBaseCentre[1], aOuterSeptumHeight ]
    ad1 = [ -iaGrooveDerivative, 0.0, 0.0 ]
    dx = vx[elementsCountAlongTrackSurface]
    dd1 = [ -d for d in vd2[elementsCountAlongTrackSurface]]
    # fudge factor
    px, pd1 = interp.sampleCubicHermiteCurves([ ax, dx ], [ ad1, dd1 ], elementsCountOut = 2, lengthFractionStart = 0.4, arcLengthDerivatives = True)[0:2]
    nx = [ ax, [ px[1][0], px[1][1], aOuterHeight ] ]
    nd1 = interp.smoothCubicHermiteDerivativesLine(nx, [ ad1, [ pd1[1][0], pd1[1][1], 0.0 ] ], fixStartDerivative = True, fixEndDirection = True)
    cx = nx[1]
    cd1 = nd1[1]
    # bx = in-between point to get more curvature near septum
    xi = 0.4
    bx = interp.interpolateCubicHermite(ax, ad1, cx, cd1, xi)
    bd1 = interp.interpolateCubicHermiteDerivative(ax, ad1, cx, cd1, xi)
    rx, rd1 = interp.sampleCubicHermiteCurves([ ax, bx, cx, dx ], [ ad1, bd1, cd1, dd1 ], elementsCountOut = elementsCountAlongTrackSurface, arcLengthDerivatives = True)[0:2]

    # get track surface points on arcs from posterior on septum end, to anterior on outer left
    nx = []
    nd1 = []
    nd2 = []
    for na in range(elementsCountAlongTrackSurface):
        np = -1 - na
        # sample arch from double cubic from posterior through ridge to anterior
        ax = [ vx[np], rx[na], vx[na] ]
        ad1 = [ vd2[np], [ rd1[na][1], -rd1[na][0], 0.0 ], [ -d for d in vd2[na]] ]
        ad1 = interp.smoothCubicHermiteDerivativesLine(ax, ad1, fixStartDirection = True, fixEndDirection = True)
        lx, ld1, le, lxi = interp.sampleCubicHermiteCurves(ax, ad1, elementsCountAcrossTrackSurface)[0:4]
        ld2 = interp.interpolateSampleLinear([ [ -d for d in vd1[np] ], rd1[na], vd1[na] ], le, lxi)
        nx += lx
        nd1 += ld1
        nd2 += ld2  # to be smoothed when all rows assembled.
    # add last point at tip of atrium, d1 all zero, d2 rotating around
    d1 = vd1[elementsCountAlongTrackSurface]
    d2 = vector.setMagnitude(vd2[elementsCountAlongTrackSurface], vector.magnitude(d1))
    for n1 in range(elementsCountAcrossTrackSurface + 1):
        nx.append(copy.deepcopy(vx[elementsCountAlongTrackSurface]))
        nd1.append([ 0.0, 0.0, 0.0 ])
        radiansAround = math.pi*n1/elementsCountAcrossTrackSurface
        wt1 = -math.cos(radiansAround)
        wt2 = -math.sin(radiansAround)
        nd2.append([ (wt1*d1[c] + wt2*d2[c]) for c in range(3) ])

    # smooth nd2 around each side, away from base edges
    for n1 in range(1, elementsCountAcrossTrackSurface):
        sx = []
        sd2 = []
        for n2 in range(elementsCountAlongTrackSurface + 1):
            n = n2*(elementsCountAcrossTrackSurface + 1) + n1
            sx.append(nx[n])
            sd2.append(nd2[n])
        sd2 = interp.smoothCubicHermiteDerivativesLine(sx, sd2, fixStartDirection = True, fixEndDirection = True)
        for n2 in range(elementsCountAcrossTrackSurface + 1):
            n = n2*(elementsCountAcrossTrackSurface + 1) + n1
            nd2[n] = sd2[n2]

    return TrackSurface(elementsCountAcrossTrackSurface, elementsCountAlongTrackSurface, nx, nd1, nd2)
