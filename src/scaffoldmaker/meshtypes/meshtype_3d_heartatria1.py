"""
Generates a 3-D heart atria model, suitable for attachment to the
3-D Heart Ventricles with Base 1.
"""

from __future__ import division

import copy
import math

from opencmiss.maths.vectorops import add, cross, dot, magnitude, mult, normalize, sub
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field, FieldGroup
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName, \
    findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.derivativemoothing import DerivativeScalingMode, DerivativeSmoothing
from scaffoldmaker.utils.eft_utils import createEftElementSurfaceLayer, remapEftLocalNodes, remapEftNodeValueLabel, \
    scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import getApproximateEllipsePerimeter, getCircleProjectionAxes, \
    getEllipseAngleFromVector, getEllipseArcLength, getEllipseRadiansToX, updateEllipseAngleByArcLength, \
    createCirclePoints
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, calculate_surface_axes


class MeshType_3d_heartatria1(Scaffold_base):
    '''
    3-D heart atria model, suitable for attachment to the 3-D Heart Ventricles with Base 1.
    '''

    lpvOstiumDefaultScaffoldPackages = {
        'LPV Human 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 2,
                'Number of elements across common' : 2,
                'Number of elements around ostium' : 12,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.19,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.16,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.11,
                'Vessel wall thickness' : 0.009,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 0.0,
                'Vessel angle 2 degrees' : 10.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'LPV Pig 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 1,
                'Number of elements across common' : 2,
                'Number of elements around ostium' : 12,
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
                'Vessel angle 1 spread degrees' : 0.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'LPV Rat 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 3,
                'Number of elements across common' : 2,
                'Number of elements around ostium' : 10,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.16,
                'Ostium length' : 0.08,
                'Ostium wall thickness' : 0.04,
                'Ostium inter-vessel distance' : 0.12,
                'Ostium inter-vessel height' : 0.01,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.5,
                'Vessel inner diameter' : 0.09,
                'Vessel wall thickness' : 0.01,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 0.0,
                'Vessel angle 2 degrees' : 0.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'LPV Rat 2' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 3,
                'Number of elements across common' : 2,
                'Number of elements around ostium' : 10,
                'Number of elements along' : 2,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.16,
                'Ostium length' : 0.1,
                'Ostium wall thickness' : 0.04,
                'Ostium inter-vessel distance' : 0.12,
                'Ostium inter-vessel height' : 0.01,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.5,
                'Vessel inner diameter' : 0.09,
                'Vessel wall thickness' : 0.01,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 0.0,
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
                'Number of elements around ostium' : 12,
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
                'Vessel angle 1 spread degrees' : 0.0,
                'Vessel angle 2 degrees' : -10.0,
                'Use linear through vessel wall' : True,
                }
            } ),
        'RPV Pig 1' : ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings' : {
                'Number of vessels' : 1,
                'Number of elements across common' : 2,
                'Number of elements around ostium' : 12,
                'Number of elements along' : 1,
                'Number of elements through wall' : 1,
                'Unit scale' : 1.0,
                'Outlet' : False,
                'Ostium diameter' : 0.26,
                'Ostium length' : 0.04,
                'Ostium wall thickness' : 0.02,
                'Ostium inter-vessel distance' : 0.16,
                'Ostium inter-vessel height' : 0.0,
                'Use linear through ostium wall' : False,
                'Vessel end length factor' : 1.0,
                'Vessel inner diameter' : 0.19,
                'Vessel wall thickness' : 0.013,
                'Vessel angle 1 degrees' : 0.0,
                'Vessel angle 1 spread degrees' : 0.0,
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
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        isHuman = 'Human' in parameterSetName
        isMouse = 'Mouse' in parameterSetName
        isPig = 'Pig' in parameterSetName
        isRat = 'Rat' in parameterSetName
        if isPig:
            lpvOstium = cls.lpvOstiumDefaultScaffoldPackages['LPV Pig 1']
            rpvOstium = cls.rpvOstiumDefaultScaffoldPackages['RPV Pig 1']
        elif isMouse or isRat:
            lpvOstium = cls.lpvOstiumDefaultScaffoldPackages['LPV Rat 1']
            rpvOstium = cls.rpvOstiumDefaultScaffoldPackages['RPV Human 1']
        else:
            lpvOstium = cls.lpvOstiumDefaultScaffoldPackages['LPV Human 1']
            rpvOstium = cls.rpvOstiumDefaultScaffoldPackages['RPV Human 1']
        # to avoid circular include, repeat inter-dependent ventriclesbase options here:
        ventriclesbaseOptions = {}
        if isMouse or isRat:
            ventriclesbaseOptions['LV outlet inner diameter'] = 0.21
            ventriclesbaseOptions['LV outlet wall thickness'] =  0.02
        else:
            ventriclesbaseOptions['LV outlet inner diameter'] = 0.28
            ventriclesbaseOptions['LV outlet wall thickness'] = 0.022
        options = {}
        options['Number of elements along atrial appendages'] = 2
        options['Number of elements along vena cava inlet'] = 3
        options['Number of elements around atrial septum'] = 3
        options['Number of elements around left atrium free wall'] = 8
        options['Number of elements around right atrium free wall'] = 6
        options['Number of elements over atria'] = 8
        options['Number of elements radial pulmonary vein annuli'] = 2
        options['Unit scale'] = 1.0
        options['Atria base inner major axis length'] = 0.47
        options['Atria base inner minor axis length'] = 0.41
        options['Atria major axis rotation degrees'] = 40.0
        options['Atria outer height'] = 0.45
        options['Atrial septum height'] = 0.25
        options['Atrial septum length'] = 0.25
        options['Atrial septum thickness'] = 0.07
        options['Atrial vestibule outer height'] = 0.05
        options['Fossa ovalis height'] = 0.1
        options['Fossa ovalis length'] = 0.14
        options['Fossa ovalis thickness'] = 0.07
        options['Fossa ovalis midpoint height'] = 0.15
        options['Left atrium venous free wall thickness'] = 0.02
        options['Right atrium venous free wall thickness'] = 0.015
        options['Crista terminalis thickness'] = 0.03
        options['Atrial base wall thickness'] = 0.06
        options['Atrial base slope degrees'] = 30.0
        options['Aorta outer plus diameter'] = ventriclesbaseOptions['LV outlet inner diameter'] + 2.0*ventriclesbaseOptions['LV outlet wall thickness']
        options['Atrial base front incline degrees'] = 15.0
        options['Atrial base back incline degrees'] = 20.0
        options['Atrial base side incline degrees'] = 20.0
        options['Atria venous anterior over'] = 0.7
        options['Atria venous midpoint over'] = 0.41
        options['Left atrium venous midpoint left'] = 0.5
        options['Right atrium venous right'] = 0.45
        options['Left atrial appendage angle axial degrees'] = 0.0
        options['Left atrial appendage angle left degrees'] = 20.0
        options['Left atrial appendage angle up degrees'] = -60.0
        options['Left atrial appendage arc length'] = 0.3
        options['Left atrial appendage arc radius'] = 0.15
        options['Left atrial appendage base length'] = 0.3
        options['Left atrial appendage left'] = 0.85
        options['Left atrial appendage midpoint left'] = 0.45
        options['Left atrial appendage midpoint over'] = 0.95
        options['Left atrial appendage wall thickness'] = 0.025
        options['Left atrial appendage wedge angle degrees'] = 60.0
        options['Right atrial appendage angle axial degrees'] = 30.0
        options['Right atrial appendage angle left degrees'] = 40.0
        options['Right atrial appendage angle up degrees'] = -20.0
        options['Right atrial appendage arc length'] = 0.5
        options['Right atrial appendage arc radius'] = 0.17
        options['Right atrial appendage base length'] = 0.3
        options['Right atrial appendage midpoint right'] = 0.47
        options['Right atrial appendage midpoint over'] = 0.92
        options['Right atrial appendage pouch right'] = 0.9
        options['Right atrial appendage wall thickness'] = 0.025
        options['Right atrial appendage wedge angle degrees'] = 90.0
        options['Common left-right pulmonary vein ostium'] = False
        options['Left pulmonary vein ostium'] = copy.deepcopy(lpvOstium)
        options['Left pulmonary vein ostium angle degrees'] = 65.0
        options['Left pulmonary vein ostium position left'] = 0.64
        options['Left pulmonary vein ostium position over'] = 0.47
        options['Right pulmonary vein ostium'] = copy.deepcopy(rpvOstium)
        options['Right pulmonary vein ostium angle degrees'] = 80.0
        options['Right pulmonary vein ostium position left'] = 0.2
        options['Right pulmonary vein ostium position over'] = 0.4
        options['Inferior vena cava inlet position over'] = 0.18
        options['Inferior vena cava inlet position right'] = 0.25
        options['Inferior vena cava inlet angle left degrees'] = 0.0
        options['Inferior vena cava inlet angle over degrees'] = -30.0
        options['Inferior vena cava inlet derivative factor'] = 1.0
        options['Inferior vena cava inlet length'] = 0.1
        options['Inferior vena cava inlet inner diameter'] = 0.22
        options['Inferior vena cava inlet wall thickness'] = 0.015
        options['Superior vena cava inlet position over'] = 0.65
        options['Superior vena cava inlet position right'] = 0.22
        options['Superior vena cava inlet angle left degrees'] = -15.0
        options['Superior vena cava inlet angle over degrees'] = 15.0
        options['Superior vena cava inlet derivative factor'] = 1.0
        options['Superior vena cava inlet length'] = 0.1
        options['Superior vena cava inlet inner diameter'] = 0.18
        options['Superior vena cava inlet wall thickness'] = 0.015
        options['Define epicardium layer'] = False
        options['Epicardium layer minimum thickness'] = 0.01
        options['Refine'] = False
        options['Refine number of elements surface'] = 4
        options['Refine number of elements through wall'] = 1
        options['Refine number of elements through epicardium layer'] = 1
        options['Use cross derivatives'] = False

        if isMouse or isRat:
            options['Atria base inner major axis length'] = 0.32
            options['Atria base inner minor axis length'] = 0.26
            options['Atria major axis rotation degrees'] = 30.0
            options['Atria outer height'] = 0.4
            options['Atria venous anterior over'] = 0.75
            options['Atria venous midpoint over'] = 0.45
            options['Atrial base back incline degrees'] = 0.0
            options['Atrial base front incline degrees'] = 15.0
            options['Atrial base side incline degrees'] = 20.0
            options['Atrial base slope degrees'] = 20.0
            options['Atrial base wall thickness'] = 0.04
            options['Atrial septum height'] = 0.25
            options['Atrial septum length'] = 0.14
            options['Atrial septum thickness'] = 0.05
            options['Atrial vestibule outer height'] = 0.06
            options['Common left-right pulmonary vein ostium'] = True
            options['Crista terminalis thickness'] = 0.03
            options['Fossa ovalis height'] = 0.07
            options['Fossa ovalis length'] = 0.07
            options['Fossa ovalis midpoint height'] = 0.15
            options['Fossa ovalis thickness'] = options['Atrial septum thickness']
            options['Inferior vena cava inlet angle left degrees'] = 0.0
            options['Inferior vena cava inlet angle over degrees'] = -20.0
            options['Inferior vena cava inlet derivative factor'] = 2.0
            options['Inferior vena cava inlet inner diameter'] = 0.16
            options['Inferior vena cava inlet length'] = 0.2
            options['Inferior vena cava inlet position over'] = 0.2
            options['Inferior vena cava inlet position right'] = 0.28
            options['Inferior vena cava inlet wall thickness'] = 0.012
            options['Left atrial appendage angle axial degrees'] = -10.0
            options['Left atrial appendage angle left degrees'] = 0.0
            options['Left atrial appendage angle up degrees'] = -60.0
            options['Left atrial appendage arc length'] = 0.5
            options['Left atrial appendage arc radius'] = 0.35
            options['Left atrial appendage base length'] = 0.35
            options['Left atrial appendage left'] = 0.99
            options['Left atrial appendage midpoint left'] = 0.78
            options['Left atrial appendage midpoint over'] = 1.0
            options['Left atrial appendage wall thickness'] = 0.025
            options['Left atrial appendage wedge angle degrees'] = 80.0
            options['Left atrium venous free wall thickness'] = 0.04
            options['Left atrium venous midpoint left'] = 0.45
            options['Left pulmonary vein ostium angle degrees'] = 50.0
            options['Left pulmonary vein ostium position left'] = 0.3
            options['Left pulmonary vein ostium position over'] = 0.36
            options['Right atrial appendage angle axial degrees'] = 35.0
            options['Right atrial appendage angle left degrees'] = 20.0
            options['Right atrial appendage angle up degrees'] = 0.0
            options['Right atrial appendage arc length'] = 0.6
            options['Right atrial appendage arc radius'] = 0.35
            options['Right atrial appendage base length'] = 0.35
            options['Right atrial appendage midpoint right'] = 0.45
            options['Right atrial appendage midpoint over'] = 0.75
            options['Right atrial appendage pouch right'] = 0.95
            options['Right atrial appendage wall thickness'] = 0.025
            options['Right atrial appendage wedge angle degrees'] = 90.0
            options['Right atrium venous free wall thickness'] = 0.03
            options['Right atrium venous right'] = 0.5
            options['Superior vena cava inlet angle left degrees'] = 0.0
            options['Superior vena cava inlet angle over degrees'] = -5.0
            options['Superior vena cava inlet derivative factor'] = 2.0
            options['Superior vena cava inlet inner diameter'] = 0.15
            options['Superior vena cava inlet length'] = 0.2
            options['Superior vena cava inlet position over'] = 0.62
            options['Superior vena cava inlet position right'] = 0.25
            options['Superior vena cava inlet wall thickness'] = 0.012
        elif isPig:
            options['Atrial base side incline degrees'] = 0.0
            options['Left atrial appendage angle axial degrees'] = -10.0
            options['Left atrial appendage angle left degrees'] = 20.0
            options['Left atrial appendage angle up degrees'] = -60.0
            options['Left atrial appendage arc length'] = 0.4
            options['Left atrial appendage arc radius'] = 0.3
            options['Left atrial appendage base length'] = 0.35
            options['Left atrial appendage left'] = 0.9
            options['Left atrial appendage midpoint left'] = 0.55
            options['Left atrial appendage midpoint over'] = 1.0
            options['Left atrial appendage wedge angle degrees'] = 50.0
            options['Left pulmonary vein ostium angle degrees'] = 65.0
            options['Left pulmonary vein ostium position left'] = 0.67
            options['Left pulmonary vein ostium position over'] = 0.42
            options['Right atrial appendage angle axial degrees'] = 10.0
            options['Right atrial appendage angle left degrees'] = -20.0
            options['Right atrial appendage angle up degrees'] = -10.0
            options['Right atrial appendage arc length'] = 0.5
            options['Right atrial appendage arc radius'] = 0.25
            options['Right atrial appendage base length'] = 0.25
            options['Right atrial appendage midpoint right'] = 0.55
            options['Right atrial appendage wedge angle degrees'] = 60.0
            options['Right pulmonary vein ostium angle degrees'] = 80.0
            options['Right pulmonary vein ostium position left'] = 0.22
            options['Right pulmonary vein ostium position over'] = 0.33
            options['Inferior vena cava inlet position over'] = 0.22
            options['Inferior vena cava inlet position right'] = 0.25
            options['Inferior vena cava inlet angle left degrees'] = 0.0
            options['Inferior vena cava inlet angle over degrees'] = 0.0
            options['Inferior vena cava inlet derivative factor'] = 0.5
            options['Superior vena cava inlet position over'] = 0.6
            options['Superior vena cava inlet position right'] = 0.25
            options['Superior vena cava inlet angle left degrees'] = 0.0
            options['Superior vena cava inlet angle over degrees'] = -10.0
        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along atrial appendages',
            'Number of elements along vena cava inlet',
            'Number of elements around atrial septum',
            'Number of elements around left atrium free wall',
            'Number of elements around right atrium free wall',
            'Number of elements over atria',
            'Number of elements radial pulmonary vein annuli',
            'Unit scale',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria major axis rotation degrees',
            'Atria outer height',
            'Atrial septum height',
            'Atrial septum length',
            'Atrial septum thickness',
            'Atrial vestibule outer height',
            'Fossa ovalis height',
            'Fossa ovalis length',
            'Fossa ovalis midpoint height',
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
            'Atria venous anterior over',
            'Atria venous midpoint over',
            'Left atrium venous midpoint left',
            'Right atrium venous right',
            'Left atrial appendage angle axial degrees',
            'Left atrial appendage angle left degrees',
            'Left atrial appendage angle up degrees',
            'Left atrial appendage arc length',
            'Left atrial appendage arc radius',
            'Left atrial appendage base length',
            'Left atrial appendage left',
            'Left atrial appendage midpoint left',
            'Left atrial appendage midpoint over',
            'Left atrial appendage wall thickness',
            'Left atrial appendage wedge angle degrees',
            'Right atrial appendage angle axial degrees',
            'Right atrial appendage angle left degrees',
            'Right atrial appendage angle up degrees',
            'Right atrial appendage arc length',
            'Right atrial appendage arc radius',
            'Right atrial appendage base length',
            'Right atrial appendage midpoint right',
            'Right atrial appendage midpoint over',
            'Right atrial appendage pouch right',
            'Right atrial appendage wall thickness',
            'Right atrial appendage wedge angle degrees',
            'Common left-right pulmonary vein ostium',
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
            'Define epicardium layer',
            'Epicardium layer minimum thickness',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall',
            'Refine number of elements through epicardium layer'
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

    @classmethod
    def checkOptions(cls, options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        for key in  [
            'Number of elements along atrial appendages',
            'Number of elements along vena cava inlet',
            'Number of elements radial pulmonary vein annuli'
            ]:
            if options[key] < 1:
                options[key] = 1
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
            if options[key] <= 6:
                options[key] = 6
            elif options[key] >= 10:
                options[key] = 10
            else:
                options[key] = 8
        for key in [
            'Unit scale',
            'Atria base inner major axis length',
            'Atria base inner minor axis length',
            'Atria outer height',
            'Atrial septum height',
            'Atrial septum length',
            'Atrial septum thickness',
            'Atrial vestibule outer height',
            'Fossa ovalis height',
            'Fossa ovalis length',
            'Fossa ovalis thickness',
            'Fossa ovalis midpoint height',
            'Left atrium venous free wall thickness',
            'Right atrium venous free wall thickness',
            'Crista terminalis thickness',
            'Atrial base wall thickness',
            'Atrial base slope degrees',
            'Left atrial appendage arc length',
            'Left atrial appendage base length',
            'Left atrial appendage wall thickness',
            'Right atrial appendage arc length',
            'Right atrial appendage base length',
            'Right atrial appendage wall thickness',
            'Inferior vena cava inlet derivative factor',
            'Inferior vena cava inlet length',
            'Inferior vena cava inlet inner diameter',
            'Inferior vena cava inlet wall thickness',
            'Superior vena cava inlet derivative factor',
            'Superior vena cava inlet length',
            'Superior vena cava inlet inner diameter',
            'Superior vena cava inlet wall thickness',
            'Epicardium layer minimum thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        for key in [
            'Atria venous anterior over',
            'Atria venous midpoint over',
            'Left atrium venous midpoint left',
            'Right atrium venous right',
            'Left atrial appendage left',
            'Left atrial appendage midpoint left',
            'Left atrial appendage midpoint over',
            'Right atrial appendage midpoint right',
            'Right atrial appendage midpoint over',
            'Left pulmonary vein ostium position left',
            'Left pulmonary vein ostium position over',
            'Right pulmonary vein ostium position left',
            'Right pulmonary vein ostium position over',
            'Inferior vena cava inlet position over',
            'Inferior vena cava inlet position right',
            'Superior vena cava inlet position over',
            'Superior vena cava inlet position right']:
            if options[key] < 0.0:
                options[key] = 0.0
            elif options[key] > 1.0:
                options[key] = 1.0
        for key in [
            'Right atrial appendage pouch right']:
            if options[key] < 0.0:
                options[key] = 0.0
            elif options[key] > 2.0:
                options[key] = 2.0
        if options['Aorta outer plus diameter'] < options['Atrial septum thickness']:
            options['Aorta outer plus diameter'] = options['Atrial septum thickness']
        for key in [
            'Atria major axis rotation degrees']:
            if options[key] < -75.0:
                options[key] = -75.0
            elif options[key] > 75.0:
                options[key] = 75.0
        for key in [
            'Left atrial appendage arc radius',
            'Right atrial appendage arc radius']:
            if options[key] < 0.1*math.fabs(options['Atria base inner minor axis length']):
                options[key] = 0.1*math.fabs(options['Atria base inner minor axis length'])
                dependentChanges = True
            elif options[key] > 1000.0*math.fabs(options['Atria base inner major axis length']):
                options[key] = 1000.0*math.fabs(options['Atria base inner major axis length'])
                dependentChanges = True
        for key in [
            'Left atrial appendage wedge angle degrees',
            'Right atrial appendage wedge angle degrees']:
            if options[key] < 5.0:
                options[key] = 5.0
            elif options[key] > 150.0:
                options[key] = 150.0
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        cls.updateSubScaffoldOptions(options)
        return dependentChanges

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        '''
        Update lpv, rpv ostium sub-scaffold options which depend on parent options.
        '''
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountOverAtria = options['Number of elements over atria']
        unitScale = options['Unit scale']
        laVenousFreeWallThickness = options['Left atrium venous free wall thickness']  # not scaled by unitScale
        commonLeftRightPvOstium = options['Common left-right pulmonary vein ostium']
        lpvOstium = options['Left pulmonary vein ostium']
        rpvOstium = options['Right pulmonary vein ostium']
        elementsCountAroundLpvOstium, elementsCountAroundRpvOstium = \
            getLeftAtriumPulmonaryVeinOstiaElementsCounts(elementsCountAroundLeftAtriumFreeWall, elementsCountOverAtria, commonLeftRightPvOstium)
        lpvOstiumSettings = lpvOstium.getScaffoldSettings()
        lpvOstiumSettings['Number of elements around ostium'] = elementsCountAroundLpvOstium
        lpvOstiumSettings['Unit scale'] = unitScale
        lpvOstiumSettings['Ostium wall thickness'] = laVenousFreeWallThickness
        lpvOstiumSettings['Outlet'] = False
        lpvOstiumSettings['Use linear through ostium wall'] = False
        rpvOstiumSettings = rpvOstium.getScaffoldSettings()
        rpvOstiumSettings['Number of elements around ostium'] = elementsCountAroundRpvOstium
        rpvOstiumSettings['Unit scale'] = unitScale
        rpvOstiumSettings['Ostium wall thickness'] = laVenousFreeWallThickness
        rpvOstiumSettings['Outlet'] = False
        rpvOstiumSettings['Use linear through ostium wall'] = False


    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        cls.updateSubScaffoldOptions(options)

        elementsCountAlongAtrialAppendages = options['Number of elements along atrial appendages']
        elementsCountAlongVCInlet = options['Number of elements along vena cava inlet']
        elementsCountAroundAtrialSeptum = options['Number of elements around atrial septum']
        elementsCountAroundLeftAtriumFreeWall = options['Number of elements around left atrium free wall']
        elementsCountAroundLeftAtrium = elementsCountAroundLeftAtriumFreeWall + elementsCountAroundAtrialSeptum
        elementsCountAroundRightAtriumFreeWall = options['Number of elements around right atrium free wall']
        elementsCountOverAtria = options['Number of elements over atria']
        elementsCountRadialPVAnnuli = options['Number of elements radial pulmonary vein annuli']
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
        aVestibuleOuterHeight = unitScale*options['Atrial vestibule outer height']
        foMidpointZ = unitScale*options['Fossa ovalis midpoint height']
        foMagZ = unitScale*0.5*options['Fossa ovalis height']
        foMagY = unitScale*0.5*options['Fossa ovalis length']
        foThickness = unitScale*options['Fossa ovalis thickness']
        laVenousFreeWallThickness = unitScale*options['Left atrium venous free wall thickness']
        raVenousFreeWallThickness = unitScale*options['Right atrium venous free wall thickness']
        cristaTerminalisThickness = unitScale*options['Crista terminalis thickness']
        aBaseWallThickness = unitScale*options['Atrial base wall thickness']
        aBaseSlopeRadians = math.radians(options['Atrial base slope degrees'])
        aVenousAnteriorOver = options['Atria venous anterior over']
        aVenousMidpointOver = options['Atria venous midpoint over']
        laVenousMidpointLeft = options['Left atrium venous midpoint left']
        raVenousRight = options['Right atrium venous right']
        laaAngleAxialRadians = math.radians(options['Left atrial appendage angle axial degrees'])
        laaAngleLeftRadians = math.radians(options['Left atrial appendage angle left degrees'])
        laaAngleUpradians = math.radians(options['Left atrial appendage angle up degrees'])
        laaArcLength = unitScale*options['Left atrial appendage arc length']
        laaArcRadius = unitScale*options['Left atrial appendage arc radius']
        laaBaseLength = unitScale*options['Left atrial appendage base length']
        laaLeft = options['Left atrial appendage left']
        laaMidpointLeft = options['Left atrial appendage midpoint left']
        laaMidpointOver = options['Left atrial appendage midpoint over']
        laaWallThickness = unitScale*options['Left atrial appendage wall thickness']
        laaWedgeAngleRadians = math.radians(options['Left atrial appendage wedge angle degrees'])
        raaAngleAxialRadians = math.radians(options['Right atrial appendage angle axial degrees'])
        raaAngleLeftRadians = math.radians(options['Right atrial appendage angle left degrees'])
        raaAngleUpradians = math.radians(options['Right atrial appendage angle up degrees'])
        raaArcLength = unitScale*options['Right atrial appendage arc length']
        raaArcRadius = unitScale*options['Right atrial appendage arc radius']
        raaBaseLength = unitScale*options['Right atrial appendage base length']
        raaMidpointRight = options['Right atrial appendage midpoint right']
        raaMidpointOver = options['Right atrial appendage midpoint over']
        raaPouchRight = options['Right atrial appendage pouch right']
        raaWallThickness = unitScale*options['Right atrial appendage wall thickness']
        raaWedgeAngleRadians = math.radians(options['Right atrial appendage wedge angle degrees'])
        commonLeftRightPvOstium = options['Common left-right pulmonary vein ostium']
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
        defineEpicardiumLayer = options['Define epicardium layer']
        epicardiumLayerMinimumThickness = unitScale*options['Epicardium layer minimum thickness']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        heartGroup = AnnotationGroup(region, get_heart_term("heart"))
        lamGroup = AnnotationGroup(region, get_heart_term("left atrium myocardium"))
        ramGroup = AnnotationGroup(region, get_heart_term("right atrium myocardium"))
        aSeptumGroup = AnnotationGroup(region, get_heart_term("interatrial septum"))
        fossaGroup = AnnotationGroup(region, get_heart_term("fossa ovalis"))
        laaGroup = AnnotationGroup(region, get_heart_term("left auricle"))
        raaGroup = AnnotationGroup(region, get_heart_term("right auricle"))
        annotationGroups = [heartGroup, lamGroup, ramGroup, aSeptumGroup, fossaGroup, laaGroup, raaGroup]
        if defineEpicardiumLayer:
            epicardiumGroup = AnnotationGroup(region, get_heart_term("epicardium"))
            annotationGroups.append(epicardiumGroup)

        lpvOstiumSettings = lpvOstium.getScaffoldSettings()
        lpvCount = lpvOstiumSettings['Number of vessels']
        if commonLeftRightPvOstium:
            # use only lpv:
            if lpvCount == 1:
                pvGroup = AnnotationGroup(region, get_heart_term("pulmonary vein"))
                lpvGroups = [ pvGroup ]
            else:
                lpvGroup = AnnotationGroup(region, get_heart_term("left pulmonary vein"))
                rpvGroup = AnnotationGroup(region, get_heart_term("right pulmonary vein"))
                if lpvCount == 2:
                    lpvGroups = [ lpvGroup, rpvGroup ]
                else:
                    mpvGroup = AnnotationGroup(region, get_heart_term("middle pulmonary vein"))
                    lpvGroups = [ lpvGroup, mpvGroup, rpvGroup ]
            annotationGroups += lpvGroups
        else:  # separate left and right pulmonary vein ostia
            lpvGroup = AnnotationGroup(region, get_heart_term("left pulmonary vein"))
            lpvGroups = [ lpvGroup ]*lpvCount
            annotationGroups.append(lpvGroup)
            rpvOstiumSettings = rpvOstium.getScaffoldSettings()
            rpvCount = rpvOstiumSettings['Number of vessels']
            rpvGroup = AnnotationGroup(region, get_heart_term("right pulmonary vein" ))
            rpvGroups = [ rpvGroup ]*rpvCount
            annotationGroups.append(rpvGroup)

        ivcInletGroup = AnnotationGroup(region, get_heart_term("inferior vena cava inlet"))
        svcInletGroup = AnnotationGroup(region, get_heart_term("superior vena cava inlet"))
        ivcGroup = AnnotationGroup(region, get_heart_term("inferior vena cava"))
        svcGroup = AnnotationGroup(region, get_heart_term("superior vena cava"))
        laeVenousMidpointGroup =\
            AnnotationGroup(region, get_heart_term("left atrium epicardium venous midpoint"), isMarker=True)
        raeVenousMidpointGroup =\
            AnnotationGroup(region, get_heart_term("right atrium epicardium venous midpoint"), isMarker=True)
        # av boundary nodes are put in left and right fibrous ring groups so they can be found by heart1
        lFibrousRingGroup = AnnotationGroup(region, get_heart_term("left fibrous ring"))
        rFibrousRingGroup = AnnotationGroup(region, get_heart_term("right fibrous ring"))
        annotationGroups += [laeVenousMidpointGroup, ivcGroup, ivcInletGroup, raeVenousMidpointGroup,
                             svcGroup, svcInletGroup, lFibrousRingGroup, rFibrousRingGroup]

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

        mesh = fm.findMeshByDimension(3)
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)

        heartMeshGroup = heartGroup.getMeshGroup(mesh)
        lamMeshGroup = lamGroup.getMeshGroup(mesh)
        ramMeshGroup = ramGroup.getMeshGroup(mesh)
        aSeptumMeshGroup = aSeptumGroup.getMeshGroup(mesh)
        fossaMeshGroup = fossaGroup.getMeshGroup(mesh)
        laaMeshGroup = laaGroup.getMeshGroup(mesh)
        raaMeshGroup = raaGroup.getMeshGroup(mesh)
        lpvMeshGroups = [ inletGroup.getMeshGroup(mesh) for inletGroup in lpvGroups ]
        if not commonLeftRightPvOstium:
            rpvMeshGroups = [ inletGroup.getMeshGroup(mesh) for inletGroup in rpvGroups ]
        ivcInletMeshGroup = ivcInletGroup.getMeshGroup(mesh)
        svcInletMeshGroup = svcInletGroup.getMeshGroup(mesh)
        ivcMeshGroup = ivcGroup.getMeshGroup(mesh)
        svcMeshGroup = svcGroup.getMeshGroup(mesh)

        # get elements count over atria, around left and right free wall base, and around ostia
        # note elementsCountOverAtriaCoronarySinus is assumed to be 1
        elementsCountOverAtriaCoronarySinus, \
        elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior, \
        elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior \
            = getOverAtriaElementsCounts(elementsCountOverAtria)
        elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftAtrialAppendageBase, elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV \
            = getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall)
        elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightAtrialAppendagePlainBase, \
            elementsCountAroundRightAtrialAppendagePouchBase, elementsCountAroundRightAtriumAorta \
            = getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall)
        elementsCountAroundLpvOstium, elementsCountAroundRpvOstium = \
            getLeftAtriumPulmonaryVeinOstiaElementsCounts(elementsCountAroundLeftAtriumFreeWall, elementsCountOverAtria, commonLeftRightPvOstium)
        if commonLeftRightPvOstium:
            elementsCountOverSideLeftAtriumLPV = elementsCountAroundLeftAtriumLPV + 1
        else:
            elementsCountOverSideLeftAtriumLPV = elementsCountAroundLeftAtriumLPV

        # GRC fudge factors:
        aOuterSeptumHeight = 1.2*aSeptumHeight
        iaGrooveDerivative = 1.0*aSeptumThickness

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
                laaLeft, laVenousMidpointLeft, raVenousRight, raaPouchRight, elementsCountAroundTrackSurface)

        laTrackSurface = getAtriumTrackSurface(elementsCountAroundTrackSurface, elementsCountAcrossTrackSurface,
            ltBaseOuterx, ltBaseOuterd1, ltBaseOuterd2, aSeptumBaseCentre, aOuterHeight, aOuterSeptumHeight, iaGrooveDerivative)
        raTrackSurface = laTrackSurface.createMirrorX()

        # need to create pulmonary vein ostia early because other derivatives are smoothed to fit them

        # create left pulmonary vein ostium
        lpvOstiumPosition = laTrackSurface.createPositionProportion(lpvOstiumPositionOver, lpvOstiumPositionLeft)
        # get absolute direction on surface corresponding to chosen angle
        cx, cd1, cd2 = laTrackSurface.evaluateCoordinates(lpvOstiumPosition, derivatives = True)
        td1, td2, td3 = calculate_surface_axes(cd1, cd2, [ 0.0, 1.0, 0.0 ])
        zAngleRadians = math.atan2(td1[0], -td2[0])
        #print('zAngleRadians',zAngleRadians)
        cosAngle = math.cos(zAngleRadians + lpvOstiumAngleRadians)
        sinAngle = math.sin(zAngleRadians + lpvOstiumAngleRadians)
        lpvOstiumDirection = [ (cosAngle*-td2[c] + sinAngle*td1[c]) for c in range(3) ]
        nodeIdentifier, elementIdentifier, (lpvox, lpvod1, lpvod2, lpvod3, lpvoNodeId, lpvoPositions) = \
            generateOstiumMesh(region, lpvOstiumSettings, laTrackSurface, lpvOstiumPosition, lpvOstiumDirection, nodeIdentifier, elementIdentifier,
                               vesselMeshGroups=[ [ heartMeshGroup, meshGroup ] for meshGroup in lpvMeshGroups ], ostiumMeshGroups=[ heartMeshGroup, lamMeshGroup ])

        if not commonLeftRightPvOstium:
            # create right pulmonary vein ostium
            rpvOstiumPosition = laTrackSurface.createPositionProportion(rpvOstiumPositionOver, rpvOstiumPositionLeft)
            # get absolute direction on surface corresponding to chosen angle
            cx, cd1, cd2 = laTrackSurface.evaluateCoordinates(rpvOstiumPosition, derivatives = True)
            td1, td2, td3 = calculate_surface_axes(cd1, cd2, [ 0.0, 1.0, 0.0 ])
            zAngleRadians = math.atan2(td1[0], -td2[0])
            #print('zAngleRadians',zAngleRadians)
            cosAngle = math.cos(zAngleRadians + rpvOstiumAngleRadians)
            sinAngle = math.sin(zAngleRadians + rpvOstiumAngleRadians)
            rpvOstiumDirection = [ (cosAngle*-td2[c] + sinAngle*td1[c]) for c in range(3) ]
            nodeIdentifier, elementIdentifier, (rpvox, rpvod1, rpvod2, rpvod3, rpvoNodeId, rpvoPositions) = \
                generateOstiumMesh(region, rpvOstiumSettings, laTrackSurface, rpvOstiumPosition, rpvOstiumDirection, nodeIdentifier, elementIdentifier,
                                   vesselMeshGroups=[ [ heartMeshGroup, meshGroup ] for meshGroup in rpvMeshGroups ], ostiumMeshGroups=[ heartMeshGroup, lamMeshGroup ])

        # get points over interatrial septum on exterior groove
        agn1Mid = elementsCountOverRightAtriumNonVenousAnterior + elementsCountOverRightAtriumVenous//2
        # 1. go up vestibule height on posterior
        nx = laTrackSurface.nx [:laTrackSurface.elementsCount1 + 1]
        nd = laTrackSurface.nd1[:laTrackSurface.elementsCount1 + 1]
        csx, csd2, e, xi = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd, aVestibuleOuterHeight)
        lagcsProportion = (e + xi)/laTrackSurface.elementsCount1
        agLength = sum(interp.getCubicHermiteArcLength(nx[e], nd[e], nx[e + 1], nd[e + 1]) for e in range(laTrackSurface.elementsCount1))
        agx  = [labx [1][0]]
        agd1 = [labd1[1][0]]
        agd2 = [labd2[1][0]]
        agd3 = [labd3[1][0]]
        # add lengths over groove from anterior to posterior, intersecting aVenousAnteriorOver, aVenousMidpointOver, lagcsProportion
        count1 = elementsCountOverLeftAtriumNonVenousAnterior
        ragProportionLengths1 = [(1.0 - aVenousAnteriorOver) / count1] * count1
        count2 = elementsCountOverLeftAtriumVenous // 2
        ragProportionLengths2 = [(aVenousAnteriorOver - aVenousMidpointOver) / count2] * count2
        count3 = elementsCountOverLeftAtriumVenous // 2 + elementsCountOverLeftAtriumNonVenousPosterior \
                 - elementsCountOverAtriaCoronarySinus
        ragProportionLengths3 = [(aVenousMidpointOver - lagcsProportion) / count3] * count3
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
            x, d1t, d2t = raTrackSurface.evaluateCoordinates(trackPosition, derivatives = True)
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
            agd2.append(vector.setMagnitude(d1t, agLength*0.5*(ragProportionLengths[e] + ragProportionLengths[e + 1])))
            agd3.append(normalize(cross(d1t, d2t)))
            if e in (0, elementsCountOverAtria - 2):
                # old calculation used for vestibule top
                agd1.append([ -(f1*d1s + f3*d1f), 0.0, 0.0 ])
                for c in range(3):
                    agd3[-1][c] *= raVenousFreeWallThickness  # kludge
            else:
                # deep in interatrial groove
                # make d1 unit tangents on raTrackSurface, normal to d2, scale with d3 later when making septum nodes asd1
                agd1.append(normalize(d2t))
            # store unit normal to raTrackSurface in d3, scale with d1 later
        ragProportions.append(1.0)
        laVenousLimitPosterior = ragProportionLengths[-1] + ragProportionLengths[-2]
        # Get heights of elements on aorta up interatrial groove, for building venous limit curves
        aoHeight1 = ragProportionLengths[0]*agLength
        # fix vestibule top d2 magnitude
        agd2[-1] = vector.setMagnitude(agd2[-1], aVestibuleOuterHeight)
        # smooth d2 up to vestibule top
        agd2 = interp.smoothCubicHermiteDerivativesLine(agx, agd2, fixAllDirections = True, fixEndDerivative = True)  # , magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        # GRC , magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        # reverse derivative on posterior vestibule top
        agd1[-1] = [ -d for d in agd1[-1] ]
        agd2[-1] = [ -d for d in agd2[-1] ]
        # add posterior crux point
        agx .append(labx [1][elementsCountAroundLeftAtriumFreeWall])
        agd1.append(labd1[1][elementsCountAroundLeftAtriumFreeWall])
        agd2.append(vector.setMagnitude(labd2[1][elementsCountAroundLeftAtriumFreeWall], aVestibuleOuterHeight))
        agd3.append(labd3[1][elementsCountAroundLeftAtriumFreeWall])
        # copy derivatives to labd2[1], rabd2[1]
        rabd2[1][elementsCountAroundRightAtriumFreeWall] = labd2[1][0] = agd2[0]
        rabd2[1][0] = labd2[1][elementsCountAroundLeftAtriumFreeWall] = agd2[-1]

        # start getting points on interatrial septum next to vestibule height, then septum "arch"
        # need these to get fossa angles
        halffoThickness = 0.5*foThickness
        halfaSeptumThickness = 0.5*aSeptumThickness
        xia = 0.35  # GRC fudge factor
        coronarySinusHeightAnterior = (1.0 - xia)*aVestibuleOuterHeight + xia*foMidpointZ
        aSeptumPosteriorY = aSeptumBaseCentre[1] - 0.5*aSeptumLength
        aSeptumAnteriorY = aSeptumBaseCentre[1] + 0.5*aSeptumLength
        x1 = [ 0.0, aSeptumPosteriorY, aVestibuleOuterHeight ]
        d1 = [ 0.0, aSeptumLength, 0.0 ]
        x2 = [ 0.0, aSeptumAnteriorY, coronarySinusHeightAnterior ]
        d2 = interp.interpolateHermiteLagrangeDerivative(x1, d1, x2, 1.0)  # GRC was d1
        asx, asd1 = interp.sampleCubicHermiteCurves([ x1, x2 ], [ d1, d2 ], elementsCountAroundAtrialSeptum, arcLengthDerivatives = True)[0:2]

        # get fossa ovalis points at centre and around
        elementsCountAroundFossa = elementsCountOverAtria + elementsCountAroundAtrialSeptum - 2
        fossaPerimeterLength = getApproximateEllipsePerimeter(foMagY, foMagZ)
        estElementSizeAroundFossa = fossaPerimeterLength/elementsCountAroundFossa
        fossaInnerDerivativeRatio = 1.0  # GRC fudge factor
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

        # complete getting points on interatrial septum at vestibule top, then septum "arch"
        archMagY = 0.5*aSeptumLength
        archMagZ = aSeptumHeight - foMidpointZ
        halfArchEllipsePerimeter = 0.5*getApproximateEllipsePerimeter(archMagY, archMagZ)
        archLength = (2.0*foMidpointZ - aVestibuleOuterHeight - coronarySinusHeightAnterior) + halfArchEllipsePerimeter
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
            # now calculate interatrial groove d1, d3 to work with the above septum collapsed elements
            # raTrackSurface tangent is currently in agd1[ng]
            asd = [-d for d in add(d1, d2)]
            agd_norm = normalize(cross(cross(agd2[ng], agd1[ng]), agd2[ng]))
            agd = mult(agd_norm, magnitude(asd))
            agd = interp.smoothCubicHermiteDerivativesLine(
                [asx[-1], agx[ng]], [asd, agd], fixStartDerivative=True, fixEndDirection=True)[1]
            mag_agd = magnitude(agd)
            if mag_agd < raVenousFreeWallThickness:
                agd = mult(agd_norm, raVenousFreeWallThickness)
            agd1[ng] = [-agd[0], 0.0, 0.0]
            agd3[ng] = [0.0, agd[1], agd[2]]
            if ng not in (2, elementsCountOverAtria // 2, elementsCountOverArch):
                # add a component of agd2 to agd1 ease building adjacent annuli
                sign = (-1.0 if (ng > (elementsCountOverAtria // 2)) else 1.0)
                scale = sign * magnitude(agd1[ng]) / magnitude(agd2[ng])
                agd1[ng] = add(agd1[ng], mult(agd2[ng], scale))

        asd2s = interp.smoothCubicHermiteDerivativesLine(
            asx[elementsCountAroundAtrialSeptum:], asd2[elementsCountAroundAtrialSeptum:],
            fixAllDirections = True, fixStartDerivative=True, fixEndDerivative=True,
            magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        asd2 = asd2[:elementsCountAroundAtrialSeptum] + asd2s
        #asd2 = interp.smoothCubicHermiteDerivativesLoop(asx, asd2, fixAllDirections = True, magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
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

        # get points on external vestibule top, around entire free wall of both atria
        # not all of these will become nodes, but they are always used to set base derivatives
        # left atrium
        lavtx  = [ agx [1] ]
        lavtd1 = [ agd1[1] ]
        lavtd2 = [ agd2[1] ]
        lavtd3 = [ agd3[1] ]
        lan1Aorta = elementsCountAroundLeftAtriumAorta
        lan1Mid = elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftAtrialAppendageBase
        lan1MidVenous = lan1Mid + elementsCountAroundLeftAtriumLPV
        lavtProportions = [ [ 1.0 - ragProportions[1], 0.0 ] ]
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
            trackDistance1 = aoHeight1 if onAorta else aVestibuleOuterHeight
            position = laTrackSurface.trackVector(startPosition, direction, trackDistance1)
            lavtProportions.append(laTrackSurface.getProportion(position))
            x, d1, d2 = laTrackSurface.evaluateCoordinates(position, derivatives = True)
            ax1, ax2, ax3 = calculate_surface_axes(d1, d2, direction)
            lavtx .append(x)
            lavtd1.append(vector.setMagnitude(ax2, -vector.magnitude(labd1[1][n1])))
            lavtd2.append(vector.setMagnitude(ax1, trackDistance1))
            # if not commonLeftRightPvOstium and \
            #         (n1 == (elementsCountAroundLeftAtriumFreeWall - elementsCountAroundLeftAtriumRPV - 1)):
            #     # derivative needs a tweak for LPV annulus
            #     lavtd2[-1] = sub(lavtd2[-1], vector.setMagnitude(lavtd1[-1], trackDistance1))
            lavtd3.append(vector.setMagnitude(ax3, laVenousFreeWallThickness))
            # fix d2 on outer base
            labd2[1][n1] = vector.setMagnitude(labd2[1][n1], trackDistance1)
        # add end points and smooth d1
        lavtx .append(agx [-2])
        lavtd1.append(agd1[-2])
        lavtd2.append(agd2[-2])
        lavtd3.append(agd3[-2])
        lavtProportions.append( [ 1.0 - ragProportions[-2], 0.0 ] )
        lavtd1 = interp.smoothCubicHermiteDerivativesLine(lavtx, lavtd1, fixAllDirections = True, fixStartDerivative = True, fixEndDerivative = True)
        # get inner points
        lavtx  = [ [agx [0]], lavtx  ]
        lavtd1 = [ [agd1[0]], lavtd1 ]
        lavtd2 = [ [agd2[0]], lavtd2 ]
        lavtd3 = [ [agd3[0]], lavtd3 ]
        for n1 in range(1, elementsCountAroundLeftAtriumFreeWall):
            x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(lavtx[1], lavtd1[1], lavtd2[1], n1, -laVenousFreeWallThickness)
            # do same upwards to get proper value of d2
            nx  = [ labx [1][n1], lavtx [1][n1] ]
            nd1 = [ labd1[1][n1], lavtd1[1][n1] ]
            nd2 = [ labd2[1][n1], lavtd2[1][n1] ]
            _, d2, _, _ = interp.projectHermiteCurvesThroughWall(nx, nd2, [ [ -d for d in d1 ] for d1 in nd1 ], 1, -laVenousFreeWallThickness)
            lavtx [0].append(x)
            lavtd1[0].append(d1)
            lavtd2[0].append(d2)
            lavtd3[0].append(d3)
            # fix d2 on inner base
            labd2[0][n1] = interp.interpolateLagrangeHermiteDerivative(labx[0][n1], x, d2, 0.0)
        lavtx [0].append(agx [-1])
        lavtd1[0].append(agd1[-1])
        lavtd2[0].append(agd2[-1])
        lavtd3[0].append(agd3[-1])
        # right atrium
        ravtx  = [ agx [-2] ]
        ravtd1 = [ agd1[-2] ]
        ravtd2 = [ agd2[-2] ]
        ravtd3 = [ agd3[-2] ]
        ravtProportions = [ [ ragProportions[-2], 0.0 ] ]
        rabProportions = [ [ 1.0, 0.0 ] ]
        ran1Ctp = elementsCountAroundRightAtriumPosteriorVenous
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
            rabProportions.append(raTrackSurface.getProportion(startPosition))
            onAorta = n1 == (elementsCountAroundRightAtriumFreeWall - 1)
            direction = [ 0.0, 0.0, 1.0 ] if onAorta else vector.normalise(rabd2[1][n1])
            trackDistance1 = aoHeight1 if onAorta else aVestibuleOuterHeight
            position = raTrackSurface.trackVector(startPosition, direction, trackDistance1)
            ravtProportions.append(raTrackSurface.getProportion(position))
            x, d1, d2 = raTrackSurface.evaluateCoordinates(position, derivatives = True)
            ax1, ax2, ax3 = calculate_surface_axes(d1, d2, direction)
            ravtx .append(x)
            ravtd1.append(vector.setMagnitude(ax2, -vector.magnitude(rabd1[1][n1])))
            ravtd2.append(vector.setMagnitude(ax1, trackDistance1))
            ravtd3.append(vector.setMagnitude(ax3, raVenousFreeWallThickness))
            # fix d2 on outer base
            rabd2[1][n1] = vector.setMagnitude(rabd2[1][n1], trackDistance1)
        # add end points and smooth d1
        ravtx .append(agx [1])
        ravtd1.append(agd1[1])
        ravtd2.append(agd2[1])
        ravtd3.append(agd3[1])
        ravtProportions.append([ ragProportions[1], 0.0 ])
        rabProportions.append([ 0.0, 0.0 ])
        # get inner points
        ravtx  = [ [agx [-1]], ravtx  ]
        ravtd1 = [ [agd1[-1]], ravtd1 ]
        ravtd2 = [ [agd2[-1]], ravtd2 ]
        ravtd3 = [ [agd3[-1]], ravtd3 ]
        for n1 in range(1, elementsCountAroundRightAtriumFreeWall):
            x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(ravtx[1], ravtd1[1], ravtd2[1], n1,
                -(raVenousFreeWallThickness if (n1 < elementsCountAroundRightAtriumPosteriorVenous) else raaWallThickness))
            # do same upwards to get proper value of d2
            nx  = [ rabx [1][n1], ravtx [1][n1] ]
            nd1 = [ rabd1[1][n1], ravtd1[1][n1] ]
            nd2 = [ rabd2[1][n1], ravtd2[1][n1] ]
            _, d2, _, _ = interp.projectHermiteCurvesThroughWall(nx, nd2, [ [ -d for d in d1 ] for d1 in nd1 ], 1, -raVenousFreeWallThickness)
            ravtx [0].append(x)
            ravtd1[0].append(d1)
            ravtd2[0].append(d2)
            ravtd3[0].append(d3)
            # fix d2 on inner base
            rabd2[0][n1] = interp.interpolateLagrangeHermiteDerivative(rabx[0][n1], x, d2, 0.0)
        ravtx [0].append(agx [0])
        ravtd1[0].append(agd1[0])
        ravtd2[0].append(agd2[0])
        ravtd3[0].append(agd3[0])

        # get points on left atrium over appendage, from aorta to laa end on vestibule top
        endDerivative = [ -d for d in lavtd2[1][lan1Mid] ]
        # for not commonLeftRightPvOstium need a reliable location for laml - sample 3 across and use point 1
        laoax, laoad1, laoad2, laoad3, laoaProportions = laTrackSurface.createHermiteCurvePoints(
            lavtProportions[1][0], lavtProportions[1][1],
            lavtProportions[lan1Mid][0], lavtProportions[lan1Mid][1],
            elementsCount=(elementsCountAroundLeftAtriumRPV + elementsCountOverSideLeftAtriumLPV) if commonLeftRightPvOstium else 3,
            derivativeStart=lavtd2[1][1],
            derivativeEnd=endDerivative)
        if (not commonLeftRightPvOstium) and (elementsCountOverSideLeftAtriumLPV != 2):
            # resample from the laml line and merge
            _laoax, _laoad1, _laoad2, _laoad3, _laoaProportions = laTrackSurface.createHermiteCurvePoints(
                laoaProportions[1][0], laoaProportions[1][1],
                lavtProportions[lan1Mid][0], lavtProportions[lan1Mid][1],
                elementsCount=elementsCountOverSideLeftAtriumLPV,
                derivativeStart=laoad1[1],
                derivativeEnd=endDerivative)
            laoax  = laoax [:1] + _laoax
            laoad1 = laoad1[:1] + _laoad1
            laoad2 = laoad2[:1] + _laoad2
            laoad3 = laoad3[:1] + _laoad3
            laoaProportions = laoaProportions[:1] + _laoaProportions
        # d2 magnitudes are the same as d1 by default; want these to be even instead
        # transition between start and end side derivatives d1
        mag0 = vector.magnitude(lavtd1[1][1])
        magn = vector.magnitude(lavtd1[1][lan1Mid])
        # could blend accurately with arc lengths; try indexes:
        laoaCount = len(laoaProportions)
        for n1 in range(1, laoaCount - 1):
            xi = n1 / laoaCount
            laoad2[n1] = vector.setMagnitude(laoad2[n1], (1.0 - xi)*mag0 + xi*magn)
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
            laoax [n3][ 0] = lavtx [n3][1]
            laoad1[n3][ 0] = lavtd1[n3][1]
            laoad2[n3][ 0] = lavtd2[n3][1]
            laoad3[n3][ 0] = lavtd3[n3][1]
            laoax [n3][-1] = lavtx [n3][lan1Mid]
            laoad1[n3][-1] = lavtd1[n3][lan1Mid]
            laoad2[n3][-1] = lavtd2[n3][lan1Mid]
            laoad3[n3][-1] = lavtd3[n3][lan1Mid]

        if not commonLeftRightPvOstium:
            # get left atrium venous mid line from 3rd point on laoa to point on lavt at venous midpoint
            # find and pass through midpoint between left and right PVs
            n1lpv = 0
            n1rpv = elementsCountOverLeftAtriumVenous + elementsCountAroundLeftAtriumRPV
            n1End = elementsCountAroundLeftAtriumFreeWall - elementsCountAroundLeftAtriumRPV
            rpvoProportion1, rpvoProportion2 = laTrackSurface.getProportion(rpvoPositions[n1rpv])
            lpvoProportion1, lpvoProportion2 = laTrackSurface.getProportion(lpvoPositions[n1lpv])
            mpx, mpd1, mpd2, _, mpProportions = laTrackSurface.createHermiteCurvePoints(rpvoProportion1, rpvoProportion2, lpvoProportion1, lpvoProportion2,
                elementsCount = 2, derivativeStart = rpvod2[1][n1rpv], derivativeEnd = [ -d for d in lpvod2[1][n1lpv] ])
            # scale mid derivative 2 to be mean of d1 in LPV, RPV
            d2mag = 0.5*vector.magnitude(lpvod1[1][n1lpv]) + 0.5*vector.magnitude(rpvod1[1][n1rpv])
            mpd2[1] = vector.setMagnitude(mpd2[1], d2mag)
            lamlx, lamld2, lamld1, lamld3, lamlProportions = laTrackSurface.createHermiteCurvePoints(
                laoaProportions[1][0], laoaProportions[1][1],
                mpProportions[1][0], mpProportions[1][1],
                elementsCount = elementsCountOverLeftAtriumVenous//2 + 1,
                derivativeStart = laoad2[1][1],
                derivativeEnd = mpd2[1])
            _lamlx, _lamld2, _lamld1, _lamld3, _lamlProportions = laTrackSurface.createHermiteCurvePoints(
                mpProportions[1][0], mpProportions[1][1],
                lavtProportions[n1End][0], lavtProportions[n1End][1],
                elementsCount = elementsCountOverLeftAtriumVenous//2 + 1,
                derivativeStart = mpd2[1],
                derivativeEnd = [ -d for d in lavtd2[0][n1End]])
            lamlx  += _lamlx [1:]
            lamld1 += _lamld1[1:]
            lamld2 += _lamld2 [1:]
            lamld3 += _lamld3 [1:]
            lamlProportions += _lamlProportions[1:]
            # reverse d1
            lamld1 = [ [ -d for d in d1 ] for d1 in lamld1 ]
            # smooth d2
            lamld2 = interp.smoothCubicHermiteDerivativesLine(lamlx, lamld2, fixAllDirections=True, fixStartDerivative=True, fixEndDerivative=True)
            # give d1 the same magnitude as the smoothed d2
            for n in range(1, len(lamld1)):
                lamld1[n] = vector.setMagnitude(lamld1[n], vector.magnitude(lamld2[n]))
            # get inner points
            lamlx  = [ [ None ], lamlx  ]
            lamld1 = [ [ None ], lamld1 ]
            lamld2 = [ [ None ], lamld2 ]
            lamld3 = [ [ None ], lamld3 ]
            lamlProportions += _lamlProportions[1:]
            for n2 in range(1, elementsCountOverLeftAtriumVenous + 3):
                x, d2, d1, d3 = interp.projectHermiteCurvesThroughWall(lamlx[1], lamld2[1], [ [ -d for d in d1 ] for d1 in lamld1[1] ], n2, -laVenousFreeWallThickness)
                lamlx [0].append(x)
                lamld1[0].append([ -d for d in d1 ])
                lamld2[0].append(d2)
                lamld3[0].append(d3)
                lamld3[1][n2] = d3
            # substitute known start and end coordinates
            for n3 in range(2):
                lamlx [n3][ 0] = laoax [n3][1]
                lamld1[n3][ 0] = laoad1[n3][1]
                lamld2[n3][ 0] = laoad2[n3][1]
                lamld3[n3][ 0] = laoad3[n3][1]
                lamlx [n3][-1] = lavtx [n3][n1End]
                lamld1[n3][-1] = lavtd1[n3][n1End]
                lamld2[n3][-1] = lavtd2[n3][n1End]
                lamld3[n3][-1] = lavtd3[n3][n1End]

        if not commonLeftRightPvOstium:
            # get points on row above left atrium venous anterior, from interatrial septum to laml[1]
            # sample points up to venous midpoint, between RPV and laoa
            agn1va = elementsCountOverLeftAtriumNonVenousAnterior
            lavbx  = [ agx [agn1va] ]
            lavbd1 = [ agd1[agn1va] ]
            lavbd2 = [ agd2[agn1va] ]
            lavbd3 = [ agd3[agn1va] ]
            lavbd2inner = [ None ]
            lavbProportions = [ [ aVenousAnteriorOver, 0.0 ] ]
            for n1 in range(elementsCountAroundLeftAtriumRPV - 1):
                n1rpv = -elementsCountOverLeftAtriumVenous//2 - n1 - 1
                rpvoProportion1, rpvoProportion2 = laTrackSurface.getProportion(rpvoPositions[n1rpv])
                startDerivative = [ (1.5*laoad2[1][n1][c] - 0.5*laoad1[1][n1][c]) for c in range(3) ] if (n1 == 0) else laoad2[1][n1]  # GRC fudge factor
                _x, _d2, _d1, _d3, _proportions = laTrackSurface.createHermiteCurvePoints(
                    laoaProportions[n1][0], laoaProportions[n1][1],
                    rpvoProportion1, rpvoProportion2,
                    elementsCount=1 + elementsCountRadialPVAnnuli,
                    derivativeStart=startDerivative,
                    derivativeEnd=[ -d for d in rpvod2[1][n1rpv] ])
                lavbx .append(_x [1])
                lavbd1.append([ -d for d in _d1[1] ])
                lavbd2.append(_d2[1])
                lavbd3.append(_d3[1])
                lavbProportions.append(_proportions[1])
                # get precise inner d2
                _d2inner = interp.projectHermiteCurvesThroughWall(_x, _d2, _d1, 1, -laVenousFreeWallThickness)[1]
                lavbd2inner.append(_d2inner)
            # use end values from laml[1][1], with d1 at 30 degree angle
            cos30 = math.cos(math.pi/6.0)
            sin30 = math.sin(math.pi/6.0)
            lavbx .append(lamlx [1][1])
            lavbd1.append([ (cos30*lamld1[1][1][c] + sin30*lamld2[1][1][c]) for c in range(3) ])
            lavbd2.append(lamld2[1][1])
            lavbd3.append(lamld3[1][1])
            lavbd2inner.append(lamld2[0][1])
            lavbProportions.append(lamlProportions[1])
            # smooth d1:
            lavbd1 = interp.smoothCubicHermiteDerivativesLine(lavbx, lavbd1, fixAllDirections = True,
                fixStartDerivative=True, fixEndDerivative=True, magnitudeScalingMode=interp.DerivativeScalingMode.ARITHMETIC_MEAN)
            # get inner points
            asn1va = elementsCountAroundAtrialSeptum - 1 + elementsCountOverLeftAtriumNonVenousAnterior
            lavbx  = [ [ asx [0][asn1va] ], lavbx  ]
            lavbd1 = [ [ asd1[0][asn1va] ], lavbd1 ]
            lavbd2 = [ [ asd2[0][asn1va] ], lavbd2 ]
            lavbd3 = [ [ asd3[0][asn1va] ], lavbd3 ]
            for n1 in range(1, len(lavbx[1])):
                x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(lavbx[1], lavbd1[1], lavbd2[1], n1, -laVenousFreeWallThickness)
                lavbx [0].append(x)
                lavbd1[0].append(d1)
                lavbd2[0].append(lavbd2inner[n1])
                lavbd3[0].append(d3)
                lavbd3[1][n1] = d3
            # copy end d1 values back to lamld1
            lamld1[0][1] = lavbd1[0][-1]
            lamld1[1][1] = lavbd1[1][-1]

        agn1vp = elementsCountOverLeftAtriumNonVenousAnterior + elementsCountOverLeftAtriumVenous
        if not commonLeftRightPvOstium:
            # get points on row above left atrium venous posterior, from interatrial septum to laml[-2]
            # sample points up to venous midpoint, between RPV and vestibule top
            lavqx  = [ agx [agn1vp] ]
            lavqd1 = [ agd1[agn1vp] ]
            lavqd2 = [ agd2[agn1vp] ]
            lavqd3 = [ agd3[agn1vp] ]
            lavqd2inner = [ None ]
            lavqProportions = [ [ laVenousLimitPosterior, 0.0 ] ]
            for n1 in range(1, elementsCountAroundLeftAtriumRPV):
                n1rpv = elementsCountOverLeftAtriumVenous//2 + n1
                n1cs = elementsCountAroundLeftAtriumFreeWall - n1
                rpvoProportion1, rpvoProportion2 = laTrackSurface.getProportion(rpvoPositions[n1rpv])
                _x, _d2, _d1, _d3, _proportions = laTrackSurface.createHermiteCurvePoints(
                    rpvoProportion1, rpvoProportion2,
                    lavtProportions[n1cs][0], lavtProportions[n1cs][1],
                    elementsCount = 1 + elementsCountRadialPVAnnuli,
                    derivativeStart = rpvod2[1][n1rpv], derivativeEnd = [ -d for d in lavtd2[1][n1cs] ])
                lavqx .append(_x [-2])
                lavqd1.append([ -d for d in _d1[-2] ])
                lavqd2.append(_d2[-2])
                lavqd3.append(_d3[-2])
                lavqProportions.append(_proportions[-2])
                # get precise inner d2
                _d2inner  = interp.projectHermiteCurvesThroughWall(_x, _d2, _d1, len(_x) - 2, -laVenousFreeWallThickness)[1]
                lavqd2inner.append(_d2inner)
            # use end values from laml[1][-2], with d1 at 30 degree angle
            cos30 = math.cos(math.pi/6.0)
            sin30 = math.sin(math.pi/6.0)
            lavqx .append(lamlx [1][-2])
            lavqd1.append([ (cos30*lamld1[1][-2][c] - sin30*lamld2[1][-2][c]) for c in range(3) ])
            lavqd2.append(lamld2[1][-2])
            lavqd3.append(lamld3[1][-2])
            lavqd2inner.append(lamld2[0][-1])
            lavqProportions.append(lamlProportions[-2])
            # smooth d1:
            lavqd1 = interp.smoothCubicHermiteDerivativesLine(lavqx, lavqd1, fixAllDirections = True,
                fixStartDerivative=True, fixEndDerivative=True, magnitudeScalingMode=interp.DerivativeScalingMode.ARITHMETIC_MEAN)
            # get inner points
            asn1vp = -1
            lavqx  = [ [ asx [0][asn1vp] ], lavqx  ]
            lavqd1 = [ [ asd1[0][asn1vp] ], lavqd1 ]
            lavqd2 = [ [ asd2[0][asn1vp] ], lavqd2 ]
            lavqd3 = [ [ asd3[0][asn1vp] ], lavqd3 ]
            for n1 in range(1, len(lavqx[1])):
                x, d1, _, d3 = interp.projectHermiteCurvesThroughWall(lavqx[1], lavqd1[1], lavqd2[1], n1, -laVenousFreeWallThickness)
                lavqx [0].append(x)
                lavqd1[0].append(d1)
                lavqd2[0].append(lavqd2inner[n1])
                lavqd3[0].append(d3)
                lavqd3[1][n1] = d3
            # copy end d1 values back to lamld1
            lamld1[0][-2] = lavqd1[0][-1]
            lamld1[1][-2] = lavqd1[1][-1]

        # get points on right atrium along crista terminalis from aorta to posterior venous limit
        xi = (1.0 - aVenousMidpointOver - ravtProportions[ran1Aorta][0])/(ravtProportions[ran1Ctp][0] - ravtProportions[ran1Aorta][0])
        rctmpProportion1 = 1.0 - aVenousMidpointOver
        rctmpProportion2 = raVenousRight  # GRC was (failing) (1.0 - xi)*ravtProportions[ran1Aorta][1] + xi*ravtProportions[ran1Ctp][1]
        elementsCountOverCristaTerminalisAnterior = elementsCountOverRightAtriumVenous//2 + 1
        elementsCountOverCristaTerminalisPosterior = elementsCountOverRightAtriumVenous//2
        # trick to get lower derivative at midpoint: sample one element higher
        _, _d2 = raTrackSurface.createHermiteCurvePoints(
            rctmpProportion1, rctmpProportion2, ravtProportions[ran1Ctp][0], ravtProportions[ran1Ctp][1],
            elementsCount = elementsCountOverCristaTerminalisPosterior + 1,
            derivativeStart = None,
            derivativeEnd = [ -d for d in ravtd2[1][ran1Ctp] ])[0:2]
        useDerivativeStart = _d2[0]
        _ractx, _ractd2, _ractd1, _ractd3, _ractProportions = raTrackSurface.createHermiteCurvePoints(
            rctmpProportion1, rctmpProportion2, ravtProportions[ran1Ctp][0], ravtProportions[ran1Ctp][1],
            elementsCount = elementsCountOverCristaTerminalisPosterior,
            derivativeStart = useDerivativeStart,
            derivativeEnd = [ -d for d in ravtd2[1][ran1Ctp] ])
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
        # get right atrium posterior venous width around vestibule top to crista terminalis, to scale set ct derivative 1 below
        raVenousWidth = sum(interp.getCubicHermiteArcLength(ravtx[1][e], ravtd1[1][e], ravtx[1][e + 1], ravtd1[1][e + 1]) for e in range(elementsCountAroundRightAtriumPosteriorVenous))
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
        for n3 in range(2):
            # overwrite venous right x, d1 on vestibule top
            ravtx [n3][ran1Ctp] = ractx [n3][-1]
            d1mag = min(vector.magnitude(ravtd1[n3][ran1Ctp]), 1.0*vector.magnitude(ractd1[n3][-1]))  # GRC fudge factor
            ravtd1[n3][ran1Ctp] = vector.setMagnitude(ravtd1[n3][ran1Ctp], d1mag)
        # substitute known start and end coordinates
        ractx [0][ 0] = asx [1][elementsCountAroundAtrialSeptum]
        ractd1[0][ 0] = asd1[1][elementsCountAroundAtrialSeptum]
        ractd2[0][ 0] = asd2[1][elementsCountAroundAtrialSeptum]
        ractd3[0][ 0] = asd3[1][elementsCountAroundAtrialSeptum]
        ractx [1][ 0] = agx [1]
        ractd1[1][ 0] = agd1[1]
        ractd2[1][ 0] = agd2[1]
        ractd3[1][ 0] = agd3[1]
        for n3 in range(2):
            ractx [n3][-1] = ravtx [n3][ran1Ctp]
            ractd1[n3][-1] = ravtd1[n3][ran1Ctp]
            ractd2[n3][-1] = ravtd2[n3][ran1Ctp]
            ractd3[n3][-1] = ravtd3[n3][ran1Ctp]

        # get points on right atrium ridge midway between inferior and superior vena cavae from crista terminalis to interatrial groove
        # minimum of 3 points over top of venous component
        elementsCountOverSideRightAtriumVC = max(elementsCountAroundRightAtriumPosteriorVenous, 3)
        ravmx, ravmd1, ravmd2, ravmd3 = raTrackSurface.createHermiteCurvePoints(
            ractProportions[elementsCountOverCristaTerminalisAnterior][0], ractProportions[elementsCountOverCristaTerminalisAnterior][1],
            1.0 - aVenousMidpointOver, 0.0,
            elementsCount = elementsCountOverSideRightAtriumVC,
            derivativeStart = [d for d in ractd1[1][elementsCountOverCristaTerminalisAnterior]],
            derivativeEnd = sub(agd1[agn1Mid], agd3[agn1Mid]))[0:4]
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
            # copy d1 back to crista terminalis
            ractd1[n3][elementsCountOverCristaTerminalisAnterior] = ravmd1[n3][0]
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

        # get points over right atrial appendage from anterior end of crista terminalis to end of pouch on vestibule top
        assert elementsCountAroundRightAtrialAppendagePlainBase == 2
        elementsCountOverSideRightAtriumPouch = elementsCountOverAtria//2 + elementsCountOverCristaTerminalisAnterior - 3
        ran1raap = ran1Ctp + elementsCountAroundRightAtrialAppendagePlainBase
        nc = 2
        raapx, raapd1, raapd2, raapd3, raapProportions = raTrackSurface.createHermiteCurvePoints(
            ractProportions[nc][0], ractProportions[nc][1],
            ravtProportions[ran1raap][0], ravtProportions[ran1raap][1],
            elementsCount = elementsCountOverSideRightAtriumPouch,
            derivativeStart = [ (0.25*ractd2[1][nc][c] - 1.5*ractd1[1][nc][c]) for c in range(3) ],
            #derivativeStart = [ -d for d in ractd1[1][nc] ],
            derivativeEnd = [ -d for d in ravtd2[1][ran1raap] ])
        # get inner points
        raapx  = [ [ None ], raapx  ]
        raapd1 = [ [ None ], raapd1 ]
        raapd2 = [ [ None ], raapd2 ]
        raapd3 = [ [ None ], raapd3 ]
        for n in range(1, len(raapx[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(raapx[1], raapd1[1], raapd2[1], n, -raaWallThickness)
            raapx [0].append(x)
            raapd1[0].append(d1)
            raapd2[0].append(d2)
            raapd3[0].append(d3)
            raapd3[1][n] = d3
        # substitute known end coordinates
        for n3 in range(2):
            raapx [n3][ 0]  = ractx [n3][nc]
            raapd1[n3][ 0]  = ractd1[n3][nc]
            raapd2[n3][ 0]  = ractd2[n3][nc]
            raapd3[n3][ 0]  = ractd3[n3][nc]
            raapx [n3][-1]  = ravtx [n3][ran1raap]
            raapd1[n3][-1]  = ravtd1[n3][ran1raap]
            raapd2[n3][-1]  = ravtd2[n3][ran1raap]
            raapd3[n3][-1]  = ravtd3[n3][ran1raap]
        # get second row between raap and crista terminalis: raaq
        ran1raaq = ran1Ctp + elementsCountAroundRightAtrialAppendagePlainBase - 1
        nc = 2
        raaqx, raaqd1, raaqd2, raaqd3, raaqProportions = raTrackSurface.createHermiteCurvePoints(
            ractProportions[nc][0], ractProportions[nc][1],
            ravtProportions[ran1raaq][0], ravtProportions[ran1raaq][1],
            elementsCount = elementsCountOverAtria//2 + elementsCountOverCristaTerminalisAnterior - 3,
            #derivativeStart = [ (ractd2[1][nc][c] - 0.5*ractd1[1][nc][c]) for c in range(3) ],
            derivativeStart = [ (0.5*ractd2[1][nc][c] - ractd1[1][nc][c]) for c in range(3) ],
            derivativeEnd = [ -d for d in ravtd2[1][ran1raaq] ])
        # get inner points
        raaqx  = [ [ None ], raaqx  ]
        raaqd1 = [ [ None ], raaqd1 ]
        raaqd2 = [ [ None ], raaqd2 ]
        raaqd3 = [ [ None ], raaqd3 ]
        for n in range(1, len(raaqx[1])):
            x, d1, d2, d3 = interp.projectHermiteCurvesThroughWall(raaqx[1], raaqd1[1], raaqd2[1], n, -raaWallThickness)
            raaqx [0].append(x)
            raaqd1[0].append(d1)
            raaqd2[0].append(d2)
            raaqd3[0].append(d3)
            raaqd3[1][n] = d3
        # substitute known end coordinates
        for n3 in range(2):
            raaqx [n3][ 0]  = ractx [n3][nc]
            raaqd1[n3][ 0]  = ractd1[n3][nc]
            raaqd2[n3][ 0]  = ractd2[n3][nc]
            raaqd3[n3][ 0]  = ractd3[n3][nc]
            raaqx [n3][-1]  = ravtx [n3][ran1raaq]
            raaqd1[n3][-1]  = ravtd1[n3][ran1raaq]
            raaqd2[n3][-1]  = ravtd2[n3][ran1raaq]
            raaqd3[n3][-1]  = ravtd3[n3][ran1raaq]
        # smooth d2 between raap and raaq
        for n1 in range(2, len(raapx[1]) - 1):
            nx  = [ raaqx [1][n1], raapx [1][n1] ]
            nd1 = [ raaqd2[1][n1], raapd2[1][n1] ]
            nd2 = [ [ -d for d in d1 ] for d1 in [ raaqd1[1][n1], raapd1[1][n1] ] ]
            raaqd2[1][n1], raapd2[1][n1] = nd1 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixAllDirections = True)
            # get inner derivatives
            _, raaqd2[0][n1], _, _ = interp.projectHermiteCurvesThroughWall(nx, nd1, nd2, 0, -raaWallThickness)
            _, raapd2[0][n1], _, _ = interp.projectHermiteCurvesThroughWall(nx, nd1, nd2, 1, -raaWallThickness)
        # make point in centre of triangle at top of raaq, for 3 quad element junction
        raaqProportions[1] = [ (ractProportions[nc + 1][i] + raapProportions[1][i] + raaqProportions[2][i])/3.0 for i in range(2) ]
        position = raTrackSurface.createPositionProportion(raaqProportions[1][0], raaqProportions[1][1])
        raaqx[1][1], d1, d2 = raTrackSurface.evaluateCoordinates(position, derivatives = True)
        # calculate derivative 1 there to fit nearest point on raaq
        raaqd1[1][1] = interp.interpolateLagrangeHermiteDerivative(raaqx[1][1], raaqx[1][2], raaqd1[1][2], 0.0)
        raaqd1[1][1] = vector.setMagnitude(calculate_surface_axes(d1, d2, raaqd1[1][1])[0], vector.magnitude(raaqd1[1][1]))
        raaqd1[1][1] = interp.smoothCubicHermiteDerivativesLine(raaqx[1][1:3], raaqd1[1][1:3], fixAllDirections = True, fixEndDerivative = True)[0]
        # calculate derivative 2 there and on nearest point on raap
        # raapd2[1][1] magnitude needs to be set to fit distance between nodes:
        d2mag = math.sqrt(sum((raapx[1][1][c] - raaqx[1][1][c])*(raapx[1][1][c] - raaqx[1][1][c]) for c in range(3)))
        raapd2[1][1] = vector.setMagnitude(raapd2[1][1], d2mag)  # must reduce this otherwise smooth will converge wrongly
        raaqd2[1][1] = interp.interpolateLagrangeHermiteDerivative(raaqx[1][1], raapx[1][1], raapd2[1][1], 0.0)
        raaqd2[1][1] = vector.setMagnitude(calculate_surface_axes(d1, d2, raaqd2[1][1])[0], vector.magnitude(raaqd2[1][1]))
        raaqd2[1][1], raapd2[1][1] = interp.smoothCubicHermiteDerivativesLine([ raaqx[1][1], raapx[1][1] ], [ raaqd2[1][1], raapd2[1][1] ], fixAllDirections = True)
        # get inner coordinates and derivatives
        raaqx[0][1], raaqd1[0][1], _, raaqd3[0][1] = interp.projectHermiteCurvesThroughWall(raaqx[1][1:3], raaqd1[1][1:3], raaqd2[1][1:3], 0, -raaWallThickness)
        raaqd3[1][1] = raaqd3[0][1]
        nx  = [ raaqx [1][1], raapx [1][1] ]
        nd1 = [ raaqd2[1][1], raapd2[1][1] ]
        nd2 = [ [ -d for d in d1 ] for d1 in [ raaqd1[1][1], raapd1[1][1] ] ]
        _, raaqd2[0][1], _, _ = interp.projectHermiteCurvesThroughWall(nx, nd1, nd2, 0, -raaWallThickness)
        _, raapd2[0][1], _, _ = interp.projectHermiteCurvesThroughWall(nx, nd1, nd2, 1, -raaWallThickness)

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

        # create septum nodes, along vestibule top and over arch
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

        # create left atrium vestibule top nodes
        # start and end with common nodes on interatrial groove or septum arch
        lavtNodeId = [ [ asNodeId[0][elementsCountAroundAtrialSeptum] ], [ agNodeId[1] ] ]
        for n3 in range(2):
            for n1 in range(1, len(lavtx[n3]) - 1):
                if elementsCountAroundLeftAtriumAorta < n1 < lan1Mid:
                    # left atrial appendage
                    lavtNodeId[n3].append(None)
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                lavtNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lavtx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lavtd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lavtd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lavtd3[n3][n1])
                nodeIdentifier += 1
        lavtNodeId[0].append(asNodeId[0][0])
        lavtNodeId[1].append(agNodeId[-2])

        # create right atrium vestibule top nodes
        ravtNodeId = [ [ asNodeId[1][0] ], [ agNodeId[-2] ] ]
        for n3 in range(2):
            for n1 in range(1, len(ravtx[n3]) - 1):
                if ran1raap < n1 <= ran1Aorta:
                    # right atrial appendage; aorta substituted from ract later
                    ravtNodeId[n3].append(None)
                    continue
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                ravtNodeId[n3].append(nodeIdentifier)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ravtx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ravtd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ravtd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, ravtd3[n3][n1])
                nodeIdentifier += 1
        ravtNodeId[0].append(asNodeId[1][elementsCountAroundAtrialSeptum])
        ravtNodeId[1].append(agNodeId[1])

        # create nodes on left atrium over appendage
        laoaNodeId = [ [ lavtNodeId[0][1] ], [ lavtNodeId[1][1] ] ]
        for n3 in range(2):
            for n1 in range(1, len(laoax[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, laoax [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, laoad1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, laoad2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, laoad3[n3][n1])
                laoaNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            laoaNodeId[n3].append(lavtNodeId[n3][lan1Mid])

        if not commonLeftRightPvOstium:
            # create left atrium venous midline nodes
            lamlNodeId = [ [], [] ]
            for n3 in range(2):
                lamlNodeId[n3].append(laoaNodeId[n3][2])
                for n2 in range(1, len(lamlx[n3]) - 1):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, lamlx [n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, lamld1[n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, lamld2[n3][n2])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, lamld3[n3][n2])
                    lamlNodeId[n3].append(nodeIdentifier)
                    nodeIdentifier += 1
                lamlNodeId[n3].append(lavtNodeId[n3][-elementsCountAroundLeftAtriumRPV - 1])

        if not commonLeftRightPvOstium:
            # create nodes on row above left atrium venous anterior to laml[1]
            lavbNodeId = [ [ asNodeId[0][asn1va] ], [ agNodeId[elementsCountOverLeftAtriumNonVenousAnterior] ] ]
            n1lpv = -elementsCountOverLeftAtriumVenous//2
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
                lavbNodeId[n3].append(lamlNodeId[n3][1])

        # create nodes on row above left atrium venous posterior to laml[-2]
        if not commonLeftRightPvOstium:
            lavqNodeId = [ [ asNodeId[0][-1] ], [ agNodeId[elementsCountOverLeftAtriumNonVenousAnterior + elementsCountOverLeftAtriumVenous ] ] ]
            n1lpv = elementsCountOverLeftAtriumVenous//2
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
                lavqNodeId[n3].append(lamlNodeId[n3][-2])

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
            ractNodeId[n3].append(ravtNodeId[n3][ran1Ctp])
            # use second ract node as aorta vestibule top node
            ravtNodeId[n3][ran1Aorta] = ractNodeId[n3][1]

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

        # create right atrial appendage plain boundary nodes
        raapNodeId = [ [], [] ]
        for n3 in range(2):
            raapNodeId[n3].append(ractNodeId[n3][2])
            for n1 in range(1, len(raapx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, raapx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, raapd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, raapd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, raapd3[n3][n1])
                raapNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            raapNodeId[n3].append(ravtNodeId[n3][ran1raap])
        # and middle row on plain appendage
        raaqNodeId = [ [], [] ]
        for n3 in range(2):
            raaqNodeId[n3].append(ractNodeId[n3][2])
            for n1 in range(1, len(raaqx[n3]) - 1):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, raaqx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, raaqd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, raaqd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, raaqd3[n3][n1])
                raaqNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
            raaqNodeId[n3].append(ravtNodeId[n3][ran1raaq])

        if False:
            # create lt nodes:
            for n1 in range(len(ltBaseOuterx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ltBaseOuterx [n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ltBaseOuterd1[n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ltBaseOuterd2[n1])
                nodeIdentifier += 1

        if defineEpicardiumLayer:
            # epicardial fat pad in RAGP region -- track surface bridging interatrial groove
            fpx = []
            fpd1 = []
            fpd2 = []
            proportionAcross = 0.2
            for i in range(2):
                if i == 0:
                    nx, nd2, nd1, nd3, nProportions = laTrackSurface.createHermiteCurvePoints(
                        0.0, proportionAcross, 1.0, proportionAcross, elementsCountAcrossTrackSurface)
                else:
                    nx, nd2, nd1, nd3 = raTrackSurface.createHermiteCurvePoints(
                        1.0, proportionAcross, 0.0, proportionAcross, elementsCountAcrossTrackSurface)[0:4]
                if nx:
                    nx = [add(nx[i], mult(nd3[i], epicardiumLayerMinimumThickness)) for i in range(len(nx))]
                    nd1 = [[-c for c in d] for d in nd1]
                fpx.append(nx)
                fpd1.append(nd1)
                fpd2.append(nd2)

            # put into single arrays cycling left to right fastest, smoothing each d1 row
            nx = []
            nd1 = []
            nd2 = []
            for n in range(elementsCountAcrossTrackSurface + 1):
                scale = interp.computeCubicHermiteDerivativeScaling(fpx[0][n], fpd1[0][n], fpx[1][n], fpd1[1][n])
                for i in range(2):
                    nx.append(fpx[i][n])
                    fpd1[i][n] = mult(fpd1[i][n], scale)
                    nd1.append(fpd1[i][n])
                    nd2.append(fpd2[i][n])

            fpTrackSurface = TrackSurface(1, elementsCountAcrossTrackSurface, nx, nd1, nd2)

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
        drawFpTrackSurface = False
        if drawFpTrackSurface:
            # create track surface nodes:
            fpTrackSurfaceFirstNodeIdentifier = nodeIdentifier
            for n in range((fpTrackSurface.elementsCount2 + 1)*(fpTrackSurface.elementsCount1 + 1)):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, fpTrackSurface.nx[n] )
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, fpTrackSurface.nd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, fpTrackSurface.nd2[n])
                nodeIdentifier += 1

        #################
        # Create elements
        #################

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # left atrium free wall elements to vestibule top, starting at cfb / anterior interatrial sulcus
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
                labNodeId[0][e1], labNodeId[0][e1 + 1], lavtNodeId[0][e1], lavtNodeId[0][e1 + 1],
                labNodeId[1][e1], labNodeId[1][e1 + 1], lavtNodeId[1][e1], lavtNodeId[1][e1 + 1]]
            if None in nids:
                continue  # left atrial appendage
            scalefactors = None
            meshGroups = [ heartMeshGroup, lamMeshGroup ]
            if e1 == -1:
                # cfb/anterior interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = rabNodeId[0][-elementsCountAroundAtrialSeptum]
                nids[2] = ravtNodeId[0][-1]
                nids.pop(6)
                nids.pop(4)
                meshGroups += [ ramMeshGroup ]
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
            if len(nids) == 6:
                aSeptumMeshGroup.addElement(element)

        # right atrium free wall elements to vestibule top, starting at crux / posterior interatrial sulcus
        for e1 in range(-1, elementsCountAroundRightAtriumFreeWall):
            eft1 = eft
            elementtemplate1 = elementtemplate
            nids = [
                rabNodeId[0][e1], rabNodeId[0][e1 + 1], ravtNodeId[0][e1], ravtNodeId[0][e1 + 1],
                rabNodeId[1][e1], rabNodeId[1][e1 + 1], ravtNodeId[1][e1], ravtNodeId[1][e1 + 1]]
            if None in nids:
                continue  # right atrial appendage
            scalefactors = None
            meshGroups = [ heartMeshGroup, ramMeshGroup ]
            # Anderson definition of right atrial appendage starts at crista terminalis:
            #if (e1 >= elementsCountAroundRightAtriumPosteriorVenous) and (e1 < elementsCountAroundRightAtriumFreeWall - elementsCountAroundRightAtriumAorta - 1):
            #    meshGroups += [ raaMeshGroup ]
            if e1 == -1:
                # crux/posterior interatrial groove straddles left and right atria, collapsed to 6 node wedge
                nids[0] = labNodeId[0][elementsCountAroundLeftAtriumFreeWall]
                nids[2] = lavtNodeId[0][elementsCountAroundLeftAtriumFreeWall]
                nids.pop(6)
                nids.pop(4)
                meshGroups += [ lamMeshGroup ]
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
            if len(nids) == 6:
                aSeptumMeshGroup.addElement(element)

        if commonLeftRightPvOstium:
            # left atrium row above vestibule beside aorta 
            meshGroups = [ heartMeshGroup, lamMeshGroup ]
            for e1 in range(elementsCountAroundLeftAtriumAorta):
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [ lavtNodeId[0][e1], lavtNodeId[0][e1 + 1], laoaNodeId[0][e1], laoaNodeId[0][e1 + 1],
                         lavtNodeId[1][e1], lavtNodeId[1][e1 + 1], laoaNodeId[1][e1], laoaNodeId[1][e1 + 1] ]
                scalefactors = None
                if e1 == 0:
                    nids[2] = asNodeId[0][elementsCountAroundAtrialSeptum + 1]
                    nids[6] = agNodeId[2]
                    # general linear map d3 adjacent to interatrial groove
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 4, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                    scalefactors = [ -1.0 ]
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                #print('create element laao', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
        else:
            # left atrium extra row between appendage and RPV, anterior
            meshGroups = [ heartMeshGroup, lamMeshGroup ]
            for e1 in range(elementsCountAroundLeftAtriumRPV):
                eft1 = eft
                elementtemplate1 = elementtemplate
                nids = [ lavtNodeId[0][0] if (e1 == 0) else laoaNodeId[0][e1 - 1], laoaNodeId[0][e1], lavbNodeId[0][e1], lavbNodeId[0][e1 + 1],
                         lavtNodeId[1][0] if (e1 == 0) else laoaNodeId[1][e1 - 1], laoaNodeId[1][e1], lavbNodeId[1][e1], lavbNodeId[1][e1 + 1] ]
                scalefactors = None
                if e1 == 0:
                    # general linear map d3 adjacent to interatrial groove
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 5, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    scalefactors = [ -1.0 ]
                elif e1 == 1:
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, []) ])
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, []) ])
                    scalefactors = [ -1.0 ]
                if eft != eft1:
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                #print('create element lavb', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        if not commonLeftRightPvOstium:
            # left atrium first row of elements above vestibule, posterior
            meshGroups = [ heartMeshGroup, lamMeshGroup ]
            scalefactors = [ -1.0 ]
            elementsCount = len(lavqx[1]) - 1 if commonLeftRightPvOstium else elementsCountAroundLeftAtriumRPV
            for e1 in range(elementsCount):
                nc = elementsCountAroundLeftAtriumFreeWall - e1
                nids = [ lavqNodeId[0][e1], lavqNodeId[0][e1 + 1], lavtNodeId[0][nc], lavtNodeId[0][nc - 1],
                            lavqNodeId[1][e1], lavqNodeId[1][e1 + 1], lavtNodeId[1][nc], lavtNodeId[1][nc - 1] ]
                eft1 = tricubichermite.createEftNoCrossDerivatives()
                setEftScaleFactorIds(eft1, [1], [])
                scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS1 ], [ 1 ])
                scaleEftNodeValueLabels(eft1, [ 3, 4, 7, 8 ], [ Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                if e1 == 0:
                    # general linear map d3 adjacent to collapsed inter-atrial groove
                    remapEftNodeValueLabel(eft1, [ 1, 7 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                    remapEftNodeValueLabel(eft1, [ 3, 5 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, []) ])
                elif commonLeftRightPvOstium and (e1 == (elementsCount - 1)):
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1]) ])
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, []) ])
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                result3 = element.setScaleFactors(eft1, scalefactors)
                #print('create element lavq', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        # create atrial septum base row of elements
        meshGroups = [ heartMeshGroup, lamMeshGroup, ramMeshGroup, aSeptumMeshGroup ]
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
        meshGroups = [ heartMeshGroup, lamMeshGroup, ramMeshGroup, aSeptumMeshGroup, fossaMeshGroup ]
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
        meshGroups = [ heartMeshGroup, lamMeshGroup, ramMeshGroup, aSeptumMeshGroup ]
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
        meshGroups = [ heartMeshGroup, lamMeshGroup, ramMeshGroup, aSeptumMeshGroup ]
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
        if commonLeftRightPvOstium:
            # set points clockwise from interatrial groove at anterior venous limit
            # insert at indexes such that 0 is the midpoint on interatrial groove
            ix = interp.getNearestPointIndex(lpvox[1], agx[agn1Mid]) - elementsCountOverLeftAtriumVenous//2
            if ix > 0:
                ix -= elementsCountAroundLpvOstium
            # down interatrial groove from anterior venous limit, including both corners
            for n1 in range(elementsCountOverLeftAtriumVenous + 2):
                ns = n1 - elementsCountOverLeftAtriumVenous - 1
                ng = elementsCountOverLeftAtriumNonVenousAnterior + n1
                lpvax [0][ix] = asx [0][ns]
                lpvad1[0][ix] = asd1[0][ns]
                lpvad2[0][ix] = asd2[0][ns]
                lpvad3[0][ix] = asd3[0][ns]
                lpvaNodeId[0][ix] = asNodeId[0][ns]
                lpvax [1][ix] = agx [ng]
                lpvad1[1][ix] = agd1[ng]
                lpvad2[1][ix] = agd2[ng]
                lpvad3[1][ix] = agd3[ng]
                lpvaNodeId[1][ix] = agNodeId[ng]
                if n1 == 0:
                    lpvaDerivativesMap[0][ix] = ( (-1, 0, 0), (-1, -1, 0), (1, 0, 1), (0, 1, 0 ) )
                    lpvaDerivativesMap[1][ix] = ( (-1, 0, 0), (-1, -1, 0), (-1, 0, 1), (0, 1, 0 ) )
                elif n1 == elementsCountOverLeftAtriumVenous:
                    lpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 1, 0), (1, 0, 1) )
                    lpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 1, -1), (-1, 0, 1) )
                elif n1 == (elementsCountOverLeftAtriumVenous + 1):
                    lpvaDerivativesMap[0][ix] = ( (0, -1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 0 ) )
                    lpvaDerivativesMap[1][ix] = ( (0, -1, 0), (1, -1, -1), (-1, 0, 1), (-1, 0, 0 ) )
                else:
                    lpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 0, 0), (1, 0, 1) )
                    lpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 0, -1), (-1, 0, 1) )
                ix += 1
            # left around posterior vestibule top lavt
            for n1 in range(1, elementsCountAroundLeftAtriumRPV + elementsCountAroundLeftAtriumLPV):
                nc = elementsCountAroundLeftAtriumFreeWall - n1
                for n3 in range(2):
                    lpvax [n3][ix] = lavtx [n3][nc]
                    lpvad1[n3][ix] = lavtd1[n3][nc]
                    lpvad2[n3][ix] = lavtd2[n3][nc]
                    lpvad3[n3][ix] = lavtd3[n3][nc]
                    lpvaNodeId[n3][ix] = lavtNodeId[n3][nc]
                    lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                ix += 1
            # up left atrium over appendage laoa to interatrial groove
            for n1 in range(laoaCount - 1):
                no = -1 - n1
                for n3 in range(2):
                    lpvax [n3][ix] = laoax [n3][no]
                    lpvad1[n3][ix] = laoad1[n3][no]
                    lpvad2[n3][ix] = laoad2[n3][no]
                    lpvad3[n3][ix] = laoad3[n3][no]
                    lpvaNodeId[n3][ix] = laoaNodeId[n3][no]
                    if n1 == 0:
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0 ) )
                    elif n1 == (laoaCount - 2):
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0 ) )
                    else:
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                ix += 1
        else:
            # set points clockwise from venous midpoint at anterior venous limit
            # insert at indexes such that 0 is one past the midpoint on venous midline
            ix = -elementsCountOverLeftAtriumVenous//2
            # down left atrium venous midpoint line
            for n1 in range(1, elementsCountOverLeftAtriumVenous + 2):
                for n3 in range(2):
                    lpvax [n3][ix] = lamlx [n3][n1]
                    lpvad1[n3][ix] = lamld1[n3][n1]
                    lpvad2[n3][ix] = lamld2[n3][n1]
                    lpvad3[n3][ix] = lamld3[n3][n1]
                    lpvaNodeId[n3][ix] = lamlNodeId[n3][n1]
                    lpvaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
                ix += 1
            # left around left vestibule top, including 2 corners
            for n1 in range(elementsCountAroundLeftAtriumLPV + 1):
                nc = lan1MidVenous - n1
                for n3 in range(2):
                    lpvax [n3][ix] = lavtx [n3][nc]
                    lpvad1[n3][ix] = lavtd1[n3][nc]
                    lpvad2[n3][ix] = lavtd2[n3][nc]
                    lpvad3[n3][ix] = lavtd3[n3][nc]
                    lpvaNodeId[n3][ix] = lavtNodeId[n3][nc]
                    if n1 == 0:
                        lpvaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, -1, 0), None, (-1, 0, 0 ) )
                    elif n1 == elementsCountAroundLeftAtriumLPV:
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0 ) )
                    else:
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                ix += 1
            # up left atrium venous side line
            for n1 in range(elementsCountOverSideLeftAtriumLPV):
                no = -2 - n1
                for n3 in range(2):
                    lpvax [n3][ix] = laoax [n3][no]
                    lpvad1[n3][ix] = laoad1[n3][no]
                    lpvad2[n3][ix] = laoad2[n3][no]
                    lpvad3[n3][ix] = laoad3[n3][no]
                    lpvaNodeId[n3][ix] = laoaNodeId[n3][no]
                    if n1 == (elementsCountOverSideLeftAtriumLPV - 1):
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (-1, -1, 0), None, (0, 1, 0) )
                    else:
                        lpvaDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
                ix += 1
        #print('lpvaNodeId[0]',lpvaNodeId[0])
        #print('lpvaNodeId[1]',lpvaNodeId[1])
        #print('lpvaDerivativesMap[0]',lpvaDerivativesMap[0])
        #print('lpvaDerivativesMap[1]',lpvaDerivativesMap[1])
            #if commonLeftRightPvOstium:
        lpvoProportions = [ list(laTrackSurface.getProportion(lpvoPositions[i])) for i in range(elementsCountAroundLpvOstium) ]
        lpvaProportions = [ list(laTrackSurface.getProportion(laTrackSurface.findNearestPosition(lpvax[1][i], startPosition=lpvoPositions[i])))
                            for i in range(elementsCountAroundLpvOstium) ]
        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            lpvox, lpvod1, lpvod2, lpvod3, lpvoNodeId, None,
            lpvax, lpvad1, lpvad2, lpvad3, lpvaNodeId, lpvaDerivativesMap,
            maxEndThickness=laVenousFreeWallThickness,
            elementsCountRadial = elementsCountRadialPVAnnuli, meshGroups = [ heartMeshGroup, lamMeshGroup ],
            tracksurface=laTrackSurface, startProportions=lpvoProportions, endProportions=lpvaProportions,
            rescaleStartDerivatives=True, rescaleEndDerivatives=True, sampleBlend=0.0, fixMinimumStart=True)

        # left atrium epicardium venous midpoint marker point
        if commonLeftRightPvOstium:
            laevmElementId = elementIdentifier - elementsCountAroundLpvOstium*elementsCountRadialPVAnnuli + (elementsCountAroundLpvOstium // 2);
            if lpvOstiumSettings['Number of vessels'] == 3:
                laevmElementId += 1
            laevmXi = [ 0.0, 0.0, 1.0 ]
        else:
            laevmElementId = elementIdentifier - elementsCountAroundLpvOstium
            laevmXi = [ 0.0, 1.0, 1.0 ]
        laevmElement = mesh.findElementByIdentifier(laevmElementId)
        markerNode = laeVenousMidpointGroup.createMarkerNode(nodeIdentifier, element=laevmElement, xi=laevmXi)
        nodeIdentifier = markerNode.getIdentifier() + 1
        for group in [ heartGroup, lamGroup ]:
            group.getNodesetGroup(nodes).addNode(markerNode)

        if not commonLeftRightPvOstium:
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
                    rpvaDerivativesMap[1][ix] = ( (-1, 0, 0), (-1, -1, -1), (-1, 0, 1), (0, 1, 0 ) )
                elif n1 == elementsCountOverLeftAtriumVenous:
                    rpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 1, 0), (1, 0, 1), (1, 0, 0 ) )
                    rpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 1, -1), (-1, 0, 1), (1, 0, 0 ) )
                else:
                    rpvaDerivativesMap[0][ix] = ( (0, 1, 0), (-1, 0, 0), (1, 0, 1) )
                    rpvaDerivativesMap[1][ix] = ( (0, 1, 0), (-1, 0, -1), (-1, 0, 1) )
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
            # up left atrium venous midline, including both corners
            for n1 in range(elementsCountOverLeftAtriumVenous + 1):
                nm = elementsCountOverLeftAtriumVenous + 1 - n1
                for n3 in range(2):
                    rpvax [n3][ix] = lamlx [n3][nm]
                    rpvad1[n3][ix] = lamld1[n3][nm]
                    rpvad2[n3][ix] = lamld2[n3][nm]
                    rpvad3[n3][ix] = lamld3[n3][nm]
                    rpvaNodeId[n3][ix] = lamlNodeId[n3][nm]
                    if n1 == 0:
                        rpvaDerivativesMap[n3][ix] = ( (1, 0, 0), (1, 1, 0), None, (0, -1, 0) )
                    elif n1 == elementsCountOverLeftAtriumVenous:
                        rpvaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, -1, 0), None, (-1, 0, 0) )
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
            rpvoProportions = [ list(laTrackSurface.getProportion(rpvoPositions[i])) for i in range(elementsCountAroundRpvOstium) ]
            rpvaProportions = [ list(laTrackSurface.getProportion(laTrackSurface.findNearestPosition(rpvax[1][i], startPosition=rpvoPositions[i])))
                                for i in range(elementsCountAroundRpvOstium) ]
            nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
                nodes, mesh, nodeIdentifier, elementIdentifier,
                rpvox, rpvod1, rpvod2, rpvod3, rpvoNodeId, None,
                rpvax, rpvad1, rpvad2, rpvad3, rpvaNodeId, rpvaDerivativesMap,
                maxEndThickness=laVenousFreeWallThickness,
                elementsCountRadial = elementsCountRadialPVAnnuli, meshGroups = [ heartMeshGroup, lamMeshGroup ],
                tracksurface=laTrackSurface, startProportions=rpvoProportions, endProportions=rpvaProportions,
                rescaleStartDerivatives=True, rescaleEndDerivatives=True, sampleBlend=0.0)

        # create inferior and superior vena cavae inlets
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
            vcx, vd1, vd2, vd3 = getCircleProjectionAxes(ocx, ocd1, ocd2, ocd3, vcLength, vcAngle1Radians, vcAngle2Radians)
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
            # vc annuli
            vcax  = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad1 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad2 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcad3 = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcaNodeId = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            vcaDerivativesMap = [ [ None ]*elementsCountAroundVC, [ None ]*elementsCountAroundVC ]
            if v == 0:  # ivc
                # set points clockwise starting above crux of heart
                ix = 0
                # up interatrial groove to venous midpoint, including both corners
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
                        # on vestibule top, d1 and d2 are reversed
                        vcaDerivativesMap[0][ix] = ( (-1, 0, 0), (-1, -1, 0), (1, 0, 1), (0, 1, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (-1, 0, 0), (-1, -1, 0), (-1, 0, 1), (0, 1, 0 ) )
                    elif n1 == 0:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, -1, 0), (-1, 0, 1), (-1, 0, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, -1, -1), (1, 0, 1), (-1, 0, 1 ) )
                    else:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, 0, 0), (-1, 0, 1) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, 0, -1), (1, 0, 1) )
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
                # back over crista terminalis to vestibule top including corners
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
                            # on vestibule top, d1 and d2 are reversed
                            vcaDerivativesMap[n3][ix] = ( (0, -1, 0), (1, -1, 0), None, (-1, 0, 0) )
                        else:
                            vcaDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
                    ix += 1
                # left around vestibule top on venous posterior
                for n1 in range(1, elementsCountAroundRightAtriumPosteriorVenous):
                    nc = elementsCountAroundRightAtriumPosteriorVenous - n1
                    #print('v',v,'ix', ix, 'n1', n1, 'nc', nc)
                    for n3 in range(2):
                        vcax [n3][ix] = ravtx [n3][nc]
                        vcad1[n3][ix] = ravtd1[n3][nc]
                        vcad2[n3][ix] = ravtd2[n3][nc]
                        vcad3[n3][ix] = ravtd3[n3][nc]
                        vcaNodeId[n3][ix] = ravtNodeId[n3][nc]
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
                        vcaDerivativesMap[1][ix] = ( (1, 0, -1), (1, 1, -1), (1, 0, 1), (0, -1, 0 ) )
                    elif n1 == 0:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, -1, 0), (-1, 0, 1), (-1, 0, 0 ) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, 0 ) )
                    else:
                        vcaDerivativesMap[0][ix] = ( (0, -1, 0), (1, 0, 0), (-1, 0, 1) )
                        vcaDerivativesMap[1][ix] = ( (0, -1, 0), (1, 0, -1), (1, 0, 1) )
                    ix += 1
                # back over crista terminalis to vestibule top including corners
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

            # add rows of VC inlets to appropriate mesh groups
            rowMeshGroups = []
            vcMeshGroup = ivcMeshGroup if (v == 0) else svcMeshGroup
            vcInletMeshGroup = ivcInletMeshGroup if (v == 0) else svcInletMeshGroup
            if elementsCountAlongVCInlet == 1:
                rowMeshGroups = [ [ heartMeshGroup, vcMeshGroup, vcInletMeshGroup, ramMeshGroup] ]
            else:
                rowMeshGroups = []
                for i in range(elementsCountAlongVCInlet):
                    xi = (i + 1)/elementsCountAlongVCInlet
                    meshGroups = [ heartMeshGroup ]
                    if xi < 0.67:
                        meshGroups +=  [ vcMeshGroup ]
                    if xi > 0.51:
                        meshGroups +=  [ vcInletMeshGroup, ramMeshGroup ]
                    rowMeshGroups.append(meshGroups)
            nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
                nodes, mesh, nodeIdentifier, elementIdentifier,
                vcvx, vcvd1, vcvd2, None, None, None,
                vcax, vcad1, vcad2, vcad3, vcaNodeId, vcaDerivativesMap,
                maxEndThickness = 1.5*raVenousFreeWallThickness,
                elementsCountRadial = elementsCountAlongVCInlet,
                meshGroups = rowMeshGroups, rescaleEndDerivatives=True, fixMinimumStart=True)

            if v == 0:  # ivc
                # right atrium epicardium venous midpoint marker point
                raevmElementId = elementIdentifier - elementsCountAroundVC + elementsCountOverRightAtriumVenous//2 + elementsCountOverSideRightAtriumVC//2
                raevmXi = [ 0.5 if (elementsCountOverSideRightAtriumVC % 2) else 0.0, 1.0, 1.0 ]
                raevmElement = mesh.findElementByIdentifier(raevmElementId)
                markerNode = raeVenousMidpointGroup.createMarkerNode(nodeIdentifier, element=raevmElement, xi=raevmXi)
                nodeIdentifier = markerNode.getIdentifier() + 1
                for group in [heartGroup, ramGroup]:
                    group.getNodesetGroup(nodes).addNode(markerNode)

        # create left atrial appendage
        position = laTrackSurface.createPositionProportion(laaMidpointOver, laaMidpointLeft)
        laamx, d1, d2 = laTrackSurface.evaluateCoordinates(position, derivatives = True)
        # force d2 to be vertical, d3, d1 to be horizontal
        laamd2 = [ 0.0, 0.0, 1.0 ]
        laamd3 = vector.normalise(vector.crossproduct3(d2, laamd2))
        laamd1 = vector.crossproduct3(laamd2, laamd3)
        if False:
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, laamx )
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, laamd1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, laamd2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, laamd3)
            nodeIdentifier += 1
        elementsCountAroundLaa = elementsCountAroundLeftAtrialAppendageBase + len(laoaProportions) + 1
        # get start points, nodes, derivative maps
        laasx  = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laasd1 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laasd2 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laasd3 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laasNodeId = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laasDerivativesMap = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        # set points anticlockwise around base first, starting at aorta
        ixStart = 1 + elementsCountAroundLaa%2  # position in final array for aorta base
        ixRotation = round(elementsCountAroundLaa*0.5*laaAngleAxialRadians/math.pi)  # rotate indexes to align with axial angle
        ix = (ixStart - ixRotation) % elementsCountAroundLaa - elementsCountAroundLaa  # works for negative values as modulo is always non-negative in python
        #print('laa ixStart',ixStart,'ixRotation',ixRotation,'ix',ix)
        if commonLeftRightPvOstium:
            nv = 1
            for n3 in range(2):
                laasx [n3][ix] = lavtx [n3][nv]
                laasd1[n3][ix] = lavtd1[n3][nv]
                laasd2[n3][ix] = lavtd2[n3][nv]
                laasd3[n3][ix] = lavtd3[n3][nv]
                laasNodeId[n3][ix] = lavtNodeId[n3][nv]
                laasDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 0, 0), None)
            ix += 1
        # left along base
        for n1 in range(elementsCountAroundLeftAtrialAppendageBase + 1):
            nb = n1 + elementsCountAroundLeftAtriumAorta
            for n3 in range(2):
                laasx [n3][ix] = labx [n3][nb]
                laasd1[n3][ix] = labd1[n3][nb]
                laasd2[n3][ix] = labd2[n3][nb]
                laasd3[n3][ix] = labd3[n3][nb]
                laasNodeId[n3][ix] = labNodeId[n3][nb]
                if n1 == 0:
                    laasDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 1, 0), None, None )
                elif n1 == elementsCountAroundLeftAtrialAppendageBase:
                    laasDerivativesMap[n3][ix] = ( None, (-1, 1, 0), None, (0, 1, 0 ) )
                else:
                    laasDerivativesMap[n3][ix] = ( None, None, None )
            ix += 1
        # left over appendage laoa
        for n1 in range(laoaCount):
            no = -1 - n1
            for n3 in range(2):
                laasx [n3][ix] = laoax [n3][no]
                laasd1[n3][ix] = laoad1[n3][no]
                laasd2[n3][ix] = laoad2[n3][no]
                laasd3[n3][ix] = laoad3[n3][no]
                laasNodeId[n3][ix] = laoaNodeId[n3][no]
                if n1 == 0:
                    laasDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, 0, 0), None )
                elif n1 == (laoaCount - 1):
                    laasDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 0, 0), None )
                else:
                    laasDerivativesMap[n3][ix] = ( (-1, 0, 0), (0, -1, 0), None )
            ix += 1
        #print('laasNodeId[0]',laasNodeId[0])
        #print('laasNodeId[1]',laasNodeId[1])
        #print('laasDerivativesMap[0]',laasDerivativesMap[0])
        #print('laasDerivativesMap[1]',laasDerivativesMap[1])
        # get end points, nodes, derivative maps, expanding from wedge
        laawx, laawd1, laawd2, laawd3, elementsCountAcrossLaaWedge, laawPointsMap, laaeDerivativesMap = \
            getAtrialAppendageWedgePoints(laamx, laamd1, laamd2, laamd3, laaAngleLeftRadians, laaAngleUpradians, laaAngleAxialRadians, laaBaseLength,
                elementsCountAroundLaa, elementsCountAlongAtrialAppendages, laaArcLength, laaArcRadius, laaWallThickness, laaWedgeAngleRadians)
        # create laa wedge nodes:
        laawNodeId = [ [], [] ]
        for n3 in range(2):
            for n1 in range(len(laawx[n3])):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, laawx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, laawd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, laawd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, laawd3[n3][n1])
                laawNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
        laaex  = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laaed1 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laaed2 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laaed3 = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        laaeNodeId = [ [ None ]*elementsCountAroundLaa, [ None ]*elementsCountAroundLaa ]
        for n1 in range(elementsCountAroundLaa):
            nw = laawPointsMap[n1]
            for n3 in range(2):
                laaex [n3][n1] = laawx [n3][nw]
                laaed1[n3][n1] = laawd1[n3][nw]
                laaed2[n3][n1] = laawd2[n3][nw]
                laaed3[n3][n1] = laawd3[n3][nw]
                laaeNodeId[n3][n1] = laawNodeId[n3][nw]
        #print('laaeNodeId[0]',laaeNodeId[0])
        #print('laaeNodeId[1]',laaeNodeId[1])
        #print('laaeDerivativesMap[0]',laaeDerivativesMap[0])
        #print('laaeDerivativesMap[1]',laaeDerivativesMap[1])
        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            laasx, laasd1, laasd2, laasd3, laasNodeId, laasDerivativesMap,
            laaex, laaed1, laaed2, laaed3, laaeNodeId, laaeDerivativesMap,
            forceMidLinearXi3 = True, forceEndLinearXi3 = True,
            maxStartThickness = laaWallThickness,
            elementsCountRadial = elementsCountAlongAtrialAppendages,
            meshGroups = [ heartMeshGroup, lamMeshGroup, laaMeshGroup ])

        # create right atrium plain elements
        # Anderson considers these part of the right atrial appendage:
        #meshGroups = [ ramMeshGroup, raaMeshGroup ]
        meshGroups = [ heartMeshGroup, ramMeshGroup ]
        for e2 in range(2):
            for e1 in range(elementsCountOverSideRightAtriumPouch):
                eft1 = eft
                elementtemplate1 = elementtemplate
                scalefactors = None

                nc = 2 + e1
                if e2 == 0:
                    nids = [ ractNodeId[0][nc], ractNodeId[0][nc + 1], raaqNodeId[0][e1], raaqNodeId[0][e1 + 1],
                             ractNodeId[1][nc], ractNodeId[1][nc + 1], raaqNodeId[1][e1], raaqNodeId[1][e1 + 1] ]
                else:
                    nids = [ raaqNodeId[0][e1], raaqNodeId[0][e1 + 1], raapNodeId[0][e1], raapNodeId[0][e1 + 1],
                             raaqNodeId[1][e1], raaqNodeId[1][e1 + 1], raapNodeId[1][e1], raapNodeId[1][e1 + 1] ]
                if e1 == 0:
                    if e2 == 0:
                        continue
                    nids[0] = ractNodeId[0][nc + 1]
                    nids[4] = ractNodeId[1][nc + 1]
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 2, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    #remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                    remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                elif (e2 == 1) and (e1 < (elementsCountOverSideRightAtriumPouch - 1)):
                    pass  # regular elements
                else:
                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                    setEftScaleFactorIds(eft1, [1], [])
                    scalefactors = [ -1.0 ]
                    if e1 < (elementsCountOverSideRightAtriumPouch - 1):
                        if e2 == 0:
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 2, 5, 6 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                            if e1 == 1:
                                remapEftNodeValueLabel(eft1, [ 3, 7 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
                    else:
                        if e2 == 0:
                            remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                            remapEftNodeValueLabel(eft1, [ 1, 5 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [] ) ])
                        remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [1] ) ])
                        remapEftNodeValueLabel(eft1, [ 2, 4, 6, 8 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
                if eft1 is not eft:
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                #print('create element raa plain', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

        # create right atrial appendage 'pouch'
        position = raTrackSurface.createPositionProportion(1.0 - raaMidpointOver, raaMidpointRight)
        raamx, d1, d2 = raTrackSurface.evaluateCoordinates(position, derivatives = True)
        # force d2 to be vertical, d3, d1 to be horizontal
        raamd2 = [ 0.0, 0.0, 1.0 ]
        raamd3 = vector.normalise(vector.crossproduct3(raamd2, d2))
        raamd1 = vector.crossproduct3(raamd2, raamd3)
        if False:
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, raamx )
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, raamd1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, raamd2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, raamd3)
            nodeIdentifier += 1
        elementsCountAroundRaa = elementsCountAroundRightAtrialAppendagePlainBase + elementsCountAroundRightAtrialAppendagePouchBase + elementsCountOverAtria - 2
        #print('elementsCountAroundRaa', elementsCountAroundRaa)
        # get start points, nodes, derivative maps
        raasx  = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raasd1 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raasd2 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raasd3 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raasNodeId = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raasDerivativesMap = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        # set points anticlockwise around base first from aorta
        ixStart = 2 + elementsCountAroundRaa % 2  # position in final array for first base on raap
        ixRotation = round(elementsCountAroundRaa*0.5*raaAngleAxialRadians/math.pi)  # rotate indexes to align with axial angle
        ix = (ixStart - ixRotation) % elementsCountAroundRaa - elementsCountAroundRaa  # works for negative values as modulo is always non-negative in python
        #print('raa ixStart',ixStart,'ixRotation',ixRotation,'ix',ix)
        # left/anticlockwise along base
        for n1 in range(elementsCountAroundRightAtrialAppendagePouchBase + 1):
            nb = elementsCountAroundRightAtriumPosteriorVenous + elementsCountAroundRightAtrialAppendagePlainBase + n1
            for n3 in range(2):
                raasx [n3][ix] = rabx [n3][nb]
                raasd1[n3][ix] = rabd1[n3][nb]
                raasd2[n3][ix] = rabd2[n3][nb]
                raasd3[n3][ix] = rabd3[n3][nb]
                raasNodeId[n3][ix] = rabNodeId[n3][nb]
                if n1 == 0:
                    raasDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 1, 0), None, None )
                elif n1 == elementsCountAroundRightAtrialAppendagePouchBase:
                    raasDerivativesMap[n3][ix] = ( None, (-1, 1, 0), None, (0, 1, 0 ) )
                else:
                    raasDerivativesMap[n3][ix] = ( None, None, None )
            ix += 1
        # back over crista terminalis ract
        for n1 in range(2):
            nc = 1 + n1
            for n3 in range(2):
                raasx [n3][ix] = ractx [n3][nc]
                raasd1[n3][ix] = ractd1[n3][nc]
                raasd2[n3][ix] = ractd2[n3][nc]
                raasd3[n3][ix] = ractd3[n3][nc]
                raasNodeId[n3][ix] = ractNodeId[n3][nc]
                if n1 == 0:
                    raasDerivativesMap[n3][ix] = ( None, (-1, 0, 0), None, (0, 1, 0) )  # compare ( None, (-1, 1, 0), None, (0, 1, 0) )
                else: # n1 == 1:
                    #raasDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, -1, 0), None, (-1, 1, 0) )  # compare ( (0, 1, 0), (-1, 0, 0), None, (-1, 1, 0) )
                    raasDerivativesMap[n3][ix] = ( (0, 1, 0), (-1, -1, 0), None, (-1, 0, 0) )
            ix += 1
        # down right atrial appendage pouch limit, raap
        n1Last = len(raapx[1]) - 1
        for n1 in range(1, n1Last + 1):
            for n3 in range(2):
                raasx [n3][ix] = raapx [n3][n1]
                raasd1[n3][ix] = raapd1[n3][n1]
                raasd2[n3][ix] = raapd2[n3][n1]
                raasd3[n3][ix] = raapd3[n3][n1]
                raasNodeId[n3][ix] = raapNodeId[n3][n1]
                if n1 == n1Last:
                    raasDerivativesMap[n3][ix] = ( (0, -1, 0), (1, 0, 0), None )
                else:
                    raasDerivativesMap[n3][ix] = ( (1, 0, 0), (0, 1, 0), None )
            ix += 1
        #print('raasNodeId[0]',raasNodeId[0])
        #print('raasNodeId[1]',raasNodeId[1])
        #print('raasDerivativesMap[0]',raasDerivativesMap[0])
        #print('raasDerivativesMap[1]',raasDerivativesMap[1])
        # get end points, nodes, derivative maps, expanding from wedge
        raawx, raawd1, raawd2, raawd3, elementsCountAcrossRaaWedge, raawPointsMap, raaeDerivativesMap = \
            getAtrialAppendageWedgePoints(raamx, raamd1, raamd2, raamd3, raaAngleLeftRadians, raaAngleUpradians, raaAngleAxialRadians, raaBaseLength,
                elementsCountAroundRaa, elementsCountAlongAtrialAppendages, raaArcLength, raaArcRadius, raaWallThickness, raaWedgeAngleRadians)
        # create raa wedge nodes:
        raawNodeId = [ [], [] ]
        for n3 in range(2):
            for n1 in range(len(raawx[n3])):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, raawx [n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, raawd1[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, raawd2[n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, raawd3[n3][n1])
                raawNodeId[n3].append(nodeIdentifier)
                nodeIdentifier += 1
        raaex  = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raaed1 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raaed2 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raaed3 = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        raaeNodeId = [ [ None ]*elementsCountAroundRaa, [ None ]*elementsCountAroundRaa ]
        for n1 in range(elementsCountAroundRaa):
            nw = raawPointsMap[n1]
            for n3 in range(2):
                raaex [n3][n1] = raawx [n3][nw]
                raaed1[n3][n1] = raawd1[n3][nw]
                raaed2[n3][n1] = raawd2[n3][nw]
                raaed3[n3][n1] = raawd3[n3][nw]
                raaeNodeId[n3][n1] = raawNodeId[n3][nw]
        #print('raaeNodeId[0]',raaeNodeId[0])
        #print('raaeNodeId[1]',raaeNodeId[1])
        #print('raaeDerivativesMap[0]',raaeDerivativesMap[0])
        #print('raaeDerivativesMap[1]',raaeDerivativesMap[1])
        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            raasx, raasd1, raasd2, raasd3, raasNodeId, raasDerivativesMap,
            raaex, raaed1, raaed2, raaed3, raaeNodeId, raaeDerivativesMap,
            forceMidLinearXi3 = True, forceEndLinearXi3 = True,
            maxStartThickness = raaWallThickness,
            elementsCountRadial = elementsCountAlongAtrialAppendages,
            meshGroups = [ heartMeshGroup, ramMeshGroup, raaMeshGroup ])

        if defineEpicardiumLayer:
            # project epicardial points over atria to build fat pad
            epiGroup = fm.createFieldGroup()
            epiMesh = epiGroup.createFieldElementGroup(mesh).getMeshGroup()
            is_a = fm.createFieldOr(lamGroup.getFieldElementGroup(mesh), ramGroup.getFieldElementGroup(mesh))
            is_aa = fm.createFieldOr(laaGroup.getFieldElementGroup(mesh), raaGroup.getFieldElementGroup(mesh))
            is_not_epi = fm.createFieldNot(fm.createFieldOr(is_aa, aSeptumGroup.getFieldElementGroup(mesh)))
            is_a_epi = fm.createFieldAnd(is_a, is_not_epi)
            epiMesh.addElementsConditional(is_a_epi)
            # print("epiMesh.getSize()", epiMesh.getSize())
            epiNodes = epiGroup.createFieldNodeGroup(nodes).getNodesetGroup()
            # add nodes on xi3=1 of epiMesh to epiNodes
            epiElementIdentifiers = []
            elementIterator = epiMesh.createElementiterator()
            epiElement = elementIterator.next()
            while epiElement.isValid():
                epiElementIdentifier = epiElement.getIdentifier()
                epiEft = epiElement.getElementfieldtemplate(coordinates, -1)
                # only implemented for regular cube elements
                nodeCount = epiEft.getNumberOfLocalNodes()
                if nodeCount == 8:
                    for ln in range(5, nodeCount + 1):
                        epiNode = epiElement.getNode(epiEft, ln)
                        epiNodes.addNode(epiNode)
                    epiElementIdentifiers.append(epiElementIdentifier)
                else:
                    print("Non-cube element ID on epicardium", epiElementIdentifier, "#nodes", nodeCount)
                epiElement = elementIterator.next()
            # print("epiNodes.getSize()", epiNodes.getSize())
            # project nodes above epicardial surface or nearest point on fatpad tracksurface if further away
            epiFatPadNodeIdentifiersMap = {}
            # make blank map first since can't iterate over nodes while creating them
            nodeIterator = epiNodes.createNodeiterator()
            epiNode = nodeIterator.next()
            while epiNode.isValid():
                epiNodeIdentifier = epiNode.getIdentifier()
                epiFatPadNodeIdentifiersMap[epiNodeIdentifier] = epiNodeIdentifier
                epiNode = nodeIterator.next()
            fpGroup = fm.createFieldGroup()
            fpGroup.setName("fp")
            fpNodes = fpGroup.createFieldNodeGroup(nodes).getNodesetGroup()
            bridgeGroup = fm.createFieldGroup()
            bridgeGroup.setName("ia_bridge")
            bridgeNodes = bridgeGroup.createFieldNodeGroup(nodes).getNodesetGroup()
            bridgeNodeTangents = {}

            # parameters for shifting centre of fatpad to spread elements around PV, VC inlets
            # have ivcPositionOver, svcPositionOver
            rpvPositionOver = lpvOstiumPositionOver if commonLeftRightPvOstium else rpvOstiumPositionOver
            # GRC fudge factors
            svcSpread = ivcSpread = ivcPositionOver
            rpvSpread = ivcSpread
            ivcMag = 0.15
            svcMag = 0.0
            rpvMag = 0.15

            for n in range(elementsCountAcrossTrackSurface + 1):
                xi_ivc = max(1.0, math.fabs(ivcPositionOver - nProportions[n][1]))
                xi_svc = max(1.0, math.fabs(ivcPositionOver - nProportions[n][1]))

            for epiNodeIdentifier in epiFatPadNodeIdentifiersMap.keys():
                epiNode = nodes.findNodeByIdentifier(epiNodeIdentifier)
                cache.setNode(epiNode)
                result, epix = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                result, epid1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                result, epid2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                result, epid3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                if result != RESULT_OK:
                    epid3 = vector.crossproduct3(epid1, epid2)
                fatx = add(epix, vector.setMagnitude(epid3, epicardiumLayerMinimumThickness))
                epifx = None

                epifPosition = fpTrackSurface.findNearestPosition(epix, startPosition=None)
                if not (((epifPosition.e1 == 0) and (epifPosition.xi1 < 0.0001)) or
                        ((epifPosition.e1 == (fpTrackSurface.elementsCount1 - 1)) and (0.9999 < epifPosition.xi1))):
                    # shift proportion 1 around inlet
                    xi_centre = 1.0 - (2.0 * math.fabs(0.5 - epifPosition.xi1))
                    proportion1, proportion2 = fpTrackSurface.getProportion(epifPosition)
                    # scale to spread out centre
                    # proportion1Scaled = interp.interpolateCubicHermite([0.0], [0.5], [0.5], [0.6], xi_centre)[0]
                    # if proportion1 < 0.5:
                    #     proportion1 = proportion1Scaled
                    # else:
                    #     proportion1 = 1.0 - proportion1Scaled
                    proportion1Shift = 0.0
                    xi_ivc = math.fabs((ivcPositionOver - proportion2) / ivcSpread)
                    if xi_ivc < 1.0:
                        proportion1Shift -= interp.interpolateCubicHermite([ivcMag], [0.0], [0.0], [0.0], xi_ivc)[0]
                    xi_svc = math.fabs((svcPositionOver - proportion2) / svcSpread)
                    if xi_svc < 1.0:
                        proportion1Shift -= interp.interpolateCubicHermite([svcMag], [0.0], [0.0], [0.0], xi_svc)[0]
                    xi_rpv = math.fabs((rpvPositionOver - proportion2) / rpvSpread)
                    if xi_rpv < 1.0:
                        proportion1Shift += interp.interpolateCubicHermite([rpvMag], [0.0], [0.0], [0.0], xi_rpv)[0]
                    proportion1Shift *= xi_centre
                    # proportion1Shift = interp.interpolateCubicHermite(
                    #     [0.0], [0.0], [proportion1Shift], [0.0], xi_centre)[0]
                    proportion1 += proportion1Shift

                    # convert shifted proportions to shifted position on fpTrackSurface
                    epifPosition2 = fpTrackSurface.createPositionProportion(proportion1, proportion2)
                    epifx, epifd1, epifd2 = fpTrackSurface.evaluateCoordinates(epifPosition2, derivatives=True)
                    delta_epi = sub(epifx, epix)
                    # epifx must be above the epicardium surface
                    # and at least epicardiumLayerMinimumThickness away from epix
                    if (dot(delta_epi, epid3) > 0.0) and (magnitude(delta_epi) >= epicardiumLayerMinimumThickness):
                        epifNormal = normalize(cross(epifd1, epifd2))
                        # epix must be under the fatpad
                        if dot(delta_epi, epifNormal) > 0.0:
                            fatx = epifx

                node = nodes.createNode(nodeIdentifier, nodetemplateLinearS3)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, fatx)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, epid1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, epid2)
                # interatrial groove nodes at crux and vestibule top need to keep their original values
                if epiNodeIdentifier not in agNodeId[-2:] + [ractNodeId[1][-1]]:
                    fpNodes.addNode(node)
                    if fatx is epifx:
                        bridgeNodes.addNode(node)
                        bridgeNodeTangents[nodeIdentifier] = (epifd1, epifd2, epifNormal)
                epiFatPadNodeIdentifiersMap[epiNodeIdentifier] = nodeIdentifier
                nodeIdentifier += 1

            # create fatpad elements
            fatEft = bicubichermitelinear.createEftNoCrossDerivatives()
            fatElementtemplate = mesh.createElementtemplate()
            fatElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            fatElementtemplate.defineField(coordinates, -1, fatEft)
            elementtemplateX = mesh.createElementtemplate()
            elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

            epicardiumMeshGroup = epicardiumGroup.getMeshGroup(mesh)
            meshGroups = [heartMeshGroup, epicardiumMeshGroup]
            elementtemplate = mesh.createElementtemplate()
            # iterate over list of identifiers since can't iterate over mesh while modifying it
            for epiElementIdentifier in epiElementIdentifiers:
                epiElement = mesh.findElementByIdentifier(epiElementIdentifier)
                epiEft = epiElement.getElementfieldtemplate(coordinates, -1)
                elementtemplate1 = fatElementtemplate
                eft1, scalefactors = createEftElementSurfaceLayer(epiElement, epiEft, bicubichermitelinear, fatEft,
                                                                  removeNodeValueLabel=Node.VALUE_LABEL_D_DS3)
                nids = []
                for ln in range(5, nodeCount + 1):
                    epiNode = epiElement.getNode(epiEft, ln)
                    nids.append(epiNode.getIdentifier())
                nids += [epiFatPadNodeIdentifiersMap[nid] for nid in nids]
                if eft1 is not fatEft:
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result2 = element.setNodesByIdentifier(eft1, nids)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                # print('create element fat', element.isValid(), elementIdentifier, result2, result3, nids)
                elementIdentifier += 1

                for meshGroup in meshGroups:
                    meshGroup.addElement(element)

            # smooth bridge nodes
            smoothing = DerivativeSmoothing(region, coordinates, selectionGroupName="ia_bridge",
                                            scalingMode=DerivativeScalingMode.HARMONIC_MEAN)
            smoothing.smooth(updateDirections=True)
            #project derivatives onto fpTrackSurface
            for nid, tangents in bridgeNodeTangents.items():
                sd1, sd2, sd3 = tangents
                node = nodes.findNodeByIdentifier(nid)
                cache.setNode(node)
                result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                newd1 = sub(d1, mult(sd3, dot(d1, sd3)))
                newd2 = sub(d2, mult(sd3, dot(d2, sd3)))
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, newd1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, newd2)
            # smooth over all new nodes without changing directions
            smoothing = DerivativeSmoothing(region, coordinates, selectionGroupName="fp",
                                            scalingMode=DerivativeScalingMode.HARMONIC_MEAN)
            smoothing.smooth(updateDirections=False)

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
        if drawFpTrackSurface:
            mesh2d = fm.findMeshByDimension(2)
            bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft2d = mesh2d.createElementfieldtemplate(bicubicHermiteBasis)
            # remove cross derivative 12
            for n in range(4):
                r = eft2d.setFunctionNumberOfTerms(n * 4 + 4, 0)
            elementtemplate2d = mesh2d.createElementtemplate()
            elementtemplate2d.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            elementtemplate2d.defineField(coordinates, -1, eft2d)
            nodesCount1 = fpTrackSurface.elementsCount1 + 1
            for e2 in range(fpTrackSurface.elementsCount2):
                for e1 in range(fpTrackSurface.elementsCount1):
                    element = mesh2d.createElement(-1, elementtemplate2d)  # since on 2-D mesh
                    nid1 = fpTrackSurfaceFirstNodeIdentifier + e2 * nodesCount1 + e1
                    element.setNodesByIdentifier(eft2d, [nid1, nid1 + 1, nid1 + nodesCount1, nid1 + nodesCount1 + 1])

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
        refineElementsCountThroughEpicardiumLayer =\
            options['Refine number of elements through epicardium layer']
        defineEpicardiumLayer = options['Define epicardium layer']

        sourceFm = meshrefinement._sourceFm
        annotationGroups = meshrefinement._sourceAnnotationGroups
        laGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left atrium myocardium"))
        laElementGroupField = laGroup.getFieldElementGroup(meshrefinement._sourceMesh)
        raGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right atrium myocardium"))
        raElementGroupField = raGroup.getFieldElementGroup(meshrefinement._sourceMesh)
        aSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interatrial septum"))
        aSeptumMeshGroup = aSeptumGroup.getMeshGroup(meshrefinement._sourceMesh)
        epicardiumGroup = None
        epicardiumMeshGroup = None
        if defineEpicardiumLayer:
            epicardiumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("epicardium"))
            epicardiumMeshGroup = epicardiumGroup.getMeshGroup(meshrefinement._sourceMesh)
        coordinates = findOrCreateFieldCoordinates(meshrefinement._sourceFm)

        # last atria element is last element in following group:
        lastGroup = epicardiumGroup
        if not lastGroup:
            lastGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right auricle"))
        lastMeshGroup = lastGroup.getMeshGroup(meshrefinement._sourceMesh)
        lastElementIdentifier = -1
        elementIter = lastMeshGroup.createElementiterator()
        element = elementIter.next()
        while element.isValid():
            lastElementIdentifier = element.getIdentifier()
            element = elementIter.next()

        cache = sourceFm.createFieldcache()
        element = meshrefinement._sourceElementiterator.next()
        wedgeElementCount = 0
        while element.isValid():
            elementIdentifier = element.getIdentifier()
            refineElements1 = refineElementsCountSurface
            refineElements2 = refineElementsCountSurface
            refineElements3 = refineElementsCountThroughWall
            cache.setElement(element)
            if aSeptumMeshGroup.containsElement(element):
                eft = element.getElementfieldtemplate(coordinates, 1)
                if eft.getNumberOfLocalNodes() == 6:
                    wedgeElementCount += 1
                    # the first two around the base are collapsed on 1-3, remainder on 2-3
                    if wedgeElementCount <= 2:
                        refineElements1 = refineElementsCountThroughWall
                    else:
                        refineElements2 = refineElementsCountThroughWall
            elif epicardiumGroup and epicardiumMeshGroup.containsElement(element):
                refineElements3 = refineElementsCountThroughEpicardiumLayer
            meshrefinement.refineElementCubeStandard3d(element, refineElements1, refineElements2, refineElements3)
            if elementIdentifier == lastElementIdentifier:
                return  # finish on last so can continue elsewhere
            element = meshrefinement._sourceElementiterator.next()


    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        # create endocardium and epicardium groups
        defineEpicardiumLayer = options['Define epicardium layer']
        fm = region.getFieldmodule()
        lamGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left atrium myocardium"))
        ramGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right atrium myocardium"))
        aSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interatrial septum"))
        laaGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left auricle"))
        raaGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right auricle"))
        lpvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left pulmonary vein"))
        # middle pulmonary vein is only present in rodents:
        mpvGroup = findAnnotationGroupByName(annotationGroups, "middle pulmonary vein")
        rpvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right pulmonary vein"))
        ivcGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("inferior vena cava"))
        svcGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("superior vena cava"))
        # following will already be defined if defineEpicardiumLayer is true
        epiGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("epicardium"))

        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_lam = lamGroup.getFieldElementGroup(mesh2d)
        is_ram = ramGroup.getFieldElementGroup(mesh2d)
        is_lam_endo = fm.createFieldAnd(is_lam, is_exterior_face_xi3_0)
        is_ram_endo = fm.createFieldOr(fm.createFieldAnd(fm.createFieldAnd(is_ram, is_exterior_face_xi3_0),
                                                        fm.createFieldNot(is_lam_endo)),
                                      fm.createFieldAnd(aSeptumGroup.getFieldElementGroup(mesh2d),
                                                        is_exterior_face_xi3_1))
        is_laa = laaGroup.getFieldElementGroup(mesh2d)
        is_raa = raaGroup.getFieldElementGroup(mesh2d)
        is_laa_endo = fm.createFieldAnd(is_laa, is_exterior_face_xi3_0)
        is_raa_endo = fm.createFieldAnd(is_raa, is_exterior_face_xi3_0)
        is_laa_epi = fm.createFieldAnd(laaGroup.getFieldElementGroup(mesh2d), is_exterior_face_xi3_1)
        is_raa_epi = fm.createFieldAnd(raaGroup.getFieldElementGroup(mesh2d), is_exterior_face_xi3_1)
        # is_myocardium = fm.createFieldOr(is_lam, is_ram)
        is_ext_xi3_1_and_not_septum = fm.createFieldAnd(
            is_exterior_face_xi3_1, fm.createFieldNot(aSeptumGroup.getFieldElementGroup(mesh2d)))
        is_os_lam = fm.createFieldAnd(is_lam, is_ext_xi3_1_and_not_septum)
        is_os_ram = fm.createFieldAnd(is_ram, is_ext_xi3_1_and_not_septum)
        is_epi = epiGroup.getFieldElementGroup(mesh2d)

        # luminal surfaces of endocardium of left/right atrium
        lslaEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of left atrium"))
        lslaEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_lam_endo)
        lsraEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of right atrium"))
        lsraEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_ram_endo)
        # endocardium groups are defined identically to luminal surfaces at scaffold scale
        laEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("left atrium endocardium"))
        laEndoGroup.getMeshGroup(mesh2d).addElementsConditional(lslaEndoGroup.getFieldElementGroup(mesh2d))
        raEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("right atrium endocardium"))
        raEndoGroup.getMeshGroup(mesh2d).addElementsConditional(lsraEndoGroup.getFieldElementGroup(mesh2d))

        laaEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("endocardium of left auricle"))
        laaEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_laa_endo)
        raaEndoGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("endocardium of right auricle"))
        raaEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_raa_endo)
        laaEpiGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("epicardium of left auricle"))
        laaEpiGroup.getMeshGroup(mesh2d).addElementsConditional(is_laa_epi)
        raaEpiGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("epicardium of right auricle"))
        raaEpiGroup.getMeshGroup(mesh2d).addElementsConditional(is_raa_epi)

        oslamGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("outer surface of myocardium of left atrium"))
        oslamGroup.getMeshGroup(mesh2d).addElementsConditional(is_os_lam)
        osramGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("outer surface of myocardium of right atrium"))
        osramGroup.getMeshGroup(mesh2d).addElementsConditional(is_os_ram)
        if defineEpicardiumLayer:
            oslamGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(is_lam, is_epi))
            osramGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(is_ram, is_epi))
        osmGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("outer surface of myocardium"))
        osmGroup.getMeshGroup(mesh2d).addElementsConditional(oslamGroup.getFieldElementGroup(mesh2d))
        osmGroup.getMeshGroup(mesh2d).addElementsConditional(osramGroup.getFieldElementGroup(mesh2d))
        if defineEpicardiumLayer:
            # future: limit to atria once ventricles have epicardium layer
            is_os_epi = fm.createFieldAnd(is_epi, is_exterior_face_xi3_1)
            osEpiGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region,get_heart_term("outer surface of epicardium"))
            osEpiGroup.getMeshGroup(mesh2d).addElementsConditional(is_os_epi)
            # note: epiGroup only contains 3-D elements in this case
        else:
            # if no volumetric epicardium group, add outer surface of atrial myocardium
            epiGroup.getMeshGroup(mesh2d).addElementsConditional(oslamGroup.getFieldElementGroup(mesh2d))
            epiGroup.getMeshGroup(mesh2d).addElementsConditional(osramGroup.getFieldElementGroup(mesh2d))

        lslpvGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of left pulmonary vein"))
        lsrpvGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of right pulmonary vein"))
        is_lpv = lpvGroup.getFieldElementGroup(mesh2d)
        is_mpv = mpvGroup.getFieldElementGroup(mesh2d) if mpvGroup else None
        is_rpv = rpvGroup.getFieldElementGroup(mesh2d)
        lslpvGroup.getMeshGroup(mesh2d).addElementsConditional(
            fm.createFieldAnd(is_exterior_face_xi3_0, is_lpv))
        if mpvGroup:
            lsmpvGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_heart_term("luminal surface of middle pulmonary vein"))
            lsmpvGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldAnd(is_exterior_face_xi3_0, is_mpv))
        lsrpvGroup.getMeshGroup(mesh2d).addElementsConditional(
            fm.createFieldAnd(is_exterior_face_xi3_0, is_rpv))

        lsivcGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of inferior vena cava"))
        lssvcGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, get_heart_term("luminal surface of superior vena cava"))
        is_ivc = ivcGroup.getFieldElementGroup(mesh2d)
        is_svc = svcGroup.getFieldElementGroup(mesh2d)
        lsivcGroup.getMeshGroup(mesh2d).addElementsConditional(
            fm.createFieldAnd(is_exterior_face_xi3_0, is_ivc))
        lssvcGroup.getMeshGroup(mesh2d).addElementsConditional(
            fm.createFieldAnd(is_exterior_face_xi3_0, is_svc))

        lFibrousRingGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left fibrous ring"))
        rFibrousRingGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right fibrous ring"))
        if (lFibrousRingGroup.getDimension() <= 0) or (rFibrousRingGroup.getDimension() <= 0):
            is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
            is_pv = fm.createFieldOr(is_lpv, is_rpv)
            if mpvGroup:
                is_pv = fm.createFieldOr(is_pv, is_mpv)
            lFibrousRingGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldAnd(is_lam, fm.createFieldAnd(is_exterior_face_xi2_0, fm.createFieldNot(is_pv))))
            is_vc = fm.createFieldOr(is_ivc, is_svc)
            rFibrousRingGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldAnd(is_ram, fm.createFieldAnd(is_exterior_face_xi2_0, fm.createFieldNot(is_vc))))

def getLeftAtriumPulmonaryVeinOstiaElementsCounts(elementsCountAroundLeftAtriumFreeWall, elementsCountOverAtria, commonLeftRightPvOstium):
    '''
    Return numbers of elements around left and right pulmonary vein ostia.
    If commonLeftRightPvOstium, value is returned in elementsCountAroundLpvOstium.
    :return: elementsCountAroundLpvOstium, elementsCountAroundRpvOstium
    '''
    elementsCountOverAtriaCoronarySinus, \
    elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior, \
    elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior \
        = getOverAtriaElementsCounts(elementsCountOverAtria)
    elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftAtrialAppendageBase, elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV \
        = getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall)
    elementsCountAroundRpvOstium = 2*(elementsCountOverLeftAtriumVenous + elementsCountAroundLeftAtriumRPV)
    if commonLeftRightPvOstium:
        elementsCountAroundLpvOstium = elementsCountOverLeftAtriumVenous + 2*(elementsCountAroundLeftAtriumLPV + elementsCountAroundLeftAtriumRPV) + 2
    else:
        elementsCountAroundLpvOstium = elementsCountOverLeftAtriumVenous + 2 + 2*elementsCountAroundLeftAtriumLPV
    return elementsCountAroundLpvOstium, elementsCountAroundRpvOstium


def getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall):
    '''
    Get the number of elements in each section of the left atrium free wall.
    :param elementsCountAroundLeftAtriumFreeWall: Valid range 6-10.
    :return: elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftAtrialAppendageBase,
        elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV
    '''
    assert 6 <= elementsCountAroundLeftAtriumFreeWall <= 10, \
        'getLeftAtriumBaseFreewallElementsCounts: elements count out of range: ' + str(elementsCountAroundLeftAtriumFreeWall)
    elementsCountAroundLeftAtriumAorta = 1
    elementsCountAroundLeftAtriumVP = (elementsCountAroundLeftAtriumFreeWall + 1)//2
    elementsCountAroundLeftAtrialAppendageBase = elementsCountAroundLeftAtriumFreeWall - elementsCountAroundLeftAtriumVP - elementsCountAroundLeftAtriumAorta
    elementsCountAroundLeftAtriumRPV = elementsCountAroundLeftAtriumVP//2
    elementsCountAroundLeftAtriumLPV = elementsCountAroundLeftAtriumVP - elementsCountAroundLeftAtriumRPV
    return elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftAtrialAppendageBase, \
        elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV


def getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall):
    '''
    Get the number of elements in each section of the right atrium free wall.
    :param elementsCountAroundRightAtriumFreeWall: Valid range 6-10.
    :return: elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightAtrialAppendagePlainBase
        elementsCountAroundRightAtrialAppendagePouchBase, elementsCountAroundRightAtriumAorta
    '''
    assert 6 <= elementsCountAroundRightAtriumFreeWall <= 10, \
        'getRightAtriumBaseFreewallElementsCounts: elements count out of range: ' + str(elementsCountAroundRightAtriumFreeWall)
    elementsCountAroundRightAtriumAorta = 1
    elementsCountAroundRightAtriumPosteriorVenous = elementsCountAroundRightAtriumFreeWall//4
    elementsCountAroundRightAtrialAppendagePlainBase = 2
    elementsCountAroundRightAtrialAppendagePouchBase = (elementsCountAroundRightAtriumFreeWall - elementsCountAroundRightAtrialAppendagePlainBase
        - elementsCountAroundRightAtriumPosteriorVenous - elementsCountAroundRightAtriumAorta)
    return elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightAtrialAppendagePlainBase, \
        elementsCountAroundRightAtrialAppendagePouchBase, elementsCountAroundRightAtriumAorta


def getOverAtriaElementsCounts(elementsCountOverAtria):
    '''
    Get the number of elements in each section over the atria, at the outer
    interatrial groove.
    :param elementsCountOverAtria: Valid values 6 or 8.
    :return: elementsCountOverAtriaCoronarySinus,
        elementsCountOverLeftAtriumNonVenousAnterior, elementsCountOverLeftAtriumVenous, elementsCountOverLeftAtriumNonVenousPosterior,
        elementsCountOverRightAtriumNonVenousAnterior, elementsCountOverRightAtriumVenous, elementsCountOverRightAtriumNonVenousPosterior
    '''
    assert elementsCountOverAtria in [ 6, 8, 10 ], \
        'getOverAtriaElementsCounts: elements count not 6, 8 or 10: ' + str(elementsCountOverAtria)
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
    laaLeft, laVenousMidpointLeft, raVenousRight, raaPouchRight, elementsCountAroundTrackSurface):
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
    aSeptumCfbSideElementLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laSeptumAnteriorRadians, laCfbLeftRadians)
    aCfbSideDerivativeLength = aSeptumCfbSideElementLength  # GRC was 2.0*aSeptumCfbSideElementLength - atrialSeptumOuterElementLength, but too high
    aRemainingLength = atrialPerimeterLength - aSeptumCfbSideElementLength - elementsCountAroundAtrialSeptum*atrialSeptumOuterElementLength

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
    ltFreeWallElementLength = (aRemainingLength - 0.5*(atrialSeptumOuterElementLength + aCfbSideDerivativeLength))/(elementsCountAroundLeftAtriumFreeWallFixed - 2)
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
        elementLength = 0.5*(aCfbSideDerivativeLength + ltFreeWallElementLength) if (n1 == 0) else ltFreeWallElementLength
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
    xi = 0.85  # GRC fudge factor
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
    er = 0.5*(2.0 - laVenousMidpointLeft)*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    laVenousMidpointX = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[0]
    laVenousMidpointRadians = getEllipseRadiansToX(axOuter, bxOuter, laVenousMidpointX - laCentreX, backRadians)
    #print('laVenousMidpointLeft', laVenousMidpointLeft, 'e', e, 'xi', xi, 'x', laVenousMidpointX, 'radians', laVenousMidpointRadians)
    # note ra points these are computed on the left atrium and mirrored at the end
    er = 0.5*(2.0 - raVenousRight)*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    raVenousRightX = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[0]
    raVenousRightRadians = getEllipseRadiansToX(axOuter, bxOuter, raVenousRightX - laCentreX, backRadians)
    #print('raVenousRight', raVenousRight, 'e', e, 'xi', xi, 'x', raVenousRightX, 'radians', raVenousRightRadians)
    er = 0.5*raaPouchRight*elementsCountAroundTrackSurface
    e = int(er)
    xi = er - e
    raaPouchRightY = interp.interpolateCubicHermite(vx[e], vd1[e], vx[e + 1], vd1[e + 1], xi)[1]
    raaPouchRightRadians = getEllipseRadiansToX(ayOuter, byOuter, raaPouchRightY - laCentreY, sideRadians)
    #print('raaPouchRight', raaPouchRight, 'e', e, 'xi', xi, 'y', raaPouchRightY, 'radians', raaPouchRightRadians)

    # get numbers of elements and lengths of sections of left atrium (outer)
    elementsCountAroundLeftAtriumAorta, elementsCountAroundLeftAtrialAppendageBase, elementsCountAroundLeftAtriumLPV, elementsCountAroundLeftAtriumRPV = \
        getLeftAtriumBaseFreewallElementsCounts(elementsCountAroundLeftAtriumFreeWall)
    laaLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laCfbLeftRadians, laaEndRadians)
    laaLeftLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laaEndRadians, laVenousMidpointRadians)
    laVenousRightLength = aRemainingLength - laaLeftLength - laaLength

    # get element lengths/derivatives at edges of each section and transition element sizes between
    laaEndDerivative = laVenousMidpointDerivative = laaLeftLength/elementsCountAroundLeftAtriumLPV
    laaElementLengths = interp.sampleCubicElementLengths(laaLength, elementsCountAroundLeftAtrialAppendageBase, startDerivative = aCfbSideDerivativeLength, endDerivative = laaEndDerivative)
    lvlElementLengths = interp.sampleCubicElementLengths(laaLeftLength, elementsCountAroundLeftAtriumLPV, startDerivative = laaEndDerivative, endDerivative = laVenousMidpointDerivative)
    lvrElementLengths = interp.sampleCubicElementLengths(laVenousRightLength, elementsCountAroundLeftAtriumRPV, startDerivative = laVenousMidpointDerivative, endDerivative = atrialSeptumOuterElementLength)

    # get radians of nodes around left atrium, starting at cfb
    elementsCountAroundLeftAtrium = elementsCountAroundLeftAtriumFreeWall + elementsCountAroundAtrialSeptum
    laRadians = []
    radiansAround = laSeptumAnteriorRadians
    for n1 in range(elementsCountAroundLeftAtrium):
        laRadians.append(radiansAround)
        if n1 == 0:
            elementLength = aSeptumCfbSideElementLength
        elif n1 < (elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftAtrialAppendageBase):
            elementLength = laaElementLengths[n1 - elementsCountAroundLeftAtriumAorta]
        elif n1 < (elementsCountAroundLeftAtriumAorta + elementsCountAroundLeftAtrialAppendageBase + elementsCountAroundLeftAtriumLPV):
            elementLength = lvlElementLengths[n1 - elementsCountAroundLeftAtriumAorta - elementsCountAroundLeftAtrialAppendageBase]
        elif n1 == (elementsCountAroundLeftAtriumFreeWall - 1):
            radiansAround = laSeptumPosteriorRadians + 2.0*math.pi
            continue
        elif n1 < elementsCountAroundLeftAtriumFreeWall:
            elementLength = lvrElementLengths[n1 - elementsCountAroundLeftAtriumAorta - elementsCountAroundLeftAtrialAppendageBase - elementsCountAroundLeftAtriumLPV]
        else:
            radiansAround = updateEllipseAngleByArcLength(aBaseInnerMajorMag, aBaseInnerMinorMag, radiansAround, atrialSeptumInnerElementLength)
            continue
        radiansAround = updateEllipseAngleByArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, radiansAround, elementLength)

    #print('la radiansAround start', laSeptumAnteriorRadians, 'end', radiansAround - 2.0*math.pi)
    #print('laRadians', laRadians)

    # get numbers of elements and lengths of sections of right atrium (outer)
    elementsCountAroundRightAtriumPosteriorVenous, elementsCountAroundRightAtrialAppendagePlainBase, \
        elementsCountAroundRightAtrialAppendagePouchBase, elementsCountAroundRightAtriumAorta \
        = getRightAtriumBaseFreewallElementsCounts(elementsCountAroundRightAtriumFreeWall)
    elementsCountAroundRightAtrialAppendageBase = elementsCountAroundRightAtrialAppendagePlainBase + elementsCountAroundRightAtrialAppendagePouchBase
    raaPouchLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, laCfbLeftRadians, raaPouchRightRadians)
    raaPlainLength = getEllipseArcLength(aBaseOuterMajorMag, aBaseOuterMinorMag, raaPouchRightRadians, raVenousRightRadians)
    raPosteriorVenousLength = aRemainingLength - raaPlainLength - raaPouchLength

    # get element lengths/derivatives at edges of each section and transition element sizes between
    #raPosteriorVenousLimitDerivative = (raaLength - 0.5*aCfbSideDerivativeLength)/(elementsCountAroundRightAtrialAppendageBase - 0.5)
    raPosteriorVenousLimitDerivative = raaPlainDerivative = raaPlainLength/elementsCountAroundRightAtrialAppendagePlainBase
    raaPouchElementLengths = interp.sampleCubicElementLengths(raaPouchLength, elementsCountAroundRightAtrialAppendagePouchBase, startDerivative = aCfbSideDerivativeLength, endDerivative = raaPlainDerivative)
    raaPlainElementLengths = interp.sampleCubicElementLengths(raaPlainLength, elementsCountAroundRightAtrialAppendagePlainBase, startDerivative = raaPlainDerivative, endDerivative = raPosteriorVenousLimitDerivative)
    ravElementLengths = interp.sampleCubicElementLengths(raPosteriorVenousLength, elementsCountAroundRightAtriumPosteriorVenous, startDerivative = raPosteriorVenousLimitDerivative, endDerivative = atrialSeptumOuterElementLength)

    # get radians of nodes around right atrium (computed on left and mirrored at the end), starting at cfb
    elementsCountAroundRightAtrium = elementsCountAroundRightAtriumFreeWall + elementsCountAroundAtrialSeptum
    raRadians = []
    radiansAround = laSeptumAnteriorRadians
    for n1 in range(elementsCountAroundRightAtrium):
        raRadians.append(radiansAround)
        if n1 == 0:
            elementLength = aSeptumCfbSideElementLength
        elif n1 < (elementsCountAroundRightAtriumAorta + elementsCountAroundRightAtrialAppendagePouchBase):
            elementLength = raaPouchElementLengths[n1 - elementsCountAroundRightAtriumAorta]
        elif n1 < (elementsCountAroundRightAtriumAorta + elementsCountAroundRightAtrialAppendageBase):
            elementLength = raaPlainElementLengths[n1 - elementsCountAroundRightAtriumAorta - elementsCountAroundRightAtrialAppendagePouchBase]
        elif n1 == (elementsCountAroundRightAtriumFreeWall - 1):
            radiansAround = laSeptumPosteriorRadians + 2.0*math.pi
            continue
        elif n1 < elementsCountAroundRightAtriumFreeWall:
            elementLength = ravElementLengths[n1 - elementsCountAroundRightAtriumAorta - elementsCountAroundRightAtrialAppendageBase]
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
    ad1 = [ -iaGrooveDerivative, 0.0, math.tan(math.pi / 3) * iaGrooveDerivative ]
    dx = vx[elementsCountAlongTrackSurface]
    dd1 = [ -d for d in vd2[elementsCountAlongTrackSurface]]
    # fudge factor
    px, pd1 = interp.sampleCubicHermiteCurves([ ax, dx ], [ ad1, dd1 ], elementsCountOut = 2, lengthFractionStart = 0.6, arcLengthDerivatives = True)[0:2]
    nx = [ ax, [ px[1][0], px[1][1], aOuterHeight ] ]
    nd1 = interp.smoothCubicHermiteDerivativesLine(nx, [ ad1, [ pd1[1][0], pd1[1][1], 0.0 ] ], fixStartDerivative = True, fixEndDirection = True)
    cx = nx[1]
    cd1 = nd1[1]
    # recalculate interatrial groove dervivative ad1
    #ad1 = interp.interpolateLagrangeHermiteDerivative(ax, cx, cd1, 0.0)
    nx = [ax, cx, dx]
    nd1 = interp.smoothCubicHermiteDerivativesLine(nx, [ ad1, cd1, dd1 ], fixAllDirections=True,
                                                   magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    rx, rd1 = interp.sampleCubicHermiteCurves([ax, cx, dx], nd1,
                                              elementsCountOut=elementsCountAlongTrackSurface,
                                              arcLengthDerivatives=True)[0:2]

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
    # add last point at central outer point of atrium, d1 all zero, d2 rotating around
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
        sd2 = interp.smoothCubicHermiteDerivativesLine(
            sx, sd2, fixStartDirection=True if n1 in (1, elementsCountAcrossTrackSurface - 1) else False,
            fixEndDirection=True)
        for n2 in range(elementsCountAcrossTrackSurface + 1):
            n = n2*(elementsCountAcrossTrackSurface + 1) + n1
            nd2[n] = sd2[n2]

    return TrackSurface(elementsCountAcrossTrackSurface, elementsCountAlongTrackSurface, nx, nd1, nd2)


def getAtrialAppendageWedgePoints(basex, based1, based2, based3, angle1radians, angle2radians, angle3Radians, baseLength,
        elementsCountAroundAppendage, elementsCountRadial, arcLength, arcRadius, wallThickness, wedgeAngleRadians):
    '''
    Get points on wedge midline at end of atrial appendage.
    :param basex, based1, based2, based3: Base point on atrium to project from.
    :param angle1radians, angle2radians: Rotation of direction toward d1, d2.
    :param baseLength: Distance to project out from base.
    :param elementsCountAcrossWedge: Elements around arc, one less than numbers of points out.
    :param elementsCountAlongAppendage: Number of radial elements along appendage, used to
    determine end derivatives.
    :param arcLength: Length of outer arc around atrial appendage wedge.
    :param arcRadius: Radius of outer arc. Must be < pi
    :param wallThickness: Atrial appendage wall thickness
    :param wedgeAngleRadians: Angle from base to top at end of wedge.
    :return: aawx[n3][n1], aawd1[n3][n1], aawd2[n3][n1], aawd3[n3][n1], 
        elementsCountAcrossWedge, wedgePointsMap, wedgeDerivativesMap.
    where:
        n3 is range 2 inner, outer; n1 is 0-elementsCountAcrossWedge
        wedgePointsMap maps points count around appendage --> node index across wedge
        derivativesMap is for passing to createAnnulusMesh3d.
    n is number of points around arc.
    '''
    elementsCountAcrossWedge = (elementsCountAroundAppendage - 4)//2
    # wedge centre:
    wcx, wd1, wd2, wd3 = getCircleProjectionAxes(basex, based1, based2, based3, baseLength, angle1radians, angle2radians, angle3Radians)
    # arc centre:
    acx = [ (wcx[c] - arcRadius*wd3[c]) for c in range(3) ]
    wedgeLength = wallThickness/math.tan(0.5*wedgeAngleRadians)
    cosHalfWedgeAngleRadians = math.cos(0.5*wedgeAngleRadians)
    sinHalfWedgeAngleRadians = math.sin(0.5*wedgeAngleRadians)
    aawx  = [ [], [] ]
    aawd1 = [ [], [] ]
    aawd2 = [ [], [] ]
    aawd3 = [ [], [] ]
    arcRadians = arcLength/arcRadius
    elementLengthOuterRadial = baseLength/elementsCountRadial
    for n2 in range(2):
        if n2 == 0:
            radius = arcRadius - wedgeLength
            derivativeMag = elementLengthOuterRadial - wedgeLength
        else:
            radius = arcRadius
            derivativeMag = elementLengthOuterRadial
        elementLengthArc = (arcLength/elementsCountAcrossWedge)*(radius/arcRadius)
        for n1 in range(0, elementsCountAcrossWedge + 1):
            angleRadians = arcRadians*(n1/elementsCountAcrossWedge - 0.5)
            cosAngleRadians = math.cos(angleRadians)
            sinAngleRadians = math.sin(angleRadians)
            x  = [ (acx[c] + radius*(cosAngleRadians*wd3[c] + sinAngleRadians*wd1[c])) for c in range(3) ]
            d1 = [ elementLengthArc*(cosAngleRadians*wd1[c] - sinAngleRadians*wd3[c]) for c in range(3) ]
            d2 = [ derivativeMag*( cosHalfWedgeAngleRadians*(-sinAngleRadians*wd1[c] - cosAngleRadians*wd3[c]) + sinHalfWedgeAngleRadians*wd2[c]) for c in range(3) ]
            d3 = [ derivativeMag*(-cosHalfWedgeAngleRadians*(-sinAngleRadians*wd1[c] - cosAngleRadians*wd3[c]) + sinHalfWedgeAngleRadians*wd2[c]) for c in range(3) ]
            #d3 = [ derivativeMag*(cosAngleRadians*wd3[c] + sinAngleRadians*wd1[c]) for c in range(3) ]  # flat d3, if using full WedgeAngleRadians
            aawx [n2].append(x )
            aawd1[n2].append(d1)
            aawd2[n2].append(d2)
            aawd3[n2].append(d3)

    wedgePointsMap = [ 0 ]*3 + list(range(1, elementsCountAcrossWedge)) + [ elementsCountAcrossWedge ]*3 + list(range(elementsCountAcrossWedge - 1, 0, -1))
    wedgeDerivativesMap = ( [ ( ( -1,  0,  0 ), (  0, -1,  0 ), None, (  0,  0,  0 ) ),
                              ( (  0,  0,  0 ), ( +1,  0,  0 ), None, (  0,  0,  0 ) ),  # inside corner, collapsed
                              ( (  0,  0,  0 ), (  0,  0, +1 ), None, None         ) ]
                          + [ ( None          , (  0,  0, +1 ), None ) ]*(elementsCountAcrossWedge - 1)  # bottom
                          + [ ( None          , (  0,  0, +1 ), None, (  0,  0,  0 ) ),
                              ( (  0,  0,  0 ), ( -1,  0,  0 ), None, (  0,  0,  0 ) ),  # inside corner, collapsed
                              ( (  0,  0,  0 ), (  0, -1,  0 ), None, ( -1,  0,  0 ) ) ]
                          + [ ( ( -1,  0,  0 ), (  0, -1,  0 ), None ) ]*(elementsCountAcrossWedge - 1) )  # top
    if (elementsCountAroundAppendage%2) == 1:
        wedgePointsMap.insert(1, 0)
        wedgeDerivativesMap.insert(1, wedgeDerivativesMap[1])
    #print('wedgePointsMap',wedgePointsMap)
    #print('wedgeDerivativesMap',wedgeDerivativesMap)
    return aawx, aawd1, aawd2, aawd3, elementsCountAcrossWedge, wedgePointsMap, [ wedgeDerivativesMap, wedgeDerivativesMap ]
