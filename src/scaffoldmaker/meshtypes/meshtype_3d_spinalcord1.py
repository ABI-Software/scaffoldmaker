"""
Generates a spinal cord scaffold using the solid cylinder elements,
"""

from __future__ import division
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.annotation.cns_terms import get_spinalcord_term
from opencmiss.zinc.node import Node


class MeshType_3d_spinalcord1(Scaffold_base):
    """
Generates a spinal cord mesh using a solid cylinder
    """
    centralPathDefaultScaffoldPackages = {
        'Spinal cord 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 5.0,
                'Number of elements': 5
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[-3.1101e+01, 5.1038e+01, -3.4868e+01], [-1.4473e+01, -8.2693e+00, -1.2407e+01],
                     [-4.2508e-01, 1.7548e+00, -6.7373e-01], [-2.3520e-01, -2.3481e-01, 2.4575e-01],
                     [1.5634e+00, -2.5602e-01, -1.6532e+00], [3.0517e-01, 1.3668e+00, 2.5293e+00]],
                    [[-3.8901e+01, 4.4771e+01, -4.7690e+01], [5.4571e+00, -3.1608e+00, -1.6864e+01],
                     [-5.0153e-01, 1.4671e+00, -4.3727e-01], [8.2260e-02, -3.4051e-01, 2.2716e-01],
                     [1.3947e+00, 5.7896e-01, 3.4282e-01], [-6.4267e-01, 3.0315e-01, 1.4628e+00]],
                    [[-2.6488e+01, 4.4886e+01, -6.0429e+01], [1.8955e+01, 3.4567e+00, -4.5915e+00],
                     [-2.4840e-01, 1.0697e+00, -2.2012e-01], [2.6155e-01, 3.6440e-02, 1.4932e-01],
                     [2.4185e-01, 3.0958e-01, 1.2315e+00], [-4.6277e-02, -2.5595e-01, 4.1616e-01]],
                    [[-4.1370e+00, 4.5254e+01, -7.0292e+01], [1.5429e+01, -2.7169e+00, -2.5934e+01],
                     [2.4610e-02, 1.6965e+00, -1.6309e-01], [8.9196e-02, -2.2273e-02, 9.7179e-02],
                     [1.7014e+00, 7.1902e-02, 1.0046e+00], [5.3281e-01, -9.8096e-02, -6.1819e-01]],
                    [[-6.6576e+00, 3.9713e+01, -1.0278e+02], [-5.4766e+00, -1.5383e+00, -4.2455e+01],
                     [-1.3103e-01, 8.0965e-01, -1.2434e-02], [-6.5550e-02, -7.5651e-01, 9.8784e-02],
                     [9.9977e-01, 1.5973e-01, -1.3475e-01], [-6.7179e-01, 2.5517e-02, -6.8272e-01]],
                    [[-1.7441e+01, 3.6668e+01, -1.5334e+02], [-1.3066e+01, -1.8613e+00, -4.2825e+01],
                     [-5.8584e-02, 2.5282e-01, 6.8850e-03], [2.1046e-01, -3.5719e-01, -6.0145e-02],
                     [3.7369e-01, 8.9800e-02, -1.1791e-01], [-5.8040e-01, -1.6538e-01, 7.1640e-01]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Spinal Cord 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Spinal cord 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements along': 12,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1
        }
        if parameterSetName == 'Default':
            parameterSetName = 'Rat 1'
        options['Base parameter set'] = parameterSetName
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements across shell',
            'Number of elements along',
            'Lower half',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
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
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1
        if options['Number of elements along'] < 1:
            options['Number of elements along'] = 1
        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isRat = 'Rat' in parameterSetName

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        mesh = fm.findMeshByDimension(3)

        cervicalGroup = AnnotationGroup(region, get_spinalcord_term("Cervical spinal cord"))
        thoracicGroup = AnnotationGroup(region, get_spinalcord_term("Thoracic spinal cord"))
        lumbarGroup = AnnotationGroup(region, get_spinalcord_term("Lumbar spinal cord"))
        sacralGroup = AnnotationGroup(region, get_spinalcord_term("Sacral spinal cord"))
        caudalGroup = AnnotationGroup(region, get_spinalcord_term("caudal segment of spinal cord"))

        # Base cylinder mesh
        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell, [0.0, 0.0, 0.0],
                            cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # Making the groups for different parts of the spinal cord
        # cervical
        cervicalMeshGroup = cervicalGroup.getMeshGroup(mesh)
        cervicalElements = elementsCountAlong // 6*(elementsCountAcrossMajor * elementsCountAcrossMinor - 4)
        for e in range(1, cervicalElements + 1):
            element = mesh.findElementByIdentifier(e)
            cervicalMeshGroup.addElement(element)

        # thoracic
        thoracicMeshGroup = thoracicGroup.getMeshGroup(mesh)
        thoracicElements = elementsCountAlong // 2*(elementsCountAcrossMajor * elementsCountAcrossMinor - 4)
        for e in range(cervicalElements + 1, thoracicElements + 1):
            element = mesh.findElementByIdentifier(e)
            thoracicMeshGroup.addElement(element)

        # lumbar
        lumbarMeshGroup = lumbarGroup.getMeshGroup(mesh)
        lumbarElements = elementsCountAlong // 4*(elementsCountAcrossMajor * elementsCountAcrossMinor - 4)
        for e in range(thoracicElements + 1, thoracicElements + lumbarElements + 1):
            element = mesh.findElementByIdentifier(e)
            lumbarMeshGroup.addElement(element)

        # sacral
        sacralMeshGroup = sacralGroup.getMeshGroup(mesh)
        sacralElements = elementsCountAlong // 6*(elementsCountAcrossMajor * elementsCountAcrossMinor - 4)
        for e in range(thoracicElements + lumbarElements + 1, thoracicElements + lumbarElements + sacralElements + 1):
            element = mesh.findElementByIdentifier(e)
            sacralMeshGroup.addElement(element)

        # caudal
        caudalMeshGroup = caudalGroup.getMeshGroup(mesh)
        caudalElements = elementsCountAlong*(elementsCountAcrossMajor * elementsCountAcrossMinor - 4)
        for e in range(caudalElements - (elementsCountAcrossMajor * elementsCountAcrossMinor - 4) + 1, caudalElements + 1):
            element = mesh.findElementByIdentifier(e)
            caudalMeshGroup.addElement(element)

        annotationGroup = [cervicalGroup, thoracicGroup, lumbarGroup, sacralGroup, caudalGroup]
        return annotationGroup

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)
