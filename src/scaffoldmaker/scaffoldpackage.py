"""
Class packaging a scaffold type, options and modifications.
Supports serialisation to/from JSON. Can be used as a scaffold option.
"""

import copy
import math
from opencmiss.utils.zinc.field import createFieldEulerAnglesRotationMatrix
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.utils.maths import vectorops
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

class ScaffoldPackage:
    '''
    Class packaging a scaffold type, options and modifications.
    '''

    def __init__(self, scaffoldType, dct={}, defaultParameterSetName='Default'):
        '''
        :param scaffoldType: A scaffold type derived from Scaffold_base.
        :param dct: Dictionary containing other scaffold settings. Key names and meanings:
            scaffoldSettings: The options dict for the scaffold, or None to generate defaults.
            meshEdits: A Zinc model file as a string e.g. containing edited node parameters, or None.
        :param defaultParameterSetName: Parameter set name from scaffoldType to get defaults from.
        '''
        #print('ScaffoldPackage.__init__',dct)
        assert issubclass(scaffoldType, Scaffold_base), 'ScaffoldPackage:  Invalid scaffold type'
        self._scaffoldType = scaffoldType
        # merge with defaults to ensure new options for scaffold type are present
        self._scaffoldSettings = scaffoldType.getDefaultOptions(defaultParameterSetName)
        scaffoldSettings = dct.get('scaffoldSettings')
        if scaffoldSettings:
            # remove obsolete options? If so, deepcopy first?
            self._scaffoldSettings.update(scaffoldSettings)
        # note rotation is stored in degrees
        self._rotation = dct.get('rotation')
        if not self._rotation:
            self._rotation = [ 0.0, 0.0, 0.0 ]
        self._scale = dct.get('scale')
        if not self._scale:
            self._scale = [ 1.0, 1.0, 1.0 ]
        self._translation = dct.get('translation')
        if not self._translation:
            self._translation = [ 0.0, 0.0, 0.0 ]
        self._meshEdits = copy.deepcopy(dct.get('meshEdits'))

    def deepcopy(self, other):
        '''
        Deep copy contents from another object.
        '''
        self._scaffoldType = other._scaffoldType
        self._scaffoldSettings = copy.deepcopy(other._scaffoldSettings)
        self._meshEdits = copy.deepcopy(other._meshEdits)

    def __eq__(self, other):
        '''
        Need equality operator to determine if custom options are in use.
        '''
        if isinstance(other, ScaffoldPackage):
            return (self._scaffoldType == other._scaffoldType) \
                and (self._scaffoldSettings == other._scaffoldSettings) \
                and (self._meshEdits == other._meshEdits)
        return NotImplemented

    def toDict(self):
        '''
        Encodes object into a dictionary for JSON serialisation.
        Key names are described in __init__().
        Scaffold type is encoded separately.
        :return: Dictionary containing object encoding.
        '''
        dct = {
            'rotation' : self._rotation,
            'scaffoldTypeName' : self._scaffoldType.getName(),
            'scaffoldSettings' : self._scaffoldSettings,
            'scale' : self._scale,
            'translation' : self._translation
            }
        if self._meshEdits:
            dct['meshEdits'] = self._meshEdits
        return dct

    def getMeshEdits(self):
        return self._meshEdits

    def setMeshEdits(self, meshEdits):
        self._meshEdits = meshEdits

    def getScaffoldSettings(self):
        return self._scaffoldSettings

    def getScaffoldType(self):
        return self._scaffoldType

    def getRotation(self):
        return self._rotation

    def setRotation(self, rotation):
        '''
        :param rotation: list of 3 rotation angles about z, y', x'', in degrees.
        :return: True if value changed, otherwise False.
        '''
        assert len(rotation) == 3
        if self._rotation != rotation:
            self._rotation = rotation
            return True
        return False

    def getScale(self):
        return self._scale

    def setScale(self, scale):
        '''
        :param scale: list of 3 scales for x, y, z.
        :return: True if value changed, otherwise False.
        '''
        assert len(scale) == 3
        if self._scale != scale:
            self._scale = scale
            return True
        return False

    def getTranslation(self):
        return self._translation

    def setTranslation(self, translation):
        '''
        :param translation: list of 3 translation values for x, y, z.
        :return: True if value changed, otherwise False.
        '''
        assert len(translation) == 3
        if self._translation != translation:
            self._translation = translation
            return True
        return False

    def getTransformationMatrix(self):
        '''
        :return: 4x4 row-major transformation matrix with first index down rows, second across columns,
        suitable for multiplication p' = Mp where p = [ x, y, Z, h ], or None if identity.
        '''
        # apply transformation in order: scale then rotation then translation
        if not all((v == 0.0) for v in self._rotation):
            rotationMatrix = vectorops.eulerToRotationMatrix3([ deg*math.pi/180.0 for deg in self._rotation ])
            return [
                [ rotationMatrix[0][0]*self._scale[0], rotationMatrix[0][1]*self._scale[1], rotationMatrix[0][2]*self._scale[2], self._translation[0] ],
                [ rotationMatrix[1][0]*self._scale[0], rotationMatrix[1][1]*self._scale[1], rotationMatrix[1][2]*self._scale[2], self._translation[1] ],
                [ rotationMatrix[2][0]*self._scale[0], rotationMatrix[2][1]*self._scale[1], rotationMatrix[2][2]*self._scale[2], self._translation[2] ],
                [ 0.0, 0.0, 0.0, 1.0 ] ]
        if not (all((v == 1.0) for v in self._scale) and all((v == 0.0) for v in self._translation)):
            return [
                [ self._scale[0], 0.0, 0.0, self._translation[0] ],
                [ 0.0, self._scale[1], 0.0, self._translation[1] ],
                [ 0.0, 0.0, self._scale[2], self._translation[2] ],
                [ 0.0, 0.0, 0.0, 1.0 ] ]
        return None

    def applyTransformation(self, region):
        '''
        If rotation, scale or transformation are set, transform node coordinates.
        '''
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName('coordinates').castFiniteElement()
        if not coordinates.isValid():
            print('Warning: ScaffoldPackage.applyTransformation: Missing coordinates field')
            return
        with ChangeManager(fieldmodule):
            componentsCount = coordinates.getNumberOfComponents()
            if componentsCount < 3:
                # pad with zeros
                coordinates = fieldmodule.createFieldConcatenate([ coordinates ] + [ fieldmodule.createFieldConstant([ 0.0 ]*(3 - componentsCount)) ])
            newCoordinates = coordinates
            # apply scale first so variable scaling in x, y, z doesn't interplay with rotation
            if not all((v == 1.0) for v in self._scale):
                #print("applyTransformation: apply scale", self._scale)
                newCoordinates = newCoordinates*fieldmodule.createFieldConstant(self._scale)
            if not all((v == 0.0) for v in self._rotation):
                #print("applyTransformation: apply rotation", self._rotation)
                newCoordinates = fieldmodule.createFieldMatrixMultiply(3,
                    createFieldEulerAnglesRotationMatrix(fieldmodule, fieldmodule.createFieldConstant([ deg*math.pi/180.0 for deg in self._rotation ])), newCoordinates)
            if not all((v == 0.0) for v in self._translation):
                #print("applyTransformation: apply translation", self._translation)
                newCoordinates = newCoordinates + fieldmodule.createFieldConstant(self._translation)
            # be sure to delete temporary fields and fieldassignment to reduce messages
            if newCoordinates is not coordinates:
                fieldassignment = coordinates.createFieldassignment(newCoordinates)
                fieldassignment.assign()
                del fieldassignment
            del newCoordinates
            del coordinates

    def generate(self, region, applyTransformation=True):
        '''
        :param applyTransformation: If True (default) apply scale, rotation and translation to
        node coordinates. Specify False if client will transform, e.g. with graphics transformations.
        '''
        #print('\nScaffoldPackage.generate: ', self.toDict())
        annotationGroups = self._scaffoldType.generateMesh(region, self._scaffoldSettings)
        if self._meshEdits:
            # apply mesh edits, a Zinc-readable model file containing node edits
            # Note: these are untransformed coordinates
            sir = region.createStreaminformationRegion()
            srm = sir.createStreamresourceMemoryBuffer(self._meshEdits)
            region.read(sir)
        if applyTransformation:
            self.applyTransformation(region)
        return annotationGroups
