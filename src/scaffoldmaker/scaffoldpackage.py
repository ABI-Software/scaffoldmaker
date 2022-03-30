"""
Class packaging a scaffold type, options and modifications.
Supports serialisation to/from JSON. Can be used as a scaffold option.
"""

import copy
import math

from opencmiss.maths.vectorops import euler_to_rotation_matrix
from opencmiss.utils.zinc.field import createFieldEulerAnglesRotationMatrix
from opencmiss.utils.zinc.finiteelement import get_maximum_node_identifier
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findAnnotationGroupByName
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
        rotation =  dct.get('rotation')
        self._rotation = copy.deepcopy(rotation) if rotation else [ 0.0, 0.0, 0.0 ]
        scale = dct.get('scale')
        self._scale = copy.deepcopy(scale) if scale else [ 1.0, 1.0, 1.0 ]
        translation = dct.get('translation')
        self._translation = copy.deepcopy(translation) if translation else [ 0.0, 0.0, 0.0 ]
        meshEdits = dct.get('meshEdits')
        if meshEdits:
            # ensure stored as bytes to match what Zinc creates
            if isinstance(meshEdits, str):
                meshEdits = bytes(meshEdits, 'utf-8')
            else:
                meshEdits = copy.deepcopy(meshEdits)
        self._meshEdits = meshEdits
        self._isGenerated = False  # set to True when generate() is called
        # annotation groups automatically created by scaffold script = set in generate()
        self._autoAnnotationGroups = []
        # read user AnnotationGroups list in dict form:
        userAnnotationGroupsDict = dct.get('userAnnotationGroups')
        # serialised form of user annotation groups, read from serialisation before generate(), updated before writing
        self._userAnnotationGroupsDict = copy.deepcopy(userAnnotationGroupsDict) if userAnnotationGroupsDict else []
        # can only have the actual user annotation groups once generate() is called
        self._userAnnotationGroups = []
        # region is set in generate(); can only instantiate user AnnotationGroups then
        self._region = None
        self._nextNodeIdentifier = 1

    def __eq__(self, other):
        '''
        Need equality operator to determine if custom options are in use.
        '''
        if isinstance(other, ScaffoldPackage):
            return (self._scaffoldType == other._scaffoldType) \
                and (self._scaffoldSettings == other._scaffoldSettings) \
                and (self._rotation == other._rotation) \
                and (self._scale == other._scale) \
                and (self._translation == other._translation) \
                and (self._meshEdits == other._meshEdits) \
                and (self._userAnnotationGroupsDict == other._userAnnotationGroupsDict)
        return NotImplemented

    def __deepcopy__(self, memo):
        '''
        Deep copies object in deserialised, pre-generated form.
        '''
        return ScaffoldPackage(self._scaffoldType, self.toDict())

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
        self.updateUserAnnotationGroups()
        if self._userAnnotationGroupsDict:
            dct['userAnnotationGroups'] = self._userAnnotationGroupsDict
        return dct

    def updateUserAnnotationGroups(self):
        '''
        Ensure user annotation groups are present in serialised form (dict).
        Only done if scaffold has been generated.
        :param force: If not True,
        '''
        if self._isGenerated:
            self._userAnnotationGroupsDict = [ annotationGroup.toDict() for annotationGroup in self._userAnnotationGroups ]

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
        suitable for multiplication p' = Mp where p = [ x, y, z, h ], or None if identity.
        '''
        # apply transformation in order: scale then rotation then translation
        if not all((v == 0.0) for v in self._rotation):
            rotationMatrix = euler_to_rotation_matrix([ deg*math.pi/180.0 for deg in self._rotation ])
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

    def applyTransformation(self):
        '''
        If rotation, scale or transformation are set, transform node coordinates.
        Only call after generate().
        :return: True if a non-identity transformation has been applied, False if not.
        '''
        assert self._region
        fieldmodule = self._region.getFieldmodule()
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
            doApply = newCoordinates is not coordinates
            if doApply:
                fieldassignment = coordinates.createFieldassignment(newCoordinates)
                fieldassignment.assign()
                del fieldassignment
            del newCoordinates
            del coordinates
        return doApply

    def generate(self, region, applyTransformation=True):
        '''
        Generate the finite element scaffold and define annotation groups.
        :param applyTransformation: If True (default) apply scale, rotation and translation to
        node coordinates. Specify False if client will transform, e.g. with graphics transformations.
        '''
        self._region = region
        with ChangeManager(region.getFieldmodule()):
            self._autoAnnotationGroups = self._scaffoldType.generateMesh(region, self._scaffoldSettings)
            # need next node identifier for creating user-defined marker points
            nodes = region.getFieldmodule().findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            self._nextNodeIdentifier = get_maximum_node_identifier(nodes) + 1
            if self._meshEdits:
                # apply mesh edits, a Zinc-readable model file containing node edits
                # Note: these are untransformed coordinates
                sir = region.createStreaminformationRegion()
                srm = sir.createStreamresourceMemoryBuffer(self._meshEdits)
                region.read(sir)
            # define user AnnotationGroups from serialised Dict
            self._userAnnotationGroups = [ AnnotationGroup.fromDict(dct, self._region) for dct in self._userAnnotationGroupsDict ]
            self._isGenerated = True
            if applyTransformation:
                self.applyTransformation()

    def getNextNodeIdentifier(self):
        """
        :return: First node identifier to try using after generating. Used for user-defined marker points.
        """
        return self._nextNodeIdentifier

    def getAnnotationGroups(self):
        '''
        Empty until after call to generate().
        :return: Alphabetically sorted list of annotation groups.
        '''
        return sorted(self._autoAnnotationGroups + self._userAnnotationGroups, key=AnnotationGroup.getName)

    def findAnnotationGroupByName(self, name):
        '''
        Invalid until after call to generate().
        :return: Annotation group with the given name or None.
        '''
        return findAnnotationGroupByName(self._autoAnnotationGroups + self._userAnnotationGroups, name)

    def createUserAnnotationGroup(self, term=None):
        '''
        Create a new, empty user annotation group.
        Only call after generate().
        :param term: Identifier for anatomical term, a tuple of (name, id), or None
        to generate a unique name 'group#' with ID "None".
        If a term is supplied, both its name and id must be strings, the name
        must be unique i.e. unused in the list of annotation groups.
        The id should be a unique string, or use "None" if unknown.
        e.g. ('heart', 'FMA:7088') or None.
        :return: New AnnotationGroup.
        '''
        assert self._region
        if term is not None:
            assert isinstance(term, tuple) and (len(term) == 2) and all(isinstance(s, str) for s in term),\
                "Invalid annotation term " + str(term)
            assert not self.findAnnotationGroupByName(term[0]), "Annotation term " + str(term) + " name is in use"
            useTerm = term
        else:
            number = 1
            while True:
                name = "group" + str(number)
                if not self.findAnnotationGroupByName(name):
                    break
                number += 1
            useTerm = (name, "None")
        annotationGroup = AnnotationGroup(self._region, useTerm)
        self._userAnnotationGroups.append(annotationGroup)
        return annotationGroup

    def deleteAnnotationGroup(self, annotationGroup):
        '''
        Delete the annotation group. Must be a user annotation group.
        :return: True on success, otherwise False
        '''
        if annotationGroup and self.isUserAnnotationGroup(annotationGroup):
            annotationGroup.clear()
            self._userAnnotationGroups.remove(annotationGroup)
            return True
        return False

    def isUserAnnotationGroup(self, annotationGroup):
        '''
        Invalid until after call to generate().
        :return: True if annotationGroup is user-created and editable.
        '''
        return annotationGroup in self._userAnnotationGroups
