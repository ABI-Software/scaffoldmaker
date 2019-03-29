"""
Class packaging a scaffold type, options and modifications.
Supports serialisation to/from JSON. Can be used as a scaffold option.
"""

from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

class ScaffoldPackage:
    '''
    Class packaging a scaffold type, options and modifications.
    '''

    def __init__(self, scaffoldType, dct={}):
        '''
        :param scaffoldType: A scaffold type derived from Scaffold_base.
        :param dct: Dictionary containing other scaffold settings. Key names and meanings:
            scaffoldSettings: The options dict for the scaffold, or None to generate defaults.
            meshEdits: A Zinc model file as a string e.g. containing edited node parameters, or None.
        '''
        #print('ScaffoldPackage.__init__',dct)
        assert issubclass(scaffoldType, Scaffold_base), 'ScaffoldPackage:  Invalid scaffold type'
        self._scaffoldType = scaffoldType
        # merge with defaults to ensure new options for scaffold type are present
        self._scaffoldSettings = scaffoldType.getDefaultOptions()
        scaffoldSettings = dct.get('scaffoldSettings')
        if scaffoldSettings:
            # remove obsolete options?
            self._scaffoldSettings.update(scaffoldSettings)
        self._meshEdits = dct.get('meshEdits')

    def __eq__(self, other):
        '''
        Need equality operator to determine if custom options are in use.
        '''
        if isinstance(other, ScaffoldPackage):
            return (self._scaffoldType == other._scaffoldType) \
                and (self._scaffoldSettings == other._scaffoldSettings) \
                and (self._meshEdits == other._meshEdits)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, ScaffoldPackage):
            return (self._scaffoldType != other._scaffoldType) \
                or (self._scaffoldSettings != other._scaffoldSettings) \
                or (self._meshEdits != other._meshEdits)
        return NotImplemented

    def toDict(self):
        '''
        Encodes object into a dictionary for JSON serialisation.
        Key names are described in __init__().
        Scaffold type is encoded separately.
        :return: Dictionary containing object encoding.
        '''
        dct = {
            'scaffoldTypeName' : self._scaffoldType.getName(),
            'scaffoldSettings' : self._scaffoldSettings
            }
        if self._meshEdits:
            dct['meshEdits'] = self._meshEdits
        return dct

    def setMeshEdits(self, meshEdits):
        self._meshEdits = meshEdits

    def getScaffoldSettings(self):
        return self._scaffoldSettings

    def setScaffoldSettings(self, scaffoldSettings):
        self._scaffoldSettings = scaffoldSettings

    def getScaffoldType(self):
        return self._scaffoldType

    def setScaffoldType(self, scaffoldType):
        if scaffoldType != self._scaffoldType:
            self.__init__(scaffoldType)

    def generate(self, region):
        return self._scaffoldType.generateMesh(region, self._scaffoldSettings)
