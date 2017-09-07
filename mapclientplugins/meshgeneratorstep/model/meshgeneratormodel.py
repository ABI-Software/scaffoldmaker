'''
Created on Aug 29, 2017

@author: Richard Christie
'''

import os, sys
import json

from opencmiss.zinc.context import Context
from opencmiss.zinc.status import OK as ZINC_OK
from opencmiss.zinc.field import Field
from opencmiss.zinc.glyph import Glyph
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_2d_plate1 import MeshType_2d_plate1
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_2d_sphere1 import MeshType_2d_sphere1
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_2d_tube1 import MeshType_2d_tube1
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_box1 import MeshType_3d_box1
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_tube1 import MeshType_3d_tube1

class MeshGeneratorModel(object):
    '''
    Framework for generating meshes of a number of types, with mesh type specific options
    '''

    def __init__(self, location):
        '''
        Constructor
        '''
        self._location = location
        self._context = Context("MeshGenerator")
        tess = self._context.getTessellationmodule().getDefaultTessellation()
        tess.setRefinementFactors(12)
        self._sceneChangeCallback = None
        # set up standard materials and glyphs so we can use them elsewhere
        self._materialmodule = self._context.getMaterialmodule()
        self._materialmodule.defineStandardMaterials()
        glyphmodule = self._context.getGlyphmodule()
        glyphmodule.defineStandardGlyphs()
        self._deleteElementRanges = []
        self._settings = {
            'meshTypeName' : '',
            'meshTypeOptions' : { },
            'identifierOffset' : 0,
            'deleteElementRanges' : '',
            'displayNodeNumbers' : True,
            'displayElementNumbers' : True
        }
        self._discoverAllMeshTypes()
        self._loadSettings()
        self._generateMesh()

    def _discoverAllMeshTypes(self):
        self._meshTypes = [
            MeshType_2d_plate1,
            MeshType_2d_sphere1,
            MeshType_2d_tube1,
            MeshType_3d_box1,
            MeshType_3d_sphereshell1,
            MeshType_3d_tube1
            ]
        self._currentMeshType = MeshType_3d_box1
        self._settings['meshTypeName'] = self._currentMeshType.getName()
        self._settings['meshTypeOptions'] = self._currentMeshType.getDefaultOptions()

    def getAllMeshTypeNames(self):
        meshTypeNames = []
        for meshType in self._meshTypes:
            meshTypeNames.append(meshType.getName())
        return meshTypeNames

    def getMeshTypeName(self):
        return self._settings['meshTypeName']

    def _getMeshTypeByName(self, name):
        for meshType in self._meshTypes:
            if meshType.getName() == name:
                return meshType
        return None

    def setMeshTypeByName(self, name):
        meshType = self._getMeshTypeByName(name)
        if meshType is not None:
            if meshType != self._currentMeshType:
                self._currentMeshType = meshType
                self._settings['meshTypeName'] = self._currentMeshType.getName()
                self._settings['meshTypeOptions'] = self._currentMeshType.getDefaultOptions()
                self._generateMesh()

    def getMeshTypeOrderedOptionNames(self):
        return self._currentMeshType.getOrderedOptionNames()

    def getMeshTypeOption(self, key):
        return self._settings['meshTypeOptions'][key]

    def setMeshTypeOption(self, key, value):
        oldValue = self._settings['meshTypeOptions'][key]
        # print('setMeshTypeOption: key ', key, ' value ', str(value))
        newValue = None
        try:
            if type(oldValue) is bool:
                newValue = bool(value)
            elif type(oldValue) is int:
                newValue = int(value)
            elif type(oldValue) is float:
                newValue = float(value)
            elif type(oldValue) is str:
                newValue = str(value)
            else:
                newValue = value
        except:
            print('setMeshTypeOption: Invalid value')
            return
        self._settings['meshTypeOptions'][key] = newValue
        self._currentMeshType.checkOptions(self._settings['meshTypeOptions'])
        # print('final value = ', self._settings['meshTypeOptions'][key])
        if self._settings['meshTypeOptions'][key] != oldValue:
            self._generateMesh()

    def getDeleteElementsRangesText(self):
        return self._settings['deleteElementRanges']

    def _parseDeleteElementsRangesText(self, elementRangesTextIn):
        '''
        :return: True if ranges changed, otherwise False
        '''
        elementRanges = []
        for elementRangeText in elementRangesTextIn.split(','):
            try:
                elementRangeEnds = elementRangeText.split('-')
                elementRangeStart = int(elementRangeEnds[0])
                if len(elementRangeEnds) > 1:
                    elementRangeStop = int(elementRangeEnds[1])
                else:
                    elementRangeStop = elementRangeStart
                if elementRangeStop >= elementRangeStart:
                    elementRanges.append([elementRangeStart, elementRangeStop])
                else:
                    elementRanges.append([elementRangeStop, elementRangeStart])
            except:
                pass
        elementRangesText = ''
        first = True
        for elementRange in elementRanges:
            if first:
                first = False
            else:
                elementRangesText += ','
            elementRangesText += str(elementRange[0])
            if elementRange[1] != elementRange[0]:
                elementRangesText += '-' + str(elementRange[1])
        changed = self._deleteElementRanges != elementRanges
        self._deleteElementRanges = elementRanges
        self._settings['deleteElementRanges'] = elementRangesText
        return changed

    def setDeleteElementsRangesText(self, elementRangesTextIn):
        if self._parseDeleteElementsRangesText(elementRangesTextIn):
            self._generateMesh()

    def getContext(self):
        return self._context

    def getRegion(self):
        return self._region

    def registerSceneChangeCallback(self, sceneChangeCallback):
        self._sceneChangeCallback = sceneChangeCallback

    def getScene(self):
        return self._region.getScene()

    def _loadSettings(self):
        try:
            with open(self._location + '-settings.json', 'r') as f:
                self._settings.update(json.loads(f.read()))
            self._currentMeshType = self._getMeshTypeByName(self._settings['meshTypeName'])
            self._parseDeleteElementsRangesText(self._settings['deleteElementRanges'])
        except:
            pass  # no settings saved yet

    def _saveSettings(self):
        with open(self._location + '-settings.json', 'w') as f:
            f.write(json.dumps(self._settings, default=lambda o: o.__dict__, sort_keys=True, indent=4))

    def _generateMesh(self):
        self._region = self._context.createRegion()
        fm = self._region.getFieldmodule()
        fm.beginChange()
        self._currentMeshType.generateMesh(self._region, self._settings['meshTypeOptions'])
        for dimension in range(3,0,-1):
            mesh = fm.findMeshByDimension(dimension)
            if mesh.getSize() > 0:
                break
        if len(self._deleteElementRanges) > 0:
            deleteElementIdentifiers = []
            elementIter = mesh.createElementiterator()
            element = elementIter.next()
            while element.isValid():
                identifier = element.getIdentifier()
                for deleteElementRange in self._deleteElementRanges:
                    if (identifier >= deleteElementRange[0]) and (identifier <= deleteElementRange[1]):
                        deleteElementIdentifiers.append(identifier)
                element = elementIter.next()
            #print('delete elements ', deleteElementIdentifiers)
            for identifier in deleteElementIdentifiers:
                element = mesh.findElementByIdentifier(identifier)
                mesh.destroyElement(element)
            del element
            # destroy all orphaned nodes
            nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            #size1 = nodes.getSize()
            nodes.destroyAllNodes()
            #size2 = nodes.getSize()
            #print('deleted', size1 - size2, 'nodes')
        fm.defineAllFaces()
        #for dimension in range(3,0,-1):
        #    mesh = fm.findMeshByDimension(dimension)
        #    print(dimension, '-D element count: ', mesh.getSize())
        fm.endChange()
        # TODO: offset identifiers of elements, faces, lines, nodes
        # TODO: delete elements and any orphaned faces, lines and nodes
        coordinates = fm.findFieldByName('coordinates')
        cmiss_number = fm.findFieldByName('cmiss_number')
        # make graphics
        scene = self._region.getScene()
        scene.beginChange()
        axes = scene.createGraphicsPoints()
        pointattr = axes.getGraphicspointattributes()
        pointattr.setGlyphShapeType(Glyph.SHAPE_TYPE_AXES_XYZ)
        pointattr.setBaseSize([1.0,1.0,1.0])
        axes.setMaterial(self._materialmodule.findMaterialByName('grey50'))
        lines = scene.createGraphicsLines()
        lines.setCoordinateField(coordinates)
        nodeNumbers = scene.createGraphicsPoints()
        nodeNumbers.setFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeNumbers.setCoordinateField(coordinates)
        pointattr = nodeNumbers.getGraphicspointattributes()
        pointattr.setLabelField(cmiss_number)
        pointattr.setGlyphShapeType(Glyph.SHAPE_TYPE_NONE)
        nodeNumbers.setMaterial(self._materialmodule.findMaterialByName('green'))
        elementNumbers = scene.createGraphicsPoints()
        elementNumbers.setFieldDomainType(Field.DOMAIN_TYPE_MESH_HIGHEST_DIMENSION)
        elementNumbers.setCoordinateField(coordinates)
        pointattr = elementNumbers.getGraphicspointattributes()
        pointattr.setLabelField(cmiss_number)
        pointattr.setGlyphShapeType(Glyph.SHAPE_TYPE_NONE)
        elementNumbers.setMaterial(self._materialmodule.findMaterialByName('cyan'))
        scene.endChange()
        if self._sceneChangeCallback is not None:
            self._sceneChangeCallback()

    def getOutputModelFilename(self):
        return self._location + '.ex2'

    def _writeModel(self):
        self._region.writeFile(self.getOutputModelFilename())

    def done(self):
        self._saveSettings()
        self._writeModel()
