
'''
Class for exporting a Scaffold from Zinc to legacy vtk text format.

@author: Richard Christie
'''

import io

from scaffoldmaker.utils.zinc_utils import getElementNodeIdentifiersCube8Node
from opencmiss.zinc.field import Field


class ExportVtk:
    '''
    Class for exporting a Scaffold from Zinc to legacy vtk text format.
    Limited to writing only 3-D hexahedral elements. Assumes all nodes have field defined.
    '''

    def __init__(self, region, description, annotationGroups = None):
        '''
        :param region: Region containing finite element model to export.
        :param description: Single line text description up to 256 characters.
        :param annotationGroups: Optional list of AnnotationGroup for model.
        '''
        self._region = region
        self._fieldmodule = self._region.getFieldmodule()
        self._coordinates = self._fieldmodule.findFieldByName('coordinates')
        self._description = description
        self._annotationGroups = annotationGroups if annotationGroups else []


    def _write(self, outstream):
        assert isinstance(outstream, io.TextIOBase), 'ExportVtk.write:  Invalid outstream argument'
        outstream.write('# vtk DataFile Version 2.0\n')
        outstream.write(self._description + '\n')
        outstream.write('ASCII\n')
        outstream.write('DATASET UNSTRUCTURED_GRID\n')
        nodeIdentifierToIndex = {}  # map needed since vtk points are zero index based, i.e. have no identifier
        coordinates = self._coordinates
        coordinatesCount = coordinates.getNumberOfComponents()
        cache = self._fieldmodule.createFieldcache()

        nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        pointCount = nodes.getSize()
        outstream.write('POINTS ' + str(pointCount) + ' double\n')
        nodeIter = nodes.createNodeiterator()
        node = nodeIter.next()
        index = 0
        while node.isValid():
            nodeIdentifierToIndex[node.getIdentifier()] = index
            cache.setNode(node)
            result, x = self._coordinates.evaluateReal(cache, coordinatesCount)
            first = True
            for s in x:
                if not first:
                    outstream.write(' ')
                outstream.write(str(s))
                first = False
            outstream.write('\n')
            index += 1
            node = nodeIter.next()

        # following assumes all hex elements
        mesh = self._fieldmodule.findMeshByDimension(3)
        cellCount = mesh.getSize()
        cellListSize = 9*cellCount
        outstream.write('CELLS ' + str(cellCount) + ' ' + str(cellListSize) + '\n')
        elementIter = mesh.createElementiterator()
        element = elementIter.next()
        # use vtk winding
        vtkHex8Indexing = [ 0, 1, 3, 2, 4, 5, 7, 6 ]
        while element.isValid():
            eft = element.getElementfieldtemplate(coordinates, -1)  # assumes all components same
            nodeIdentifiersHex8 = getElementNodeIdentifiersCube8Node(element, eft)
            outstream.write('8')
            for localIndex in vtkHex8Indexing:
                index = nodeIdentifierToIndex[nodeIdentifiersHex8[localIndex]]
                outstream.write(' ' + str(index))
            outstream.write('\n')
            element = elementIter.next()
        outstream.write('CELL_TYPES ' + str(cellCount) + '\n')
        for i in range(cellCount - 1):
            outstream.write('12 ')
        outstream.write('12\n')

        if self._annotationGroups:
            outstream.write('CELL_DATA ' + str(cellCount) + '\n')
            for annotationGroup in self._annotationGroups:
                safeName = annotationGroup.getName().replace(' ', '_')
                outstream.write('SCALARS ' + safeName + ' int 1\n')
                outstream.write('LOOKUP_TABLE default\n') 
                meshGroup = annotationGroup.getMeshGroup(mesh)
                elementIter = mesh.createElementiterator()
                element = elementIter.next()
                while element.isValid():
                    outstream.write('1 ' if meshGroup.containsElement(element) else '0 ')
                    element = elementIter.next()
                outstream.write('\n')


    def writeFile(self, filename):
        '''
        Export to legacy vtk file/
        '''
        with open(filename, 'w') as outstream:
            self._write(outstream)
