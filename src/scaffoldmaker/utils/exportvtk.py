'''
Class for exporting a Scaffold from Zinc to legacy vtk text format.
'''

import io
from sys import version_info
from scaffoldmaker.utils import zinc_utils
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
        for dimension in range(3, 1, -1):
            self._mesh = self._fieldmodule.findMeshByDimension(dimension)
            if self._mesh.getSize() > 0:
                break
        self._coordinates = self._fieldmodule.findFieldByName('coordinates')
        self._description = description
        self._annotationGroups = annotationGroups if annotationGroups else []


    def _write(self, outstream):
        if version_info.major > 2:
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

        # following assumes all hex (3-D) or all quad (2-D) elements
        if self._mesh.getDimension() == 2:
            localNodeCount = 4
            vtkIndexing = [ 0, 1, 3, 2 ]
            cellTypeString = '9'
        else:
            localNodeCount = 8
            vtkIndexing = [ 0, 1, 3, 2, 4, 5, 7, 6 ]
            cellTypeString = '12'
        localNodeCountStr = str(localNodeCount)
        cellCount = self._mesh.getSize()
        cellListSize = (1 + localNodeCount)*cellCount
        outstream.write('CELLS ' + str(cellCount) + ' ' + str(cellListSize) + '\n')
        elementIter = self._mesh.createElementiterator()
        element = elementIter.next()
        while element.isValid():
            eft = element.getElementfieldtemplate(coordinates, -1)  # assumes all components same
            if localNodeCount == 4:
                nodeIdentifiers = zinc_utils.getElementNodeIdentifiers4Node(element, eft)
            else:
                nodeIdentifiers = zinc_utils.getElementNodeIdentifiers8Node(element, eft)
            outstream.write(localNodeCountStr)
            for localIndex in vtkIndexing:
                index = nodeIdentifierToIndex[nodeIdentifiers[localIndex]]
                outstream.write(' ' + str(index))
            outstream.write('\n')
            element = elementIter.next()
        outstream.write('CELL_TYPES ' + str(cellCount) + '\n')
        for i in range(cellCount - 1):
            outstream.write(cellTypeString + ' ')
        outstream.write(cellTypeString + '\n')

        if self._annotationGroups:
            outstream.write('CELL_DATA ' + str(cellCount) + '\n')
            for annotationGroup in self._annotationGroups:
                safeName = annotationGroup.getName().replace(' ', '_')
                outstream.write('SCALARS ' + safeName + ' int 1\n')
                outstream.write('LOOKUP_TABLE default\n') 
                meshGroup = annotationGroup.getMeshGroup(self._mesh)
                elementIter = self._mesh.createElementiterator()
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
