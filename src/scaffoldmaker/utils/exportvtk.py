'''
Class for exporting a Scaffold from Zinc to legacy vtk text format.
'''

import io
import os
import sys
from sys import version_info

from opencmiss.utils.zinc.finiteelement import getElementNodeIdentifiersBasisOrder
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field
from opencmiss.zinc.result import RESULT_OK


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
        self._nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._coordinates = self._fieldmodule.findFieldByName('coordinates').castFiniteElement()
        self._description = description
        self._annotationGroups = annotationGroups if annotationGroups else []
        self._markerNodes = None
        markerGroup = self._fieldmodule.findFieldByName("marker")
        if markerGroup.isValid():
            markerGroup = markerGroup.castGroup()
            markerNodeGroup = markerGroup.getFieldNodeGroup(self._nodes)
            if markerNodeGroup.isValid():
                self._markerNodes = markerNodeGroup.getNodesetGroup()

    def _write(self, outstream):
        if version_info.major > 2:
          assert isinstance(outstream, io.TextIOBase), 'ExportVtk.write:  Invalid outstream argument'
        outstream.write('# vtk DataFile Version 2.0\n')
        outstream.write(self._description + '\n')
        outstream.write('ASCII\n')
        outstream.write('DATASET UNSTRUCTURED_GRID\n')
        nodeIdentifierToIndex = {}  # map needed since vtk points are zero index based, i.e. have no identifier
        coordinatesCount = self._coordinates.getNumberOfComponents()
        cache = self._fieldmodule.createFieldcache()

        # exclude marker nodes from output
        pointCount = self._nodes.getSize()
        if self._markerNodes:
            pointCount -= self._markerNodes.getSize()
        outstream.write('POINTS ' + str(pointCount) + ' double\n')
        nodeIter = self._nodes.createNodeiterator()
        node = nodeIter.next()
        index = 0
        while node.isValid():
            if not (self._markerNodes and self._markerNodes.containsNode(node)):
                nodeIdentifierToIndex[node.getIdentifier()] = index
                cache.setNode(node)
                result, x = self._coordinates.evaluateReal(cache, coordinatesCount)
                if result != RESULT_OK:
                    print("Coordinates not found for node", node.getIdentifier())
                    x = [ 0.0 ]*coordinatesCount
                if coordinatesCount < 3:
                    for c in range(coordinatesCount - 1, 3):
                        x.append(0.0)
                outstream.write(" ".join(str(s) for s in x) + "\n")
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
            eft = element.getElementfieldtemplate(self._coordinates, -1)  # assumes all components same
            nodeIdentifiers = getElementNodeIdentifiersBasisOrder(element, eft)
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

        # use cell data for annotation groups containing elements of mesh dimension
        # use point data for lower dimensional annotation groups
        pointAnnotationGroups = []
        cellAnnotationGroups = []
        for annotationGroup in self._annotationGroups:
            if annotationGroup.hasMeshGroup(self._mesh):
                cellAnnotationGroups.append(annotationGroup)
            elif annotationGroup.hasNodesetGroup(self._nodes):
                pointAnnotationGroups.append(annotationGroup)

        if cellAnnotationGroups:
            outstream.write('CELL_DATA ' + str(cellCount) + '\n')
            for annotationGroup in cellAnnotationGroups:
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

        if pointAnnotationGroups:
            outstream.write('POINT_DATA ' + str(pointCount) + '\n')
            for annotationGroup in pointAnnotationGroups:
                safeName = annotationGroup.getName().replace(' ', '_')
                outstream.write('SCALARS ' + safeName + ' int 1\n')
                outstream.write('LOOKUP_TABLE default\n') 
                nodesetGroup = annotationGroup.getNodesetGroup(self._nodes)
                nodeIter = self._nodes.createNodeiterator()
                node = nodeIter.next()
                while node.isValid():
                    if not (self._markerNodes and self._markerNodes.containsNode(node)):
                        outstream.write('1 ' if nodesetGroup.containsNode(node) else '0 ')
                    node = nodeIter.next()
                outstream.write('\n')


    def _writeMarkers(self, outstream):
        coordinatesCount = self._coordinates.getNumberOfComponents()
        cache = self._fieldmodule.createFieldcache()
        markerLocation = self._fieldmodule.findFieldByName("marker_location")
        markerName = self._fieldmodule.findFieldByName("marker_name")
        if markerLocation.isValid() and markerName.isValid():
            with ChangeManager(self._fieldmodule):
                markerCoordinates = self._fieldmodule.createFieldEmbedded(self._coordinates, markerLocation)
                nodeIter = self._markerNodes.createNodeiterator()
                node = nodeIter.next()
                while node.isValid():
                    cache.setNode(node)
                    result, x = markerCoordinates.evaluateReal(cache, coordinatesCount)
                    if result == RESULT_OK:
                        name = markerName.evaluateString(cache)
                        outstream.write(",".join(str(s) for s in x) + "," + name + "\n")
                    node = nodeIter.next()
                del markerCoordinates


    def writeFile(self, filename):
        '''
        Export to legacy vtk file.
        '''
        try:
            with open(filename, 'w') as outstream:
                self._write(outstream)
            if self._markerNodes and (self._markerNodes.getSize() > 0):
                markerFilename = os.path.splitext(filename)[0] + "_marker.csv"
                with open(markerFilename, 'w') as outstream:
                    self._writeMarkers(outstream)
        except Exception as e:
            print("Failed to write VTK file", filename, file=sys.stderr);
