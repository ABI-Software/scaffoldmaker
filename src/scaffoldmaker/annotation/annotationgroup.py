"""
Describes subdomains of a scaffold with attached names and terms.
"""

from opencmiss.zinc.field import FieldGroup

class AnnotationGroup(object):
    '''
    Describes subdomains of a scaffold with attached names and terms.
    '''

    def __init__(self, region, name, FMANumber, lyphID):
        '''
        :param region: The Zinc region the AnnotationGroup is to be made for.
        :param name: The name of an annotation group e.g. common medical term.
        :param FMANumber: The FMA Number of the group.
        :param lyphID: The Apinatomy Lyph ID for the group.
        '''
        self._name = name
        self._FMANumber = FMANumber
        self._lyphID = lyphID
        fm = region.getFieldmodule()
        field = fm.findFieldByName(name)
        if field.isValid():
            self._group = field.castGroup()
            assert self._group.isValid(), 'AnnotationGroup found existing non-group field called ' + name
        else:
            # assume client is calling between fm.begin/endChange()
            self._group = fm.createFieldGroup()
            self._group.setName(name)
            self._group.setManaged(True)

    def getName(self):
        return self._name

    def getFMANumber(self):
        return self._FMANumber

    def getLyphID(self):
        return self._lyphID

    def getGroup(self):
        return self._group

    def getFieldElementGroup(self, mesh):
        '''
        :param mesh: The Zinc mesh to manage a sub group of.
        :return: The Zinc element group field for mesh in this AnnotationGroup.
        '''
        elementGroup = self._group.getFieldElementGroup(mesh)
        if not elementGroup.isValid():
            elementGroup = self._group.createFieldElementGroup(mesh)
        return elementGroup

    def getFieldNodeGroup(self, nodeset):
        '''
        :param nodeset: The Zinc nodeset to manage a sub group of.
        :return: The Zinc node group field for nodeset in this AnnotationGroup.
        '''
        nodeGroup = self._group.getFieldNodeGroup(nodeset)
        if not nodeGroup.isValid():
            nodeGroup = self._group.createFieldNodeGroup(nodeset)
        return nodeGroup

    def getMeshGroup(self, mesh):
        '''
        :param mesh: The Zinc mesh to manage a sub group of.
        :return: The Zinc meshGroup for adding elements of mesh in this AnnotationGroup.
        '''
        return self.getFieldElementGroup(mesh).getMeshGroup()

    def getNodesetGroup(self, nodeset):
        '''
        :param nodeset: The Zinc nodeset to manage a sub group of.
        :return: The Zinc nodesetGroup for adding nodes from nodeset in this AnnotationGroup.
        '''
        return self.getFieldNodeGroup(nodeset).getNodesetGroup()

    def addSubelements(self):
        '''
        Call after group is complete and faces have been defined to add faces and
        nodes for elements in group to related subgroups.
        '''
        self._group.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
        fm = self._group.getFieldmodule()
        for dimension in range(1, 4):
            mesh = fm.findMeshByDimension(dimension)
            elementGroup = self._group.getFieldElementGroup(mesh)
            if elementGroup.isValid():
                meshGroup = elementGroup.getMeshGroup()
                #print('Mesh group:', self._name, ', size', meshGroup.getSize())
                meshGroup.addElementsConditional(elementGroup)  # use FieldElementGroup as conditional field


def findAnnotationGroupByName(annotationGroups, name):
    '''
    :param annotationGroups: list(AnnotationGroup)
    :param name: Name of group.
    :return: AnnotationGroup or None.
    '''
    for annotationGroup in annotationGroups:
        if annotationGroup._name == name:
            return annotationGroup
    return None
