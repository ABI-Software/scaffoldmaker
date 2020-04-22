"""
Describes subdomains of a scaffold with attached names and terms.
"""

from opencmiss.zinc.field import FieldGroup

class AnnotationGroup(object):
    '''
    Describes subdomains of a scaffold with attached names and terms.
    '''

    def __init__(self, region, term):
        '''
        :param region: The Zinc region the AnnotationGroup is to be made for.
        :param term: Identifier for anatomical term, currently a tuple of name, id.
        e.g. ("heart", "FMA:7088")
        '''
        self._name = term[0]
        self._id = term[1]
        fm = region.getFieldmodule()
        field = fm.findFieldByName(self._name)
        if field.isValid():
            self._group = field.castGroup()
            assert self._group.isValid(), 'AnnotationGroup found existing non-group field called ' + self._name
        else:
            # assume client is calling between fm.begin/endChange()
            self._group = fm.createFieldGroup()
            self._group.setName(self._name)
            self._group.setManaged(True)

    def getName(self):
        return self._name

    def getId(self):
        return self._id

    def getFMANumber(self):
        """
        :return: Integer FMA number or None.
        """
        if self._id and (self.id[:4] == "FMA:"):
            return int(self._id[4:])
        return None

    def getTerm(self):
        return ( self._name, self._id )

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

    def hasMeshGroup(self, mesh):
        '''
        :param mesh: The Zinc mesh to query a sub group of.
        :return: True if MeshGroup for mesh exists and is not empty, otherwise False.
        '''
        elementGroup = self._group.getFieldElementGroup(mesh)
        return elementGroup.isValid() and (elementGroup.getMeshGroup().getSize() > 0)

    def getNodesetGroup(self, nodeset):
        '''
        :param nodeset: The Zinc nodeset to manage a sub group of.
        :return: The Zinc nodesetGroup for adding nodes from nodeset in this AnnotationGroup.
        '''
        return self.getFieldNodeGroup(nodeset).getNodesetGroup()

    def hasNodesetGroup(self, nodeset):
        '''
        :param nodeset: The Zinc nodeset to query a sub group of.
        :return: True if NodesetGroup for nodeset exists and is not empty, otherwise False.
        '''
        nodeGroup = self._group.getFieldNodeGroup(nodeset)
        return nodeGroup.isValid() and (nodeGroup.getNodesetGroup().getSize() > 0)

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


def findAnnotationGroupByName(annotationGroups: list, name: str):
    '''
    Find existing annotation group for name.
    :param annotationGroups: list(AnnotationGroup)
    :param name: Name of group.
    :return: AnnotationGroup or None if not found.
    '''
    for annotationGroup in annotationGroups:
        if annotationGroup._name == name:
            return annotationGroup
    return None


def findOrCreateAnnotationGroupForTerm(annotationGroups: list, region, term) -> AnnotationGroup:
    '''
    Find existing annotation group for term, or create it for region if not gound.
    If annotation group created here, append it to annotationGroups.
    :param annotationGroups: list(AnnotationGroup)
    :param region: Zinc region to create group for.
    :param term: Identifier for anatomical term, currently a tuple of name, id.
    :return: AnnotationGroup.
    '''
    name = term[0]
    annotationGroup = findAnnotationGroupByName(annotationGroups, name)
    if annotationGroup:
        assert annotationGroup._id == term[1], "Annotation group '" + name + "' id '" + term[1] + "' does not match existing id '" + annotationGroup._id + "'"
    else:
        annotationGroup = AnnotationGroup(region, term)
        annotationGroups.append(annotationGroup)
    return annotationGroup


def getAnnotationGroupForTerm(annotationGroups: list, term) -> AnnotationGroup:
    '''
    Get existing annotation group for term. Raise exception if not found.
    :param annotationGroups: list(AnnotationGroup)
    :param term: Identifier for anatomical term, currently a tuple of name, id.
    :return: AnnotationGroup.
    '''
    name = term[0]
    annotationGroup = findAnnotationGroupByName(annotationGroups, name)
    if annotationGroup:
        assert annotationGroup._id == term[1], "Annotation group '" + name + "' id '" + term[1] + "' does not match existing id '" + annotationGroup._id + "'"
        return annotationGroup
    raise NameError("Annotation group '" + name + "' not found.")


def mergeAnnotationGroups(*annotationGroupsIn):
    '''
    Merge the supplied sequence of list(annotationGroups) to a single list,
    without duplicates.
    :param annotationGroupsIn: Variable number of list(AnnotationGroup) to merge.
     Groups must be for the same region.
    :return: Merged list(AnnotationGroup)
    '''
    annotationGroups = []
    for agroups in annotationGroupsIn:
        for agroup in agroups:
            if not findAnnotationGroupByName(annotationGroups, agroup._name):
                annotationGroups.append(agroup)
    return annotationGroups
