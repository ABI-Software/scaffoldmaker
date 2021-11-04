"""
Describes subdomains of a scaffold with attached names and terms.
"""

from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.field import Field, FieldGroup
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.utils.zinc_utils import group_get_highest_dimension, \
    identifier_ranges_from_string, identifier_ranges_to_string, \
    mesh_group_add_identifier_ranges, mesh_group_to_identifier_ranges, \
    nodeset_group_add_identifier_ranges, nodeset_group_to_identifier_ranges


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
        fieldmodule = region.getFieldmodule()
        field = fieldmodule.findFieldByName(self._name)
        if field.isValid():
            self._group = field.castGroup()
            assert self._group.isValid(), 'AnnotationGroup found existing non-group field called ' + self._name
        else:
            with ChangeManager(fieldmodule):
                self._group = fieldmodule.createFieldGroup()
                self._group.setName(self._name)
                self._group.setManaged(True)

    def toDict(self):
        '''
        Encodes object into a dictionary for JSON serialisation.
        Used only for user-defined annotation groups.
        :return: Dictionary containing object encoding.
        '''
        # get identifier ranges from highest dimension domain in group
        dimension = self.getDimension()
        fieldmodule = self._group.getFieldmodule()
        if dimension > 0:
            mesh = fieldmodule.findMeshByDimension(dimension)
            meshGroup = self._group.getFieldElementGroup(mesh).getMeshGroup()
            identifierRanges = mesh_group_to_identifier_ranges(meshGroup)
        else:
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            nodesetGroup = self._group.getFieldNodeGroup(nodes).getNodesetGroup()
            if nodesetGroup.isValid():
                identifierRanges = nodeset_group_to_identifier_ranges(nodesetGroup)
            else:
                identifierRanges = []
        dct = {
            '_AnnotationGroup' : True,
            'name' : self._name,
            'ontId' : self._id,
            'dimension' : dimension,
            'identifierRanges' : identifier_ranges_to_string(identifierRanges)
            }
        return dct

    @classmethod
    def fromDict(cls, dct, region):
        '''
        Instantiate from dict. See toDict()
        :param region: Zinc region.
        :return: AnnotationGroup
        '''
        assert dct['_AnnotationGroup']
        name = dct['name']
        ontId = dct['ontId']
        dimension = dct['dimension']
        identifierRangesString = dct['identifierRanges']
        identifierRanges = identifier_ranges_from_string(identifierRangesString)
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            annotationGroup = cls(region, (name, ontId))
            annotationGroup._group.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
            if dimension > 0:
                meshGroup = annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(dimension))
                mesh_group_add_identifier_ranges(meshGroup, identifierRanges)
            else:
                nodesetGroup = annotationGroup.getNodesetGroup(fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES))
                nodeset_group_add_identifier_ranges(nodesetGroup, identifierRanges)
        return annotationGroup

    def getName(self):
        return self._name

    def setName(self, name):
        '''
        Client must ensure name is unique for all annotation groups.
        First tries to rename zinc group field; if that fails, it won't rename group
        as the name is already in use.
        :return:  True on success, otherwise False
        '''
        fieldmodule = self._group.getFieldmodule()
        # use ChangeManager so multiple name changes are atomic
        with ChangeManager(fieldmodule):
            if RESULT_OK == self._group.setName(name):
                # workaround for zinc issue: must rename subelement groups
                for dimension in range(3, 0, -1):
                    mesh = fieldmodule.findMeshByDimension(dimension)
                    elementGroup = self._group.getFieldElementGroup(mesh)
                    if elementGroup.isValid():
                        elementGroup.setName(name + '.' + mesh.getName())
                nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
                nodeGroup = self._group.getFieldNodeGroup(nodes)
                if nodeGroup.isValid():
                    nodeGroup.setName(name + '.' + nodes.getName())
                self._name = name
                return True
        return False

    def getId(self):
        return self._id

    def setId(self, id):
        '''
        Client must ensure id is unique for all annotation groups.
        :return:  True on success, otherwise False
        '''
        self._id = id
        return True

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

    def getDimension(self):
        '''
        Get dimension 3, 2 or 1 of mesh which group is annotating, 0 if nodes or -1 if empty.
        '''
        return group_get_highest_dimension(self._group)

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
    Find existing annotation group for term, or create it for region if not found.
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
