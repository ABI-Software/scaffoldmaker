"""
Describes subdomains of a scaffold with attached names and terms.
"""
import sys

from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.utils.zinc.field import find_or_create_field_coordinates, find_or_create_field_group, \
    find_or_create_field_stored_mesh_location, find_or_create_field_stored_string
from opencmiss.zinc.element import Element, Mesh
from opencmiss.zinc.field import Field, FieldFiniteElement, FieldGroup, FieldStoredMeshLocation, FieldStoredString
from opencmiss.zinc.fieldmodule import Fieldmodule
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.utils.zinc_utils import get_highest_dimension_mesh, get_next_unused_node_identifier, \
    group_get_highest_dimension, identifier_ranges_from_string, identifier_ranges_to_string, \
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
        self._group = find_or_create_field_group(region.getFieldmodule(), self._name, managed=False)
        assert self._group.getName() == self._name, \
            'AnnotationGroup found existing non-group field called ' + self._name
        self._isMarker = False
        self._markerGroup = None  # held while marker point exists
        self._markerMaterialCoordinatesField = None

    def clear(self):
        """
        Manual clean up of user defined annotation group. See use in ScaffoldPackage.
        If this is a marker, destroys the marker node.
        Finally clears group
        """
        fieldmodule = self._group.getFieldmodule()
        with ChangeManager(fieldmodule):
            if self._isMarker:
                nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
                markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
                nodes.destroyNode(markerNode)
                self._isMarker = False
                self._markerGroup = None
                self._markerMaterialCoordinatesField = None
            self._group.clear()

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
        if self._isMarker:
            dct['marker'] = self._markerToDict()
        return dct

    def _markerToDict(self):
        """
        If the node is in the marker group, return a dict mapping "marker_location" to [elementId, xi] and if another
        entry is present it is a map from the materialCoordinatesField name to the coordinates value. Eg.
        {
            "heart coordinates": [11.5, 16.2, 7.6],
            "marker_location": [5, [0.0, 0.0, 1.0]]
        }
        :return: A python dict as above or None if not a marker point.
        """
        assert self._isMarker
        markerDct = {}
        if self._markerMaterialCoordinatesField:
            markerDct[self._markerMaterialCoordinatesField.getName()] = self.getMarkerMaterialCoordinates()[1]
        element, xi = self.getMarkerLocation()
        markerDct['marker_location'] = [element.getIdentifier(), xi]
        return markerDct

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
        markerDct = dct.get('marker')
        identifierRanges = identifier_ranges_from_string(identifierRangesString)
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            annotationGroup = cls(region, (name, ontId))
            annotationGroup._group.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
            if dimension > 0:
                meshGroup = annotationGroup.getMeshGroup(fieldmodule.findMeshByDimension(dimension))
                mesh_group_add_identifier_ranges(meshGroup, identifierRanges)
            else:
                if markerDct:
                    nodeIdentifier = int(identifierRanges[0][0])
                    annotationGroup._markerFromDict(nodeIdentifier, markerDct)
                nodesetGroup = annotationGroup.getNodesetGroup(
                    fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES))
                nodeset_group_add_identifier_ranges(nodesetGroup, identifierRanges)
        return annotationGroup

    def _markerFromDict(self, nodeIdentifier, markerDct: dict):
        """
        Define a new marker point node with nodeIdentifier and fields defined as in markerDct.
        Warn but do nothing if node with that identifier exists.
        :param nodeIdentifier:
        :param markerDct: A dict mapping names of fields to values, e.g.
        {
            "heart coordinates": [11.5, 16.2, 7.6],
            "marker_location": [5, [0.0, 0.0, 1.0]]
        }
        """
        assert not self._isMarker
        fieldmodule = self._group.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        node = nodes.findNodeByIdentifier(nodeIdentifier)
        if node.isValid():
            print("Error: Read annotation group " + self._name + ".  "\
                  "Marker point node " + str(nodeIdentifier) + " already exists", file=sys.stderr)
            return
        mesh = get_highest_dimension_mesh(fieldmodule)
        if not mesh:
            print("Error: Read annotation group " + self._name + ".  Marker cannot be made for empty mesh",
                  file=sys.stderr)
            return
        with ChangeManager(fieldmodule):
            materialCoordinatesField = None
            materialCoordinates = None
            element = None
            xi = None
            for key in markerDct.keys():
                if key != "marker_location":
                    fieldName = key
                    materialCoordinatesField = fieldmodule.findFieldByName(fieldName).castFiniteElement()
                    materialCoordinates = markerDct[fieldName]
                    if materialCoordinatesField.isValid() and isinstance(materialCoordinates, list):
                        break;
                    print("Error: Read annotation group " + self._name + ".  " \
                          "Invalid marker material coordinates field " + fieldName + " or value.", file=sys.stderr)
                    materialCoordinatesField = None
                    materialCoordinates = None
                    break
            if not materialCoordinatesField:
                element_xi = markerDct.get("marker_location")
                if not element_xi:
                    print("Error: Read annotation group " + self._name + ".  " \
                          "Marker missing marker_location entry", file=sys.stderr)
                    return
                elementIdentifier, xi = element_xi
                element = mesh.findElementByIdentifier(elementIdentifier)
                if not element.isValid():
                    print("Error: Read annotation group " + self._name + ".  " \
                          "Marker element " + str(elementIdentifier) + " not found", file=sys.stderr)
                    return
            self.createMarkerNode(nodeIdentifier, materialCoordinatesField, materialCoordinates, element, xi)

    def isMarker(self):
        """
        Query if this annotation group is a created marker node.
        """
        return self._isMarker

    def createMarkerNode(self, startNodeIdentifier=1,
                         materialCoordinatesField: FieldFiniteElement=None, materialCoordinates=None,
                         element=None, xi=[0.0, 0.0, 0.0]):
        """
        Convert annotation group into a marker point annotation.
        Important: annotation group must currently be empty, and elements must exist.
        The marker node is added to the marker group in addition to this group.
        Raises an exception if marker creation cannot be achieved.
        :param startNodeIdentifier: First unused node identifier >= this may use for marker node. Incremented until
        an unused node identifier is found for the marker node.
        :param materialCoordinatesField: Material coordinates field to define location of marker point in.
        Must be a finite element type field for which isTypeCoordinates() is True, with up to 3 components,
        and at least as many components as the highest mesh dimension.
        Only one of materialCoordinatesField or element may be specified.
        :param materialCoordinates: The coordinates to assign to the materialCoordinatesField, if field supplied.
        The element:xi coordinates are computed from it.
        :param element: Optional element to set initial location from; must be in highest dimension mesh.
        Only one of materialCoordinatesField or element may be specified.
        If neither is specified the first element in the highest dimension mesh is used.
        :param xi: If element is supplied or first is assumed, the xi coordinates to embed at.
        :return: Zinc Node representing marker point, with either default location in first element or that supplied.
        """
        assert not self._isMarker
        assert self._group.isEmpty(), "Can only create marker in empty annotation group"
        fieldmodule = self._group.getFieldmodule()
        mesh = get_highest_dimension_mesh(fieldmodule)
        assert mesh, "Can only create marker point if there is a mesh with elements"
        assert not (element and materialCoordinatesField)
        if materialCoordinatesField:
            coordinatesCount = materialCoordinatesField.getNumberOfComponents()
            assert materialCoordinates and (len(materialCoordinates) >= coordinatesCount) and \
                   materialCoordinatesField.castFiniteElement().isValid() and \
                   materialCoordinatesField.isTypeCoordinate() and (mesh.getDimension() <= coordinatesCount <= 3)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = get_next_unused_node_identifier(nodes, startNodeIdentifier)
        with ChangeManager(fieldmodule):
            markerGroup = getAnnotationMarkerGroup(fieldmodule)
            markerNodeGroup = markerGroup.getFieldNodeGroup(nodes)
            if not markerNodeGroup.isValid():
                markerNodeGroup = markerGroup.createFieldNodeGroup(nodes)
            markerNodesetGroup = markerNodeGroup.getNodesetGroup()
            markerLocation = getAnnotationMarkerLocationField(fieldmodule, mesh)
            markerName = getAnnotationMarkerNameField(fieldmodule)
            nodetemplate = nodes.createNodetemplate()
            if materialCoordinatesField:
                assert RESULT_OK == nodetemplate.defineField(materialCoordinatesField)
            assert RESULT_OK == nodetemplate.defineField(markerLocation)
            assert RESULT_OK == nodetemplate.defineField(markerName)
            markerNode = nodes.createNode(nodeIdentifier, nodetemplate)
            assert RESULT_OK == self.getNodesetGroup(nodes).addNode(markerNode)
            assert RESULT_OK == markerNodesetGroup.addNode(markerNode)
            self._isMarker = True
            self._markerGroup = markerGroup  # maintain reference so group is not destroyed
            # assign marker name to be same as this group's name. This needs to be maintained. See setName()
            fieldcache = fieldmodule.createFieldcache()
            fieldcache.setNode(markerNode)
            markerName.assignString(fieldcache, self._name)
            if materialCoordinatesField:
                self.setMarkerMaterialCoordinates(materialCoordinatesField, materialCoordinates)
            else:
                if not element:
                    element = mesh.createElementiterator().next()
                self.setMarkerLocation(element, xi)
        return markerNode

    def getMarkerLocation(self):
        """
        If the annotation is a created marker point, get its element:xi location.
        :return: Zinc Element, xi (list of float)
        """
        assert self._isMarker
        fieldmodule = self._group.getFieldmodule()
        mesh = get_highest_dimension_mesh(fieldmodule)
        markerLocation = getAnnotationMarkerLocationField(fieldmodule, mesh)
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
        fieldcache.setNode(markerNode)
        element, xi = markerLocation.evaluateMeshLocation(fieldcache, mesh.getDimension())
        if not isinstance(xi, list):
            xi = [xi]  # workaround for Zinc 1-D xi being a plain float
        return element, xi

    def setMarkerLocation(self, element: Element, xi: list):
        """
        If the annotation group is a created marker point, set its element:xi location.
        If the marker also has a materialCoordinatesField, it is updated to the value
        at the supplied element:xi location.
        :param element: Element to set location in
        :param xi: xi coordinates (list of float)
        """
        assert self._isMarker
        fieldmodule = self._group.getFieldmodule()
        mesh = get_highest_dimension_mesh(fieldmodule)
        assert mesh.containsElement(element), "Invalid element, not in highest dimension mesh"
        assert isinstance(xi, list) and (len(xi) >= mesh.getDimension())
        with ChangeManager(fieldmodule):
            markerLocation = getAnnotationMarkerLocationField(fieldmodule, mesh)
            fieldcache = fieldmodule.createFieldcache()
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
            fieldcache.setNode(markerNode)
            markerLocation.assignMeshLocation(fieldcache, element, xi)
            if self._markerMaterialCoordinatesField:
                fieldcache.setMeshLocation(element, xi)
                result, materialCoordinates = self._markerMaterialCoordinatesField.evaluateReal(
                    fieldcache, self._markerMaterialCoordinatesField.getNumberOfComponents())
                fieldcache.setNode(markerNode)
                self._markerMaterialCoordinatesField.assignReal(fieldcache, materialCoordinates)

    def getMarkerMaterialCoordinates(self):
        """
        If the annotation is a created marker point, get its material coordinates field and location, if set.
        :return: materialCoordinatesField, materialCoordinates or None, None
        """
        assert self._isMarker
        if not self._markerMaterialCoordinatesField:
            return None, None
        fieldmodule = self._group.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
        fieldcache.setNode(markerNode)
        result, materialCoordinates = self._markerMaterialCoordinatesField.evaluateReal(
            fieldcache, self._markerMaterialCoordinatesField.getNumberOfComponents())
        return self._markerMaterialCoordinatesField, materialCoordinates

    def setMarkerMaterialCoordinates(self, materialCoordinatesField, materialCoordinates=None):
        """
        Also updates the marker location when this is assigned, forcing it to be within the mesh.
        The material coordinates are then recalculated to be in the mesh.
        Some approximations may occur if the point is outside the mesh - user beware.
        :param materialCoordinatesField: Material coordinates field to set or change to, or None to remove
        material coordinates field so marker point only has an element:xi location. Must be defined on elements
        of highest dimension mesh.
        :param materialCoordinates: If None, evaluate materialCoordinatesField at current marker location.
        """
        assert self._isMarker
        if not (self._markerMaterialCoordinatesField or materialCoordinatesField):
            return
        fieldmodule = self._group.getFieldmodule()
        mesh = get_highest_dimension_mesh(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
        with ChangeManager(fieldmodule):
            oldMaterialCoordinatesField = self._markerMaterialCoordinatesField
            if materialCoordinatesField != oldMaterialCoordinatesField:
                nodetemplate = nodes.createNodetemplate()
                if self._markerMaterialCoordinatesField:
                    assert RESULT_OK == nodetemplate.undefineField(self._markerMaterialCoordinatesField)
                if materialCoordinatesField:
                    assert RESULT_OK == nodetemplate.defineField(materialCoordinatesField)
                assert RESULT_OK == markerNode.merge(nodetemplate)
                self._markerMaterialCoordinatesField = materialCoordinatesField
            if materialCoordinatesField:
                if materialCoordinates is not None:
                    # find nearest location in mesh
                    constCoordinates = fieldmodule.createFieldConstant(materialCoordinates)
                    findMeshLocation = fieldmodule.createFieldFindMeshLocation(
                        constCoordinates, materialCoordinatesField, mesh)
                    fieldcache = fieldmodule.createFieldcache()
                    fieldcache.setNode(markerNode)
                    element, xi = findMeshLocation.evaluateMeshLocation(fieldcache, mesh.getDimension())
                    if not element.isValid():
                        self.setMarkerMaterialCoordinates(oldMaterialCoordinatesField)
                        print("AnnotationGroup setMarkerMaterialCoordinates.  Field is not defined on mesh. Reverting.")
                        return
                    if not isinstance(xi, list):
                        xi = [xi]  # workaround for Zinc 1-D xi being a plain float
                    del findMeshLocation
                    del constCoordinates
                else:
                    element, xi = self.getMarkerLocation()
                # use this function to reassign material coordinates to be within mesh
                self.setMarkerLocation(element, xi)

    def getName(self):
        return self._name

    def setName(self, name):
        '''
        Client must ensure name is unique for all annotation groups.
        First tries to rename zinc group field; if that fails, it won't rename group
        as the name is already in use.
        :return:  True on success, otherwise False
        '''
        if self._name == name:
            return True
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
                if self._isMarker:
                    # assign marker name to be same as this group's name. This needs to be maintained
                    markerName = getAnnotationMarkerNameField(fieldmodule)
                    fieldcache = fieldmodule.createFieldcache()
                    markerNode = self.getNodesetGroup(nodes).createNodeiterator().next()
                    fieldcache.setNode(markerNode)
                    markerName.assignString(fieldcache, self._name)
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
        assert annotationGroup._id == term[1], "Annotation group '" + name + "' id '" + term[1]\
                                               + "' does not match existing id '" + annotationGroup._id + "'"
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
        assert annotationGroup._id == term[1], "Annotation group '" + name + "' id '" + term[1]\
                                               + "' does not match existing id '" + annotationGroup._id + "'"
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

def getAnnotationMarkerGroup(fieldmodule: Fieldmodule) -> FieldGroup:
    """
    Find or create the standard Zinc Group which marker points are created in.
    Clients should use this method rather than finding by standard name "marker".
    :param fieldmodule: Zinc Fieldmodule to find/create group in.
    :return: FieldGroup
    """
    return find_or_create_field_group(fieldmodule, "marker", managed=False)

def getAnnotationMarkerLocationField(fieldmodule: Fieldmodule, mesh: Mesh) -> FieldStoredMeshLocation:
    """
    Find or create the standard Zinc Field used to store marker mesh locations in the highest dimension mesh.
    Clients should use this method rather than finding by standard name "marker_location".
    :param fieldmodule: Zinc Fieldmodule to find/create field in.
    :param mesh: The highest dimension mesh in region.
    :return: FieldStoredMeshLocation
    """
    return find_or_create_field_stored_mesh_location(fieldmodule, mesh, name="marker_location", managed=False)

def getAnnotationMarkerNameField(fieldmodule: Fieldmodule) -> FieldStoredString:
    """
    Find or create the standard Zinc Field used to store marker names.
    Clients should use this method rather than finding by standard name "marker_name".
    :param fieldmodule: Zinc Fieldmodule to find/create group in.
    :return: FieldStoredString
    """
    return find_or_create_field_stored_string(fieldmodule, name="marker_name", managed=False)


