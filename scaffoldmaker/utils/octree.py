'''
Octree for searching for objects by coordinates
Created on April 4, 2018

@author: Richard Christie
'''

import copy
import math

class Octree:
    '''
    Octree for searching for objects by coordinates
    '''

    def __init__(self, minimums, maximums, tolerance = None):
        '''
        :param minimums: List of 3 minimum coordinate values. Caller to include any edge allowance.
        :param maximums: List of 3 maximum coordinate values. Caller to include any edge allowance.
        :param tolerance: If supplied, tolerance to use, or None to compute as 1.0E-6*diagonal.
        '''
        self._dimension = 3
        self._maxObjects = 20
        assert len(minimums) == self._dimension, 'Octree minimums is invalid length'
        assert len(maximums) == self._dimension, 'Octree maximums is invalid length'
        if tolerance is None:
            self._tolerance = 1.0E-6*math.sqrt(sum(((maximums[i] - minimums[i])*(maximums[i] - minimums[i])) for i in range(self._dimension)))
        else:
            self._tolerance = tolerance
        self._minimums = copy.deepcopy(minimums)
        self._maximums = copy.deepcopy(maximums)
        # Octree is either leaf with _coordinatesObjects, or has 2**self._dimension children
        self._coordinatesObjects = []
        # exactly 2^self._dimension children, cycling in lowest x index fastest
        self._children = None
        #if tolerance is None:
        #print('Octree:', self._minimums, self._maximums, self._tolerance)


    def _isWithinBounds(self, x):
        '''
        :return: True if x is within bounds of minimums and maximums.
        '''
        for i in range(self._dimension):
            if (x[i] < self._minimums[i]) or (x[i] > self._maximums[i]):
                return False
        return True


    def _isWithinBoundsPlusTolerance(self, x):
        '''
        :return: True if x is within bounds of minimums and maximums plus tolerance.
        '''
        for i in range(self._dimension):
            if (x[i] < (self._minimums[i] - self._tolerance)) or \
                (x[i] > (self._maximums[i] + self._tolerance)):
                return False
        return True


    def _findObjectByCoordinates(self, x):
        '''
        Find closest existing object with |x - ox| < tolerance.
        :param x: 3 coordinates in a list.
        :return: nearest distance, nearest object or None, None if none found.
        '''
        if not self._isWithinBoundsPlusTolerance(x):
            return None, None
        nearestDistance = None
        nearestObject = None
        if self._coordinatesObjects is not None:
            for coordinatesObject in self._coordinatesObjects:
                distance = math.sqrt(sum((x[i] - coordinatesObject[0][i])*(x[i] - coordinatesObject[0][i]) for i in range(self._dimension)))
                if (distance < self._tolerance) and ((nearestDistance is None) or (distance < nearestDistance)):
                    nearestDistance = distance
                    nearestObject = coordinatesObject[1]
        else:
            for child in self._children:
                distance, obj = child._findObjectByCoordinates(x)
                if (distance is not None) and ((nearestDistance is None) or (distance < nearestDistance)):
                    nearestDistance = distance
                    nearestObject = obj
        return nearestDistance, nearestObject


    def findObjectByCoordinates(self, x):
        '''
        Find closest existing object with |x - ox| < tolerance.
        :param x: 3 coordinates in a list.
        :return: nearest object or None if not found.
        '''
        nearestDistance, nearestObject = self._findObjectByCoordinates(x)
        return nearestObject


    def addObjectAtCoordinates(self, x, obj):
        '''
        Add object at coordianates to octree.
        Caller must have received None result for findObjectByCoordinates() first!
        Assumes caller has verified x is within range of Octree.
        :param x: 3 coordinates in a list.
        :param obj: object to store with coordinates.
        '''
        if self._coordinatesObjects is not None:
            if len(self._coordinatesObjects) < self._maxObjects:
                self._coordinatesObjects.append( (copy.deepcopy(x), obj) )
                return
            else:
                # subdivide and add coordinatesObjects plus new object to new children
                coordinatesObjects = self._coordinatesObjects
                self._coordinatesObjects = None
                self._children = []
                for i in range(1 << self._dimension):
                    childMinimums = copy.deepcopy(self._minimums)
                    childMaximums = copy.deepcopy(self._maximums)
                    for j in range(self._dimension):
                        if i & (1 << j):
                            childMinimums[j] = 0.5*(self._minimums[j] + self._maximums[j])
                        else:
                            childMaximums[j] = 0.5*(self._minimums[j] + self._maximums[j])
                    child = Octree(childMinimums, childMaximums, self._tolerance)
                    self._children.append(child)
                    # add any coordinatesObjects to child, if in range
                    addedCoordinatesObjects = []
                    for coordinatesObject in coordinatesObjects:
                        if child._isWithinBounds(coordinatesObject[0]):
                            #print('+ add ', coordinatesObject[1])
                            child.addObjectAtCoordinates(coordinatesObject[0], coordinatesObject[1])
                            addedCoordinatesObjects.append(coordinatesObject)
                    coordinatesObjects = [ v for v in coordinatesObjects if v not in addedCoordinatesObjects ]
                assert (len(coordinatesObjects) == 0), 'Octree.addObjectAtCoordinates failed to add existing objects'
        # add the new object to the first child it fits in
        for child in self._children:
            if child._isWithinBounds(x):
                child.addObjectAtCoordinates(x, obj)
                return
