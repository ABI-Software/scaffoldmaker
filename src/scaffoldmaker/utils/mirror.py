'''
Utilities for getting a mirror image of a scaffold.
'''

from __future__ import division

from scaffoldmaker.utils import vector


class Mirror:
    '''
    Utilities for getting a mirror image of a scaffold.
    '''
    def __init__(self, mirrorPlane):
        '''
        :param mirrorPlane: Plane ax+by+cd=d defined as a list p=[a,b,c,d].
        '''
        self._plane = mirrorPlane
        self._planeUnitNormalVector = self.getPlaneNormalVector()
        self._magNormal = vector.magnitude(mirrorPlane[:3])
        self._planeDParam = mirrorPlane[3]

    def getPointDistanceFromPlane(self, x):
        """
        Get distance from point x to plane p. Returns negative distance, if point x is on the negative side of the plane.
        :param x: Point coordinates.
        :return: Distance.
        """
        n = self._planeUnitNormalVector
        ma = self._magNormal
        d = self._planeDParam
        return vector.dotproduct(x, n) - d / ma

    def getPlaneNormalVector(self):
        """
        Get plane normal vector
        :return: Unit normal vector.
        """
        p = self._plane
        return vector.normalise([p[0], p[1], p[2]])

    def pointProjectionToPlane(self, x):
        """
        Find the projection of point x onto plane p
        :param x:  Point coordinates.
        :return: x0, the point projection onto the plane
        """
        D = self.getPointDistanceFromPlane(x)
        n = self._planeUnitNormalVector
        return [x[c] - D*n[c] for c in range(3)]

    def mirrorImageOfPoint(self, x):
        """
        Find the mirror image of point x on plane p
        :param x:  Point coordinates.
        :return: xm, mirror image of point x.
        """
        D = self.getPointDistanceFromPlane(x)
        n = self._planeUnitNormalVector
        return [x[c] - 2*D*n[c] for c in range(3)]

    def mirrorVector(self,v):
        """
        Find the mirror image of vector v on plane p
        :param v: Vector.
        :return: image vector.
        """
        n = self._planeUnitNormalVector
        return [(v[c] - 2 * vector.dotproduct(v, n) * n[c]) for c in range(3)]

    def reverseMirrorVector(self,v):
        """
        Find the reverse mirror image of vector v on plane p
        :param v: Vector.
        :return: image vector.
        """
        n = self._planeUnitNormalVector
        return [-v[c] + 2 * vector.dotproduct(v, n) * n[c] for c in range(3)]
