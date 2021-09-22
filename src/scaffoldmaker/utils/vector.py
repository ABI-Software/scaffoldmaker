'''
Utility functions for vectors.
'''

import math

def crossproduct3(a, b):
    '''
    :return: vector 3-D cross product of a and b
    '''
    return [ a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] ]

def dotproduct(a, b):
    '''
    :return: vector dot (inner) product of a and b
    '''
    return sum(a[i]*b[i] for i in range(len(a)))

def magnitude(v):
    '''
    return: scalar magnitude of vector v
    '''
    return math.sqrt(sum(c*c for c in v))

def normalise(v):
    '''
    :return: vector v normalised to unit length
    '''
    mag = math.sqrt(sum(c*c for c in v))
    return [ c/mag for c in v ]

def setMagnitude(v, mag):
    '''
    return: Vector v with magnitude set to mag.
    '''
    scale = mag/math.sqrt(sum(c*c for c in v))
    return [ c*scale for c in v ]

def addVectors(vectors, scalars = None):
    '''
    returns s1*v1+s2*v2+... where scalars = [s1, s2, ...] and vectors=[v1, v2, ...].
    :return: Resultant vector
    '''
    if not scalars:
        scalars = [1]*len(vectors)
    else:
        assert len(vectors) == len(scalars)

    resultant = [0, 0, 0]
    for i in range(len(vectors)):
        resultant = [resultant[c] + scalars[i] * vectors[i][c] for c in range(3)]
    return resultant


def scalarProjection(v1, v2):
    """
    :return: Scalar projection of v1 onto v2.
    """
    return dotproduct(v1, normalise(v2))


def vectorProjection(v1, v2):
    """
    Calculate vector projection of v1 on v2
    :return: A projection vector.
    """
    s1 = scalarProjection(v1, v2)
    return scaleVector(normalise(v2), s1)


def vectorRejection(v1, v2):
    """
    Calculate vector rejection of v1 on v2
    :return: A rejection vector.
    """
    v1p = vectorProjection(v1, v2)
    return addVectors([v1, v1p], [1.0, -1.0])


def scaleVector(v, s):
    """
    Calculate s * v
    :param v: Vector.
    :param s: Scalar.
    :return:
    """
    return [s * c for c in v]


def parallelVectors(v1, v2):
    """
    :return: True if the vectors are parallel.
    """
    assert (len(v2) == len(v1)), 'Vectors lengths are not the same.'
    TOL = 1.0e-6/2.0
    if magnitude(crossproduct3(v1, v2)) < TOL * (magnitude(v1)+magnitude(v2)):
        return True
    return False


def angleBetweenVectors(v1, v2):
    """
    :return: Angle between vectors v1 and v2 in radians
    """
    return math.acos(dotproduct(normalise(v1), normalise(v2)))


def rotateVectorAroundVector(v, k, a):
    """
    Rotate vector v, by an angle a (right-hand rule) in radians around vector k.
    :return: rotated vector.
    """
    k = normalise(k)
    vperp = addVectors([v, crossproduct3(k, v)], [math.cos(a), math.sin(a)])
    vparal = scaleVector(k, dotproduct(k, v)*(1 - math.cos(a)))
    return addVectors([vperp, vparal])
