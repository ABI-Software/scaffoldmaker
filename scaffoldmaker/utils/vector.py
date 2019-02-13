'''
Utility functions for vectors.
Created on Apr 13, 2018

@author: Richard Christie
@edit: Mahyar Osanlouy
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

def diagonal(m):
    '''
    :param m: matrix
    :return: Diagonal elements of m in a list.
    '''
    width, height = len(m[0]), len(m)

    def diag(sx, sy):
        for x, y in zip(range(sx, height), range(sy, width)):
            yield m[x][y]
    for sx in range(height):
        yield list(diag(sx, 0))
    for sy in range(1, width):
        yield list(diag(0, sy))

def identity(n):
    '''
    :param n: Number of rows (and columns) in n x n output.
    :return: The identity matrix
    '''
    return [[0] * i + [1] + [0] * (n - i - 1) for i in range(n)]

def outerProduct(a, b):
    '''

    :param a: First input vector
    :param b: Second input vector
    :return: The outer product of a and b.
    '''
    return b[:, None] * a