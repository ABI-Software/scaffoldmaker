'''
Utility class for representing surfaces on which features can be located.
'''

from __future__ import division
import copy
import math
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector


class TrackSurface:
    '''
    A surface description on which positions can be stored and tracked for
    specified directions and distances for location surface features.
    Currently represented by a lattice of elementsCount1 x elementsCount2
    square elements with bicubic Hermite interpolation but zero cross derivatives.
    '''

    def __init__(self, elementsCount1, elementsCount2, nx, nd1, nd2):
        '''
        Creates a TrackSurface with a lattice of elementsCount1*elementsCount2
        elements, with nodes and elements varying across direction 1 fastest.
        :param elementsCount1, elementsCount2: Number of elements in directions 1 and 2.
        :param nx, nd1, nd2: List of (elementsCount2 + 1)*(elementsCount1 + 1) node
        coordinates, derivatives 1 and derivatives 2. Cross derivatives are zero.
        All coordinates are 3 component.
        '''
        self.elementsCount1 = elementsCount1
        self.elementsCount2 = elementsCount2
        self.nx = nx
        self.nd1 = nd1
        self.nd2 = nd2

    def createPositionProportion(self, proportion1, proportion2):
        '''
        Return position on surface for proportions across directions 1 and 2.
        :param proportion1, proportion2: Proportions across directions 1 and 2,
        each varying from 0.0 to 1.0, with elements equal sized.
        :return: TrackSurfacePosition
        '''
        assert (proportion1 >= 0.0) and (proportion1 <= 1.0) and (proportion2 >= 0.0) and (proportion2 <= 1.0), 'Proportion out of range'
        pe1 = proportion1*self.elementsCount1
        if pe1 < self.elementsCount1:
            e1 = int(pe1)
            xi1 = pe1 - e1
        else:
            e1 = self.elementsCount1
            xi1 = 1.0
        pe2 = proportion2*self.elementsCount2
        if pe2 < self.elementsCount2:
            e2 = int(pe2)
            xi2 = pe2 - e2
        else:
            e2 = self.elementsCount2
            xi2 = 1.0
        return TrackSurfacePosition(e1, e2, xi1, xi2)

    def evaluateCoordinates(self, position, derivatives=False):
        '''
        Evaluate coordinates on surface at position, and optionally
        derivative w.r.t. xi1 and xi2.
        :param position: A valid TrackSurfacePosition.
        :return: If derivatives is False: coordinates [ x, y, z];
        If derivatives is True: coordinates, derivative1, derivative2.
        '''
        n1 = position.e2*(self.elementsCount1 + 1) + position.e1
        n2 = n1 + 1
        n3 = n1 + self.elementsCount1 + 1
        n4 = n3 + 1
        ni = [ n1, n2, n3, n4 ]
        f1x1, f1d1, f1x2, f1d2 = interp.getCubicHermiteBasis(position.xi1)
        f2x1, f2d1, f2x2, f2d2 = interp.getCubicHermiteBasis(position.xi2)
        fx  = [ f1x1*f2x1, f1x2*f2x1, f1x1*f2x2, f1x2*f2x2 ]
        fd1 = [ f1d1*f2x1, f1d2*f2x1, f1d1*f2x2, f1d2*f2x2 ]
        fd2 = [ f1x1*f2d1, f1x2*f2d1, f1x1*f2d2, f1x2*f2d2 ]
        nx = self.nx
        nd1 = self.nd1
        nd2 = self.nd2
        coordinates = []
        for c in range(3):
            x = 0.0
            for ln in range(4):
                gn = ni[ln]
                x += fx[ln]*nx[gn][c] + fd1[ln]*nd1[gn][c] + fd2[ln]*nd2[gn][c]
            coordinates.append(x)
        if not derivatives:
            return coordinates
        df1x1, df1d1, df1x2, df1d2 = interp.getCubicHermiteBasisDerivatives(position.xi1)
        d1fx  = [ df1x1*f2x1, df1x2*f2x1, df1x1*f2x2, df1x2*f2x2 ]
        d1fd1 = [ df1d1*f2x1, df1d2*f2x1, df1d1*f2x2, df1d2*f2x2 ]
        d1fd2 = [ df1x1*f2d1, df1x2*f2d1, df1x1*f2d2, df1x2*f2d2 ]
        df2x1, df2d1, df2x2, df2d2 = interp.getCubicHermiteBasisDerivatives(position.xi2)
        d2fx  = [ f1x1*df2x1, f1x2*df2x1, f1x1*df2x2, f1x2*df2x2 ]
        d2fd1 = [ f1d1*df2x1, f1d2*df2x1, f1d1*df2x2, f1d2*df2x2 ]
        d2fd2 = [ f1x1*df2d1, f1x2*df2d1, f1x1*df2d2, f1x2*df2d2 ]
        derivative1 = []
        derivative2 = []
        for c in range(3):
            d1 = 0.0
            d2 = 0.0
            for ln in range(4):
                gn = ni[ln]
                d1 += d1fx[ln]*nx[gn][c] + d1fd1[ln]*nd1[gn][c] + d1fd2[ln]*nd2[gn][c]
                d2 += d2fx[ln]*nx[gn][c] + d2fd1[ln]*nd1[gn][c] + d2fd2[ln]*nd2[gn][c]
            derivative1.append(d1)
            derivative2.append(d2)
        return coordinates, derivative1, derivative2

    def trackVector(self, startPosition, direction, trackDistance):
        '''
        Track from startPosition the given distance in the vector direction.
        Approximate, uses improved Euler method (mean of original & final gradient).
        :param startPosition: TrackSurfacePosition
        :param direction: 3-D vector (x, y, z) to track along. Projected onto surface.
        :param trackDistance: Distance to track along. Can be negative.
        :return: Final TrackSurfacePosition
        '''
        useDirection = direction
        useTrackDistance = trackDistance
        if trackDistance < 0.0:
            useDirection = [ -d for d in direction ]
            useTrackDistance = -trackDistance
        position = copy.deepcopy(startPosition)
        distance = 0.0
        distanceLimit = 0.9999*useTrackDistance
        max_mag_dxi = 0.02  # target/maximum magnitude of xi increment

        while distance < useTrackDistance:
            xi1 = position.xi1
            xi2 = position.xi2
            ax, ad1, ad2 = self.evaluateCoordinates(position, derivatives = True)
            adelta_xi1, adelta_xi2 = calculate_surface_delta_xi(ad1, ad2, useDirection)
            scale = max_mag_dxi/math.sqrt(adelta_xi1*adelta_xi1 + adelta_xi2*adelta_xi2)
            adxi1 = dxi1 = scale*adelta_xi1
            adxi2 = dxi2 = scale*adelta_xi2
            for i in range(1):  # can try more in future
                # can go slightly outside element to get predictor/correct position
                position.xi1 = xi1 + dxi1
                position.xi2 = xi2 + dxi2
                bx, bd1, bd2 = self.evaluateCoordinates(position, derivatives = True)
                bdelta_xi1, bdelta_xi2 = calculate_surface_delta_xi(bd1, bd2, useDirection)
                # use mean of start and end derivatives
                delta_xi1 = 0.5*(adelta_xi1 + bdelta_xi1)
                delta_xi2 = 0.5*(adelta_xi2 + bdelta_xi2)
                scale = max_mag_dxi/math.sqrt(delta_xi1*delta_xi1 + delta_xi2*delta_xi2)
                dxi1 = scale*delta_xi1
                dxi2 = scale*delta_xi2
            bxi1, bxi2, proportion, faceNumber = increment_xi_on_square(xi1, xi2, dxi1, dxi2)
            position.xi1 = bxi1
            position.xi2 = bxi2
            #print(distance, '-->', position)
            bx, bd1, bd2 = self.evaluateCoordinates(position, derivatives = True)
            bdelta_xi1, bdelta_xi2 = calculate_surface_delta_xi(bd1, bd2, useDirection)
            scale = max_mag_dxi/math.sqrt(bdelta_xi1*bdelta_xi1 + bdelta_xi2*bdelta_xi2)
            bdxi1 = scale*bdelta_xi1
            bdxi2 = scale*bdelta_xi2
            ad = [ proportion*(adxi1*ad1[c] + adxi2*ad2[c]) for c in range(3) ]
            bd = [ proportion*(bdxi1*bd1[c] + bdxi2*bd2[c]) for c in range(3) ]
            # GRC check:
            #arcLength = interp.computeCubicHermiteArcLength(ax, ad, bx, bd, rescaleDerivatives = True)
            arcLength = interp.getCubicHermiteArcLength(ax, ad, bx, bd)
            #print('scales', vector.magnitude([ (bx[c] - ax[c]) for c in range(3) ]), vector.magnitude(ad), vector.magnitude(bd), 'arc length', arcLength)
            if (distance + arcLength) >= distanceLimit:
                # limit to useTrackDistance, approximately, and finish
                r = proportion*(useTrackDistance - distance)/arcLength
                position.xi1 = xi1 + r*dxi1
                position.xi2 = xi2 + r*dxi2
                #print(distance, '-->', position, '(final)')
                break
            if arcLength == 0.0:
                print('TrackSurface.trackVector. No increment at', position, 'final distance', distance, 'of', useTrackDistance)
                break
            distance += arcLength
            if faceNumber:
                onBoundary = False
                if faceNumber == 1:  # xi1 == 0.0
                    if position.e1 > 0:
                        position.e1 -= 1
                        position.xi1 = 1.0
                    else:
                        onBoundary = True
                if faceNumber == 2:  # xi1 == 1.0
                    if position.e1 < (self.elementsCount1 - 2):
                        position.e1 += 1
                        position.xi1 = 0.0
                    else:
                        onBoundary = True
                if faceNumber == 3:  # xi2 == 0.0
                    if position.e2 > 0:
                        position.e2 -= 1
                        position.xi2 = 1.0
                    else:
                        onBoundary = True
                if faceNumber == 4:  # xi2 == 1.0
                    if position.e2 < (self.elementsCount2 - 2):
                        position.e2 += 1
                        position.xi2 = 0.0
                    else:
                        onBoundary = True
                if onBoundary:
                    write('TrackSurface.trackVector:  End on boundary at', position)
                    break
                #print('  cross face', faceNumber, 'new position', position)
        return position


class TrackSurfacePosition:
    '''
    Position on TrackSurface. Create with createPosition~ method of TrackSurface.
    '''

    def __init__(self, e1, e2, xi1, xi2):
        '''
        :param e1, e2: Element index in directions 1 and 2, starting at 0.
        :param xi1, xi2: 2-D element chart coordinates 0.0 <= xiN <= 1.0
        '''
        self.e1 = e1
        self.e2 = e2
        self.xi1 = xi1
        self.xi2 = xi1

    def __str__(self):
        return 'element (' + str(self.e1) + ',' + str(self.e2) + ') xi (' + str(self.xi1) + ',' + str(self.xi2) + ')'

def calculate_surface_delta_xi(d1, d2, direction):
    '''
    Calculate dxi1, dxi2 in 3-D vector direction.
    :param d1, d2: Derivatives of coordinate w.r.t. xi1, xi2.
    :param direction: 3-D vector.
    :return: delta_xi1, delta_xi2
    '''
    # overdetermined system (3-D vector, 2-D surface), solve least squares
    # A transpose A x = A transpose b
    a = [ [ 0.0, 0.0 ], [ 0.0, 0.0 ]]
    b = [ 0.0, 0.0 ]
    dx_dxi = [ d1, d2 ]
    for i in range(2):
        for j in range(2):
            for k in range(3):
                a[i][j] += dx_dxi[i][k]*dx_dxi[j][k]
        for k in range(3):
            b[i] += dx_dxi[i][k]*direction[k];
    # 2x2 matrix inverse
    deta = a[0][0]*a[1][1] - a[0][1]*a[1][0]
    inva = [ [ a[1][1]/deta, -a[0][1]/deta ], [ -a[1][0]/deta, a[0][0]/deta ] ]
    delta_xi1 = inva[0][0]*b[0] + inva[0][1]*b[1]
    delta_xi2 = inva[1][0]*b[0] + inva[1][1]*b[1]
    return delta_xi1, delta_xi2

def calculate_surface_axes(d1, d2, direction):
    '''
    :return: Vectors ax1, ax2, ax3; ax1 in-plane in 3-D vector direction,
    ax2 in-plane normal to a and ax3 normal to the surface plane.
    Vectors all have unit magnitude.
    '''
    ax3 = vector.normalise(vector.crossproduct3(d1, d2))
    delta_xi1, delta_xi2 = calculate_surface_delta_xi(d1, d2, direction)
    ax1 = vector.normalise([ delta_xi1*d1[c] + delta_xi2*d2[c] for c in range(3) ])
    ax2 = vector.normalise(vector.crossproduct3(ax3, ax1))
    return ax1, ax2, ax3

def increment_xi_on_square(xi1, xi2, dxi1, dxi2):
    '''
    Increment xi1, xi2 by dxi1, dxi2 limited to square element bounds on [0,1].
    :return: New xi1, xi2, proportion, face number 1-4 or None if within boundary.
    Proportion is 1.0 if increment is within element, otherwise equal to proportion
    of dxi1, dxi2 incremented to hit the boundary. If so limited to the boundary,
    the face number is returned as an integer:
    1 is xi1=0.0, 2 is xi1=1.0, 3 is xi2=0.0, 4 is xi2=1.0
    '''
    nxi1 = xi1 + dxi1
    nxi2 = xi2 + dxi2
    proportion = 1.0
    faceNumber = None
    if (nxi1 < 0.0) or (nxi1 > 1.0) or (nxi2 < 0.0) or (nxi2 > 1.0):
        # come back in direction of dxi1, dxi2 to first boundary
        if (nxi1 < 0.0) and (dxi1 < 0.0):
            thisProportion = -xi1/dxi1
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 1
                nxi1 = 0.0
                nix2 = xi2 + proportion*dxi2
        if (nxi1 > 1.0) and (dxi1 > 0.0):
            thisProportion = (1.0 - xi1)/dxi1
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 2
                nxi1 = 1.0
                nix2 = xi2 + proportion*dxi2
        if (nxi2 < 0.0) and (dxi2 < 0.0):
            thisProportion = -xi2/dxi2
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 3
                nxi1 = xi1 + proportion*dxi1
                nxi2 = 0.0
        if (nxi2 > 1.0) and (dxi2 > 0.0):
            thisProportion = (1.0 - xi2)/dxi2
            if thisProportion < proportion:
                proportion = thisProportion
                faceNumber = 4
                nxi1 = xi1 + proportion*dxi1
                nxi2 = 1.0
        #print('  increment to face', faceNumber, 'xi (' + str(xi1) + ',' + str(xi2) + ')', 'nxi (' + str(nxi1) + ',' + str(nxi2) + ')')
    return nxi1, nxi2, proportion, faceNumber
