#!/usr/bin/env python3
"""
Local Phase Coherence Implementation
"""

import os, sys, cv2
import numpy as np
from SteerablePyramid import SteerablePyramid
import argparse
import progressbar

DEBUG_PRINT = 0  # for debug purposes

# Building root path
ROOT_PATH=""
for s in os.path.abspath(__file__).split('/'):
    ROOT_PATH+=s+'/'
    if s=='Wave_project':
        break


def get_w(N, d=1.0, s1=1.0, mode='lin'):
    """
        Computation of w vector used for the LPC strength computation.
        w is the solution of the linear system :

                        Aw = b
        where A (N+1)x(N+1) and b (N+1) both depend on real coefficients determined by d.
        Those real coefficients, stored in a vector, can be either computed with a linear
        mode such as s = [s1, s1 + d, s1 + 2d, ..., s1 + (N-1)d] or a log mode s.t
                     s = [d^0, d, d^2, ..., d^(N-1)]

        Parameters
        ----------
        N           Integer
                    number of scales of the steerable pyramid
        d           float, optional
                    step value for s vector computation
        s1          float, optional
                    choice of strength of the finest scale, default = 1
        mode        string, {'lin', 'log'}
                    mode choice for scales strength, default = 'lin'

        Returns
        -------
        numpy.ndarray
                    w vector of length N
    """
    # vector containing s2, ..., sN values for system resolution
    _s = np.zeros(N-1)
    if (mode=='lin'):
        _s = np.array([1+i*d for i in range(1,N)])
    elif (mode=='log'):
        _s = np.array([d**i for i in range(1, N)])
    if (DEBUG_PRINT):
        print("s vector :\n", _s)
    # A matrix construction
    A = np.zeros((N+1, N+1))
    A[:-2,:-2] = np.identity(N-1)
    A[-2,:-2] = 1
    A[-1,:-2] = 1/_s
    A[:-2,-2] = 0.5
    A[:-2,-1] = 0.5*1/_s.T
    if (DEBUG_PRINT):
        print("A matrix \n", A)
    # b vector construction
    b = -1 * np.ones(N+1)
    b[:-2] = 0
    if (DEBUG_PRINT):
        print("b vector \n", b)
    res = np.ones(N)
    res[1:] = np.linalg.solve(A, b)[:-2]
    return res


def compute_LPC_strength(steer, j, k, d=1, s1=1, mode='lin'):
    """
        Computation of Local Phase Coherence strength

        Parameters
        ----------
        steer       SteerablePyramid
                    input steerable Pyramid
        j           integer
                    chosen orientation
        k           integer
                    spatial location

        Returns
        -------
        float
                    S_{LPC}^{j,k}
    """
    _w = get_w(steer.N, d=d, s1=s1, mode=mode)
    _phi = []
    ind = k
    for i in range(steer.N):
        #print(ind, steer.BND[i][j]['s'].size)
        _cijk = steer.BND[i][j]['s'].flatten()[ind]
        _phi.append(np.angle(_cijk))
        ind >>= 2
    _phi = np.array(_phi)
    return np.cos(np.dot(_w, _phi))



def compute_spatial_LPC(steer, k, C):
    """
        Computation of Local Phsae Coherence coefficient.

        Parameters
        ----------
        steer       SteerablePyramid
                    input steerable pyramid
        k           integer
                    spatial location
        C           float
                    constant for stabilisation purposes

        Returns
        -------
        float
                    S_{LPC}^{k}
    """
    c1jk = np.absolute(np.array([steer.BND[0][j]['s'].flatten()[k] for j in range(steer.K)]))
    phi = np.array([compute_LPC_strength(steer,j,k) for j in range(steer.K)])
    return np.sum(c1jk * phi) / (np.sum(c1jk) + C)



def compute_LPC_map(steer, C):
    """
        Computation of LPC map.

        Parameters
        ----------
        steer       SteerablePyramid
                    input steerable pyramid
        C           float > 0
                    constant for stabilisation purposes

        Returns
        -------
        numpy.ndarray
                    map of LPC coefficient
    """
    res = np.zeros(steer.IMAGE_ARRAY.size)
    res = res.flatten()
    bar = progressbar.ProgressBar(max_value=res.shape[0])
    for k in range(res.shape[0]):
        res[k] = compute_spatial_LPC(steer, k, C)
        bar.update(k)
    return res.reshape(steer.IMAGE_ARRAY.shape[:2])


def compute_LPC_index(steer, C, beta):
    """
        Computation of LPC index.

        Parameters
        ----------
        steer           SteerablePyramid
                        input steerable pyramid
        C               float > 0
                        constant for stabilisation purposes
        beta            float > 0
                        constant of weighting LPC coefficients

        Returns
        -------
        float in [0, 1]
                        LPC-based sharpness index of steer.IMAGE_ARRAY
    """
    # TODO

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description="Wavelets project 2018-2019 : Sharpness index based on Local Phase Coherence")
    parser.add_argument("--input_file", '-i', help='Input File Name (file must be in images/ folder)')
    parser.add_argument('--depth', '-N', default=3, type=int, help='Depth of Pyramid. Integer')
    parser.add_argument('--orientation', '-M', default=8, type=int, help='Number of orientations of pyramid. Integer')
    parser.add_argument('--constant', '-C', default=2, type=float, help='Constant for LPC coefficients computation. Float')
    parser.add_argument('--beta', '-B', default=0.0001, type=float, help='Constant for LPC coefficients weighting. Float')
    parser.add_argument('--sampled', '-K', default=1.0, type=float, help='Rate of extracted values from LPC map. Float')
    parser.add_argument('--verbose', '-v', default=1, type=int, help='verbose')
    args = parser.parse_args()

    # Grayscale mode image reading
    image = cv2.imread(ROOT_PATH+"/images/"+args.input_file, 0)
    assert image is not None
    yres, xres = image.shape[:2]

    # Steerable pyramid building
    image_name = args.input_file.split('.')[0]
    steer = SteerablePyramid(image, xres, yres, args.depth, args.orientation, image_name, ROOT_PATH+"output", args.verbose)
    steer.create_pyramids()

    # Computation of LPC map
    LPC_map = compute_LPC_map(steer, args.constant)

    print(LPC_map)

    cv2.imwrite(ROOT_PATH+"output/{}_LPC_map.png".format(image_name), np.absolute(250*LPC_map))

    # Computation of LPC score
