#!/usr/bin/env python3
"""
Local Phase Coherence Implementation
"""

import os, sys, cv2
import numpy as np
from ComplexGaborFilters import ComplexGabor
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


def compute_LPC_strength(gabor, j, k, d=1, s1=1, mode='lin'):
    """
        Computation of Local Phase Coherence strength

        Parameters
        ----------
        gabor       ComplexGabor
                    complex Gabor decomposition
        j           integer
                    chosen orientation
        k           integer
                    spatial location

        Returns
        -------
        float
                    S_{LPC}^{j,k}
    """
    _w = get_w(gabor.nb_scales, d=d, s1=s1, mode=mode)
    _phi = []
    for i in range(gabor.nb_scales):
        #print(ind, gabor.BND[i][j]['s'].size)
        _cijk = gabor.bands[i][j].flatten()[k]
        _phi.append(np.angle(_cijk))
    _phi = np.array(_phi)
    return np.cos(np.dot(_w, _phi))



def compute_spatial_LPC(gabor, k, C):
    """
        Computation of Local Phsae Coherence coefficient.

        Parameters
        ----------
        gabor       ComplexGabor
                    input complex Gabor filtering
        k           integer
                    spatial location
        C           float
                    constant for stabilisation purposes

        Returns
        -------
        float
                    S_{LPC}^{k}
    """
    c1jk = np.absolute(np.array([gabor.bands[0][j].flatten()[k] for j in range(gabor.nb_theta)]))
    phi = np.array([compute_LPC_strength(gabor,j,k) for j in range(gabor.nb_theta)])
    return np.sum(c1jk * phi) / (np.sum(c1jk) + C)



def compute_LPC_map(gabor, C):
    """
        Computation of LPC map.

        Parameters
        ----------
        gabor       ComplexGabor
                    input complex Gabor filtering
        C           float > 0
                    constant for stabilisation purposes

        Returns
        -------
        numpy.ndarray
                    map of LPC coefficient
    """
    res = np.zeros(gabor.image.shape[:2])
    res = res.flatten()
    print("\n===== LCP map computation =====")
    bar = progressbar.ProgressBar(max_value=res.shape[0])
    for k in range(res.shape[0]):
        res[k] = compute_spatial_LPC(gabor, k, C)
        bar.update(k+1)
    return res.reshape(gabor.image.shape[:2])


def compute_LPC_index(gabor, C, beta, K=1.0):
    """
        Computation of LPC index.

        Parameters
        ----------
        gabor           ComplexGabor
                        input steerable pyramid
        C               float > 0
                        constant for stabilisation purposes
        beta            float > 0
                        constant of weighting LPC coefficients
        K               float in [0, 1.0]n optional
                        rate of chosen coefficients from LCP map, default = 1.0
                        (i.e all the coefficients are chosen by default)

        Returns
        -------
        float in [0, 1]
                        LPC-based sharpness index of gabor.IMAGE_ARRAY
    """
    print("\n===== LCP index computation =====")
    LPC_map = compute_LPC_map(gabor, C)
    if (gabor.verbose):
        cv2.imwrite(ROOT_PATH+"output/{}_LPCmap_N_{}_M_{}_C_{}_B_{}.png".format(gabor.image_name, gabor.nb_scales, gabor.nb_theta, C, beta), np.absolute(255*LPC_map))
    nb_coeffs = int(K * LPC_map.size)
    sorted_LPC = np.flip(np.sort(LPC_map.flatten())[:nb_coeffs], 0)
    weights = np.exp(np.array([-(k-1)/(nb_coeffs-1)/beta for k in range(nb_coeffs)]))
    return np.dot(sorted_LPC, weights) / np.sum(weights)

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description="Wavelets Ensimag project 2018-2019 : Sharpness index based on Local Phase Coherence")
    parser.add_argument("--input_file", '-i', help='Input File Name (file must be in images/ folder)')
    parser.add_argument('--depth', '-N', default=3, type=int, help='Depth of Pyramid. Integer')
    parser.add_argument('--orientation', '-M', default=8, type=int, help='Number of orientations of pyramid. Integer')
    parser.add_argument('--constant', '-C', default=2, type=float, help='Constant for LPC coefficients computation. Float')
    parser.add_argument('--beta', '-B', default=0.0001, type=float, help='Constant for LPC coefficients weighting. Float')
    parser.add_argument('--sampled', '-K', default=1.0, type=float, help='Rate of extracted values from LPC map. Float')
    parser.add_argument('--verbose', '-v', default=1, type=int, help='verbose')
    args = parser.parse_args()

    # Comple Gabor filtering
    image_path = args.input_file
    nb_scales = args.depth
    nb_theta = args.orientation
    gabor = ComplexGabor(image_path, nb_scales, nb_theta, verbose=args.verbose)

    # Computation of LPC index
    S = compute_LPC_index(gabor, args.constant, args.beta)
    print("\n===== Sharpness Value: ", S)
