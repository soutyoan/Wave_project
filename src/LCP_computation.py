#!/usr/bin/env python3
"""
Local Phase Coherence Implementation
"""

import os, sys, cv2
import numpy as np
from SteerablePyramid import SteerablePyramid
import argparse

DEBUG_PRINT = 1  # for debug purposes

# Building root path
ROOT_PATH=""
for s in os.path.abspath(__file__).split('/'):
    ROOT_PATH+=s+'/'
    if s=='Wave_project':
        break



def compute_LPC_strength(steer, j, k):
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
    # TODO

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
    # TODO
    return 0.0

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
    # TODO
    

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
    parser.add_argument('--depth', '-N', default=3, type=int, help='Depth of Pyramid. Integer in {3, 4, 5}')
    parser.add_argument('--orientation', '-M', default=8, type=int, help='Number of orientations of pyramid. Integer')
    parser.add_argument('--constant', '-C', default=2, type=float, help='Constant for LPC coefficients computation. Float')
    parser.add_argument('--beta', '-B', default=0.0001, type=float, help='Constant for LPC coefficients weighting. Float')
    parser.add_argument('--sampled', '-K', default=1.0, type=float, help='Rate of extracted values from LPC map. Float')
    parser.add_argument('--verbose', '-v', default=1, type=int, help='verbose')
    args = parser.parse_args()

    print(args)

    # Grayscale mode image reading
    image = cv2.imread(ROOT_PATH+"/images/"+args.input_file, 0)
    assert image is not None
    yres, xres = image.shape[:2]

    # Steerable pyramid building
    image_name = args.input_file.split('.')[0]
    steer = SteerablePyramid(image, xres, yres, args.depth, args.orientation, image_name, ROOT_PATH+"/output", args.verbose)

    # Computation of LPC map
    LPC_map = compute_LPC_map(steer, C)

    if steer.verbose:
        cv2.imwrite(ROOT_PATH+"output/{}_LPC_map.png".format(image_name), np.absolute(LPC_map))

    # Computation of
