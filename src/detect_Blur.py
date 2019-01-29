#!/usr/bin/env python3
"""
Blurring appearence detection using LPC index
"""

import os, sys, cv2
import numpy as np
from ComplexGaborFilters import ComplexGabor
from LPC_computation import *
import argparse
import progressbar

DEBUG_PRINT = 0  # for debug purposes

# Building root path
ROOT_PATH=""
for s in os.path.abspath(__file__).split('/'):
    ROOT_PATH+=s+'/'
    if s=='Wave_project':
        break

def compare_LPC(image_name, sig, beta=0.01):
    """
        Displaying the LPC index of input image and blurred image with Gaussian blur.

        Parameters
        ----------
        image_name          string
                            input image name

        sig                 int
                            input sigma of Gaussian blur

    """
    gabor = ComplexGabor(image_name, 3, 8)
    s1 = compute_LPC_index(gabor, 2, beta)
    gabor.image = cv2.blur(gabor.image, (sig, sig))
    gabor.compute_bands()
    s2 = compute_LPC_index(gabor, 2, beta)
    return s1, s2

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description="Wavelets Ensimag project 2018-2019 : blur detection with LPC sharpness index")
    parser.add_argument("--input_file", '-i', help='Input File Name (file must be in images/ folder)')
    parser.add_argument("--sigma", '-sig', type=int, help='Value of sigma of Gaussian blur')
    parser.add_argument('--beta', '-B', default=0.01, type=float, help='Constant for LPC coefficients weighting. Float')
    parser.add_argument('--verbose', '-v', default=1, type=int, help='verbose')
    args = parser.parse_args()

    image_name = args.input_file
    sig = args.sigma
    beta = args.beta
    s1, s2 = compare_LPC(image_name, sig, beta=beta)

    print("\n===== Sharpness index of original and blurred images : ", s1, s2)
