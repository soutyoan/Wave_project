#!/usr/bin/env python3
"""
Local Phase Coherence Parameters testing with images/small_lena.jpg
"""

import os, sys, cv2
import numpy as np
from ../src/LPC_computation import *
from ../src/ComplexGaborFilters import *
from ../src/detect_Blur import *
import argparse
import progressbar
import matplotlib.pyplot as plt

DEBUG_PRINT = 0  # for debug purposes

# Building root path
ROOT_PATH=""
for s in os.path.abspath(__file__).split('/'):
    ROOT_PATH+=s+'/'
    if s=='Wave_project':
        break

OUTPUT_PATH = ROOT_PATH+"/output/"


def test_N():
    """
        Graph displaying of evolution of LPC index regarding number of scales N
    """
    n_values = np.array([n for n in range(1,11)])
    s_values = np.zeros(10)
    for _ind, n in enumerate(n_values):
        gabor = ComplexGabor("small_lena.jpg", n, 8)
        s_values[_ind] = compute_LPC_index(gabor, 2, 0.0001)
    plt.figure("Test parameter N")
    plt.xlabel("Number of scales N")
    plt.ylabel("S_LPC value")
    plt.plot(n_values, s_values, 'r+', ms=4)
    plt.savefig(OUTPUT_PATH+"LPC_test_N.png")

def test_M():
    """
        Graph displaying of evolution of LPC index regarding number of orientations M
    """
    m_values = np.array([n for n in range(1,11)])
    s_values = np.zeros(10)
    for _ind, m in enumerate(m_values):
        gabor = ComplexGabor("small_lena.jpg", 3, m)
        s_values[_ind] = compute_LPC_index(gabor, 2, 0.0001)
    plt.figure("Test parameter M")
    plt.xlabel("Number of orientations M")
    plt.ylabel("S_LPC value")
    plt.plot(m_values, s_values, 'r+', ms=4)
    plt.savefig(OUTPUT_PATH+"LPC_test_M.png")

def test_B():
    """
        Graph displaying of evolution of LPC index regarding beta parameter
    """
    b_values = np.array([10**-i for i in range(5,-1,-1)])
    s_values = np.zeros(b_values.shape)
    gabor = ComplexGabor("small_lena.jpg", 3, 8)
    for _ind, beta in enumerate(b_values):
        s_values[_ind] = compute_LPC_index(gabor, 2, beta)
    plt.figure("Test parameter B")
    plt.xlabel("Beta value")
    plt.ylabel("S_LPC value")
    plt.plot(b_values, s_values, 'r+', ms=4)
    plt.savefig(OUTPUT_PATH+"LPC_test_B.png")

def test_blur():
    """
        Graph displaying evolution of LPC index regarding blurring effect
    """
    sig_values = np.array([sig for sig in range(11)])
    b_values = np.array([10**-i for i in range(5,-1,-1)])
    s_values = np.zeros((b_values.shape, sig_values.shape))
    gabor = ComplexGabor("small_lena.jpg", 3, 8)
    image_original = np.copy(gabor.image)
    plt.figure("Test blur")
    for _indb, beta in enumerate(b_values):
        for _ind, sig in enumerate(sig_values):
            if (sig!=0):
                gabor.image = cv2.blur(image_original, (sig, sig))
                gabor.compute_bands()
            s_values[_indb, _ind] = compute_LPC_index(gabor, 2, beta)
        plt.plot(sig_values, s_values[_indb], "+", ms=4, label="beta = {}".format(beta))
    plt.legend()
    plt.savefig(OUTPUT_PATH+"LPC_test_blur.png")

if (__name__=="__main__"):
    test_N()
    test_M()
    test_B()
    test_blur()
