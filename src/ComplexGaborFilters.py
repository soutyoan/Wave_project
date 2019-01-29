#!/usr/bin/env python3
"""
Implementation of complex Gabor filters and bands
"""
import numpy as  np
import os
import cv2
from scipy.ndimage.filters import gaussian_filter1d

# Building root path
ROOT_PATH=""
for s in os.path.abspath(__file__).split('/'):
    ROOT_PATH+=s+'/'
    if s=='Wave_project':
        break

IMG_PATH = ROOT_PATH+"/images/"

OUTPUT_PATH = ROOT_PATH+"/output/"

class ComplexGabor():
    """
    complex Gabor filters
    """

    def __init__(self, image_path, nb_scales, nb_theta):
        """
            Initializer class.

            Parameters
            ----------
            image_path      string
                            input image name stored in IMG_PATH

            nb_scales       integer
                            number of wanted scales

            nb_theta        integer
                            number of wished orientations

        """
        self.image = cv2.imread(IMG_PATH+image_path, 0)
        assert self.image is not None
        self.nb_scales = nb_scales
        self.nb_theta = nb_theta
        self.frequencies = np.linspace(1.0, float(nb_scales), num=nb_scales)

    def compute_bands(self):
        """
            Subbands computations. We do not process with downsampling method

            Returns
            -------
            np.ndarray
                    F of shape (nb_scales, nb_theta, self.image.shape) containing
                    all complex subbands
        """
        H, W = self.image.shape[:2]
        N = self.nb_theta
        O = self.nb_scales
        f = self.image
        F = np.zeros((O, N, H, W))

        for i in range(O):
            w_i = self.frequencies[i]
            sig_i = 2*np.pi/w_i
            N_h = N//2
            for k in range(N):
                theta_k = k*np.pi/N
                # computation of J_w_theta_sig (x,y)
                F[i, k, : , :] = self.compute_F(w_i, theta_k, sig_i)
        return F

    def compute_F(self, w_i, theta_k, sig_i):
                """
                    Computation of filter band F_w_theta_sig

                    Parameters
                    ----------
                    w_i         float
                                frequency value
                    theta_k     Float
                                orientation (in rad.)
                    sig_i       float
                                scale value

                    Returns
                    -------
                    np.ndarray
                                Shape : (H, W) of J_w_theta_sig
                """
                f = self.image
                H, W = f.shape[:2]
                J = self.compute_J(w_i, theta_k, sig_i)
                F = np.empty((H, W), dtype=complex)
                cos_w_theta = w_i*np.cos(theta_k)
                sin_w_theta = w_i*np.sin(theta_k)
                for x in range(W):
                    f_cr = np.copy(np.real(J[:, x]).flatten()) * np.array([np.cos(l*sin_w_theta) for l in range(H)])
                    f_sr = np.copy(np.real(J[:, x]).flatten()) * np.array([np.sin(l*sin_w_theta) for l in range(H)])
                    f_ci = np.copy(np.imag(J[:, x]).flatten()) * np.array([np.cos(l*sin_w_theta) for l in range(H)])
                    f_si = np.copy(np.imag(J[:, x]).flatten()) * np.array([np.sin(l*sin_w_theta) for l in range(H)])
                    if (theta_k <= np.pi / 2):
                        fcr_fsi = gaussian_filter1d(f_cr+f_si, sig_i)
                        fsr_fci = gaussian_filter1d(f_sr-f_ci, sig_i)
                    else:
                        fcr_fsi = gaussian_filter1d(f_cr-f_si, sig_i)
                        fsr_fci = gaussian_filter1d(f_sr+f_ci, sig_i)
                    _temp = fcr_fsi + 1j * fsr_fci
                    F[:, x] = _temp
                return F

    def compute_J(self, w_i, theta_k, sig_i):
                """
                    Computation of J_w_theta_sig

                    Parameters
                    ----------
                    w_i         float
                                frequency value
                    theta_k     Float
                                orientation (in rad.)
                    sig_i       float
                                scale value

                    Returns
                    -------
                    np.ndarray
                                Shape : (H, W) of J_w_theta_sig
                """
                if (theta_k > np.pi / 2):
                    return np.conjugate(self.compute_J(w_i, np.pi - theta_k, sig_i))

                f = self.image
                H, W = f.shape[:2]
                J = np.empty((H, W), dtype=complex)
                cos_w_theta = w_i*np.cos(theta_k)
                sin_w_theta = w_i*np.sin(theta_k)
                for y in range(H):
                    f_c = np.copy(f[y, :]) * np.array([np.cos(k*cos_w_theta) for k in range(W)])
                    f_s = np.copy(f[y, :]) * np.array([np.sin(k*cos_w_theta) for k in range(W)])
                    smoothed_fc = gaussian_filter1d(f_c, sig_i)
                    smoothed_fs = gaussian_filter1d(f_s, sig_i)
                    J[y, :] = smoothed_fc + smoothed_fs + 1j * (smoothed_fs - smoothed_fc)
                return J

if (__name__=="__main__"):
    image_path = "small_lena.jpg"
    nb_scales = 5
    nb_theta = 8
    gabor = ComplexGabor(image_path, nb_scales, nb_theta)
    F = gabor.compute_bands()
    for i in range(nb_scales):
        for k in range(nb_theta):
            cv2.imwrite(OUTPUT_PATH+"subband_{}_{}.jpg".format(i, k), np.absolute(F[i, k]))
