from SteerablePyramid import *
import cv2
import numpy as np

image = cv2.imread("../images/lena.jpg", 0)

yres, xres = image.shape[:2]

steer = SteerablePyramid(image, xres, yres, 3, 8, "lena", "../output", False)
cv2.imwrite("../output/low_filt_py.png", steer.H0_FILT)
