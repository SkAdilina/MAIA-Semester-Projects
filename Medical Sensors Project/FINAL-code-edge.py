"""pip install numpy
pip install matplotlib
pip install PIL
pip install opencv-python
pip install scikit-image
pip install argparse
pip install imutils
"""

import os
dirname = os.path.dirname(__file__)

import numpy as np
import cv2
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

from skimage.feature import peak_local_max
from skimage.morphology import watershed
from scipy import ndimage
import argparse
import imutils

for filename in os.listdir(os.path.join(dirname, 'prediction')):

  # Importing Image
  img = cv2.imread((os.path.join(dirname, 'prediction', filename)), 1)
  print("Working on {0}.......".format(filename))
  #serial = 1
  #img = cv2.imread(os.path.join(dirname, 'prediction', 'patient{0}-pred-fig-2.png'.format(serial)),1)
  plt.axis('off')
  plt.imshow(img)

  #######################################################################
  ######################## ROI (Region of Interest) #####################
  #######################################################################

  #img (y1,y2 and x1,x2)
  roi_img = img[135:290, 265:410]
  plt.axis('off')
  plt.imshow(roi_img)
  #plt.show()

  gray = cv2.cvtColor(roi_img, cv2.COLOR_BGR2GRAY) #erosion 
  blur = cv2.GaussianBlur(gray,(7,7),5)

  plt.imshow(cv2.cvtColor(blur ,cv2.COLOR_BGR2RGB)) # load image using cv2....and do processing.
  plt.axis('off')
  #plt.show() # as opencv loads in BGR format by default, we want to show it in RGB
  
  #######################################################################
  ### Otsu
  ### https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_thresholding/py_thresholding.html
  #######################################################################
  
  # gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
  thresh = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
  plt.imshow(cv2.cvtColor(thresh ,cv2.COLOR_BGR2RGB)) # load image using cv2....and do processing.
  plt.axis('off')
  #plt.show() # as opencv loads in BGR format by default, we want to show it in RGB.

  
  #######################################################################
  ### Erosion (image noise) 
  ### https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_morphological_ops/py_morphological_ops.html"""
  #######################################################################

  kernel = np.ones((3,3),np.uint8)
  erosion = cv2.erode(thresh,kernel,iterations = 7)
  #print(img)
  plt.axis('off')
  plt.imshow(erosion)
  #plt.show()

  #######################################################################
  ### Edge detection : https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_canny/py_canny.html
  #######################################################################

  edges = cv2.Canny(erosion,70,300) #(100,200)
  #edges = cv2.Laplacian(erosion,cv2.CV_64F)
  plt.axis('off')
  plt.imshow(edges)
  #plt.show()

  #######################################################################
  ################# Euclidian-distance transform ########################
  #######################################################################

  # compute the exact Euclidean distance from every binary pixel to the nearest zero pixel, then find peaks in this distance map
  D = ndimage.distance_transform_edt(edges)
  localMax = peak_local_max(D, indices=False, min_distance=45, labels=thresh)

  ##############################################################################
  ### Watershed 
  ### https://www.pyimagesearch.com/2015/11/02/watershed-opencv/
  ##############################################################################

  # perform a connected component analysis on the local peaks, using 8-connectivity, then appl#y the Watershed algorithm
  markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0]
  labels = watershed(-D, markers, mask=thresh)
  #print("[INFO] {} unique segments found".format(len(np.unique(labels)) - 1))

  # loop over the unique labels returned by the Watershed algorithm
  for label in np.unique(labels):
    if label == 0 :
      continue
    # otherwise, allocate memory for the label region and draw it on the mask
    mask = np.zeros(gray.shape, dtype="uint8")
    mask[labels == label] = 255
    
    #resized = cv2.resize(mask, img)

    # detect contours in the mask and grab the largest one
    cnts = cv2.findContours(mask, cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)  
    cnts = imutils.grab_contours(cnts)
    c = max(cnts, key=cv2.contourArea)
    

  edges = cv2.Canny(mask, 100, 200) #100,200
  #edges = cv2.Laplacian(mask,cv2.CV_64F)
  plt.axis('off')
  plt.imshow(edges)
  plt.imshow(roi_img, alpha=0.5)
  #plt.savefig(os.path.join(dirname, 'edge', 'patient{0}_e_segments.png'.format(serial)), transparent=True) 
  plt.savefig(os.path.join(dirname, 'all', filename), transparent=True) 
  #plt.show()
  