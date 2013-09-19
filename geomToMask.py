#!/usr/bin/env python

#Written by Duane Loh.
#First distributed on 6th May 2013.
#
# Usage:
#	./geomToMask.py

import os
import numpy as N
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as PL
import time
from myModules import geo 

#Import geometry file -- change the path name if you'd like
relPath = os.path.abspath(".")
#geometryFile = "/reg/d/psdm/cxi/cxi73013/scratch/duaneloh/config/geometry/CxiDs1-geometry.h5"		#Config folder by default is set one directory up
geometryFile = "/Users/duaneloh/Dropbox/CurrentProjects/L73013/config/geometry/CxiDs1-geometry.h5"		
writeDir = "/Users/duaneloh/Dropbox/CurrentProjects/L73013/config/mask/"

#Use h5py module to open geometry file, read (x,y,z) coordinates into [x,y,z] arrays, then close the file
# z-coordinate array isn't used but shown here for reference
f = h5py.File(geometryFile, "r")
[x, y, z] = [N.array(f['x']), N.array(f['y']), N.array(f['z'])]  #Python syntax for performing three assignments sequentially
f.close()

#Shift the array such that average pixel position is at (0,0)
x -= x.mean()
y -= y.mean()

#Compute radial distance of pixels from (0,0) position
rad = N.sqrt(x*x + y*y)

#Specify default parameters about the detector
wavelengthInAngs = 1.693
detectorDist = 137.E-3
pixSize = 0.0001099


##########################################################################
# Start of lines to modify if you'd like a different mask
##########################################################################
#Mask to keep only pixels between spatial frequencies q1-q2
#[q1,q2] = [0., 1.]		#Units in inverse Angstroms
[q1,q2] = [1.15, 1.43]
[q3, q4] = [1.73, 1.86]
[q5, q6] = [2.1, 2.5]
ang1 = 2.*N.arcsin( (wavelengthInAngs*q1) / (4.*N.pi) )
ang2 = 2.*N.arcsin( (wavelengthInAngs*q2) / (4.*N.pi) )
ang3 = 2.*N.arcsin( (wavelengthInAngs*q3) / (4.*N.pi) )
ang4 = 2.*N.arcsin( (wavelengthInAngs*q4) / (4.*N.pi) )
ang5 = 2.*N.arcsin( (wavelengthInAngs*q5) / (4.*N.pi) )
ang6 = 2.*N.arcsin( (wavelengthInAngs*q6) / (4.*N.pi) )

lmin = detectorDist * N.tan(ang1)
lmax = detectorDist * N.tan(ang2)
lmin3 = detectorDist * N.tan(ang3)
lmax4 = detectorDist * N.tan(ang4)
lmin5 = detectorDist * N.tan(ang5)
lmax6 = detectorDist * N.tan(ang6)


#Mask out pixels to the left and right of the center by 300 pixels
#You might have to change these number from run to run
mask = ((rad < lmax)*(rad > lmin) + (rad < lmax4)*(rad > lmin3) + (rad < lmax6)*(rad > lmin5))*((y > 300*pixSize)+(y < -300*pixSize))

mask2 = N.ones((370,388)).astype('int')
#mask2[:,:240] = 0
mask2[:,184:234] = -1

##########################################################################
# End of lines to modify if you'd like a different mask
##########################################################################

#Shift the detector pixels back to positive memory addresses (added a 3 pixel padding on the left and bottom)
x -= x.min() - 3.*pixSize
y -= y.min() - 3.*pixSize

#Recast the coordinates as float types for Cythonic module
xf = x.astype('float')
yf = y.astype('float')
maskf = mask.astype('float')

f = h5py.File('cxi73013-r0039.cxi','r')
d = (f[u'entry_1/instrument_1/detector_1/data'][0]).astype('float')
f.close()

#Timed interpolation of intensities
startTime = time.clock()
arr = geo.interpolateGeometry(xf, yf, maskf, pixSize)
#arr = geo.interpolateGeometry(xf, yf, d, pixSize)
endTime = time.clock()
print "Took %f seconds"%(endTime - startTime)

#Create directory if it's present
#if not os.path.exists(writeDir):
#	os.makedirs(writeDir)
#
##Write viewing mask to file
#outFile = h5py.File(writeDir + "/L73013_viewingMask_FrontDetector.h5", "w")
#g1 = outFile.create_group('data')
#g1.create_dataset('data', data=mask)
#outFile.close()
#
##TODO: Change back detector mask to include common mode pixels.
#outFile = h5py.File(writeDir + "/L73013_viewingMask_BackDetector.h5", "w")
#g1 = outFile.create_group('data')
#g1.create_dataset('data', data=mask2)
#outFile.close()
