#!/usr/bin/env python

#Written by Duane Loh.
#First distributed on 6th May 2013.
#
# Usage:
#    ./viewAverages.py  --help

#Import common modules
import os
import sys
import string
import re
from optparse import OptionParser
import glob as G
import datetime

#Import extra modules
import numpy as N
import h5py
import matplotlib as M
import matplotlib.pyplot as P

######################################################################################################
# Read and parse command line input flags
######################################################################################################
parser = OptionParser()
parser.add_option("-b", "--back", action="store", type="string", dest="background", help="NEEDED: Run number used to compute background statistics", metavar="10 (or 0010)", default="")
parser.add_option("-f", "--backforce", action="store_true", dest="forceComputeBackground", help="OPTIONAL: Force a computation of background statistics", metavar="Force recompute", default=False)
parser.add_option("-n", "--nFrames", action="store", type="int", dest="nFrames", help="OPTIONAL: Number of frames to process", metavar="numFrames", default=-1)
parser.add_option("-d", "--detector", action="store", type="int", dest="detID", help="OPTIONAL: Which detector to process\n -1 for both, 1 for front only, 2 for back only", metavar="2", default=2)
parser.add_option("-P", "--partition", action="store", type="float", dest="partition", help="OPTIONAL: Partition data stream at this fraction", metavar="fraction", default=-1)
parser.add_option("-B", "--bootstrapnum", action="store", type="int", dest="bootstrapNum", help="OPTIONAL: Number of samples to draw in bootstrapping", metavar="numsamples", default=1000)

parser.add_option("-D", "--donotsaveresults", action="store_true", dest="doNotSave", help="OPTIONAL: Do not save intermediate results for background or data runs", metavar="False", default=False)

parser.add_option("-M", "--max", action="store", type="float", dest="maxVal", help="OPTIONAL: Initialize plots to this maximum value", metavar="1000.", default=-1)
parser.add_option("-m", "--min", action="store", type="float", dest="minVal", help="OPTIONAL: Initialize plots to this minimum value", metavar="0.", default=-1)
parser.add_option("-t", "--tag", action="store", type="string", dest="fTag", help="OPTIONAL: File tag to prepend output", metavar="FILETAG", default="")
parser.add_option("-s", "--scale", action="store", type="float", dest="imgScale", help="OPTIONAL: Scaling factor of image", metavar="1.0", default=1.)
#TODO: Put on viewer option?


(op, args) = parser.parse_args()

#Import modules we wrote
from myModules import imagingClass as IM
from myModules import geo 
from myModules import stat
from myModules import cxi_unstable as cxi

############################################################################################
# Set directories for necessary files.
############################################################################################
print '..........................\n'
#expDir = "/reg/d/psdm/cxi/cxi73013/"
#expH5Dir  = expDir + "/scratch/duaneloh/cxidb/"
#expH5Dir  = expDir + "/scratch/common/no_bg_correction/"
#writeDir  = expDir + "/scratch/duaneloh/analysis/"
#maskFName1 = expDir + "/scratch/duaneloh/config/mask/L73013_viewingMask_FrontDetector.h5"
#geomFName1 = expDir + "/scratch/duaneloh/config/geometry/CxiDs1-geometry.h5"
#maskFName2 = expDir + "/scratch/duaneloh/config/mask/L73013_viewingMask_BackDetector.h5"
#geomFName2 = expDir + "???"

formattedBgNum = "r%04d"%int(op.background)
#bgDir = expH5Dir + "/" + formattedBgNum + "/"
#bgH5Glob = bgDir + "/*.cxi"
#bgH5Glob = expH5Dir + "*" + formattedBgNum + ".cxi"
#bgHistName = writeDir + formattedBgNum + "_bg.h5"

expDir = "/Volumes/HD2/cxi73013/"
expH5Dir = expDir + "cxidb/"
writeDir = expDir + "analysis/"+formattedBgNum+"/"
bgH5Glob = expH5Dir + "*" + formattedBgNum + ".cxi"
bgHistName = writeDir + formattedBgNum + "_bg.h5"
maskFName1 = expDir + "/config/mask/L73013_viewingMask_FrontDetector.h5"
geomFName1 = expDir + "/config/geometry/CxiDs1-geometry.h5"
maskFName2 = expDir + "/config/mask/L73013_viewingMask_BackDetector.h5"

############################################################################################
# Set some default experimental values
############################################################################################
pixSize = 0.0001099			#detector pixel size in meters
dDist1 = (162)*1E-3			#sample detector distance in meters
dDist2 = (500)*1E-3			#TODO: What is the back detector distance?
angAvgMax = 300.			#ignore angular average values beyond this
defaultWavelength = 1.693	#In angstroms

#TODO: Make mask for back detector
[rr, rc] = [1480, 1552]		#Size of raw array
[ar, ac] = [1736, 1736]		#Size of assembled array
detector1 = cxi.cxiDetector(rr, rc, ar, ac, pixSize, dDist1, defaultWavelength, geomFileLoc=geomFName1, maskFileLoc=maskFName1)

[rr2, rc2] = [370, 388]		#Size of raw array
[ar2, ac2] = [500, 500]		#Size of assembled
detector2 = cxi.cxiDetector(rr2, rc2, ar2, ac2, pixSize, dDist2, defaultWavelength, maskFileLoc=maskFName2)

bg = cxi.cxidbBGRun(expH5Dir, formattedBgNum, bgH5Glob, writeDir, bgHistName, detector1, detector2, bgHistMin=1200., bgHistMax=2000., bgHistMin2=1200., bgHistMax2=3000.,bgHistBinWidth2=4.)

############################################################################################
# Read background run file
############################################################################################
if os.path.isfile(bg.outputFN) and not op.forceComputeBackground:
	print "%s found."%(bg.outputFN)
	print "Will use it instead of recomputing background."
	print "Reading file.................."
	bg.readStatsFromFile(bg.outputFN, detID=op.detID)
	bg.computeCDF(detID=op.detID)
	bg.computeGaussianDiv(detID=op.detID)
else:
	if op.forceComputeBackground:
		print "Recomputing background statistics"
	else:
		print "Collecting background statistics."
	#bg.collectBasicStats(maxFNum=op.nFrames, detID=op.detID)
	bg.collectStatistics(op.nFrames, detID=op.detID)	#input arg: max. frames to scan; -1 scans all frames
	bg.computeCDF(detID=op.detID)
	bg.computeGaussianDiv(detID=op.detID)

	if not op.doNotSave:
		print "Saving statistics to file: %s"%(bg.outputFN)
		bg.saveStatsToFile(detID=op.detID)

if op.partition > 0:
	bg.partition(maxFNum=op.nFrames, detID=op.detID, part=op.partition, bootstrapNum=op.bootstrapNum)


