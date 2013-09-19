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
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="NEEDED: Run number you wish to view", metavar="1 (or 0001)", default="")
parser.add_option("-b", "--back", action="store", type="string", dest="background", help="NEEDED: Run number used to compute background statistics", metavar="10 (or 0010)", default="")
parser.add_option("-D", "--donotsaveresults", action="store_true", dest="doNotSave", help="OPTIONAL: Do not save intermediate results for background or data runs", metavar="False", default=False)
parser.add_option("-d", "--detector", action="store", type="int", dest="detID", help="OPTIONAL: Which detector to process\n -1 for both, 1 for front only, 2 for back only", metavar="-1", default=2)

parser.add_option("-f", "--runforce", action="store_true", dest="forceComputeRun", help="OPTIONAL: Force a compution of data run's significance", metavar="Force recompute", default=False)
parser.add_option("-n", "--nFrames", action="store", type="int", dest="nFrames", help="OPTIONAL: Number of frames to process", metavar="numFrames", default=-1)

parser.add_option("-M", "--max", action="store", type="float", dest="maxVal", help="OPTIONAL: Initialize plots to this maximum value", metavar="maxVal", default=-1)
parser.add_option("-m", "--min", action="store", type="float", dest="minVal", help="OPTIONAL: Initialize plots to this minimum value", metavar="minVal", default=-1)
parser.add_option("-t", "--tag", action="store", type="string", dest="fTag", help="OPTIONAL: File tag to prepend output", metavar="FILETAG", default="")
parser.add_option("-s", "--scale", action="store", type="float", dest="imgScale", help="OPTIONAL: Scaling factor of image", metavar="1.0", default=1.)
parser.add_option("-k", "--mask", action="store_true", dest="mask", help="OPTIONAL: Mask pixels using viewingMask", metavar="mask", default=False)

(op, args) = parser.parse_args()

#Import modules we wrote
from myModules import imagingClass as IM
from myModules import geo 
from myModules import stat
from myModules import cxi_unstable as cxi

######################################################################################################
# Set directories for necessary files.
######################################################################################################
print '..........................\n'
#expDir = "/reg/d/psdm/cxi/cxi73013/"
#expH5Dir  = expDir + "/scratch/duaneloh/cxidb/"
#expH5Dir  = expDir + "/scratch/common/no_bg_correction/"
#writeDir  = expDir + "/scratch/duaneloh/analysis/"
#maskFName1 = expDir + "/scratch/duaneloh/config/mask/L73013_viewingMask_FrontDetector.h5"
#geomFName1 = expDir + "/scratch/duaneloh/config/geometry/CxiDs1-geometry.h5"
#maskFName2 = expDir + "/scratch/duaneloh/config/mask/L73013_viewingMask_BackDetector.h5"
#geomFName2 = expDir + "???"

expDir = "/Volumes/HD2/cxi73013/"
maskFName1 = expDir + "/config/mask/L73013_viewingMask_FrontDetector.h5"
geomFName1 = expDir + "/config/geometry/CxiDs1-geometry.h5"
maskFName2 = expDir + "/config/mask/L73013_viewingMask_BackDetector.h5"



######################################################################################################
# Set some default experimental values
######################################################################################################
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

######################################################################################################
# Set directories for BG and Data runs
######################################################################################################

expH5Dir = expDir + "cxidb/"
formattedBgNum = "r%04d"%int(op.background)
writeDir = expDir + "analysis/"+formattedBgNum+"/"
#bgDir = expH5Dir + "/" + formattedBgNum + "/"
#bgH5Glob = bgDir + "/*.cxi"
bgH5Glob = expH5Dir + "*" + formattedBgNum + ".cxi"
bgHistName = writeDir + formattedBgNum + "_bg.h5"
bg = cxi.cxidbBGRun(expH5Dir, formattedBgNum, bgH5Glob, writeDir, bgHistName, detector1, detector2)


formattedRunNum = "r%04d"%int(op.runNumber)
#runDir = expH5Dir + "/" + formattedRunNum + "/"
#runH5Glob = runDir + "/*.cxi"
runH5Glob = expH5Dir + "*" + formattedRunNum + ".cxi"
avgSignificanceHistName = writeDir + formattedRunNum + "_data_" + formattedBgNum + "_bg.h5"

dataRun = cxi.cxidbDataRun(expH5Dir, formattedRunNum, runH5Glob, writeDir, avgSignificanceHistName, detector1, detector2)

######################################################################################################
# Read background run file
# input: bgH5Glob, bgHistName
# output: avgBGWL, bgHist, bgCDF
# output_to_file: bgHistName 
######################################################################################################
print '..........................\n'
if os.path.isfile(bgHistName):
	print "%s found."%(bgHistName)
	print "Will use it instead of recomputing background."
	print "Reading file.................."
	bg.readStatsFromFile(bgHistName, detID=op.detID)
	bg.computeCDF(detID=op.detID)
	bg.computeGaussianDiv(detID=op.detID)
	bg.makeStatMasks(detID=op.detID)
else:
	print "Unable to find %s. Aborting"%(bgHistName)
	sys.exit(1)

######################################################################################################
# Read data run file, interpolate for viewing
# input: runH5Glob
# output: wLArr&avgWL, avgSignificance,  
######################################################################################################
if os.path.isfile(dataRun.outputFN) and not op.forceComputeRun:
	print "%s found."%dataRun.outputFN
	print "Will use it instead of recomputing data run."
	dataRun.readStatsFromFile(dataRun.outputFN, detID=op.detID)
else:
	if op.forceComputeRun:
		print "Re-testing run against background"
	else:
		print "Testing data against background statistics....."
	dataRun.collectStatistics(bg, op.nFrames, detID=op.detID)
	if not op.doNotSave:
		print "Saving data to file............................"
		dataRun.saveStatsToFile(detID=op.detID)	
