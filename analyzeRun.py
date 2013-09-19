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
parser.add_option("-d", "--detector", action="store", type="int", dest="detID", help="OPTIONAL: Which detector to process\n -1 for both, 1 for front only, 2 for back only", metavar="-1", default=-1)

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
from myModules import cxi

######################################################################################################
# Set directories for necessary files.
######################################################################################################
print '..........................\n'
expDir = "/reg/d/psdm/cxi/cxi73013/"
expH5Dir  = expDir + "/scratch/common/no_bg_correction/"
writeDir  = expDir + "/scratch/duaneloh/analysis/"
maskFName = expDir + "/scratch/duaneloh/config/mask/viewingMask.h5"
geomFName = expDir + "/scratch/duaneloh/config/geometry/CxiDs1-geometry.h5"

######################################################################################################
# Set some default experimental values
######################################################################################################
pixSize = 0.0001099			#detector pixel size in meters
dDist = (162)*1E-3		#sample detector distance in meters
angAvgMax = 300.			#ignore angular average values beyond this
defaultWavelength = 1.693	#In angstroms

histRange = 100
histBinWidth = 10.
[rr, rc] = [1480, 1552]		#Size of raw array
[ar, ac] = [1736, 1736]		#Size of assembled array
currRun = cxi.cxiExpSetup(rr, rc, ar, ac, pixSize, dDist, defaultWavelength, geomFileLoc=geomFName, maskFileLoc=maskFName)
if not currRun.hasMask:
	op.mask = False
viewer = cxi.cxiViewAssemble(currRun, writeDir)
backViewer = cxi.cxiViewRaw(currRun, writeDir)
######################################################################################################
# Set directories for BG and Data runs
######################################################################################################
formattedBgNum = "r%04d"%int(op.background)
bgDir = expH5Dir + "/" + formattedBgNum + "/"
bgH5Glob = bgDir + "/*.cxi"
bgHistName = writeDir + formattedBgNum + "_bg.h5"
bg = cxi.cxidbBGRun(expH5Dir, formattedBgNum, bgH5Glob, writeDir, bgHistName, currRun, bgHistMin=1500., bgHistMax=1800., bgHistMin2=1200., bgHistMax2=1700.)

formattedRunNum = "r%04d"%int(op.runNumber)
runDir = expH5Dir + "/" + formattedRunNum + "/"
runH5Glob = runDir + "/*.cxi"
avgSignificanceHistName = writeDir + formattedRunNum + "_data_" + formattedBgNum + "_bg.h5"
#avgSignificanceHistName = writeDir + formattedRunNum + "_data.h5"
dataRun = cxi.cxidbDataRun(expH5Dir, formattedRunNum, runH5Glob, writeDir, avgSignificanceHistName, currRun)

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

######################################################################################################
# Display processed, assembled data
######################################################################################################

def viewFrontDetector(tag=-1, defaultADU=20., in_cmin=op.minVal, in_cmax=op.maxVal, in_imgScale=op.imgScale):

	if tag == -1:
		d = dataRun.avgSig.reshape(dataRun.rr, dataRun.rc)
	else:
		d = dataRun.sigADU[tag].reshape(dataRun.rr, dataRun.rc)/defaultADU
	viewer.draw(d, dataRun.avgWL, formattedRunNum, formattedRunNum+op.fTag, plotRings=False, cmax=in_cmax, cmin=in_cmin, imgScale=in_imgScale)


def viewBackDetector(tag=-1, defaultADU=20., in_cmin=op.minVal, in_cmax=op.maxVal, in_imgScale=op.imgScale):

	if tag == -1:
		d = dataRun.avgSig2.reshape(dataRun.rr2, dataRun.rc2)
	else:
		d = dataRun.sigADU2[tag].reshape(dataRun.rr2, dataRun.rc2)/defaultADU
	backViewer.draw(d, dataRun.avgWL, formattedRunNum, formattedRunNum+op.fTag, plotRings=False, cmax=in_cmax, cmin=in_cmin, imgScale=in_imgScale)


