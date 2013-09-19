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
parser.add_option("-f", "--backforce", action="store_true", dest="forceComputeBackground", help="OPTIONAL: Force a compution of background statistics", metavar="Force recompute", default=False)
parser.add_option("-n", "--nFrames", action="store", type="int", dest="nFrames", help="OPTIONAL: Number of frames to process", metavar="numFrames", default=-1)
parser.add_option("-d", "--detector", action="store", type="int", dest="detID", help="OPTIONAL: Which detector to process\n -1 for both, 1 for front only, 2 for back only", metavar="-1", default=-1)


parser.add_option("-D", "--donotsaveresults", action="store_true", dest="doNotSave", help="OPTIONAL: Do not save intermediate results for background or data runs", metavar="False", default=False)

parser.add_option("-M", "--max", action="store", type="float", dest="maxVal", help="OPTIONAL: Initialize plots to this maximum value", metavar="1000.", default=-1)
parser.add_option("-m", "--min", action="store", type="float", dest="minVal", help="OPTIONAL: Initialize plots to this minimum value", metavar="0.", default=-1)
parser.add_option("-t", "--tag", action="store", type="string", dest="fTag", help="OPTIONAL: File tag to prepend output", metavar="FILETAG", default="")
parser.add_option("-s", "--scale", action="store", type="float", dest="imgScale", help="OPTIONAL: Scaling factor of image", metavar="1.0", default=1.)
parser.add_option("-k", "--mask", action="store_true", dest="mask", help="OPTIONAL: Mask pixels using viewingMask", metavar="mask", default=False)
parser.add_option("-K", "--makeMask", action="store_true", dest="reMask", help="OPTIONAL: replace viewingMask using background stats.", metavar="makeMask", default=False)

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
######################################################################################################
# Set directories for BG and Data runs
######################################################################################################
formattedBgNum = "r%04d"%int(op.background) 
bgDir = expH5Dir + "/" + formattedBgNum + "/"
bgH5Glob = bgDir + "/*.cxi"
bgHistName = writeDir + formattedBgNum + "_bg.h5"
bg = cxi.cxidbBGRun(expH5Dir, formattedBgNum, bgH5Glob, writeDir, bgHistName, currRun, bgHistMin=1500., bgHistMax=1800., bgHistMin2=1200., bgHistMax2=1700.)

######################################################################################################
# Read background run file
# input: bgH5Glob, bgHistName
# output: avgBGWL, bgHist, bgCDF
# output_to_file: bgHistName 
######################################################################################################
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
	bg.collectStatistics(op.nFrames, detID=op.detID)	#input arg: max. frames to scan; -1 scans all frames
	bg.computeCDF(detID=op.detID)
	bg.computeGaussianDiv(detID=op.detID)
	
	if not op.doNotSave:
		print "Saving statistics to file: %s"%(bg.outputFN)
		bg.saveStatsToFile(detID=op.detID)

#bg.viewBgStats()
#bg.pruneForSignificantMasks()
#toPlot = (bg.sigMasks[:,2] - (bg.sigMasks[:,3]*(bg.sigMasks[:,2]>0))).reshape(bg.avgBg.shape)

#if op.mask:
#	viewer.draw(toPlot*(bg.rMF.reshape(bg.avgBg.shape)), N.mean(bg.bgWLArr), formattedBgNum, formattedBgNum+op.fTag, plotRings=True, cmax=op.maxVal, cmin=op.minVal)
#else:
#	viewer.draw(toPlot, N.mean(bg.bgWLArr), formattedBgNum, formattedBgNum+op.fTag, plotRings=True, cmax=op.maxVal, cmin=op.minVal)


