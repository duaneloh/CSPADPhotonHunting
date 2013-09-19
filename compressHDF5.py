#!/usr/bin/env python

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

######################################################################################################
# Read and parse command line input flags
######################################################################################################
parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNum", help="NEEDED: Run number used to compress", metavar="10 (or 0010)", default="")
parser.add_option("-n", "--numFrames", action="store", type="int", dest="numFrames", help="OPTIONAL: Number of frames to compy", metavar="1000 ", default=-1)
(op,args) = parser.parse_args()


expTag = "cxi73013"
expDir = "/reg/d/psdm/cxi/" + expTag + "/"
h5SourceDir = expDir + "hdf5/"
h5DestDir = expDir + "scratch/duaneloh/cxidb/"
dataField = 'Configure:0000/Run:0000/CalibCycle:0000/CsPad2x2::ElementV1/CxiSc2.0:Cspad2x2.0/data'
ebeamField = 'Configure:0000/Run:0000/CalibCycle:0000/Bld::BldDataEBeamV3/EBeam/data'

formattedRunNum = "r%04d"%int(op.runNum)
h5SourceFile = h5SourceDir +expTag+"-"+formattedRunNum+".h5"
h5DestFile = h5DestDir + expTag+"-"+formattedRunNum+".cxi"


fsource = h5py.File(h5SourceFile, 'r')
fdest = h5py.File(h5DestFile, 'w')

numFramesToCopy = fsource[dataField].shape[0]
if (op.numFrames == -1):
	print "Getting detector information for %d frames"%numFramesToCopy
elif (op.numFrames > numFramesToCopy):
	print "Specified frames to copy (%d) exceeded actual number of frames (%d)"%(op.numFrames, numFramesToCopy)
	print "Will copy actual number instead %d"%numFramesToCopy
else:
	numFramesToCopy = op.numFrames

print "Extracting............"
backDetectorData = fsource[dataField]
waveLengthInfo = fsource[ebeamField][:numFramesToCopy]
outputWavelengthInA = N.zeros(numFramesToCopy)

gdest = fdest.create_group("entry_1/instrument_1/detector_2")
outputBackDetectorData = gdest.create_dataset("data", (numFramesToCopy, 370, 388), "int16")
print "Writing frames........"

for nn in range(numFramesToCopy):
	tempD = backDetectorData[nn]
	outputBackDetectorData[nn,:185,:] = tempD[:,:,0]
	outputBackDetectorData[nn,185:,:] = tempD[:,:,1]
	print "%05d of %05d frames"%(nn,numFramesToCopy)
fdest['entry_1/instrument_1/detector_2/data'].attrs['numEvents'] = [numFramesToCopy]

print "Writing wavelengths........"
for nn,wlVal in enumerate(waveLengthInfo):
	#print wlVal
	#print len(wlVal)
	peakCurrent = wlVal[7]
	DL2energyGeV = 0.001*wlVal[2]
	LTUwakeLoss = 0.0016293*peakCurrent
	SRlossPerSegment = 0.63*DL2energyGeV
	wakeLossPerSegment = 0.0003*peakCurrent
	energyLossPerSegment = SRlossPerSegment + wakeLossPerSegment
	energyProfile = DL2energyGeV - 0.001*LTUwakeLoss - 0.0005*energyLossPerSegment
	photonEnergyeV = 44.42*energyProfile*energyProfile
	outputWavelengthInA[nn] = 12398.42/photonEnergyeV
fdest.create_dataset("LCLS/photon_wavelength_A", data=outputWavelengthInA)

#[uDMask, fBeamCharge, fEbeamL3Energy, fEbeamLTUPosX, fEbeamLTUPosY, fEbeamLTUAngX, fEbeamLTUAngY, fEbeamPkCurrBC2, fEbeamEnergyBC2, fEbeamPkCurrBC1, fEbeamEnergyBC1] = waveLengthInfo

print "Copying done"

fsource.close()
fdest.close()


