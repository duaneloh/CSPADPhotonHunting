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

import imagingClass as IM
import geo 
import stat

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as PL
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


def cxiViewAssembled(inData, inWL, detSetup, filename="fn", writeDir="", runTag="runxxx", plotRings=False, cmax=2000., cmin=0., imgScale=1., angAvgMax=1000.):

	d = geo.interpolateGeometry(detSetup.x, detSetup.y, inData.astype('float'), detSetup.pixSize)
	aa = geo.angAve2D(inData.astype('double'), detSetup.radInt)
	qArr = (detSetup.qRaw / inWL).astype('double')
	aq = geo.angAve2D(qArr, self.radInt)
	currImg = IM.img_class(d, inWL, detSetup.dDist, angAvg=[aq[(aa > 0)*(aa < angAvgMax)],aa[(aa > 0)*(aa < angAvgMax)]], writeDir=writeDir, filename=filename, runTag=runTag, plotRings=plotRings, cmax=cmax, cmin=cmin, imgScale=imgScale)
	currImg.draw_img()



def cxiViewUnassembled(inData, filename="fn", writeDir="", runTag="runxxx", plotRings=False, cmax=2000., cmin=0., imgScale=1., angAvgMax=1000.):
	
	currImg = IM.img_class(inData, "", "", writeDir=writeDir, filename=filename, runTag=runTag, plotRings=False, cmax=cmax, cmin=cmin, imgScale=imgScale)
	currImg.draw_img()



class cxiDetector(object):
	def __init__(self, rr, rc, ar, ac, pixSize, dDist, defaultWavelength, geomFileLoc=None, maskFileLoc=None,cmNorm=1450.):
	
		self.defaultWavelength = defaultWavelength
		
		[self.rr, self.rc] = [rr, rc]
		[self.ar, self.ac] = [ar, ac]
		[self.xleftPad, self.yleftPad] = [3., 3.]
		[self.x, self.y, self.z] = [None, None, None]
		[self.rad, self.radInt, self.qRaw] = [None, None, None]
		self.rMF = N.ones(self.rr*self.rc).astype('int')
		[self.pixSize, self.dDist] = [pixSize, dDist]
		self.cmNorm = cmNorm
		
		self.hasMask = False

		self.maskFileLoc = maskFileLoc
		if self.maskFileLoc is not None:
			try:
				print "Reading mask %s:" % self.maskFileLoc
				mask = h5py.File(self.maskFileLoc, "r")
				rM = (mask['data/data'].value).astype('int')
				mask.close()
				self.rMF = rM.flatten()
				self.hasMask = True
				print ".............. succeeded!\n"
			except:
				print ".............. failed!\n"
		
		self.geomFileLoc = geomFileLoc
		if self.geomFileLoc is not None:
			try:
				print "Reading geometry for detector 1: %s" % self.geomFileLoc
				geom = h5py.File(self.geomFileLoc, "r")		#Open geometry H5 file for reading
				[self.x, self.y, self.z] = [N.array(geom['x']), N.array(geom['y']), N.array(geom['z'])]
				self.rr = self.x.shape[0]
				self.rc = self.x.shape[1]
				geom.close()
				self.x -= self.x.mean()
				self.y -= self.y.mean()
				self.rad = N.sqrt(self.x*self.x + self.y*self.y)
				self.radInt = (self.rad / self.pixSize).astype('int')
				self.qRaw = 2.*N.pi * 2. * N.sin(0.5*N.arctan2(self.rad, self.dDist))
				self.x -= self.x.min() - (self.xleftPad * self.pixSize)
				self.y -= self.y.min() - (self.yleftPad * self.pixSize)
				rM = self.rMF.reshape(self.rr,-1)
				if self.hasMask:
					self.radInt = (self.radInt*(rM == 1) - (rM != 1)).astype('int')
				print ".............. succeeded!\n"
			except:
				print ".............. failed!\n"
				


class cxidbRun(object):
	def __init__(self, dataDir, runTag, dataGlob, writeDir, outputFN):
		self.dataDir = dataDir
		self.runTag	= runTag
		self.dataGlob = dataGlob
		self.writeDir = writeDir
		self.outputFN = outputFN



class cxidbBGRun(cxidbRun):

	def __init__(self, dataDir, runTag, dataGlob, writeDir, outputFN, det1Setup, det2Setup, cmNorm2=1400., bgHistMin=1500., bgHistMax=2000., bgHistBinWidth=1., bgHistMin2=1000., bgHistMax2=3000., bgHistBinWidth2=2.):
		cxidbRun.__init__(self, dataDir, runTag, dataGlob, writeDir, outputFN)
		self.det1Setup = det1Setup
		self.det2Setup = det2Setup
		self.bgCfXIDB = ""
		self.bgWLArr = []
		self.pulseEnergy_eV = None
		self.fcounter = 0.
		self.runTag = runTag
		self.statLabels = ["avg.", "var.", "K.L. Div", "Gaussian norm.", "Chi-square"]
		
		self.rr = self.det1Setup.rr
		self.rc = self.det1Setup.rc
		self.rMF = det1Setup.rMF
		self.numPixsInBgHist = int(N.sum(self.rMF!=0))
		self.bgHistMin = bgHistMin
		self.bgHistMax = bgHistMax
		self.bgHistBinWidth = bgHistBinWidth
		self.bgHist = N.zeros((self.numPixsInBgHist, ((self.bgHistMax - self.bgHistMin)/self.bgHistBinWidth) + 1))
		self.avgBg = N.zeros(self.rr*self.rc)
		self.bgCDF = []
		self.bgStat = []
		self.sigMasksF = None
		self.numFrames = 0
		self.histStats = N.zeros((self.numPixsInBgHist,6))

		self.rr2 = self.det2Setup.rr
		self.rc2 = self.det2Setup.rc
		self.rMF2 = det2Setup.rMF.astype('int')
		self.numPixsInBgHist2 = int(N.sum(self.rMF2!=0))
		self.bgHistMin2 = bgHistMin2
		self.bgHistMax2 = bgHistMax2
		self.bgHistBinWidth2 = bgHistBinWidth2
		self.bgHist2 = N.zeros((self.numPixsInBgHist2, ((self.bgHistMax2 - self.bgHistMin2)/self.bgHistBinWidth2) + 1))
		self.avgBg2 = N.zeros(self.rr2*self.rc2)
		self.bgCDF2 = []
		self.bgStat2 = []
		self.sigMasksF2 = None
		self.numFrames2 = 0
		self.histStats2 = N.zeros((self.numPixsInBgHist2, 6))
		
		self.maxFNum = 0
		self.likelyBG = []
		self.likelyBGScores = []
		self.boostrapScores = []
		self.likelyHits = []
		self.likelyHitScores = []
			

	#TODO: Combine IO functions into general object
	def collectBasicStats(self, maxFNum=-1, detID=2, fileSelector=[]):
		self.bgCXIDB = G.glob(self.dataGlob)[0]
		print "Background CXIDB file found: %s\n"%(self.bgCXIDB)
		f = h5py.File(self.bgCXIDB,"r")
		if (detID == 1) or (detID == -1):
			try:
				self.histStats *= 0.
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"
		
		if (detID == 2) or (detID == -1):
			try:
				self.histStats2 *= 0.
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]
			except:
				print "Unable to read detector_2"
		
		if fileSelector == []:
			if (detID == 1) and ((maxFNum == -1) or (self.numFrames < maxFNum)):
				maxFNum = self.numFrames
			elif (detID == 2) and ((maxFNum == -1) or (self.numFrames2 < maxFNum)):
				maxFNum = self.numFrames2
			elif (detID == -1):
				if (self.numFrames < self.numFrames2):
					if self.numFrames < maxFNum:
						maxFNum = self.numFrames
				else:
					if self.numFrames2 < maxFNum:
						maxFNum = self.numFrames2
			
			self.maxFNum = maxFNum
			toVisit = range(maxFNum)
		else:
			toVisit = fileSelector
			maxFNum = len(toVisit)
		
		startTime = datetime.datetime.now()
		for self.fcounter in toVisit:
			if (detID == 1) or (detID == -1):
				rawBgF = data[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF, self.rMF, self.det1Setup.cmNorm)
				self.avgBg += rawBgF
				stat.updatePixValStats(rawBgF, self.rMF, self.histStats)
			if (detID == 2) or (detID == -1):
				rawBgF2 = data2[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF2, self.rMF2, self.det2Setup.cmNorm)
				self.avgBg2 += rawBgF2
				stat.updatePixValStats(rawBgF2, self.rMF2, self.histStats2)
			if (self.fcounter%10 == 0):
				print "Done with BG %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)
		
		
		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing BG files took %f seconds each"%((differenceTime.total_seconds())/maxFNum)
		
		print "Normalizing statistics....."
		if (detID == 1) or (detID == -1):
			mm = self.histStats[:,4]>0.
			mmc = 1.-mm
			av = (self.histStats[:,2]/(self.histStats[:,4]*mm + mmc))
			self.histStats[:,5] = (self.histStats[:,3]/(self.histStats[:,4]*mm + mmc)) - av*av
			
		if (detID == 2) or (detID == -1):
			mm2 = self.histStats2[:,4]>0.
			mmc2 = 1.-mm2
			av2 = (self.histStats2[:,2]/(self.histStats2[:,4]*mm2 + mmc2))
			self.histStats2[:,5] = (self.histStats2[:,3]/(self.histStats2[:,4]*mm2 + mmc2)) - av2*av2
		print "..done!"
		f.close()
	
		
	def collectStatistics(self, maxFNum=-1, detID=2, fileSelector=[]):
	
		self.bgCXIDB = G.glob(self.dataGlob)[0]
		print "Background CXIDB file found: %s\n"%(self.bgCXIDB)
		
		f = h5py.File(self.bgCXIDB,"r")
		try:
			self.bgWLArr = f['LCLS/photon_wavelength_A'].value
		except:
			print "Unable to read photon wavelengths"
			
		if (len(self.bgWLArr) > 0):
			numBadWLEntries = 0
			for nn,cWL in enumerate(self.bgWLArr):
				if (cWL == float('Inf') or cWL == -float('Inf') or cWL == float('NaN') or cWL == 0.):
					self.bgWLArr[nn] = self.det1Setup.defaultWavelength
					numBadWLEntries += 1
			print "%d of %d wavelength entries were undefined and replaced with default"%(numBadWLEntries, len(self.bgWLArr))
		try:
			self.pulseEnergy_eV = f['LCLS/photon_energy_eV'].value
		except:
			print "Unable to read pulse energies"
			
		if (detID == 1) or (detID == -1):
			try:
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"

		if (detID == 2) or (detID == -1):
			try:
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]				
			except:
				print "Unable to read detector_2"

		if fileSelector == []:
			if (detID == 1) and ((maxFNum == -1) or (self.numFrames < maxFNum)):
				maxFNum = self.numFrames
			elif (detID == 2) and ((maxFNum == -1) or (self.numFrames2 < maxFNum)):
				maxFNum = self.numFrames2
			elif (detID == -1):
				if (self.numFrames < self.numFrames2):
					if self.numFrames < maxFNum:
						maxFNum = self.numFrames
				else:
					if self.numFrames2 < maxFNum:
						maxFNum = self.numFrames2

			self.maxFNum = maxFNum
			toVisit = range(maxFNum)
		else:
			toVisit = N.array(fileSelector)
			toVisit.sort()
			maxFNum = len(toVisit)

		startTime = datetime.datetime.now()
		for nn,self.fcounter in enumerate(toVisit):
			if (detID == 1) or (detID == -1):
				rawBgF = data[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF, self.rMF, self.det1Setup.cmNorm)
				self.avgBg += rawBgF
				stat.updatePixValHist(rawBgF, self.rMF, self.bgHist, self.bgHistMin, self.bgHistMax, self.bgHistBinWidth)
			if (detID == 2) or (detID == -1):
				rawBgF2 = data2[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF2, self.rMF2, self.det2Setup.cmNorm)
				self.avgBg2 += rawBgF2
				stat.updatePixValHist(rawBgF2, self.rMF2, self.bgHist2, self.bgHistMin2, self.bgHistMax2, self.bgHistBinWidth2)
			if (nn%10 == 0):
				print "Done with BG %s %04d of %04d." % (self.runTag, nn, maxFNum)


		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing BG files took %f seconds each"%((differenceTime.total_seconds())/maxFNum)

		print "Normalizing statistics....."
		if (detID == 1) or (detID == -1):
			self.avgBg /= 1.0*maxFNum
			stat.normalizePixValHist(self.bgHist, 1.E-10)
		if (detID == 2) or (detID == -1):
			self.avgBg2 /= 1.0*maxFNum
			stat.normalizePixValHist(self.bgHist2, 1.E-10)
		self.avgBgWL = N.array(self.bgWLArr).mean()
		print "..done!"
		f.close()


	def computeCDF(self, detID=2):
		if (detID == 1) or (detID == -1):
			print "Computing c.d.f. for data from detector 1."
			self.bgCDF = stat.computeCDFFromValHist(self.bgHist)
			print "Done!"
		if (detID == 2) or (detID == -1):
			print "Computing c.d.f. for data from detector 2."
			self.bgCDF2 = stat.computeCDFFromValHist(self.bgHist2)
			print "Done!"
					
		
	def computeGaussianDiv(self, detID=2):
		if (detID == 1) or (detID == -1):
			print "Computing GaussianDiv for detector 1."
			self.bgStat = stat.computeStatFromValHist(self.bgHist, self.bgHistMin, self.bgHistBinWidth)
			print "Done!"
	
		if (detID == 2) or (detID == -1):
			print "Computing GaussianDiv for detector 2."
			self.bgStat2 = stat.computeStatFromValHist(self.bgHist2, self.bgHistMin2, self.bgHistBinWidth2)
			print "Done!"


	def selfTest(self, maxFNum=-1, detID=2, bootstrapNum=-1, fileSelector=[]):
	
		self.bgCXIDB = G.glob(self.dataGlob)[0]
		print "Background CXIDB file found: %s\n"%(self.bgCXIDB)
		
		f = h5py.File(self.bgCXIDB,"r")
		if (detID == 1):
			try:
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"
				
		if (detID == 2):
			try:
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]
			except:
				print "Unable to read detector_2"

		if fileSelector == []:
			if (detID == 1) and ((maxFNum == -1) or (self.numFrames < maxFNum)):
				maxFNum = self.numFrames
			elif (detID == 2) and ((maxFNum == -1) or (self.numFrames2 < maxFNum)):
				maxFNum = self.numFrames2
			toVisit = range(maxFNum)
		else:
			toVisit = fileSelector
			maxFNum = len(toVisit)

		if (bootstrapNum > 0):
			toVisit = N.random.choice(toVisit, bootstrapNum)
			maxFNum = bootstrapNum

		results = N.zeros(maxFNum)
		toVisit.sort()
		startTime = datetime.datetime.now()
		for nn,self.fcounter in enumerate(toVisit):
			if (detID == 1):
				rawBgF = data[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF, self.rMF, self.det1Setup.cmNorm)
				results[nn] = stat.computeLogLikelihood(rawBgF, self.rMF, self.bgHist, self.bgHistMin, self.bgHistBinWidth)

			if (detID == 2):
				rawBgF2 = data2[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawBgF2, self.rMF2, self.det2Setup.cmNorm)
				results[nn] = stat.computeLogLikelihood(rawBgF2, self.rMF2, self.bgHist2, self.bgHistMin2, self.bgHistBinWidth2)
			if (nn%10 == 0):
				print "Done with BG %s %04d of %04d." % (self.runTag, nn, maxFNum)

		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing BG files took %f seconds each"%((differenceTime.total_seconds())/maxFNum)
		f.close()
		
		return N.array(results)
		
		
	def partition(self, maxFNum=-1, detID=2, part=0.1, bootstrapNum=1000):
		print "Identifying most likely background in frames"
		self.likelyBGScores = self.selfTest(maxFNum=maxFNum, detID=detID)
		ord = self.likelyBGScores.argsort()
		cPos = N.floor(part*len(self.likelyBGScores))
		self.likelyBG = ord[cPos:]
		self.likelyHits = ord[:cPos]
				
		print "Recomputing statistics for most likely background"
		self.collectStatistics(maxFNum=maxFNum, detID=detID, fileSelector=self.likelyBG)
		self.computeCDF(detID=detID)
		self.computeGaussianDiv(detID=detID)
		
		print "Bootstrapping to calibrate background-background comparisons"
		self.boostrapScores = self.selfTest(maxFNum=maxFNum, detID=detID, bootstrapNum=bootstrapNum, fileSelector=self.likelyBG)
		
		print "Comparing likely background against likely hits"
		self.likelyHitScores = self.selfTest(maxFNum=maxFNum, detID=detID, fileSelector=self.likelyHits)
		
		return None


	def avgSlices(self, fileSelector, detID=2):
		
		toVisit = N.sort(fileSelector)
		maxFNum = len(toVisit)

		self.bgCXIDB = G.glob(self.dataGlob)[0]
		print "Background CXIDB file found: %s\n"%(self.bgCXIDB)
		f = h5py.File(self.bgCXIDB,"r")

		startTime = datetime.datetime.now()		
		if (detID == 1):
			data = f['entry_1/instrument_1/detector_1/data']
			(self.numFrames, rr, rc) = data.shape
			if (data.attrs["numEvents"][0] != self.numFrames):
				print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
				print "Will default to attributed value instead"
				self.numFrames = data.attrs["numEvents"][0]

			avg = N.zeros((rr,rc))
			for self.fcounter in toVisit:
				dataF = data[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(dataF, self.rMF, self.det1Setup.cmNorm)
				avg += dataF.reshape(self.det1Setup.rr,-1)
				print "Processed %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)
			avg /= 1.0*maxFNum

		if (detID == 2):
			data2 = f['entry_1/instrument_1/detector_2/data']
			(self.numFrames2, rr2, rc2) = data2.shape
			if (data2.attrs["numEvents"][0] != self.numFrames2):
				print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
				print "Will default to attributed value instead"
				self.numFrames2 = data2.attrs["numEvents"][0]

			avg = N.zeros((rr2,rc2))
			for self.fcounter in toVisit:
				dataF2 = data2[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(dataF2, self.rMF2, self.det2Setup.cmNorm)
				avg += dataF2.reshape(self.det2Setup.rr,-1)
				print "Processed %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)
			avg /= 1.0*maxFNum
		
		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing %s files took %f seconds each"%(self.runTag, (differenceTime.total_seconds())/maxFNum)
		f.close()
		
		return avg


	def makeStatMasks(self, detID=2):
		
		if (detID == 1) or (detID == -1):
			self.sigMasksF = stat.makeStatMasks(self.bgStat, self.rMF)
		
		if (detID == 2) or (detID == -1):
			self.sigMasksF2 = stat.makeStatMasks(self.bgStat2, self.rMF2)
			

	def viewSinglePixelHistogram(self, pixNum, detID=2):
		fig = plt.figure()
		gs = gridspec.GridSpec(4,4, height_ratios=[1,1,1,1])
		gs.update(top=0.9, bottom=0.1, hspace=0.25)
		ax = fig.add_subplot(gs[:,:])
		ax.set_xlabel("ADU")
		ax.set_ylabel("Probability(ADU)")

		if (detID == 1):
			ax.set_title("Distribution for pixel %d:\navg:%.4f, var:%.4f, KLDivGaussian:%.4f"%(pixNum, self.bgStat[pixNum,0], self.bgStat[pixNum,1], self.bgStat[pixNum,2]))
			ax.set_yscale('log')
			plt.plot(self.bgHistMin+self.bgHistBinWidth*N.arange(self.bgHist.shape[1]), self.bgHist[pixNum])

		if (detID == 2):
			ax.set_title("Distribution for pixel %d:\navg:%.4f, var:%.4f, KLDivGaussian:%.4f"%(pixNum, self.bgStat2[pixNum,0], self.bgStat2[pixNum,1], self.bgStat2[pixNum,2]))
			ax.set_yscale('log')
			plt.plot(self.bgHistMin2+self.bgHistBinWidth2*N.arange(self.bgHist2.shape[1]), self.bgHist2[pixNum])

		plt.show()


	def viewSelectedPixelStatistics(self, kmin, kmax, detID=2, statNum=2, fPlotGamma=0.2, imgScale=1., ADULow=15, ADUHigh=30, logScale=False, ymin=0., ymax=1.):
		fig = plt.figure(figsize=(imgScale*15.5, imgScale*9.5), dpi=80)
		gs = gridspec.GridSpec(4,3)
		gs.update(top=0.95, bottom=0.08, hspace=0.25, wspace=0.25)
							
		if (detID == 1):
			if kmin < 0:
				raise ValueError("kmin, the minimum ranked pixel plotted, has to be greater than or equal zero.")
			if kmax > self.bgHist.shape[0]:
				raise ValueError("kmax, the maximum ranked pixel plotted, has to be less than %d"%self.bgHist.shape[0])
			kl = self.bgStat[:,statNum].copy()
			ordering = self.bgStat[:,statNum].argsort()
			selectedPixHist = self.bgHist[ordering[kmin:kmax]]
			xx = self.bgHistMin+self.bgHistBinWidth*N.arange(self.bgHist.shape[1])
			
		elif (detID == 2):
			if kmin < 0:
				raise ValueError("kmin, the minimum ranked pixel plotted, has to be greater than or equal zero.")
			if kmax > self.bgHist2.shape[0]:
				raise ValueError("kmax, the maximum ranked pixel plotted, has to be less than %d"%self.bgHist2.shape[0])
			kl = self.bgStat2[:,statNum].copy()
			ordering = self.bgStat2[:,statNum].argsort()
			selectedPixHist = self.bgHist2[ordering[kmin:kmax]]
			xx = self.bgHistMin2+self.bgHistBinWidth2*N.arange(self.bgHist2.shape[1])

		selectedPixHistF = N.array([N.abs(N.fft.ifftn(N.abs(N.fft.fftn(item))**2)) for item in selectedPixHist])[:,:int(ADUHigh*2)]
		if (detID == 1):
			xxf = N.arange(0, self.bgHistBinWidth*len(selectedPixHistF[0]), self.bgHistBinWidth)
		if (detID == 2):
			xxf = N.arange(0, self.bgHistBinWidth2*len(selectedPixHistF[0]), self.bgHistBinWidth2)
		sortingArray = N.array([N.max(item[ADULow:ADUHigh]) for item in selectedPixHistF])
		sortingArrayArg = sortingArray.argsort()
		selectedPixHist = selectedPixHist[sortingArrayArg]
		selectedPixHistF = selectedPixHistF[sortingArrayArg]		
		
		kl.sort()
		klInMask = fig.add_subplot(gs[:,0])
		klInMask.set_xlabel("pixel (sorted by Div.)")
		klInMask.set_ylabel("%s"%self.statLabels[statNum])
		#klInMask.set_yscale('log')
		klInMask.set_title("%s"%self.statLabels[statNum])
		plt.plot([kmin,kmin],[0,1], "k-", lw=2)
		plt.plot([kmax,kmax],[0,1], "k-", lw=2)
		for label in klInMask.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		klInMask.set_ylim(ymin,ymax)
		plt.plot(kl)
		
		pltAspect=selectedPixHist.shape[1]/(1.0*selectedPixHist.shape[0])
		pltAspectAutoCorr = selectedPixHistF.shape[1]/(1.0*selectedPixHistF.shape[0])
		firstLine = N.floor(0.25*selectedPixHist.shape[0])
		secondLine = N.floor(0.75*selectedPixHist.shape[0])

		plot1 = fig.add_subplot(gs[:2,1])
		plot1.set_xlabel("num. of histogram bins")
		plot1.set_ylabel("pixel")
		plot1.set_title("Stacked histograms for pixels\n ranked %d--%d of %d"%(kmin, kmax, len(ordering)))
		plt.plot([0, selectedPixHist.shape[1]],[firstLine, firstLine], "w-", lw=2)
		plt.plot([0, selectedPixHist.shape[1]],[secondLine, secondLine], "w-", lw=2)
		plt.imshow(selectedPixHist**fPlotGamma, aspect=pltAspect)
		
		plotf = fig.add_subplot(gs[2:,1])
		plotf.set_xlabel("num. of histogram bins")
		plotf.set_title("Stacked Autocorr. for pixels\nranked %d--%d of %d"%(kmin, kmax, len(ordering)))
		plt.plot([0, selectedPixHist.shape[1]],[firstLine, firstLine], "w-", lw=2)
		plt.plot([0, selectedPixHist.shape[1]],[secondLine, secondLine], "w-", lw=2)
		plt.imshow(selectedPixHistF**fPlotGamma, aspect=pltAspectAutoCorr)

		ax1 = fig.add_subplot(gs[0,2])
		ax1.set_xlabel("ADU")
		ax1.set_ylabel("Probability(ADU)")
		for label in ax1.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		if logScale:
			ax1.set_yscale('log')
		plt.plot(xx, selectedPixHist[firstLine])

		ax2 = fig.add_subplot(gs[1,2])
		ax2.set_xlabel("ADU")
		ax2.set_ylabel("Probability(ADU)")
		for label in ax2.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		if logScale:
			ax2.set_yscale('log')
		plt.plot(xx, selectedPixHist[secondLine])

		ax3 = fig.add_subplot(gs[2,2])
		ax3.set_xlabel("ADU")
		ax3.set_ylabel("Autocorr. of P(ADU)")
		for label in ax3.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		if logScale:
			ax3.set_yscale('log')
		plt.plot(xxf, selectedPixHistF[firstLine])
		
		ax4 = fig.add_subplot(gs[3,2])
		ax4.set_xlabel("ADU")
		ax4.set_ylabel("Autocorr. of P(ADU)")
		for label in ax4.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		if logScale:
			ax4.set_yscale('log')
		plt.plot(xxf, selectedPixHistF[secondLine])

		plt.show()


	def viewUnassembledSigMasks(self, sigMaskNum, detID=2, runTag="", cmin=0., cmax=2000., imgScale=1.0, ADU=1.):
		
		avgWL = N.mean(self.bgWLArr)

		if (detID == 1):
			if self.sigMasksF is None:
				self.makeStatMasks(detID=1)
			toPlot = self.sigMasksF[sigMaskNum].reshape(self.det1Setup.rr, self.det1Setup.rc) / ADU
			cxiViewUnassembled(toPlot, writeDir=self.writeDir, filename="%s %s"%(self.runTag, self.statLabels[sigMaskNum]), runTag=runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)

		elif (detID == 2):
			if self.sigMasksF2 is None:
				self.makeStatMasks(detID=2)
			toPlot = self.sigMasksF2[:,sigMaskNum].reshape(self.det2Setup.rr, self.det2Setup.rc) / ADU
			cxiViewUnassembled(toPlot, writeDir=self.writeDir, filename="%s %s det. 2"%(self.runTag, self.statLabels[sigMaskNum]), runTag=runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)

		else:
			print "Invalid detector ID %d given."%(detID)


	def viewAssembledSigMasks(self, sigMaskNum, detID=2, runTag="", cmin=0., cmax=2000., imgScale=1.0, ADU=1., plotRings=False):
		
		avgWL = N.mean(self.bgWLArr)
		
		if (detID == 1):
			if self.sigMasksF is None:
				self.makeStatMasks(detID=1)
			toPlot = self.sigMasksF[sigMaskNum].reshape(self.det1Setup.rr, self.det1Setup.rc) / ADU
			cxiViewAssembled(toPlot, avgWL, self.det1Setup, filename="%s %s"%(self.runTag, self.statLabels[sigMaskNum]), writeDir=self.writeDir, runTag=runTag, plotRings=plotRings, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		elif (detID == 2):
			if self.sigMasksF2 is None:
				self.makeStatMasks(detID=2)
			toPlot = self.sigMasksF2[:,sigMaskNum].reshape(self.det2Setup.rr, self.det2Setup.rc) / ADU
			cxiViewAssembled(toPlot, avgWL, self.det2Setup, filename="%s %s det. 2"%(self.runTag, self.statLabels[sigMaskNum]), writeDir=self.writeDir, runTag=runTag, plotRings=plotRings, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		else:
			print "Invalid detector ID %d given."%(detID)


	def viewUnassembledAveragedSlices(self, fileSelector, detID=2, cmin=0., cmax=2000., imgScale=1.0, ADU=1.):
		avgs = self.avgSlices(fileSelector, detID) / ADU
		if (detID == 1):
			cxiViewUnassembled(avgs, writeDir=self.writeDir, filename=self.runTag, runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		elif (detID == 2):
			cxiViewUnassembled(avgs, writeDir=self.writeDir, filename=self.runTag+"slices_2", runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		else:
			print "Invalid detector ID %d given."%(detID)


	def viewAssembledAveragedSlices(self, fileSelector, detID=2,  cmin=0., cmax=2000., imgScale=1.0, ADU=1., plotRings=False):
		avgs = self.avgSlices(fileSelector, detID) / ADU
		avgWL = N.mean(self.bgWLArr)
		if (detID == 1):
			cxiViewAssembled(avgs, avgWL, self.det1Setup, writeDir=self.writeDir, filename=self.outputFN, runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		elif (detID == 2):
			cxiViewAssembled(avgs, avgWL, self.det2Setup, writeDir=self.writeDir, filename=self.outputFN+"slices_2", runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		else:
			print "Invalid detector ID %d given."%(detID)


	def saveStatsToFile(self, detID=2):
		if not os.path.exists(self.writeDir):
			os.makedirs(self.writeDir)
		outFile = h5py.File(self.outputFN, "w")
		g1 = outFile.create_group('data')
		g1.create_dataset('bgWLArr', data=self.bgWLArr)
		#g1.create_dataset('pulseEnergy_eV', data=self.pulseEnergy_eV)
		
		if (detID == 1) or (detID == -1):
			g1.create_dataset('avgBg', data=self.avgBg)
			g1.create_dataset('rawMaskFlat', data=self.rMF)
			g1.create_dataset('histMax', data=self.bgHistMax)
			g1.create_dataset('histMin', data=self.bgHistMin)
			g1.create_dataset('histBinWidth', data=self.bgHistBinWidth)
			g1.create_dataset('bgHist', data=self.bgHist)
			g1.create_dataset('bgStat', data=self.bgStat)
			g1.create_dataset('rr', data=self.rr)
			g1.create_dataset('rc', data=self.rc)
		
		if (detID == 2) or (detID == -1):
			g1.create_dataset('avgBg2', data=self.avgBg2)
			g1.create_dataset('rawMaskFlat2', data=self.rMF2)
			g1.create_dataset('histMax2', data=self.bgHistMax2)
			g1.create_dataset('histMin2', data=self.bgHistMin2)
			g1.create_dataset('histBinWidth2', data=self.bgHistBinWidth2)
			g1.create_dataset('bgHist2', data=self.bgHist2)
			g1.create_dataset('bgStat2', data=self.bgStat2)
			g1.create_dataset('rr2', data=self.rr2)
			g1.create_dataset('rc2', data=self.rc2)
		
		outFile.close()


	def readStatsFromFile(self, fileLoc, detID=2):
		try:
			inFile = h5py.File(fileLoc, "r")
			self.bgWLArr = inFile['data/bgWLArr'].value
			#self.pulseEnergy_eV = inFile['data/pulseEnergy_eV'].value
			
			if (detID == 1) or (detID == -1):
				try:
					self.avgBg = inFile['data/avgBg'].value
					self.rMF = inFile['data/rawMaskFlat'].value
					self.bgHistMax = inFile['data/histMax'].value
					self.bgHistMin = inFile['data/histMin'].value
					self.bgHistBinWidth = inFile['data/histBinWidth'].value
					self.bgHist = inFile['data/bgHist'].value
					self.bgStat = inFile['data/bgStat'].value
					self.rr = inFile['data/rr'].value
					self.rc = inFile['data/rc'].value
				except:
					print "Unable to read data for detector 1"
			
			if(detID == 2) or (detID == -1):
				try:
					self.avgBg2 = inFile['data/avgBg2'].value
					self.rMF2 = inFile['data/rawMaskFlat2'].value
					self.bgHistMax2 = inFile['data/histMax2'].value
					self.bgHistMin2 = inFile['data/histMin2'].value
					self.bgHistBinWidth2 = inFile['data/histBinWidth2'].value
					self.bgHist2 = inFile['data/bgHist2'].value
					self.bgStat2 = inFile['data/bgStat2'].value
					self.rr2 = inFile['data/rr2'].value
					self.rc2 = inFile['data/rc2'].value
				except:
					print "Unable to read data for detector 2"
		
		except:
			print "Unable to read %s."%(fileLoc)
		inFile.close()


	#TODO:Save bootstrapping results
	def saveBootStrappingResults(self):
	
		pass



class cxidbDataRun(cxidbRun):

	def __init__(self, dataDir, runTag, dataGlob, writeDir, outputFN, det1Setup, det2Setup, sigHistMin=1200., sigHistMax=1800., sigHistBinWidth=1.0, sigHistMin2=1000., sigHistMax2=2000., sigHistBinWidth2=1.0, sigLvls=N.array([0.99,0.95,0.90])):
		cxidbRun.__init__(self, dataDir, runTag, dataGlob, writeDir, outputFN)
		self.det1Setup = det1Setup
		self.det2Setup = det2Setup
		self.wLArr = []
		self.fcounter = 0.
		self.pulseEnergy_eV = None
		self.dataCXIDB = ""
		self.sigLvls = sigLvls
		self.avgWL = self.det1Setup.defaultWavelength
		self.runTag = runTag

		[self.rr, self.rc] = det1Setup.rr, det1Setup.rc
		self.avgSig = N.zeros((self.rr*self.rc))
		self.sigADUF = N.zeros((len(self.sigLvls), self.rr*self.rc))
		self.sigADUFTally = []
		#self.numPixsInSigHist = (det1Setup.rMF!=0).sum()
		#self.sigHistMin = sigHistMin
		#self.sigHistMax = sigHistMax
		#self.sigHistBinWidth = sigHistBinWidth
		#self.sigHist = N.zeros((self.numPixsInSigHist, (self.sigHistMax - self.sigHistMin)/self.sigHistBinWidth + 1))
		self.sigEvents = []
		self.loglikelihood = []
		self.numFrames = 0
		
		[self.rr2, self.rc2] = det2Setup.rr, det2Setup.rc
		self.avgSig2 = N.zeros((self.rr2*self.rc2))
		self.sigADUF2 = N.zeros((len(self.sigLvls), self.rr2*self.rc2))
		self.sigADUFTally2 = []
		#self.numPixsInSigHist2 = (det2Setup.rMF!=0).sum()
		#self.sigHistMin2 = sigHistMin2
		#self.sigHistMax2 = sigHistMax2
		#self.sigHistBinWidth2 = sigHistBinWidth2
		#self.sigHist2 = N.zeros((self.numPixsInSigHist2, (self.sigHistMax2 - self.sigHistMin2)/self.sigHistBinWidth2 + 1))
		self.diag2 = N.zeros((self.rr2*self.rc2))
		self.diag2_b = N.zeros((self.rr2*self.rc2))
		self.diag2_c = []
		self.diag2_d = []
		
		self.sigEvents2 = []
		self.loglikelihood2 = []
		self.numFrames2 = 0
				

	def collectStatistics(self, bg, maxFNum=-1, detID=2, fileSelector=[], histSigSelector=0):

		self.dataCXIDB = G.glob(self.dataGlob)[0]
		print "Data CXIDB file found: %s\n"%(self.dataCXIDB)
		
		f = h5py.File(self.dataCXIDB,"r")
		try:
			self.wLArr = f['LCLS/photon_wavelength_A'].value
		except:
			print "Unable to read photon wavelengths"
		
		if (len(self.wLArr) > 0):
			numBadWLEntries = 0
			for nn,cWL in enumerate(self.wLArr):
				if (cWL == float('Inf') or cWL == -float('Inf') or cWL == float('NaN') or cWL == 0.):
					self.wLArr[nn] = self.det1Setup.defaultWavelength
					numBadWLEntries += 1
			print "%d of %d wavelength entries were undefined and replaced with default"%(numBadWLEntries, len(self.wLArr))
		
		try:
			self.pulseEnergy_eV = f['LCLS/photon_energy_eV'].value
		except:
			print "Unable to read pulse energies"

		if (detID == 1) or (detID == -1):
			try:
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"
				
		if (detID == 2) or (detID ==-1):
			try:
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]
			except:
				print "Unable to read detector_2"


		
		if fileSelector == []:
			if (detID == 1) and ((maxFNum == -1) or (self.numFrames < maxFNum)):
				maxFNum = self.numFrames
			elif (detID == 2) and ((maxFNum == -1) or (self.numFrames2 < maxFNum)):
				maxFNum = self.numFrames2
			elif (detID == -1):
				if (self.numFrames < self.numFrames2):
					if self.numFrames < maxFNum:
						maxFNum = self.numFrames
				else:
					if self.numFrames2 < maxFNum:
						maxFNum = self.numFrames2
			toVisit = range(maxFNum)
		else:
			toVisit = fileSelector
			maxFNum = len(toVisit)

		#HACK:
		self.diag2_c = N.zeros(len(toVisit))
		self.diag2_d = N.zeros(len(toVisit))
		
		self.loglikelihood = N.zeros(maxFNum)
		self.sigEvents = N.zeros((maxFNum, len(self.sigLvls)))
		self.sigADUFTally = N.zeros((maxFNum, len(self.sigLvls)))
		
		self.loglikelihood2 = N.zeros(maxFNum)
		self.sigEvents2 = N.zeros((maxFNum, len(self.sigLvls)))
		self.sigADUFTally2 = N.zeros((maxFNum, len(self.sigLvls)))
		
		startTime = datetime.datetime.now()
		#HACK:
		maskF2 = (bg.rMF2 == -1)
		numMaskF2 = 1.0/(maskF2.sum())
		for self.fcounter in toVisit:
			if (detID == 1) or (detID == -1):
				rawDataF = data[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawDataF, bg.rMF, self.det1Setup.cmNorm)
				processedRawDataF = stat.computeTailSignificanceFromCDF(rawDataF, bg.rMF, bg.bgCDF, bg.bgHistMin, bg.bgHistBinWidth)
				self.loglikelihood[self.fcounter] = stat.computeLogLikelihood(rawDataF, bg.rMF, bg.bgHist, bg.bgHistMin, bg.bgHistBinWidth)
				subtractedDataF = (rawDataF - bg.avgBg)
				
				#toHist = (processedRawDataF>=self.sigLvls[histSigSelector])*subtractedDataF
				#stat.updatePixValHist(toHist, bg.rMF, bg.bgHist, bg.bgHistMin, bg.bgHistMax, bg.bgHistBinWidth)
				self.avgSig += processedRawDataF
				for n,s in enumerate(self.sigLvls):
					tmp = (processedRawDataF >= s)
					self.sigEvents[self.fcounter, n] += tmp.sum()
					tmp *= subtractedDataF
					self.sigADUF[n] += tmp
					self.sigADUFTally[self.fcounter, n] += tmp.sum()

			if (detID == 2) or (detID == -1):
				rawDataF2 = data2[self.fcounter].flatten().astype('float')
				stat.commonModeCorrection(rawDataF2, bg.rMF2, self.det2Setup.cmNorm)
				#HACK: This doesn't always produce positive values...
				subtractedDataF2 = (rawDataF2 - bg.avgBg2)
				
				processedRawDataF2 = stat.computeTailSignificanceFromCDF(rawDataF2, bg.rMF2, bg.bgCDF2, bg.bgHistMin2, bg.bgHistBinWidth2)
				self.loglikelihood2[self.fcounter] = stat.computeLogLikelihood(rawDataF2, bg.rMF2, bg.bgHist2, bg.bgHistMin2, bg.bgHistBinWidth2)
				self.diag2_c[self.fcounter] = ((processedRawDataF2>0.995)*(maskF2)*subtractedDataF2).min()
				self.diag2_d[self.fcounter] = ((processedRawDataF2>0.995)*(1-maskF2)*subtractedDataF2).min()
				#TODO: Figure out why this bg.avgBg2 larger than bg.bgStat[:,0] by ~1.0 in all pixels?
				self.diag2 += processedRawDataF2>self.sigLvls[0]
				self.diag2_b += subtractedDataF2
				#toHist2 = (processedRawDataF2>=self.sigLvls[histSigSelector])*subtractedDataF2
				#stat.updatePixValHist(toHist2, bg.rMF2, bg.bgHist2, bg.bgHistMin2, bg.bgHistMax2, bg.bgHistBinWidth2)
				self.avgSig2 += (processedRawDataF2)
				for n,s in enumerate(self.sigLvls):
					tmp = (processedRawDataF2 >= s)
					self.sigEvents2[self.fcounter, n] += tmp.sum()
					tmp *= subtractedDataF2
					self.sigADUF2[n] += tmp
					self.sigADUFTally2[self.fcounter, n] += tmp.sum()
			if (self.fcounter%10 == 0):
				print "Done with Data %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)
		
		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing Data files took %f seconds each"%((differenceTime.total_seconds())/maxFNum)
		
		if (detID == 1) or (detID == -1):
			self.avgSig /= 1. * maxFNum
			#stat.normalizePixValHist(self.sigHist, 1.E-10)
		if (detID == 2) or (detID == -1):
			self.avgSig2 /= 1. * maxFNum
			#stat.normalizePixValHist(self.sigHist2, 1.E-10)
		self.avgWL = N.array(self.wLArr).mean()
		f.close()

	
	


	def viewPattern(self, bg, fileSelector, detID=2, sigLvlNum=0, runTag="", cmax=2000., cmin=0., imgScale=1.0, ADU=1.):
		self.dataCXIDB = G.glob(self.dataGlob)[0]
		print "Data CXIDB file found: %s\n"%(self.dataCXIDB)
		
		f = h5py.File(self.dataCXIDB,"r")
		
		if (detID == 1) or (detID == -1):
			try:
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"
		
		if (detID == 2) or (detID ==-1):
			try:
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]
			except:
				print "Unable to read detector_2"

		if (detID == 1) or (detID == -1):
			rawDataF = data[fileSelector].flatten().astype('float')
			stat.commonModeCorrection(rawDataF, bg.rMF, self.det1Setup.cmNorm)
			processedRawDataF = stat.computeTailSignificanceFromCDF(rawDataF, bg.rMF, bg.bgCDF, bg.bgHistMin, bg.bgHistBinWidth)
			subtractedDataF = (rawDataF - bg.avgBg)
			toPlot = ((processedRawDataF >= self.sigLvls[sigLvlNum])*subtractedDataF).reshape(self.det1Setup.rr,-1)

		if (detID == 2) or (detID == -1):
			rawDataF2 = data2[fileSelector].flatten().astype('float')
			stat.commonModeCorrection(rawDataF2, bg.rMF2, self.det2Setup.cmNorm)
			processedRawDataF2 = stat.computeTailSignificanceFromCDF(rawDataF2, bg.rMF2, bg.bgCDF2, bg.bgHistMin2, bg.bgHistBinWidth2)
			subtractedDataF2 = (rawDataF2 - bg.avgBg2)
			toPlot = ((processedRawDataF2 >= self.sigLvls[sigLvlNum])*subtractedDataF2).reshape(self.det2Setup.rr,-1)

		f.close()
		
		cxiViewUnassembled(toPlot/ADU, writeDir=self.writeDir, filename="%s_%s_sig_frame_%d"%(self.runTag, self.sigLvls[sigLvlNum], fileSelector), runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		return None


	def viewAveragePattern(self, bg, fileSelector, detID=2, sigLvlNum=0, bgKLDivFilter=0.02, runTag="", cmax=2000., cmin=0., imgScale=1.0, ADU=1.):
		self.dataCXIDB = G.glob(self.dataGlob)[0]
		print "Data CXIDB file found: %s\n"%(self.dataCXIDB)
		
		f = h5py.File(self.dataCXIDB,"r")
		
		if (detID == 1) or (detID == -1):
			try:
				data = f['entry_1/instrument_1/detector_1/data']
				(self.numFrames, rr, rc) = data.shape
				if (data.attrs["numEvents"][0] != self.numFrames):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames = data.attrs["numEvents"][0]
			except:
				print "Unable to read detector_1"
		
		if (detID == 2) or (detID ==-1):
			try:
				data2 = f['entry_1/instrument_1/detector_2/data']
				(self.numFrames2, rr2, rc2) = data2.shape
				if (data2.attrs["numEvents"][0] != self.numFrames2):
					print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
					print "Will default to attributed value instead"
					self.numFrames2 = data2.attrs["numEvents"][0]
			except:
				print "Unable to read detector_2"
		toVisit = N.array(fileSelector)
		toVisit.sort()
				
		if (detID == 1):
			toPlot = N.zeros((self.det1Setup.rr,self.det1Setup.rc))
		elif (detID == 2):
			toPlot = N.zeros((self.det2Setup.rr,self.det2Setup.rc))
		
		for fnum in toVisit:
			if (detID == 1) or (detID == -1):
				rawDataF = data[fnum].flatten().astype('float')
				stat.commonModeCorrection(rawDataF, bg.rMF, self.det1Setup.cmNorm)
				processedRawDataF = stat.computeTailSignificanceFromCDF(rawDataF, bg.rMF, bg.bgCDF, bg.bgHistMin, bg.bgHistBinWidth)
				subtractedDataF = (rawDataF - bg.avgBg)
				toPlot += ((processedRawDataF >= self.sigLvls[sigLvlNum])*subtractedDataF).reshape(self.det1Setup.rr,-1)
			
			if (detID == 2) or (detID == -1):
				rawDataF2 = data2[fnum].flatten().astype('float')
				stat.commonModeCorrection(rawDataF2, bg.rMF2, self.det2Setup.cmNorm)
				processedRawDataF2 = stat.computeTailSignificanceFromCDF(rawDataF2, bg.rMF2, bg.bgCDF2, bg.bgHistMin2, bg.bgHistBinWidth2)
				subtractedDataF2 = (rawDataF2 - bg.avgBg2)
				toPlot += ((processedRawDataF2 >= self.sigLvls[sigLvlNum])*subtractedDataF2).reshape(self.det2Setup.rr,-1)
		
		f.close()
		
		if (detID == 1):
			toPlot *= (bg.bgStat[:,2]>bgKLDivFilter).reshape(self.det1Setup.rr,-1)
		elif (detID == 2):
			toPlot *= (bg.bgStat2[:,2]>bgKLDivFilter).reshape(self.det2Setup.rr,-1)

		
		cxiViewUnassembled(toPlot/ADU, writeDir=self.writeDir, filename="%s_%s_sig_frames"%(self.runTag, self.sigLvls[sigLvlNum]), runTag=self.runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)

		return None


	def plotBackDetectorStats(self, phADU=20., imgScale=1.):
		numCols = len(self.sigLvls)
		if len(self.sigLvls) < 3:
			numCols = 3
		fig = plt.figure(figsize=(imgScale*16.5, imgScale*9.5), dpi=80)
		gs = gridspec.GridSpec(3, numCols)
		gs.update(top=0.95, bottom=0.08, hspace=0.3, wspace=0.25)
				
		ll = self.loglikelihood2.copy()
		ll -= ll.mean()
		for n,s in enumerate(self.sigLvls):
			sd = self.sigEvents2[:,n].copy()
			(y, x) = PL.histogram(sd, bins=100)
			ax = fig.add_subplot(gs[1,n])
			ax.text(.9, .9, "Hist. of %.2f sig pix./frame"%(s), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
			ax.set_xlabel("# sig. pixels")
			ax.set_ylabel("num. frames")
			for label in ax.xaxis.get_ticklabels():
				label.set_rotation(-30)
				label.set_fontsize(10)
			plt.plot(x[1:], y)

			sd = self.sigADUFTally2[:,n].copy()
			(y, x) = PL.histogram(sd, bins=100)
			ax = fig.add_subplot(gs[0,n])
			ax.text(.9, .9, "Hist. of %.2f sig. ADUs/frame"%(s), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
			ax.set_xlabel("Total sig. ADUs")
			ax.set_ylabel("num. frames")
			for label in ax.xaxis.get_ticklabels():
				label.set_rotation(-30)
				label.set_fontsize(10)
			plt.plot(x[1:], y)

			if n < 2:
				ax = fig.add_subplot(gs[2,n])
				ax.text(.92, .92, "Num of sig. pix. vs Log-likelihood"%(s), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
				ax.set_ylabel("Frame Log-likelihood")
				ax.set_xlabel("Num of %.2f sig. pix./frame"%(s))
				for label in ax.xaxis.get_ticklabels():
					label.set_rotation(-30)
					label.set_fontsize(10)
				plt.plot(self.sigEvents2[:,n], ll, 'r.')
		
		(y, x) = PL.histogram(ll, bins=100)
		ax = fig.add_subplot(gs[2,2])
		ax.set_xlabel("Log-likelihood")
		ax.set_ylabel("Num. frames")
		ax.set_yscale('log')
		for label in ax.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		plt.plot(x[1:], y)

		plt.show()



	def viewUnassembledSigMasks(self, sigLvlNum, detID=2, runTag="", cmin=0., cmax=2000., imgScale=1.0, ADU=1.):
		
		avgWL = N.mean(self.wLArr)
		
		if (detID == 1):
			toPlot = self.sigADUF[sigLvlNum].reshape(self.det1Setup.rr, self.det1Setup.rc) / ADU
			cxiViewUnassembled(toPlot, writeDir=self.writeDir, filename="%s_%s"%(self.runTag, self.sigLvls[sigLvlNum]), runTag=runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		elif (detID == 2):
			toPlot = self.sigADUF2[sigLvlNum].reshape(self.det2Setup.rr, self.det2Setup.rc) / ADU
			cxiViewUnassembled(toPlot, writeDir=self.writeDir, filename="%s_%s_sig_det2"%(self.runTag, self.sigLvls[sigLvlNum]), runTag=runTag, cmax=cmax, cmin=cmin, imgScale=imgScale)
		
		else:
			print "Invalid detector ID %d given."%(detID)


	def saveStatsToFile(self, detID=2):
		if not os.path.exists(self.writeDir):
			os.makedirs(self.writeDir)
		
		outFile = h5py.File(self.outputFN, "w")
		g1 = outFile.create_group('data')
		g1.create_dataset('wLArr', data=self.wLArr)
		#g1.create_dataset('pulseEnergy_eV', data=self.pulseEnergy_eV)
		g1.create_dataset('sigLvls', data=self.sigLvls)
		
		if (detID == 1) or (detID == -1):
			try:
				g1.create_dataset('avgSignificance', data=self.avgSig)
				g1.create_dataset('sigADUF', data=self.sigADUF)
				g1.create_dataset('sigADUFTally', data=self.sigADUFTally)
				g1.create_dataset('sigEvents', data=self.sigEvents)
				g1.create_dataset('loglikelihood', data=self.loglikelihood)
				g1.create_dataset('rr', data=self.rr)
				g1.create_dataset('rc', data=self.rc)
				#g1.create_dataset('sigHist', data=self.sigHist)
				#g1.create_dataset('sigHistMin', data=self.sigHistMin)
				#g1.create_dataset('sigHistMax', data=self.sigHistMax)
				#g1.create_dataset('sigHistBinWidth', data=self.sigHistBinWidth)
			except:
				print "Unable to write data for detector 1"
		
		if (detID == 2) or (detID == -1):
			try:
				g1.create_dataset('avgSignificance2', data=self.avgSig2)
				g1.create_dataset('sigADUF2', data=self.sigADUF2)
				g1.create_dataset('sigADUFTally2', data=self.sigADUFTally2)
				g1.create_dataset('sigEvents2', data=self.sigEvents2)		
				g1.create_dataset('loglikelihood2', data=self.loglikelihood2)
				g1.create_dataset('rr2', data=self.rr2)
				g1.create_dataset('rc2', data=self.rc2)
				#g1.create_dataset('sigHist2', data=self.sigHist2)
				#g1.create_dataset('sigHistMin2', data=self.sigHistMin2)
				#g1.create_dataset('sigHistMax2', data=self.sigHistMax2)
				#g1.create_dataset('sigHistBinWidth2', data=self.sigHistBinWidth2)
			except:
				print "Unable to write data for detector 2"
		
		outFile.close()
		

	def readStatsFromFile(self, fileLoc, detID=2):
		
		inFile = h5py.File(fileLoc, "r")
		self.wLArr = inFile['data/wLArr'].value
		#self.pulseEnergy_eV = inFile['data/pulseEnergy_eV'].value

		self.sigLvls = inFile['data/sigLvls'].value

		if (detID == 1) or (detID == -1):
			try:
				self.avgSig = inFile['data/avgSignificance'].value
				self.sigADUF = inFile['data/sigADUF'].value
				self.sigADUFTally = inFile['data/sigADUFTally'].value
				self.sigEvents = inFile['data/sigEvents'].value
				self.loglikelihood = inFile['data/loglikelihood'].value
				self.rr = inFile['data/rr'].value
				self.rc = inFile['data/rc'].value
				#self.sigHist = inFile['data/sigHist'].value
				#self.sigHistMin = inFile['data/sigHistMin'].value
				#self.sigHistMax = inFile['data/sigHistMax'].value
				#self.sigHistBinWidth = inFile['data/sigHistBinWidth'].value
				
			except:
				print "Unable to read data for detector 1"
				
		if (detID == 2) or (detID == -1):
			try:
				self.avgSig2 = inFile['data/avgSignificance2'].value
				self.sigADUF2 = inFile['data/sigADUF2'].value
				self.sigADUFTally2 = inFile['data/sigADUFTally2'].value
				self.sigEvents2 = inFile['data/sigEvents2'].value
				self.loglikelihood2 = inFile['data/loglikelihood2'].value
				self.rr2 = inFile['data/rr2'].value
				self.rc2 = inFile['data/rc2'].value
				#self.sigHist2 = inFile['data/sigHist2'].value
				#self.sigHistMin2 = inFile['data/sigHistMin2'].value
				#self.sigHistMax2 = inFile['data/sigHistMax2'].value
				#self.sigHistBinWidth2 = inFile['data/sigHistBinWidth2'].value
			except:
				print "Unable to read data for detector 2"


		inFile.close()



