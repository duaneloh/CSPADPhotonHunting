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


class cxiViewRaw():

	def __init__(self, ex, writeDir):
		self.dDist = ex.dDist
		self.writeDir = writeDir
	
	def draw(self, inData, inWL, runTag, inFn, plotRings=False, cmax=2000., cmin=0., imgScale=1.):
		currImg = IM.img_class(inData, inWL, self.dDist, writeDir=self.writeDir, filename=inFn, runTag=runTag, plotRings=plotRings, cmax=cmax, cmin=cmin, imgScale=imgScale)
		currImg.draw_img()



class cxiViewAssemble():

	def __init__(self, ex, writeDir, angAvgMax=1000.):
		self.x = ex.x.astype('float')
		self.y = ex.y.astype('float')
		self.pixSize = ex.pixSize
		self.radInt = ex.radInt
		self.qRaw = ex.qRaw
		self.dDist = ex.dDist
		self.angAvgMax = angAvgMax
		self.writeDir = writeDir
	
	def draw(self, inData, inWL, runTag, inFn, plotRings=True, cmax=2000., cmin=0., imgScale=1.):
		d = geo.interpolateGeometry(self.x, self.y, inData.astype('float'), self.pixSize)
		aa = geo.angAve2D(inData.astype('double'), self.radInt)
		qArr = (self.qRaw / inWL).astype('double')
		aq = geo.angAve2D(qArr, self.radInt)
		currImg = IM.img_class(d, inWL, self.dDist, angAvg=[aq[(aa > 0)*(aa < self.angAvgMax)],aa[(aa > 0)*(aa < self.angAvgMax)]], writeDir=self.writeDir, filename=inFn, runTag=runTag, plotRings=plotRings, cmax=cmax, cmin=cmin, imgScale=imgScale)
		currImg.draw_img()



class cxiExpSetup(object):
	def __init__(self, rr, rc, ar, ac, pixSize, dDist, defaultWavelength, rr2=370, rc2=388, geomFileLoc=None, maskFileLoc=None):
		[self.rr, self.rc] = [rr, rc]
		[self.ar, self.ac] = [ar, ac]
		[self.rr2, self.rc2] = [rr2, rc2]
		[self.pixSize, self.dDist, self.defaultWavelength] = [pixSize, dDist, defaultWavelength]
		[self.xleftPad, self.yleftPad] = [3., 3.]
		[self.x, self.y, self.z] = [None, None, None]
		[self.rad, self.radInt, self.qRaw] = [None, None, None]
		
		self.rM = N.ones((self.rr, self.rc))
		self.rMF = N.ones(self.rr*self.rc)
		
		self.rM2 = N.ones((self.rr2, self.rc2))
		self.rMF2 = N.ones(self.rr2*self.rc2)
		
		self.hasMask = False
		self.hasGeometry = False

		self.maskFileLoc = maskFileLoc
		if self.maskFileLoc is not None:
			try:
				print "Reading mask %s:" % self.maskFileLoc
				mask = h5py.File(self.maskFileLoc, "r")
				self.rM = (mask['data/data'].value).astype('int')
				self.rM2 = (mask['data/mask2'].value).astype('int')
				mask.close()
				self.rMF = self.rM.flatten()
				self.rMF2 = self.rM2.flatten()
				self.hasMask = True
				print ".............. succeeded!\n"
			except:
				print ".............. failed!\n"
		
		self.geomFileLoc = geomFileLoc
		if self.geomFileLoc is not None:
			try:
				print "Reading geometry %s:" % self.geomFileLoc
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
				self.hasGeometry = True
				if self.hasMask:
					self.radInt = (self.radInt*self.rM - (1-self.rM)).astype('int')
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

	def __init__(self, dataDir, runTag, dataGlob, writeDir, outputFN, expSetup, bgHistMin=1500., bgHistMax=2000., bgHistBinWidth=1., bgHistMin2=0., bgHistMax2=2000., bgHistBinWidth2=1.):
		cxidbRun.__init__(self, dataDir, runTag, dataGlob, writeDir, outputFN)
		self.ex = expSetup
		self.bgCfXIDB = ""
		self.bgWLArr = []
		self.pulseEnergy_eV = None
		self.fcounter = 0.
		self.runTag = runTag

		self.rr = self.ex.rr
		self.rc = self.ex.rc
		self.rMF = expSetup.rMF
		self.numPixsInBgHist = int(N.sum(self.rMF))
		self.bgHistMin = bgHistMin
		self.bgHistMax = bgHistMax
		self.bgHist = N.zeros((self.numPixsInBgHist, self.bgHistMax - self.bgHistMin + 1)).astype('float')
		self.bgHistBinWidth = bgHistBinWidth
		self.avgBg = N.zeros((self.rr, self.rc))
		self.bgCDF = []
		self.bgKLGaussian = []
		self.sigMasks = None
		self.numFrames = 0

		self.rr2 = self.ex.rr2
		self.rc2 = self.ex.rc2
		self.rMF2 = expSetup.rMF2.astype('int')
		self.numPixsInBgHist2 = int(N.sum(self.rMF2))
		self.bgHistMin2 = bgHistMin2
		self.bgHistMax2 = bgHistMax2
		self.bgHist2 = N.zeros((self.numPixsInBgHist2, self.bgHistMax2 - self.bgHistMin2 + 1)).astype('float')
		self.bgHistBinWidth2 = bgHistBinWidth2
		self.avgBg2 = N.zeros((self.rr2, self.rc2))
		self.bgCDF2 = []
		self.bgKLGaussian2 = []
		self.sigMasks2 = None
		self.numFrames2 = 0			
			
		
	def collectStatistics(self, maxFNum=-1, detID=-1):
	
		self.bgCXIDB = G.glob(self.dataGlob)[0]
		print "Background CXIDB file found: %s\n"%(self.bgCXIDB)
		
		f = h5py.File(self.bgCXIDB,"r")
		try:
			self.bgWLArr = f['LCLS/photon_wavelength_A'].value
		except:
			print "Unable to read photon wavelengths"
			
		if (len(self.bgWLArr) > 0):
			for nn,cWL in enumerate(self.bgWLArr):
				if (cWL == float('Inf') or cWL == -float('Inf') or cWL == float('NaN') or cWL == 0.):
					self.bgWLArr[nn] = self.ex.defaultWavelength						
		try:
			self.pulseEnergy_eV = f['LCLS/photon_energy_eV'].value
		except:
			print "Unable to read pulse energies"

		try:
			data = f['entry_1/instrument_1/detector_1/data']
			(self.numFrames, rr, rc) = data.shape
			if (data.attrs["numEvents"][0] != self.numFrames):
				print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames, data.attrs["numEvents"][0])
				print "Will default to attributed value instead"
				self.numFrames = data.attrs["numEvents"][0]
		except:
			print "Unable to read detector_1"

		if (rr != self.rr) or (rc != self.rc):
			raise ValueError("Shape of data in detector1 in %s is different from input experimental geometry."%(self.bgCXIDB))

		try:
			data2 = f['entry_1/instrument_1/detector_2/data']
			(self.numFrames2, rr2, rc2) = data2.shape
			if (data2.attrs["numEvents"][0] != self.numFrames2):
				print "Data shape (%05d) different from attributed number of frames (%05d)" %(self.numFrames2, data2.attrs["numEvents"][0])
				print "Will default to attributed value instead"
				self.numFrames2 = data2.attrs["numEvents"][0]
		except:
			print "Unable to read detector_2"

			if (rr2 != self.rr2) or (rc2 != self.rc2):
				raise ValueError("Shape of data (%d, %d) in detector2 in %s is different from expected (%d, %d)."%(rr2, rc2, self.bgCXIDB, self.rr2, self.rc2))

		if (maxFNum == -1) or (maxFNum > self.numFrames) or (maxFNum > self.numFrames2):
			if (self.numFrames < self.numFrames2):
				maxFNum = self.numFrames
			else:
				maxFNum = self.numFrames2


		startTime = datetime.datetime.now()
		for self.fcounter in range(maxFNum):
			if (detID == 1) or (detID == -1):
				rawData = data[self.fcounter]
				self.avgBg += rawData
				rawBgF = (rawData.flatten()).astype('int')
				stat.binDetectorIntADUs(rawBgF, self.bgHistBinWidth)
				stat.updatePixValHist(rawBgF, self.rMF, self.bgHist, self.bgHistMin, self.bgHistMax)
			if (detID == 2) or (detID == -1):
				rawData2 = data2[self.fcounter]
				rawBgF2 = (rawData2.flatten()).astype('int')
				self.avgBg2 += rawData2
				stat.updatePixValHist(rawBgF2, self.rMF2, self.bgHist2, self.bgHistMin2, self.bgHistMax2)
				print "Done with BG %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)


		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing BG files took %f seconds each"%((differenceTime.total_seconds())/self.fcounter)

		print "Normalizing statistics....."
		self.avgBg /= self.fcounter
		self.avgBg2 /= self.fcounter
		self.avgBgWL = N.array(self.bgWLArr).mean()
		stat.normalizePixValHist(self.bgHist, 0.5/self.numPixsInBgHist)
		stat.normalizePixValHist(self.bgHist2, 0.5/self.numPixsInBgHist2)
		print "..done!"
		f.close()


	def computeCDF(self, detID=-1):
		if (detID == 1) or (detID == -1):
			print "Computing c.d.f. for data from detector 1."
			self.bgCDF = stat.computeCDFFromValHist(self.bgHist)
			print "Done!"
		if (detID == 2) or (detID == -1):
			print "Computing c.d.f. for data from detector 2."
			self.bgCDF2 = stat.computeCDFFromValHist(self.bgHist2)
			print "Done!"
					
		
	def computeGaussianDiv(self, detID=-1):
		if (detID == 1) or (detID == -1):
			print "Computing GaussianDiv for detector 1."
			self.bgKLGaussian = stat.computeKLDivFromValHist(self.bgHist, self.bgHistMin, self.bgHistBinWidth)
			print "Done!"
	
		if (detID == 2) or (detID == -1):
			print "Computing GaussianDiv for detector 2."
			self.bgKLGaussian2 = stat.computeKLDivFromValHist(self.bgHist2, self.bgHistMin2, self.bgHistBinWidth2)
			print "Done!"


	def saveStatsToFile(self, detID=-1):
		if not os.path.exists(self.writeDir):
			os.makedirs(self.writeDir)
		outFile = h5py.File(self.outputFN, "w")
		g1 = outFile.create_group('data')
		g1.create_dataset('bgWLArr', data=self.bgWLArr)
		g1.create_dataset('pulseEnergy_eV', data=self.pulseEnergy_eV)

		if (detID == -1) or (detID == 1):
			try:
				g1.create_dataset('avgBg', data=self.avgBg)
				g1.create_dataset('rawMaskFlat', data=self.rMF)
				g1.create_dataset('histMax', data=self.bgHistMax)
				g1.create_dataset('histMin', data=self.bgHistMin)
				g1.create_dataset('histBinWidth', data=self.bgHistBinWidth)
				g1.create_dataset('bgHist', data=self.bgHist)
				g1.create_dataset('bgKLGaussian', data=self.bgKLGaussian)
				g1.create_dataset('rr', data=self.rr)
				g1.create_dataset('rc', data=self.rc)
			except:
				print "Unable to write data for detector 1"

		if (detID == -1) or (detID == 2):
			try:
				g1.create_dataset('avgBg2', data=self.avgBg2)
				g1.create_dataset('rawMaskFlat2', data=self.rMF2)
				g1.create_dataset('histMax2', data=self.bgHistMax2)
				g1.create_dataset('histMin2', data=self.bgHistMin2)
				g1.create_dataset('histBinWidth2', data=self.bgHistBinWidth2)
				g1.create_dataset('bgHist2', data=self.bgHist2)
				g1.create_dataset('bgKLGaussian2', data=self.bgKLGaussian2)
				g1.create_dataset('rr2', data=self.rr2)
				g1.create_dataset('rc2', data=self.rc2)
			except:
				print "Unable to write data for detector 2"
		
		outFile.close()


	def readStatsFromFile(self, fileLoc, detID=-1):
	
		try:
			inFile = h5py.File(fileLoc, "r")
			self.bgWLArr = inFile['data/bgWLArr'].value
			self.pulseEnergy_eV = inFile['data/pulseEnergy_eV'].value

			if (detID == -1) or (detID == 1):
				try:
					self.avgBg = inFile['data/avgBg'].value
					self.rMF = inFile['data/rawMaskFlat'].value
					self.bgHistMax = inFile['data/histMax'].value
					self.bgHistMin = inFile['data/histMin'].value
					self.bgHistBinWidth = inFile['data/histBinWidth'].value
					self.bgHist = inFile['data/bgHist'].value
					self.bgKLGaussian = inFile['data/bgKLGaussian'].value
					self.rr = inFile['data/rr'].value
					self.rc = inFile['data/rc'].value
				except:
					print "Unable to read data for detector 1"

			if(detID == -1) or (detID == 2):
				try:
					self.avgBg2 = inFile['data/avgBg2'].value
					self.rMF2 = inFile['data/rawMaskFlat2'].value
					self.bgHistMax2 = inFile['data/histMax2'].value
					self.bgHistMin2 = inFile['data/histMin2'].value
					self.bgHistBinWidth2 = inFile['data/histBinWidth2'].value
					self.bgHist2 = inFile['data/bgHist2'].value
					self.bgKLGaussian2 = inFile['data/bgKLGaussian2'].value
					self.rr2 = inFile['data/rr2'].value
					self.rc2 = inFile['data/rc2'].value
				except:
					print "Unable to read data for detector 2"

		except:
			print "Unable to read %s."%(fileLoc)
		inFile.close()


	def plotSinglePixelHistogram(self, pixNum, detID=1):
		fig = plt.figure()
		gs = gridspec.GridSpec(4,4, height_ratios=[1,1,1,1])
		gs.update(top=0.9, bottom=0.1, hspace=0.25)
		ax = fig.add_subplot(gs[:,:])
		ax.set_xlabel("ADU")
		ax.set_ylabel("Probability(ADU)")

		if (detID == 1):
			ax.set_title("Distribution for pixel %d:\navg:%.4f, var:%.4f, KLDivGaussian:%.4f"%(pixNum, self.bgKLGaussian[pixNum,0], self.bgKLGaussian[pixNum,1], self.bgKLGaussian[pixNum,2]))
			plt.plot(self.bgHistMin+self.bgHistBinWidth*N.arange(self.bgHist.shape[1]), self.bgHist[pixNum])

		if (detID == 2):
			ax.set_title("Distribution for pixel %d:\navg:%.4f, var:%.4f, KLDivGaussian:%.4f"%(pixNum, self.bgKLGaussian2[pixNum,0], self.bgKLGaussian2[pixNum,1], self.bgKLGaussian2[pixNum,2]))
			plt.plot(self.bgHistMin2+self.bgHistBinWidth2*N.arange(self.bgHist2.shape[1]), self.bgHist2[pixNum])

		plt.show()


	def plotPhotonCountingByKLDiv(self, kmin, kmax, kselect=2, fPlotGamma=0.2, imgScale=1., aduLow=15, aduHigh=30, detID=1):
		fig = plt.figure(figsize=(imgScale*15.5, imgScale*9.5), dpi=80)
		gs = gridspec.GridSpec(4,3)
		gs.update(top=0.95, bottom=0.08, hspace=0.25, wspace=0.25)
							
		if (detID == 1):
			if kmin < 0:
				raise ValueError("kmin, the minimum ranked pixel plotted, has to be greater than or equal zero.")
			if kmax > self.bgHist.shape[0]:
				raise ValueError("kmax, the maximum ranked pixel plotted, has to be less than %d"%self.bgHist.shape[0])
			kl = self.bgKLGaussian[:,kselect].copy()
			ordering = self.bgKLGaussian[:,kselect].argsort()
			selectedPixHist = self.bgHist[ordering[kmin:kmax]]
			xx = self.bgHistMin+self.bgHistBinWidth*N.arange(self.bgHist.shape[1])
			
		elif (detID == 2):
			if kmin < 0:
				raise ValueError("kmin, the minimum ranked pixel plotted, has to be greater than or equal zero.")
			if kmax > self.bgHist2.shape[0]:
				raise ValueError("kmax, the maximum ranked pixel plotted, has to be less than %d"%self.bgHist2.shape[0])
			kl = self.bgKLGaussian2[:,kselect].copy()
			ordering = self.bgKLGaussian2[:,kselect].argsort()
			selectedPixHist = self.bgHist2[ordering[kmin:kmax]]
			xx = self.bgHistMin2+self.bgHistBinWidth2*N.arange(self.bgHist2.shape[1])

		selectedPixHistF = N.array([N.abs(N.fft.ifftn(N.abs(N.fft.fftn(item))**2)) for item in selectedPixHist])[:,:N.floor(0.5*selectedPixHist.shape[1])]
		sortingArray = N.array([N.max(item[aduLow:aduHigh]) for item in selectedPixHistF])
		sortingArrayArg = sortingArray.argsort()
		selectedPixHist = selectedPixHist[sortingArrayArg]
		selectedPixHistF = selectedPixHistF[sortingArrayArg]		
		
		kl.sort()
		klInMask = fig.add_subplot(gs[:,0])
		klInMask.set_xlabel("pixel (sorted by ADU)")
		klInMask.set_ylabel("KL Div. type %d"%kselect)
		klInMask.set_title("KL Div. type %d"%kselect)
		plt.plot([kmin,kmin],[0,2], "k-", lw=2)
		plt.plot([kmax,kmax],[0,2], "k-", lw=2)
		for label in klInMask.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		plt.plot(kl)
		
		pltAspect=selectedPixHist.shape[1]/(1.0*selectedPixHist.shape[0])
		firstLine = N.floor(0.25*selectedPixHist.shape[0])
		secondLine = N.floor(0.75*selectedPixHist.shape[0])

		plot1 = fig.add_subplot(gs[:2,1])
		plot1.set_xlabel("ADU histogram")
		plot1.set_ylabel("pixel")
		plot1.set_title("Stacked histograms for pixels\n ranked %d--%d of %d"%(kmin, kmax, len(ordering)))
		plt.plot([0, selectedPixHist.shape[1]],[firstLine, firstLine], "w-", lw=2)
		plt.plot([0, selectedPixHist.shape[1]],[secondLine, secondLine], "w-", lw=2)
		plt.imshow(selectedPixHist**fPlotGamma, aspect=pltAspect)
		
		plotf = fig.add_subplot(gs[2:,1])
		plotf.set_xlabel("ADU histogram")
		plotf.set_title("Stacked Autocorr. for pixels\nranked %d--%d of %d"%(kmin, kmax, len(ordering)))
		plt.plot([0, selectedPixHist.shape[1]],[firstLine, firstLine], "w-", lw=2)
		plt.plot([0, selectedPixHist.shape[1]],[secondLine, secondLine], "w-", lw=2)
		plt.imshow(selectedPixHistF**fPlotGamma, aspect=0.5*pltAspect)

		ax1 = fig.add_subplot(gs[0,2])
		ax1.set_xlabel("ADU")
		ax1.set_ylabel("Probability(ADU)")
		plt.plot(xx, selectedPixHist[firstLine])

		ax2 = fig.add_subplot(gs[1,2])
		ax2.set_xlabel("ADU")
		ax2.set_ylabel("Probability(ADU)")
		plt.plot(xx, selectedPixHist[secondLine])

		ax3 = fig.add_subplot(gs[2,2])
		ax3.set_xlabel("ADU")
		ax3.set_ylabel("Autocorr. of P(ADU)")
		plt.plot(selectedPixHistF[firstLine])
		
		ax4 = fig.add_subplot(gs[3,2])
		ax4.set_xlabel("ADU")
		ax4.set_ylabel("Autocorr. of P(ADU)")
		plt.plot(selectedPixHistF[secondLine])

		plt.show()

		
	def makeStatMasks(self, detID=-1):
	
		if (detID == 1) or (detID == -1):
			self.sigMasks = stat.makeStatMasks(self.bgKLGaussian, self.rMF)
			
		if (detID == 2) or (detID == -1):
			self.sigMasks2 = stat.makeStatMasks(self.bgKLGaussian2, self.rMF2)



class cxidbDataRun(cxidbRun):

	def __init__(self, dataDir, runTag, dataGlob, writeDir, outputFN, expSetup, sigLvls=N.array([0.99,0.95,0.90])):
		cxidbRun.__init__(self, dataDir, runTag, dataGlob, writeDir, outputFN)
		self.ex = expSetup
		self.wLArr = []
		self.fcounter = 0.
		self.pulseEnergy_eV = None
		self.dataCXIDB = ""
		self.sigLvls = sigLvls
		self.avgWL = self.ex.defaultWavelength
		self.runTag = runTag

		[self.rr, self.rc] = expSetup.rr, expSetup.rc
		self.avgSig = N.zeros((self.rr*self.rc))
		self.sigADU = N.zeros((len(self.sigLvls), self.rr*self.rc))
		self.sigADUTally = []
		self.sigEvents = []
		self.loglikelihood = []
		self.numFrames = 0
		
		[self.rr2, self.rc2] = expSetup.rr2, expSetup.rc2
		self.avgSig2 = N.zeros((self.rr2*self.rc2))
		self.sigADU2 = N.zeros((len(self.sigLvls), self.rr2*self.rc2))
		self.sigADUTally2 = []
		self.sigEvents2 = []
		self.loglikelihood2 = []
		self.numFrames2 = 0
				

	def collectStatistics(self, bg, maxFNum=-1, detID=-1):

		self.dataCXIDB = G.glob(self.dataGlob)[0]
		print "Data CXIDB file found: %s\n"%(self.dataCXIDB)
		
		f = h5py.File(self.dataCXIDB,"r")
		try:
			self.wLArr = f['LCLS/photon_wavelength_A'].value
		except:
			print "Unable to read photon wavelengths"
		if (len(self.wLArr) > 0):
			for nn,cWL in enumerate(self.wLArr):
				if (cWL == float('Inf') or cWL == -float('Inf') or cWL == float('NaN') or cWL == 0.):
					self.wLArr[nn] = self.ex.defaultWavelength
		
		try:
			self.pulseEnergy_eV = f['LCLS/photon_energy_eV'].value
		except:
			print "Unable to read pulse energies"

		try:
			data = f['entry_1/instrument_1/detector_1/data']
		except:
			print "Unable to read detector_1"
		
		try:
			data2 = f['entry_1/instrument_1/detector_2/data']
		except:
			print "Unable to read detector_2"
		
		(self.numFrames, rr, rc) = data.shape
		if (rr != self.rr) or (rc != self.rc):
			raise ValueError("Shape of data in detector1 in %s is different from input experimental geometry."%(self.dataCXIDB))

		(self.numFrames2, rr2, rc2) = data2.shape
		if (rr2 != self.rr2) or (rc2 != self.rc2):
			raise ValueError("Shape of data (%d, %d) in detector2 in %s is different from expected (%d, %d)."%(rr2, rc2, self.dataCXIDB, self.rr2, self.rc2))

		
		if (maxFNum == -1) or (maxFNum > self.numFrames) or (maxFNum > self.numFrames2):
			if (self.numFrames < self.numFrames2):
				maxFNum = self.numFrames
			else:
				maxFNum = self.numFrames2

		self.loglikelihood = N.zeros(maxFNum)
		self.sigEvents = N.zeros((maxFNum, len(self.sigLvls)))
		self.sigADUTally = N.zeros((maxFNum, len(self.sigLvls)))
		
		self.loglikelihood2 = N.zeros(maxFNum)
		self.sigEvents2 = N.zeros((maxFNum, len(self.sigLvls)))
		self.sigADUTally2 = N.zeros((maxFNum, len(self.sigLvls)))
		
		startTime = datetime.datetime.now()
		for self.fcounter in range(maxFNum):
			if (detID == 1) or ((detID == -1)):
				rawDataF = data[self.fcounter].flatten().astype('int')
				stat.binDetectorIntADUs(rawDataF, bg.bgHistBinWidth)
				processedRawDataF = stat.computeTailSignificanceFromCDF(rawDataF, bg.rMF, bg.bgCDF, bg.bgHistMin, bg.bgHistBinWidth)
				self.loglikelihood[self.fcounter] = stat.computeLogLikelihood(rawDataF, bg.rMF, bg.bgHist, bg.bgHistMin, bg.bgHistBinWidth)
				subtractedDataF = (rawDataF - bg.avgBg.flatten())
							
				self.avgSig += processedRawDataF
				for n,s in enumerate(self.sigLvls):
					tmp = (processedRawDataF >= s)
					self.sigEvents[self.fcounter, n] += tmp.sum()
					tmp *= subtractedDataF
					self.sigADU[n] += tmp
					self.sigADUTally[self.fcounter, n] += tmp.sum()

			if (detID == 2) or ((detID == -1)):
				rawDataF2 = data2[self.fcounter].flatten().astype('int')
				stat.binDetectorIntADUs(rawDataF2, bg.bgHistBinWidth2)
				processedRawDataF2 = stat.computeTailSignificanceFromCDF(rawDataF2, bg.rMF2, bg.bgCDF2, bg.bgHistMin2, bg.bgHistBinWidth2)
				self.loglikelihood2[self.fcounter] = stat.computeLogLikelihood(rawDataF2, bg.rMF2, bg.bgHist2, bg.bgHistMin2, bg.bgHistBinWidth2)
				subtractedDataF2 = (rawDataF2 - bg.avgBg2.flatten())
				
				self.avgSig2 += processedRawDataF2
				for n,s in enumerate(self.sigLvls):
					tmp = (processedRawDataF2 >= s)
					self.sigEvents2[self.fcounter, n] += tmp.sum()
					tmp *= subtractedDataF2
					self.sigADU2[n] += tmp
					self.sigADUTally2[self.fcounter, n] += tmp.sum()

			print "Done with Data %s %04d of %04d." % (self.runTag, self.fcounter, maxFNum)
		
		endTime = datetime.datetime.now()
		differenceTime = endTime - startTime
		print "Processing Data files took %f seconds each"%((differenceTime.total_seconds())/self.fcounter)
		
		self.avgSig /= 1. * maxFNum
		self.avgSig2 /= 1. * maxFNum
		self.avgWL = N.array(self.wLArr).mean()
		f.close()


	def plotBackDetectorStats(self, phADU=20., imgScale=1.):
		numCols = len(self.sigLvls)
		fig = plt.figure(figsize=(imgScale*16.5, imgScale*9.5), dpi=80)
		gs = gridspec.GridSpec(3, numCols)
		gs.update(top=0.95, bottom=0.08, hspace=0.3, wspace=0.25)
				
		ll = self.loglikelihood2.copy()
		ll -= ll.mean()
		for n,s in enumerate(self.sigLvls):
			sd = self.sigEvents2[:,n].copy()
			(y, x) = PL.histogram(sd, bins=100)
			ax = fig.add_subplot(gs[0,n])
			ax.set_title("Hist. of %.2f sig pix./frame"%(s))
			ax.set_xlabel("# sig. pixels")
			ax.set_ylabel("num. frames")
			plt.plot(x[1:], y)

			sd = self.sigADUTally2[:,n].copy()
			(y, x) = PL.histogram(sd, bins=100)
			ax = fig.add_subplot(gs[1,n])
			ax.set_title("Hist. of %.2f sig. ADUs/frame"%(s))
			ax.set_xlabel("Total sig. ADUs")
			ax.set_ylabel("num. frames")
			plt.plot(x[1:], y)

			if n < 2:
				ax = fig.add_subplot(gs[2,n])
				ax.set_ylabel("Frame Log-likelihood")
				ax.set_xlabel("Sum of %.2f sig. ADU"%(s))
				for label in ax.xaxis.get_ticklabels():
					label.set_rotation(-30)
					label.set_fontsize(10)
				plt.plot(self.sigADUTally2[:,n], ll, 'r.')
		
		(y, x) = PL.histogram(ll, bins=100)
		ax = fig.add_subplot(gs[2,2])
		ax.set_xlabel("Frame num. (by log-likelihood)")
		ax.set_ylabel("Log-likelihood")
		for label in ax.xaxis.get_ticklabels():
			label.set_rotation(-30)
			label.set_fontsize(10)
		plt.plot(x[1:], y)

		plt.show()


	def saveStatsToFile(self, detID=-1):
		if not os.path.exists(self.writeDir):
			os.makedirs(self.writeDir)
		
		outFile = h5py.File(self.outputFN, "w")
		g1 = outFile.create_group('data')
		g1.create_dataset('wLArr', data=self.wLArr)
		g1.create_dataset('pulseEnergy_eV', data=self.pulseEnergy_eV)
		g1.create_dataset('sigLvls', data=self.sigLvls)
		
		try:
			g1.create_dataset('avgSignificance', data=self.avgSig)
			g1.create_dataset('sigADU', data=self.sigADU)
			g1.create_dataset('sigADUTally', data=self.sigADUTally)
			g1.create_dataset('sigEvents', data=self.sigEvents)
			g1.create_dataset('loglikelihood', data=self.loglikelihood)
			g1.create_dataset('rr', data=self.rr)
			g1.create_dataset('rc', data=self.rc)
		except:
			print "Unable to write data for detector 1"
		
		try:
			g1.create_dataset('avgSignificance2', data=self.avgSig2)
			g1.create_dataset('sigADU2', data=self.sigADU2)
			g1.create_dataset('sigADUTally2', data=self.sigADUTally2)
			g1.create_dataset('sigEvents2', data=self.sigEvents2)		
			g1.create_dataset('loglikelihood2', data=self.loglikelihood2)
			g1.create_dataset('rr2', data=self.rr2)
			g1.create_dataset('rc2', data=self.rc2)
		except:
			print "Unable to write data for detector 2"
		
		outFile.close()
		

	def readStatsFromFile(self, fileLoc, detID=-1):
		
		inFile = h5py.File(fileLoc, "r")
		self.wLArr = inFile['data/wLArr'].value
		self.pulseEnergy_eV = inFile['data/pulseEnergy_eV'].value

		self.sigLvls = inFile['data/sigLvls'].value

		try:
			self.avgSig = inFile['data/avgSignificance'].value
			self.sigADU = inFile['data/sigADU'].value
			self.sigADUTally = inFile['data/sigADUTally'].value
			self.sigEvents = inFile['data/sigEvents'].value
			self.loglikelihood = inFile['data/loglikelihood'].value
			self.rr = inFile['data/rr'].value
			self.rc = inFile['data/rc'].value
		except:
			print "Unable to read data for detector 1"

		try:
			self.avgSig2 = inFile['data/avgSignificance2'].value
			self.sigADU2 = inFile['data/sigADU2'].value
			self.sigADUTally2 = inFile['data/sigADUTally2'].value
			self.sigEvents2 = inFile['data/sigEvents2'].value
			self.loglikelihood2 = inFile['data/loglikelihood2'].value
			self.rr2 = inFile['data/rr2'].value
			self.rc2 = inFile['data/rc2'].value
		except:
			print "Unable to read data for detector 2"


		inFile.close()



