#!/usr/bin/env python

import os
import sys
import string
import re
from optparse import OptionParser
import rings

import numpy as N
import h5py

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as PL
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, waveLength, detectorDistance, angAvg=[], filename="plot", writeDir="", runTag="", plotRings=False, cmax=-1, cmin=-1, imgScale=1., adj=[0,0]):
		self.inarr = inarr*(inarr>0)
		for i in range(len(inarr)):
			self.inarr[i] = self.inarr[i][::-1]
		self.filename = filename
		self.write_dir = writeDir
		self.runTag = runTag
		if (cmax == -1):
			self.cmax = self.inarr.max()
		else:
			self.cmax = cmax 

		if (cmin == -1):
			self.cmin = self.inarr.min()
		else:
			self.cmin = cmin 

		self.hist = PL.histogram(self.inarr, bins=300)
		self.imgScale = imgScale
		self.wavelength = waveLength
		self.detectorDistance = detectorDistance
		self.HIceQ = {}
		self.plotRings = plotRings
		self.angAvg = angAvg
		self.adj = adj
		

	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(self.write_dir + self.runTag):
				os.makedirs(self.write_dir + self.runTag)
			pngtag = self.write_dir + self.runTag + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			plt.savefig(pngtag)
		if event.key == 'r':
			colmin, colmax = self.orglims
			plt.clim(colmin, colmax)
			plt.draw()


	def on_click(self, event):
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin, colmax = self.orglims
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			plt.clim(colmin, colmax)
			plt.draw()
				

	def draw_img(self):
	
		print "Right-click on colorbar to set maximum scale."
		print "Left-click on colorbar to set minimum scale."
		print "Center-click on colorbar (or press 'r') to reset color scale."
		print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
		print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
		print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."
		
		fig = plt.figure(figsize=(self.imgScale*9.5, self.imgScale*12.5), dpi=80)
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)

		
		if (len(self.angAvg) > 0):
			gs = gridspec.GridSpec(3,1,height_ratios=[4.5,1,1])
			gs.update(top=0.95, bottom=0.08, hspace=0.25)
			angAvgCanvas = fig.add_subplot(gs[2])
			angAvgCanvas.text(.02, .98, "Angular average", horizontalalignment='left', verticalalignment='top', transform=angAvgCanvas.transAxes)
			angAvgCanvas.set_xlabel("q(inv Angs.)")
			angAvgCanvas.set_ylabel("I(q)")
			if (len(self.angAvg) == 2):
				plt.plot(self.angAvg[0], self.angAvg[1])
			else:
				plt.plot(self.angAvg)
		else:
			gs = gridspec.GridSpec(2,1, height_ratios=[4.5,1])
			gs.update(top=0.95, bottom=0.08, hspace=0.25)

		histCanvas = fig.add_subplot(gs[1])
		histCanvas.text(.02, .98, "Histogram of pixel values", horizontalalignment='left', verticalalignment='top', transform=histCanvas.transAxes)
		histCanvas.set_xlabel("pixel value")
		histCanvas.set_ylabel("num. of pixels")
		histCanvas.set_yscale('log')
		plt.plot(PL.movavg(self.hist[1],2)[2:], self.hist[0][2:])


		canvas = fig.add_subplot(gs[0])
		canvas.set_title(self.filename)
		plt.rc('image',origin='lower')
		self.axes = plt.imshow(self.inarr, vmax = self.cmax, vmin = self.cmin)
		divider = make_axes_locatable(canvas)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		self.colbar = plt.colorbar(self.axes, cax=cax)
		self.orglims = self.axes.get_clim()


		if(self.plotRings):
			for i,j in rings.iceHInvAngQ.iteritems():
				self.HIceQ[i] = rings.get_pix_from_invAngsQ_and_detectorDist(self.runTag, j, self.detectorDistance, wavelengthInAngs=self.wavelength)

			(rTemp,cTemp) = self.inarr.shape
			rTemp += self.adj[0]
			cTemp += self.adj[1]
			#Approximate, center
			for i,j in self.HIceQ.iteritems():
				circ = PL.Circle((rTemp/2, cTemp/2), radius=j)
				circ.set_fill(False)
				circ.set_edgecolor('k')
				canvas.add_patch(circ)

		plt.show()

