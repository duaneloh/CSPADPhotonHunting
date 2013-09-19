import numpy as N
import h5py
import sys
import os
import re

pixSize=110.E-6
nominalWavelengthInAngs=1.693

iceHInvAngQ={'100':1.611, '002':1.717, '101':1.848, '102':2.353, '110':2.793, '103':3.035, '200':3.222, '112':3.272, '201':3.324}
#iceHInvAngQ={'Al111': 2.687, 'Al200':3.103, 'Al220':4.388, 'Al311':5.145}
#iceHInvAngQ={'Decane001':2*N.pi/13.53, 'Decane002': 2*N.pi/6.716, 'Decane010':2*N.pi/4.575, 'Decane01-1':2*N.pi/4.370, 'Decane100':2*N.pi/4.010, 'Decane01-2':2*N.pi/3.833, 'Decane10-1':2*N.pi/3.709, 'Decane1-10':2*N.pi/3.506, 'Decane004':2*N.pi/3.347}

def get_pix_from_invAngsQ_and_detectorDist(run_tag, invAngsQ, detDist, wavelengthInAngs=nominalWavelengthInAngs):
	temp = 2*N.arcsin(0.5*invAngsQ*wavelengthInAngs/(2*N.pi))
	return detDist*N.tan(temp)/pixSize
