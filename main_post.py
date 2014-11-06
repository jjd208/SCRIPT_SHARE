# A function for calculating where the energy goes in corrosion signals

from scipy import signal
from scipy import ndimage
import numpy as np
import math as m
import cmath
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import externalScriptsPost as pP
from copy import copy, deepcopy

f = open('utheta-0.txt','r')
uthetaf = f.readlines()
f.close()
lux = len(uthetaf[1].split())

r_o = 0.06+0.007
d = 0.003
downSample = 1

num_files = 1
nodez_num = 52		# Number of monitoring points to import
orders = 5		# Number of circumferential orders you want to look at
dx = 0.0225		# Spacing of monitoring points
dt = 2e-06		# Time step of results file

fstart = 25000		# Frequency range to look at
fend = 35000
fpadding = 10		# Padding to apply in frequency and wavenumber
kpadding = 10
order = 2		# Order of interpolation for disperse curves

aRsS = np.zeros((1,orders,3,lux*fpadding))
aTsS = np.zeros((1,orders,3,lux*fpadding))

for jfj in range(0,1):
	# First look in theta direction
	filenum = jfj
	[utr,uti,t] = pP.extractModes(r_o,d,num_files,nodez_num,orders,downSample,0,filenum) 
	[a,b,c] = np.shape(utr)
	F = np.fft.fftfreq(n=c*fpadding,d=dt)
	[aRs,aTs,aSent] = pP.sepModes(utr,uti,dx,dt,fstart,fend,F,kpadding,fpadding,order)
	[aSent,aRs,aTs] = pP.normaliseDisplacements(aRs,aTs,aSent,fstart,fend,F,0,orders,order)
	aSentF = deepcopy(aSent)
	aTsF = deepcopy(aTs)
	aRsF = deepcopy(aRs)

	# Then look at radial
	[utr,uti,t] = pP.extractModes(r_o,d,num_files,nodez_num,orders,downSample,1,filenum) 
	[aRs,aTs,aSent] = pP.sepModes(utr,uti,dx,dt,fstart,fend,F,kpadding,fpadding,order)
	[aSent,aRs,aTs] = pP.normaliseDisplacements(aRs,aTs,aSent,fstart,fend,F,1,orders,order)
	aSentF[1,0,:] = aSent[1,0,:]
	aTsF[1,0,:] = aTs[1,0,:]
	aRsF[1,0,:] = aTs[1,0,:]

	# Then look at axial
	'''[utr,uti,t] = pP.extractModes(r_o,d,num_files,nodez_num,orders,downSample,2,filenum) 
	[aRs,aTs,aSent] = pP.sepModes(utr,uti,dx,dt,fstart,fend,F,kpadding,fpadding,order)
	[aSent,aRs,aTs] = pP.normaliseDisplacements(aRs,aTs,aSent,fstart,fend,F,2,orders,order)
	[aRs, aTs, aSent] = pP.squareSignals(aRs,aTs,aSent)
	aSentF[1,2,:] = aSent[1,2,:]
	aTsF[1,2,:] = aTs[1,2,:]
	aRsF[1,2,:] = aTs[1,2,:]'''
		
	[aRsF, aTsF, aSentF] = pP.squareSignals(aRsF,aTsF,aSentF)
	[aRsF,aTsF] = pP.normaliseBySent(abs(aSentF),abs(aRsF),abs(aTsF))

	aRsS[jfj,:,:,:] = aRsF[:,:,:]
	aTsS[jfj,:,:,:] = aTsF[:,:,:]

aRsM = np.mean(aRsS,axis=0)
aTsM = np.mean(aTsS,axis=0)

pP.plotModes(F,aRsM,aTsM,fstart,fend)
[reflect,transmit,total] = pP.sumSignals(aRsM,aTsM)

plt.figure()
plt.plot(F,total)
plt.xlim(fstart,fend)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Energy balance')
plt.ylim(0.8,1.1)

plt.savefig('total.eps', bbox_inches='tight')
plt.show()
