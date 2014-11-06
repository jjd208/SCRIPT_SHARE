# A script for generating a gaussian random surface
# Created by Jacob Dobson
# Created on 26th June 2014
# Last modified 03 November 2014

import sys
import dispy
import os
import csv
import numpy as np
import math as m
import datetime
import random as rdom

#################################################
# Set the parameters
r_i = 0.06			# Inner radius of pipe
LenX = 2*m.pi*r_i
LenZ = 1.0			# Length of corrosion section (m)
ds = 0.003			# Spacing between points on corrosion surface (m)
rms = 1.0e-03			# Surface RMS
correlationLength = 0.05	# Surface correlation length

for kk in range(0,1):		# Each loop generates a new surface

        WghtCorrLen = (correlationLength/(m.sqrt(2)))
        M = int(m.floor((4*WghtCorrLen)/(ds*m.sqrt(2)))+1)
        N = M
        w = np.empty((2*M+1,2*N+1))
        ci = 0
        for ii in range(-M,M):
                ci += 1
                cj = 0
                for jj in range(-N,N):
                        cj += 1
                        w[ci,cj]=m.exp(-1*(((np.sqrt(((ii*ds)**2)+((jj*ds)**2)))**2)/(WghtCorrLen**2)))

	# Simple (imperfect) error catching for cases when w isnan
        [a,b] = np.shape(w)
        for ii in range(a):
                for jj in range(b):
                        if m.isnan(w[ii,jj])==True:
                                w[ii,jj] = 0

        w = w/np.sum(w)
        mm = int(m.floor((round((LenX/ds)*10)/10))+1)
        nn = int(m.floor((round((LenZ/ds)*10)/10))+1)
	NumelSurf=mm*nn
        NumelW=((2*M)+1)*((2*N)+1)
        Numelvi= int((mm+2*M)*(nn+2*N))
        v = []
        while len(v) < Numelvi:
                s = rdom.gauss(0,1)
                v.append(s)
        v = np.reshape(v,(mm+2*M,nn+2*N))
        h = np.empty((mm,nn))
        ci = -1
        for ii in range(M,mm+M):
                ci += 1
                cj = -1
	        for jj in range(N,nn+N):
        	        cj += 1
			vlist = v[ii-M:ii+M+1,jj-N:jj+N+1]
                        h[ci,cj] = np.sum(np.multiply(vlist,w))

	# Then, correct the rms
	[a,b] = np.shape(h)
	hs = np.empty((a,b))
	for ii in range(0,a):
		for jj in range(0,b):
			hs[ii,jj] = h[ii,jj]**2

	hrms = np.sqrt(np.mean(hs))

	for ii in range(0,a):
		for jj in range(0,b):
			h[ii,jj] = h[ii,jj]*(rms/hrms)

	# Reduce the tails of the surface to avoid extreme values
	benchm = np.mean(np.mean(h))
	for ii in range(0,a):
		for jj in range(0,b):
			if h[ii,jj] < benchm-(3*rms):
				h[ii,jj] = benchm-(3*rms)
			elif h[ii,jj] > benchm+(3*rms):
				h[ii,jj] = benchm+(3*rms)

        np.savetxt("hvalues"+str(kk)+".dat",h,fmt="%1.6f")


