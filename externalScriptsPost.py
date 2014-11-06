from scipy import signal
from scipy import ndimage
import numpy as np
import math as m
import cmath
import time
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

def extractModes(r_o,d,num_files,nodez_num,orders,downSample,key,filenum):

	f = open('utheta-'+str(filenum)+'.txt','r')
	uthetaf = f.readlines()
	f.close()

	f = open('t.txt','r')
	timef = f.readlines()
	f.close()

	timep = [0]*(len(timef)-1)
	for ii in range(0,len(timef)-1):
		timep[ii] = float(timef[ii])

	lux = len(uthetaf[1].split())

	f = open('nodeLocs.txt','r')
	ntmp = f.readlines()
	f.close()

	nodeLocs = [1]*len(ntmp)
	for ii in range(0,len(ntmp)):
	     nodeLocs[ii] = float(ntmp[ii])

	node_x = nodeLocs[::9]
	node_y = nodeLocs[1::9]
	node_z = nodeLocs[2::9]
	node_theta = np.arctan2(node_y,node_x)
	node_r = np.hypot(node_x,node_y)

	# Load in displacements
	if key == 0:
		f = open('utheta-'+str(filenum)+'.txt','r')
	if key == 1:
		f = open('ur-'+str(filenum)+'.txt','r')
	if key == 2:
		f = open('uz-'+str(filenum)+'.txt','r')
	uthetaf = f.readlines()
	f.close()

	uthetaT = np.zeros((len(uthetaf),lux))

	for ii in range(0,len(uthetaf)):
	     uthetaT[ii,:] = uthetaf[ii].split()

	uzu = np.unique(node_z)
	uzu.sort()
	print(uzu)
	s = np.array(np.zeros((nodez_num,len(node_z)/len(uzu))),dtype=int)

	for jj in range(0,nodez_num):
		cnt = 0
		for ii in range(0,len(node_z)):
			if node_z[ii] == uzu[jj]:
				s[jj,cnt] = int(ii)
				cnt += 1

	[a,b] = np.shape(s)
	uthetaR = np.zeros((nodez_num,b,lux))	
	for jj in range(0,a):
		for kk in range(0,b):
			uthetaR[jj,kk,:] = uthetaT[s[jj,kk],:]

	#######################################################
	# Next, separate into the different cicumferential 
	# Start with the reflected only first
	utr = np.zeros((nodez_num,orders,lux))
	uti = np.zeros((nodez_num,orders,lux))

	print('Entering bottleneck')


	for gg in range(0,nodez_num):
		utr[gg,0,:] = np.sum(uthetaR[gg,:,:],axis=0)

		for cc in range(1,orders):
		       uthetar = np.zeros((b,lux))
		       uthetai = np.zeros((b,lux))
		       mur = np.zeros((b,lux))
		       mui = np.zeros((b,lux))

		       m1 = np.exp(cc*1j*node_theta[s[gg,:]])
		       m1r = m1.real
		       m1i = m1.imag
		       for ii in range(0,lux):
				mur[:,ii] = m1r	
				mui[:,ii] = m1i

       		       uthetar = np.multiply(mur,uthetaR[gg,:,:])
		       uthetai = np.multiply(mui,uthetaR[gg,:,:])
			       
		       utr[gg,cc,:] = np.sum(uthetar, axis=0)    
		       uti[gg,cc,:] = np.sum(uthetai, axis=0)

	print('Exitiing bottleneck')
	return utr, uti, timep

def sepModes(utr,uti,dx,dt,fstart,fend,F,kpadding,fpadding,order):

	##########################
	# Structure of the data:
	# aRs = reflected components
	# aTs = transmitted components
	# rows = different circumferential orders
	# columns = different mode indices
	
	[a,b,c] = np.shape(utr)
	
	aRs = np.array(np.zeros((b,3,c*fpadding)),dtype=complex)
	aTs = np.array(np.zeros((b,3,c*fpadding)),dtype=complex)
	aSent = np.array(np.zeros((b,3,c*fpadding)),dtype=complex)

	for ii in range(0,b):	# For each order
		uR = utr[0:(a/2),ii,:] + 1j*uti[0:(a/2),ii,:]
		uT = utr[(a/2):,ii,:] + 1j*uti[(a/2):,ii,:]

		# Perform 2D-FFT
		aR = np.fft.fft2(uR,s=(a*(kpadding/2),c*fpadding))
		aT = np.fft.fft2(uT,s=(a*(kpadding/2),c*fpadding))
		fy = np.fft.fftfreq(n=a*(kpadding/2),d=dx)
		fx = np.fft.fftfreq(n=c*fpadding,d=dt)		
		stdDev = fy[1]-fy[0]
	
		# Then separate different modes
		if ii == 0:	#T(0,1)

			T01rng = []
			L02rng = []
			for cig in range(0,2):
				if cig == 0:	
					f = open('DISPERSE/T01k.txt','r')
				elif cig == 1:
					f = open('DISPERSE/L02k.txt','r')
				t = f.readlines()
				f.close()

				xi = [0]*(len(t)-2)
				yi = [0]*(len(t)-2)
				for bb in range(0,len(t)-2):
				        av = t[bb].split('\t')
				        xi[bb] = 1e6*float(av[0])
				        yi[bb] = 1e3*(float(av[1]))

				s = InterpolatedUnivariateSpline(xi, yi, k=order)
				kval = s(fx)
	
				for jlc in range(0,c*fpadding):
					if fx[jlc] >fstart and fx[jlc] < fend:
						k = kval[jlc]
						if k < 0:
							k = 0
						rng = []
						for fts in range(0,len(fy)):
							if fy[fts] > k-stdDev and fy[fts] < k+stdDev:
								rng.append(fts)
								break
						if cig == 0 and rng:
							T01rng.append(rng[0])
						elif cig ==0 and not rng:
							T01rng.append(T01rng[-1])
						elif cig == 1 and rng:
							L02rng.append(rng[0])
						elif cig == 1 and not rng:
							L02rng.append(L02rng[-1])			
				# Now divide up the space ....
			cnt = 0
			for vn in range(0,c*fpadding):
				if fx[vn] > fstart and fx[vn] < fend:
					barrierOne = int((T01rng[cnt] + L02rng[cnt]) / 2.0)
					aRs[ii,0,vn] =  (np.sqrt(np.sum((abs(aR[0:a*(kpadding/4),vn]))**2)))
					aSent[ii,0,vn] =  (np.sqrt(np.sum((abs(aR[a*(kpadding/4):-1,vn]))**2)))
					aTs[ii,0,vn] =  np.sqrt(np.sum((abs(aT[a*(kpadding/4):-1,vn]))**2))
					cnt += 1

		elif ii == 1:   #F(1,*)
			F11rng = []
			F12rng = []
			F13rng = []

			for cig in range(0,3):
				if cig == 0:	
					f = open('DISPERSE/F11k.txt','r')
				elif cig == 1:
					f = open('DISPERSE/F12k.txt','r')
				elif cig == 2:
					f = open('DISPERSE/F13k.txt','r')
				t = f.readlines()
				f.close()

				xi = [0]*(len(t)-2)
				yi = [0]*(len(t)-2)
				for bb in range(0,len(t)-2):
				        av = t[bb].split('\t')
				        xi[bb] = 1e6*float(av[0])
				        yi[bb] = 1e3*(float(av[1]))

				s = InterpolatedUnivariateSpline(xi, yi, k=order)
				kval = s(fx)

				for jlc in range(0,c*fpadding):
					if fx[jlc] > fstart and fx[jlc] < fend:
						k = kval[jlc]
						if k < 0:
							k = 0
						rng = []
						for fts in range(0,len(fy)):
							if fy[fts] > k-stdDev and fy[fts] < k+stdDev:
								rng.append(fts)
								break
						if cig == 0 and rng:
							F11rng.append(rng[0])
						elif cig ==0 and not rng:
							F11rng.append(F11rng[-1])
						elif cig == 1 and rng:
							F12rng.append(rng[0])
						elif cig == 1 and not rng:
							F12rng.append(F12rng[-1])			
						elif cig == 2 and rng:										
							F13rng.append(rng[0])
						elif cig == 2 and not rng:										
							F13rng.append(F13rng[-1])
				# Now divide up the space ....
			cnt = 0
			jh = []
			yc = []
			jd = []
			for vn in range(0,c*fpadding):
				if fx[vn] > fstart and fx[vn] < fend:
					barrierTwo = int((F11rng[cnt] + F12rng[cnt]) / 2.0)
					barrierOne = int((F12rng[cnt] + F13rng[cnt]) / 2.0)
					aRs[ii,2,vn] =  (np.sqrt(np.sum((abs(aR[0:barrierOne,vn]))**2)))
					aSent[ii,2,vn] =  (np.sqrt(np.sum((abs(aR[-barrierOne:-1,vn]))**2)))
					aTs[ii,2,vn] =  (np.sqrt(np.sum((abs(aT[-barrierOne:-1,vn]))**2)))
					aRs[ii,1,vn] =  (np.sqrt(np.sum((abs(aR[barrierOne:barrierTwo,vn]))**2)))
					aSent[ii,1,vn] =  (np.sqrt(np.sum((abs(aR[-barrierTwo:-barrierOne,vn]))**2)))
					aTs[ii,1,vn] = (np.sqrt(np.sum((abs(aT[-barrierTwo:-barrierOne,vn]))**2)))
					aRs[ii,0,vn] =  (np.sqrt(np.sum((abs(aR[barrierTwo:a*(kpadding/4),vn]))**2)))
					aSent[ii,0,vn] =  (np.sqrt(np.sum((abs(aR[a*(kpadding/4):-barrierTwo,vn]))**2)))
					aTs[ii,0,vn] =  np.sqrt(np.sum((abs(aT[a*(kpadding/4):-barrierTwo,vn]))**2))
					cnt += 1
			#######################################################################################
		elif ii == 2:   #F(2,*)
	
			F22rng = []
			F23rng = []
			kstart = []

			for cig in range(0,2):
				if cig == 0:	
					f = open('DISPERSE/F22k.txt','r')
				elif cig == 1:
					f = open('DISPERSE/F23k.txt','r')
				t = f.readlines()
				f.close()

				xi = [0]*(len(t)-2)
				yi = [0]*(len(t)-2)
				for bb in range(0,len(t)-2):
				        av = t[bb].split('\t')
				        xi[bb] = 1e6*float(av[0])
				        yi[bb] = 1e3*float(av[1])
				kstart.append(xi[0])
				s = InterpolatedUnivariateSpline(xi, yi, k=order)
				kval = s(fx)

				for jlc in range(0,c*fpadding):
					if fx[jlc] > fstart and fx[jlc] < fend:
						k = kval[jlc]
						rng = []
						for fts in range(0,len(fy)):
							if fy[fts] > k-stdDev and fy[fts] < k+stdDev:
								rng.append(fts)
						if cig == 0 and rng:
							F22rng.append(rng[0])
						elif cig ==0 and not rng:
							F22rng.append(F22rng[-1])
						elif cig == 1 and rng:
							F23rng.append(rng[0])
						elif cig == 1 and not rng:
							F23rng.append(F23rng[-1])			

			# Now divide up the space ....
			cnt = 0
			for vn in range(0,c*fpadding):
				if fx[vn] > fstart and fx[vn] < fend:
					if fx[vn] < kstart[1]:	# No f23
						aRs[ii,1,vn] =  0
						aTs[ii,1,vn] =  0
						aRs[ii,0,vn] =  np.sqrt(np.sum((abs(aR[0:a*(kpadding/4),vn]))**2))
						aTs[ii,0,vn] =  np.sqrt(np.sum((abs(aT[a*(kpadding/4):-1,vn]))**2))
						cnt += 1

					elif fx[vn] < kstart[0]: # No f22 and f23
						aRs[ii,1,vn] =  0
						aTs[ii,1,vn] =  0
						aRs[ii,0,vn] =  0
						aTs[ii,0,vn] =  0
						cnt += 1
					else:
						barrierOne = int((F22rng[cnt] + F23rng[cnt]) / 2.0)
						aRs[ii,1,vn] =  np.sqrt(np.sum((abs(aR[0:barrierOne,vn]))**2))
						aTs[ii,1,vn] =  np.sqrt(np.sum((abs(aT[-barrierOne:-1,vn]))**2))
						aRs[ii,0,vn] =  np.sqrt(np.sum((abs(aR[barrierOne:a*(kpadding/4),vn]))**2))
						aTs[ii,0,vn] =  np.sqrt(np.sum((abs(aT[a*(kpadding/4):-barrierOne,vn]))**2))
						cnt += 1

		if ii == 3:	#F(3,*)
			aRs[ii,0,:] = np.sqrt(np.sum((abs(aR[0:(a*kpadding)/4.0,:]))**2,axis=0))
			aTs[ii,0,:] = np.sqrt(np.sum((abs(aT[(a*kpadding)/4.0:,:]))**2,axis=0))
			for sr in range(0,len(fx)):
				if fx[sr] < 24000:
					aRs[ii,0,sr] = 0
					aTs[ii,0,sr] = 0

		if ii == 4:	#F(4,2)
			aRs[ii,0,:] = np.sqrt(np.sum((abs(aR[0:(a*kpadding)/4.0,:]))**2,axis=0))
			aTs[ii,0,:] = np.sqrt(np.sum((abs(aT[(a*kpadding)/4.0:,:]))**2,axis=0))
			for sr in range(0,len(fx)):
				if fx[sr] < 31500:
					aRs[ii,0,sr] = 0
					aTs[ii,0,sr] = 0

	for ii in range(0,len(fx)):		
		if fx[ii] < fstart or fx[ii] > fend:
			aRs[:,:,ii] = 0
			aTs[:,:,ii] = 0
			aSent[:,:,ii] = 0
	return aRs, aTs, aSent
	
def normaliseDisplacements(aRs,aTs,aSent,fstart,fend,F,key,orders,order):
	
	if orders == 2:
		guy = 5
	elif orders == 3:
		guy = 7
	elif orders == 4:
		guy = 8
	elif orders == 5:
		guy = 9

	[a,b,c] = np.shape(aRs)
	
	for cfj in range(0,guy):
		print(cfj)
		if cfj == 0 and key == 0:
			f = open('DISPERSE/THETA/T01.txt','r')
		if cfj == 1 and key == 0:
			f = open('DISPERSE/THETA/L02.txt','r')
		elif cfj == 2 and key == 0:
			f = open('DISPERSE/THETA/F11.txt','r')
		elif cfj == 3 and key == 0:
			f = open('DISPERSE/THETA/F12.txt','r')
		elif cfj == 4 and key == 0:
			f = open('DISPERSE/THETA/F13.txt','r')
		elif cfj == 5 and key == 0:
			f = open('DISPERSE/THETA/F22.txt','r')
		elif cfj == 6 and key == 0:
			f = open('DISPERSE/THETA/F23.txt','r')
		elif cfj == 7 and key == 0:
			f = open('DISPERSE/THETA/F32.txt','r')
		elif cfj == 8 and key == 0:
			f = open('DISPERSE/THETA/F42.txt','r')
		elif cfj == 0 and key == 1:
			f = open('DISPERSE/RADIAL/T01.txt','r')
		elif cfj == 1 and key == 1:
			f = open('DISPERSE/RADIAL/L02.txt','r')
		elif cfj == 2 and key == 1:
			f = open('DISPERSE/RADIAL/F11.txt','r')
		elif cfj == 3 and key == 1:
			f = open('DISPERSE/RADIAL/F12.txt','r')
		elif cfj == 4 and key == 1:
			f = open('DISPERSE/RADIAL/F13.txt','r')
		elif cfj == 5 and key == 1:
			f = open('DISPERSE/RADIAL/F22.txt','r')
		elif cfj == 6 and key == 1:
			f = open('DISPERSE/RADIAL/F23.txt','r')
		elif cfj == 7 and key == 1:
			f = open('DISPERSE/RADIAL/F32.txt','r')
		elif cfj == 8 and key == 1:
			f = open('DISPERSE/THETA/F42.txt','r')
		elif cfj == 0 and key == 2:
			f = open('DISPERSE/AXIAL/T01.txt','r')
		elif cfj == 1 and key == 2:
			f = open('DISPERSE/AXIAL/L02.txt','r')
		elif cfj == 2 and key == 2:
			f = open('DISPERSE/AXIAL/F11.txt','r')
		elif cfj == 3 and key == 2:
			f = open('DISPERSE/AXIAL/F12.txt','r')
		elif cfj == 4 and key == 2:
			f = open('DISPERSE/AXIAL/F13.txt','r')
		elif cfj == 5 and key == 2:
			f = open('DISPERSE/AXIAL/F22.txt','r')
		elif cfj == 6 and key == 2:
			f = open('DISPERSE/AXIAL/F23.txt','r')
		elif cfj == 7 and key == 2:
			f = open('DISPERSE/AXIAL/F32.txt','r')
		elif cfj == 8 and key == 2:
			f = open('DISPERSE/THETA/F42.txt','r')

		t = f.readlines()
		f.close()

		xi = [0]*(len(t)-2)
		yi = [0]*(len(t)-2)
		for ii in range(0,len(t)-2):
		        a = t[ii].split('\t')
		        xi[ii] = 1e6*float(a[0])
		        yi[ii] = abs(float(a[1]))

		s = InterpolatedUnivariateSpline(xi, yi, k=order)
		mval = s(F)

		# Adjust all T(0,1)
		for ii in range(0,c):
			if cfj == 0:
				aSent[0,0,ii] = aSent[0,0,ii] / mval[ii]
				aRs[0,0,ii] =  aRs[0,0,ii] / mval[ii]
				aTs[0,0,ii] =  aTs[0,0,ii] / mval[ii]
			if cfj == 1:
				aSent[0,1,ii] = aSent[0,1,ii] / mval[ii]
				aRs[0,1,ii] =  aRs[0,1,ii] / mval[ii]
				aTs[0,1,ii] =  aTs[0,1,ii] / mval[ii]
			elif cfj == 2:
				aSent[1,0,ii] =  aSent[1,0,ii] / mval[ii]
				aRs[1,0,ii] =  aRs[1,0,ii] / mval[ii]
				aTs[1,0,ii] =  aTs[1,0,ii] / mval[ii]
			elif cfj == 3:
				aSent[1,1,ii] =  aSent[1,1,ii] / mval[ii]
				aRs[1,1,ii] =  aRs[1,1,ii] / mval[ii]
				aTs[1,1,ii] =  aTs[1,1,ii] / mval[ii]
			elif cfj == 4:
				aSent[1,2,ii] =  aSent[1,2,ii] / mval[ii]
				aRs[1,2,ii] =  aRs[1,2,ii] / mval[ii]
				aTs[1,2,ii] =  aTs[1,2,ii] / mval[ii]
			elif cfj == 5:
				aRs[2,0,ii] =  aRs[2,0,ii] / mval[ii]
				aTs[2,0,ii] =  aTs[2,0,ii] / mval[ii]
			elif cfj == 6:
				aRs[2,1,ii] =  aRs[2,1,ii] / mval[ii]
				aTs[2,1,ii] =  aTs[2,1,ii] / mval[ii]
			elif cfj == 7:
				aRs[3,0,ii] =  aRs[3,0,ii] / mval[ii]
				aTs[3,0,ii] =  aTs[3,0,ii] / mval[ii]
			elif cfj == 8:
				aRs[4,0,ii] =  aRs[4,0,ii] / mval[ii]
				aTs[4,0,ii] =  aTs[4,0,ii] / mval[ii]

	return aSent, aRs, aTs
		
def normaliseBySent(aSent,aRs,aTs):
	
	[a,b,c] = np.shape(aRs)
	
	for ii in range(0,a):
		for jj in range(0,b):
			aRs[ii,jj,:] = aRs[ii,jj,:] / aSent[0,0,:]
			aTs[ii,jj,:] = aTs[ii,jj,:] / aSent[0,0,:]
	return aRs, aTs 

def squareSignals(aRs,aTs,aSent):
	
	[a,b,c] = np.shape(aRs)

	for ii in range(0,a):
		for jj in range(0,b):
			aRs[ii,jj,:] = aRs[ii,jj,:]**2
			aTs[ii,jj,:] = aTs[ii,jj,:]**2
			aSent[ii,jj,:] = aSent[ii,jj,:]**2

	return aRs, aTs, aSent

def sumSignals(aRs,aTs):

	a = len(aRs[0,0,:])
	reflect = np.zeros((a,1))
	transmit = np.zeros((a,1))
	total = np.zeros((a,1))
	for ii in range(0,a):
		reflect[ii] = aRs[0,0,ii] + aRs[1,0,ii] + aRs[1,1,ii] + aRs[1,2,ii] + aRs[2,0,ii] + aRs[2,1,ii] + aRs[3,0,ii] + aRs[0,1,ii] + aRs[4,0,ii]
		transmit[ii] = aTs[0,0,ii] + aTs[1,0,ii] + aTs[1,1,ii] + aTs[1,2,ii] + aTs[2,0,ii] + aTs[2,1,ii] + aTs[3,0,ii] + aTs[0,1,ii] + aTs[4,0,ii]
		total[ii] = reflect[ii] + transmit[ii]

	return reflect,transmit,total

def plotModes(F,aRs,aTs,fstart,fend):

        plt.figure()
	ax = plt.subplot(1,1,1)
        ax.plot(F,aTs[0,0,:],label='T(0,1) transmitted')
        ax.plot(F,aRs[0,0,:],label='T(0,1) reflected')
        ax.plot(F,aTs[0,1,:],label='L(0,2) transmitted')
        ax.plot(F,aRs[0,1,:],label='L(0,2) reflected')
        ax.plot(F,aTs[1,0,:],'--',label='F(1,1) transmitted')
        ax.plot(F,aRs[1,0,:],'--',label='F(1,1) reflected')
        ax.plot(F,aTs[1,1,:],'--',label='F(1,2) transmitted')
	ax.plot(F,aRs[1,1,:],'--',label='F(1,2) reflected')
        ax.plot(F,aTs[1,2,:],'--',label='F(1,3) transmitted')
        ax.plot(F,aRs[1,2,:],'--',label='F(1,3) reflected')
        ax.plot(F,aTs[2,0,:],'.-',label='F(2,2) transmitted')
        ax.plot(F,aRs[2,0,:],'.-',label='F(2,2) reflected')
        ax.plot(F,aTs[2,1,:],'.-',label='F(2,3) transmitted')
        ax.plot(F,aRs[2,1,:],'.-',label='F(2,3) reflected')
        ax.plot(F,aTs[3,0,:],'.-',label='F(3,2) transmitted')
        ax.plot(F,aRs[3,0,:],'.-',label='F(3,2) reflected')
        ax.plot(F,aTs[4,0,:],':',label='F(4,2) transmitted')
        ax.plot(F,aRs[4,0,:],':',label='F(4,2) reflected')
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.xlim(fstart,fend)
#	plt.ylim(0,1.1)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Energy spectrum for different modes')
	plt.savefig('waves.eps', bbox_inches='tight')

def smoothSignals(aSent,aRs,aTs,num_sig,sig_std):
	
	sg = signal.gaussian(num_sig,sig_std)
	[a,b,c] = np.shape(aRs)

	for ii in range(0,a):
		for jj in range(0,b):
			aSent[ii,jj,:] = ndimage.filters.convolve1d(aSent[ii,jj,:],sg/sg.sum())
			aRs[ii,jj,:] = ndimage.filters.convolve1d(aRs[ii,jj,:],sg/sg.sum())
			aTs[ii,jj,:] = ndimage.filters.convolve1d(aTs[ii,jj,:],sg/sg.sum())

	return aSent, aRs, aTs

def divideByGVelocity(aRs,aTs,aSent,fstart,fend,F,key,orders,order):
	
	if orders == 2:
		guy = 5
	elif orders == 3:
		guy = 7
	elif orders == 4:
		guy = 8
	elif orders == 5:
		guy = 9

	[a,b,c] = np.shape(aRs)

	filenames = 'GROUP_VEL/'	

	for cfj in range(0,guy):
		if cfj == 0 and key == 0:
			f = open(filenames+'T01.txt','r')
		if cfj == 1 and key == 0:
			f = open(filenames+'L02.txt','r')
		elif cfj == 2 and key == 0:
			f = open(filenames+'F11.txt','r')
		elif cfj == 3 and key == 0:
			f = open(filenames+'F12.txt','r')
		elif cfj == 4 and key == 0:
			f = open(filenames+'F13.txt','r')
		elif cfj == 5 and key == 0:
			f = open(filenames+'F22.txt','r')
		elif cfj == 6 and key == 0:
			f = open(filenames+'F23.txt','r')
		elif cfj == 7 and key == 0:
			f = open(filenames+'F32.txt','r')
		elif cfj == 8 and key == 0:
			f = open(filenames+'F42.txt','r')

		t = f.readlines()
		f.close()

		xi = [0]*(len(t)-2)
		yi = [0]*(len(t)-2)
		for ii in range(0,len(t)-2):
		        a = t[ii].split('\t')
		        xi[ii] = 1e6*float(a[0])
		        yi[ii] = abs(float(a[1]))

		s = InterpolatedUnivariateSpline(xi, yi, k=order)
		mval = s(F)

		# Adjust all T(0,1)
		for ii in range(0,c):
			if cfj == 0:
				aSent[0,0,ii] = aSent[0,0,ii] / mval[ii]
				aRs[0,0,ii] =  aRs[0,0,ii] / mval[ii]
				aTs[0,0,ii] =  aTs[0,0,ii] / mval[ii]
			if cfj == 1:
				aSent[0,1,ii] = aSent[0,1,ii] / mval[ii]
				aRs[0,1,ii] =  aRs[0,1,ii] / mval[ii]
				aTs[0,1,ii] =  aTs[0,1,ii] / mval[ii]
			elif cfj == 2:
				aSent[1,0,ii] =  aSent[1,0,ii] / mval[ii]
				aRs[1,0,ii] =  aRs[1,0,ii] / mval[ii]
				aTs[1,0,ii] =  aTs[1,0,ii] / mval[ii]
			elif cfj == 3:
				aSent[1,1,ii] =  aSent[1,1,ii] / mval[ii]
				aRs[1,1,ii] =  aRs[1,1,ii] / mval[ii]
				aTs[1,1,ii] =  aTs[1,1,ii] / mval[ii]
			elif cfj == 4:
				aSent[1,2,ii] =  aSent[1,2,ii] / mval[ii]
				aRs[1,2,ii] =  aRs[1,2,ii] / mval[ii]
				aTs[1,2,ii] =  aTs[1,2,ii] / mval[ii]
			elif cfj == 5:
				aRs[2,0,ii] =  aRs[2,0,ii] / mval[ii]
				aTs[2,0,ii] =  aTs[2,0,ii] / mval[ii]
			elif cfj == 6:
				aRs[2,1,ii] =  aRs[2,1,ii] / mval[ii]
				aTs[2,1,ii] =  aTs[2,1,ii] / mval[ii]
			elif cfj == 7:
				aRs[3,0,ii] =  aRs[3,0,ii] / mval[ii]
				aTs[3,0,ii] =  aTs[3,0,ii] / mval[ii]
			elif cfj == 8:
				aRs[4,0,ii] =  aRs[4,0,ii] / mval[ii]
				aTs[4,0,ii] =  aTs[4,0,ii] / mval[ii]

	return aSent, aRs, aTs
