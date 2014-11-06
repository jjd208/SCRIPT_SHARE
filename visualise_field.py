from scipy import signal
from scipy import ndimage
import numpy as np
import math as m
import cmath
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.interpolate 

def cart2pol(node_x,node_y,node_z):
        n = len(node_x)
        node_r = [0]*(n)
        node_theta = [0]*(n)

        for i in range(0,n):
                node_r[i] = m.hypot(node_x[i],node_y[i])
                node_theta[i] = np.arctan2(node_y[i],node_x[i])

        return(node_r,node_theta)


def unwrapCyl2Flat(node_x,node_y,node_z,r_i):

        [node_r,node_theta]=cart2pol(node_x,node_y,0)

        # Then convert back to X,Y
	n1 = len(node_r)
        node_conv_x = np.empty(n1)
        node_conv_y = np.empty(n1)
        node_conv_z = np.empty(n1)

        for i in range(0,n1):
                node_conv_x[i] = (2*m.pi*node_r[i])*(node_theta[i]/(2*m.pi))
                node_conv_y[i] = node_r[i] - r_i
                node_conv_z[i] = node_z[i]

        return (node_conv_x,node_conv_y,node_conv_z)

def makeCal():
	filenum = 0
        f = open('utheta-'+str(filenum)+'.txt','r')
        uthetaf = f.readlines()
        f.close()

        lux = len(uthetaf[1].split())

        f = open('nodeLocs.txt','r')
        ntmp = f.readlines()
        f.close()

        nodeLocs = [1]*len(ntmp)
        for ii in range(0,len(ntmp)):
             nodeLocs[ii] = float(ntmp[ii])

        node_x = nodeLocs[::3]
        node_y = nodeLocs[1::3]
        node_z = nodeLocs[2::3]
	r_i = 0.06
        node_theta = np.arctan2(node_y,node_x)
	key = 0
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


	[node_conv_x,node_conv_y,node_conv_z] = unwrapCyl2Flat(node_x,node_y,node_z,r_i)	
	d = 0.02
	# Get those nodes on the outer surface
	cnt = []
	for ii in range(0,len(node_conv_z)):
		if node_conv_y[ii] >= (0.007 - d/2.0) and node_conv_y[ii] <= (0.007 + d/2.0):
			cnt.append(ii)	

	disp_x = np.zeros((len(cnt)))
	disp_z = np.zeros((len(cnt)))
	disp_ut = np.zeros((len(cnt),lux))

	for ii in range(0,len(cnt)):
		disp_x[ii] = node_conv_x[cnt[ii]]
		disp_z[ii] = node_conv_z[cnt[ii]]
		disp_ut[ii] = uthetaT[cnt[ii],:]

	return [disp_x,disp_z,disp_ut,lux]

[disp_x,disp_y,disp_ut,lux] = makeCal()

# Generate grid for plotting
x = np.linspace(-0.2,0.2,100)
y = np.linspace(0,3.0,601)

xi,yi = np.meshgrid(x,y)

xI = xi.flatten()
yI = yi.flatten()

fig = plt.figure()
ax = plt.axes(xlim=(np.min(disp_x), np.max(disp_x)), ylim=(np.min(disp_y), np.max(disp_y)))  
plt.xlabel(r'x')
plt.ylabel(r'y')

# animation function
def animate(i):
    z = scipy.interpolate.griddata( (disp_x,disp_y) , disp_ut[:,5*i], (xI,yI), method='linear')
    zi = np.reshape(z,(601,100))
    v = np.linspace(0,5000,100)
    cont = plt.contourf(x, y, zi, v)
    return cont  

anim = animation.FuncAnimation(fig, animate, frames=40)
anim.save('test.mp4',fps=10)
