import sys
import dispy
import os
import numpy as np
import math as m
import datetime
import random as rdom
import externalScripts as eS
from scipy.interpolate import griddata

#################################################

if __name__ == '__main__':

        print(sys.version)

        r_i = 0.06		# Inner radius
        t = 0.007		# Wall thickness
        r_o = r_i+t
        dz = 3.4		# Model length
        cl = 0.005		# Characteristic element length
 
        filename = 'cylinder_1.geo'
        eS.makeGeo(filename,r_o,r_i,dz,cl)
        os.system('gmsh -3 '+filename+' -rand 1e-9 -smooth 1')

	# Material Properties
	E0 = 204e9
	nu = 0.29
	rho = 7900
	d = cl 	# Element size

	# Source location ...
	rs = 0.06+0.007
	zs = 0.32

	# Monitor location ...
	rm = 0.06 + 0.007
	zm1 = 0.40
	zm2 = 2.30
	monN = 40		# Number of monitoring points at either end
	zspc = 0.02		# Spacing of monitoring points (m)

        # Make the corrosion file and the geo file
	for kk in range(0,1):

		f = open('hvalues.dat','r')
		tmp = f.readlines()
		f.close()

		tmp2 = tmp[0].split(' ')
		h = np.empty((len(tmp),len(tmp2)))

		for ii in range(0,len(tmp)):
		        h[ii,:] = tmp[ii].split(' ')

		# Merge the mesh and the corrosion profile
		[node_num,node_x,node_y,node_z,n1,n2] = eS.getNodeData('cylinder_1.msh')
		[node_conv_x,node_conv_y,node_conv_z] = eS.unwrapCyl2Flat(node_x,node_y,node_z,n1,n2,r_i)
		[top,srf_num] = eS.findTopNodes(n2,n1,node_conv_y,cl)
		[a,b] = np.shape(h)
		
		[Xm,Zm] = eS.findMeshSurfNode(srf_num,top,node_conv_x,node_conv_z)
		[Xc,Zc] = eS.findCorrSurfNode(a,b,cl)

		# So, we have two grids, M and C. We want to interpolate
		#  grid C onto grid M

		for ii in range(0,len(Zc)):
			Zc[ii] = Zc[ii]+1.25		# Offset specifies axial location at which the corrosion starts
		[Xc_mesh,Zc_mesh] = np.meshgrid(Xc,Zc)
		Xc_flat = Xc_mesh.flatten()
		Zc_flat = Zc_mesh.flatten()
		h_flat = h.flatten()
		yi = griddata((Xc_flat,Zc_flat),h_flat,(Xm,Zm),method='nearest')	# Could do more sophsticated interpolation if needed

		# Change node radius
		[node_r,node_theta] = eS.cart2pol(node_x,node_y,node_z)

		node_r_new = node_r

		for ii in range(0,len(top)):
			if yi[ii] > 0:			# Only want to take material away, not add any
		        	node_r_new[top[ii]] = node_r[top[ii]]+(yi[ii])

		[node_x_final,node_y_final] = eS.pol2cart(node_r_new,node_theta,node_z)
		# Write to the mesh file
		eS.writeNewMesh('cylinder_1.msh','cylinder_2.msh',node_x_final,node_y_final,node_z,n1,n2)
		[node_num,node_x,node_y,node_z,n1,n2] = eS.getNodeData('cylinder_2.msh')
		[node_x_new,node_y_new,node_z_new,node_num_new,node_special] = eS.newNodeCoord(node_x,node_y,node_z,node_num,dz)
		[eStart,eEnd] = eS.writeNewElements('cylinder_2.msh',node_special,node_num_new)
		eS.removeCentralNodes('cylinder_2.msh',n2-1,n1,eStart,eEnd-3)

		############################################
		path_home = os.getcwd()
		# Get nodal information
		inputfile = open('cylinder_2.msh','r')
		text = inputfile.readlines()
		inputfile.close()

		f_out = 0
		for lin in text:
		      f_out += 1
		      if '$Node' in lin:
			  n1 = f_out
		      if '$EndNodes' in lin:
			  n2 = f_out

		nomat = np.empty((n2-n1-2,4))
		for i in range(0,n2-n1-2):
			nomat[i,:] = text[n1+1+i].split(' ')

		# Get element information
		f_out = 0
		for lin in text:
		      f_out += 1
		      if '$Elements' in lin:
			  et = f_out
		      if '$EndElements' in lin:
			  e2 = f_out
		f_out = 0
		for i in range(et,e2):
			f_out += 1
			tmp = text[i].split(' ')
			if len(tmp)==9:
			        e1 = f_out+et
			        break

		elmat = np.ndarray(shape=(e2-e1-2,5),dtype=int)
		cnt = 0
		for i in range(0,e2-e1-2):
			temp = text[e1+i].split(' ')
			cnt += 1
			elmat[i,0] = cnt
			elmat[i,1] = int(temp[5])
			elmat[i,2] = int(temp[6])
			elmat[i,3] = int(temp[7])
	       		elmat[i,4] = int(temp[8])

		nonum = len(nomat[:,0])
		elnum = len(elmat[:,0])

		#########################################
		# Check the volume of all the elements to flag any elements which have
		# become very distorted due to disruption of corrosion
		volm = [0]*elnum
		n = np.empty((4,3))
		V1 = [0]*3
		V2 = [0]*3
		V3 = [0]*3
		nodes = [0]*4
		cntv = 0

		for i in range(0,elnum):        # cycle through all the elements
			nodes[0] = elmat[i,1]-1 # Get the nodes for each element
			nodes[1] = elmat[i,2]-1
			nodes[2] = elmat[i,3]-1
	       		nodes[3] = elmat[i,4]-1
			n[0,0] = nomat[nodes[0],1]
			n[0,1] = nomat[nodes[0],2]
			n[0,2] = nomat[nodes[0],3]
			n[1,0] = nomat[nodes[1],1]
			n[1,1] = nomat[nodes[1],2]
			n[1,2] = nomat[nodes[1],3]
			n[2,0] = nomat[nodes[2],1]
			n[2,1] = nomat[nodes[2],2]
			n[2,2] = nomat[nodes[2],3]
			n[3,0] = nomat[nodes[3],1]
			n[3,1] = nomat[nodes[3],2]
			n[3,2] = nomat[nodes[3],3]

		# calculte the three vectors V1,V2,V3
			V1[0] = n[1,0] - n[0,0]         # X
			V2[0] = n[2,0] - n[0,0]
			V3[0] = n[3,0] - n[0,0]
			V1[1] = n[1,1] - n[0,1]         # Y
			V2[1] = n[2,1] - n[0,1]         # Y
			V3[1] = n[3,1] - n[0,1]         # Y
			V1[2] = n[1,2] - n[0,2]         # Z
			V2[2] = n[2,2] - n[0,2]         # Z
			V3[2] = n[3,2] - n[0,2]         # Z

		# perform the volume calculation
			v321 = V3[0]*V2[1]*V1[2]
			v231 = V2[0]*V3[1]*V1[2]
			v312 = V3[0]*V1[1]*V2[2]
			v132 = V1[0]*V3[1]*V2[2]
			v213 = V2[0]*V1[1]*V3[2]
			v123 = V1[0]*V2[1]*V3[2]
			volm[i] = (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123)

		# If the volume is negative, or very small, delete the element
			if (volm[i] < 4e-15):
		       	        del text[e1+i-cntv]
			        cntv += 1

		print('A total of %d elements were removed.....' %cntv)
	
		# Renumber the elements...
		cs = text[e1].split(' ')
		cnt = int(cs[0])
		for i in range(0,(e2-e1-cntv-1)):
			temp = text[e1+i].split(' ')
			temp[0] = str(cnt)
			text[e1+i] = temp[0]+' '+temp[1]+' '+temp[2]+' '+temp[3]+' '+temp[4]+' '+temp[5]+' '+temp[6]+' '+temp[7]+' '+temp[8]
			cnt += 1
		En = int(text[et]) - cntv
		text[et] = str(En)+'\n'
		outputfile = open('cylinder_4.msh','w')
		outputfile.writelines(text)
		outputfile.close()

		############################################
		# Now we start writing all of this to the .inp file
		# wave speed for bulk wave
		CT = m.sqrt(0.5*E0/(rho*(1+nu)))
		CL = m.sqrt(E0*(1-nu)/(rho*(1+nu)*(1-2*nu)))

		fcen = 28e3
		L = CT/fcen
		k = 2*m.pi*fcen/CL
		d = cl 	# Element size
		# Absorbing layers?
		Nlayers = int(m.ceil(2*L/d))

		[dt_step, T, deltat] = [1e-08, 4e-3, 1e-6]	# [time step, total time, how often to output to hist file]

		# number of dof
		n_disp = 3
		n_stress = 6

		eltype = 'C3D4'
		nodenum = 4	# number of nodes per element

		############################################
		fin_node = 'AllNodes.dat'
		fin_elem = 'AllElements.dat'
		finpsig = 'TB_6cycle_28kHz.dat'

		path_home = os.getcwd()
		# Get nodal information
		inputfile = open('cylinder_4.msh','r')
		text = inputfile.readlines()
		inputfile.close()
		f_out = 0
		for lin in text:
		      f_out += 1
		      if '$Node' in lin:
			  n1 = f_out
		      if '$EndNodes' in lin:
			  n2 = f_out

		nomat = np.empty((n2-n1-2,4))
		for i in range(0,n2-n1-2):
			nomat[i,:] = text[n1+1+i].split(' ')

		# Get element information
		f_out = 0
		for lin in text:
		      f_out += 1
		      if '$Elements' in lin:
			  et = f_out
		      if '$EndElements' in lin:
			  e2 = f_out

		f_out = 0
		for i in range(et,e2):
			f_out += 1
			tmp = text[i].split(' ')
			if len(tmp)==9:
				e1 = f_out+et
				break

		elmat = np.ndarray(shape=(e2-e1,5),dtype=int)
		cnt = 1
		for i in range(0,e2-e1):
			tmp = text[e1+i-1].split(' ')
			elmat[i,0] = cnt
			cnt += 1
			elmat[i,1] = int(tmp[5])
			elmat[i,2] = int(tmp[6])
			elmat[i,3] = int(tmp[7])
			elmat[i,4] = int(tmp[8])

		eS.checkForUnconnected(nomat,elmat)

		fnodes = open('AllNodes.dat','w')
		np.savetxt(fnodes,nomat,fmt='%d, %f, %f, %f')
		fnodes.close()

		felements = open('AllElements.dat','w')
		np.savetxt(felements,elmat,fmt='%d',delimiter=',')
		felements.close()

		nonum = np.arange(len(nomat[:,0]))+1
		elnum = elmat[:,0]

		################################################
		# Coordinates of the center of the elements
		zmax = max(nomat[:,3])
		zmin = min(nomat[:,3])

		# coordinates of the center of the elements
		xc=nomat[elmat[:,1:]-1,1].mean(axis=1)
		yc=nomat[elmat[:,1:]-1,2].mean(axis=1)
		zc=nomat[elmat[:,1:]-1,3].mean(axis=1)

		elayer_abs=[]
		for i in range(Nlayers):
		    z0=zmin+ Nlayers*d - (i+1)*d
		    z1=zmax- Nlayers*d + (i+1)*d

		    elayer=elnum[((zc>=z0)&(zc<z0+d)) | ((zc>z1-d)&(zc<=z1))]
		    filename='layer'+str(i+1)+'.dat'	    
	            np.savetxt(filename,elayer,fmt='%9.0f')
		    elayer_abs=np.hstack((elayer_abs,elayer))

		# main domain
		elayerZero=np.setdiff1d(elnum,np.unique(elayer_abs))
		filename='layer0.dat'
		np.savetxt(filename,elayerZero,fmt='%9.0f')

		# Find the source nodes
		source = nonum[ ((np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) <= rs+(d/4)) &
			 (np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) >= rs-(10.0*d)) &
			 (nomat[:,3]>zs-(d/2)) & (nomat[:,3]<zs+(d/2))) ]

		###############################################
		# First monitor location
		monitor1 = nonum[ ((np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) <= rm+(d/4)) &
			 (np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) >= rm-(d/2.0)) &
			 (nomat[:,3]>zm1-(d/2.0))&(nomat[:,3]<zm1+(d/2.0))) ]

		zm1 = zm1+zspc

		for cc in range(0,monN):
			monitort = nonum[ ((np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) <= rm+(d/4)) &
	                         (np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) >= rm-(d/2.0)) &
	                         (nomat[:,3]>zm1-(d/2.0))&(nomat[:,3]<zm1+(d/2.0))) ]
			monitor1 = np.append(monitor1,monitort)
			zm1 = zm1+zspc			

		###############################################
		# Second monitor location
		monitor2 = nonum[ ((np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) <= rm+(d/4)) &
			 (np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) >= rm-(d/2.0)) &
			 (nomat[:,3]>zm2-(d/2.0))&(nomat[:,3]<zm2+(d/2.0))) ]

		zm2 = zm2+zspc

		for cc in range(0,monN):
			monitort = nonum[ ((np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) <= rm+(d/4)) &
	                         (np.sqrt(nomat[:,1]**2 + nomat[:,2]**2) >= rm-(d/2.0)) &
	                         (nomat[:,3]>zm2-(d/2.0))&(nomat[:,3]<zm2+(d/2.0))) ]
			monitor2 = np.append(monitor2,monitort)
			zm2 = zm2+zspc			

		monitor = np.append(monitor1, monitor2)
		np.savetxt('SOURCENODES.dat',source,fmt='%9.0f')
		np.savetxt('MONITORNODES.dat',monitor,fmt='%9.0f')

		X = nomat[source,1]
		Y = nomat[source,2]
		Z = nomat[source,3]
		[R, theta] = eS.cart2pol(X,Y,Z)

		A = np.empty((2,len(X)))

		for ii in range(0,len(X)):
			A[0,ii] = -m.sin(theta[ii])
			A[1,ii] = m.cos(theta[ii])

		fid = open('EXCITE.dat','w')
		fid.write('*Cload,amplitude=Torsional'+'\n')

		for ii in range(0,len(X)):
			fid.write('Part-1-1.'+str(source[ii])+',1,'+str(A[0,ii])+'\n')

		for ii in range(0,len(X)):
			fid.write('Part-1-1.'+str(source[ii])+',2,'+str(A[1,ii])+'\n')

		fid.close()
		################################################
		f = open('job-0.inp','w')
		f.write('*Heading'+'\n');
		f.write('SI units (kg, m, s, N)'+'\n');
		f.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO'+'\n');

		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*Part, name=Part-1'+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*INCLUDE, INPUT=ELSET4LAYERS.DAT'+'\n');

		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*end part'+'\n');
		f.write('**'+'\n');

		f.write('**------------------------------------------------------------'+'\n');
		f.write('** ASSEMBLY'+'\n')
		f.write('**------------------------------------------------------------'+'\n')

		f.write('*Assembly, name=Assembly'+'\n')
		f.write('**'+'\n')

		f.write('*Instance, name=Part-1-1, part=Part-1'+'\n')
		f.write('**'+'\n')
		f.write('*Node, input='+fin_node+'\n')
		f.write('**'+'\n')
		f.write('**'+'\n');
		f.write('*Element, type='+eltype+', input='+fin_elem+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*End Instance'+'\n');


		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*nset,nset=monitor,instance=Part-1-1,input=MONITORNODES.dat'+'\n');
		f.write('*nset,nset=SOURCE,instance=Part-1-1,input=SOURCENODES.dat'+'\n');
		f.write('**'+'\n');

		f.write('*End Assembly'+'\n');

		f.write('**'+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('**------------------------------------------------------------'+'\n')
		f.write('**MATERIALS'+'\n');
		f.write('**------------------------------------------------------------'+'\n')

		f.write('**'+'\n');
		f.write('*INCLUDE, INPUT=MATERIALS.DAT'+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');


		f.write('**------------------------------------------------------------'+'\n')
		f.write('**Amplitude definition'+'\n');
		f.write('**------------------------------------------------------------'+'\n')

		f.write('**'+'\n');
		f.write('*AMPLITUDE,TIME=TOTAL TIME,NAME=Torsional, INPUT='+finpsig+',DEFINITION=TABULAR'+'\n');
		f.write('**'+'\n');


		f.write('**------------------------------------------------------------'+'\n')
		f.write('**Step definition'+'\n');
		f.write('**------------------------------------------------------------'+'\n')

		f.write('**'+'\n');
		f.write('*Step, name=ApplyBC, nlgeom=NO'+'\n'); 
		f.write('*Dynamic, Explicit,direct user control'+'\n');
		f.write(str(dt_step)+', '+str(T)+'\n');
		f.write('*Bulk Viscosity'+'\n');
		f.write('0,0'+'\n');

		f.write('**------------------------------------------------------------'+'\n')
		f.write('**Apply Load'+'\n');
		f.write('**------------------------------------------------------------'+'\n')

		f.write('*INCLUDE, INPUT=EXCITE.dat'+'\n');

		f.write('**'+'\n');

		f.write('**'+'\n');

		f.write('**'+'\n');

		f.write('**------------------------------------------------------------'+'\n')
		f.write('**Output Requests'+'\n');
		f.write('**------------------------------------------------------------'+'\n')
		f.write('**'+'\n');
		f.write('*Restart, write, number interval=1, time marks=NO'+'\n');
		f.write('**'+'\n');

		f.write('**------------------------------------------------------------'+'\n')
		f.write('** FIELD output'+'\n');
		f.write('**------------------------------------------------------------'+'\n')
		f.write('**'+'\n');
		f.write('*Output,field, number interval=40'+'\n');
		f.write('*Node Output'+'\n');
		f.write('U'+'\n');

		f.write('**------------------------------------------------------------'+'\n')
		f.write('** HISTORY output'+'\n');
		f.write('**------------------------------------------------------------'+'\n')
		f.write('**'+'\n');
		f.write('*Output,history, frequency='+str(int(round(deltat/dt_step)))+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*Node Output, nset=monitor'+'\n');
		f.write('U1,U2,U3'+'\n');
		f.write('**'+'\n');
		f.write('**'+'\n');

		f.write('*End Step'+'\n');

		f.close();

		###############################################
		f=open('MATERIALS.DAT','w');

		f.write('*parameter'+'\n')
		f.write('K=3e5'+'\n');

		for ii in range(Nlayers):

		    jj=((ii+1.0)/(Nlayers*1.0))**3;

		    f.write('alpha'+str(ii)+'=K*'+str(jj)+'\n');

		f.write('**'+'\n');
		f.write('**'+'\n');

		# abaqus material property definition commands
		E_decay = 0.01;
		alphamax = -m.log(0.01)/(k*L)
		E = np.arange(Nlayers)*0
		Lx = np.arange(0,Nlayers*d,d)

		f.write('*MATERIAL,NAME=expmaterial'+'\n');
		f.write('*ELASTIC, TYPE=isotropic'+'\n');
		f.write(str(E0)+','+str(nu)+'\n');
		f.write('*DENSITY'+'\n');
		f.write(str(rho)+'\n');


		for i in range(Nlayers):
		    jj = ((i+1.0)/(Nlayers*1.0))**3
		    alphaE = alphamax*jj

		    E = E0*m.exp(-alphaE*k*Lx[i])    
		    f.write('*MATERIAL,NAME=abs_mat'+str(i+1)+'\n');
		    f.write('*ELASTIC, TYPE=isotropic'+'\n');
		    f.write(str(E)+','+str(nu)+'\n');
		    f.write('*DENSITY'+'\n');
		    f.write(str(rho)+'\n');
		    f.write('*DAMPING, ALPHA=<alpha'+str(i)+'>'+'\n');

		f.close()

		###############################################
		f=open('ELSET4LAYERS.DAT','w')

		f.write('** layer0 is the real domain'+'\n')
		for c8 in range(Nlayers+1):
		    f.write('*ELSET, ELSET=layer'+str(c8)+',input='+'layer'+str(c8)+'.dat'+'\n')
		#*SOLID section definition for abaqus

		Nlayers=Nlayers	
		for i in range (Nlayers+1):
		    if i==0:
		        # for real domain
		        f.write('*SOLID SECTION,ELSET=layer'+str(i)+',MATERIAL=expmaterial'+'\n')
		    else:
		        f.write('*SOLID SECTION,ELSET=layer'+str(i)+',MATERIAL=abs_mat'+str(i)+'\n')
		f.close()
