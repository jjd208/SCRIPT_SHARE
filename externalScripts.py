import sys
import dispy
import os
import csv
import numpy as np
import math as m
import datetime
import random as rdom

#########################################################

def makeGeo(filename,r_o,r_i,dz,cl):

        f=open(filename,'w')

        f.write('/* File created:'+'\n')
        f.write(str(datetime.datetime.now())+'*/'+'\n')
        f.write('cl = '+str(cl)+'; \n')
        f.write('r_o = '+str(r_o)+'; \n')
        f.write('r_i = '+str(r_i)+'; \n')
        # Write the points
        f.write('Point(0) = {0,0,0,cl};\n')
        f.write('Point(1) = {r_o,0,0,cl};\n')
        f.write('Point(2) = {0,r_o,0,cl};\n')
        f.write('Point(3) = {-r_o,0,0,cl};\n')
        f.write('Point(4) = {0,-r_o,0,cl};\n')
        f.write('Point(5) = {r_i,0,0,cl};\n')
        f.write('Point(6) = {0,r_i,0,cl};\n')
        f.write('Point(7) = {-r_i,0,0,cl};\n')
        f.write('Point(8) = {0,-r_i,0,cl};\n')

        # Write the lines
        f.write('Circle(9) = {1,0,2};\n')
        f.write('Circle(10) = {2,0,3};\n')
        f.write('Circle(11) = {3,0,4};\n')
        f.write('Circle(12) = {4,0,1};\n')
        f.write('Circle(13) = {5,0,6};\n')
        f.write('Circle(14) = {6,0,7};\n')
        f.write('Circle(15) = {7,0,8};\n')
        f.write('Circle(16) = {8,0,5};\n')

	nmb = np.ceil(dz/cl);
        # Tidy up
        f.write('Line Loop(17) = {9,10,11,12}; \n')
        f.write('Line Loop(18) = {13,14,15,16}; \n')
        f.write('Plane Surface(19) = {17,18}; \n')
        f.write('Extrude {0,0,'+str(dz)+'} {Surface{19}; Layers{'+str(nmb)+'}; }; \n')
        f.write('Mesh.CharacteristicLengthMin = cl;\n')
        f.write('Mesh.Optimize = 1;\n')
        f.close
###################################################

def getNodeData(filename):
        inputfile = open(filename,'r')
        text = inputfile.readlines()
        inputfile.close()
        f_out = 0
        for lin in text:
                f_out += 1
                if '$Node' in lin:
                        n1 = f_out
                if '$EndNodes' in lin:
                        n2 = f_out

        nodes = text[n1+1:n2]
        node_num = np.empty(n2-n1-2)
        node_x = np.empty(n2-n1-2)
        node_y = np.empty(n2-n1-2)
        node_z = np.empty(n2-n1-2)

        for i in range(0,n2-n1-2):
                a = text[n1+1+i].split(' ')
                node_num[i] = float(a[0])
                node_x[i] = float(a[1])
                node_y[i] = float(a[2])
                node_z[i] = float(a[3])
        return (node_num,node_x,node_y,node_z,n1,n2)                

###################################################

def unwrapCyl2Flat(node_x,node_y,node_z,n1,n2,r_i):

        [node_r,node_theta]=cart2pol(node_x,node_y,0)

        # Then convert back to X,Y
        node_conv_x = np.empty(n2-n1-2)
        node_conv_y = np.empty(n2-n1-2)
        node_conv_z = np.empty(n2-n1-2)

        for i in range(0,n2-n1-2):
                node_conv_x[i] = (2*m.pi*node_r[i])*(node_theta[i]/(2*m.pi))
                node_conv_y[i] = node_r[i] - r_i
                node_conv_z[i] = node_z[i]

        return (node_conv_x,node_conv_y,node_conv_z)                

####################################################

def findTopNodes(n2,n1,node_conv_y,cl):
        top = []
        for i in range(0,n2-n1-2):
                if node_conv_y[i] <= cl/4.0 and node_conv_y[i] >= -cl/4.0:
                               top.append(i)
        return (top,len(top))

####################################################

def findMeshSurfNode(srf_num,top,node_conv_x,node_conv_z):
        Xm = [0]*srf_num
        Zm = [0]*srf_num
        for i in range(0,srf_num):
                Xm[i] = node_conv_x[top[i]]
                Zm[i] = node_conv_z[top[i]]
        return (Xm,Zm)

###################################################

def findCorrSurfNode(a,b,ds):
        Xc = np.arange(0,(a*ds),ds)-((a*ds)/2)
        Zc = np.arange(0,(b*ds),ds)
        return (Xc[0:a],Zc[0:b])        

##################################################

def cart2pol(node_x,node_y,node_z):
        n = len(node_x)
        node_r = [0]*(n)
        node_theta = [0]*(n)

        for i in range(0,n):
                node_r[i] = m.hypot(node_x[i],node_y[i])
                node_theta[i] = np.arctan2(node_y[i],node_x[i]) 

        return(node_r,node_theta)

def pol2cart(node_r,node_theta,node_z):
        n = len(node_r)
        node_x = np.empty(n)
        node_y = np.empty(n)

        for i in range(0,n):
                node_x[i] = node_r[i]*m.cos(node_theta[i])
                node_y[i] = node_r[i]*m.sin(node_theta[i])

        return(node_x,node_y)

#################################################

def cSurfHeight(s,node_r,a,b,srf_num,top,h):
        jacob = list(node_r)
        for j in range(0,a):
                       for k in range(0,b):
                        ii = int(s[j,k])
                        jacob[top[ii]] = node_r[top[ii]] + h[j,k]
        return(jacob)

#################################################

def writeNewMesh(input_filename,output_filename,node_x,node_y,node_z,n1,n2):
        inputfile = open(input_filename,'r')
        text = inputfile.readlines()
        inputfile.close()

        for i in range(0, n2-n1-2):
                text[n1+1+i] = str(i+1) +' '+str(node_x[i])+' '+str(node_y[i])+' '+str(node_z[i])+'\n'
        outputfile = open(output_filename,'w')
        outputfile.writelines(text)
        outputfile.close()

def nSmallest(a,n):
        a = a.tolist()
        for i in range(0,n):
                d = a.index(min(a))
                a[d] = 1000
        return min(a)

def writeS(filename,s,a,b):
        outf = open(filename,'w')
        for i in range(0,a):
                for j in range(0,b):
                        outf.write('{}'.format(s[i,j]))
                        outf.write('\t')        
                outf.write('\n')

#################################################
#################################################
def newNodeCoord(node_x,node_y,node_z,node_num,dz):
        endD = node_num[len(node_num)-1]
        dist_z = [1]*len(node_z)

        node_special = []

        for ii in range(0,len(node_z)):
                dist_z[ii] = dz-node_z[ii]
                if dist_z[ii] == 0:
                        node_special.append(node_num[ii])

        node_new_len = len(node_z)-len(node_special)
        node_x_new = [1]*node_new_len
        node_y_new = [1]*node_new_len
        node_z_new = [1]*node_new_len
        node_num_new = [0]*len(node_z)
        cnt = 0
        for ii in range(0,len(node_z)):
                if node_num[ii] in node_special:
                        node_num_new[ii] = node_num[ii]
                else:
                        node_num_new[ii] = cnt + len(node_z)+1
                        node_x_new[cnt] = node_x[ii]
                        node_y_new[cnt] = node_y[ii]
                        node_z_new[cnt] = dz+dist_z[ii]
                        cnt += 1
	return(node_x_new,node_y_new,node_z_new,node_num_new,node_special)

def writeNewNodes(in_file,out_file,node_x,node_y,node_z,n2,n1,numNew):

        inputfile = open(in_file,'r')
        text = inputfile.readlines()
        inputfile.close()
        textlen = len(text)
        text = [0]*(textlen+numNew)

        inputfile = open(in_file,'r')
        textn = inputfile.readlines()
        inputfile.close()

        for ii in range(0,n2-1):
                text[ii] = textn[ii]

        # Correct the number of nodes
        nodeNums = int(text[n1])
        nodeNums2 = nodeNums + numNew
        text[n1] = str(nodeNums2)+'\n'

        for ii in range(0,numNew):      # Insert the new nodes
                text[n2+ii-1] = str(nodeNums+1+ii)+' '+str(node_x[ii])+' '+str(node_y[ii])+' '+str(node_z[ii])+'\n'

        # Then insert the rest
        for ii in range(0,textlen-n2+1):
                text[ii+n2+numNew-1] = textn[n2+ii-1]

        for ii in range(0,len(text)):
                if type(text[ii]) != str:
                        print(ii)
        outputfile = open(out_file,'w')
        outputfile.writelines(text)
        outputfile.close()
        return (nodeNums)

def writeNewElements(file,node_special,node_num_new):
	# nodeNum = number of added nodes
        inputfile = open(file,'r')
        text = inputfile.readlines()
        inputfile.close()
	textlen = len(text)

	# Get element information
	f_out = 0
	for lin in text:
      		f_out += 1
      		if '$Elements' in lin:
          		eStart = f_out
      		if '$EndElements' in lin:
          		eEnd = f_out

	return(eStart,eEnd)

def removeCentralNodes(file,n2,n1,eStart,eEnd):

	inputfile = open(file,'r')
	text = inputfile.readlines()
	inputfile.close()
	
	s = [] # Stores the index of the deleted node

	for ii in range(n1+1,n2):
		tmp = text[ii].split(' ')
		if float(tmp[1])==0 and float(tmp[2])==0:
			s.append(ii)	
	s.reverse()
	print(len(s))

	for jj in range(0,len(s)):
		del text[s[jj]]
	s.reverse()

	for ii in range(0,len(s)):	# Convert from line number
		s[ii] = s[ii]-4		# to node number
	textn = list(text)
	print(s)

	textn[n1] = str(int(text[n1])-len(s))+'\n'

	cnt = 1
	for ii in range(n1+1,n2-len(s)):
		tmp = text[ii].split(' ')
		tmp[0] = cnt
		textn[ii] = str(tmp[0])+' '+tmp[1]+' '+tmp[2]+' '+tmp[3]
		cnt += 1

	# Then adjust the element numbers
	for ii in range(eStart-1,eEnd):
		tmp = text[ii].split(' ')
		if len(tmp)==6:	# Check last two nodes
			tmp[4] = int(tmp[4])
			tmp[5] = int(tmp[5])
			for jj in range(0,len(s)):
				if tmp[4] >= s[jj]:
					tmp[4] -= 1
				if tmp[5] >= s[jj]:
					tmp[5] -= 1
			textn[ii] = tmp[0]+' '+tmp[1]+' '+tmp[2]+' '+tmp[3]+' '+str(tmp[4])+' '+str(tmp[5])+'\n' 			
		elif len(tmp)==7:
			tmp[5] = int(tmp[5])
			tmp[6] = int(tmp[6])
			for jj in range(0,len(s)):
				if tmp[5] >= s[jj]:
					tmp[5] -= 1
				if tmp[6] >= s[jj]:
					tmp[6] -= 1
			textn[ii] = tmp[0]+' '+tmp[1]+' '+tmp[2]+' '+tmp[3]+' '+tmp[4]+' '+str(tmp[5])+' '+str(tmp[6])+'\n' 			

		elif len(tmp)==8:
			tmp[5] = int(tmp[5])
			tmp[6] = int(tmp[6])
			tmp[7] = int(tmp[7])
			for jj in range(0,len(s)):
				if tmp[5] >= s[jj]:
					tmp[5] -= 1
				if tmp[6] >= s[jj]:
					tmp[6] -= 1
				if tmp[7] >= s[jj]:
					tmp[7] -= 1
			textn[ii] = tmp[0]+' '+tmp[1]+' '+tmp[2]+' '+tmp[3]+' '+tmp[4]+' '+str(tmp[5])+' '+str(tmp[6])+' '+str(tmp[7])+'\n' 			

		elif len(tmp)==9:
			tmp[5] = int(tmp[5])
			tmp[6] = int(tmp[6])
			tmp[7] = int(tmp[7])
			tmp[8] = int(tmp[8])
			
			for kk in range(5,9):
				tmp[kk] -= 1
				if tmp[kk] >= s[1]:
					tmp[kk] -= 1

			textn[ii] = tmp[0]+' '+tmp[1]+' '+tmp[2]+' '+tmp[3]+' '+tmp[4]+' '+str(tmp[5])+' '+str(tmp[6])+' '+str(tmp[7])+' '+str(tmp[8])+'\n' 			

        inputfile = open(file,'w')
        text = inputfile.writelines(textn)
        inputfile.close()

def checkForUnconnected(nomat,elmat):
        L1 = set(nomat[:,0])
        L2 = set(elmat[:,1])
        L3 = set(elmat[:,2])
        L4 = set(elmat[:,3])
        L5 = set(elmat[:,4])
        print(L1.difference(L2,L3,L4,L5))
