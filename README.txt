/home/jjd208/VISUALHi everyone,

This folder contains the following scripts:
	- main_pre.py
	- externalScripts.py
	- generate_corrosion.py
	- main_post.py
	- externalScriptsPost.py
	- main_matlab.m

Which I will try to explain a bit below. To generate the file for input to POGO you need to run:

> python2 generate_corrosion.py
> mv hvalues0.dat hvalues.dat
> python2 main_pre.py

To then do the post-processing you run:

> matlab main_matlab.m
> python2 main_post.py


generate_corrosion.py -->
	
	This is a script for generating a gaussian random surface of prescribed correlation length and rms height.
	The script outputs hvalue.dat files which are an input to my main_pre.py script

main_pre.py -->
	
	This script is for creating/meshing a model and then setting up the model with source, monitoring points and absorbing region.
	This will need to be extenively rewritten for your own geometries, but should act as a useful start point.

externalScripts.py -->

	Simply stores functions which are called by main_pre.py. I created this file to stop main_pre.py getting too crowded

main_post.py -->
	
	This is the script for post-processing of the model data. It takes the files from MATLAB and separates out the different modes. It also does 
	the energy calculation. It makes use of data in the folder DISPERSE but please note that this is specific to my geometry and you MUST 
 	create your own DISPERSE directory (email me if you want help with this).

externalScriptsPost.py -->
	
	Simply stores functions which are called by main_post.py

main_matlab.m -->

	A MATLAB script for extracting data from the pogo-hist files and exporting as *.txt files. Note it converts the cartesian displacement data 
	into a circumferential and radial component.

visualise_field.py -->
	
	For those interested in visualising the field output, this file takes the field data and plots it as a contour map. The data must first be 
	converted into txt files using a script similar to main_matlab.m
