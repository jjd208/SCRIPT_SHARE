function [ ntMeas, dtMeas, nodeNums, nodeDofs, nodeLocs, histTraces ] = loadPogoHist( fileName )
%loadPogoHist - load history data from Pogo FE
%
% [ ntMeas, dtMeas, nodeNums, nodeDofs, nodeLocs, histTraces ] = loadPogoHist( fileName )
%
%fileName - the file name
%ntMeas - the number of measurement times
%dtMeas - the time spacing between measurement times (s)
%nodeNums - the numbers of the nodes at which measurementns are given
%nodeDofs - the degree of freedom for each measurement
%nodeLocs - location of each node; dimension fast, node number slow
%histTraces - the measurements. Time is fast, node number slow
%
% Written by P. Huthwaite, September 2012

%Pogo - a finite element package to simulate elastic wave propagation on the GPU
%Matlab tools, used to read and post-process the files produced with Pogo.
%Copyright (C) 2013 Peter Huthwaite
%
%If you find Pogo useful in your academic work, please cite the relevant papers;
%information on our latest papers is available at <http://pogo.peterhuthwaite.com/>.
%
%Pogo is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%Pogo is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with Pogo.  If not, see <http://www.gnu.org/licenses/>.

    fid = fopen(fileName,'rb');
    if (fid == -1) 
        disp('File could not be opened.')
        return;
    end

    text = fread(fid, 20, '*char')';
    
    if strcmp(text, '%pogo-hist1.0')
        disp('File is wrong format.')
        return
    end
    
    prec = fread(fid, 1, 'int32');
    if prec ~= 4 && prec ~= 8
        disp('Precision incorrect.')
        return 
    end
    nDims = fread(fid, 1, 'int32');
    nMeas = fread(fid, 1, 'int32');
    ntMeas = fread(fid, 1, 'int32');
    
    if prec == 4
        dtMeas = fread(fid, 1, 'float32');
    else
        dtMeas = fread(fid, 1, 'float64');
    end
    nodeNums = zeros(nMeas,1);
    nodeDofs = zeros(nMeas,1);
    nodeLocs = zeros(nDims,nMeas);
    histTraces = zeros(ntMeas,nMeas);
    
    for cnt = 1:nMeas
        nodeNums(cnt) = fread(fid, 1, 'int32');
        nodeDofs(cnt) = fread(fid, 1, 'int32');
        if prec == 4
            nodeLocs(:,cnt) = fread(fid, nDims, 'float32');
            histTraces(:,cnt) = fread(fid, ntMeas, 'float32');
        else
            nodeLocs(:,cnt) = fread(fid, nDims, 'float64');
            histTraces(:,cnt) = fread(fid, ntMeas, 'float64');
        end
    end
    
    
    fclose(fid);


end

