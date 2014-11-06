% Pogo - a finite element package to simulate elastic wave propagation on the GPU
% Copyright (C) 2013 Peter Huthwaite
% 
% This file is part of the Matlab post-processing tools made available with Pogo.
% 
% If you find Pogo useful in your academic work, please cite the relevant papers;
% information on our latest papers is available at <http://pogo.peterhuthwaite.com/>.
% 
% Pogo is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Pogo is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Pogo.  If not, see <http://www.gnu.org/licenses/>.
%

clear

for kk = 0
	fileBase = strcat('job-',num2str(kk));
	[ntMeas,dtMeas,nodeNums,nodeDofs,nodeLocs,histTraces] = loadPogoHist(sprintf('%s.pogo-hist',fileBase));

	downSample = 1
	tP = (1:ntMeas)*dtMeas;
	ux = 1e12*histTraces(1:downSample:end,1:3:end).';
	uy = 1e12*histTraces(1:downSample:end,2:3:end).';
	uz = 1e12*histTraces(1:downSample:end,3:3:end).';

	[a,b] = size(ux);

	fid = fopen('nodeNum.txt','w');
	fprintf(fid,'%f\n',nodeNums);
	fclose(fid)

	fid = fopen('nodeLocs.txt','w');
	fprintf(fid,'%f\n',nodeLocs);
	fclose(fid)

	fid = fopen('t.txt','w');
	fprintf(fid,'%.9f\n',tP(downSample:downSample:end));
	fclose(fid);
	
	node_x = nodeLocs(1:9:end);
	node_y = nodeLocs(2:9:end);

	for ii = 1:length(node_x)
		node_theta(ii) = atan2(node_y(ii),node_x(ii));
	end

	for ii = 1:a
		for jj = 1:b
			utheta(ii,jj) = ux(ii,jj)*sin(node_theta(ii))-(uy(ii,jj)*cos(node_theta(ii)));
			ur(ii,jj) = ux(ii,jj)*cos(node_theta(ii))+(uy(ii,jj)*sin(node_theta(ii)));
		end
	end

	yas = 0

	flen = strcat('utheta-',num2str(yas),'.txt')
	fid = fopen(flen,'w');
	for ii = 1:a
		fprintf(fid,'%.5f ',utheta(ii,:));
		fprintf(fid,'\n');
	end
	fclose(fid);

	flen = strcat('ur-',num2str(yas),'.txt')
	fid = fopen(flen,'w');
	for ii = 1:a
		fprintf(fid,'%.5f ',ur(ii,:));
		fprintf(fid,'\n');
	end
	fclose(fid);

	flen = strcat('uz-',num2str(yas),'.txt')
	fid = fopen(flen,'w');
	for ii = 1:a
		fprintf(fid,'%.5f ',uz(ii,:));
		fprintf(fid,'\n');
	end
	fclose(fid);


end
