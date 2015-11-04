clear all
close all
%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

% load the data
load ../data/surfstain_3D.mat
imsegm = imsegm(:,:,19);

% Smoothing
prm.smoothim.method = 'dirced';

% Ridge filtering
prm.filterridges = 1;

% Segmentation
prm.classifycells.convexarea = 0.50;
prm.classifycells.convexperim = 0.45;
[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurf(imsegm,20,100,'prm',prm);

cellsegm.show(imsegm,1);title('Raw image');axis off;
cellsegm.show(minima,2);title('Markers');axis off;
cellsegm.show(wat,3);title('Watershed image');axis off;
cellsegm.show(cellbw,4);title('Cell segmentation');axis off;
