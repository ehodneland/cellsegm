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

clear all
close all

% load the data
load ../data/surfstain_and_nucleus_3D.mat

% No smoothing to save time
prm.smoothim.method = 'none';

% No ridge filtering
prm.filterridges = 0;

% No illumination
prm.illum = 0;

% Lower threshold for nucleus markers
prm.getminima.nucleus.segmct.thrs.th = 0.5;

% split markers
prm.getminima.nucleus.segmct.split = 1;
prm.getminima.nucleus.segmct.splitth = 1;

% Subtract the nucleus channel from the surface staining to reduce the
% cross talk effect. 
imsegm1 = imsegm;
imsegm = smooth3(imsegm) - smooth3(imnucl);

[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurf(imsegm,3,100,'imnucleus',imnucl,'prm',prm);

plane = 6;
cellsegm.show(imsegm(:,:,plane),1);title('Surface stain');axis off;
cellsegm.show(imnucl(:,:,plane),2);title('Nucleus stain');axis off;
cellsegm.show(minima(:,:,plane),3);title('Markers');axis off;
cellsegm.show(wat(:,:,plane),4);title('Watershed image');axis off;
cellsegm.show(cellbw(:,:,plane),5);title('Cell segmentation');axis off;
