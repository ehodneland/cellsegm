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
load ../data/surfstain_and_manual_3D.mat

% no ridge filtering
prm.filterridges = 0;

% directional coherence enhancing diffusion
prm.smoothim.method = 'dirced';

% method for minima
prm.getminima.method = 'manual';

% method for classification
prm.classifycells.method = 'minimacell';

% segmentation
[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurf(imsegm,5,50,...
    'prm',prm,'minima',minima,'minimacell',minimacell);

plane = 5;
show(imsegm(:,:,plane),1);
title('Raw image');axis off;
show(imsegmout(:,:,plane),2);
title('Smoothed raw image');axis off;
show(maxprojimage(minima),3);
title('All markers (maximum projection)');axis off;
show(maxprojimage(minimacell),4);
title('Markers cells (maximum projection)');axis off;
show(wat(:,:,plane),5);
title('Watershed image');axis off;
show(cellbw(:,:,plane),6);
title('Cell segmentation');axis off;
