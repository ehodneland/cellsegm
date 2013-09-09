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
load ../data/surfstain_3D.mat

% No smoothing to save time
prm.smoothim.method = 'none';

% Ridge filtering
prm.filterridges = 1;

% specifying cell quantities
prm.classifycells.convexarea = 0.5;
prm.classifycells.convexperim = 0.40;
prm.classifycells.intincell = 1.30;

% minima level around the middle of cells
prm.getminima.level = 0.45;

% segmentation
[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurf(imsegm,20,100,'prm',prm);
