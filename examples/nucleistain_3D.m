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
load ../data/nucleistain_3D.mat

% for visualization
plane = 8;

%
% Segmentation by iterative thresholding
%

prm.method = 'thrs';
prm.thrs.th = 1.4;
[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,3,30,'prm',prm);

cellsegm.show(imsegm(:,:,plane),1);
title('Hoechst staining');axis off;
cellsegm.show(cellbw1(:,:,plane),2);
title('Cell segmentation, iterative thresholding, without splitting');axis off;

% Add splitting of cells. Can do this separately to have better control
prm.splitth = 1.8;
% cells above this threshold are split (all cells here)
n = prmout.minvolvox;
h = [0.5 0.5 1.5];
cellbw2 = cellsegm.splitcells(cellbw1,prm.splitth,n,h);
cellsegm.show(cellbw2(:,:,plane),3);
title('Cell segmentation, iterative thresholding, with splitting');axis off;
