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
imsegm = imsegm(:,:,7);

%
% Segmentation by adaptive thresholding
%
% 
prm.method = 'adth';
prm.adth.th = 0.2;
prm.smoothim.method = 'eed';
prm.smoothim.eed.kappa = 0.1;
[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,5,50,'prm',prm);
cellsegm.show(imsegmout,1);title('Raw image');axis off;
cellsegm.show(cellbw1,1);title('Cell segmentation by ADTH');axis off;

% improving the results by splitting of cells
splitth = 2;
plane = 1;
% cells above this threshold are split (all cells here)
n = prmout.minvolvox;
h = [0.5 0.5 1.5];
cellbw2 = cellsegm.splitcells(cellbw1,splitth,n,h);
cellsegm.show(cellbw2,2);title('Cell segmentation by ADTH with splitting');axis off;

%
% Segmentation by iterative thresholding
%

prm.method = 'thrs';
prm.thrsth = 1.2;
prm.split = 0;
prm.smoothim.method = 'eed';
prm.smoothim.eed.kappa = 0.1;
[cellbw3,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,5,50,'prm',prm);
cellsegm.show(cellbw3,3);title('Cell segmentation by THRS');axis off;

% improving the results by splitting of cells
splitth = 1;
plane = 1;
cellbw4 = cellsegm.splitcells(cellbw3,splitth,n,h);
cellsegm.show(cellbw4,4);title('Cell segmentation by THRS with splitting');axis off;