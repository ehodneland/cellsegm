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

clear all;
close all;

% segmentation channel
prm.segmch = 1;

% automated method to find minima
prm.segmsurf.getminima.method = 'automated';

% adaptive thresholding for finding minima
prm.segmsurf.getminima.automated.adth = 0.01;

% filter ridges
prm.segmsurf.filterridges = 0;

% no smoothing to save time
% prm.segmsurfsmoothimmethod = 'none';
prm.smoothim.method = 'dirced';

% start segmentation at plane 1
prm.segmstart = 1;

% use gpu if Jacket is installed, otherwise set to 0
prm.gpu = 0;

% classification thresholds
prm.segmsurf.classifycells.convexperim = 0.15;
prm.segmsurf.classifycells.intincell = 60;
prm.segmsurf.classifycells.intborder = 120;

% if there is no backbround we must define the background level
prm.segmsurf.classifycells.meanintbck = 1;

% use the strongest signal plane for finding minima
% prm.segmsurf.getminima.level = 'strong';
prm.segmsurf.getminima.automated.level = 0.1;

% run cell segmentation with no parameter file
cellsegm.cellsegmentation('../data/tissue',3,4,1,Inf,0.1,1.7,'','prm',prm);
