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

% the name of the folders where the data is, as a cell array
name = {'../data/condition1','../data/condition2'};

% start segmentation at plane 3
prmin.segmstart = [3 3;3 3];

% segmentation
cellsegm.cellsegmentation(name,1,2,1,100,2,80,'prmfilenucleus','prm',prmin);

% see the results
a = pwd;
cd(name{1});
cellsegm.viewsegm(1,2,1,2);
cd(a);
