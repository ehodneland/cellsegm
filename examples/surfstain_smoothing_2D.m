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

load ../data/surfstain_3D.mat;
imsegm = imsegm(150:400,50:300,15);

% smoothing with coherence enhancing diffusion
imsm1 = cellsegm.smoothim(imsegm,'ced');

% smoothing with directional coherence enhancing diffusion
imsm2 = cellsegm.smoothim(imsegm,'dirced');

% smoothing with Gaussian smoothing
imsm3 = cellsegm.smoothim(imsegm,'gaussian');

show(imsegm,1);axis off;axis image;title('Raw image');
show(imsm1,2);axis off;axis image;title('Coherence enhancing diffusion');
show(imsm2,3);axis off;axis image;title('Directional coherence enhancing diffusion');
show(imsm3,4);axis off;axis image;title('Gaussian smoothing');


