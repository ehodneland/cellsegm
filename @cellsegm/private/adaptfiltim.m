%   ADAPTFILTIIM Performs adaptive filtering of image
%
%   ADAPTFILTIIM(I,RAD,D,H) Adaptive filtering on I, RAD defines the filter size 
%   of average filter in relation to the units of the voxel size. D is 
%   the additive value above background value, can be set to 0.02-0.2. H is 
%   the voxelsize  
%   Image should be scaled to [0 1] before entering the algorithm.
%
% 
% 
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
%
function [thim] = adaptfiltim(varargin)

im = varargin{1};
rad = varargin{2};
d = varargin{3};
h = varargin{4};

im = scale(im);

dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;

% radius in voxels
rad = round(rad./h);

% at least one voxel
rad(rad < 1) = 1;

% make threshold image
if dim(3) == 1
    rad = rad(1:2);
    g = fspecial('average',rad);    
else 
    n = prod(rad);
    g = (1/n)*ones(rad(1),rad(2),rad(3));     
end
th = imfilter(im,g,'replicate');

% multiply thresholds above image mean values
th = th + d;

%
% Create ridges
%

% threshold image at different scales to get ridges
thim = gt(im,th);

