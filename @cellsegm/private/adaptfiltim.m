%   ADAPTFILTIdim(1) Performs adaptive filtering
%
%   ADAPTFILTIdim(1)(Idim(1),RAD,D) Adaptive filtering on Idim(1), RAD defines the filter size 
%   of average filter. D is 
%   the additive value above background value, can be set to 0.02-0.2,  
%   Image should be scaled to [0 1] before entering the algorithm.
%
%   ADAPTFILT(Idim(1),RAD,D,H) includes the pixel size as well.
% 
%   Ex adaptfiltim(im,30,0.05);
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
if nargin == 3
    h = [1 1 1];
elseif nargin == 4
    h = varargin{4};
else
    error('Wrong number of inputs to ADAPTFILT');
end;

im = scale(im);
hx = h(1);
hy = h(2);
hz = h(3);

dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;
% z dimension based on the resolution
radz = floor(rad*(hx/hz));

% make threshold image
if dim(3) == 1
    th = imfilter(im,fspecial('average',rad),'replicate');
else 
    g = 1/(rad*rad*radz)*ones(rad,rad,radz); 
    th = imfilter(im,g,'replicate');
end
    

% multiply thresholds above image mean values
th = th + d;

%
% Create ridges
%

% threshold image at different scales to get ridges
thim = gt(im,th);

