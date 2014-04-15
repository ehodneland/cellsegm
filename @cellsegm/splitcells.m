function [cellbw] = splitcells(cellbw,splitth,splitvolvox,h)
% SPLITCELLS Splitting segmented binary image of cells 
%
%   CELLBW = SPLITCELLS(CELLBW,SPLITTH,MINVOLVOX,H) splitting the cells 
%   in CELLBW and the splitting threshold in SPLITTH with the voxel volume 
%   MINVOLVOX. SPLITTH is the second argument in
%   IMEXTENDEDMAX. H is the voxel size. 
%   Returning the splitted cells in the binary image CELLBW.
%
%   SPLITCELLS is based upon the distance function and local maxima for 
%   splitting, similar to 
%
%   "Combining intensity, edge and shape information for 2d and 3D
%   segmentaton fo cell nuclei in tissue sections, Wahlby, Journal of
%   microscopy, Vol 215, pp 67-76, 2004". 
%
%
%   See also segmct
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
vis = 0;

msg = ['This is ' upper(mfilename) ' splitting objects'];
disp(msg);

dim = size(cellbw);
ndim = numel(dim);
if ndim == 2
    dim = [dim 1];
end;

if ndim == 2
    conn = 4;
else
    conn = 6;
end;

% make sure that the splitting threshold is integer valued
splitvolvox = round(splitvolvox);

% extract the large regions what you want to split
cellbwin = cellbw;
cellbwlarge = bwareaopen(cellbw,splitvolvox,6); 
cellbwsmall = cellbwin - cellbwlarge;
cellbw = cellbwlarge;

% % extract one plane to work on
% if ~isequal(plane,'all')
%     cellbw = cellbw(:,:,plane);
% end;

cellbw = cellbw > 0;
% load ball1;se = getball(ball,1,1);
% cellbw = imdilate(cellbw,se);

[faser,L] = bwlabeln(cellbw);
msg = ['Number of objects due to splitting: ' int2str(L)];
disp(msg);
if L == 0
    warning('Number of objects due to splitting is too low, you may want to decrease the volume threshold');
    cellbw = cellbwin;
    return;
end;

% make double for interpolation
cellbw = double(cellbw);

% convert to high resolution (hr) image for bwdist
minh = min(h);
hhr = [minh,minh,minh];
fov = dim .* h;
dimhr = fov./hhr;

% resize to high resolution
cellbwhr = cellsegm.imresize3d(cellbw,dimhr,'linear');
cellbwhr = cellbwhr > 0;

% distance function
distimhr = bwdist(cellbwhr == 0);

% some smooothing of the distances
clear prmin;
prmin.gaussian.stdev = 3;
prmin.gaussian.diameter = 9;
distimhr = cellsegm.smoothim(distimhr,'gaussian','prm',prmin);

% resize back
distim = cellsegm.imresize3d(distimhr,dim,'linear');

% distim = imcomplement(distim);
% distim(cellbw == 0) = -Inf;
% wat = watershed(distim);
% showall(cellbw,distim,wat,wat > 0)

vis = 0;
if vis    
    disp('After distance function')
    showall(cellbw,distim)
end;

% maximum points
maximumbw = imextendedmax(distim,splitth);
% remove single voxels
load ball1;se = getball(ball,1,1);
maximumbw = imopen(maximumbw,se);

% background
if vis
    disp('After imextendedmax');
    showall(cellbw,distim,maximumbw);
end;

% find the Euclidean distances between the maxima
distim = distimage(maximumbw,conn);

% find the division lines
division = absgradient(distim,[1 1 1]);
division = division > 0;

% only within the large regions
division = division.*cellbw;
cellbw(division == 1) = 0;

% large and small regions
cellbw = cellbwsmall | cellbw;

% vis = 1
if vis
    disp('After distim')
    showall(cellbw,distim,maximumbw)
end;

% if there are artifacts from lines
load ball1;se = getball(ball,1,1);
cellbw = imopen(cellbw,se);

if vis
    disp('Final')
    showall(cellbwin,cellbw)
end;
