function [cellbw] = splitcells(cellbw,splitth,plane)
% SPLITCELLS Splitting segmented binary image of cells 
%
%   CELLBW = SPLITCELLS(CELLBW,SPLITTH,PLANE) splitting the cells in CELLBW and the
%   splitting threshold in SPLITTH. SPLITTH is the second argument in
%   IMEXTENDEDMAX. PLANE is the plane where the distance 
%   function in 2D is computed, although the code runs for 3D data. 
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

dim = size(cellbw);
ndim = numel(dim);
if ndim == 2
    dim = [dim 1];
end;

if ndim == 2
    conn = 8;
else
    conn = 18;
end;

cellbwin = cellbw;
cellbw = cellbw(:,:,plane);

% input volumes
volin = bwsize(cellbw,conn);
meanvol = mean(volin);
minvol = min(volin);

% distance function must be done in 2D since it othwerwise becomes very
% strange
distim = bwdist(imcomplement(cellbw));

if vis    
    cellbwin = cellbw;
    disp('After distance function')
    showall(cellbw,distim)
end;

% maximum points
maximumbw = imextendedmax(distim,splitth);

if vis
    disp('After imextendedmax');
    showall(cellbw,maximumbw);
end;

% division lines
division = distimage(maximumbw,conn);
if vis
    disp('After distim2')
    showall(cellbw,division)
end;


% find lines
val = unique(division(maximumbw == 1));
perimall = zeros(dim(1:2));
for i = 1 : length(val)
    reghere = division == val(i);
   
   perim = bwperim(reghere);        
   perimall = perimall + perim;
   
end;
perimall = perimall > 0;
cellbw = cellbwin;
for i = 1 : dim(3)
    cellbwhere = cellbw(:,:,i);    
    cellbwhere(perimall>0) = 0;
    cellbw(:,:,i) = cellbwhere;
end;

% volume
[vol,faser] = bwsize(cellbw,conn);
ind = find(vol < 0.5*minvol);
for i = 1 : numel(ind)    
    cellbw(faser == ind(i)) = 0;
end;

if vis
    disp('Final')
    showall(cellbwin,cellbw)
end;
