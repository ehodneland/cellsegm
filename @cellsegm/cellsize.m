function [minvol,minvolvox,maxvol,maxvolvox] = cellsize(minvol,maxvol,h,just,O)
% CELLSIZE Computing the modified cell size in reduced data sets
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

if sum(h == 0) > 0
    error('There is zero voxel size')
end;
minr = ((3*minvol)/(4*pi))^(1/3);
% maxr = ((3*maxvol)/(4*pi))^(1/3);
% maxz = maxr/h(3);
minz = minr/h(3);
if O < minz
    % shrinink a sphere to a circle
    r = ((3*maxvol)/(4*pi))^(1/3);
    maxvol = O^just*pi*r^2;
    r = ((3*minvol)/(4*pi))^(1/3);
    minvol = O^just*pi*r^2;
end;
voxelvol = prod(h);
minvolvox = round(minvol/voxelvol);
maxvolvox = round(maxvol/voxelvol);
