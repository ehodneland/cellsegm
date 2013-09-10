%   BWRANGE Find the z-range of a binary component
%
%   [RANGE,MINZ,MAXZ] = BWRANGE(BW) Finding range of (one) connection in BW 
%   Returning the range of the connection in RANGE, and minimum and maximum z
%   planes in MINZ and MAXZ
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
function [range,minZ,maxZ] = bwrange(connhere)

[k(:,1) k(:,2) k(:,3)] = ind2sub(size(connhere),find(connhere));    
range = max(k(:,3)) - min(k(:,3)) + 1;
minZ = min(k(:,3));
maxZ = max(k(:,3));