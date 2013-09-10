%   BOUNDREG Finding the smallest rectangle around a binary object
%
%   [REG,BOX] = BOUNDREG(REG,PXY,PZ) Finding the smalles rectangle 
%   around the region in REG, 
%   PXY outside the smallest rectangle in XY plane, and PZ is the number of 
%   pixels outside the smalles rectangel. If PXY = PZ = 0, then it is close to the
%   region
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
function [reg,box] = boundreg(reg,pXY,pZ)

if isempty(find(reg,1))    
    box = [];
end;

[M N O] = size(reg);

[k(:,1) k(:,2) k(:,3)] = ind2sub(size(reg),find(reg));
minX = min(k(:,1));maxX = max(k(:,1));
minY = min(k(:,2));maxY = max(k(:,2));
minZ = min(k(:,3));maxZ = max(k(:,3));
box(1,1) = max(1,minX-pXY); box(1,2) = min(M,maxX+pXY); 
box(2,1) = max(1,minY-pXY); box(2,2) = min(N,maxY+pXY); 
box(3,1) = max(1,minZ-pZ);  box(3,2) = min(O,maxZ+pZ); 

reg = reg(box(1,1):box(1,2), box(2,1):box(2,2), box(3,1):box(3,2));
