function [faser] = resegmentwat(faser,im)
% RESEGMENTWAT Resegment after merging
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

val = unique(faser(faser > 0));
numval = length(val);
minima = zeros(size(im));
load ball2;se = getball(ball,2,1);
for i = 1 : numval
    % this region
    reghere = eq(faser,val(i));

    % erode it
    reghere = imerode(reghere,se);
    
    % assign to minima
    minima(reghere) = 1;
end;
% minima = gt(faser,0);

% load ball2;se = getball(ball,2,1);
% minima = imerode(minima,se);
iimpose = imimposemin(im,minima);

% return variablees
faser = watershed(iimpose);