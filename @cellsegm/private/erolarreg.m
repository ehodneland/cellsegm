function [BW] = erolarreg(BW,p,th)
% EROLARREG Erode large regions in image
%
%   EROLARREG(BW,P,TH) Eroding binary large regions above 
%   TH pixels using disc-like  structuring element with radius p.
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



[faser,L] = bwlabeln(BW);
large = zeros(size(BW));

for j = 1 : L
    regHere = eq(faser,j);
    ind = find(regHere);
    range = bwrange(regHere);
    numPix = size(ind,1)/range;
    if numPix > th
        large(ind) = 1;
        BW(ind) = 0;
    end;    
end;


% erosion
name = ['ball' int2str(p)];
load(name);se = getball(ball,p,1);
large = imerode(large,se);

% combine small and large
BW = logical(BW + large);
