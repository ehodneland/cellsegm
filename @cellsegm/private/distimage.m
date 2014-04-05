function [a] = distimage(minima,conn)
% DISTIMAGE Partition the image according to the Euclidian distances
%
%   A = DISTIMAGE(MINIMA,CONN) Partition the image according to the 
%   Euclidian distances from the markers in MINIMA to other pixels in image.
%   A pixel belongs to a certain phase if the distance to that marker is the
%   smallest of distances to all markers. CONN is the connectivity that is
%   used in labelling the minima
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


if isempty(find(minima,1))
    a = zeros(size(minima));
    return;
end;

[M N O] = size(minima);

% label the minima
[faser,L] = bwlabeln(minima,conn);

di = zeros([M N O L]);
for i = 1 : L
    regHere = eq(faser,i);
    di(:,:,:,i) = bwdist(regHere);
end;

[minVal,a] = min(di,[],4);
a = squeeze(a);
