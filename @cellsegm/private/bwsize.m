%   BWSIZE Find the size of binary segments
%
%   [DIM,FASER] = BWSIZE(BW,CONN) finds the size of the disconnected binary
%   segments inside BW and returns the number of pixels in DIM and the
%   piecewise constant image FASER connected to the label in DIM
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
function [dim,faser] = bwsize(BW,conn)

[faser,L] = bwlabeln(BW,conn);
val = faser(faser > 0);
valu = unique(val);
n = numel(valu);
dim = zeros(L,1);
for i = 1 : n
    dim(i) = nnz(val == valu(i));        
end;
