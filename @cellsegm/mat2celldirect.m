function [A] = mat2celldirect(a)
% MAT2CELLDIRECT Convert a matrix into a cell matrix
% A = CELL2MATDIRECT(A) Converts the matrix A into a cel matrix which is
% returned
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



% if already a cell then dont do anything
if iscell(a)
    A = a;
    return;
end;

[M N O] = size(a);
A = cell(M,N,O);
for i = 1 : M
    for j = 1 : N
        for k = 1 : O
            A{i,j,k} = a(i,j,k);
        end;
    end;
end;
