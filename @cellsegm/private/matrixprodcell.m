function [ab] = matrixprodcell(a,b)
% MATRIXPRODCELL computes matrixproduct A*B of cell arrays
%
% AB = MATRIXPRODCELL(A,B) computes matrix product of two cell arrays A and
% B
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
dima = size(a);
dimb = size(b);
dimim = size(a{1,1});

ab = cell(dima(1),dimb(2));
for i = 1 : dima(1)
    for j = 1 : dimb(2)
        ab{i,j} = zeros(dimim);
        for k = 1 : dimb(1)            
            ab{i,j} = ab{i,j} + a{i,k} .* b{k,j};            
        end;
    end;
end;

    




