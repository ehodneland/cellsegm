function [V,D] = eigcell(A)
% EIGCELL Computing eigenvalues and vectors for a cell array with a matrix
% in each cell
%
% Returning the eigenvectors V and the eigenvalues D as cell arrays
%
% From Panju, Iterative methods for computing eigenvalues and eigenvectors
% NB only square matrices
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


[m n] = size(A);

if iscell(A)
    dimim = size(A{1,1});
else
    dimim = size(A);
end;
if iscell(A)
    Q = eyecell(n,dimim);
else
    Q = eye(n,n);
end;

% can increase the number of iterations if necessary
for i = 1 : 15
    
    % compute QR factorization
    [q,r] = qrcell(A); 

    % flip around
    if iscell(A)
        A = matrixprodcell(r,q);
        Q = matrixprodcell(Q,q);
    else
        A = r*q;
        Q = Q*q;        
    end;
     
end;
D = A;
V = Q;

if iscell(A)
    for i = 1 : m
        for j = 1 : n
            D{i,j} = D{i,j} .* gt(abs(D{i,j}),1e-5);
            % can get nan for zero gradients in R in coh enh diffusion,
            % means zero on the diagonal here??         
            ind = isnan(V{i,j});
            V{i,j}(ind) = 0;
            ind = isnan(D{i,j});
            D{i,j}(ind) = 0;
        end;
    end;
else
    D = D .* gt(abs(D),1e-5);
end;
D = flipud(fliplr(D));
V = fliplr(V);



