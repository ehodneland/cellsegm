function [Q,R] = qrcell(A)
% QRCELL QR factorization of cell array with a matrix in each cell array.
%
% From Golub, Matrix computations, page 224 and 210
% Working both for matrices and cell arrays. 
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

[m n] = size(A);
if iscell(A)
    dimim = size(A{1,1});
else
    dimim = size(A);
end;

if ~iscell(A)

    Q = eye(m,n);
    H = cell(n,1);
    for j = 1 : n
    
        % compuate for subpart of A
        [v,beta(j)] = house(A(j:m,j));        
    
        % Householder transformation
        I = eye(m-j+1);
        H{j} = I-beta(j)*(v*v');

        % This becomes R iteratively on the upper diagonal
        A(j:m,j:n) = H{j}*A(j:m,j:n);
        if j < m
            A(j+1:m,j) = v(2:m-j+1);
        end;
        % Hfull = eye(m,n);
        % Hfull(j:n,j:n) = H;
        % Q = Q*Hfull;
    
    end;
    R = triu(A);
    % Back substitution to find Q
    % We dont need to mutliply all matrices at all steps since they are
    % partly identity matrices
    for j = n : -1 : 1
        Q(j:n,j:n) = H{j}*Q(j:n,j:n);        
    end;
    
else
    Q = eyecell(m,dimim);
    H = cell(n,1);
    for j = 1 : n
    
        % compute Householder vector for column
        [v,beta] = housecell(A(j:m,j));        

%         if ~isempty(find(isnan(v{1})))
%             showall(v{1,1},isnan(v{1,1}))
%         end
        
        % Householder transformation
        I = eyecell(m-j+1,dimim);
        H{j} = outerprodcell(v,v);
        H{j} = multconstcell(H{j},-beta);
        H{j} = sumcell(I,H{j});

        % This becomes R iteratively on the upper diagonal
        A(j:m,j:n) = matrixprodcell(H{j},A(j:m,j:n));
        if j < m
            A(j+1:m,j) = v(2:m-j+1);
        end;
        
        % Hfull = eyecell(m,dimim);
        % Hfull(j:n,j:n) = H;
        % Q = matrixprodcell(Q,Hfull);
    

    end;
    R = triucell(A);
    
    % Back substitution to find Q
    for j = n : -1 : 1
        Q(j:n,j:n) = matrixprodcell(H{j},Q(j:n,j:n));        
    end;
    
end;


% ------------------------------

function [v,beta] = housecell(x)

dim = size(x{1,1});

n = numel(x);
if n == 1
    v = x;
    beta = zeros(dim);
    return;
end;

% sig = x(2:n)'*x(2:n);
% for i = 1 : 3
%     x{i}
% end;
sig = innerprodcell(x(2:n),x(2:n));

one = cell(1,1);
one{1} = ones(dim);
v = [one ; x(2:end)];

% beta = zeros(dim);
% ind = ne(sig,0);
mu = sqrt(x{1}.^2 + sig);
v{1} = zeros(dim);
ind1 = x{1} <= 0;
v{1}(ind1 == 1) = x{1}(ind1 == 1) - mu(ind1 == 1);
v{1}(ind1 == 0) = -sig(ind1 == 0)./(x{1}(ind1 == 0) + mu(ind1 == 0)); 

% added extra
v{1}(isnan(v{1})) = 0;

beta = 2*(v{1}.^2)./(sig + v{1}.^2);
v = divideconstcell(v,v{1});
ind = eq(sig,0);
beta(ind) = 0;

% added extra
beta(isnan(beta)) = 0;

% if ~isempty(find(isnan(beta)))
%     showall(beta,isnan(beta))
% end

%----------------------------------

function [v,beta] = house(x)

n = numel(x);
if n == 1
    v = x;
    beta = 0;
    return;
end;
sig = x(2:n)'*x(2:n);
v = [1 ; x(2:end)];

if sig == 0
    beta = 0;    
else
    mu = sqrt(x(1)^2 + sig);

    if x(1) <= 0
        v(1) = x(1) - mu;
    else
        v(1) = -sig./(x(1) + mu);
    end;
    beta = 2*(v(1)^2)./(sig + v(1)^2);
    v = v./v(1);    
end;




