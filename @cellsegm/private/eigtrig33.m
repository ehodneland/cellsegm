function [V,D] = eigtrig33(A)
% EIGTRIG33 Computing eigenvalues and vectors for a 3x3 cell array A with a 
% matrix in each cell, using a trigonometric method to solve the
% characteristic equation.
%
% Returning the eigenvectors V and the eigenvalues D as cell arrays
%
% Taken from http://en.wikipedia.org/wiki/Eigenvalue_algorithm. Inspiration
% can also be found in Platero et al, "Analytic formulation for 3D
% diffusion tensor".
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

%
% Compute eigenvalues
%

if iscell(A)
    [V,D] = compcell(A);
else
    [V,D] = compreal(A);
end;


%-----------------------------------------


function [V,D] = compcell(A)

dim = size(A{1,1});
D = eyecell(3,dim);
I = eyecell(3,dim);
p1 = A{1,2}.^2 + A{1,3}.^2 + A{2,3}.^2;

% q is the mean trace
q = (A{1,1} + A{2,2} + A{3,3})/3;

% p2 = tr((A-qI))^2/6
p2 = (A{1,1} - q).^2 + (A{2,2} - q).^2 + (A{3,3} - q).^2 + 2 * p1;
p = sqrt(p2 / 6);

% I is the identity matrix   
%B = (1 / p) * (A - q * I);
B = multconstcell(I,-q);
B = sumcell(A,B);
B = multconstcell(B,1./(p+eps));

%r = det(B) / 2;
r = detcell33(B)/2;


% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
% if (r <= -1) 
%    phi = pi / 3;
% elseif (r >= 1)
%    phi = 0;
% else
%    phi = acos(r) / 3;
% end
ind1 = r <= -1;
ind2 = r >= 1;
phi = acos(r) / 3;
phi(ind1) = pi/3;
phi(ind2) = 0;

% the eigenvalues satisfy eig3 <= eig2 <= eig1
D{3,3} = q + 2 * p .* cos(phi);
D{1,1} = q + 2 * p .* cos(phi + (2*pi/3));
% since trace(A) = eig1 + eig2 + eig3
D{2,2} = 3 * q - D{1,1} - D{3,3};


%
% Compute eigenvectors
%
% v1 = (A - D(2,2)*I)*(A - D(3,3)*I);
% v2 = (A - D(1,1)*I)*(A - D(3,3)*I);
% v3 = (A - D(1,1)*I)*(A - D(2,2)*I);
a1 = multconstcell(I,-D{2,2});
a1 = sumcell(A,a1);
a2 = multconstcell(I,-D{3,3});
a2 = sumcell(A,a2);
v{1} = matrixprodcell(a1,a2);

a1 = multconstcell(I,-D{1,1});
a1 = sumcell(A,a1);
a2 = multconstcell(I,-D{3,3});
a2 = sumcell(A,a2);
v{2} = matrixprodcell(a1,a2);

a1 = multconstcell(I,-D{1,1});
a1 = sumcell(A,a1);
a2 = multconstcell(I,-D{2,2});
a2 = sumcell(A,a2);
v{3} = matrixprodcell(a1,a2);

% sum the columns
%v1 = sum(v1,2);
%v2 = sum(v2,2);
%v3 = sum(v3,2);
for j = 1 : 3
    for k = 1 : 3
        V{k,j} = v{j}{k,1} + v{j}{k,2} + v{j}{k,3};
    end
end;

% normalize the columns
%v1 = normc(v1);
%v2 = normc(v2);
%v3 = normc(v3);
for i = 1 : 3
    s = 0;
    for j = 1 : 3
        s = s + V{j,i}.^2;
    end;
    s = sqrt(s);
    for j = 1 : 3
        V{j,i} = V{j,i}./(s + eps);
    end;
end;


%----------------------------------------

function [V,D] = compreal(A)

D = eye(3,3);
I = eye(3,3);
p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2;
if (p1 == 0) 
    % A is diagonal.
    D = A;
else
    % q is the mean trace
    q = trace(A)/3;
    % p2 = tr((A-qI))^2/6
    p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1;
    p = sqrt(p2 / 6);
    % I is the identity matrix   
    B = (1 / p) * (A - q * I);
    r = det(B) / 2;
 
    
    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    % but computation error can leave it slightly outside this range.
    if (r <= -1) 
       phi = pi / 3;
    elseif (r >= 1)
       phi = 0;
    else
       phi = acos(r) / 3;
    end
 
    % the eigenvalues satisfy eig3 <= eig2 <= eig1
    D(3,3) = q + 2 * p * cos(phi);
    D(1,1) = q + 2 * p * cos(phi + (2*pi/3));
    % since trace(A) = eig1 + eig2 + eig3
    D(2,2) = 3 * q - D(1,1) - D(3,3);
end

%
% Compute eigenvectors
%

v1 = (A - D(2,2)*I)*(A - D(3,3)*I);
v2 = (A - D(1,1)*I)*(A - D(3,3)*I);
v3 = (A - D(1,1)*I)*(A - D(2,2)*I);

% eigenvector is the space spanned by all the columns!
v1 = sum(v1,2);
v2 = sum(v2,2);
v3 = sum(v3,2);

v1 = normc(v1);
v2 = normc(v2);
v3 = normc(v3);
V = [v1,v2,v3];

