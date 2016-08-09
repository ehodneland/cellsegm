% GAUSSIAN Make a Gaussian filter
%
%   GAUSSIAN(DIM,SIGMA) Make a Gaussian filter of dimension DIM. 
%   The standard deviation is SIGMA. 
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
function [g] = gaussian(dim,sigma)


h = [1,1,1];
[x,minx,maxx] = cellcenteredgrid(dim,h);
% Center the grid around middle of filter
middim = dim/2;
for i = 1 : numel(x)
    x{i} = x{i} - middim(i)*h(i);
end;


% mu is zero, we center the distribution
if dim(1) > 1 && dim(2) == 1
    n = 2*pi*sigma^2;
    n = sqrt(n);
    g = (1/n)*exp(-(1/2)*(x{1}.^2)/(sigma^2));
elseif dim(3) == 1
    sigma = diag(sigma);
    detsigma = det(sigma);
    sigmainv = inv(sigma);
    n = detsigma*(2*pi)^2;
    g = (1/n)*exp(-(1/2)*(sigmainv(1,1)*x{1}.^2 + sigmainv(2,2)*x{2}.^2));
elseif dim(3) > 1
    sigma = diag(sigma);
    detsigma = det(sigma);
    sigmainv = inv(sigma);
    n = detsigma*(2*pi)^2;
    g = (1/n)*exp(-(1/2)*(sigmainv(1,1)*x{1}.^2 + sigmainv(2,2)*x{2}.^2 + sigmainv(3,3)*x{3}.^2));
else
    error('Wrong option for DIMZ')
end;
% normalize to 1 due to cut off of the filter, discrete effect!
g = g/sum(g(:));
