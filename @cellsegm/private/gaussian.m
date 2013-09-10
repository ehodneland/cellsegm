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



dim = round(dim);
midx = dim(1)/2;
midy = dim(2)/2;
midz = dim(3)/2;


[x y z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
x = x - midx-0.5;
y = y - midy-0.5;
z = z - midz-0.5;

if dim(3) == 1
    g = (1/(2*pi*sigma^2))*exp(-(x.^2 + y.^2)/(2*sigma.^2));
elseif dim(3) > 1
    g = (1/(2*pi*sigma^2)^(3/2))*exp(-(x.^2 + y.^2 + z.^2)/(2*sigma.^2));
else
    error('Wrong option for DIMZ')
end;
