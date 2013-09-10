% DERIVCENTR3 Finding derivatives
%
%   [UX UY UZ] = DERIVCENTR3(U,HX,HY,HZ) Finding the derivatives of U 
%   using a 3 neighbourhood and central differences.
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
function [ux,uy,uz] = derivcentr3(u,hx,hy,hz)

dim = size(u);
ndim = numel(dim);
ux = (u([2:end end],:,:,:) - u([1 1:end-1],:,:,:))/(2*hx);
uy = (u(:,[2:end end],:,:) - u(:,[1 1:end-1],:,:))/(2*hy);
if ndim >= 3
    uz = (u(:,:,[2:end end],:) - u(:,:,[1 1:end-1],:))/(2*hz);
else
    uz = zeros(dim);
end;

