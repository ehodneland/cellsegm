% DERIVFORW3 Finding derivatives using forward difference
%
%   [UX UY UZ] = DERIVFORW3(U,HX,HY,HZ) Finding the derivatives of U using a 
%   3 neighbourhood and forward differences.
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
function [ux,uy,uz] = derivforw3(u,hx,hy,hz)

[M N O] = size(u);
ux = (-u([3:end end end],:,:,:) + 4*u([2:end end],:,:,:) - 3*u)/(2*hx);
uy = (-u(:,[3:end end end],:,:) + 4*u(:,[2:end end],:,:) - 3*u)/(2*hy);
if O < 3
    uz = zeros(M,N,O);
else
    uz = (-u(:,:,[3:end end end],:) + 4*u(:,:,[2:end end],:) - 3*u)/(2*hz);
end;

