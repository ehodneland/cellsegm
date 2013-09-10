% DERIV2CENTR5 Computing second derivatives in image
%
%   [UXX UYY UZZ UXY UYZ UXZ] = DERIV2CENTR5(U,HX,HY,HZ) Finding the second 
%   derivatives of U using a 5 neighbourhood and central differences.
%
%   Ex : [uxx,uyy,uzz,uxy,uyz,uxz] = deriv2centr5(u,1,1,1)
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
function [uxx,uyy,uzz,uxy,uyz,uxz] = deriv2centr5(u,hx,hy,hz)



[M N O] = size(u);

% central
uxx = (-u([3:end end end],:,:,:) + 16*u([2:end end],:,:,:) - 30*u + ...
      16*u([1 1:end-1],:,:,:) - u([1 1 1:end-2],:,:,:))/(12*hx^2);
uyy = (-u(:,[3:end end end],:,:) + 16*u(:,[2:end end],:,:) - 30*u + ...
      16*u(:,[1 1:end-1],:,:) - u(:,[1 1 1:end-2],:,:))/(12*hy^2);
  
if O < 3
    uzz = zeros(M,N,O);
else
    uzz = (-u(:,:,[3:end end end],:) + 16*u(:,:,[2:end end],:) - 30*u + ...
        16*u(:,:,[1 1:end-1],:) - u(:,:,[1 1 1:end-2],:))/(12*hz^2);
end;

%
% mixed
%
[ux,uy,uz] = derivforw3(u,hx,hy,hz);


% backward y
uxy = (3*ux - 4*ux(:,[1 1:end-1],:,:) + ux(:,[1 1 1:end-2],:,:))/(2*hy);
% backward x
uyx = (3*uy - 4*uy([1 1:end-1],:,:,:) + uy([1 1 1:end-2],:,:,:))/(2*hx);
uxy = (uxy + uyx)/2;

% backward z
if O < 3
    uyz = zeros(M,N,O);
else
    uyz = (3*uy - 4*uy(:,:,[1 1:end-1],:) + uy(:,:,[1 1 1:end-2],:))/(2*hz);
    uzy = (3*uz - 4*uz(:,[1 1:end-1],:,:) + uz(:,[1 1 1:end-2],:,:))/(2*hy);
    uyz = (uyz + uzy)/2;
end;


% backward z
if O < 3
    uxz = zeros(M,N,O);    
else
    uxz = (3*ux - 4*ux(:,:,[1 1:end-1],:) + ux(:,:,[1 1 1:end-2],:))/(2*hz);
    uzx = (3*uz - 4*uz([1 1:end-1],:,:,:) + uz([1 1 1:end-2],:,:,:))/(2*hz);
    uxz = (uxz + uzx)/2;
end;

