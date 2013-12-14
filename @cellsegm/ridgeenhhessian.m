function [v] = ridgeenhhessian(varargin)
% Ridge enhancement using Hessian matrix
% V = RIDGEENHSHESSIAN(U,H) produces a ridge enhancement of the image I using
% the Hessian matrix. H = [HX,HY,HZ] is the pixel relations.
%
% V = RIDGEENHSHESSIAN(..,name,'optoin') also adds the filter
% radius 'rad' or/and the standard deviation 'stdev' of the Gaussian smoothing
%
% Ex: filtim = ridgeenhhessian(filtim,[1 1 3]);
%
% 
% Ex: Gaussian filter radius 3
% filtim = ridgeenhhessian(filtim,[1 1 3],'rad',3);
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


u = varargin{1};
prm.h = varargin{2};
% default values
% was at 3
prm.rad = 3;
prm.radz = 1;
prm.stdev = 2;
if nargin > 2
    for i = 3 : 2 : nargin
        name = varargin{i};
        var = varargin{i+1};
        switch name
            case 'rad'
                prm.rad = var;
            case 'stdev'
                prm.stdev = var;
            otherwise
                error('Wrong option')        
        end;            
    end;    
end;
prm.h = rowvect(prm.h);

msg = ['This is RIDGENHHESSIAN using settings'];
disp(msg);
printstructscreen(prm);

[M N O] = size(u);

u = scale(u);

dim = [prm.rad prm.rad prm.rad]./prm.h;
g1 = gaussian(dim,prm.stdev);
% dim = [7 7 1];
% g2 = gaussian(dim,3);
u = imfilter(u,g1,'replicate');

% showall(u)

% derivatives second order
[uxx,uyy,uzz,uxy,uyz,uxz] = deriv2centr5(u,prm.h(1),prm.h(2),prm.h(3));

% intiell smoothing to avoid rapid oscillations
% stdDev = 3;
% dim = 9;
uxx = imfilter(uxx,g1,'replicate');
uyy = imfilter(uyy,g1,'replicate');
uzz = imfilter(uzz,g1,'replicate');
uxy = imfilter(uxy,g1,'replicate');
uxz = imfilter(uxz,g1,'replicate');
uyz = imfilter(uyz,g1,'replicate');

   
E = zeros(M,N,O,3);
for i = 1 : M
    for j = 1 : N               
        for k = 1 : O
            t = ([uxx(i,j,k) uxy(i,j,k) uxz(i,j,k) ;
                  uxy(i,j,k) uyy(i,j,k) uyz(i,j,k) ;
                  uxz(i,j,k) uyz(i,j,k) uzz(i,j,k)]);
            
            E(i,j,k,:) = eig(t);
        end;
    end;
end;

% ridge in 3D (also for 2D ok)
v = -E(:,:,:,1) - E(:,:,:,2).^2 - E(:,:,:,3).^2;

