function [filtim] = dircohenh(varargin)
% DIRCOHENH Structural smoothing of image
%
%   DIRCOHENH(IM,D,H) Smoothes the image IM using a square filter with
%   diameter D and voxel size H which must be a 1x3 array
% 
%   DIRCOHENH(IM,D,H,GPU) specifies whether to use GPU (GPU=1) (by Jacket) or 
%   not (GPU=0)
%
%   Ex: filtim = structsmooth(im,13,[1 1 3]);
% 
%   Literature:
%   High-throughput Anaysis of Multispectral Imagse of Breast Cancer Tisue.
%   Umesh Adiga et al. IEEE Transactions on image processing, Vol 15, No 8,
%   August 2006
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

im = varargin{1};
prm.d = varargin{2};
prm.h = varargin{3};
prm.gpu = 0;
if nargin == 4
    prm.gpu = varargin{4};
end;

msg = ['This is ' upper(mfilename) ' using settings'];
disp(msg);
printstructscreen(prm);

% dimxy = 3;
dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;
h = prm.h;

% number of elements to remove before one takes the Olympic average. Set to
% 1-->3
numremele = 2;


%
% Make the filters
%
stepxy = 30;
% stepxy = 30;
% stepz = 45;
% stepz = 30;
deghere = (0:stepxy:(180-stepxy))';
numdeghere = length(deghere);
deg2D = [deghere repmat(90,numdeghere,1)];


% deghere  = (0 : stepxy : (360-stepxy))';
% numdeghere = length(deghere);
% deg3D = [deg2D ; ...
%          0 0 ; ...
%          deghere repmat(stepz,numdeghere,1)];

% this is quicker, not 100% correct but works as good
deg3D = [deg2D; ...
         0 0];
     
if dim(3) == 1
    deg = deg2D;
else
    deg = deg3D;
end;


% the half filter size
p = floor(prm.d/2);
   
% make filter
c = makefilter(deg,p);

% number of directions
numdir = size(c,2);
 
if prm.gpu
    % for Jacket speedup
    im = gsingle(im);
end;

% inisization
filtim = zeros(dim);
for i = 1 : numdir
    
    % coordinates of points in this direction after rotation of filter
    chere = c(i);

    % filter image for max value
    [maxim,sumim,numpoints] = maxfilt3(im,chere,numremele,dim,h,prm.gpu);

    % structural filtering   
    filtim = max(filtim,(sumim - sum(maxim,4))/ (numpoints-numremele));

end;

if prm.gpu
    % gather from GPU
    filtim = double(filtim);
end;

% ----------------------------------------------------------

function [maxim,sumim,numpoints] = maxfilt3(im,c,numremele,dim,h,gpu)

% % relative coordinate where we want the image values
% cx = c.x;
% cy = c.y;
% cz = c.z;

% grid
if dim(3) == 1
%     [x y] = ndgrid(h(1):h(1):h(1)*dim(1),h(2):h(2):h(2)*dim(2));
    [x, y] = ndgrid(1:dim(1),1:dim(2));
    z = ones(dim(1:2));
else    
%     [x y z] = ndgrid(h(1):h(1):h(1)*dim(1),h(2):h(2):h(2)*dim(2),h(3):h(3):h(3)*dim(3));
    [x, y, z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
end;

if gpu
    x = gsingle(x);
    y = gsingle(y);
    if dim(3) > 1
        z = gsingle(z);
    end;
end;


maxim = zeros([dim numremele]);
sumim = zeros(dim);
numpoints = length(c.x);
for i = 1 : numpoints

    % the coordinate where we want to find the values
    xhere = x + c.x(i);
    yhere = y + c.y(i);
    zhere = z + c.z(i);

    % NB!!!!!!!!!
    % Must not do this for z direction, then it goes wrong.
    % Do NOOOOOT replacestr any Inf or NaN by image values, then the filter
    % does not work anymore!!!!!!
%     xhere(xhere < h(1)) = h(1);xhere(xhere > dim(1)*h(1)) = dim(1)*h(1);
%     yhere(yhere < h(1)) = h(1);yhere(yhere > dim(2)*h(1)) = dim(2)*h(2);
    xhere(xhere < 1) = 1;xhere(xhere > dim(1)) = dim(1);
    yhere(yhere < 1) = 1;yhere(yhere > dim(2)) = dim(2);
    
    
    if dim(3) == 1
        if gpu
            imhere = interp2(im,yhere,xhere);        
        else
            imhere = interp2(y,x,im,yhere,xhere,'cubic',-Inf);
        end;
        
    else        
        if gpu
            imhere = interp3(im,yhere,xhere,zhere);
        else
            imhere = interp3(y,x,z,im,yhere,xhere,zhere,'cubic',-Inf);
        end;
    end;
    
    
    % sort previous max image and 
    sortarray = maxim;
    sortarray(:,:,:,numremele+1) = imhere;
    sorted = sort(sortarray,4,'descend');

    % max image
    maxim(:,:,:,1:numremele) = sorted(:,:,:,1:numremele);    
    sumim = sumim + imhere;

    
end;
%       showall(maxim(:,:,:,1),maxim(:,:,:,2),sumim)
    
%-------------------------------------------

function [c] = makefilter(deg,p)

% step = round(p/2);
% r = -p:step:p;
% this is a bit tricky.....??? is this optimal???
r = linspace(-p,p,7);

% degrees, sphaerical coordinates
% theta : angle to positive x-axis 
% phi   : angle to positive z-axis 
% r     : distance along vector to point
%
for i = 1 : size(deg,1)

    % the degrees
    theta = deg(i,1);
    phi = deg(i,2);        

    % the positions
    c(i).x = r*sind(phi)*cosd(theta);
    c(i).y = r*sind(phi)*sind(theta);
    c(i).z = r*cosd(phi);

end;

%--------------------------------------------

