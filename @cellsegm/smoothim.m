function [im] = smoothim(varargin)
% SMOOTHIM  3D smoothing of images.
%   IM = SMOOTHIM(IM,METHOD) performs 3D smoothing of the iamge IM using
%   method METHOD. 
%   
%   Possible values for METHOD:
%   'ced'       : Coherence enhancing diffusion
%   'dirced'    : Directional coherence enhancing diffusion 
%   'eed'       : Edge enhancing diffusion
%   'gaussian'  : Gaussian smoothing
%   'none'      : No smoothing
%
%   IM = SMOOTHIM(IM,METHOD,'PRM',PRM) enables the parameter struct PRM
%   for command line based settings. Possible settings in PRM:
%
%   PRM.PLANEWISE : either planewise (1, default) or full 3D (0)
%   PRM.*DT       : controls the time step in method 1 and 3 (0.1, default)
%   PRM.*MAXNITER : maximum number of iterations in CED
%                   and EED (100, default) 
%   PRM.*DIAMETER : controls the diameter of the filter in method DIRCED and GAUSSIAN (13, default) 
%   PRM.*STDEV    : controls the standard deviation of the filter in method
%                   GAUSSSIAN (3, default)
%   PRM.*KAPPA    : controls the edge sensitivity in method CED and EED. This 
%                   is sensitiv to the scalar image values and a 
%                   normalization of the image is recommended 
%                   (0.001 for method CED, 10 for method EED, default)
%   The star * here means that it can be replaced by one of the mentioned 
%   methods, for instance PRM.CED.DT is the time step for method 'ced'.
%
%   Example
%   -------
%
%   load data/surfstain_3D.mat
%   imsegm = imsegm(:,:,15);
%   smimsegm = cellsegm.smoothim(imsegm,'ced');
%   show(imsegm,1);title('Raw image')
%   show(smimsegm,2);title('Smoothed image')
%
%   See also segmct, classifycells, getminima
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


im = varargin{1};
method = varargin{2};

% planewise
prm.planewise = 1;

% time step
prm.eed.dt = 0.1;
prm.ced.dt = 0.1;

% number of iterations
prm.ced.maxniter = 100;
prm.eed.maxniter = 100;

% no gpu needs Jacket!
prm.gpu = 0;

% kappa sensitivity to edges
prm.ced.kappa = 0.0001;
prm.eed.kappa = 10;

% diameter of filter
prm.dirced.diameter = 13;

% std of filter and diameter
prm.gaussian.diameter = 5;
prm.gaussian.stdev = 2;

% voxel size
prm.h = [0.5 0.5 1.5];

for i = 3 : 2 : nargin
    namehere = varargin{i};
    
    switch(namehere)
        case('prm')
            prmin = varargin{i+1};
            % merge input parameters
            [prm] = mergeinputpar(prm,prmin);
    end;
end;

if isequal(method,'none')
    msg = ['No smoothing in SMOOTHIM'];
    disp(msg);
    return;
end;

msg = ['This is ' upper(mfilename) ' using settings'];
disp(msg);
printstructscreen(prm);

if ~ismember(prm.planewise,[0 1])
    error('Wrong option in PLANEWISE');
end;

dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;

% if prm.gpu
%     im = gsingle(im);
% end;

if isequal(method,'ced')
    % coherence enhancing diffusion
    msg = 'Using coherence enhancing diffusion';
    disp(msg);
    if prm.planewise
        for i = 1 : dim(3)
            imhere = im(:,:,i);
            im(:,:,i) = cellsegm.cohenhdiff(imhere,prm.ced.dt,prm.ced.maxniter,prm.ced.kappa,prm.h);       
        end;
    else
        % true 3D
        im = cellsegm.cohenhdiff(im,prm.ced.dt,prm.ced.maxniter,prm.ced.kappa,prm.h,'opt','num');       
        error('Wrong option to SMPRM.SMOOTHDIM')        
    end;
    
    
elseif isequal(method,'dirced')
    % Directional coherence enhancing diffusion    
    msg = 'Using directional coherence enhancement filtering';
    disp(msg);
    if prm.planewise
        for i = 1 : dim(3)
            imhere = im(:,:,i);
            im(:,:,i) = cellsegm.dircohenh(imhere,prm.dirced.diameter,prm.h,prm.gpu);
        end;
    else
        % true 3D
         im = cellsegm.dircohenh(im,prm.dirced.diameter,prm.h);
    end;

elseif isequal(method,'eed')
    msg = 'Using edge enhancing diffusion';
    disp(msg);
    if prm.planewise
        for i = 1 : dim(3)
            im(:,:,i) = cellsegm.edgeenhdiff(im(:,:,i),prm.eed.dt,prm.eed.maxniter,prm.eed.kappa,prm.h); 
        end;
    else
        % edge enhancing diffusion
        im = cellsegm.edgeenhdiff(im,prm.eed.dt,prm.eed.maxniter,prm.eed.kappa,prm.h);
    end;
    
elseif isequal(method,'gaussian')
    msg = ['Using Gaussian smoothing'];
    disp(msg);
    
    msg = ['Using settings'];
    disp(msg);    
    printstructscreen(prm.gaussian);
    
    if prm.planewise
        filt = fspecial('gaussian',prm.gaussian.diameter,prm.gaussian.stdev);                
        for i = 1 : dim(3)
            im(:,:,i) = imfilter(im(:,:,i),filt,'replicate');
        end;
    else
        % edge enhancing diffusion
        diam = prm.gaussian.diameter*ones(1,3)./prm.h(1:3);
        im = smooth3(im,'gaussian',diam);
    end;
else        
    error('Wrong option to SMOOTHIM, terminating')
end;

% if prm.gpu
%     im = double(im);
% end;