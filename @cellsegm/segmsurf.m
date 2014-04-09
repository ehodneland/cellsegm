function [cellbw,wat,imsegm,minima,minimacell,info] = segmsurf(varargin)
% SEGMSURFWAT 3D segmentation of surface stained cells
%
%   [CELLBW,WAT,IMSEGM,MINIMA,MINIMACELL,INFO] = SEGMSURFWAT(IM,MINV,MAXV) 
%   Finding cells in a surface stained image IM. MINV and MAXV are the lower 
%   and upper expected volume of a cell in thousand cubic microns.
%   Returning the classified cells in CELLBW, the watershed image WAT, the
%   image used for segmentation IMSEGM (after smoothing),
%   the markers MINIMA, the cell markers MINIMACELL, 
%   information of parameters and classification descisions in INFO.
%  
%   [...] = SEGMSURFWAT(...,NAME,OPTION)
%   specifies an option in NAME for segmentation that can be:
%   PRMFILE         : Path to parameter file
%   MINIMA          : All markers, for instance from manual drawings
%   MINIMACELL      : Markers for cells, for instance from manual drawings
%   IMNUCLEUS       : Nucleus image for defining markers
%   PRM             : General parameters for SEGMSURFWAT
%   PRM can have no prefixes, or can have prefixes:
%   PRM.SMOOTHIM      -> Smoothing parameters (for SMOOTHIM)
%   PRM.CLASSIFYCELLS -> Cell classification parameters (for CLASSIFYCELLS)
%   PRM.GETMINIMA     -> Marker parameters for finding minima (for GETMINIMA)
%
%   Example 1
%   -------
%   
%   A segmentation with cells with volume of between 3000 and 50000 mcm3:
%
%   load data/surfstain_3D.mat
%   imsegm = imsegm(:,:,19);
%   [cellbw,wat,imsegmout,minima,minimacell,info] = cellsegm.segmsurf(imsegm,3,100);
%   show(imsegmout,1);title('Image used for segmentation')
%   show(wat,2);title('Watershed image');
%   show(cellbw,3);title('Segmented cells');
%
%   Example 2
%   -------
%
%   Specifying a particular method for smoothing. 
%
%   prmin.smoothim.method = 'ced';
%   load data/surfstain_3D.mat
%   imsegm = imsegm(:,:,19);
%   [cellbw,wat,imsegmout,minima,minimacell,info] = cellsegm.segmsurf(imsegm,3,100,'prm',prmin);
%   show(imsegmout,1);title('Image used for segmentation')
%   show(wat,2);title('Watershed image');
%   show(cellbw,3);title('Segmented cells');
%   
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
%

printmsg(['This is ' upper(mfilename) ' for segmentation of cells']);



% image for segmentation
imsegm = varargin{1};
% cell volumes in cubic microns
minvolfull = 1000*varargin{2};
maxvolfull = 1000*varargin{3};
imsegmini = imsegm;
prmfile = [];
minimacell = [];
minima = [];
imnucl = [];

% for visualilzation    
vis = 0;
visfinal = 0;

%
% Standard parameters
%

% voxel size
prm.h = [0.5 0.5 1.5];

% adjust for not fully circular shape of 3D cells
prm.just = 0.9;

% filter ridges
prm.filterridges = 0;

% GPU, requires Jacket for matlab!
prm.gpu = 0;

% correct illumination
prm.illum = 0;
prm.illumdiameter = 25;


%
% Minima parameters
%

% level
prm.getminima.level = 0.30;

% method
prm.getminima.method = 'automated';


%
% Smoothing parameters
%

% directional coherence enhancing diffusion
prm.smoothim.method = 'dirced';

% 2D or 3D smoothing
prm.smoothim.dim = 2;

% number of iterations of smoothing in CED
prm.smoothim.ced.maxniter = 100;


% dimension image
dim = size(imsegm);
if numel(dim) == 2
    dim = [dim 1];
end;

% 
% Classification parameters
%
prm.classifycells.method = 'threshold';


%
% Watershed parameters
%

% To test intensity of wat lines; merging
prm.merge = 0;

% methods for merging
prm.mergefragments.optlog = 2;
prm.mergefragments.optint = 2;

% The signicance of watershed lines, stronger if larger, done on ridgeim!!
% if the raw image, then use much smaller, around 1
% set to 1.05-1.10
% prm.watthint = 1.06;
prm.mergefragments.int = 1.3;

% The conexity measure for the merging, we should not allow a very much
% smaller convexity after merging
% thconv = 1 if no change in convexity
% thconv > 1 if increase in convexity
% thconv < 1 if decrease in convexity
prm.mergefragments.conv = 1.0;


if maxvolfull < minvolfull
    error([mfilename ': The upper cell volume is lower than the lower, check the input'])
end;

% read from parameter file?
for i = 4 : 2 : nargin
    namehere = varargin{4};
    varhere = varargin{5};
    switch(namehere)
        case('prmfile')
            prmfile = varhere;            
        otherwise
            % nothing
    end;
end;

%
% Load the parameter file if one is given
%
if ~isempty(prmfile)
    isfolder = ~isempty(strfind(prmfile,'/'));
    [a b c] = fileparts(prmfile);
    if isfolder        
        A = pwd;
        cd(a);        
    end;    
    msg = sprintf('Reading from parameter file %s',prmfile);
    disp(msg);
    cmd = [b c];
    prmin = eval(cmd);
    
    if isfolder
        cd(A);
    end;
    % merge the structs by overriding the existing values
    prm = mergestruct(prm,prmin);
end;


%
% The command line has the highest priority, merge the variables here
%
for i = 4 : 2 : nargin
    
    % this variable
    namehere = varargin{i};
    varhere = varargin{i+1};
    
    switch namehere
        case 'minima'
            minima = varhere;
        case('prm')
            prmin = varhere;
            prm = mergeinputpar(prm,prmin);            
        case('minimacell')
            minimacell = varhere;
        case('imnucleus')
            imnucl = varhere;           
        otherwise
            msg = ['Wrong option ' namehere];
            error(msg)                
    end;
end;

prm.dim = size(imsegm);

% if you give in a nucleus image you want to use it!
if ~isempty(imnucl)
    prm.getminima.method = 'nucleus';
    prm.classifycells.method = 'minimacell';
end;

% test for credibility!
if isequal(prm.getminima.method,'manual')
    if isempty(minima) && isempty(minimacell)
        error('You must specify either minima eller minimacell')
    end;
end;
if isequal(prm.getminima.method,'nucleus')
    if isempty(imnucl)
        error('You must specify a nucleus image');
    end;
end;
%
% Computed settings
%

% get the ADJUSTED 3D cell size, in case its not a full 3D volume. 
% NB: This is not always linearly scalable, try to avoid using reduced 3D.
[prm.minvol,prm.minvolvox,prm.maxvol,prm.maxvolvox] = ...
    cellsegm.cellsize(minvolfull,maxvolfull,prm.h,prm.just,dim(3));

% give the voxel sizes in thousand, since all programs are tuned for this
prm.classifycells.minvolfull = minvolfull/1000;
prm.classifycells.maxvolfull = maxvolfull/1000;
prm.getminima.minvolfull = minvolfull/1000;
prm.getminima.maxvolfull = maxvolfull/1000;
prm.watminvolfull = minvolfull/1000;
prm.watmaxvolfull = maxvolfull/1000;
prm.minvolfull = minvolfull;
prm.maxvolfull = maxvolfull;

disp('Using settings:');
printstructscreen(prm);

% pixel size
prm.getminima.h = prm.h;
prm.classifycells.h = prm.h;
prm.smoothim.h = prm.h;

% gpu for smoothing
prm.smoothim.gpu = prm.gpu;

%--------------------------------------------%
%-----------Cell segmentation start----------%
%--------------------------------------------%


%
% Correcting for uneven illumination
%
if prm.illum
    msg = ['Correcting uneven illumination'];
    disp(msg);
    se = strel('disk',prm.illumdiameter);
    a = imopen(imsegm,se);
    imsegm = imsegm - a;    
    if isequal(prm.getminima.method,'nucleus')    
        a = imopen(imsegm,se);
        imnucl = imnucl - a;        
    end;
    clear a;
end;

%
% Smoothing
%
try
    prmin = prm.smoothim;
catch
    
end;
imsegm = cellsegm.smoothim(imsegm,prm.smoothim.method,'prm',prmin);

% 
% Ridge enhancement
% 
if prm.filterridges == 1
    msg = ['Ridge enhancement'];
    disp(msg);
    imsegm = cellsegm.ridgeenhhessian(imsegm,prm.h);           
else    
    msg = 'No ridge enhancement';
    disp(msg);
end;

%
% Finding minima
%

createminima = 0;
if isequal(prm.getminima.method,'manual')
    if isempty(minimacell) || isempty(minima)
        createminima = 1;        
    end;    
elseif isequal(prm.getminima.method,'nucleus')
    createminima = 1;
elseif isequal(prm.getminima.method,'automated')
    createminima = 1;
else
    error('Wrong option in PRM.GETMINIMA.METHOD');
end;

% we create minima
if createminima
    % find the minima
    try
        prmin = prm.getminima;
    catch
    
    end;
    if isequal(prm.getminima.method,'nucleus')
        [minima,minimacell,prmout] = cellsegm.getminima(imsegm,minimacell,prmin,imnucl);
    else        
        [minima,minimacell,prmout] = cellsegm.getminima(imsegm,minimacell,prmin);
    end;    
else
    disp('Not creating any additional minima, only using manual minima');
end;

% showall(minima,minimacell,imsegm)

if vis == 1
    showall(imsegm,minima)
end;

%
% 3D watershed segmentation
%

msg = ['3D watershed segmentation'];
disp(msg);
Iimpose = imimposemin(imsegm,minima);
wat = watershed(Iimpose);
clear Iimpose
if prm.merge == 1
    % Test the significance of watershed lines, the new version
    % optlog = 1 : merge the strong lines (for CT)
    % optlog = 2 : merge the weak lines (for WGA)
    % optint = 1 : use the region based merging
    % optint = 2 : use the snake based merging
    msg = ['Merging watershed regions'];
    disp(msg);
    wat = cellsegm.mergefragments(wat,imsegm,prm.mergefragments.int,prm.mergefragments.conv,prm.mergefragments.optlog,prm.mergefragments.optint);
    
    % re segment after fooling around with the image
    msg = ['Resegment watershed'];
    disp(msg);
    wat = cellsegm.resegmentwat(wat,imsegm);
    
end;


if vis == 1
    show(minima,10)
    show(wat,11)
end;


%
% Classify cells
%
try
    prmin = prm.classifycells;
catch
    
end;
if ~isempty(minimacell)
    [cellbw,info.classifycells] = cellsegm.classifycells(wat,imsegmini,prmin,minimacell);
else
    [cellbw,info.classifycells] = cellsegm.classifycells(wat,imsegmini,prmin);
end;


if vis == 1 || visfinal == 1
    if dim(3) == 1
        show(imsegm,1);title('Image');
        show(wat,10);title('Wat');
        show(minima,11);title('Minima')
        show(minimacell,12);title('Minimacell')
        show(cellbw,13);title('Final segmentation');
        pause
    else       
        showall(imsegm,cellbw,wat,minima)
    end;
    
end;

% return variables
info.prm = prm;


        