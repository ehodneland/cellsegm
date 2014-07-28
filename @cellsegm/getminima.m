function [minima,minimacell,prm] = getminima(varargin)
% GETMINIMA Find markers for cell segmentation
%
%   [MINIMA,MINIMACELL,PRM] = GETMINIMA(IM,MINIMACELL,PRM) 
%   Finds the markers for image IM, cell markers
%   MINIMACELL, and settings in struct PRM. Let PRM = [] for default
%   settings and let MINIMACELL = [] for no cell markers. If MINIMACELL is
%   non-empty and prm.method = 'automatic', the method only finds the
%   background markers.
%
%   [...] = GETMINIMA(...,IMNUCL) also read the nucleus image to find 
%   markers. For this option to be used, set prm.method = 'nucleus'. 
%
%   General fields for PRM
%
%   PRM.H (default = [0.5 0.5 1.5])            
%   Voxelsize 1x3 vector
% 
%   PRM.MINVOL  (default = 5)
%   Lower cell volume in 1000 mcm3.
%
%   PRM.MAXVOL (default = 50);
%   Upper cell volume in 1000 mcm3.
% 
%   PRM.METHOD defines the choice of method for finding markers. There are
%   three options available:
%
%   PRM.METHOD = 'nucleus'
%   The markers are created from the nucleus channel
%
%   PRM.METHOD = 'manual'
%   Cell markers are given from manual delinear or other sources. This
%   option requires a non-empty binary image for the cell
%   markers, and only the background markers are found 
%   by the program.
%
%   PRM.METHOD = 'automated' (default)
%   The markers are created automatically
%
%   Specific settings for PRM.METHOD = 'automated':
%
%   PRM.AUTOMATED.ADFILTTH (default = 0.01)
%   Threshold in adaptive filtering. Default 0.01
%
%   Specific fields in PRM for PRM.METHOD = 'nucleus'
%   PRM.NUCLEUS*      : Any option to SEGMCT for segmentation of the
%                       nuclei. See SEGMCT for options. Default settings
%                       are PRM.NUCLEUS.METHOD = 'itth', PRM.NUCLEUS.SPLIT =
%                       0, PRM.NUCLEUS.ITTH = 1.
%
%   PRM.LEVEL (default = 0.3)
%   Normalized plane level 0 -> 1 to make markers. 
%                       If PRM.LEVEL = 'ALL' all planes 
%                       are used to create make markers. 
%                       If PRM.LEVEL = 'STRONG' the brightest plane is used
%                       as a marker image
%
%
%
%   See also segmct, classifycells, segmsurfwat
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
msg = ['This is ' upper(mfilename) ' finding minima'];
disp(msg);


imsegm = varargin{1};
minimacell = varargin{2};
prmin = varargin{3};
imnucl = [];

if nargin == 4
    imnucl = varargin{4};
end;
dim = size(imsegm);
if numel(dim) == 2
    dim = [dim 1];
end;

vis = 0;

%
% Minima parameters
%

% minimum and maximum volume
prm.minvolfull = 3;
prm.maxvolfull = 50;

% I change this 20120320 by adding 0.5 in front!
% prm.volfactor = 0.50;

% pixelsize
prm.h = [0.5 0.5 1.5];

% NB may not work if this is too high, flooding goes out
prm.rangemin = 0.4;

% how quickly to decrease the minima
prm.dimminball = 2;

% method
prm.method = 'automated';

% normalized level to make minima
prm.level = 0.30;

% adjust for 3D size
prm.just = 0.9;

% the threshold for the adaptive filter to find the minima(0.01-0.03)
prm.automated.adth = 0.01;
% the dimension of adaptive filter
prm.automated.diameter = 20;

% settings for segmct in nucleus segmentation
prm.nucleus.method = 'thrs';
prm.nucleus.thrs.th = 1.0;

% the planes for making markers, specified further down
prm.plane = [];

if numel(dim) == 2
    dim = [dim 1];
end;
prm.dim = dim;

% assign the input variables
if ~isempty(prmin)
    prm = mergeinputpar(prm,prmin);
end;
prm.minvolfull = prm.minvolfull*1000;
prm.maxvolfull = prm.maxvolfull*1000;

% voxel volume
prm.voxelvol = prod(prm.h);

% adjust cell volume
[prm.minvol,prm.minvolvox,prm.maxvol,prm.maxvolvox] = cellsegm.cellsize(prm.minvolfull,prm.maxvolfull,prm.h,prm.just,dim(3));
if isequal(prm.level,'all')
    prm.plane = 1 : dim(3);
elseif isequal(prm.level,'strong')
    % find the plane of strongest intensities. Especially useful for
    % tissue samples
    val = zeros(dim(3),1);
    for i = 1 : dim(3)
        imhere = imsegm(:,:,i);
        val(i) = mean(imhere(:));
    end;
    [maxval,ind] = max(val);
    prm.plane = ind;
else
    prm.plane = max(round(prm.level * dim(3)),1);
    prm.plane = min(dim(3),prm.plane);
    prm.plane = max(1,prm.plane);
end;

% adjust for not using all planes for finding minima
nplane = numel(prm.plane);
if ~isequal(prm.level,'all')    
    factor = nplane/dim(3);
    prm.minvol = round(prm.minvol*factor);
    prm.minvolvox = round(prm.minvolvox*factor);
    prm.maxvol = round(prm.maxvol*factor);    
    prm.maxvolvox = round(prm.maxvolvox*factor);    
end;

msg = ['Using settings'];
disp(msg);
printstructscreen(prm);

msg = ['Minima plane ' int2str(prm.plane)];
disp(msg);


%
% Finding minima
%

if isequal(prm.method,'manual')
    msg = ['Using MANUAL method for cell markers'];
    disp(msg);
    
    % get background
    minimabck = getbackground(imsegm,minimacell,prm);

    % add background and cell
    minima = gt(minimabck + minimacell,0);
    
elseif isequal(prm.method,'nucleus')  
    
    msg = ['Using NUCLEUS method for markers'];
    disp(msg);
    
    % add nucleus information    
    try
        prmin = prm.nucleus;
    catch
    
    end;
    prmin.h = prm.h;
    prmin.plane = prm.plane;
    prmin.minvolfull = prm.minvolfull;
    prmin.maxvolfull = prm.maxvolfull;
    minimacell = nucleus(imnucl,prmin);

    % get background minima    
    minimabck = getbackground(imsegm,minimacell,prm);

    % add background and cell
    minima = gt(minimabck + minimacell,0);

    if vis
        'After remove small-final'
        showall(imhere,minimabck)
    end;

%         showall(minima,minimabck)
elseif isequal(prm.method,'automated')
    
    msg = ['Using AUTOMATED method for markers'];
    disp(msg);

    try
        prmin = prm.automated;
    catch
    
    end;    
    prmin.h = prm.h;
    prmin.plane = prm.plane;
    prmin.minvolvox = prm.minvolvox;
    prmin.maxvolvox = prm.maxvolvox;
    if isempty(minimacell)        
        [minima,prm] = automated(imsegm,prmin);
    else
        msg = ['Only finding background since cells are known'];
        disp(msg);
        minimabck = getbackground(imsegm,minimacell,prmin);
        minima = minimacell;
        minima(minimabck == 1) = 1;
    end;
        
else
    error([mfilename ': Wrong option for prm.method']);
end;

%--------------------------------------------------

function [im] = repeatplane(im,start,erorad)

% dim = size(im);
% add = start-1;
% dimnew = dim;
% dimnew(3) = dimnew(3) + add;
% imnew = zeros(dimnew);
% imnew(:,:,start:end) = im;
dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;
name = ['ball' int2str(erorad)];
load(name);se = getball(ball,erorad,1);
imhere = im(:,:,start);
% down
for i = start-1:-1:max(start-5,1)    
    imhere = imerode(imhere,se);
    im(:,:,i) = imhere;
end;
imhere = im(:,:,start);
% up
for i = start+1:1:min(start+5,dim(3))    
    imhere = imerode(imhere,se);
    im(:,:,i) = imhere;
end;

% showall(im)


%-----------------------

function [minimabck] = getbackground(imsegm,minimacell,prm)


vis = 0;


dim = size(imsegm);
if numel(dim) == 2
    dim = [dim 1];
end;

% find background
imhere = scale(imsegm);    
th = 0.01;
totvol = numel(imsegm);
while 1        
    minimabck = lt(imhere,th);
    vol = nnz(minimabck);
    relvol = vol/totvol;
    if relvol > 0.10
        msg = ['Threshold for background: ' num2str(th)];
        disp(msg);
        break;
    end;
    th = th + 0.01;
end;

if vis
    'After threshold'
    showall(imhere,minimabck)
end;

% close with one
load ball1;se = getball(ball,1,1);
minimabck = imclose(minimabck,se);

% erode to remove from nucleus markers to avoid overlap
load ball2;se = getball(ball,2,1);
minimabck = imerode(minimabck,se);

% % remove maximum projection of the nucleus markers
% m = maxprojimage(minimacell);
% load ball4;se = getball(ball,4,1);
% m = imdilate(m,se);
% for i = 1 : dim(3)
%     ima = minimabck(:,:,i);
%     ima(m == 1) = 0;    
%     minimabck(:,:,i) = ima;
% end;

% open to disconnect to cells
load ball4;se = getball(ball,4,1);
minimabck = imopen(minimabck,se);

if vis
    'After close and open'
    showall(imhere,minimabck,minimacell)
end;

% remove small background minima
% no less than maximum cell volume in each background area
th = round(prm.maxvolvox);
for i = 1 : dim(3)
    minimabck(:,:,i) = bwareaopen(minimabck(:,:,i),th);
end;
% minimabck = bwareaopen(minimabck,round(prm.minvolvox/8));

if vis
    'After remove small'
    showall(imhere,minimabck,minimacell)
end;

% remove lower ones to avoid overlap with nucleus markers
[range,minz,maxz] = bwrange(minimacell);
if maxz < dim(3)
    minimabck(:,:,1:maxz) = 0;
end;
if vis
    'After removing lower ones'
    showall(imhere,minimabck,minimacell)
end;

% remove all markers with overlap to nucleus markers (minimacell)
% [faser,L] = bwlabeln(minimabck);
% overlap = faser(minimacell == 1);
% overlap = minimabck .* minimacell;
r = prm.h(1)/prm.h(3);
dxy = 10;
dz = round(dxy*r);
load ball10;se = getball(ball,dxy,dz);
dilminimacell = imdilate(minimacell,se);
minimabck(dilminimacell == 1) = 0;

% % took this away since it removed all background if there is only one voxel
% % overlap!
% overlap(overlap == 0) = [];
% overlap = unique(overlap);
% for i = 1 : numel(overlap)
%     minimabck(faser == overlap(i)) = 0;
% end;

if vis
    'After remove overlap'
    showall(imhere,minimabck,minimacell)
end;


% % minimabck = imfill(minimabck,'holes');
% minimabck = fillholes(minimabck,0,prm.minvolvox,4);
% load ball2;se = getball(ball,2,1);
% minimabck = imopen(minimabck,se);    
% minimabck = imclose(minimabck,se);  
% 
% if vis
%     'After fill,open.close'
%     showall(imhere,minimabck)
% end;


% remove small background minima, if some new small ones were created in
% previous steps
minimabck = bwareaopen(minimabck,round(prm.minvolvox/2));

if vis
    'Final'
    showall(imhere,minimabck)
end;

%---------------------------------------------------

function [minima] = nucleus(imnucl,prm)

vis = 0;

dim = size(imnucl);

% find minima in the nucldus staining using cell segmentation of CT cells!
% msg = ['In FINDMINIMANUCLEUS using plane ' int2str(prm.plane)];
% disp(msg);

imhere = imnucl(:,:,prm.plane);

% segmct

try
    prmin = prm.segmct;
catch
    
end;
prmin.h = prm.h;
[minimahere,wat] = cellsegm.segmct(imhere,0.1*prm.minvolfull/1000,prm.maxvolfull/1000,'prm',prmin);

if vis
    'after segmct'
    showall(imnucl(:,:,prm.plane),minimahere)
end;


if vis
    'after open'
    showall(imnucl(:,:,prm.plane),minimahere)
end;


minima = zeros(dim);
minima(:,:,prm.plane) = minimahere;

% repeat the plane downward to improve the segmentation
if numel(prm.plane) == 1
    minima = repeatplane(minima,prm.plane,2);
end;

if vis
    'after remove small'
    showall(imnucl(:,:,prm.plane),minimahere)
end;


%-------------------------------------------------------


function [minima,prm] = automated(im,prm)


vis = 0;

dim = size(im);


imhere = im(:,:,prm.plane);

% connectivity
conn = 4;    

msg = ['This is method AUTOMATED finding minima'];
disp(msg);

msg = ['Using settings'];
disp(msg);
printstructscreen(prm);


% thresholding
im = scale(im);
thim = logical(adaptfiltim(imhere,prm.diameter,prm.adth,prm.h));
if vis == 1
    disp('Before remove small')
    showall(imhere,thim)
end;

% remove at border
thim = constborder(thim,5,0);

% remove very small
thhere = max(3,round(0.2*prm.minvolvox));
thim = bwareaopen(thim,thhere,conn);
if vis == 1
    'after remove small and fix border'
    showall(imhere,thim)
end;

%
% find the small regions that can be filled, and keep those as thim
%
disp('  Closing and filling');
filled = zeros(size(thim));
step = 2;
high = 10;

% make constant around to also include the cells at border. 
% May also create problems in the classification due to half cells 
% NB this may also create some large minima in small images that make out the whole image!!!!

thim = constborder(thim,5,1);
thimOld = thim;


% iterative closing
for i = 0 : step : high
    
    thim = thimOld;
    if i > 0                
        name = ['ball' int2str(i)];
        load(name);se = getball(ball,i,1);
        thim = imclose(thim,se);
    end;

    % fill each region
    [faser,L] = bwlabeln(thim,conn);  
        
    for j = 1 : L
        regHere = eq(faser,j);    

        % reduce to save time
        [regHereRed,box] = boundreg(regHere,1,0);
        [m n o] = size(regHereRed);
    
        
        % fill                
        filledRed = logical(zeros(size(regHereRed)));                
        for  k = 1 : o
            filledRed(:,:,k) = logical(imfill(regHereRed(:,:,k),'holes') - regHereRed(:,:,k));                    
        end;

        filledRed = bwareaopen(filledRed,prm.minvolvox,conn);
        filledRed = filledRed - bwareaopen(filledRed,prm.maxvolvox,conn);              
        filledRed = logical(filledRed);

        % reduced old filled regions
        filledRedOld = filled(box(1,1):box(1,2),box(2,1):box(2,2),box(3,1):box(3,2));        
        filledRed = removeoverlap(filledRed,filledRedOld,conn);
        
        % put back
        filled(box(1,1):box(1,2),box(2,1):box(2,2),box(3,1):box(3,2)) = ...
            gt(filled(box(1,1):box(1,2),box(2,1):box(2,2),box(3,1):box(3,2)) + filledRed,0);       
        
    end;
%     showall(filled)
end;

% the minima becomes the filled regions
minima = filled;

if vis == 1
    'After closing and filling'
    showall(imhere,minima)
end


% MUST open large regions!
% estimate the diameter of a sphere
dimopen = 2*((3*prm.minvolvox)/(4*pi))^(1/3);
dimopen = round(dimopen/6);
dimopen = max(dimopen,1);
% dimopen = 10;
name = ['ball' int2str(dimopen)];
load(name);se  = getball(ball,dimopen,1);
[faser,L] = bwlabeln(minima);
minima = zeros(size(minima));
for i = 1 : L
    reghere = faser == i;
%     vol = sum(reghere(:));
    reghere = imopen(reghere,se);    
    minima(reghere == 1) = 1;
end;

% must erode again after opeining to take apart from each other minima we
% have opened
% Cannot erode if minvolvox is too small, then you remove everything

thhere = max(10,prm.minvolvox*2);
minima = erolarreg(minima,1,thhere);

if vis == 1
    'After erode and open'
    showall(imhere,minima)
end
 
% must remove small after open
minima = bwareaopenrange(minima,round(0.1*prm.minvolvox),conn);

% clear border
minima = logical(constborder(minima,5,0));

%
% Close each region separately
%

[faser,L] = bwlabeln(minima);
for j = 1 : L
    reghere = eq(faser,j);
    load ball3;se = getball(ball,3,1);
    reghere = imclose(reghere,se);
    minima(reghere) = 1;
end;


% fill small holes 
minima = fillholes(minima,0,prm.minvolvox,conn);

if vis == 1
    'After closing each region separately'
    showall(imhere,minima)
end


% get background minima
minimabck = getbackground(imhere,minima,prm);

% add to minima
minima(minimabck == 1) = 1; 

% make the minima over many planes if we use only one plane for minima
% NB this must come after the making of background, otherwise the
% background will flood the cells! 
% took this away 20120727
% minima = makeMinimaRange(minima,prm,prm.plane,dim(3));

% make 3D again
m = minima;
minima = zeros(dim);
minima(:,:,prm.plane) = m;

% repeat the plane downward to improve the segmentation
% NB only downward since this is often where the staining is missing
if numel(prm.plane) == 1
    minima = repeatplane(minima,prm.plane,2);
end;

% 3D again!!
if vis == 1  
    'After fill holes in background and clear border'
    showall(im,minima)
end

% remove small as final step
minima = bwareaopen(minima,prm.minvolvox);

if vis == 1
    'Final'
    showall(im,minima)
end

%--------------------------------------------

function [filledred] = removeoverlap(filledred,filledredold,conn)


[m n o] = size(filledred);

% loop over the reduced filled objects to find overlap, dont want
% overlap
for i = 1 : o
    % new region
    reg = filledred(:,:,i);
    
    % old region
    regold = filledredold(:,:,i);
    
    % count the objects in this plane of new region
    [faser,L] = bwlabeln(reg,conn);
    
    % check the overlap in this plane
    overlap = logical(reg.*regold);
    valoverlap = unique(faser(overlap));
    
    % if no overlap than continue to next object, we do nothing
    if isempty(valoverlap)
        continue;
    end;
        
    % remove the overlap regions
    for j = 1 : length(valoverlap)
        reg(faser == valoverlap(j)) = 0;
    end;
    
    % put back the info
    filledred(:,:,i) = reg;
    
end;

%----------------------------------------------------        
