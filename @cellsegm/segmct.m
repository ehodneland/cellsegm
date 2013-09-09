function [cellbw,wat,imsegm,prmout] = segmct(varargin)
% SEGMCT Segmentation of cytosplasmically stained cells and nuclei.
%
%   [CELLBW,WAT,IMSEGM,PRMOUT] = SEGMCT(IM,MINV,MAXV) Segmentation of 
%   cytoplasmically stained image IM. Can also be used for Hoechst stained nuclei.
%   MINV and MAXV are the expected minimum and maximum object volumes.
%   The cells are given as disconnected components in a binary image CELLBW, 
%   and labelled in WAT. 
%
%   [...] = SEGMCT(...,'PRM',PRM) gives an option to
%   specify a set of parameters.
%
%   PRM.METHOD (default = 'ITTH')
%   There are several methods available, specified by PRM.METHOD
%
%       1. ADAPTIVE THRESHOLDING, PRM.METHOD = 'ADTH'. Parameters to set are 
%       PRM.ADTHADTH (default = 0.20) and PRM.ADTHFILTRAD (default = 40).
%       NB: This method is robust but slow for 3D data.
%
%       2. ITERATIVE THRESHOLDING (default), PRM.METHOD = 'ITTH'. Iterative 
%       thresholding of the image until the largest region has the
%       volume of the smallest expected cell volume. This method is suitable
%       for stained nuclei. Parameters to set are PRM.ITTH.ITTH (default = 0.8) 
%       which is the lowest threshold to allow. Note that the threshold is not 
%       absolute but a multiple of the threshold arising using GRAYTHRESH with 
%       no arguments.
%   
%
%   PRM.SPLIT (default = 0)
%   Controls splitting of cells that are incorrectly connected. Using
%   distance function and local minima. Can be useful if the cells (nuclei)
%   are incorrectly connected. Can also be run separately. 
%
%   PRM.SPLITTH (default = 1)
%   The amount of splitting is controlled by PRM.SPLITTH
%   where lower values gives stronger splitting. PRM.SPLITTH is the second 
%   input argument to IMEXTENDEDMAX in MATLAB.
%
%   Example
%   -------
%
%   examples/nucleistain_2D;
%
%
%   Example
%   -------
%
%   examples/nucleistain_3D;
%
%
%   See also segmsurf, splitcells, getminima
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



msg = ['This is ' upper(mfilename) ' for segmentation of cytoplasmically stained cells'];
disp(msg);

vis = 0;

im = varargin{1};
minvolfull = 1000*varargin{2};
maxvolfull = 1000*varargin{3};

% default method
prm.method = 'thrs';

% gaussian smoothin
prm.smoothim.method = 'gaussian';

% illumination correction
prm.illum = 0;

%
% Iterative threshold
%

% threshold
prm.thrs.th = 1.2;

%
% Adative thresholding
%

% threshold
prm.adth.adth = 0.20;

% filter radius in microns
prm.adth.filtrad = 40;


% %
% % Gradient method
% % 
% prm.gradientadth = 0.2;
% prm.gradientthrsth = 2;


%
% Splitting
%
prm.splitth = 1;

% splitting of cells
prm.split = 0;

prm.h = [0.5 0.5 1.5];
prm.just = 0.9;
dim = size(im);
ndim = numel(dim);
if numel(dim) == 2
    dim = [dim 1];
end;

% is there a paramter file given?
for i = 4 : 2 : nargin
    namehere = varargin{i};
    
    switch(namehere)
        case('prm')
            prmin = varargin{i+1};
            % merge input parameters
            prm = mergeinputpar(prm,prmin);            
    end;
end;
disp('Segmentation by SEGMCT.M');
[M N O] = size(im);

if O == 1 
    conn = 8;
else
    conn = 26;
end;


% voxelvol
prm.voxelvol = prod(prm.h);

% expected number of pixels in a cell per plane
% [prm.minvol prm.maxvol] = cellsize(diamlow,diamhigh,1,1,0.9,O);
[prm.minvol,prm.minvolvox,prm.maxvol,prm.maxvolvox] = ...
cellsegm.cellsize(minvolfull,maxvolfull,prm.h,prm.just,dim(3));
prm.minvolfull = minvolfull;
prm.maxvolfull = maxvolfull;

msg = ['This is SEGMCT using settings'];
disp(msg);
printstructscreen(prm);

% scale image
im = scale(im);

%
% Correcting for uneven illumination
%
if prm.illum
    msg = ['Correcting uneven illumination'];
    disp(msg);
    se = strel('disk',25);
    im = imtophat(im,se);
end;

%
% Smoothing
%
try
    prmin = prm.smoothim;
catch
    
end;
prmin.h = prm.h;
im = cellsegm.smoothim(im,prm.smoothim.method,'prm',prmin);

% Almost 
% High througput analysis of multispectral images of breast cancer
% tissue, Umesh Adiga, IEEE, Vol8, No 8, 2006
msg = ['Using segmentation method ' upper(prm.method) ' in SEGMCT'];
disp(msg);
if isequal(prm.method,'thrs')
    % Segmentation by iterative thresholding until largest region is larger
    % than expected cell volume
    [cellbw,wat] = segmthrs(im,prm,conn);
elseif isequal(prm.method,'gradient')
    % Almost 
    % Combining intensity, edge and shape information for 2d and 3D
    % segmentaton fo cell nuclei in tissue sections, Wahlby, Journal of
    % microscopy, Vol 215, pp 67-76, 2004
    % This one is quite good for Tanjas NRK cells!
    [cellbw,wat] = segmgrad(im,prm);
elseif isequal(prm.method,'adth')
    % adpative thresholding
    [cellbw,wat] = segmadth(im,prm,conn);
else
    error('No valid method given');
end;

% split cells?
if prm.split
    msg = ['Splitting cells using distance function'];
    disp(msg);
    plane = ceil(dim(3)/3);
    [cellbw] = cellsegm.splitcells(cellbw,prm.splitth,plane);
   
    % relabel
    [wat,L] = bwlabeln(cellbw);
end;

% remove the large cells
% Lin = Lout;
cellbw = cellbw - bwareaopen(cellbw,prm.maxvolvox,6);

% to return
imsegm = im;
prmout = prm;

[wat,L] = bwlabeln(cellbw);
msg = ['Number of objects left: ' int2str(L)];
disp(msg);

% [faser,Lout] = bwlabeln(cellbw);
% msg = ['Removed ' int2str(Lin-Lout) ' regions due to large size'];
% disp(msg);


% show(im,1);show(cellbw,2)

%--------------------------------------------------------------------

function [cellbw,wat] = segmadth(im,prm,conn)

% adaptive thresholding
d = round(prm.adth.filtrad/mean(prm.h(1:2)));
msg = ['Adaptive filtering with filter radius ' num2str(d) ' and threshold ' num2str(prm.adth.adth)];
disp(msg);
cellbw = adaptfiltim(im,d,prm.adth.adth);

% remove the small cells
disp('Remove small parts that are not cells')
[faser,Lin] = bwlabeln(cellbw,conn);
cellbw = bwareaopen(cellbw,round(prm.minvolvox),6);
[faser,Lout] = bwlabeln(cellbw,conn);
msg = ['Removed ' int2str(Lin-Lout) ' regions due to small size'];
disp(msg);

% % remove the large cells
% Lin = Lout;
% cellbw = cellbw - bwareaopen(cellbw,prm.maxvolvox,6);
% [faser,Lout] = bwlabeln(cellbw,conn);
% 
% msg = ['Removed ' int2str(Lin-Lout) ' regions due to large size'];
% disp(msg);

% % split the cells
% cellbw = splitcells(cellbw,conn);

% fill holes
holes = imfill(cellbw,'holes') - cellbw;
holes = holes - bwareaopen(holes,prm.minvolvox,conn);
cellbw(holes == 1) = 1;

% erode and open to disconnects
load ball2;se = getball(ball,2,1);
cellbw = imopen(cellbw,se);
load ball1;se = getball(ball,1,1);
cellbw = imerode(cellbw,se);

% return
[wat,L] = bwlabeln(cellbw);


%-----------------------------------------------------------------


function [cellbw,wat] = segmgrad(im,prm)

dim = size(im);
ndim = numel(dim);
if numel(dim) == 2
    dim = [dim 1];
end;


if ndim == 2
    [dx dy] = gradient(im);
    absgrad = sqrt(dx.^2 + dy.^2);
elseif ndim == 3
    [dx dy dz] = gradient(im);
    absgrad = sqrt(dx.^2 + dy.^2 + dz.^2);
else
    error('Wrong dimension');
end;

% get minimacell using one of the internal routines
conn = 6;

% not too low threshold, then you 
prmin.thrsth = prm.gradientthrsth;
prmin.minvolvox = prm.minvolvox;
prmin.maxvolvox = prm.maxvolvox;
msg = ['Finding markers for a watershed segmentation'];
disp(msg);
[minimacell,wat] = segmthrs(im,prmin,conn);

% use segmsurfwat to segment further
prmin.filterridges = 0;

% smoothing has already been done
prmin.smoothimmethod = 'none';
prmin.getminimamethod = 'manual';
prmin.classifycellsmethod = 'minimacell';
prmin.illum = 0;
[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurfwat(absgrad,prm.minvolfull/1000,prm.maxvolfull/1000,...
    'minimacell',minimacell,'prm',prmin);

%-----------------------------------------------------------------
function [cellbw,wat] = segmdist(im,prm,conn,opDiam)

[M N O] = size(im);

% filt = fspecial('gaussian',13,7);
% im = imfilter(im,filt,'replicate');

test = 0;

% 
% Make the segmentaiton essemtially, foreground seeds
%

% whos
% pause
% must have this for the uint8 thing!!
im = 255*scale(im);
im = uint8(im);
% NB, lavere gir høyere!!!
cellbw = zeros(size(im));
level = prm.otsuth*graythresh(im);
level = max(0,level);
level = min(level,1);
for i = 1 : O
    cellbw(:,:,i) = im2bw(im(:,:,i),level);
end;
im = double(im);

if test
    'graythresh'
    showall(im,cellbw)
end;

% remove small parts
name = ['ball' int2str(opDiam)];
load(name);se = getball(ball,opDiam,1);
cellbw = imopen(cellbw,se);

% remove small parts
disp('Remove small and big parts')
cellbw = bwareaopen(cellbw,prm.minvol,conn);    
% cellbw = cellbw - bwareaopen(cellbw,thpixhigh,conn);    
% showall(im,cellbw)
if test
    'small'
    showall(im,cellbw)
end;

dist = bwdist(imcomplement(cellbw));   
filt = fspecial('gaussian',13,9);
dist = imfilter(dist,filt,'replicate');
dist = imcomplement(dist);
if test
    'dist'
    showall(im,cellbw,dist)
end;

% make minima here!!!
minima = imregionalmin(dist,conn);
load  ball4;se = getball(ball,4,O);
minima = imdilate(minima,se);
minima = imclose(minima,se);
minima(cellbw == 0) = 1;

% make the "ridge"
dist(cellbw == 0) = min(dist(:));
% for the VECT
dist = scale(dist);

% % make a "fake" 3D minima!!!
% mid = round(O/2);
% if round(O/2) == O/2
%     vect = [mid:-1:1 1:1:mid];
% else
%     vect = [mid:-1:2 1 2:1:mid];
% end;
% vect = vect.^2;
% vect = 0.02*scale(vect);
% % dist = scale(dist);
% % vect = vect * 0.05;
% disp(sprintf('Weighting the planes to create a 3D minima using vect = %s',num2str(vect)));
% % add the weights
% for i = 1 : O
%     dist(:,:,i) = dist(:,:,i) + vect(i);
% end;

% %
% % To combine minima we make markers!!
% %
% minima = imregionalmin(dist,conn);
% % minima = eroLarMin(minima
% load  ball5;se = getball(ball,5,O);
% minima = imdilate(minima,se);
% minima = imclose(minima,se);

% showall(im,dist,minima)

% this watershed captures the boundaries found in graythresh, but also
% divides the connected cells!!!
Iimpose = imimposemin(dist,minima);
wat = watershed(Iimpose);
% save test dist wat im
% showall(im,dist,minima,wat)
if test
    'wat'
    showall(im,minima,wat)
end;


% cellbwold = cellbw;
% watold = wat;
%
% Merge
%
% if we make it smaller we merge weaker boundaries, i.e. if too little
% merging make it smaller, if too much merging (we loose tru regions to
% background) make it bigger
thint = 0.86;
% smaller than zero means that the region must become more convex after
% merging
thconv = 0.87;
optlog = 1;
optint = 1;
wat = mergefragments(wat,im,thint,thconv,optlog,optint);

% re segment after fooling around with the image
wat = resegmentwat(wat,dist);

% make cells from watershed
disp('Classify cells')
cellbw = classify(wat,im,prm.minvol,prm.maxvol);

% showall(im,cellbw,dist,minima,wat,cellbwold,watold)
% pause

%------------------------------------------------------------

function cellbw = classify(wat,im,thpixlow,thpixhigh)

meanall = mean(im(:));
L = max(wat(:));

cellbw = zeros(size(im));
for i = 1 : L
    reghere = eq(wat,i);
   
    % mean intensities in this region
    meanint = mean(im(reghere));
    
    % the reatio to the rest of the image
    ratio = meanint/meanall;
    
    % number of pixels
    numpix = length(find(reghere));
    
%     show(im,1)
%     show(reghere,2)
%     ratio

%     pause
    % was at 1.2
    th = 1.2;
    if ratio > th && numpix > thpixlow && numpix < thpixhigh
        cellbw(reghere) = 1;
    else        
        disp('Not cell,')
        disp(sprintf('Ratio intensity: %6.2f [ ''>'' %6.2f ]',ratio,th))
        disp(sprintf('Number pixels: %6.2f [%6.2f  ''>'' %6.2f ]',numpix,thpixhigh,thpixlow))
        
    end;
end;



%------------------------------------------------------------

% segmentation with iterative thresholding until regions are larger than
% expected cell volume
function [cellbw,wat] = segmthrs(im,prm,conn)

vis = 0;

% [M N O] = size(im);
dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;

im = 255*scale(im);
im = uint8(im);
% NB, lavere gir høyere!!!
cellbw = zeros(size(im));
% cellbwold = cellbw;

th = 2;
% interatively until we reach the given threshold
while 1
    level = th*graythresh(im);
    % IA dd this since hte crasched when level > 1, which can happen
    level = min(level,0.5);
    
    for i = 1 : dim(3)
        cellbw(:,:,i) = im2bw(im(:,:,i),level);
    end;

    [vol,faser] = bwsize(cellbw,conn);    
    msg = ['Thresholding at ' num2str(th)];
    disp(msg);
    
    if th <= prm.thrs.th
        msg = ['Threshold is lower than specified, terminating loop'];
        disp(msg);
        break;
    end;
    
    if max(vol) > 2*prm.maxvolvox
        msg = ['Threshold creates too large regions, terminating loop'];
        disp(msg);
        % use the last good segmentation
%         cellbw = cellbwold;
        break;
    end;

    th = th - 0.1;
   
%     cellbwold = cellbw;
end;
im = double(im);
cellbw = bwareaopen(cellbw,round(prm.minvolvox/10));

if vis == 1
    disp('After remove small parts')
    showall(im,cellbw)
end;

% connect each
load ball3;se = getball(ball,3,1);
[faser,L] = bwlabeln(cellbw);
for i = 1 : L
    reghere = faser == i;
    reghere = imclose(reghere,se);
    cellbw(reghere == 1) = 1;
end;

%
% Fill holes
%
disp('Fill holes')
for i = 1 : dim(3)
    cellbw(:,:,i) = imfill(cellbw(:,:,i),'holes');
end;

load ball1;se = getball(ball,1,1);
cellbw = imerode(cellbw,se);

% must do a opening to disconnect
% was at 1
load ball7;se = getball(ball,7,1);
cellbw = imopen(cellbw,se);



% load ball3;se = getball(ball,3,1);
% cellbw = imopen(cellbw,se);


% % erode to disconnect after open
% load ball1;se = getball(ball,1,1);
% cellbw = imerode(cellbw,se);


if vis == 1
    disp('After filling of holes and opening')
    showall(im,cellbw);
end;

% % cut the cells
% disp('Cutting cells')

% if vis == 1
%     disp('After cutting')
%     showall(im,cellbw)
% end;

% remove the small cells
disp('Remove small parts that are not cells')
[faser,Lin] = bwlabeln(cellbw);
cellbw = bwareaopen(cellbw,prm.minvolvox,6);
[faser,Lout] = bwlabeln(cellbw);

msg = ['Removed ' int2str(Lin-Lout) ' regions due to small size'];
disp(msg);

% Cannot do this before split!
% % remove the large cells
% Lin = Lout;
% cellbw = cellbw - bwareaopen(cellbw,prm.maxvolvox,6);
% [faser,Lout] = bwlabeln(cellbw);
% 
% msg = ['Removed ' int2str(Lin-Lout) ' regions due to large size'];
% disp(msg);

% % fine tune segmentation
% intth = 0.7;
% vis = 1;
% if vis == 1
%     cellbwold = cellbw;
% end;
% cellbw = dilatestrong( cellbw,im,intth,conn);
% 
% 
% if vis == 1
%     disp('After dilate strong')
%     showall(im,cellbwold,cellbw)
% end;

% mid = ceil(dim(3)/3);
% cellbw = splitcells(cellbw,1,mid);


% output watershed image
[wat,Lwat] = bwlabeln(cellbw);


if vis == 1
    disp('Final')
    showall(im,cellbw)
end;

% ----------------------------------------------

% this method is dilating one pixel at a time, dont think it makes a big
% differnece to dilatestrong
function [cellbw] = dilatestrong2(bw,im,intth,conn)

imini = im;
load ball2;se = getball(ball,2,1);
niter = 0;
while 1
    niter = niter + 1;
    bwold = bw;
    
    % label the binary regions
    [faser,L] = bwlabeln(bw,conn);
%     showall(im,bwold,bw,faser)
        
    % dilated border
%     border = zeros(size(im));
    faserth = zeros(L,1);
    for i = 1 : L
        reghere = faser == i;
        
        % the dilation
        borderhere = imdilate(reghere,se) - reghere;
        
%         % update BORDER to use it later 
%         border(borderhere == 1) = 1;
        
        % update the FASERBORDER to know which FASER the dilated pixel
        % belongs to
        faser(borderhere == 1) = i;
        
        % the threshold for this FASER
        faserth(i) = mean(im(reghere));
    end;
    
%     showall(im,bw,faser)
    % the border to sort the values
    border = gt(faser,0) - bw;

    % sort the pixels
    ind = find(border == 1);
    val = im(ind);
    [sorted,indsorted] = sort(val,1,'descend');
    
    % loop over the sorted pixels to start wiht the smallest one!!
    for i = 1 : length(val)
        indarray = indsorted(i);
        indimage = ind(indarray);
        
        % this image value
        valhere = val(indarray);
        
        % this faser value to find the connection to FASER to check for
        % threshold value
        faservalhere = faser(indimage);

        % the threshold for this FASER
        thhere = faserth(faservalhere)*intth;
        
        % neighbour values to look for crash!!
        point = zeros(size(im));
        point(indimage) = 1;
        point = imdilatefast(point,ones(5,5));
        faservalaround = faser(point == 1);
        
        % take away itself and 0, that is no crash
        faservalaround(faservalaround == faservalhere) = [];
        faservalaround(faservalaround == 0) = [];

%         if indimage == 38279
%             faserhere = faser;
%             im(indimage) = max(im(:));
%             faserhere(indimage) = max(faser(:)) + 2;
%             valhere
%             thhere
%             indimage
%             showall(im,faserhere)
%             im = imini;
%         end;
        
        % too close neighbour region
        if length(faservalaround)> 0
%             disp('Too close neighbour, continue')
            continue;
        % too weak, we stop since the array is sorted
        elseif valhere < thhere
            disp('Too weak, break');
            
            break;            
        end;
        
        % assign the value
        bw(indimage) = 1;
        
    end;
    
%     faserth
%     showall(im,bwold,bw,faser)
    
    % stop criterion
    diffbw = bw-bwold;
    numchange = sum(diffbw(:));
    disp(sprintf('Number of dilations: %i, changing pixels: %i',niter,numchange))
    if numchange == 0
        break;
    end;    
end;
% output
cellbw = bw;
cellbw = bw;

% ----------------------------------------------

function [cellbw] = dilatestrong(BW,im,intth,conn)

[M N O] = size(im);

[faser,numcells] = bwlabeln(BW,conn);

%
% Get the means to have a stop criteria
%
borderth = zeros(M,N,O);
% borderth = ones(size(im))*mean(im(:));

% to be similar as below
cellbw = gt(faser,0);   
% radius of dilation between each step
rad = 2;
name = ['ball' int2str(rad)];
load(name);se = getball(ball,rad,O);
niter = 0;
while 1
    niter = niter + 1;
    
    % dilate  to find the border of the cell
    % NOTE: BORDERFASER is the border of the cells, not the image!!!!
    perim = zeros(M,N,O);
    for i = 1 : numcells 
        reghere = eq(faser,i);
        perimhere = imdilatefast(reghere,se) - reghere;
        perim = perim + perimhere;
        
        % too keep track of the segmentation later
        faser(perimhere == 1) = i;
        
%         borderfaser(perimhere == 1) = i;
        % update the thresholds
        borderth(reghere == 1) = mean(im(reghere));
        borderth(perimhere == 1) = mean(im(reghere));
    end;

    % the places of conflict between regions
    perim(perim > 1) = 0;

%     showall(im,perim,cellbw)
        
    % remove the pixels with too low intensities
    perim = perim.*ge(im,intth*borderth);
        
    % add the perims to the cells
    cellbwold = cellbw;
    cellbw(perim == 1) = 1;

    % update the regions in FASER to keep track of the cells
    faser = faser .*cellbw;


    diffBW = cellbw-cellbwold;
    numchange = sum(diffBW(:));
    disp(sprintf('Number of dilations: %i, changing pixels: %i',niter,numchange))
    if numchange == 0
        break;
    end;
    

%     showall(im,perim,cellbw)
%     show(im,1)
%     show(cellbw,2)
%     show(cellbw,6)
%     pause

end

% showall(im,cellbw)

% -----------------------------------------------

function [cellbw] = convexmerge(cellbw)

[M N O] = size(cellbw);

[faser,L] = bwlabeln(cellbw);

for i = 1 : L
    reghere = eq(faser,i);
    
    % get the neighbours
    load ball3;se = getball(ball,3,1);
    dilreghere = imdilate(reghere,se);
    
    neigh = unique(faser(dilreghere == 1));
    neigh(neigh == i) = [];
    neigh(neigh == 0) = [];
    numneigh = length(neigh);
    
    
    % loop over the neighbrous to test convexity
    ratioboth = zeros(numneigh,1);
    for j = 1 : numneigh
        j
        
        % this neighbour
        neighhere = neigh(j);                
        regneigh = eq(faser,neighhere);

        % the combination together
        load ball3;se = getball(ball,3,1);
        regboth = imclose(regneigh + reghere,se);
        
        % the convex are of the combination of this region and the
        % neighbour
        convareaboth = calcconvarea(regboth);        
        areaboth = length(find(regboth));
        
        % the ratio
        ratioboth(j) = areaboth / convareaboth; 
        
    end;
                
    % convex area of thie region
    convareareghere = calcconvarea(reghere);
    areahere = length(find(reghere));
    ratiohere = areahere / convareareghere;
    
    
   numneigh
    show(reghere,2)
    show(faser,3)
    neigh
    ratioboth
    ratiohere
    pause
    
    
end;

% ---------------------------------------------

function [convarea] = calcconvarea(BW)

[M N O] = size(BW);
convarea = 0;
for j = 1 : O
    props = regionprops(double(BW(:,:,j)),'ConvexArea');
    convarea = convarea + props.ConvexArea;        
end;
