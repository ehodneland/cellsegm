function [] = cellsegmentation(varargin)
% CELLSEGMENTATION Cell segmentation of data stored on disk
%
%   CELLSEGMENTATION(FOLDER,STS,STE,PLS,PLE,MINV,MAXV,PRMFILE) 
%
%   Segmentation of cells in the folders specified in FOLDERS which can 
%   either be a single string or a cell array of folders. 
%   Reading all stacks from STS to STE, assuming the name STACK1,STACK2, as 
%   from READBIOFORMAT. The planes are
%   specified as start and stop in PLS and PLE. If PLE = Inf all planes are used.
%
%   MINV and MAXV are the limits for the diameter of a cell in 
%   thousand cubic microns, PRMFILE is the possibility to give a parameter file
%   for segmentation. This is very useful for reproducibility of settings
%   and results. If PRMFILE = [] no parameter file is loaded.
%
%   CELLSEGMENTATION(...,'prm',PRM)
%   adds a parameter struct with command-line defined parameters for the
%   segmentation process. The command line parameters have priority to the
%   parameter file settings.
%
%   The segmentation method must be chosen by setting PRM.METHOD.
%   Choices for PRM.METHOD:
% 
%   PRM.METHOD = 'segmsurf' (default) : Calling SEGMSURF for segmentation of 
%   surface stained cells. The marker method should be selected via 
%   PRM.SEGMSURF.GETMINIMA.METOD (see below).
%
%   PRM.METHOD = 'segmct': Calling SEGMCT for cytoplasmically stained cells. Can also
%   be used for segmentation of nuclei.
%
%   PRM.SEGMSTART (default = 1)
%   Starts segmentation at a higher plane than the bottom
%   plane. This is useful is the signal is very weak at the bottom, which
%   can mislead the segmentation.
%
%   PRM.SEGMCH (default = 1) 
%   Defines the number of the segmentation channel in the 4D image stack. 
%
%   PRM.NUCLEUSCH (default = [])
%   Defines the number of the nucleus channel in the 4D image stack. 
%   
%   PRM.SUBTRACT (default = 1)
%   Defines whether the nucleus channel is subtracted (1) or not (0) from the
%   surface stain (in case of cross talk). 
%
%   PRM.FORMAT (default = 'mat')
%   Defines the image format to be loaded and saved after processing. Can
%   be either 'mat' or 'tif'. Remember to also set
%   prm.nchannel when using option 'tif'. 
%  
%   PRM.NCHANNEL (default = 1)
%   Defines the number of channels in your data. It is only used for reading the
%   multiple tif files (with prm.format = 'tif') when the 4D image stack is 
%   sorted after each other in a 3D manner. This variable is used to know 
%   the distribution between number of channels and number of planes. 
%
%   The method for selecting markers should be specified.
%
%   PRM.SEGMSURF.GETMINIMA.METHOD = 'MANUAL': Attempts to load manually defined markers
%   which must be of the format STACK1-MASK.MAT, STACK2-MASK.MAT and in the
%   same folder as the data. The manual markers must be binary images with
%   ones at the markers and zero elsewhere, and the image must have same
%   dimension as the segmentation image. There must be two variables, MINIMA 
%   with all markers and MINIMACELL with only cell markers. 
%
%   PRM.SEGMSURF.GETMINIMA.METHOD = 'NUCLEUS': Using nucleus markers. This option
%   requires PRM.NUCLEUSCH to be defined.
%
%   PRM.SEGMSURF.GETMINIMA.METHOD = 'AUTOMATED' (default): Automated 
%   detection of cell markers from the surface staining channel.
%   The results from CELLSEGMENTATION can be visualized by using VIEWSEGM
%
%
%   Example
%   -------
%
%   examples/surfstain_and_nucleus_cellsegmentation_3D.m
%
%
%   Example (Will not run, no data is supporting the commands)
%   -------
%
%   folder = {'folder1','folder2','folder3'};    
%   cellsegmentation(folder,1,10,1,20,1,1,100,'myprmfile')
% 
%
%   See also viewsegm
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


prm = [];
prmincmd = [];
folder = varargin{1};
prm.stackstart = varargin{2};
prm.stackstop = varargin{3};
prm.planestart = varargin{4};
prm.planestop = varargin{5};
prm.minvolfull = varargin{6};
prm.maxvolfull = varargin{7};
prm.prmfile = varargin{8};
for i = 9 : 2 : nargin
    varname = varargin{i};
    switch varname
        case 'prm'
            prmincmd = varargin{i+1};
        otherwise
            error(['Wrong option to ' mfilename]);
    end;
end;

prm.stacks = prm.stackstart : prm.stackstop;
prm.nstacks = numel(prm.stacks);
prm.format = 'mat';
prm.segmstart = 1;
prm.subtract = 1;
prm.segmch = 1;
prm.nucleusch = [];
prm.method = 'segmsurf';
prm.h = [0.5 0.5 1.5];

% method for segmentation
prm.segmsurf.getminima.method = 'automated';

% requires Jacket for matlab!
prm.gpu = 0;

% the number of folders and make the folder name into a cell array
if iscell(folder)
    nfolder = length(folder);
else
    nfolder = 1;
    f = folder;
    clear folder;
    folder{1} = f;    
end;

%
% Load the parameter file if one is given
%

prmfile = prm.prmfile;
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

msg = ['This is ' mfilename ' using settings'];
disp(msg);
printstructscreen(prm);

% from command line, these have priority
prm = mergeinputpar(prm,prmincmd);
if numel(prm.segmstart) == 1
    prm.segmstart = prm.segmstart*ones(prm.nstacks,nfolder);
end;

% in prmsegm we only want scalar settings, for this image only
prm.segmstart(prm.segmstart == 0) = NaN;
prm.segmstart = prm.segmstart;
% prmsegm = rmfield(prm,'segmstart');

% % fill out the empty slots of NaN's
% d = size(prm.startplane);
% s = prm.startplane;
% prm.startplane = NaN(prm.stacks(end),nfolder);
% prm.startplane(1:d(1),1:d(2)) = s;
% d = size(prm.stopplane);
% s = prm.stopplane;
% prm.stopplane = NaN(prm.stacks(end),nfolder);
% prm.stopplane(1:d(1),1:d(2)) = s;

%
% Segmentation, loop over folders
%
for i = 1 : nfolder

    % segmentation in this folder
    segmentfolder(folder{i},prm);

end;%for


%-------------------------------------


function [] = segmentfolder(folderhere,prm)

nstacks = length(prm.stacks);
% loop over all images
for i = 1 : nstacks
    
    % the stack here
    stack = prm.stacks(i);
    name = [folderhere '/'  'stack' int2str(stack)];    
    
    
    % get the image
    [im,imsegm,imnucl,prm] = getim(name,prm,i);    
    
    if isempty(im)
        msg = ['Could not load ' name ', continuing'];
        disp(msg);
        continue;
    end;
        
    dim = size(im);
    if numel(dim) == 2
        dim = [dim 1];
    end;

    if dim(1) == 0 || dim(2) == 0 || dim(3) == 0        
        error('Wrong input, empty image');
        continue;
    end;

    if prm.segmstarti > dim(3) 
        msg = ['Segmentation plane is above stack, continuing'];
        disp(msg);
        continue;
    end;
        
    % so we dont need to repeat plane for this image
    imsegminput = imsegm;

    
    % Segmentation of surface stained cells
    if isequal(prm.method,'segmsurf')       
            
        % Segmentation of surface stained cells
        [cellbw,wat,imsegm,minima,minimacell,info] = segmsurfhere(imsegm,imnucl,prm);

    % segmentation of cytoplasmically stained cells, or nuclei only
    elseif isequal(prm.method,'segmct')
               
        % dont take the lowest planes?
        if prm.segmstart > 1            
            msg = ['Segmentation from plane ' int2str(prm.segmstart)];
            disp(msg);
            imsegm = imsegm(:,:,prm.segmstart:end);
        end;

        % segmentation     
        try
            prmin = prm.segmct;
        catch
            
        end;
           
        prmin.h = prm.h;
        [cellbw,wat] = cellsegm.segmct(imsegm,prm.minvolfull,prm.maxvolfull,'prm',prmin);        
        
        % no minimum is created
        minima = [];
    else
        error([mfilename ': Wrong option in METHOD']);
    end;
            
    
    if prm.segmstarti > 1
        msg = ['Adding lower planes'];
        disp(msg);
                
        % fix the lower planes              
        cellbw = repeatplane(cellbw,prm.segmstarti);
        wat = repeatplane(wat,prm.segmstarti);                
        imsegm = imsegminput;
        m = zeros(dim(1:3));
        m(:,:,prm.segmstarti:end) = minima;
        minima = m;        
        if ~isempty(minimacell)
            m = zeros(dim(1:3));
            m(:,:,prm.segmstarti:end) = minimacell;
            minimacell = m;
        end;
    end;
    
    wat = double(wat);    
  
    if isequal(prm.format,'mat') || isequal(prm.format,'all')
        % save variables
        name = [folderhere '/' 'stack' int2str(stack) '-segm.mat' ];
        save(name,'cellbw','wat','im','imsegm','minima','minimacell','info','-V6')    
        msg = ['Saving ' name];
        disp(msg);
    end;
    if isequal(prm.format,'tif')|| isequal(prm.format,'all')
        name = [folderhere '/' 'stack' int2str(stack) '-cellbw.tif' ];
        imwritemulttif(name,cellbw);
        name = [folderhere '/' 'stack' int2str(stack) '-wat.tif' ];
        imwritemulttif(name,wat);
        name = [folderhere '/' 'stack' int2str(stack) '-imsegm.tif' ];
        imwritemulttif(name,imsegm);
        name = [folderhere '/' 'stack' int2str(stack) '-minima.tif' ];
        imwritemulttif(name,minima);
        if ~isempty(minimacell)
            name = [folderhere '/' 'stack' int2str(stack) '-minimacell.tif' ];
            imwritemulttif(name,minimacell);
        end;
        
    end;
end%for


%---------------------------------

function [cellbw,wat,imsegm,minima,minimacell,info] = segmsurfhere(imsegm,imnucl,prm)

dim = size(imsegm);
if numel(dim) == 2
    dim = [dim 1];
end;

% nucleus channel for markers
if isequal(prm.segmsurf.getminima.method,'nucleus') 
    if prm.subtract  
        
        % subtract nucleus from WGA
        if dim(3) == 1
            filt = fspecial('gaussian',3,1);
            b1 = imfilter(imsegm,filt,'replicate');
            b2 = imfilter(imnucl,filt,'replicate');
        else
            b1 = smooth3(imsegm,'gaussian',[3 3 3]);
            b2 = smooth3(imnucl,'gaussian',[3 3 3]);
        end;

        % for watershed we use subtracted image                
        load ball2;se = getball(ball,2,1);            
        for j = 0.1:0.1:1.5
            % How much to subtract? Try iterating until values get
            % negative or loop is done
            imsegmtemp = b1 - j*b2;
            vol = imsegmtemp < 0;
            vol = imopen(vol,se);                    
            vol = sum(vol(:));
            % more than just noise?
            if vol > 0.10*numel(imsegm)                    
                break;
            end;
            imsegm = imsegmtemp;                
        end;

        clear imsegmtemp
        imsegm(imsegm < 0) = 0;
    end;
end;       


% cut the image for segmentation
if prm.segmstarti > dim(3)
    error([mfilename ': Segmentation starts outside image volume']);
end;
imsegm = imsegm(:,:,prm.segmstarti:end);

% Manual minima given? Load if they are
if isequal(prm.segmsurf.getminima.method,'minimacell')
    % load manual minima
    minimaload = [folderhere '/' 'stack' int2str(stack) '-mask.mat'];
    D = load(minimaload);            
    % we expect to find two variables: minima and minimacell
    minima = D.minima;
    minimacell = D.minimacell;

    % check for manual labels in the lower planes; we include them
    % if they exist! Thus we may change prmsegm.segmstart!!           
    for j = 1 : prm.segmstarti                    
       imhere = minima(:,:,j);
       if ~isempty(find(imhere,1));
           prm.segmstarti = j;
           break;
       end;
    end;

    % cut the minima
    minima = minima(:,:,prm.segmstarti:end,:);
    minimacell = minimacell(:,:,prm.segmstarti:end,:);

    % check for dimension of loaded minima
    dimhere3 = size(minima);
    if numel(dimhere) == 2
        dimhere3 = [dimhere3 1];
    end;
    if ~isequal(dimhere3,dim3)
        error('Wrong dimension of images');
    end;

    % fill holes if some of the objects are not filled
    for j = 1 : dim(3)
        minima(:,:,j) = imfill(minima(:,:,j),'holes');
        minimacell(:,:,j) = imfill(minimacell(:,:,j),'holes');
    end;
elseif isequal(prm.segmsurf.getminima.method,'nucleus')
    % cut nucleus image
    imnucl = imnucl(:,:,prm.segmstarti:end);
end;

msg = ['Segmentation from plane ' int2str(prm.segmstarti)];
disp(msg);

%
% Segmentation    
%

try
    prmin = prm.segmsurf;
catch
    
end;
prmin.h = prm.h;
prmin.gpu = prm.gpu;
if isequal(prmin.getminima.method,'minimacell')
    [cellbw,wat,imsegm,minima,minimacell,info] = ...
        cellsegm.segmsurf(imsegm,prm.minvolfull,prm.maxvolfull,'prm',prmin,'minima',minima,'minimacell',minimacell);                        
    clear temp;
elseif isequal(prmin.getminima.method,'nucleus')
    [cellbw,wat,imsegm,minima,minimacell,info] = ...
        cellsegm.segmsurf(imsegm,prm.minvolfull,prm.maxvolfull,'prm',prmin,'imnucleus',imnucl);                                        
elseif isequal(prmin.getminima.method,'automated')                  
    [cellbw,wat,imsegm,minima,minimacell,info] = ...
        cellsegm.segmsurf(imsegm,prm.minvolfull,prm.maxvolfull,'prm',prmin);                        
else
    error([mfilename ': Wrong minima method in GETMINIMAMETHOD']);
end;
        

%---------------------------------

function [imnew] = repeatplane(im,start)

dim = size(im);
add = start-1;
dimnew = dim;
dimnew(3) = dimnew(3) + add;
imnew = zeros(dimnew);
imnew(:,:,start:end) = im;
for i = 1 : start - 1
    imnew(:,:,i) = im(:,:,1);
end;

%-------------------------------------------------

function [im,imsegm,imnucl,prm] = getim(name,prm,i)

try

    if isequal(prm.format,'mat')
        D = load(name);
        
    elseif isequal(prm.format,'tif')
        namehere = [name '.tif'];        
        D.im = imreadmulttif(namehere);
        D.im = reordermultipletif(D.im,prm.nchannel);
    else
        error('Wrong format');
    end;
    disp(sprintf('Segmentation of %s using channel %i',name,prm.segmch));
catch ME1
    disp(sprintf('Unable to load %s',name));
    disp(ME1.message)
    im = [];
    imsegm = [];
    imnucl = [];
    return;
end;
im = D.im;

try
    % in microns
    prm.h = D.h;
catch
    % use default
    warning('Could not read proper voxel size, check your data');
end;

clear D;
dim = size(im);

% start and stop planes in this stack
if isequal(prm.planestart,'automatic')
    val = zeros(dim(3),1);
    for i = 1 : dim(3)
        imhere = im(:,:,i,prm.segmch);
        val(i,1) = mean(abs(imhere(:)));
    end;

    
    % start at maximum signal intensity, assuming thats the bottom of the
    % cells
    [~,ind] = max(val);
    planestart = ind-1;    
end;
prm.planestarti = max(planestart,1);
prm.planestopi = min(prm.planestop,dim(3));    
msg = ['Starting the stack at plane ' int2str(prm.planestarti)];
disp(msg);
msg = ['Stopping the stack at plane ' int2str(prm.planestopi)];
disp(msg);

% this becomes all planes in stack
prm.planei = prm.planestarti:prm.planestopi;

% all channels
im = im(:,:,prm.planei,:);


% the image to segment
imsegm = im(:,:,:,prm.segmch);

% the nucleus image
imnucl = [];
if ~isempty(prm.nucleusch)
    imnucl = im(:,:,:,prm.nucleusch);
end;

% where to start and stop segmentation in this stack
prm.segmstarti = max(1,prm.segmstart(i));

 
