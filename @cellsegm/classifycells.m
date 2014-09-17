function [cellbw,infocells] = classifycells(varargin)
% CLASSIFY Classify cells from user defined thresholds
%
%   [CELLBW, INFOCELLS] = CLASSIFY(WAT,IM,PRM) Classify the region in WAT from the image
%   IM. PRM is a struct array of settings used in the classification.
%   CELLBW is the cell image after classification and INFOCELLS contains
%   information about the classification.
%
%   PRM.H (default = [0.5 0.5 1.5])
%   Voxel size. 
%
%   1. PRM.METHOD = 'automated' (default)
%   Defines the various modes of classification, with
%   the following options:
%
%   PRM.MINVOLFULL    
%   Minimum volume of 3D cell. Default 3. NB This parameter has the given
%   name of "FULL" since it should be valid for a full 3D volume. If your 
%   data is truly 3D nothing is changed about it and PRM.MINVOL = PRM.MINVOLFULL. If your
%   data is 2D, this parameter is adjusted automatically towards a 2D
%   volume, and a modified parameter PRM.MINVOL is created and used for
%   classification.
%
%   PRM.MAXVOLFULL    
%   Maximum volume of 3D cell. Default 100. Also see description for
%   PRM.MINVOLFULL above for 2D/3D considerations.
%
%   PRM.INTINCELL     
%   Intensity in cell. Default 0.7.
%
%   PRM.INTBORDER     
%   Intensity of border. Default 1.20.
%
%   PRM.CONVEXAREA      
%   Convex area. Default 0.4. 
%
%   PRM.CONVEXPERIM     
%   Convex perimeter. Default 0.35.
%
%   PRM.INTBORDER and PRM.INTINCELL are relative thersholds, compared to the
%   background. A reduced set of thresholds can be defined by a cell array
%   of strings in PRM.PROPNAME. By prm.propname = 'all', all available 
%   thresholds are used.
%
%   2. PRM.METHOD = 'minimacell'
%   The variable MINIMACELL must be defined as
%   [CELLBW, INFOCELLS] = CLASSIFY(...,'MINIMACELL',MINIMACELL) 
%   Classify the region in WAT from the image im using information in
%   MINIMACELL for selecting cells. MINIMIACELL must be a binary image with
%   one white region inside each cell, and otherwise black.
% 
%
%   Ex:
%   cprm.minvolfull = 5;
%   cprm.maxvolfull = 50;
%   [cellbw,infocell] = classifycells(wat,im,cprm);
%
%   See also cellsegm.segmct, cellsegm.classifycells, cellsegm.getminima, 
%   cellsegm.segmsurfwat
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


msg = ['This is ' upper(mfilename) ' for classification of cells and background'];
disp(msg);

wat = varargin{1};
im = varargin{2};
prmin = varargin{3};
% minimacell = ones(size(im));
if nargin == 4
    % nucleus minima?
    minimacell = varargin{4}; 
end;
dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;


%
% Classification parameters
%

% default method
prm.method = 'threshold';

% pixel sizes
prm.h = [0.5 0.5 1.5];

% adjust for 3D
prm.just = 0.9;

% to cut cells
prm.cut = 0;

% minimum and maximum volume
prm.minvolfull = 5;
prm.maxvolfull = 100;

% the intensity threshold relative to the mean of the background
% was at 1.15
prm.intincell = 1.35;

% the intensity threshold on border relative to the mean of the image,
% higher more intensities are required
% was at 1.25
prm.intborder = 1.20;

% the concavity measure, fully concave at 1, set to 0.4-0.5
prm.convexarea = 0.5;
prm.convexperim = 0.35;

% a cell in a marker image
prm.iscellmarker = 1;

% merge input
prm = mergestruct(prm,prmin);
prm.minvolfull = prm.minvolfull*1000;
prm.maxvolfull = prm.maxvolfull*1000;

% voxel volume
prm.voxelvol = prod(prm.h);


% adjust cell volume
[prm.minvol,prm.minvolvox,prm.maxvol,prm.maxvolvox] = cellsegm.cellsize(prm.minvolfull,prm.maxvolfull,prm.h,prm.just,dim(3));

msg = ['Using settings'];
disp(msg);
printstructscreen(prm);

% the integer values to loop over
valwat = unique(wat(wat > 0));

if isequal(prm.method,'threshold')

    msg = ['Using thresholds for classification'];
    disp(msg);
    
    prm.propname = {'volume',...
                    'volume',...
                    'intincell',...
                    'intborder',...
                    'convexarea',...
                    'convexperim'};
    prm.thname = {'minvol',...
                'maxvol',...
                'intincell',...
                'intborder',...
                'convexarea',...
                'convexperim'};
    prm.logic = {'gt','lt','gt','gt','gt','gt'};
    
    
elseif isequal(prm.method,'minimacell')
    msg = ['Using minimacell and volume for classification'];
    disp(msg);
    prm.propname = {'volume',...                    
                    'volume',...
                    'cellmarker'};
    prm.thname = {'minvol',...
                'maxvol',...
                'iscellmarker'};

    prm.logic = {'gt','lt','eq'};
    
else
    error([mfilename ': Wrong option']);
end;

if isequal(prm.method,'threshold') || isequal(prm.method,'minimacell')
    msg = ['Thresholds used for classification:'];
    disp(msg);
    nth = numel(prm.thname);
    p = 15;
    for i = 1 : nth
        name = prm.propname{i};
        logic = prm.logic{i};
        val = prm.(prm.thname{i});
        msg = [makestr(name,p) makestr(logic,p) makestr(num2str(val),p)];
        disp(msg);
    end;
end;

% number of properties to use
prm.npropname = numel(prm.propname);

%
% Mean int background to adjust the thresholds for intensites
%

% if the image is empty
if isempty(valwat)
    cellbw = zeros(size(wat));
    infocells = [];
    return;
end;

% the number of regions
nwat = length(valwat);

% compute background value
if ~isfield(prm,'meanintbck')
    % mean intensity of background (largest region)
    [vol,faser] = bwsize(wat > 0,6);
    [maxvol,ind] = max(vol);
    bck = faser == ind;    
    prm.meanintbck = mean(im(bck));
end;

% get properties
prop = cellsegm.cellprop(im,wat,valwat,prm.propname,prm.h,prm.meanintbck);

% classify
cellbw = zeros(dim);
for i = 1 : nwat

    % this value of WAT
    valwathere = valwat(i);
    
    % this region
    reghere = eq(wat,valwathere);
            
    if isequal(prm.method,'threshold') || isequal(prm.method,'minimacell')
            
        if isequal(prm.method,'minimacell')
            % find overlap to markers
            overlap = minimacell .* reghere;            
            prop.cellmarker(i,1) = ~isempty(find(overlap,1));  
            
        end;
        
        nth = numel(prm.thname);
        dec = NaN(1,nth);
        for j = 1 : nth
            str1 = prop.(prm.propname{j})(i);      
            str2 = prm.(prm.thname{j});
            arg = [prm.logic{j} '(' num2str(str1) ',' num2str(str2) ')'];        
            dec(1,j) = eval(arg);
        end;
        infocells.iscellhere(i,:) = dec;

        % is it a cell?
        iscell = sum(dec) == nth;

        % store
        infocells.iscellfinal(i,1) = 0;
        if iscell
            infocells.iscellfinal(i,1) = 1;
        end;            
        
    end;
       
    if iscell
        cellbw(reghere) = 1;       
        infocells.iscellfinal(i,1) = 1;
    else
        infocells.iscellfinal(i,1) = 0;
    end;

    if iscell == 1
        res = 'cell';
    else
        res = 'background';
    end;
    msg = ['  Classifying object ' int2str(i) ' out of ' int2str(nwat) ' as ' res ', ' num2str(dec)];       
    disp(msg);
    
end;
% to save
infocells.prm = prm;
infocells.prop = prop;
infocells.propname = prm.propname;


%----------------------------------------------------------
