function [] = readbioformatmatlab(varargin)
% READBIOFORMATMATLAB Reads Bio-Formats like .LIF using the bfopen tool for
%   matlab. Can be downloaded from 
%   http://www.openmicroscopy.org/site/support/bio-formats5/developers/matlab-dev.html
%   The folder with the tools must be placed in the matlab path prior to
%   usage.
%
%   READBIOFORMATMATLAB(INFILE) reads Bio-Formats like .lif specified in the cell array
%   INFILE using the command line tools. If there is only one file the file 
%   name can be given as a string.
%   Data is saved as stack1, stack2, and so on, in a folder with subscript 
%   '-matlab'. Data is saved both as .tif and .mat
%
%
%   READBIOFORMATMATLAB(INFILE,'mat') save the data only as .mat format and the
%   .tif files are deleted after use.
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
inputfile = varargin{1};
format = 'all';
if nargin == 2
    format = varargin{2};
end;

if ~ismember(format,{'all','mat'})
    error('Wrong option to READBIOFORMAT');
end;

if ~iscell(inputfile)
    inputfile = {inputfile};
end;

for i = 1 : numel(inputfile)
    % file type
    [A,B,fileformat] = fileparts(inputfile{i});

    % reading pixel size and number of series
    [folder file ext] = fileparts(inputfile{i});
    namehere = [file ext];
    
    msg = ['Converting ' namehere ' in folder ' folder];
    disp(msg);
    
    % read the data
    try
        data = bfopen(inputfile{i});
    catch
        msg = ['Could not open ' inputfile{i}];
        disp(msg);
    end;
    nseries = size(data,1);
    
    msg = ['Number of series: ' int2str(nseries)];
    disp(msg);
        
    for j = 1 : nseries
        
        dim = size(data{j,1}{1,1});
        h = zeros(1,3);        

        metadata = data{j,4};
        try
            a = metadata.getPixelsPhysicalSizeX(0).getValue();
            h(1) = double(a);
        catch
            h(1) = nan;
        end;
        try
            a = metadata.getPixelsPhysicalSizeY(0).getValue();
            h(2) = double(a);
        catch
            h(2) = nan;    
        end;
        try
            a = metadata.getPixelsPhysicalSizeZ(0).getValue();
            h(3) = double(a);
        catch
            h(3) = nan;
        end;
        
        if sum(isnan(h)) == 3
            warning(['Could not read voxel size from series ' int2str(j)]);   
            h = [1,1,1];
        end;
        
        msg = ['Voxel size: ' num2str(h)];
        disp(msg);
        
        
        nimages = size(data{j,1},1);        
        imtif = zeros([dim(1:2),nimages]);
        for k = 1 : nimages            
            imhere = data{j,1}{k,1};
            imtif(:,:,k) = double(imhere);
        end;

        try
            key.X.max = metadata.getPixelsSizeX(0).getValue();
            key.Y.max = metadata.getPixelsSizeY(0).getValue();
            key.Z.max = metadata.getPixelsSizeZ(0).getValue();
            key.C.max = metadata.getPixelsSizeC(0).getValue();
            key.T.max = metadata.getPixelsSizeT(0).getValue();
        catch
            hashtable = data{j,2};
            key.X.max = hashtable.get('Global Image height');
            key.Y.max = hashtable.get('Global Image width');
            key.Z.max = 1;
            key.C.max = 1;
            key.T.max = 1;            
        end;
        
%         str = data{j,1}{1,2};
%         key = readstr(str);

        % reorder data into an array
        msg = ['Number of planes: ' int2str(key.Z.max)];
        disp(msg);

        msg = ['Number of time points: ' int2str(key.T.max)];
        disp(msg);

        msg = ['Number of channels: ' int2str(key.C.max)];
        disp(msg);        

        % save file as mat file
        foldersave = fullfile(folder,[file '-matlab']);
        [a,b,c] = mkdir(foldersave);        

        % save as tif also
        if isequal(format,'all')            
            pathsave = fullfile(foldersave,['stack' int2str(j) '.tif']);
            msg = ['Saving ' pathsave];
            disp(msg);
            imwritemulttif(pathsave,imtif);
        end;

        % reshape image
        im = zeros([dim(1:2),key.Z.max,key.T.max,key.C.max]);
        for k = 1 : nimages
            str = data{j,1}{k,2};           
            key = readstr(str);
            im(:,:,key.Z.current,key.T.current,key.C.current) = imtif(:,:,k);
        end;        
        clear imtif;
        im = squeeze(im);
                        
        % anyway save as mat format
        pathsave = fullfile(foldersave,['stack' int2str(j) '.mat']);
        msg = ['Saving ' pathsave];
        disp(msg);
        save(pathsave,'im','h','-v7.3');
        clear im;
        
    end;
    
    
end;

function [val] = readstr(str)
str = [str ';'];
key = {'C=','Z=','T='};
key2 = {'C?=','Z?=','T?='};
for i = 1 : numel(key)
    keyhere = key{i};
    str2 = str;
    ind = strfind(str,key{i});
    if isempty(ind)
        % in lsm files the '?' can show up, is it a bug??
        ind = strfind(str,key2{i});
    end;
    keyfield = keyhere(1);
    if isempty(ind)
        val.(keyfield).current = 1;
        val.(keyfield).max = 1;
        continue;
    end;
    str2 = str2(ind:end);    
    indslash = strfind(str2,'/');
    indscolon = strfind(str2,';');
    indequal = strfind(str2,'=');
    val.(keyfield).max = str2double(str2(indslash+1:indscolon-1));
    val.(keyfield).current = str2double(str2(indequal+1:indslash-1));               
end;


