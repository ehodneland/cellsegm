function [] = readbioformat(varargin)
% READBIOFORMAT Reads Bio-Formats like .LIF
%
%   READBIOFORMAT(INFILE) reads Bio-Formats like .lif specified in the cell array
%   INFILE using the command line tools. If there is only one file the file 
%   name can be given as a string.
%   SHOWINF and BFCONVERT from http://loci.wisc.edu/bio-formats/downloads.
%   Data is saved as stack1, stack2, and so on, in a folder with subscript 
%   '-matlab'. Data is saved both as .tif and .mat
%
%
%   READBIOFORMAT(INFILE,'mat') save the data only as .mat format and the
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
liffile = varargin{1};
format = 'all';
if nargin == 2
    format = varargin{2};
end;

if ~ismember(format,{'all','mat'})
    error('Wrong option to READBIOFORMAT');
end;

if ~iscell(liffile)
    liffile = {liffile};
end;

for i = 1 : numel(liffile)


    % reading pixel size and number of series
    [folder file ext] = fileparts(liffile{i});
    cdir = pwd;
    if ~isempty(folder)
        cd(folder);
    end;
    namehere = [file ext];
    
    msg = ['Converting ' namehere ' in folder ' folder];
    disp(msg);
    

    tname = tempname;
    msg = ['showinf -nopix ''' namehere ''' >> ' tname ];
    disp(msg);
    system(msg);    
    
    h = [];
    fid = fopen(tname,'r');
    c = 0;
    cchannel = 0;
    cseries = 0;
    cplane = 0;
    while 1
        str = fgetl(fid);
       
        if isequal(str,-1)
            break;
        end
        s = strfind(str,'HardwareSetting|ScannerSettingRecord|dblVoxel');
        if ~isempty(s)
            c = c + 1;
            h(c) = str2double(str(s+50:end));
            
        end;
        
        s = strfind(str,'Series count ');
        if ~isempty(s)
            cseries = cseries + 1;
            nseries(cseries) = str2double(str(s+15:end));
        end;
        
        s = strfind(str,'SizeC =');                
        if ~isempty(s)
            % some times there is a bracket after this, remove!
            s2 = strfind(str,'(');
            if ~isempty(s2)
                str = str(1:s2-1);
            end;
            cchannel = cchannel + 1;
            nch(cchannel) = str2double(str(s+7:end));
        end;
        
        s = strfind(str,'SizeZ =');
        if ~isempty(s)
            cplane = cplane + 1;
            nplane(cplane) = str2double(deblank(str(s+7:end)));
        end;

    end;

    msg = ['Number of series: ' int2str(nseries)];
    disp(msg);
    
    msg = ['Number of channels: ' int2str(nch)];
    disp(msg);

    msg = ['Number of planes: ' int2str(nplane)];
    disp(msg);

    msg = ['Voxel size: ' num2str(h)];
    disp(msg);
    
    % converting data to multiple .tif
    msg = ['bfconvert ''' namehere ''' stack%s.tif'];
    disp(msg);
    system(msg);
    
    
    % move files to matlab folder
    for j = 0 : nseries-1
        in = ['stack' int2str(j) '.tif'];
        targetfolder = [ file '-matlab'];
%         target = [targetfolder '/' 'stack' int2str(j) '.tif'];
        target = fullfile(targetfolder,['stack' int2str(j) '.tif']);
        if ~exist(targetfolder,'dir')
            tmp = mkdir(targetfolder);
        end;
        movefile(in,target);
%         msg = ['mv ''' in ''' ''' target ''''];
%         disp(msg);
%         system(msg);
        
    end;
    
    % change name to matlab numbering starting at 1
    for j = nseries-1 : -1 : 0
        in = fullfile([file '-matlab'], ['stack' int2str(j) '.tif']);
        out = fullfile([file '-matlab'], ['stack' int2str(j+1) '.tif']);
        movefile(in,out);
%         msg = ['mv ''' in '''  ''' out ''''];
%         disp(msg);
%         system(msg);        
    end;
    
    fclose(fid);    
    delete(tname);
    
    % change to microns
    h = h*1e6;
    
    for j = 1 : nseries
    
        % save as .mat
        imload = fullfile([file '-matlab'],['stack' int2str(j) '.tif']);
        msg = ['Converting ' imload ' to .mat format'];
        disp(msg);
    
        % read multiple tif
        im = imreadmulttif(imload,1);
        
        % in case the image is in the fourth position
        im = squeeze(im);
    
        % reorder the data in 4D stack
        im = reordermultipletif(im,nch(j));
    
        if ~isempty(im)
            imsave = fullfile([file '-matlab'],['stack' int2str(j) '.mat']);
            save(imsave,'im','h','-v7.3');
        end;
    end;

    % delete .tif files
    if isequal(format,'mat')
        for j = 1 : nseries
%             imdelete = [file '-matlab' '/' 'stack' int2str(j) '.tif'];
            imdelete = fullfile([file '-matlab'],['stack' int2str(j) '.tif']);
            msg  = ['Deleting ' imdelete];
            disp(msg);
            delete(imdelete);
        end;
    end;

%     % save as mat
%     if isequal(format,'all') || isequal(format,'mat')
%         for j = 1 : nseries
%             imload = [file '-matlab' '/' 'stack' int2str(j) '.tif'];
%             msg = ['Converting ' imload];
%             disp(msg);
%             
%             % read multiple tif
%             im = imreadmulttif(imload);
%             
%             % reorder the data in 4D stack
%             im = reordermultipletif(im,nch(j));
%             
%             imsave = [file '-matlab' '/' 'stack' int2str(j) '.mat'];
%             save(imsave,'im','h');
%             if isequal(format,'mat')
%                 delete(imload);
%             end;
%         end;
%     end;
    cd(cdir);
    
end;
