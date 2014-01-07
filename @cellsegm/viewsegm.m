function [] = viewsegm(varargin)
% VIEWSEGM Visual inspection of the segmentation from SEGMCELL
%
%   VIEWSEGM(START,STOP,CH1,CH2) starts at STACK[START]-SEGM.MAT
%   for visualization of the segmented cells. CH1 and CH2 are numbers 
%   defining the image channels for visualization. There are several 
%   option available on the GUI:
%
%   Up            : Up in the stack
%   Down          : Down in the stack
%   Next          : Next stack
%   Previous      : Previous stack
%   Classification: By clicking on the cells one can see why a cell was
%                   accepted or rejected.
%   Frame         : Manually enter in the pop-up window the plane number
%   Quit          : Quit VIEWSEGM
%
%   VIEWSEGM(START,STOP,CH1,CH2,'raw') reads the raw data in stack1.mat,
%   stack2.mat etc for visualization
%
%   See also segmcell, readbioformat
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

clear choice ch;
global choice 
global ch;

start = varargin{1};
stop = varargin{2};
prm.vis(1).ch = varargin{3};
prm.vis(2).ch = varargin{4};
prm.nch = 0;
prm.option = 'segm';
if nargin < 3
    error('Tou must specify an image channel');
end;
for j = 3 : nargin
    v = varargin{j};
    if isnumeric(v)
        prm.vis(j-2).ch = varargin{j};
        prm.nch = prm.nch + 1;
    elseif ischar(v)
        prm.option = v;
    end;
end;


% fig.pos = [400 50 500 500];
% fig.handle = figure('Position',fig.pos);
% control.pos = [100 50 300 500];
% control.handle = figure('Position',control.pos);

% image number
i = start;
iold = inf;
finish = 0;
plane = 1;

% get username
prm.username = char(java.lang.System.getProperty('user.name'));


% % change size if necessary
% p = get(fig.handle,'OuterPosition');

handle.fig.w = 400*prm.nch;
handle.fig.h = 700;
handle.fig.b = 100;
handle.fig.l = 400;
handle.fig.pos = [handle.fig.l handle.fig.b handle.fig.w handle.fig.h];
handle.fig.handle = figure('Position',handle.fig.pos,'KeyReleaseFcn',@cb);


% set(fig.handle,'Position',fig.pos);
% pause

handle.control.w = 400;
handle.control.h = handle.fig.h;    
handle.control.pos = [handle.fig.pos(1)-handle.control.w handle.fig.b handle.control.w handle.control.h];
handle.control.handle = figure('Position',handle.control.pos);
% set(control.handle,'Position',control.pos);
% set(fig.handle,'KeyPressFcn',@(h_obj,evt)disp(evt.Key));

handle.checkbox.w = 50;
handle.editfield.w = 50;
% figure(fig.handle);
set(handle.fig.handle,'Toolbar','figure');
set(handle.control.handle,'Toolbar','figure');

% margins and settings
handle.button.h = 20;
handle.button.w = 100;
handle.button.left = 10;

% set image scale
for j = 1 : prm.nch
    for k = 1 : 2
        prm.imscale(j).sc{k} = [];
    end;
end;

left = handle.button.left;
b = handle.control.h - 2*handle.button.h;
w = handle.button.w;
h = handle.button.h;
handle.up.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Up', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''up'';uiresume(gcbf);');    
b = b - h - 5;
handle.down.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Down', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''down'';uiresume(gcbf);');    
b = b - h - 5;
handle.next.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Next', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''next'';uiresume(gcbf);');    
b = b - h - 5;
handle.previous.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Previous', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''previous'';uiresume(gcbf);');    
b = b - h - 5;
handle.frame.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Frame number', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''framenum'';uiresume(gcbf);');    
b = b - h - 5;
handle.classification.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Classification', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''classification'';uiresume(gcbf);');    
% Delete/add cells permanently!
b = b - h - 5;
handle.cellstatus.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Cell status', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''cellstatus'';uiresume(gcbf);');    


% Image scale
b = b - h - 5;
handle.imscale.header(1).pos = [left b w h];
handle.imscale.header(1).h = uicontrol('Parent',handle.control.handle,'Style','Text','String','', ...
     'Position',handle.imscale.header(1).pos);   
dw = w;
handle.imscale.header(2).pos = [left+dw b handle.checkbox.w h];
handle.imscale.header(2).h = uicontrol('Parent',handle.control.handle,'Style','Text','String','', ...
     'Position',handle.imscale.header(2).pos);   
dw = dw + handle.checkbox.w;
handle.imscale.header(3).pos = [left+dw b handle.editfield.w handle.button.h];
handle.imscale.header(3).h = uicontrol('Parent',handle.control.handle,'Style','Text','String','Low', ...
     'Position',handle.imscale.header(3).pos);   
dw = dw + handle.editfield.w;
handle.imscale.header(4).pos = [left+dw b handle.editfield.w handle.button.h];
handle.imscale.header(4).h = uicontrol('Parent',handle.control.handle,'Style','Text','String','High', ...
     'Position',handle.imscale.header(4).pos);   
 
% controlling scale image 
c = 9;
for j = 1 : prm.nch
    b = b - h - 5;    
    handle.imscale.title(j).pos = [left b w h];
    handle.imscale.title(j).handle = uicontrol('Parent',handle.control.handle,'Style','Text','String',['Scale image ' int2str(j)], ...
         'Position',handle.imscale.title(j).pos);   
    dw = w;
    c = c + 1;
    handle.imscale.checkbox(j).pos = [left+dw b handle.checkbox.w h];
    handle.imscale.checkbox(j).handle = uicontrol('Parent',handle.control.handle,'Style','Checkbox', ...
         'Position',handle.imscale.checkbox(j).pos,'CallBack', ...
        ['global choice;choice = 1;uiresume(gcbf);']);    
    dw = dw + handle.checkbox.w;
    c = c + 1;
    handle.imscale.checkbox(j).editfield(1).pos = [left+dw b handle.editfield.w h];    
    handle.imscale.checkbox(j).editfield(1).handle = uicontrol('Parent',handle.control.handle,'Style','edit', ...
        'Position',handle.imscale.checkbox(j).editfield(1).pos,'CallBack', ...
        ['global choice;choice = 1;uiresume(gcbf);']); 
    
    dw = dw + handle.checkbox.w;
    c = c + 1;
    handle.imscale.checkbox(j).editfield(2).pos = [left+dw b handle.editfield.w h];
    handle.imscale.checkbox(j).editfield(2).handle = uicontrol('Parent',handle.control.handle,'Style','edit', ...
        'Position',handle.imscale.checkbox(j).editfield(2).pos,'CallBack', ...
        ['global choice;choice = 1;uiresume(gcbf);']);    
    
    
end;



% Remove small cells
% NB not permanently, only visually
b = b - h - 5;
dw = 0;
handle.celldim.button.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Vol(mcm^3)', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''volume'';uiresume(gcbf);');    
dw = dw + w;
dw = dw + handle.checkbox.w;
handle.celldim.editfield(1).handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','100', ...
    'Position',[left+dw b handle.editfield.w h]);    
dw = dw + handle.editfield.w;
handle.celldim.editfield(2).handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','Inf', ...
    'Position',[left+dw b handle.editfield.w h]);    


% Remove boundary cells
% NB not permanently, only visually
b = b - h - 5;
dw = 0;
handle.boundary.button.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Boundaries', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''boundary'';uiresume(gcbf);');    
dw = dw + w + handle.checkbox.w;
handle.boundary.editfield.handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','0.05', ...
    'Position',[left+dw b handle.editfield.w h]);    

% draw line
b = b - h - 5;
dw = 0;
handle.draw.button.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Draw', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''draw'';uiresume(gcbf);');    
dw = dw + w;
handle.draw.checkbox.pos = [left+dw b handle.checkbox.w h];
handle.draw.checkbox.handle = uicontrol('Parent',handle.control.handle,'Style','Checkbox', ...
         'Position',handle.draw.checkbox.pos,'CallBack', ...
        'global choice;choice = ''drawim'';uiresume(gcbf);');    
dw = dw + handle.checkbox.w;
handle.draw.editfield.pos = [left+dw b handle.editfield.w h];
handle.draw.editfield.handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','1', ...
    'Position',handle.draw.editfield.pos);    

% Quit
b = b - h - 5;
handle.quit.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Quit', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = ''quit'';uiresume(gcbf);');        
b = b - h - 5;
handle.control.title.pos = [left b 4*w h];
handle.control.title.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','Message field', ...
     'Position',handle.control.title.pos);   


% to save settings
prm.pathsettings = fullfile(['~' prm.username],'.cellsegm');
if ~isdir(prm.pathsettings)
    mkdir(prm.pathsettings);
end;
prm.pathsettings = fullfile(prm.pathsettings,'settings.mat');


% load last used settings
if exist(prm.pathsettings,'file')
    msg = ['Loading previous settings'];
    disp(msg);
    D = load(prm.pathsettings);    
    % list of handles
    for i = 1 : prm.nch
        for j = 1 : numel(D.handle.vis)
            % we have the same channel as in the settings file?
            if D.handle.vis(j).ch == prm.vis(i).ch
                try
                    % set the two edit fields
                    v = D.handle.imscale.checkbox(i).editfield(1).value;
                    set(handle.imscale.checkbox(i).editfield(1).handle,'String',v);
                catch
                    
                end;
                try
                    v = D.handle.imscale.checkbox(i).editfield(2).value;
                    set(handle.imscale.checkbox(i).editfield(2).handle,'String',v);                    
                catch
                    
                end;
                try
                    % set the checkbox
                    v = D.handle.imscale.checkbox(i).value;
                    set(handle.imscale.checkbox(i).handle,'Value',v);                    
                catch
                    
                end;
            end;
        end;
    end;
end;


while 1
    
   
    if isequal(prm.option,'segm')
         name = [ 'stack' int2str(i) '-segm.mat'];
     elseif isequal(prm.option,'raw')
         name = [ 'stack' int2str(i) '.mat'];
     end;
    

    if i > stop
        close(fig.handle);
        close(control.handle);
        return;
    end;
      
    % load new image if necessary
    if ne(i,iold)
        
        try
            D = load(name);
            disp(sprintf('Loading %s',name));
            pathstack = name;
            im = double(D.im);
            dim = size(D.im);
            if numel(dim) == 2
                dim = [dim 1];
            end;
            dim3 = dim(1:3);
            if isequal(prm.option,'segm')
                cellbw = D.cellbw;
                wat = D.wat;
                minima = D.minima;
                imsegm = D.imsegm;
                info = D.info;
%                 prmin = D.info.classifycells;
            elseif isequal(prm.option,'raw')
                cellbw = zeros(dim3);
                wat = zeros(dim3);
                minima = zeros(dim3);
                imsegm = zeros(dim3);
            end;
            prm.h = info.prm.h;
            
            % for visualization
            cellbwvis = cellbw;
            
            % fore removing cells from boundary and due to size
            boundarybw = zeros(dim);
            volbw = zeros(dim);
            
            % we have not deleted any cells
            prm.cellmod = 0;
                       
            % we have not taken away small or large cells
            prm.volbw = 0;

            % we have note removed or0 added any cells
            prm.boundarybw = 0;

            % empty draw fig handle
            handle.draw.handle = [];
            
            % number of lines
            cline = 0;
            linedata.x = [];
            linedata.y = [];
            
            prm.voxelvol = prod(info.prm.h);
            clear D;
            msg = ['Image dimension: ', num2str(dim)];
            disp(msg);

                   
        catch
            msg = ['Could not load ' name];
            disp(msg);
            if isequal(choice,'previous')
                i = i - 1;
            else isequal(choice,'next')
                i = i + 1;
            end;
            continue;
        end;
        
    end;
    iold = i;        

    
         
   % figure size in relation to screen
    scrsz = get(0,'ScreenSize');
    scrsz = scrsz(3:4);
    scrsz = scrsz-150;
    maxw = scrsz(1);
    maxh = scrsz(2);
    ratio = dim(2)/dim(1);
    if ratio > 1
        w = maxw;
        h = w/ratio;
    elseif ratio <= 1
        h = maxh;
        w = h*ratio;
    end;        
    handle.fig.w = w;
    handle.fig.h = h;

    % show images
%     prm.fig = fig;
    prm.plane = plane;
    prm.dim = dim;   
%     prm.draw = draw;
    handle = showimagecells(im,linedata,cellbwvis,wat,minima,prm,handle);

    
    % wait for user
    flagui = 0;
%     k = waitforbuttonpress;
    
%     curchar=uint8(get(gcf,'CurrentCharacter'));
%     if isequal(curchar,30);
%         plane = plane + 1;
%         choice = 'up';
%         uiresume(gcbf);
%     elseif isequal(curchar,31);
%         plane = plane - 1;
%         choice = 'down';
%         uiresume(gcbf);
%     end;
%     choice
%     waitfor(fig.handle,'CurrentCharacter');
    
%     if k == 0
        uiwait(handle.control.handle)      
        flagui = 1;
%     end;
    
    % no uicontrol input
    if flagui == 0
        % shortcuts from keyboard
        if isequal(ch,'uparrow');
            choice = 'up';
        elseif isequal(ch,'downarrow');
            choice = 'down';
        elseif isequal(ch,'rightarrow');
            choice = 'next';
        elseif isequal(ch,'leftarrow');
            choice = 'previous';
        end;
    end;
   
    
    if isequal(choice,'up')
        plane = min(plane + 1,dim(3));
    elseif isequal(choice,'down')
        plane = max(plane - 1,1);
    elseif isequal(choice,'next')
        i = i + 1;
        if prm.cellmod == 1
            savestack(pathstack,info,cellbw);        
        end;
        if ~isempty(linedata.x)
            saveline(pathstack,linedata,prm);
        end;
    elseif isequal(choice,'previous')        
        i = i - 1;        
        if prm.cellmod == 1        
            savestack(pathstack,info,cellbw);
        end
        if ~isempty(linedata.x)
            saveline(pathstack,linedata,prm);
        end;
    elseif isequal(choice,'framenum')
        planenew = inputdlg('Enter plane number: ');
        if ~isempty(planenew),
            plane = str2double(planenew);
        end;
    elseif isequal(choice,'classification')
        % user defined input
        figure(handle.fig.handle);
        [y x] = ginput(1);
        x = round(x);
        y = round(y);
        val = wat(x,y);
        a = info.classifycells;
        nth = numel(a.propname);

        A = {'Classificator','Value','Threshold','Cell'};
        for j = 1 : nth
            name = a.propname{j};
            A{j+1,1} = name;
            A{j+1,2} = num2str(a.prop.(name)(val));
            A{j+1,3} = num2str(double(a.prm.(a.prm.thname{j})));            
            A{j+1,4} = num2str(double(a.iscellhere(val,j)));            
        end;
%         A{nth+2,1} = 'All classificators must be 1 to become a cell';
        
         
        % print the string to the figure 
        b = handle.control.title.pos(2) - 2*handle.button.h;
        h = handle.button.h;

        w = handle.button.w;        
        for j = 1 : size(A,1)-1
            
            left = handle.button.left;
            for k = 1 : size(A,2)                
                handle.control.message.handle(j,k) = uicontrol('Parent',handle.control.handle,'Style','Text','String',A{j,k}, ...
                'Position',[left b w h]);   
                left = left + handle.button.w;           
            end;
            b = b - handle.button.h;
        end;
%         % last message string
%         left = handle.button.left;
%         w = handle.button.h*4;
%         handle.control.message.handle(j+1,1) = uicontrol('Parent',handle.control.handle,'Style','Text','String',A{nth+2,1},'Position',[left b w h]);
        
    elseif isequal(choice,'cellstatus')
               
        % user defined input
        figure(handle.fig.handle);
        [y x] = ginput(1);
        x = round(x);
        y = round(y);

        val = wat(x,y);
        
        % this is a boundary
        if val == 0
            continue;
        end;
        
        % new cell status, we must save the stack again
        prm.cellmod = 1;

        % remove/add the cell
        reg = wat == val;
        flag = cellbw(x,y,plane) == 1;
        if flag == 1
            cellbw(reg) = 0;
            cellbwvis(reg) = 0;
        elseif flag == 0
            cellbw(reg) = 1;
            cellbwvis(reg) = 1;
        end;
        % add a column with manual removal in the info struct
        if ~ismember('manual',info.classifycells.propname) 
            info.classifycells.propname = [info.classifycells.propname {'manual'}];
            ncells = size(info.classifycells.iscellfinal,1);
            info.classifycells.iscellhere = [info.classifycells.iscellhere ones(ncells,1)];
            info.classifycells.prm.thname = [info.classifycells.prm.thname {'ismanual'}];
            info.classifycells.prm.propname = [info.classifycells.prm.propname {'manual'}];
            info.classifycells.prm.logic = [info.classifycells.prm.propname {'eq'}];
            info.classifycells.prm.ismanual = 1;
            info.classifycells.prop.manual = zeros(ncells,1);            
        end;
        if flag == 1
            info.classifycells.iscellhere(val,end) = 0;
            info.classifycells.iscellfinal(val) = 0;
            info.classifycells.prop.manual(val) = 0;            
        elseif flag == 0
            info.classifycells.iscellhere(val,end) = 1;
            info.classifycells.iscellfinal(val) = 1;            
            info.classifycells.prop.manual(val) = 1;            
        end;
        
%         show(cellbw(:,:,1),3)
    elseif isequal(choice,'volume')
        
        % to use in boundary
        volbw = zeros(dim);
            
        if prm.volbw == 1
            cellbwvis = cellbw;
            prm.volbw = 0;
        elseif prm.volbw == 0;
            
            % small cells
            v = get(handle.celldim.editfield(1).handle,'String');
            if ~isempty(v)
                try
                    v = str2double(v);                    
                catch
                    set(handle.celldim.editfield(1).handle,'');
                end;
            end; 
            if ~isinf(v)
                a = cellbw;
                v = round(v/prm.voxelvol);            
                msg = ['Removing small cells'];
                disp(msg);                
                cellbwvis = bwareaopenrange(a,v,6);
                d = a - cellbwvis;
                volbw(d == 1) = 1;
            end;

            % large cells
            v = get(handle.celldim.editfield(2).handle,'String');
            if ~isempty(v)
                try
                    v = str2double(v);                    
                catch
                    set(handle.celldim.editfield(2).handle,'');
                end;
            end;        
            if ~isinf(v)
                a = cellbwvis;
                v = round(v/prm.voxelvol);
                msg = ['Removing large cells'];
                disp(msg);
                cellbwvis = cellbwvis - bwareaopenrange(a,v,6);
                d = a - cellbwvis;
                volbw(d == 1) = 1;
            end;
            prm.volbw = 1;
        end;        
        
        % add boundary information
        cellbwvis(boundarybw == 1) = 0;
        
    elseif isequal(choice,'boundary')
        
        % to use in size removal
        boundarybw = zeros(dim);
        
        if prm.boundarybw == 0                    
            
            a = cellbwvis;
            v = get(handle.boundary.editfield.handle,'String');
            if ~isempty(v)
                try
                    v = str2double(v);                    
                catch
                    set(handle.boundary.editfield.handle,'');
                end;
            end;    
            cellbwvis = imclearborderth(cellbw,v);   
            d = a  - cellbwvis;
            % the boundary cells we keep if we want to put them back
            boundarybw(d == 1) = 1;
            prm.boundarybw = 1;
        elseif prm.boundarybw == 1
            cellbwvis = cellbw;
            prm.boundarybw = 0;
        end;
        % add size information
        cellbwvis(volbw == 1) = 0;
        
%     elseif isequal(choice,'drawch')
%         v = get(draw.checkbox.handle,'String');
%         if v == 1
%             showimagedraw(im,draw,plane);
%         end;
        
    elseif isequal(choice,'draw')        
        figure(handle.draw.handle);
        [x, y] = ginput(2);
        cline = cline + 1;
        linedata.x{cline} = round(x);
        linedata.y{cline} = round(y);            
        linedata.z{cline} = [plane; plane];            
    elseif isequal(choice,'quit')
        if prm.cellmod == 1
            savestack(pathstack,info,cellbw);        
        end;
        if ~isempty(linedata.x)
            saveline(pathstack,linedata,prm);
        end;
        
        finish = 1;
    end;

        
    msg = ['Plane ' int2str(plane)];
    disp(msg);
  
    if finish == 1
        close(handle.fig.handle);
        close(handle.control.handle);
        if ~isempty(handle.draw.handle)
            close(handle.draw.handle);
        end;

        pathsave = prm.pathsettings;
        msg = ['Saving settings for VIEWSEGM in ' pathsave];
        disp(msg);
        msg = ['To restore default settings delete ' pathsave];
        disp(msg);
        % visualization channels
        handle.vis = prm.vis;
        save(pathsave,'handle');
        break
    end;
    
    
end;

%-----------------------------------------------


function [draw] = showimagedraw(im,draw,plane)

v = get(draw.editfield.handle,'String');
try
    v = str2double(v);
    v = round(v);    
%     scrsz = get(0,'ScreenSize');
    pos = [1000 100 1200 1200];
    if isempty(draw.handle)
        draw.handle = figure('Position',pos);
    end;
    implane = im(:,:,plane,v);
    lim = [min(implane(:)) max(implane(:))];
    draw.handle = figure(draw.handle);imshow(implane,lim);colormap(gray);axis image;axis off;
catch
    warning('Channel must be valid');
    set(draw.editfield.handle,'');
end;

%--------------------------------------------------------------------

% function key_releaseFcn
% % Press and release various key combinations in the figure.
% % Values returned by the event structure are displayed
% % in the command window.
% %
% 
% global ch
% h = figure('KeyReleaseFcn',@cb);
% annotation('textbox',[.25 .4  .5 .2],...
%   'FitHeightToText','on',...
%   'String',{'Press and release various key combination';...
%   'to see the result in the command window'});
% 
% %--------------------------------------------------

function [] = cb(src,evnt)

global ch;
% if ~isempty(evnt.Modifier)
%    for ii = 1:length(evnt.Modifier)
%       out = sprintf('Character: %c\nModifier: %s\nKey: %s\n',evnt.Character,evnt.Modifier{ii},evnt.Key);
%       disp(out)
%    end
% else
%    out = sprintf('Character: %c\nModifier: %s\nKey: %s\n',evnt.Character,'No modifier key',evnt.Key);
%    disp(out)
% end
ch = evnt.Key;
% uiresume;

%-------------------------------------------
      
function [] = saveline(pathstack,linedata,prm)

h = prm.h;
[a b c] = fileparts(pathstack);
pathsave = [a b '-line.mat'];
msg = ['Saving line data ' pathsave];
disp(msg);
save(pathsave,'linedata','h');

%------------------------------------------------------------------
        
function [] = savestack(pathstack,info,cellbw)

msg = ['Saving modified segmentation in ' pathstack];
disp(msg);
save(pathstack,'cellbw','info','-append');

%--------------------------------------------------------------------

function [handle] = showimagecells(im,linedata,cellbw,wat,minima,prm,handle)

imscale = handle.imscale;

% only use these scales if checkbox is on
for i = 1 : prm.nch
    flag = get(handle.imscale.checkbox(i).handle,'Value');
    % for saving
    handle.imscale.checkbox(i).value = flag;
    
    if flag == 1
        % get image scale values
        v = get(handle.imscale.checkbox(i).editfield(1).handle,'String'); 
        vstr = v;
        if ~isempty(v)
            try
                v = str2double(v);
                prm.imscale(i).sc{1} = v;
            catch
                vstr = '';
                prm.imscale(i).sc{1} = vstr;
                set(handle.imscale.checkbox(i).editfield(1).handle,vstr);
            end;
        end;    
        % for saving
        handle.imscale.checkbox(i).editfield(1).value = vstr;
    
        % get image scale values
        v = get(handle.imscale.checkbox(i).editfield(2).handle,'String');        
        vstr = v;
        if ~isempty(v)
            try
                v = str2double(v);
                prm.imscale(i).sc{2} = v;
            catch
                vstr = '';
                prm.imscale(i).sc{2} = vstr;
                set(handle.imscale.checkbox(i).editfield(2).handle,vstr);
            end;
        end;
        % for saving
        handle.imscale.checkbox(i).editfield(2).value = vstr;
    end;
end;

sc = cell(prm.nch,1);
plane = min(prm.plane,prm.dim(3));
plane = max(plane,1);

% these planes
for j = 1 : prm.nch
    implane{j} = im(:,:,plane,prm.vis(j).ch);
end;

l = 0;
w = 1/prm.nch-0.10/prm.nch;
b = 0.5;
h = 0.45;
dl = 1/prm.nch;

% raw images with scaling
figure(handle.fig.handle);
for j = 1 : prm.nch        
        
    pos{j} = [l b w h];
    sc{j} = [min(implane{j}(:)) max(implane{j}(:))];
    for i = 1 : 2
         if ~isempty(prm.imscale(j).sc{i})
            sc{j}(i) = prm.imscale(j).sc{i};
        end;
    end;
        
    if sc{j}(2) <= sc{j}(1)
        warning('Cannot have high limit smaller than low, setting high');
        v = sc{j}(1) + 1;
        v = round(v);
        sc{j}(2) = v;
        v = int2str(v);
        set(imscale.checkbox(j).editfield(2).handle,'String',v);
    end;
    
    figure(handle.fig.handle);
    subpl.handle(1,j) = subplot('Position',pos{j},'Parent',handle.fig.handle);imshow(implane{j},sc{j});colormap(gray);drawnow;axis off;axis image;
    if isequal(prm.vis(j).ch,'imsegm')
        axis image;title(['Channel ' prm.vis(j).ch])
    else
        axis image;title(['Channel ' int2str(prm.vis(j).ch)])
    end;
    l = l + dl;
end;

cellbwplane = double(cellbw(:,:,plane));
perimplane = bwperim(cellbwplane);
minimaplane = minima(:,:,plane);

dim = size(cellbwplane);
Dim = [dim 3];
figure(handle.fig.handle);
for j = 1 : prm.nch-1    
    
    % overlay image
    implane1 = implane{j};
    implane1(implane1 < sc{j}(1)) = sc{j}(1);
    implane1(implane1 > sc{j}(2)) = sc{j}(2);
    implane1 = scale(implane1);
    
    implane2 = implane{j};
    implane2(implane2 < sc{j}(1)) = sc{j}(1);
    implane2(implane2 > sc{j}(2)) = sc{j}(2);    
    implane2 = scale(implane2);
    
    implane3 = implane{j};
    implane3(implane3 < sc{j}(1)) = sc{j}(1);
    implane3(implane3 > sc{j}(2)) = sc{j}(2);    
    implane3 = scale(implane3);
    
    % cell boundaries as red
    implane1(perimplane == 1) = 1;

    
    % markers as blue
    implane3(minimaplane == 1) = 1;
    

    rgboverlayplane = zeros(Dim);
    rgboverlayplane(:,:,1) = implane1;
    rgboverlayplane(:,:,2) = implane2;
    rgboverlayplane(:,:,3) = implane3;
    
    rgboverlayplane = uint8(round(rgboverlayplane*255));    
    p = pos{j};
    p(2) = 0.01;
    subpl.handle(2,j) = subplot('Position',p,'Parent',handle.fig.handle);imagesc(rgboverlayplane);colormap(gray);drawnow;axis off;axis image;
    axis image;title('Overlay image, markers (blue), segmentation (red)');    
    
end;

% cellbw image
watplane = wat(:,:,plane);
cellbwplane(watplane == 0) = 0.5;
p = pos{end};
p(2) = 0.01;
figure(handle.fig.handle);
subpl.handle(2,end) = subplot('Position',p,'Parent',handle.fig.handle);imagesc(cellbwplane);colormap(gray);drawnow;axis off;axis image;
axis image;title('Found cells');

% make lines
for i = 1 : size(subpl.handle,1)
    for j = 1 : size(subpl.handle,2)
        for k = 1 : numel(linedata.x)
            x = linedata.x{k};
            y = linedata.y{k};
            z = linedata.z{k};
            if plane == z(1)
                subplot(subpl.handle(i,j));
                line(x,y,'Color','g');                      
            end;
        end;
    end;
end;

% image to draw on
v = get(handle.draw.checkbox.handle,'Value');
v = double(v);
% checkbox is on
if v == 1
    handle.draw = showimagedraw(im,handle.draw,plane);
% checkbox is off
elseif v == 0
    if ~isempty(handle.draw.handle)
        close(handle.draw.handle);
    end;
    handle.draw.handle = [];
end

if ~isempty(linedata.x)
    for k = 1 : numel(linedata.x)
        x = linedata.x{k};
        y = linedata.y{k};
        z = linedata.z{k};
        if plane == z(1)
            figure(handle.draw.handle);
            line(x,y,'Color','g');                      
        end;
    end;    
end;


%----------------------------------------------------------------

function [] = showimagetnt()

global plane vischglobal im cellbw tntcand wat overlaynum fignum external



% these planes
% imhere = 0.5*scale(im(:,:,plane,vischglobal));
imhere = (im(:,:,plane,vischglobal));

cellbwhere = double(cellbw(:,:,plane));
perimhere = bwperim(cellbwhere);


% TNT candidates
tntcandhere = double(tntcand(:,:,plane));

% to show the watershed that were not classified as cells
wathere = wat(:,:,plane);
cellbwhere(wathere == 0) = 0.5;

% the TNTs in this plane
tntcandhere = logical(tntcand(:,:,plane));

cellbwhere2 = double(cellbw(:,:,plane));
rgboverlay(:,:,1) = scale(cellbwhere2 + tntcandhere);
rgboverlay(:,:,1) = scale(cellbwhere2 + tntcandhere);
rgboverlay(:,:,2) = scale(cellbwhere2 + tntcandhere);
rgboverlay(:,:,2) = scale(cellbwhere2 + tntcandhere);
rgboverlay(:,:,3) = scale(cellbwhere2 + tntcandhere);
rgboverlay(:,:,3) = scale(cellbwhere2 + tntcandhere);

for i = 1  : length(overlaynum)
   channel = overlaynum(i);
    
   % if shall be empty
   if channel == 0
       continue;
   end;
   
   % overlay
   rgboverlay(:,:,i) = scale(rgboverlay(:,:,i) + im(:,:,plane,overlaynum(i))); 
end

% figure(10);imagesc(rgboverlay)

% rgbOverlayhere = imhere;
% rgbOverlayhere(perimhere) = 1;
% rgbOverlayhere(tntcandhere) = 1;
% rgbOverlayhere(:,:,2) = imhere;
% rgbOverlayhere(:,:,3) = imhere;

d{1} = imhere;
d{2} = cellbwhere;
d{2} = cellbwhere;
d{3} = tntcandhere;
d{4} = rgboverlay;

figure(fignum);
for i = 1 : 4
    subplot(2,2,i);imagesc(d{i});colormap(gray);drawnow;axis off
end;

% external?
for i = 1 : length(external)
    channel = external(i);
    figure(fignum + i);imagesc(d{channel});colormap(gray);drawnow;axis off;
end;


