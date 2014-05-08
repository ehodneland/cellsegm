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
% global ch;

start = varargin{1};
stop = varargin{2};
prm.vis(1).ch = varargin{3};
prm.nch = 0;
prm.option = 'segm';
if nargin < 3
    error('Tou must specify an image channel');
end;
c = 3;
ch = 0;
cp = 0;
prm.plugin.n = 0;
while 1
    v = varargin{c};
    
    if isnumeric(v)        
        ch = ch + 1;        
        prm.vis(ch).ch = varargin{c};
        prm.nch = prm.nch + 1;                
        c = c + 1;
        
    elseif isequal(v,'plugin')
        cp = cp + 1;        
        prm.plugin.name{cp} = varargin{c+1};
        prm.plugin.n = prm.plugin.n + 1;
        c = c + 2;
    elseif ischar(v)        
        prm.option = v;        
        c = c + 1;
    end;
    if c > nargin
        break;
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
handle.fig.handle = figure('Position',handle.fig.pos);


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
    'Position',[left b w h],'CallBack',{@callbackup,handle});  
b = b - h - 5;
handle.down.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Down', ...
    'Position',[left b w h],'CallBack',{@callbackdown,handle});   
b = b - h - 5;
handle.next.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Next', ...
    'Position',[left b w h],'CallBack',{@callbacknext,handle});
b = b - h - 5;
handle.previous.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Previous', ...
    'Position',[left b w h],'CallBack',{@callbackprevious,handle});
b = b - h - 5;
handle.frame.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Frame number', ...
    'Position',[left b w h],'CallBack',{@callbackframenum,handle});
b = b - h - 5;
handle.classification.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Classification', ...
    'Position',[left b w h],'CallBack',{@callbackclassification,handle});
% Delete/add cells permanently!
b = b - h - 5;
handle.cellstatus.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Cell status', ...
    'Position',[left b w h],'CallBack',{@callbackcellstatus,handle});


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
         'Position',handle.imscale.checkbox(j).pos,'CallBack',{@callbackcheckboxscale,handle});    
    dw = dw + handle.checkbox.w;
    c = c + 1;
    handle.imscale.checkbox(j).editfield(1).pos = [left+dw b handle.editfield.w h];    
    handle.imscale.checkbox(j).editfield(1).handle = uicontrol('Parent',handle.control.handle,'Style','edit', ...
        'Position',handle.imscale.checkbox(j).editfield(1).pos,'CallBack',{@callbackcheckboxscale,handle});
    
    dw = dw + handle.checkbox.w;
    c = c + 1;
    handle.imscale.checkbox(j).editfield(2).pos = [left+dw b handle.editfield.w h];
    handle.imscale.checkbox(j).editfield(2).handle = uicontrol('Parent',handle.control.handle,'Style','edit', ...
        'Position',handle.imscale.checkbox(j).editfield(2).pos,'CallBack',{@callbackcheckboxscale,handle});
        
end;

% Remove small cells
% NB not permanently, only visually
b = b - h - 5;
dw = 0;
handle.celldim.button.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Vol(mcm^3)', ...
    'Position',[left b w h],'CallBack',{@callbackvolume,handle});  
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
    'Position',[left b w h],'CallBack',{@callbackboundary,handle});
dw = dw + w + handle.checkbox.w;
handle.boundary.editfield.handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','0.05', ...
    'Position',[left+dw b handle.editfield.w h]);    

% % draw line
% b = b - h - 5;
% dw = 0;
% handle.draw.button.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Draw', ...
%     'Position',[left b w h],'CallBack',{@callbackdrawim,handle});
% dw = dw + w;
% % handle.draw.checkbox.pos = [left+dw b handle.checkbox.w h];
% % handle.draw.checkbox.handle = uicontrol('Parent',handle.control.handle,'Style','Checkbox', ...
% %          'Position',handle.draw.checkbox.pos,'CallBack',{@callbackdrawcheck,handle});
% % set(handle.draw.checkbox.handle,'Value',0);
% dw = dw + handle.checkbox.w;
% handle.draw.editfield.pos = [left+dw b handle.editfield.w h];
% handle.draw.editfield.handle = uicontrol('Parent',handle.control.handle,'Style','edit','String','1', ...
%     'Position',handle.draw.editfield.pos);    

% Quit
b = b - h - 5;
handle.quit.handle = uicontrol('Parent',handle.control.handle,'Style','PushButton','String','Quit', ...
    'Position',[left b w h],'CallBack',{@callbackquit,handle});  
b = b - h - 5;
handle.control.title.pos = [left b 4*w h];
handle.control.title.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','Message field', ...
     'Position',handle.control.title.pos);   

% show position
left = handle.button.left;

% Bottom line show coordinates of mousepointer
b = 1;
handle.cellstatus.coordx.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','cursorx', ...
    'Position',[left b w h]);
left = left + w;
handle.cellstatus.coordy.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','cursory', ...
    'Position',[left b w h]);
left = left + w;
handle.cellstatus.watval.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','watval', ...
    'Position',[left b w h]);

% Bottom line plane
b = b + h;
left = handle.button.left;
handle.up.plane.handle = uicontrol('Parent',handle.control.handle,'Style','Text','String','plane', ...
    'Position',[left b w h]);

% to save settings
if ismac
    prm.pathsettings = fullfile('Users',prm.username,'.cellsegm');
elseif isunix
    prm.pathsettings = fullfile('home',prm.username,'.cellsegm');
elseif ispc
    prm.pathsettings = fullfile('Documents and Settings',prm.username,'.cellsegm');
end;
prm.pathsettings = ['/' prm.pathsettings];

if ~isdir(prm.pathsettings)
    mkdir(prm.pathsettings);
end;
prm.pathsettings = fullfile(prm.pathsettings,'settings.mat');

% load last used settings
if exist(prm.pathsettings,'file')
    msg = ['Loading previous settings'];
    disp(msg);
    D = load(prm.pathsettings);    
    plane = D.plane;
    % list of handles
    for k = 1 : prm.nch
        for j = 1 : numel(D.handle.vis)
            % we have the same channel as in the settings file?
            if D.handle.vis(j).ch == prm.vis(k).ch
                try
                    % set the two edit fields
                    v = D.handle.imscale.checkbox(k).editfield(1).value;
                    set(handle.imscale.checkbox(i).editfield(1).handle,'String',v);
                catch
                    
                end;
                try
                    v = D.handle.imscale.checkbox(k).editfield(2).value;
                    set(handle.imscale.checkbox(k).editfield(2).handle,'String',v);                    
                catch
                    
                end;
                try
                    % set the checkbox
                    v = D.handle.imscale.checkbox(k).value;
                    set(handle.imscale.checkbox(k).handle,'Value',v);                    
                catch
                    
                end;
            end;
        end;
    end;
end;

% draw image
pos = [0 0 1200 1000];
handle.draw.handle = figure('Position',pos);
handle.draw.pos = pos;
set(handle.draw.handle,'Visible','off');
set(handle.draw.handle,'WindowButtonMotionFcn',{@drawcoordupdate,handle});
set(handle.draw.handle,'pointer','Crosshair');
set(handle.draw.handle,'WindowButtonDownFcn',{@drawcoordselect,handle});
%uiwait(h);
%set(h,'pointer','Arrow');
%uirestore(uistate);

while 1
    
    choice = '';
    
    pathbase = [ 'stack' int2str(i)];
    if isequal(prm.option,'segm')
        pathstack = [ pathbase '-segm.mat'];
    elseif isequal(prm.option,'raw')
        pathstack = [ pathbase '.mat'];
    end;
    prm.pathstack = pathstack;
    prm.pathbase = pathbase;

    if i > stop
        close(handle.fig.handle);
        close(handle.control.handle);
        return;
    end;
      
    % load new image if necessary
    if ne(i,iold)
        
        try
            D = load(pathstack);
            disp(sprintf('Loading %s',pathstack));
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
            
            % number of lines
            cline = 0;
            linedata.x = [];
            linedata.y = [];
            
            % plugin
            prm.plugin.plugindonorcell.loaded = 0;
            
            % cellstatus
            handle.cellstatus.value = [];
            
            prm.voxelvol = prod(info.prm.h);
            clear D;
            msg = ['Image dimension: ', num2str(dim)];
            disp(msg);

                   
        catch
            msg = ['Could not load ' pathstack];
            disp(msg);
            if isequal(choice,'previous')
                i = i - 1;
            else isequal(choice,'next')
                i = i + 1;
            end;
            continue;
        end;
        
        % set figure size
        set(handle.draw.handle,'Position',[500 0 dim(1) dim(2)]);
        
    end;
    iold = i;        
            
   % figure size in relation to screen
    scrsz = get(0,'ScreenSize');
    scrsz = scrsz(3:4);
    scrsz = scrsz-150;
    prm.maxw = scrsz(1);
    prm.maxh = scrsz(2);
    ratio = dim(2)/dim(1);
    if ratio > 1
        w = prm.maxw;
        h = w/ratio;
    elseif ratio <= 1
        h = prm.maxh;
        w = h*ratio;
    end;        
    handle.fig.w = w;
    handle.fig.h = h;

    % show images
    prm.plane = plane;
    prm.dim = dim;   
    [handle,prm] = showimagecells(im,linedata,cellbwvis,wat,minima,prm,handle);
    set(handle.fig.handle,'WindowKeyPressFcn',@callbackreading);
%     set(handle.fig.handle,'WindowKeyPressFcn',@callbackreading);
    
    % wait for user
%     figure(handle.fig.handle);
    uiwait(handle.fig.handle);
        
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
        callbckclassification(handle,wat,info,plane);
               
    elseif isequal(choice,'cellstatus')
               
        % change the cell status
        [info,prm,cellbw,cellbwvis] = callbckcellstatus(wat,cellbw,cellbwvis,info,prm,handle,plane,linedata,im,minima);
       
    elseif isequal(choice,'volume')
        
        % remove cells that are either below a threshold or above a
        % threshold
        [cellbwvis,boundarybw,volbw,prm] = callbckvolume(cellbwvis,boundarybw,prm,cellbw,handle);        
        
    elseif isequal(choice,'boundary')
        
        % remove cells that are attached to the boundary
        [cellbwvis,boundarybw,prm] = callbckboundary(cellbwvis,prm,cellbw,volbw,handle);
                        
%     elseif isequal(choice,'drawim')      
%         
%         % draw a line in the image
%         handle = callbckdrawim(handle,plane,im);
%         
%     elseif isequal(choice,'drawselect')
%                 
%         % draw the image to draw on
%         [handle,linedata,cline] = callbckdrawselect(handle,linedata,plane,cline);
        
    elseif isequal(choice,'quit')
        if prm.cellmod == 1
            savestack(pathstack,info,cellbw);        
        end;
        if ~isempty(linedata.x)
            saveline(pathstack,linedata,prm);
        end;
        
        finish = 1;
    end;

        
    set(handle.up.plane.handle,'String',['plane ' int2str(plane)]);
  
    if finish == 1
        close(handle.fig.handle);
        close(handle.control.handle);
        if ~isempty(handle.draw.handle)
            close(handle.draw.handle);
        end;

        pathsave = prm.pathsettings;
        msg = ['Saving settings for VIEWSEGM in ' pathsave];
        disp(msg);
        msg = ['To restore default settings, delete ' pathsave];
        disp(msg);
        % visualization channels
        handle.vis = prm.vis;
        save(pathsave,'handle','plane');
        break
    end;
    
    
end;

%-------------------------------------------
function []  = callbackcheckboxscale(src,evnt,handle)
global choice;
choice = '';
uiresume(handle.fig.handle);

%-------------------------------------------
function [] = callbackup(src,evnt,handle)
global choice;
choice = 'up';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbackdown(src,evnt,handle)
global choice;
choice = 'down';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbackprevious(src,evnt,handle)
global choice;
choice = 'previous';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbacknext(src,evnt,handle)
global choice;
choice = 'next';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbackframenum(src,evnt,handle)
global choice;
choice = 'framenum';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbackclassification(src,evnt,handle)
global choice;
choice = 'classification';
uiresume(handle.fig.handle);
%-------------------------------------------
function [] = callbackcellstatus(src,evnt,handle)
global choice;
choice = 'cellstatus';
uiresume(handle.fig.handle);
%--------------------------------------------
function []  = callbackvolume(src,evnt,handle)
global choice;
choice = 'volume';
uiresume(handle.fig.handle);
%--------------------------------------------
% function []  = callbackdrawcheck(src,evnt,handle)
% global choice;
% choice = 'drawcheck';
% uiresume(handle.fig.handle);
%--------------------------------------------
function []  = callbackdrawim(src,evnt,handle)
global choice;
choice = 'drawim';
uiresume(handle.fig.handle);
%--------------------------------------------
function []  = callbackboundary(src,evnt,handle)
global choice;
choice = 'boundary';
uiresume(handle.fig.handle);
%--------------------------------------------
function []  = callbackquit(src,evnt,handle)
global choice;
choice = 'quit';
uiresume(handle.fig.handle);
%--------------------------------------------
function [] = callbackreading(src,evnt)
global choice;
switch(evnt.Key)
    case 'uparrow'
        choice = 'up';
        uiresume(src);
    case 'downarrow'
        choice = 'down';
        uiresume(src);
    otherwise
        choice = '';
end
%---------------------------------------------
function [] = cellstatuscoordupdate(src,evnt,handle,wat,plane,cellbw)
     
handlesubpl = handle.subpl(1,1).handle;
v = get(handlesubpl,'CurrentPoint');
w2 = round(v(1,1));
w1 = round(v(1,2));
v1 = num2str(w1);
v2 = num2str(w2);
set(handle.cellstatus.coordx.handle,'String',v1);
set(handle.cellstatus.coordy.handle,'String',v2);
set(handle.cellstatus.coordx.handle,'UserData',w1);
set(handle.cellstatus.coordy.handle,'UserData',w2);
try
    value = wat(w1,w2,plane);
catch
    value = []; 
end;
% valueold = get(handle.cellstatus.watval.handle,'UserData');
valueold = get(handle.cellstatus.watval.handle,'String');
valueold = str2double(valueold);
    
% when to draw a gray cell?
draw = 0;
if ne(valueold,value)
    if ~isequal(value,0)
        draw = 1;
    end;
end;
if isnumeric(valueold) && isempty(value)
    draw = 1;
end;
if draw
    % cellbw image
    cellbwplane = cellbw(:,:,plane);           
    watplane = wat(:,:,plane);
    cellbwplane(watplane == 0) = 0.5;
    if ~isempty(value)
        cellbwplane(watplane == value) = 0.5;
    end;
    figure(handle.fig.handle);
    subplot('Position',handle.subpl(2,end).pos,'Parent',handle.fig.handle);
    imagesc(cellbwplane);axis off;colormap(gray);drawnow;axis image;
    axis image;title('Found cells');
end;

set(handle.cellstatus.watval.handle,'String',num2str(value));
% get(handle.cellstatus.watval.handle,'String')
% set(handle.cellstatus.watval.handle,'UserData',value);

 %-------------------------------------------
 
 function [] = cellstatuscoordselect(src,evnt)     
 uiresume(src);
   
 
%-----------------------------------------------

function [info,prm,cellbw,cellbwvis] = callbckcellstatus(wat,cellbw,cellbwvis,info,prm,handle,plane,linedata,im,minima)


h = handle.fig.handle;                    
uistate = uisuspend(h);

set(h,'WindowButtonMotionFcn',{@cellstatuscoordupdate,handle,wat,plane,cellbw});
set(h,'pointer','Crosshair');
set(h,'WindowButtonDownFcn',{@cellstatuscoordselect});
uiwait(h);
set(h,'pointer','Arrow');
uirestore(uistate);

x = round(get(handle.cellstatus.coordx.handle,'UserData'));
y = round(get(handle.cellstatus.coordy.handle,'UserData'));

if x < 1 || y < 1 || x > prm.dim(1) || y > prm.dim(2)
    msg = ['You hit outside the image'];
    disp(msg);
    return;
end;

val = wat(x,y,plane);

% this is a boundary of a cell, try again
if val == 0
    msg = ['You hit outside the boundary of a cell'];
    disp(msg);
    return;
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

%-------------------------------------------

% function [handle] = callbckdrawcheck(handle,im,plane)
% 
% v = get(handle.draw.checkbox.handle,'Value');
% % v
% % pause
% % v = get(handle.draw.checkbox.handle,'Value');
% % v
% if v == 1
%     showimagedraw(im,handle.draw,plane);
% else
%     set(handle.draw.handle,'Visible','off');
% end;
% if v == 0
%     set(handle.draw.checkbox.handle,'Value',1);    
% elseif v == 1
%     set(handle.draw.checkbox.handle,'Value',0);
% end;

%--------------------------------------

function [] = drawcoordupdate(src,evnt,handle)

figure(handle.draw.handle);
v = get(gca,'CurrentPoint');
w2 = round(v(1,1));
w1 = round(v(1,2));
v1 = num2str(w1);
v2 = num2str(w2);

set(handle.cellstatus.coordx.handle,'String',v1);
set(handle.cellstatus.coordy.handle,'String',v2);
set(handle.cellstatus.coordx.handle,'UserData',w1);
set(handle.cellstatus.coordy.handle,'UserData',w2);

%-------------------------------------------

function [handle] = callbckdrawim(handle,plane,im)

v = get(handle.draw.handle,'Visible');
if isequal(v,'off');
    set(handle.draw.handle,'Visible','on')
else isequal(v,'on');
    % do no more
    set(handle.draw.handle,'Visible','off');
    return;
end;

% v = get(handle.draw.handle,'Position')
% set(handle.draw.handle,'Position',v);

% v = get(handle.draw.editfield.handle,'String');
% v = str2double(v);
% implane = im(:,:,plane,v);


% handle.draw = showimagedraw(im,handle.draw,plane);
showimagedraw(im,handle.draw,plane);
 
return;

% minim = min(implane(:));
% maxim = max(implane(:));
% lim = [minim, maxim];
% figure(handle.draw.handle);
% imshow(implane,lim);colormap(gray);axis image;axis off;

% pos = get(handle.draw.handle,'Position')

% plot older lines
for k = 1 : numel(linedata.x)
    x = linedata.x{k};
    y = linedata.y{k};
    z = linedata.z{k};
    if plane == z(1)
        figure(handle.draw.handle);
%         set(handle.draw.handle,'Position',pos);
        line(y,x,'Color','g');                      
    end;
end;

uistate = uisuspend(handle.draw.handle);
% set(handle.draw.handle,'Units','pixels');
c = 0;
while 1
    
    set(handle.draw.handle,'WindowButtonMotionFcn',{@cellstatuscoordupdatedraw,handle});
    set(handle.draw.handle,'pointer','Crosshair');
    set(handle.draw.handle,'WindowButtonDownFcn',{@cellstatuscoordselect});
%     set(handle.draw.handle,'pointer','Arrow');    
    uiwait(handle.draw.handle);    
    
    v = round(get(handle.draw.button.handle,'UserData'));
    
    
    if v(1) < 1 || v(2) < 1 || v(1) > prm.dim(1) || v(2) > prm.dim(2)
        msg = ['You hit outside the image'];
        disp(msg);
        continue;
    end;
    c = c + 1;
    
    x(c,1) = v(1);y(c,1) = v(2);    
    
    if c == 2
        break;
    end;

end;

cline = cline + 1;
linedata.x{cline} = x;
linedata.y{cline} = y;            
linedata.z{cline} = [plane; plane];            
uirestore(uistate);

% set(handle.draw.handle,'Visible','off');

%---------------------------------------------

function [] = drawcoordselect(src,evnt,handle)     

global choice;
choice = 'drawselect';
uiresume(handle.fig.handle);
 
 %-------------------------------------------
 
function [handle,linedata,cline] = callbckdrawselect(handle,linedata,plane,cline)
     
figure(handle.draw.handle);
v = get(gca,'CurrentPoint');
v = v(1,1:2);
v = fliplr(v);
v = round(v);
vold = get(handle.draw.button.handle,'UserData');
if ~isempty(vold)
    cline = cline + 1;
    linedata.x{cline}(1,1) = vold(1);
    linedata.x{cline}(2,1) = v(1);
    linedata.y{cline}(1,1) = vold(2);
    linedata.y{cline}(2,1) = v(2);
    linedata.z{cline}(1,1) = plane;
    linedata.z{cline}(2,1) = plane;
    set(handle.draw.button.handle,'UserData','');    
else
    set(handle.draw.button.handle,'UserData',v);
end;

% get(handle.draw.button.handle,'UserData')

% get(handle.draw.handle,'Position')
 
%---------------------------------------------

function [cellbwvis,boundarybw,prm] = callbckboundary(cellbwvis,prm,cellbw,volbw,handle)
    
dim = prm.dim;
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

%--------------------------------------------

function [cellbwvis,boundarybw,volbw,prm] = callbckvolume(cellbwvis,boundarybw,prm,cellbw,handle)

dim = prm.dim;

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

%----------------------------------------------

function [] = callbckclassification(handle,wat,info,plane)

% user defined input
figure(handle.fig.handle);
[y, x] = ginput(1);
x = round(x);
y = round(y);
val = wat(x,y,plane);
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

%---------------------------------------------


function [] = showimagedraw(im,draw,plane)

v = get(draw.editfield.handle,'String');
% pos = get(draw.handle,'Position');

% try
    v = str2double(v);
    v = round(v);    
    implane = im(:,:,plane,v);
    lim = [min(implane(:)) max(implane(:))];
    figure(draw.handle);imshow(implane,lim);colormap(gray);axis image;axis off;
%     set(draw.handle,'Position',pos);
% catch
%     warning('Channel must be valid');
%     set(draw.editfield.handle,'String','1');
% end;

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
      
function [] = saveline(pathstack,linedata,prm)

h = prm.h;
[a,b,c] = fileparts(pathstack);
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

function [handle,prm] = showimagecells(im,linedata,cellbw,wat,minima,prm,handle)

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



%
% Show raw images in the upper row
%

% these planes
implane = cell(prm.nch,1);
for j = 1 : prm.nch
    implane{j} = im(:,:,plane,prm.vis(j).ch);
end;

l = 0;
w = 1/prm.nch-0.10/prm.nch;
b = 0.5;
h = 0.45;
dl = 1/prm.nch;

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
    handle.subpl(1,j).data = implane{j};    
    handle.subpl(2,j).sc = sc{j};
    % same scaling for lower row
    handle.subpl(1,j).sc = sc{j};
    handle.subpl(1,j).handle = subplot('Position',pos{j},'Parent',handle.fig.handle);imshow(handle.subpl(1,j).data,handle.subpl(1,j).sc);
    colormap(gray);drawnow;axis off;axis image;
    handle.subpl(1,j).pos = pos{j};
    if isequal(prm.vis(j).ch,'imsegm')
        axis image;title(['Channel ' prm.vis(j).ch])
    else
        axis image;title(['Channel ' int2str(prm.vis(j).ch)])
    end;
    l = l + dl;
end;

%
% Show the nuclei and lines of segmentation in the lower row
%

cellbwplane = double(cellbw(:,:,plane));
perimplane = bwperim(cellbwplane);
minimaplane = minima(:,:,plane);
dim = size(cellbwplane);
Dim = [dim 3];
figure(handle.fig.handle);
for j = 1 : prm.nch-1    
    
    % overlay image
    implane1 = implane{j};
    sc = handle.subpl(1,j).sc;
    implane1(implane1 < sc(1)) = sc(1);
    implane1(implane1 > sc(2)) = sc(2);
    implane1 = scale(implane1);
    
    implane2 = implane{j};
    implane2(implane2 < sc(1)) = sc(1);
    implane2(implane2 > sc(2)) = sc(2);
    implane2 = scale(implane2);
    
    implane3 = implane{j};
    implane3(implane3 < sc(1)) = sc(1);
    implane3(implane3 > sc(2)) = sc(2);
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
    handle.subpl(2,j).data = rgboverlayplane;
    sc = [0 255];
    handle.subpl(2,j).sc = sc;
    handle.subpl(2,j).pos = p;
    handle.subpl(2,j).handle = subplot('Position',p,'Parent',handle.fig.handle);imshow(handle.subpl(2,j).data,handle.subpl(2,j).sc);
    colormap(gray);drawnow;axis off;axis image;axis image;title('Overlay image, markers (blue), segmentation (red)');    
%     subpl.pos{2,j} = p;
end;

% cellbw image
watplane = wat(:,:,plane);
cellbwplane(watplane == 0) = 0.5;

%
% Show the cell segmentation
%
p = pos{end};
p(2) = 0.01;
% value = get(handle.cellstatus.watval.handle,'String');
% value = str2double(value);
% if ~isempty(value)
%     cellbwplane(watplane == value) = 0.5;
% end;
figure(handle.fig.handle);
% make RGB of this to enable colored overlays
for i = 2 : 3
    cellbwplane(:,:,i) = cellbwplane(:,:,1);
end;
handle.subpl(2,end).data = cellbwplane;
handle.subpl(2,end).sc = [0 1];
handle.subpl(2,end).pos = p;
handle.subpl(2,end).handle = subplot('Position',handle.subpl(2,end).pos,'Parent',handle.fig.handle);imshow(handle.subpl(2,end).data,handle.subpl(2,end).sc);
colormap(gray);drawnow;axis off;axis image;
axis image;title('Found cells');



%
% Show lines
%

% make lines
for i = 1 : size(handle.subpl,1)
    for j = 1 : size(handle.subpl,2)
        for k = 1 : numel(linedata.x)
            x = linedata.x{k};
            y = linedata.y{k};
            z = linedata.z{k};
            if plane == z(1)
                subplot(handle.subpl(i,j).handle);
                line(y,x,'Color','g');                      
            end;
        end;
    end;
end;

%
% Show image to draw on
%
v = get(handle.draw.handle,'Visible');
if isequal(v,'on')
    % visible
%     handle.draw = showimagedraw(im,handle.draw,plane);    
    showimagedraw(im,handle.draw,plane);    
    for k = 1 : numel(linedata.x)
        x = linedata.x{k};
        y = linedata.y{k};
        z = linedata.z{k};
        if plane == z(1)
            figure(handle.draw.handle);
            line(y,x,'Color','g');                      
        end;
    end;
end;

% Load any plugins
for i = 1 : prm.plugin.n
    if isequal(prm.plugin.name{i},'plugindonorcell')
        prm = plugindonorcell(handle,prm);
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


