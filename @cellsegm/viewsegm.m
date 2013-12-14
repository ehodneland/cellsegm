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

global choice 

start = varargin{1};
stop = varargin{2};
prm.visch1 = varargin{3};
prm.visch2 = varargin{4};
prm.option = 'segm';
if nargin == 5
    prm.option = varargin{5};
end;
% vischglobal = visch;


% fig.pos = [400 50 500 500];
% fig.handle = figure('Position',fig.pos);
% control.pos = [100 50 300 500];
% control.handle = figure('Position',control.pos);

% image number
i = start;
iold = inf;
finish = 0;
plane = 1;


% % change size if necessary
% p = get(fig.handle,'OuterPosition');

fig.w = 800;
fig.h = 700;
fig.b = 100;
fig.l = 400;
fig.pos = [fig.l fig.b fig.w fig.h];
fig.handle = figure('Position',fig.pos);
% set(fig.handle,'Position',fig.pos);

control.w = 400;
control.h = fig.h;    
control.l = fig.l;
control.h = fig.h;
control.pos = [fig.pos(1)-control.w fig.b control.w control.h];
control.handle = figure('Position',control.pos);
% set(control.handle,'Position',control.pos);


% figure(fig.handle);
set(fig.handle,'Toolbar','figure');
set(control.handle,'Toolbar','figure');

% % margins and settings
% margin.left = 50;
% margin.right = 20;
% margin.top = 10;
% margin.bottom = 5;
button.h = 30;
button.w = 100;
button.left = 10;

left = button.left;
b = control.h - 2*button.h;
w = button.w;
h = button.h;
handle.h1 = uicontrol('Parent',control.handle,'Style','PushButton','String','Up', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 1;uiresume(gcbf);');    
b = b - h - 5;
handle.h2 = uicontrol('Parent',control.handle,'Style','PushButton','String','Down', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 2;uiresume(gcbf);');    
b = b - h - 5;
handle.h3 = uicontrol('Parent',control.handle,'Style','PushButton','String','Next', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 3;uiresume(gcbf);');    
b = b - h - 5;
handle.h4 = uicontrol('Parent',control.handle,'Style','PushButton','String','Previous', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 4;uiresume(gcbf);');    
b = b - h - 5;
handle.h5 = uicontrol('Parent',control.handle,'Style','PushButton','String','Frame number', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 5;uiresume(gcbf);');    
b = b - h - 5;
handle.h7 = uicontrol('Parent',control.handle,'Style','PushButton','String','Classification', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 6;uiresume(gcbf);');    
b = b - h - 5;
handle.h6 = uicontrol('Parent',control.handle,'Style','PushButton','String','Quit', ...
    'Position',[left b w h],'CallBack', ...
    'global choice;choice = 7;uiresume(gcbf);');        
b = b - h - 5;
control.title.pos = [left b 4*w h];
control.title.handle = uicontrol('Parent',control.handle,'Style','Text','String','Message field', ...
     'Position',control.title.pos);   

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
                prmin = D.info.classifycells;
            elseif isequal(prm.option,'raw')
                cellbw = zeros(dim3);
                wat = zeros(dim3);
                minima = zeros(dim3);
                imsegm = zeros(dim3);
            end;
            clear D;
            msg = ['Image dimension: ', num2str(dim)];
            disp(msg);

                   
        catch
            msg = ['Could not load ' name];
            disp(msg);
% %             close(fig.handle);
% %             close(control.handle);
            i = i + 1;
%             if i == stop + 1
%                 finish = 1;
%             end;
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
    fig.w = w;
    fig.h = h;

%     % change size if necessary
%     p = get(fig.handle,'OuterPosition');    
% %     fig.pos = [p(1)+4 p(2)+4 w h];
%     set(fig.handle,'Position',p);
% %     control.pos = [p(1)+4-control.w p(2)+4 control.w h];
%     p = get(control.handle,'OuterPosition');    
%     set(control.handle,'Position',p);

    
    % show images
    prm.fig = fig;
    prm.plane = plane;
    prm.dim = dim;
    showimagecells(im,imsegm,cellbw,wat,minima,prm);
    
    % wait for user
    uiwait(control.handle);
    
    if choice == 1
        plane = min(plane + 1,dim(3));
    elseif choice == 2
        plane = max(plane - 1,1);
    elseif choice == 3
        i = i + 1;
    elseif choice == 4
        i = i - 1;
    elseif choice == 5
        planenew = inputdlg('Enter plane number: ');
        if ~isempty(planenew),
            plane = str2double(planenew);
        end;
    elseif choice == 6
        % user defined input
        figure(fig.handle);
        [y x] = ginput(1);
        x = round(x);
        y = round(y);
        val = wat(x,y);
        nth = numel(prmin.prm.thname);

        A = {'Classificator','Value','Threshold','Cell'};
        for j = 1 : nth
            name = prmin.prm.propname{j};
            A{j+1,1} = name;
            A{j+1,2} = num2str(prmin.prop.(name)(val));
            A{j+1,3} = num2str(double(prmin.prm.(prmin.prm.thname{j})));            
            A{j+1,4} = num2str(double(prmin.iscellhere(val,j)));            
        end;
        A{nth+2,1} = 'All classificators must be 1 to become a cell';
        
         
        % print the string to the figure 
        b = control.title.pos(2) - 2*button.h;
        h = button.h;

        w = button.w;        
        for j = 1 : size(A,1)-1
            
            left = button.left;
            for k = 1 : size(A,2)                
                control.message.handle(j,k) = uicontrol('Parent',control.handle,'Style','Text','String',A{j,k}, ...
                'Position',[left b w h]);   
                left = left + button.w;           
            end;
            b = b - button.h;
        end;
        % last message string
        left = button.left;
        w = button.h*4;
        control.message.handle(j+1,1) = uicontrol('Parent',control.handle,'Style','Text','String',A{nth+2,1},'Position',[left b w h]);
        
        
    elseif choice == 7
        finish = 1;
    end;

    
    msg = ['Plane ' int2str(plane)];
    disp(msg);
  
    if finish == 1
        close(fig.handle);
        close(control.handle);
        break
    end;
    
    
end;


%--------------------------------------------------------------------

function [] = showimagecells(im,imsegm,cellbw,wat,minima,prm)



plane = min(prm.plane,prm.dim(3));
plane = max(plane,1);

% these planes
if isequal(prm.visch1,'imsegm')
    implane = 0.5*scale(imsegm(:,:,plane));
else
    implane = 0.5*scale(im(:,:,plane,prm.visch1));
end;
% implane = im(:,:,plane,prm.visch1);

cellbwplane = double(cellbw(:,:,plane));
perimplane = bwperim(cellbwplane);
% load ball1;se = getball(ball,1,1);
% perimhere = imdilate(perimhere,se);
watplane = wat(:,:,plane);
if isequal(prm.visch2,'imsegm')
    im2plane = scale(imsegm(:,:,plane));
else
    im2plane = scale(im(:,:,plane,prm.visch2));
end;
cellbwplane(watplane == 0) = 0.5;

implanem = implane;
rgboverlayplane = implane;
rgboverlayplane(perimplane) = 1;
rgboverlayplane(:,:,2) = implane;
ind = minima(:,:,plane) == 1;
implanem(ind) = 1;
rgboverlayplane(:,:,3) = implanem;

implane = uint8(round(implane*255));
im2plane = uint8(round(im2plane*255));
rgboverlayplane = uint8(round(rgboverlayplane*255));

% figure(prm.fig.handle);
% offset = prm.button.hn;
% d = 0.02;
% h = (1-offset)/2 - 4*d;
% m = margin.leftn;
% mm = margin.rightn;

figure(prm.fig.handle);

% b = offset + d;
b = 0;
h = 0.45;
w = 0.45;
l = 0.025;
pos = [l b w h];
subplot('Position',pos);imagesc(rgboverlayplane);colormap(gray);drawnow;axis off;
axis image;title('Overlay image, markers (blue), segmentation (red)');

pos = [0.5+l b w h];
subplot('Position',pos);imagesc(cellbwplane);colormap(gray);drawnow;axis off;
axis image;title('Found cells');

b = b + 0.5;
pos = [l b w h];
% implane(implane < 100) = 0;
subplot('Position',pos);imagesc(implane);colormap(gray);drawnow;axis off;
if isequal(prm.visch1,'imsegm')
    axis image;title(['Channel ' prm.visch1])
else
    axis image;title(['Channel ' int2str(prm.visch1)])
end;

pos = [0.5+l b w h];
subplot('Position',pos);imagesc(im2plane);colormap(gray);drawnow;axis off;
if isequal(prm.visch2,'imsegm')
    axis image;title(['Channel ' prm.visch2])
else
    axis image;title(['Channel ' int2str(prm.visch2)])
end;



%----------------------------------------------------------------

function [] = showimagetnt()

global plane vischglobal im cellbw tntcand wat overlaynum fignum external
global plane vischglobal im cellbw tntcand wat overlaynum fignum external



% these planes
% imhere = 0.5*scale(im(:,:,plane,vischglobal));
imhere = (im(:,:,plane,vischglobal));

cellbwhere = double(cellbw(:,:,plane));
cellbwhere = double(cellbw(:,:,plane));
perimhere = bwperim(cellbwhere);
perimhere = bwperim(cellbwhere);


% TNT candidates
tntcandhere = double(tntcand(:,:,plane));

% to show the watershed that were not classified as cells
wathere = wat(:,:,plane);
cellbwhere(wathere == 0) = 0.5;
cellbwhere(wathere == 0) = 0.5;

% the TNTs in this plane
tntcandhere = logical(tntcand(:,:,plane));

cellbwhere2 = double(cellbw(:,:,plane));
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


