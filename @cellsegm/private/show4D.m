function [ ] = show4D( varargin )
%  SHOW4D Show 4D data
%
%   SHOW4D(IM) Show 4D data IM inteartively
%
%     =======================================================================================
%     Copyright (C) 2013  Erlend Hodneland
%     Email: erlend.hodneland@biomed.uib.no 
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     =======================================================================================
%

for i = 1 : nargin
    f{i} = varargin{i};
end;

global prm

for i = 1 : nargin
    h(i) = figure(i);
end;


dim = ones(1,4);
a = size(f{1});
n = numel(a);
dim(1:n) = a;
for i = 1 : nargin
    dimhere = size(f{i});
    try
        dim(4) = max(dim(4),dimhere(4));
    catch
    end;
end;

prm.plane = round(dim(3)/2);
prm.time = 1;
prm.cont = 0;

while 1

    
    prm.plane = min(prm.plane,dim(3));
    prm.plane = max(prm.plane,1);
    prm.time = min(prm.time,dim(4));
    prm.time = max(prm.time,1);
    msg = ['Plane: ' int2str(prm.plane) ', time: ' int2str(prm.time)];
    disp(msg);
    
    for i = 1 : nargin
        figure(h(i));
        try
            imagesc(f{i}(:,:,prm.plane,prm.time));colormap(gray);axis image;
        catch
            imagesc(f{i}(:,:,prm.plane));colormap(gray);axis image;
        end;
    end;

    choice = menu('show4D','Up','Down','Next','Previous','Quit');
    
    if choice == 1
        prm.plane = prm.plane + 1;
    elseif choice == 2
        prm.plane = prm.plane - 1;
    elseif choice == 3
        prm.time = prm.time + 1;
    elseif choice == 4
        prm.time = prm.time - 1;
    elseif choice == 5
        for i = 1 : numel(h)
            close(h(i));
        end;
        return;
    end;
    
            
    
    % reset
    prm.cont = 0;

    
end;



