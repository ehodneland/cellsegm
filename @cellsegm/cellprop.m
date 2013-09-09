function [prop] = cellprop(im,wat,cellv,propname,h,meanintbck)
% CELLPROP Computing cell properties
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
msg = ['This is ' upper(mfilename) ' computing cell properties'];
disp(msg);


voxelvol = prod(h);
for i = 1 : numel(cellv)
    
    dim = size(wat);
    if numel(dim) == 2
        dim = [dim 1];
    end;
    
    % this region
    reghere = wat == cellv(i);
    
    % middle plane region
    [range,minz,maxz] = bwrange(reghere);
    mid = round(mean([minz maxz]));        
    regpl = double(reghere(:,:,mid));
    voxvolregpl = nnz(regpl);    
    
    % perimeter
    perim = zeros(dim);
    for j = 1 : dim(3)
        perim(:,:,j) = bwperim(reghere(:,:,j));
    end;
        
    % must always compute volume
    prop.volume(i,1) = nnz(reghere)*voxelvol;
        
    name = 'intincell';
    if ismember(name,propname)
        v = im(reghere == 1);
        v = mean(v(:))/meanintbck;
        prop.intincell(i,1) = v;
    end;
    
    name = 'intborder';
    if ismember(name,propname)
        v = im(perim == 1);
        v = mean(v(:))/meanintbck;
        prop.intborder(i,1) = v;        
    end;
    
    name = 'convexarea';
    if ismember(name,propname)
        prop.convexarea(i,1) = 0;
        a = regionprops(double(regpl),name);
        if ~isempty(a)
            prop.convexarea(i,1) = prop.convexarea(i,1) + a.ConvexArea;
        end;
        prop.convexarea(i,1) = prop.convexarea(i,1)/voxvolregpl;
    end;
    
    
    name = 'convexperim';
    if ismember(name,propname)
        prop.convexperim(i,1) = 0;
        [regmax,box] = boundreg(regpl,0,0);
        a = regionprops(regmax,'ConvexImage');

        load ball1;se = getball(ball,1,1);
        regmaxconvperim = bwperim(imerode(a.ConvexImage,se));   

        % find perimeter which is convex
        regmaxperim = bwperim(regmax);

        % the convex parts of the boundary
        load ball1;se = getball(ball,1,1);
        regmaxconvperim = imdilate(regmaxconvperim,se) .* regmaxperim;

        % the ratio of convex boundary
        prop.convexperim(i,1) =  nnz(regmaxconvperim) / nnz(regmaxperim);
            
    end;
    
    
end;


