% SHOW Drawing images using colormap gray and imagesc function.
%
%   SHOW(IM) is drawing the image in IM at the middle plane
%
%   SHOW(IM,FIGNUM) is drawing the image in IM at the middle plane in the
%   figure FIGNUM
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
function [] = show(varargin)

global inactiv    

if nargin == 1    
    I = varargin{1};
    [M N O] = size(I);
    middle = round(O/2);
    if isempty(inactiv)
        figure;colormap(gray);imagesc(I(:,:,middle));axis image;drawnow
    end;
elseif nargin == 2
    I = varargin{1};
    [M N O] = size(I);    
    number = varargin{2};
    middle = round(O/2);
    if isempty(inactiv)    
        figure(number);colormap(gray);imagesc(I(:,:,middle));axis image;drawnow
    end;
end;%if    
