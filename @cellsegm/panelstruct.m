function [h] = panelstruct(varargin)
% PANELSTRUCT Make panel of images
% H = PANELSTRUCT(IM,D,H) makes panel of images 
% in structural array IM with
% spacing D and image height in pixels H
% 
% Example:
% panelstruct(im,0.005,800)
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


im = varargin{1};
d = varargin{2};
H = varargin{3};
cbar = zeros(size(im));
for i = 4  : 2 : nargin
    namehere = varargin{i};
    varhere = varargin{i+1};
    switch namehere
        case 'lim'
            lim = varhere;
        case 'colorbar'
            cbar = varhere;
    end;
end;

[m n] = size(im);
[M N O] = size(im{1});

% this is right:
t1 = (m * M);
t2 = (n * N);
% t1 = t1 + d*t1*(m+1);
% t2 = t2 + d*t2*(n+1);
ratio = t2/t1;

% stepcol = (1-(n+1)*d)/n;
% steprow = (1-(m+1)*d)/m;

stepcol = N/t2;
steprow = M/t1;

pos = [70 70 H*ratio H*1];
h = figure('position',pos);

% row index
posrow = 0;
for i = m:-1:1    
    poscol = 0;    
    % column index
    for j = 1:n        
        imhere = im{i,j};        

        pos = [poscol posrow stepcol-d steprow-d];        
        if ~exist('lim','var')
            limhere = [min(imhere(:)) max(imhere(:))]; 
        else
            limhere = lim{i,j};
        end;
        limhere = double(limhere);
        limhere(isnan(limhere)) = 0;
        dim = size(imhere);
        ndim = numel(dim);
        if ndim > 2            
            figure(h);subplot('position',pos);imagesc(imhere);axis image;axis off;
        else
            figure(h);subplot('position',pos);imshow(imhere,limhere);axis image;axis off;colormap(gray);
        end;
        if cbar(i,j) == 1
            colorbar('East');
        end;
        poscol = poscol + stepcol + d/4;
        
    end;
    posrow = posrow + steprow;
end;

set(h,'PaperPositionMode','auto') 

% F = getframe(h);
% imwrite(F.cdata,figname,format);
% saveas(h,figname,format);
% print(h,format,figname);

