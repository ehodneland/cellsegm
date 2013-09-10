% SHOWALL Show one plane awt a time for stack.
%
%   SHOWALL(IM) Shows all planes in a 3D stack succesively when hitting ENTER. 
%
%   SHOWALL(IM1,IM2,...) There can be several input images, then they are 
%   plottet in different figures.  They have to have the same number of 
%   planes. If the last argument is 'colorbar', then the images are plotted 
%   with the colorbar option. 
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
function [] = showall(varargin)

numim = 0;
col = 0;
% plane = 1;
for i = 1 : nargin    
    varhere = varargin{i};
    
    if isequal(varhere,'colorbar')
        col = 1;
        continue;
    end;
    
    f{i} = double(varargin{i});
    numim = numim + 1;
end;
        



F = f{1};
[M N O P] = size(F);

if P > 1
    t = 1;
    while 1        
        for j = 1 : numim
            F = f{j};
            show(F(:,:,:,t),j);            
        end;
        msg = ['Time point ' int2str(t)];
        disp(msg);
        
        t = input('Continue: 0, Time: timepoint, Next timepoint: Enter ');

        if isequal(t,0)
            return
        elseif isempty(t)
            t = t + 1;
        else
            t = t;
        end;                        

        if t > O
            break;
        end;

    end;
else


    niter = 1;
    while 1

        for j = 1 : numim
            F = f{j};

            if col == 0
                show(F(:,:,niter),j)        
            else
                show(F(:,:,niter),j);colorbar
            end;
        end;
        disp(sprintf('Plane %i',niter))

        c = input('Continue: 0, Plane: planenumber, Next plane: Enter ');

        if isequal(c,0)
            return
        elseif isempty(c)
            niter = niter + 1;
        else
            niter = c;
        end;                        

        if niter > O
            break;
        end;
    end;
end;
    