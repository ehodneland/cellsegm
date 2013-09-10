% CONSTBORDER Placing a constant value around the border of image
%
%   CONSTBORDER(IM,P,V) Placing value V around xy-border as far as P in 
%   image IM. Image must be 2D image, if more than three planes use 
%   constborder3D.
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
function [im] = constborder(im,p,v)

[M N O] = size(im);
blockones = zeros(M,N,O);
blockones(p+1:end-p,p+1:end-p,:) = 1;
blockones = imcomplement(blockones);
im(blockones == 1) = v;
