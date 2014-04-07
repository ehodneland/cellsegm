function [absgrad] = absgradient(im,h)
% ABSGRADIENT Computing the absolute gradient image.
%
% ABSGRAD = ABSGRADIENT(IM,H) computes the absolute gradient in 3D for
% image IM and stepsize H
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


[ux,uy,uz] = derivcentr3(im,h(1),h(2),h(3));
absgrad = sqrt(ux.^2 + uy.^2 + uz.^2);

