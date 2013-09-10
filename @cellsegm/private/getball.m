% GETBALL Get ball for morphological operations
%
%   BALL2 = GETBALL(BALL,R,O) is returning the correct size of the 
%   structuring element based on the radius of it and the number of planes.
%   BALL  : structural element
%   R     : radius of structural element
%   O     : number of planes
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
function [ball2] = getball(ball,r,O)

if O > 2*r+1
    ball2 = ball;
else
    fra = r-floor(O/2)+eq(round(O/2),O/2) +1;
    til = r+floor(O/2)+1;
    ball2 = ball(:,:,fra:til);
end;




