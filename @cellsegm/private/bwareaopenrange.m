% BWAREAOPENRANGE Removing small objects in binary image
%
% BW = BWAREAOPENRANGE(BW,TH,CONN) Removing small regions in binary 
% image BW below the threshold TH, scaled with the range of each 
% connection. CONN is he connectivity. Returning the remaining objects in
% image BW
%
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
function [bw] = bwareaopenrange(bw,th,conn)

[faser,L] = bwlabeln(bw,conn);


for j = 1 : L
    reghere = eq(faser,j);
    ind = find(reghere);
    range = bwrange(reghere);
    numpixrel = numel(ind) / range;
    
    if numpixrel < th
        bw(ind) = 0;
    end;
end;

