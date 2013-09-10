% FILLHOLES Fill holes with a given size
%
%   FILLHOLES(BW,LOWTH,HIGHTH,CONN) Fill holes in BW being within 
%   the range of LOWTH and HIGHTH and with connectivity CONN
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
function [bw] = fillholes(bw,lowth,highth,conn)

[M N O] = size(bw);

[faser,L] = bwlabeln(bw);

holes = zeros(size(bw));
holesall = holes;
for i = 1 : L
    
    % this region
    reghere = eq(faser,i);
    
    % fill this region
    for j = 1 : O        
        holes(:,:,j) = imfill((reghere(:,:,j)),conn,'holes') - reghere(:,:,j);        
    end;
    holes = logical(holes);
    holesall(holes == 1) = 1;
end

holesall = bwareaopen(holesall,round(lowth));
holesall = holesall - bwareaopen(holesall,round(highth));

% fill
bw(holesall == 1) = 1;