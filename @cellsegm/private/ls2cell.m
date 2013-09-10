% LS2CELL Write folder and file names to a cell array
%
%   LS2CELL(OPTION,PATTERN) writes file or folders matching PATTERN to 
%   a cell array. OPTION = 1 writes file names, OPTION = 2 writes folder 
%   names matching the pattern.
%   Returning a cell array with the strings as cells.
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
function [s] = ls2cell(option,pattern)


d = dir(pattern);

[A B C] = fileparts(pattern);
if ~isempty(A)
    for i = 1 : numel(d)
        d(i).name = [A '/' d(i).name];
    end;
end;

c = 1;
s = '';
for i = 1 : length(d)
    a = d(i).name;    
    if option == 1
        if exist(a,'file') == 2
            s{c,1} = a;
            c = c + 1;
        end;
    elseif option == 2
        if exist(a,'dir') == 7
            s{c,1} = a;
            c = c + 1;
        end;
    end;    
end;