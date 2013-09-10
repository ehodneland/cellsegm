function [] = printstructscreen(var)
% PRINTSTRUCTSCREEN prints a struct array to the screen
%
%    PRINTSTRUCTSCREEN(VAR) prints the struct array in VAR to the screen
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
p = 20;
f = fieldnames(var);
n = 0;
for i = 1 : length(f)
    n = max(n,length(f{i}));    
end;

for i = 1 : length(f)
    val = var.(f{i});
    if numel(val) > 20 && isnumeric(val)
        val = [num2str(val(1)) ' to ' num2str(val(end))];
    end;
    if isnumeric(val)
        val = num2str(val);
    end;

    if iscell(val) || isstruct(val)
        continue;
    end;
    try
        msg = [makestr(f{i},n+3) makestr(':',n+3) makestr(val,n+3)];
        disp(msg);
    catch
        
    end;
    
end;