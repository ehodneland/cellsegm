function [] = printcell(fid,A,thmin)
% PRINTCELL Print cell to file
% PRINTCELL(FID,A,TH) Print cell A to file fID, length of TH between
% numbers. The format needs to be string
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

[M N] = size(A);
% check the length
th = Inf(M,N);
for i = 1 : M
    for j = 1 : N
        val = A{i,j};
        if isnumeric(val)
            val = num2str(val);
        end;
        th(i,j) = max(length(val)+2,thmin);
    end;
end;

th = max(th,[],1);

for i = 1 : M
    for j = 1 : N
        val = A{i,j};
        if isnumeric(val) || islogical(val)
            if round(val) == val
                val = num2str(val);
            else
                val = num2str(val,'%6.6f');
            end;
        end;
        val = makestr(val,th(j));
        fprintf(fid,'%s',val);      
    end;
    fprintf(fid,'\n');
end;


