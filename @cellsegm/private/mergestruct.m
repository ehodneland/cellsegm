function [prmdefault] = mergestruct(prmdefault,prmin)
% MERGESTRUCT Merging two structs
%
%   PRMDEFAULT = MERGESTRUCT(PRMDEFEAULT,PRMIN) merging the two structs 
%   PRMDEFAULT and PRMIN where the fields in PRMDEFAULT are overwritten by 
%   those in PRMIN.
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

    
    
    if isempty(prmin)        
        return;
    end;
        
    % names of parameter fields
    varnames = fieldnames(prmin);
    for i = 1 : size(varnames,1)
        varnamehere = varnames{i};
        valhere = eval(['prmin.' varnamehere]);
        % we keep the given value!!
        if isequal(valhere,'default')
            continue;
        end;
        evalc(['prmdefault' '.' varnamehere '=' 'prmin' '.' varnamehere]);
    end;
