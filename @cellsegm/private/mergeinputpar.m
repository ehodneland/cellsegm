function [prm] = mergeinputpar(prm,prmin)
% MERGEINPUTPAR Merge input parameters in struct arrays
% 
%   PRM = MERGEINPUTPAR(PRM,PRMIN) merging input parameters in PRM and
%   PRMIN where PRM is the default struct. Returning the merged parameters.
%   MERGEINPUTPAR is different from MERGESTRUCT as it can go down several
%   levels in the struct array.
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

if isempty(prmin)
    return;
end;

f = fieldnames(prmin);
for i = 1 : numel(f)
    prm2 = prmin.(f{i});
    if isstruct(prm2)
        % if the field does not exist we do not overwrite and create it
        % with the new field
        if ~isfield(prm,f{i})
            prm.(f{i}) = prm2;
        else
            prm.(f{i}) = mergeinputpar(prm.(f{i}),prm2);
        end;
    else
        prm.(f{i}) = prm2;
    end;    
end;


