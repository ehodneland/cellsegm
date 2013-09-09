function [] = runexamples(option)
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

% RUNEXAMPLES Run examples in CellSegm
%
% RUNEXAMPLES(OPTION) runs selected examples as in option numbering.
%
% RUNEXAMPLES('all') runs all examples
%
if isequal(option,'all')
    for i = 1 : 10        
        run(i);
        if i < 10
            msg = ['Option ' int2str(i) ' done, press any key'];
            disp(msg);
            pause
        end;
    end;
else
    run(option);
end;


function [] = run(option)

if isequal(option,1)
    surfstain_smoothing_2D    
    
elseif isequal(option,2)  
    nucleistain_2D                           
    
elseif isequal(option,3)
    nucleistain_3D                             
    
elseif isequal(option,4)
    nucleistain_smoothing_2D                   
    
elseif isequal(option,5)
    surfstain_2D                               
    
elseif isequal(option,6)
    surfstain_3D                               
    
elseif isequal(option,7)
    surfstain_and_manual_3D                    
    
elseif isequal(option,8)
    surfstain_and_nucleus_2D                   
    
elseif isequal(option,9)
    surfstain_and_nucleus_3D                   
    
elseif isequal(option,10)
    surfstain_and_nucleus_cellsegmentation_3D  
        
end;