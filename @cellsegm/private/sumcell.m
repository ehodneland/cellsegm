% SUMCELL Summing to cell arrays A and B
% C = SUMCELL(A,B) Summing the cell arrays A and B
%
function [c] = sumcell(a,b)


if iscell(a)    
  
    [rows1,cols1] = size(a);
    [rows2,cols2] = size(b);
    if ~isequal(rows1,rows2) || ~isequal(cols1,cols2)
        error('Wrong number of elements')
    end;
    rows = rows1;
    cols = cols1;
    
    c = cell(rows,cols);
    for i=1:rows
        for j=1:cols        
            c{i,j} = a{i,j} + b{i,j};
        end
    end  
  
else
    c = a+b;
end