% MULTCONSTCELL Multiplies a cell matrix with a constant
% B = MULTCONSTCELL(A,C) multiplies the cell matrix or cell vector A with a
% constant C. Returns the multiplied cell matrix A;
%
function [A] = multconstcell(A,c)

[m n] = size(A);
for i = 1 : m
    for j = 1 : n
        A{i,j} = A{i,j}.*c;
    end;
end;