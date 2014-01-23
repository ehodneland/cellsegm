% ROWVECT Making row vector 
%
function [a] = rowvect(a)

[M N] = size(a);
if M > N;
    a = a';
end;