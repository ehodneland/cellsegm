% BWSIZE Finds the size of binary segments
%
% [DIM,FASER] = BWSIZE(BW,CONN) finds the size of the disconnected binary
% segments inside BW and returns the number of pixels in DIM and the
% piecewise constant image FASER connected to the label in DIM
%
function [dim,faser] = bwsize(BW,conn)

[faser,L] = bwlabeln(BW,conn);
val = faser(faser > 0);
valu = unique(val);
n = numel(valu);
dim = zeros(L,1);
for i = 1 : n
    dim(i) = nnz(val == valu(i));
    
%     % this region
%     reg = faser == i;
    
%     % number of pixels
%     dim(i) = nnz(reg);
    
end;
