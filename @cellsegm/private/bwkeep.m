% BWKEEP Keeping the largest regions in BW
%
% BW = BWKEEP(BW,VAL,CONN) keeps the VAL largest regions in BW using the
% connectivity in CONN
%
function [bw] = bwkeep(bw,val,conn)


[dim,faser] = bwsize(bw,conn);
bw = zeros(size(bw));
[sorted,ind] = sort(dim,'descend');

numreg = length(dim);
% val = unique(faser(faser > 0)));
for i = 1 : val
    if i > numreg
        break;
    end;
    valhere = ind(i);
    reghere = faser == valhere;

    bw(reghere == 1) = 1; 
end;