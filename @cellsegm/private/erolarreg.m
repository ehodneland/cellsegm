% EROLARREG(BW,P,TH) Eroding binary large regions above 
% TH pixels using disc-like 
% structuring element with radius p.
%
function [BW] = erolarreg(BW,p,th)


[faser,L] = bwlabeln(BW);
large = zeros(size(BW));

for j = 1 : L
    regHere = eq(faser,j);
    ind = find(regHere);
    range = bwrange(regHere);
    numPix(j) = size(ind,1)/range;
    if numPix(j) > th
        large(ind) = 1;
        BW(ind) = 0;
    end;    
end;


% erosion
name = ['ball' int2str(p)];
load(name);se = getball(ball,p,1);
large = imerode(large,se);

%show(large,20)
%show(BW,21)

% combine small and large
BW = logical(BW + large);
