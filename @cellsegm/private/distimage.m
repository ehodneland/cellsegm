% DISTIM Partition the image according to the Euclidian distances
% BOUNDIM = DISTIM(MINIMA,CONN) Partition the image according to the 
% Euclidian distances from the markers in MINIMA to other pixels in image.
% A pixel belongs to a certain phase if the distance to that marker is the
% smallest of distances to all markers. CONN is the connectivity that is
% used in labelling the minima
%
function [a] = distimage(minima,conn)

if isempty(find(minima,1))
    a = zeros(size(minima));
    return;
end;

[M N O] = size(minima);

% label the minima
[faser,L] = bwlabeln(minima,conn);

for i = 1 : L
    regHere = eq(faser,i);
%     showall(regHere)
    di(:,:,:,i) = bwdist(regHere);
end;

[minVal,a] = min(di,[],4);
a = squeeze(a);
