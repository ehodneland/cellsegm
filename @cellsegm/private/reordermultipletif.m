function [im] = reordermultipletif(imhere,nch)

dim = size(imhere);
if numel(dim) == 2
    dim = [dim 1];
end;
nplane = dim(3)/nch;

if round(nplane) ~= nplane
    warning('Wrong number of channels')
    im = [];
    return;
end;
dim = [dim(1:2) nplane nch];            
im = zeros(dim);
% reorder data
c = 0;
for k = 1 : nplane
    for l = 1 : nch
        c = c + 1;
        im(:,:,k,l) = imhere(:,:,c);
    end;
end;
