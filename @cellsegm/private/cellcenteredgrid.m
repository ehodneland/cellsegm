function [x,minx,maxx] = cellcenteredgrid(dim,h)
% CENTERGRID Find a cell-centered grid of dimension DIM and voxelsize H
%
%
ndim = length(dim);

% natural coordinates
x = cell(ndim,1);
for i = 1 : ndim
    
    % make grid
    x{i,1} = (0.5:1:dim(i)-0.5)*h(i);

end;

% grid
if ndim == 2
    [x{1,1},x{2,1}] = ndgrid(x{1},x{2});
elseif ndim == 3
    [x{1,1},x{2,1},x{3,1}] = ndgrid(x{1},x{2},x{3});
elseif ndim == 4
    [x{1,1},x{2,1},x{3,1},x{4,1}] = ndgrid(x{1},x{2},x{3},x{4});
else
    error('Wrong dimension');
end;

% max and min
maxx = NaN(ndim,1);
minx = NaN(ndim,1);

for j = 1 : ndim
    maxx(j) = max(x{j}(:));
    minx(j) = min(x{j}(:));
end;    
