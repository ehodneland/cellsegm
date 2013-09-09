% GAUSSIAN Make a Gaussian filter
% GAUSSIAN(DIM,SIGMA) Make a Gaussian filter of dimension in xy plane
% of DIM and of z dimension DIMZ. The standard deviation is SIGMA. 
%
function [g] = gaussian(dim,sigma)



dim = round(dim);
midx = dim(1)/2;
midy = dim(2)/2;
midz = dim(3)/2;


[x y z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
x = x - midx-0.5;
y = y - midy-0.5;
z = z - midz-0.5;

if dim(3) == 1
    g = (1/(2*pi*sigma^2))*exp(-(x.^2 + y.^2)/(2*sigma.^2));
elseif dim(3) > 1
    g = (1/(2*pi*sigma^2)^(3/2))*exp(-(x.^2 + y.^2 + z.^2)/(2*sigma.^2));
else
    error('Wrong option for DIMZ')
end;
