% IMRESIZE3D Resize image in 3D
% IMRESIZE3D(IM,G,METHOD) Resize image IM in 3D by the factor G. 
% If G has more than one element, it gives the image dimension of the
% rescaling.
% METHOD is the interpolation method, as in INTERP3
%
function [u] = imresize3d(im,g,method)

D = size(im);
ndim = numel(D);

% if IM is 2D
if numel(D) == 2  
    if numel(g) > 2
        g = g(1:2);
    end;
    
    try
        u = imresize(im,g,method);    
    catch
        u = imresize(im,g);    
    end;
    
    return;
end;

if length(g) == 1;
    d = round(D*g);
elseif length(g) == 2    
    d = g;
    d(3) = 1;    
elseif length(g) >= 3
    d = g(1:3);
end;
d = round(d);

%
% We must do 3D interpolation
%

% relative steplength in input image
H = ones(3,1);

% relative steplength in returned image
h = NaN(3,1);
for i = 1 : 3
    h(i) = D(i)./d(i);
end;

% the original grid
X = cell(3,1);
for i = 1 : 3
    X{i} = 0.5:H(i):D(i)-0.5;
end;
% the new grid
% NOTE: We start with h/2 so we get the center of the pixel in the right
% position; earlier we startet in 1, that is wrong!
x = cell(3,1);
for i = 1 : 3
    x{i} = h(i)/2 : h(i) : D(i)-h(i)/2;
    % NB we need to do this for when we scale UP, then the interpolated
    % values are outside the grid!
    minx = min(X{i});
    maxx = max(X{i});
    x{i}(x{i} < minx) = minx;
    x{i}(x{i} > maxx) = maxx;
end;
[x{1},x{2},x{3}] = ndgrid(x{1},x{2},x{3});
[X{1},X{2},X{3}] = ndgrid(X{1},X{2},X{3});

% imnew = im;
% for i = 1 : size(im,3)
%     imnew(:,:,i) = medfilt2(im(:,:,i),[3 3]);
% end;
% im = imnew;

% interpolation
if isequal(method,'bilinear')
    method = 'linear';
end;

gpu = 0;
if isa(im,'gsingle') || isa(im,'gdouble')
    gpu = 1;
end;
if gpu    
    for i = 1 : ndim
        x{i} = gsingle(x{i});
        X{i} = gsingle(X{i});
    end;
    u = interpgpu(X,im,x);
else
    u = interp3(X{2},X{1},X{3},im,x{2},x{1},x{3},method,Inf);    
end;

if gpu
    if isempty(find(im)) && ~isempty(find(u))
        warning('The resized image is empty');    
    end;
else
    if isempty(find(im,1)) && ~isempty(find(u,1))
        warning('The resized image is empty');    
    end;
end;


















