% STRUCTSMOOTH Structural smoothing of image
% STRSMOOTH(IM,D,HX,HY,HZ) Smoothes the image IM using a square filter with
% diameter D and stepsize STEPS which must be a 1x3 array
% 
% Also possible to run with spedifying the dimensio of the filter. This is
% important when swithching between PC12 cells nad NRK cells
% STRSMOOTH(IM,D,HX,HY,HZ,[DIMXY DIMZ STDEV]) DIMXY and DIMZ specifies the
% dimension of the filter to smooth the image by a gaussian convolution.
% STDEV is the standarad deviation in the Gaussian
%
% For NRK cells, set to DIMXY = 7, DIMZ = 3, STDEV = 5;
%
% Ex: filtim = structsmooth(im,13,1,1,3);
% 
% High-throughput Anaysis of Multispectral Imagse of Breast Cancer Tisue.
% Umesh Adiga et al. IEEE Transactions on image processing, Vol 15, No 8,
% August 2006
%
function [filtim] = dircohenh(varargin)

im = varargin{1};
d = varargin{2};
h = varargin{3};
hx = h(1);hy = h(2);hz = h(3);
gpu = 0;
if nargin == 4
    gpu = varargin{4};
end;

% dimxy = 3;
dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;

% gpu = 1;

% if gpu
%     msg = ['Using GPU code'];
%     disp(msg);
% else
%     msg = ['Using CPU code'];
%     disp(msg);    
% end;
% number of elements to remove before one takes the Olympic average. Set to
% 1-->3
numremele = 2;


%
% Make the filters
%
stepxy = 30;
% stepxy = 30;
% stepz = 45;
% stepz = 30;
deghere = (0:stepxy:(180-stepxy))';
numdeghere = length(deghere);
deg2D = [deghere repmat(90,numdeghere,1)];


% deghere  = (0 : stepxy : (360-stepxy))';
% numdeghere = length(deghere);
% deg3D = [deg2D ; ...
%          0 0 ; ...
%          deghere repmat(stepz,numdeghere,1)];

% this is quicker, not 100% correct but works as good
deg3D = [deg2D; ...
         0 0];
     
if dim(3) == 1
%     disp('Using 2D directional coherence enhancement filter')
    deg = deg2D;
else
%     disp('Using 3D directional coherence enhancement filter')
    deg = deg3D;
end;


% the half filter size
p = floor(d/2);
   
% make filter
c = makefilter(deg,p);

% number of directions
numdir = size(c,2);
 
if gpu
    % for Jacket speedup
    im = gsingle(im);
end;

% inisization
filtim = zeros(dim);
for i = 1 : numdir
    
    % coordinates of points in this direction after rotation of filter
    chere = c(i);

    % filter image for max value
    [maxim,sumim,numpoints] = maxfilt3(im,chere,numremele,dim,hx,hy,hz,gpu);

    % structural filtering   
    filtim = max(filtim,(sumim - sum(maxim,4))/ (numpoints-numremele));

end;

if gpu
    % gather from GPU
    filtim = double(filtim);
end;

% ----------------------------------------------------------

function [maxim,sumim,numpoints] = maxfilt3(im,c,numremele,dim,hx,hy,hz,gpu)

% % relative coordinate where we want the image values
% cx = c.x;
% cy = c.y;
% cz = c.z;

% grid
if dim(3) == 1
%     [x y] = ndgrid(hx:hx:hx*dim(1),hy:hy:hy*dim(2));
    [x y] = ndgrid(1:dim(1),1:dim(2));
    z = ones(dim(1:2));
else    
%     [x y z] = ndgrid(hx:hx:hx*dim(1),hy:hy:hy*dim(2),hz:hz:hz*dim(3));
    [x y z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
end;

if gpu
    x = gsingle(x);
    y = gsingle(y);
    if dim(3) > 1
        z = gsingle(z);
    end;
end;


maxim = zeros([dim numremele]);
sumim = zeros(dim);
numpoints = length(c.x);
for i = 1 : numpoints

    % the coordinate where we want to find the values
    xhere = x + c.x(i);
    yhere = y + c.y(i);
    zhere = z + c.z(i);

    % NB!!!!!!!!!
    % Must not do this for z direction, then it goes wrong.
    % Do NOOOOOT replacestr any Inf or NaN by image values, then the filter
    % does not work anymore!!!!!!
%     xhere(xhere < hx) = hx;xhere(xhere > dim(1)*hx) = dim(1)*hx;
%     yhere(yhere < hx) = hx;yhere(yhere > dim(2)*hx) = dim(2)*hy;
    xhere(xhere < 1) = 1;xhere(xhere > dim(1)) = dim(1);
    yhere(yhere < 1) = 1;yhere(yhere > dim(2)) = dim(2);
    
    
    if dim(3) == 1
        if gpu
            imhere = interp2(im,yhere,xhere);        
        else
            imhere = interp2(y,x,im,yhere,xhere,'cubic',-Inf);
        end;
        
    else        
        if gpu
            imhere = interp3(im,yhere,xhere,zhere);
        else
            imhere = interp3(y,x,z,im,yhere,xhere,zhere,'cubic',-Inf);
        end;
    end;
    
    
    % sort previous max image and 
    sortarray = maxim;
    sortarray(:,:,:,numremele+1) = imhere;
    sorted = sort(sortarray,4,'descend');

    % max image
    maxim(:,:,:,1:numremele) = sorted(:,:,:,1:numremele);    
    sumim = sumim + imhere;

    
end;
%       showall(maxim(:,:,:,1),maxim(:,:,:,2),sumim)
    
%-------------------------------------------

function [c] = makefilter(deg,p)

% step = round(p/2);
% r = -p:step:p;
% this is a bit tricky.....??? is this optimal???
r = linspace(-p,p,7);

% degrees, sphaerical coordinates
% theta : angle to positive x-axis 
% phi   : angle to positive z-axis 
% r     : distance along vector to point
%
for i = 1 : size(deg,1)

    % the degrees
    theta = deg(i,1);
    phi = deg(i,2);        

    % the positions
    c(i).x = r*sind(phi)*cosd(theta);
    c(i).y = r*sind(phi)*sind(theta);
    c(i).z = r*cosd(phi);

end;

%--------------------------------------------

