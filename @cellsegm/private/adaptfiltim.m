% ADAPTFILTIM(IM,RAD,D) Adaptive filtering on IM, RAD defines the filter size 
% of average filter, try 10. D is
% the additive value above background value, can be set to 0.02-0.2,  
% Image should be scaled to [0 1] before entering the algorithm.
% ADAPTFILT(IM,RAD,D,HX,HZ) includes the pixel resolution also
% 
% Ex adaptfiltim(im,30,0.05);
% 
function [thIm] = adaptfiltim(varargin)

im = varargin{1};
rad = varargin{2};
d = varargin{3};
if nargin == 3
    h = [1 1 1];
elseif nargin == 4
    h = varargin{4};
else
    error('Wrong number of inputs to ADAPTFILT');
end;
% im,rad,d,hx,hz

im = scale(im);
hx = h(1);
hy = h(2);
hz = h(3);

[M N O] = size(im);

% z dimension based on the resolution
radz = floor(rad*(hx/hz));

% make threshold image
if O == 1
    th = imfilter(im,fspecial('average',rad),'replicate');
else 
    g = 1/(rad*rad*radz)*ones(rad,rad,radz); 
    th = imfilter(im,g,'replicate');
end
    
% % smooth image
% p = 3;
% if O == 1
%     th = imfilter(im,fspecial('average',p),'replicate');
% else    
%     g = 1/(p*p*p)*ones(p,p,p);
%     th = imfilter(im,g,'replicate');
% end


% multiply thresholds above image mean values
th = th + d;

%
% Create ridges
%

% threshold image at different scales to get ridges
thIm = gt(im,th);

% showall(th,im,thIm)

test = 0;
if test ==1
    show(thIm,2)
    %pause
end;

