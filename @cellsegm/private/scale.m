% SCALE(IM) Scaling.
% SCALE(IM) scales the image IM between 0 and 1.
% SCALE(IM,LOW,HIGH) scales the image IM between LOW and HIGH.
%
function [im,minscale,maxscale] = scale(varargin)

im = varargin{1};
if nargin == 1
    low = 0;
    high = 1;
else
    low = varargin{2};
    high = varargin{3};
end;

if isempty(find(im,1))
    disp('Empty matrix, no scaling.')
    return;
end;

% scale to 0 and 1
minscale = min(im(:));
im = im - minscale;

maxscale = max(im(:));
if ne(maxscale,0) || ~isempty(find(im,1))  
    im = im / maxscale;
end;

% scale to high and low
im = im * (high-low);
im = im + low;

