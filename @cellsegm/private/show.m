% SHOW Drawing images using colormap gray and imagesc function.
% If 3D matrix then draw the middle plane.
%
function [] = show(varargin)

global inactiv    

if nargin == 1    
    I = varargin{1};
    [M N O] = size(I);
    middle = round(O/2);
    if isempty(inactiv)
        figure;colormap(gray);imagesc(I(:,:,middle));axis image;drawnow
    end;
elseif nargin == 2
    I = varargin{1};
    [M N O] = size(I);    
    number = varargin{2};
    middle = round(O/2);
    if isempty(inactiv)    
        figure(number);colormap(gray);imagesc(I(:,:,middle));axis image;drawnow
    end;
end;%if    
