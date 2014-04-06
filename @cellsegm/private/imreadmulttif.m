function [im] = imreadmulttif(imload,gray)
% IMREADMULTTIF Reading multiple tif file
%
% IM = IMREADMULTTIF(IMLOAD,GRAY) reads image multiple tif file IMLOAD and
% converts it to grayscale if GRAY = 1.
% 
% Returning the 3D image in IM
%
info = imfinfo(imload);
nim = numel(info);
imhere = imread(imload,1);          
a = size(imhere,3);
if gray == 1
    dim = [info(1).Height info(1).Width nim 1]; 
else
    dim = [info(1).Height info(1).Width nim a]; 
end;
im = zeros(dim);
for i = 1 : nim
    imhere = imread(imload,i);          
    
    % if loaded tif is rgb image that we want to turn
    % into grayscale
    if gray == 1 && dim(4) == 3
        
        imhere = rgb2gray(imhere);  
        im(:,:,i) = double(imhere);    
        
    % if we want to preserve the tif format as it is
    else            
        for j = 1 : dim(4)
            im(:,:,i,j) = double(imhere(:,:,j));
        end;
    end;
    
end;
