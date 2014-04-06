function [] = imwritemulttif(varargin)
% IMWRITEMULTTIF Print a matrix to a multiple tif file
%
% IMWRITEMULTTIF(NAMETIF,A) Prints the matrix A to a multiple tif file with
% the name nametif. 
%

nametif = varargin{1};
A = varargin{2};
% h = [1 1 1];
% if nargin == 3
%     h = varargin{3};
% end;
dim = size(A);

% make 0 255 and uint8
A = scale(A);
A = 255*scale(A);
A = round(A);
A = uint8(A);

for  i = 1 : dim(3)
   imhere = A(:,:,i);
   
   % write to file
   if i == 1
       imwrite(imhere,nametif,'tif','Compression','none');
   else
       imwrite(imhere,nametif,'tif','WriteMode','Append','Compression','none');
   end;
    
end;
