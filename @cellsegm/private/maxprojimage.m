function [maxImage] = maxprojimage(I)

[M N O] = size(I);
maxImage = I(:,:,1);
depImage = ones(M,N);

for i = 2 : O
    maxImageOld = maxImage;
    maxImage = max(maxImage,I(:,:,i));
    depImage(find(I(:,:,i) > maxImageOld)) = i;
end;%for

% show(maxImage,2)
% show(depImage,3)



