% [UX UY UZ] = DERIVCENTR3(U,HX,HY,HZ) Finding the derivatives of U using a 
% 3 neighbourhood and central differences.
%
function [ux,uy,uz] = derivcentr3(u,hx,hy,hz)

dim = size(u);
ndim = numel(dim);
ux = (u([2:end end],:,:,:) - u([1 1:end-1],:,:,:))/(2*hx);
uy = (u(:,[2:end end],:,:) - u(:,[1 1:end-1],:,:))/(2*hy);
if ndim >= 3
    uz = (u(:,:,[2:end end],:) - u(:,:,[1 1:end-1],:))/(2*hz);
else
    uz = zeros(dim);
end;

% 
% % forwaard differences at the negatve border
% ux(1,:,:) = (u(2,:,:) - u(:,:,1))/hx;
% % backward differences at the posiotive border
% ux(end,:,:) = (u(end,:,:) - u(end-1,:,:))/hx;
% 
% % forwaard differences at the negatve border
% uy(:,1,:) = (u(:,2,:) - u(:,1,:))/hy;
% % backward differences at the positive border
% uz(:,end,:) = (u(:,end,:) - u(:,end-1,:))/hy;
% 
% % forwaard differences at the negatve border
% uz(:,:,1) = (u(:,:,2) - u(:,:,2))/hz;
% % backward differences at the positive border
% uz(:,:,end) = (u(:,:,end) - u(:,:,end-1))/hz;

