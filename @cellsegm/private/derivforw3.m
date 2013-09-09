% [UX UY UZ] = DERIVFORW3(U,HX,HY,HZ) Finding the derivatives of U using a 
% 3 neighbourhood and forward differences.
%
function [ux,uy,uz] = derivforw3(u,hx,hy,hz)

[M N O] = size(u);
ux = (-u([3:end end end],:,:,:) + 4*u([2:end end],:,:,:) - 3*u)/(2*hx);
uy = (-u(:,[3:end end end],:,:) + 4*u(:,[2:end end],:,:) - 3*u)/(2*hy);
if O < 3
    uz = zeros(M,N,O);
else
    uz = (-u(:,:,[3:end end end],:) + 4*u(:,:,[2:end end],:) - 3*u)/(2*hz);
end;

