% [UXX UYY UZZ UXY UYZ UXZ] = DERIV2CENTR5(U,HX,HY,HZ) Finding the second 
% derivatives of U using a 5 neighbourhood and central differences.
%
% Ex : [uxx,uyy,uzz,uxy,uyz,uxz] = deriv2centr5(u,1,1,1)
%
function [uxx,uyy,uzz,uxy,uyz,uxz] = deriv2centr5(u,hx,hy,hz)



[M N O] = size(u);

% central
uxx = (-u([3:end end end],:,:,:) + 16*u([2:end end],:,:,:) - 30*u + ...
      16*u([1 1:end-1],:,:,:) - u([1 1 1:end-2],:,:,:))/(12*hx^2);
uyy = (-u(:,[3:end end end],:,:) + 16*u(:,[2:end end],:,:) - 30*u + ...
      16*u(:,[1 1:end-1],:,:) - u(:,[1 1 1:end-2],:,:))/(12*hy^2);
  
if O < 3
    uzz = zeros(M,N,O);
else
    uzz = (-u(:,:,[3:end end end],:) + 16*u(:,:,[2:end end],:) - 30*u + ...
        16*u(:,:,[1 1:end-1],:) - u(:,:,[1 1 1:end-2],:))/(12*hz^2);
end;

% showall(u,uzz)

%
% mixed
%
% [ux,uy,uz] = derivcentr5(u,hx,hy,hz);
[ux,uy,uz] = derivforw3(u,hx,hy,hz);

% use one forward step for uz
% uz = (u(:,:,[2:end end]) - u)/hz;

% backward y
uxy = (3*ux - 4*ux(:,[1 1:end-1],:,:) + ux(:,[1 1 1:end-2],:,:))/(2*hy);
% backward x
uyx = (3*uy - 4*uy([1 1:end-1],:,:,:) + uy([1 1 1:end-2],:,:,:))/(2*hx);
uxy = (uxy + uyx)/2;

% backward z
if O < 3
    uyz = zeros(M,N,O);
else
%     uyz = (uy - uy(:,:,[1 1:end-1]))/hz;
%     uzy = (uz - uz(:,[1 1:end-1],:))/hy;
    uyz = (3*uy - 4*uy(:,:,[1 1:end-1],:) + uy(:,:,[1 1 1:end-2],:))/(2*hz);
    uzy = (3*uz - 4*uz(:,[1 1:end-1],:,:) + uz(:,[1 1 1:end-2],:,:))/(2*hy);
    uyz = (uyz + uzy)/2;
end;
% backward y
% uzy = (3*uz - 4*uz(:,[1 1:end-1],:) + uz(:,[1 1 1:end-2],:))/(2*hy);


% backward z
if O < 3
    uxz = zeros(M,N,O);    
else
%     uxz = (ux - ux(:,:,[1 1:end-1]))/hz;
%     uzx = (uz - uz([1 1:end-1],:,:))/hx;
    uxz = (3*ux - 4*ux(:,:,[1 1:end-1],:) + ux(:,:,[1 1 1:end-2],:))/(2*hz);
    uzx = (3*uz - 4*uz([1 1:end-1],:,:,:) + uz([1 1 1:end-2],:,:,:))/(2*hz);
    uxz = (uxz + uzx)/2;
end;
% backward x
% uzx = (3*uz - 4*uz([1 1:end-1],:,:) + uz([1 1 1:end-2],:,:))/(2*hx);


% I took this away 20100909
% % fix borders
% p = 4;
% uxx = multborder(uxx,1,p);
% uyy = multborder(uyy,1,p);
% uxy = multborder(uxy,1,p);
% if O >= 3
%     uxz = multborder(uxz,1,p);
%     uyz = multborder(uyz,1,p);
%     uzz = multborder(uzz,1,p);
% end;
