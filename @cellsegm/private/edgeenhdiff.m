% EDGEENHDIFF Edge enhancing diffusion
%
%   EDGEENHDIFF(U,DT,MAXNITER,KAPPA) Performing edge enhancing diffusion 
%   of image U using time step DT, MAXNITER number of iterations, and 
%   conductivity KAPPA
% 
%   Ex: b = edgeenhdiff(im,0.2,100,10);show(im,1);
%
%   Literature:
%   Algorithms for non-linear diffusion, MATLAB in a literate programming
%   style, Section 4.1
%   Rein van den Boomgard
%
%
%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================
%
function [u] = edgeenhdiff( u, dt, maxniter, kappa)

prm.vis = 0;
prm.visit = 10;


% 
% NB Do not filter the original image!!!
%

prm.hx = 1;
prm.hy = 1;
prm.hz = 3;
[M N O] = size(u);

% filter image
if O == 1
    disp('2D edge enhancing diffusion')
    u = filter2d(u,kappa,dt,maxniter,prm);
elseif O >= 2
    disp('3D edge enhancing diffusion')    
    u = filter3d(u,kappa,dt,maxniter,prm);
end;

%-------------------------------------------------

function [u] = filter3d(u,kappa,dt,niter,prm)

hx = prm.hx;
hy = prm.hy;
hz = prm.hz;

filt = gaussian(3,1,1);
for i = 1 : niter
    

    [Rx Ry Rz] = derivcentr3(u,hx,hy,hz);
%     [Rx Ry Rz] = gradient(u,hx,hy,hz);
    Rx = imfilter(Rx,filt,'replicate');
    Ry = imfilter(Ry,filt,'replicate');
    Rz = imfilter(Rz,filt,'replicate');
    
    % for diffusion tensor
    Rw2 = Rx.^2 + Ry.^2 + Rz.^2;
    Rw = sqrt(Rw2);

    % diffusion coefficients
    c3 = exp( - (Rw / kappa).^2 );
    c2 = exp( - (Rw / kappa).^2 );    
%     c3 = 1/5 * c2;
    c1 = 1/5 * c3;

    %
    % To find three orthonormal vectors to the gradient vector
    %

%   http://sci.tech-archive.net/Archive/sci.math/2007-02/msg01210.html
%     Let b=[b1, b2, b3] be your vector.
%     Let B be its size:
%     B = sqrt(b1^2 + b2^2 + b3^2).
%     Set up a "mirror" vector
%     v = [b1+B, b2, b3]
%     or
%     v = [b1-B, b2, b3]
%     whichever makes the first entry bigger in magnitude,
%     and calculate a matrix
%     H = I - 2*v'*v / (v*v')
%     where v' is the transpose of v (writing v in a column), and I is the
%     3-by-3 identity matrix.
%     Then the first row of H will be a unit vector parallel to your b, and the
%     other two rows will be unit vectors orthogonal to b and orthogonal to each
%     other. The programming of this recipe requires no creative thinking, and
%     verification is mechanical, just messy.
    

    % three orthonormal vectors, from the gradient vector
%     eps = 0.001;
    r11 = 1-2.*(Rx+Rw).^2./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r21 = -2.*(Rx+Rw).*Ry./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r31 = -2.*(Rx+Rw).*Rz./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    
    r12 = -2.*(Rx+Rw).*Ry./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r22 = 1-2.*Ry.^2./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r32 = -2.*Ry.*Rz./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    
    r13 = -2.*(Rx+Rw).*Rz./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r23 = -2.*Ry.*Rz./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
    r33 =1-2.*Rz.^2./((Rx+Rw).^2+Ry.^2+Rz.^2 + eps);
 
    % diffusion matrix
    D.d11 = c1.*r11.^2+c2.*r21.^2+c3.*r31.^2;
    D.d12 = r11.*c1.*r12+r21.*c2.*r22+r31.*c3.*r32;
    D.d13 = r11.*c1.*r13+r21.*c2.*r23+r31.*c3.*r33;
    
    D.d21 = D.d12;
    D.d22 = c1.*r12.^2+c2.*r22.^2+c3.*r32.^2;
    D.d23 = r12.*c1.*r13+r22.*c2.*r23+r32.*c3.*r33;
    
    D.d31 = D.d13;
    D.d32 = D.d23;
    D.d33 = c1.*r13.^2+c2.*r23.^2+c3.*r33.^2;

    % update
    u = u + dt * tnldstep3d(u,D,prm);

    % show?
    if prm.vis == 1 && round(niter/prm.visit) == niter/prm.visit
        figure(1);
        subplot(1,2,1);imagesc(imini);colormap(gray);title('Original image')
        subplot(1,2,2);imagesc(u);colormap(gray);title('Edge enhanced image')
    end
end


if prm.vis == 1
    figure(1);
    subplot(1,2,1);imagesc(imini);colormap(gray);title('Original image')
    subplot(1,2,2);imagesc(u);colormap(gray);title('Edge enhanced image')
end

%----------------------------------------------------

function [r] = tnldstep3d(L,D,prm)
               
d11 = D.d11;
d12 = D.d12;
d13 = D.d13;

d21 = D.d21;
d22 = D.d22;
d23 = D.d23;

d31 = D.d31;
d32 = D.d32;
d33 = D.d33;

hx = prm.hx;
hy = prm.hy;
hz = prm.hz;

% function for computing div(D*grad(u))

% NB!!!!
% Use central differences for the mixed terms!!!!

% original scheme
p1 = d11.*(transim(L,1,0,0)-L)/hx;
p2 = d12.*(transim(L,0,1,0)-transim(L,0,-1,0))/(2*hy);
p3 = d13.*(transim(L,0,0,1)-transim(L,0,0,-1))/(2*hz);

p4 = d21.*(transim(L,1,0,0)-transim(L,-1,0,0))/(2*hx);
p5 = d22.*(transim(L,0,1,0)-L)/hy;
p6 = d23.*(transim(L,0,0,1)-transim(L,0,0,-1))/(2*hz);

p7 = d31.*(transim(L,1,0,0)-transim(L,-1,0,0))/(2*hx);
p8 = d32.*(transim(L,0,1,0)-transim(L,0,-1,0))/(2*hy);
p9 = d33.*(transim(L,0,0,1)-L)/hz;

r = ...
    (p1-transim(p1,-1,0,0))/hx + ...
    (transim(p2,1,0,0)-transim(p2,-1,0,0))/(2*hx) + ...
    (transim(p3,1,0,0)-transim(p3,-1,0,0))/(2*hx) + ...
    ...
    (transim(p4,0,1,0)-transim(p4,0,-1,0))/(2*hy) + ...
    (p5-transim(p5,0,-1,0))/hy + ...
    (transim(p6,0,1,0)-transim(p6,0,-1,0))/(2*hy) + ...    
    ...
    (transim(p7,0,0,1)-transim(p7,0,0,-1))/(2*hz) + ...
    (transim(p8,0,0,1)-transim(p8,0,0,-1))/(2*hz) + ...
    (p9-transim(p9,0,0,-1))/hz;


%----------------------------------------------------

function [u] = filter2d(u,kappa,dt,niter,prm)

hx = prm.hx;
hy = prm.hy;
hz = prm.hz;

% iterate
filt = fspecial('gaussian',3,1);
for i = 1 : niter

%     [Rx Ry Rz] = derivcentr3(u,hx,hy,hz);
    [Rx Ry] = gradient(u,hx,hy);
    Rx = imfilter(Rx,filt,'replicate');
    Ry = imfilter(Ry,filt,'replicate');
    
    % for diffusion tensor
    Rw2 = Rx.^2 + Ry.^2;
    Rw = sqrt(Rw2);
    c2 = exp( - (Rw / kappa).^2 );
    c1 = 1/5 * c2;
    a = (c1 .* Rx.^2 + c2 .* Ry.^2) ./ (Rw2+eps);
    b = (c2-c1) .* Rx .* Ry ./ (Rw2+eps);
    c = (c1 .* Ry.^2 + c2 .* Rx.^2) ./ (Rw2+eps);

    % update
    u = u + dt * tnldstep(u, a, b, c );

    % show?
    if prm.vis == 1 && round(niter/prm.visit) == niter/prm.visit
        figure(1);
        subplot(1,2,1);imagesc(imini);colormap(gray);title('Original image')
        subplot(1,2,2);imagesc(u);colormap(gray);title('Edge enhanced image')
    end
end


if prm.vis == 1
    figure(1);
    subplot(1,2,1);imagesc(imini);colormap(gray);title('Original image')
    subplot(1,2,2);imagesc(u);colormap(gray);title('Edge enhanced image')
end