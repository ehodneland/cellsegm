function [u] = transim(varargin)
% TRANSIM Translate image for use in derivatives
%
%   U = TRANSIM(U,DI,DJ,DK) Translates image for use in derivatives, returns
%   the translated image U. DI, DJ, DK are integer valued translations.
%
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
u = varargin{1};
di = varargin{2};
dj = varargin{3};
dk = varargin{4};
if nargin == 5
    dl = varargin{5};
end;
[M N O P] = size(u);

di = sign(di) * min(M,abs(di));
dj = sign(dj) * min(N,abs(dj));
dk = sign(dk) * min(O,abs(dk));
if nargin == 5
    dl = sign(dl) * min(P,abs(dl));
end;

if di > 0
    indi = [(di+1:M) M*ones(1,di)];
elseif di == 0
    indi = 1:M;
elseif di < 0
    indi = [ones(1,-di) 1:M+di];
end;

if dj > 0
    indj = [(dj+1:N) N*ones(1,dj)];
elseif dj == 0
    indj = 1:N;
elseif dj < 0
    indj = [ones(1,-dj) 1:N+dj];
end;

if dk > 0
    indk = [(dk+1:O) O*ones(1,dk)];
elseif dk == 0
    indk = 1:O;
elseif dk < 0
    indk = [ones(1,-dk) 1:O+dk];
end;


% the image to return
if nargin == 5

    if dl > 0
        indl = [(dl+1:P) P*ones(1,dl)];
    elseif dl == 0
        indl = 1:P;
    elseif dl < 0
        indl = [ones(1,-dl) 1:P+dl];
    end;
    u = u(indi,indj,indk,indl);
    
else
    u = u(indi,indj,indk,:);
end;