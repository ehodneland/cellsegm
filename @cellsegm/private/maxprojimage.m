% MAXPROJIMAGE Making the maximum projection
%
%   MAXIMAGE = MAXPROJIMAGE(I) makes the maximum projection along z and
%   returns it in MAXIMAGE
%
function [maximage] = maxprojimage(I)

maximage = max(I,[],4);



