% function [ball2] = getball(ball,r,O)
% Returning the correct size of the structuring element 
% based on the radius of it and the number of planes.
% BALL  : structural element
% R     : radius of structural element
% O     : number of planes
% radius of BALL must be equal to R
% 
function [ball2] = getball(ball,r,O)

if O > 2*r+1
    ball2 = ball;
else
    fra = r-floor(O/2)+eq(round(O/2),O/2) +1;
    til = r+floor(O/2)+1;
    ball2 = ball(:,:,fra:til);
end;%
% fra = r-floor(O/2)+eq(round(O/2),O/2) + eq(O,1)
% til = r + floor(O/2)+eq(O,1)
% lt(fra,1) + fra*ge(fra,1) : til*le(fra,2*r+1) + (2*r+1)*gt(til,2*r+1)
% ball2 = ball(:,:,lt(fra,1) + fra*ge(fra,1) : til*le(fra,2*r+1) + (2*r+1)*gt(til,2*r+1));




