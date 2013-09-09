% [RANGE,MINZ,MAXZ] = FINDRANGE(CONNHERE) Finding range of (one) connection CONNHERE, 
% number of
% z-planes. Returning the range of the connection in CONNHERE.
%
function [range,minZ,maxZ] = bwrange(connHere)

[k(:,1) k(:,2) k(:,3)] = ind2sub(size(connHere),find(connHere));    
range = max(k(:,3)) - min(k(:,3)) + 1;
minZ = min(k(:,3));
maxZ = max(k(:,3));