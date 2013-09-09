% BW = BWAREAOPENRANGE(BW,TH,CONN) Removing small regions in binary 
% image BW below the threshold
% TH, scaled with the range of each connection. CONN is he connectivity
% Returning the new image
%
function [bw] = bwareaopenrange(bw,th,conn)

[faser,L] = bwlabeln(bw,conn);


for j = 1 : L
    regHere = eq(faser,j);
    ind = find(regHere);
    range = bwrange(regHere);
    numPixRel = size(ind,1) / range;
%     show(regHere,3)
%     th
%     pause
    
    if numPixRel < th
        bw(ind) = 0;
    end;
end;

% showall(bw)
