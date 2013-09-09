% FILLHOLES FILL HOLES WITHIN A SPECIFIED RANGE
% FILLHOLES(BW,LOWTH,HIGHTH,CONN) Fill holes in BW being within 
% the range of LOWTH and HIGHTH.
%
function [bw] = fillholes(bw,lowth,highth,conn)

[M N O] = size(bw);

[faser,L] = bwlabeln(bw);

holes = zeros(size(bw));
holesall = holes;
for i = 1 : L
    
    % this region
    reghere = eq(faser,i);
    
    % fill this region
    for j = 1 : O        
        holes(:,:,j) = imfill((reghere(:,:,j)),conn,'holes') - reghere(:,:,j);        
    end;
    holes = logical(holes);
    holesall(holes == 1) = 1;
end

holesall = bwareaopen(holesall,round(lowth));
holesall = holesall - bwareaopen(holesall,round(highth));

% fill
bw(holesall == 1) = 1;