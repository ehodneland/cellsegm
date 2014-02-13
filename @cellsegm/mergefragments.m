function [faser] = mergefragments(faser,im,thint,thconv,optlog,optint)
% MERGEFRAGMENTS Merging regiong in segmented image
%
% FASER = MERGEFRAGMENTS (FASER,IM,THINT,THCONV,OPTLOG,OPTINT) Merge regions 
% in piecewise constant image FASER, the true image is IM. 
% NB Only considering regions with values > 0, integers. 
% If borders like in watershed image, let them be zero .
% THINT is the threshold for intensities, THCONV for convexity, and 
% optlog = 1 : merge the strong lines (for CT)
% optlog = 2 : merge the weak lines (for WGA)
% optint = 1 : use the region based merging
% optint = 2 : use the snake based merging
%
% THINT is the relativ threshold for the structure to be merged (<THINT), 
% THCONV is the convexity threshold for merging (<THCONV)
%
% Example for optint = 2: 
% [wat] = mergefragments(wat,im,1.1,1.1,2,2);
%
% Example for optint = 1: 
% [wat] = mergefragments(wat,im,1.5,1.2,2,1);
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



% dim = size(im);

if isempty(intersect(optlog,[1 2]))
    error('Wrong option for OPTLOG')
end;

if isempty(intersect(optint,[1 2]))
    error('Wrong option for OPTINT')
end;

if optlog == 1
    sign = '>';
elseif optlog == 2
    sign = '<';
end

% thresholds
th.int = thint;
th.conv = thconv;
th.vol = Inf;

msg = ['This is ' upper(mfilename) ' using parameters'];
disp(msg);
printstructscreen(th);

% no border set to 0
prm.isborder = 0;

disp(sprintf('  Merging two regions if the value is %s than the threshold',sign))
if optint == 1
    disp(sprintf('  Comparing the intensities of the watershed lines to the mean of the regions around'))
elseif optint == 2
    disp(sprintf('  Comparing the intensities of the watershed lines to a snake around'))
end;

pairs.checked = [];
pairs.meanint =  [];
pairs.ind = cell(0,1);

prm.isborder = 1;
if isempty(find(faser == 0,1))
    prm.isborder = 0;
end;

% nfaserold = Inf;
% nfaser = 0;
% while ne(nfaserold,nfaser)
    
% get info for all regions
% the regions to loop over
valfaser = unique(faser(faser > 0));
L = length(valfaser);

for i = 1 : L    
    % get parameters
    pairs = getprm(valfaser(i),im,faser,pairs,optint,prm);       
end;

disp('  The regions checked')
[pairs.checked pairs.meanint]


%
% Start to merge
%
niter = 0;
while ~isempty(pairs.meanint)
    niter = niter + 1;

    %
    % Check for intensity
    %
    
    % always pick the maximum/minimum intense border, here we find it
    if optlog == 1
        [trial.int,indext] = max(pairs.meanint);
    elseif optlog == 2
        [trial.int,indext] = min(pairs.meanint);
    else
        error('Wrong option')
    end;

    % the regions on each side
    labelhere = [pairs.checked(indext,1) pairs.checked(indext,2)];


    % the lowest and highes wataershed region to choose one of them to
    % survive
    minreg = min(labelhere);
    maxreg = max(labelhere);

    %
    % Check for convexity
    %

    % the first region
    reghere1 = eq(faser,minreg);    
    area1 = length(find(reghere1));
    convarea1 = calcconvarea(reghere1);

    % the other region
    reghere2 = eq(faser,maxreg);   
    area2 = length(find(reghere2));    
    convarea2 = calcconvarea(reghere2);

    % both regions combined
    regboth = gt(reghere1 + reghere2,0);
    load ball2;se = getball(ball,2,1);
    regboth = imclose(regboth,se);
    areaboth = length(find(regboth));

    % convex area combined
    convareaboth = calcconvarea(regboth);

    % how it works:
    % ratioconv = 1 if no change in convexity after merging
    % ratioconv > 1 if decrease in convexity after merging
    % ratioconv < 1 if increase in convexity after merging
    % merge for ratioconv < prm.conv
    a = area1;
    b = area2;
    A = convarea1;
    B = convarea2;
    ab = areaboth;
    AB = convareaboth;
    convbef = (a+b)/(A+B);
    convaft = ab/AB;
    trial.conv = convbef/convaft;

    %
    % Ratio volume
    %
    trial.vol = max(area1/area2,area2/area1);

%     ratioconv = convareaalone / convareaboth;

    msg = ['*********************************'];    
    disp(msg);
    msg = ['  Regions: ' int2str(minreg) ' and ' int2str(maxreg)];
    disp(msg);
    msg = ['  Ratio int to merge: ',num2str(trial.int) ' ' sign ' ' num2str(th.int)];
    disp(msg);
    msg = ['  Ratio convexity to merge: ',num2str(trial.conv) ' ' sign ' ' num2str(th.conv)];            
    disp(msg);
    msg = ['  Ratio volume to merge: ' num2str(trial.vol)  ' ' sign ' '  num2str(th.vol)];        
    disp(msg);
    merge = 1;

    % shall we continue to loop?
    if eval(['trial.int' sign 'th.int'])
        msg = ['  Merging due to intensity'];
    else
        msg = ['  Not merging due to intensity, terminating merging process'];
        break;
    end;
    disp(msg);

    %
    % Test convexity
    %

    if eval(['trial.conv' sign 'th.conv'])
        msg = ['  Merging due to convexity'];        
    else
        msg = ['  Not merging due to convexity'];            
        merge = 0;
    end;
    disp(msg);

    %
    % Test volume
    %
    if eval(['trial.vol' sign 'th.vol'])
        msg = ['Merging due to difference in volume'];        
    else
        msg = ['  Not merging due to difference in volume'];            
        merge = 0;
    end;
    disp(msg);

    if merge == 0
        % this border shall never be removed, take away from the queue but
        % we want to continue for next loop
        pairs = rembordinfo(pairs,indext);        
        continue;
    end;    

    %
    % Remove the border
    %      
    msg = ['  Merging ' int2str(minreg) ' and ' int2str(maxreg)];
    disp(msg);

    % this border piece
    indhere = pairs.ind{indext};

    prm.nitervis = Inf;
    if niter > prm.nitervis
        show(faser,3)
    end;

    % give the border the value of the lowest watershed label
    faser(indhere) = minreg;


    if niter > prm.nitervis
        show(faser,4)
    end;

    % THIS border is gone, so also cut the values belonging to it
    pairs = rembordinfo(pairs,indext);        

    % also remove from the queue the borders AROUND each of the regions, facing outward
    ind1 = find(pairs.checked(:,1) == minreg | pairs.checked(:,1) == maxreg);
    ind2 = find(pairs.checked(:,2) == minreg | pairs.checked(:,2) == maxreg);        
    ind = unique([ind1 ; ind2]);
    pairs = rembordinfo(pairs,ind);

    %
    % Fix HIGHEST watershed region
    %


    % give the highest the value of the lowest watershed label
    faser(faser == maxreg) = minreg;

%         % close them to remove small outliers in the border
%         reg = faser == minreg;
%         % NB must be O here!
%         load ball2;se = getball(ball,2,O);
%         reg = imclose(reg,se);
%         faser(reg) = minreg;

    %
    % Add values for newly created region
    %

    % get parameters again for the LOWEST watershed region which is a new
    % region!!
    pairs = getprm(minreg,im,faser,pairs,optint,prm);       

    if niter > prm.nitervis
        minreg
        pairs
        show(reghere1,1)
        show(reghere2,2)
        show(faser,5)
        a = zeros(size(faser));
        a(indhere) = 1;
        show(a,6)
        pause
    end;
end;

msg = ['Finally checked'];
disp(msg);
[pairs.checked pairs.meanint]

% relabel
val = unique(faser(faser > 0));
faserout = zeros(size(faser));
for i = 1 : length(val)
    faserout(faser == val(i)) = i;
end;
faser = faserout;

%----------------------------------------------------------



function [pairs] = rembordinfo(pairs,ind)

label = 1 : length(pairs.meanint);
label = setdiff(label,ind);
pairs.ind = cutcell(pairs.ind,label);
pairs.meanint(ind) = [];
pairs.checked(ind,:) = [];                

%----------------------------------------------------------------------        
        
function [pairs] = getprm(valhere,im,faser,pairs,optlog,prm)

dim = size(im);
if numel(dim) == 2
    dim = [dim 1];
end;
% this region
segm.reghere = eq(faser,valhere);

load ball1;se1 = getball(ball,1,dim(3));
load ball2;se2 = getball(ball,2,dim(3));
if prm.isborder
    segm.dilreghere = imdilatefast(segm.reghere,se2);
else
    segm.dilreghere = imdilatefast(segm.reghere,se1);
end;

% get the neighbours of the new combined region
neigh = getneigh(faser,valhere);

% number of neighbours
nneigh = length(neigh);

% must initiliaze since we may not loop at all, and need something to
% return
count = length(pairs.ind);
for i = 1 : nneigh

    % this neighbour
    neighhere = neigh(i);    
    
    % check for earlier intersection
    chkbef = ~isempty(intersect(pairs.checked,[neighhere valhere]));

    % have checked before so continue
    if chkbef
        msg = ['Checked before ' int2str(valhere) ' and ' int2str(neighhere) ', continue'];
        disp(msg);
        continue;
    end;
    
    % this neighbour    
    segm.regneigh = eq(faser,neighhere);

    % both
    segm.both = zeros(size(im));
    segm.both(segm.reghere) = 1;
    segm.both(segm.regneigh) = 1;
    
    if prm.isborder
        segm.dilregneigh = imdilatefast(segm.regneigh,se2);
    
        % the common border within the watershed image
        se = ones(3,3,3);
        % NB need to inside the dilation of both also!! Imclose itself gives
        % far away small regions
    
        segm.border = segm.dilreghere .* segm.dilregneigh .* (imclose(segm.both,se) - segm.both);
    else
         segm.dilregneigh = imdilatefast(segm.regneigh,se1);
        segm.border = segm.dilreghere .* segm.dilregneigh .* (segm.reghere  + segm.regneigh);
    end;
    
    % save the indices for the removing step
    indhere = find(segm.border);
    
    % only take if long enough
    [range,minz,maz] = bwrange(segm.border);
    relnumpixborder = length(indhere)/range;    
    
    if isempty(find(segm.border,1))
        msg = ['Empty border, continuing'];
        disp(msg);
        continue;
    end;
    
    % local indexing
    count = count + 1;

    % save the indices to return
    pairs.ind{count,1} = indhere;    

    if optlog == 1
        
        %
        % Mean intensity of border compared to the regions that are divided
        %

        % relative intensity
        b = mean(im(segm.border == 1));
        valreghere = mean(im(segm.reghere));
        valregneigh = mean(im(segm.regneigh));
        a = valreghere;
        c = valregneigh;
        
              
        pairs.meanint(count,1) = max(b/a,b/c);

    elseif optlog == 2
        %
        % Make the dilation of the border, a snake along
        %

        se21 = strel('disk',5);
        se22 = strel('disk',2);    
        segm.snake = imdilatefast(segm.border,se21) - imdilatefast(segm.border,se22);

        % only consider intensities on the current regions!
        segm.snake = segm.snake .* segm.both;    
        
        % relative intensity of border        
        pairs.meanint(count,1) = mean(im(segm.border == 1)) / mean(im(segm.snake == 1));   

    end;
    vis = 0;
    if vis == 1
    save test faser
    pairs.checked
    [valhere neighhere]
    pairs.meanint
    pairs.meanint(count)
    count
    pairs.checked
    show(im,1)
    show(segm.border,2)
    show(segm.snake,3)
    show(segm.both,4)
    show(segm.reghere,5)
    show(segm.regneigh,6)
    show(faser,7)
    pause
    end
    % the checked
    pairs.checked(count,:) = [valhere neighhere];

end;


% pairs.meanint
% show(faser,9)
% show(reghere,10)
% neigh
% pause
%---------------------------------------------------------------        
        

function [ind] = cutcell(ind,I)

indout = [];
for i = 1 : length(I)
    indout{i} = ind{I(i)};    
end;
ind = indout';



% -------------------------------------------------------

function [neigh] = getneigh(faser,i)
    
% this region
reghere = eq(faser,i);

load ball2;se = getball(ball,2,1);
dilreghere = imdilate(reghere,se);

neighall = faser(dilreghere);
neighall(neighall == i) = [];
neighall(neighall == 0) = [];
neigh = unique(neighall);
numneigh = length(neigh);
for j = 1 : numneigh
    neighhere = neigh(j);
    numborder = length(neighall == j);
%     numborder
    if numborder < 10
        neighall(neighall == neighhere) = [];
    end;
end;
neigh = unique(neighall);

%----------------------------------------------------------

function [convarea] = calcconvarea(BW)

dim = size(BW);
if numel(dim) == 2
    dim = [dim 1];
end;
convarea = 0;
for j = 1 : dim(3)
    if isempty(find(BW(:,:,j),1))
        continue;
    end;
%     try
    props = regionprops(double(BW(:,:,j)),'ConvexArea');
    convarea = convarea + props.ConvexArea;        
%     catch
%         continue;
%     end;
end;

% -------------------------------------------------------------

        
