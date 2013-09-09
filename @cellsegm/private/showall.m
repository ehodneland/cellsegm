% SHOWALL Show one plane awt a time for stack.
% SHOWALL Shows all planes in a 3D stack succesively. There can be several
% input images, then they are plottet in different figures.  They have to 
% have the same number of planes. If the last
% argument is 'colorbar', then the images are plotted with the colorbar
% option. 
%
function [] = showall(varargin)

numim = 0;
col = 0;
% plane = 1;
for i = 1 : nargin    
    varhere = varargin{i};
    
    if isequal(varhere,'colorbar')
        col = 1;
        continue;
    end;
    
    f{i} = double(varargin{i});
    numim = numim + 1;
end;
        



F = f{1};
[M N O P] = size(F);

if P > 1
    t = 1;
    while 1        
        for j = 1 : numim
            F = f{j};
            show(F(:,:,:,t),j);            
        end;
        msg = ['Time point ' int2str(t)];
        disp(msg);
        
        t = input('Continue: 0, Time: timepoint, Next timepoint: Enter ');

        if isequal(t,0)
            return
        elseif isempty(t)
            t = t + 1;
        else
            t = t;
        end;                        

        if t > O
            break;
        end;

    end;
else


    niter = 1;
    while 1

        for j = 1 : numim
            F = f{j};

            if col == 0
                show(F(:,:,niter),j)        
            else
                show(F(:,:,niter),j);colorbar
            end;
        end;
        disp(sprintf('Plane %i',niter))

        c = input('Continue: 0, Plane: planenumber, Next plane: Enter ');

        if isequal(c,0)
            return
        elseif isempty(c)
            niter = niter + 1;
        else
            niter = c;
        end;                        

        if niter > O
            break;
        end;
    end;
end;
    