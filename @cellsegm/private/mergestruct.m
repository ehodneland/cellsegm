% MERGESTRUCT Merging two structs
% PRMDEFAULT = MERGESTRUCT(PRMDEFEAULT,PRMIN) merging the two structs PRMDEFAULT and
% PRMIN where the fields in PRMDEFAULT are overwritten by those in PRMIN.
%
function [prmdefault] = mergestruct(prmdefault,prmin)
    
    
    if isempty(prmin)        
        return;
    end;
        
    % names of parameter fields
    varnames = fieldnames(prmin);
    for i = 1 : size(varnames,1)
        varnamehere = varnames{i};
        valhere = eval(['prmin.' varnamehere]);
        % we keep the given value!!
        if isequal(valhere,'default')
            continue;
        end;
        evalc(['prmdefault' '.' varnamehere '=' 'prmin' '.' varnamehere]);
    end;
