function [prm] = mergeinputpar(prm,prmin)

f = fieldnames(prmin);
for i = 1 : numel(f)
    prm2 = prmin.(f{i});
    if isstruct(prm2)
        % if the field does not exist we do not overwrite and create it
        % with the new field
        if ~isfield(prm,f{i})
            prm.(f{i}) = prm2;
        else
            prm.(f{i}) = mergeinputpar(prm.(f{i}),prm2);
        end;
    else
        prm.(f{i}) = prm2;
    end;    
end;


