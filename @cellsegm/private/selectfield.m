function [prmout] = selectfield(prm,field)

prmout = [];
f = fieldnames(prm);
n = numel(field);
for i = 1 : numel(f)
    
    
    if numel(f{i}) < n
        continue;
    end;
    % keep the field    
    if isequal(f{i}(1:n),field)
        prmout.(f{i}(n+1:end)) = prm.(f{i});
    end;
end;