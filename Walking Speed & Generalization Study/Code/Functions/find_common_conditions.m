function          [conds_union, conds_intersect, conds_elect] = find_common_conditions(allSubs, cfold)

% cfold=[folder '\'];
% allSubs=dir([cfold '*params.mat']);

%1. Store all the conditions
ns = length(allSubs);
conditions = cell(1,ns);
for s=1:ns
    cname = allSubs(s).name; 
    load([cfold '\' cname]);
    conditions{s} = adaptData.metaData.conditionName;
end

%2. Find common conditions (conditions that everybody did)
condSubset = intersect(removeEmptyCells(conditions{1}), removeEmptyCells(conditions{2}), 'stable');
for s=3:ns
    condSubset = intersect(condSubset, removeEmptyCells(conditions{s}), 'stable');
end
conds_intersect = condSubset;

%3. Find union of conditions (all the different conditions done)
condSuperset = union(removeEmptyCells(conditions{1}), removeEmptyCells(conditions{2}), 'stable');
for s=3:ns
    condSuperset = union(condSuperset, removeEmptyCells(conditions{s}), 'stable');
end
conds_union = condSuperset;

conds_elect = setdiff(conds_union, conds_intersect);

end