function [subjects_struct, miMaMatrix, conds] = fill_subject_struct(allSubs, parameters, ns, cfold)
np=length(parameters);
cname=allSubs(1).name;
load([cfold '\' cname]);
nc = length(adaptData.metaData.conditionName);
subjects_struct = cell(ns, nc, np);
conds=cell(1,ns);
for sub=1:ns
    %% Load subject data
    cname=allSubs(sub).name;
    load([cfold '\' cname]);
    conditions = adaptData.metaData.conditionName;
    nc = length(conditions);
    conds{sub} = conditions; 
    for p=1:np
        cp=parameters{p};
        for cond=1:nc
            data = adaptData.getParamInCond(cp, conditions{cond});
            subjects_struct{sub,cond,p} = data;
        end
    end
    
    
end

lengths = cellfun(@length, subjects_struct);
lengths = squeeze(lengths(:,:,1)); %All params have same length
miMaMatrix = zeros(2, nc);
miMaMatrix(1,:) = min(lengths,[],1);
miMaMatrix(2,:) = max(lengths,[],1);

end
