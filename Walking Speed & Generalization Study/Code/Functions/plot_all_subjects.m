function         plot_all_subjects(allSubs,parameters,ns,cfold)
subjects
np=length(parameters);
for p=1:np
    fh(p)=figure;
end

for sub=1:ns
    %% Load subject data
    cname=allSubs(sub).name;
    load([cfold '\' cname]);
    conditions = adaptData.metaData.conditionName;
    nc = length(conditions);
    
    for p=1:np
        cp=parameters{p};
        figure(fh(p))
        
        for cond=1:nc
            data = adaptData.getParamInCond(cp, conditions{cond});
            ax(cond) = subplot(ns, nc , sub + (cond-1))
            hold on
            plot(data,'o');
            linkaxes(ax,'y');
        end
    end
    
    
end

end
