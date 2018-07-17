function conditions = display_all_conditions(folder,display)
cfold=[folder '\'];
allSubs=dir([cfold '*params.mat']);

ns=length(allSubs);
conditions = cell(1,ns);

for s=1:ns
    cname = allSubs(s).name; 
    load([cfold '\' cname]);
    conditions{s} = adaptData.metaData.conditionName;
    
    if display
       disp(cname)
       conditions{s} 
    end
    
end

end