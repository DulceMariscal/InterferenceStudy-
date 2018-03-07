%With S009
% pathToData= 'C:\Users\dum5\OneDrive\_Shared Drive - Interference Project - Alessandro - Dulce_\Params Files\ForceParams';
clear all
close all
clc

pathToData= 'C:\Users\dum5\OneDrive\_Shared Drive - Interference Project - Alessandro - Dulce_\Params Files\subjectsData';
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0]];
%For step length asymmetry
% load('SlaAvgIdv.mat') 
% % load('stepLengthAsymAllData.mat')
% load('StridesToRemove.mat')
% param={'stepLengthAsym'};
% slaI=stepLengthAsymINT;
% slaS=stepLengthAsymSAV;


%for ExtAdaptaiton= SLA-SV
% load('ExtAdaptNorm2AllData.mat') 
% load('StridesToRemoveExtAdaptNorm2.mat')
% param={'ExtAdaptNorm2'};
% slaI=ExtAdaptINT;
% slaS=ExtAdaptSAV;

%No hip 
% load('ExtAdaptPNormAllData.mat') 
% load('StridesToRemoveExtAdaptPNorm.mat')
% param={'ExtAdaptPNorm'};
% slaI=ExtAdaptINT;
% slaS=ExtAdaptSAV;

% Shifted Data SLA 
 load('stepLengthAsymShifted.mat')
load('StridesToRemovestepLengthAsymShifted.mat')
param={'stepLengthAsymShifted'};

%Shifted Data Ext Adapt Norm2 
% load('ExtAdaptNorm2Shifted.mat')
% load('StridesToRemoveExtAdaptNorm2.mat')
% param={'ExtAdaptNorm2Shifted'};

%Shifted Data Ext Adapt PNorm
% load('ExtAdaptPNormShifted.mat')
% load('StridesToRemoveExtAdaptPNorm.mat')
% param={'ExtAdaptPNormShifted'};


groups = 2;
charGroups = {'S','I'};
indSubs = {setdiff(1:9, [2 5]), setdiff(1:8, [1 7])};
%  indSubs = {setdiff(1:9, []), setdiff(1:8, [])};

subs = {length(indSubs{1}) , length(indSubs{2})};
% subs ={9,8};

incrORdecr={'Increasing'};

cond_Inter=[2,5,6];
cond_Sav=[2,4,5];

NumbStridesToMidPert=nan(9,4);
MinPointAdaptCurve_Inter=nan(9,1);
MinPointAdaptCurve_Sav=nan(9,1);
Tau_INT=nan(8,2);
Tau_SAV=nan(9,2);
STRIDES_1_INT=nan(8,2);
STRIDES_5_INT=nan(8,2);
STRIDES_SS_INT=nan(8,2);
STRIDES_1_SAV=nan(9,2);
STRIDES_5_SAV=nan(9,2);
STRIDES_SS_SAV=nan(9,2);
Washout_5_SAV=nan(9,1);
Washout_5_INT=nan(8,1);
Washout_ALL_SAV=nan(9,1);
Washout_ALL_INT=nan(8,1);


for s=1:subs{2}
    sub=indSubs{2}(s);
%     sub=s;
    data_A1_Inter=slaI{sub,cond_Inter(1)};
    [nos_to_mid_pert_InterA1, dataFit_InterA1, expFitPars_InterA1] = find_nstrides_to_mid_pert(data_A1_Inter(minNOSToRemove(sub,1):end,1),incrORdecr);
    MinPointAdaptCurve_Inter(sub,1)=nanmin(data_A1_Inter(minNOSToRemove(sub,1):end,1));
    Tau_INT(sub,1)=expFitPars_InterA1(3);
    STRIDES_1_INT(sub,1)=data_A1_Inter(minNOSToRemove(sub,1):minNOSToRemove(sub,1),1);
    STRIDES_5_INT(sub,1)=nanmean(data_A1_Inter(minNOSToRemove(sub,1):minNOSToRemove(sub,1)+4,1));
    STRIDES_SS_INT(sub,1)=nanmean(data_A1_Inter(end-40:end-5,1));
    
    data_A2_Inter=slaI{sub,cond_Inter(2)};
    [nos_to_mid_pert_A2_Inter, dataFit_A2_Inter, expFitPars_A2_Inter] = find_nstrides_to_mid_pert(data_A2_Inter(minNOSToRemove(sub,1):end,1),incrORdecr);
    Tau_INT(sub,2)=expFitPars_A2_Inter(3);
    STRIDES_1_INT(sub,2)=data_A2_Inter(minNOSToRemove(sub,1):minNOSToRemove(sub,1),1);
    STRIDES_5_INT(sub,2)=nanmean(data_A2_Inter(minNOSToRemove(sub,1):minNOSToRemove(sub,1)+4,1));        
    STRIDES_SS_INT(sub,2)=nanmean(data_A2_Inter(end-40:end-5,1));
    
    data_A3_Inter=slaI{sub,cond_Inter(3)};
    if sub==8 && strcmp(param,'stepLengthAsym') || sub==8 && strcmp(param,'ExtAdaptNorm2Shifted') || sub==8 && strcmp(param,'ExtAdaptNorm2')||...
         sub==8 && strcmp(param,'ExtAdaptPNorm')|| sub==8 && strcmp(param,'ExtAdaptPNormShifted') || sub==8 && strcmp(param,'stepLengthAsymShifted')   
        Washout_5_INT(sub,1)=nanmean(data_A3_Inter(9:14,1));
    else
        Washout_5_INT(sub,1)=nanmean(data_A3_Inter(1:5,1));
    end
    
%     Washout_5_INT(sub,1)=nanmean(data_A3_Inter(1:5,1));
    
    Washout_ALL_INT(sub,1)=nanmean(data_A3_Inter(1:end-5,1));
    
    
    NumbStridesToMidPert(sub,1:2)=[nos_to_mid_pert_InterA1 nos_to_mid_pert_A2_Inter ];
end

Inter_Epoachs=[STRIDES_5_INT(:,1) STRIDES_SS_INT(:,1) STRIDES_5_INT(:,2) STRIDES_SS_INT(:,2) Washout_5_INT Washout_ALL_INT];

for s=1:subs{1}
    sub=indSubs{1}(s);
%     sub=s;
    data_A1_Sav=slaS{sub,cond_Sav(1)};
    [nos_to_mid_pert_A1_Sa, dataFit_A1_Sa, expFitPars_A1_Sa] = find_nstrides_to_mid_pert(data_A1_Sav(minNOSToRemove(sub,2):end,1),incrORdecr);
    MinPointAdaptCurve_Sav(sub,1)=nanmin(data_A1_Sav(minNOSToRemove(sub,2):end,1));
    Tau_SAV(sub,1)=expFitPars_A1_Sa(3);
    STRIDES_1_SAV(sub,1)=data_A1_Sav(minNOSToRemove(sub,2):minNOSToRemove(sub,2),1);
    STRIDES_5_SAV(sub,1)=nanmean(data_A1_Sav(minNOSToRemove(sub,2):minNOSToRemove(sub,2)+4,1));
    STRIDES_SS_SAV(sub,1)=nanmean(data_A1_Sav(end-40:end-5,1));
    
    
    data_A2_Sav=slaS{sub,cond_Sav(2)};
    [nos_to_mid_pert_A2_Sa, dataFit_A2_Sa, expFitPars_A2_Sa] = find_nstrides_to_mid_pert(data_A2_Sav(minNOSToRemove(sub,2):end,1),incrORdecr);
    Tau_SAV(sub,2)=expFitPars_A2_Sa(3);
    STRIDES_1_SAV(sub,2)=data_A2_Sav(minNOSToRemove(sub,2):minNOSToRemove(sub,2),1);
    STRIDES_5_SAV(sub,2)=nanmean(data_A2_Sav(minNOSToRemove(sub,2):minNOSToRemove(sub,2)+4,1));
    STRIDES_SS_SAV(sub,2)=nanmean(data_A2_Sav(end-40:end-5,1));
    
    
    data_A3_Sav=slaS{sub,cond_Sav(3)};
    Washout_5_SAV(sub,1)=nanmean(data_A3_Sav(1:5,1));
    Washout_ALL_SAV(sub,1)=nanmean(data_A3_Sav(1:end-5,1));
    
    NumbStridesToMidPert(sub,3:4)=[nos_to_mid_pert_A1_Sa nos_to_mid_pert_A2_Sa];
end

Sav_Epoachs=[STRIDES_5_SAV(:,1) STRIDES_SS_SAV(:,1) STRIDES_5_SAV(:,2) STRIDES_SS_SAV(:,2) Washout_5_SAV Washout_ALL_SAV];


Inter_Epoachs(any(isnan(Inter_Epoachs), 2), :) = [];
Sav_Epoachs(any(isnan(Sav_Epoachs), 2), :) = [];

MinPointAdaptation=nanmin([MinPointAdaptCurve_Inter; MinPointAdaptCurve_Sav]);




save([ param{1} 'MinPointAdapt.mat'], 'MinPointAdaptation')
save([ param{1} 'NumbStridesToMidPert.mat'], 'NumbStridesToMidPert')
save([ param{1} 'First&SS.mat'], 'STRIDES_1_INT','STRIDES_SS_INT','STRIDES_1_SAV','STRIDES_SS_SAV')
save([param{1} 'Epochs.mat'],'Sav_Epoachs','Inter_Epoachs')

cond=5;
Sav_Results=[];
Inter_Results=[];
for c=1:cond
Inter_Results=[Inter_Results; Inter_Epoachs(:,c) c*ones(size(Inter_Epoachs,1),1) ones(size(Inter_Epoachs,1),1)];
Sav_Results=[Sav_Results; Sav_Epoachs(:,c) c*ones(size(Sav_Epoachs,1),1) 2*ones(size(Sav_Epoachs,1),1)];
end 
stata_5_=[Inter_Results;Sav_Results];

Sav_Results=[];
Inter_Results=[];
for c=[1 2 3 4 6]
Inter_Results=[Inter_Results; Inter_Epoachs(:,c) c*ones(size(Inter_Epoachs,1),1) ones(size(Inter_Epoachs,1),1)];
Sav_Results=[Sav_Results; Sav_Epoachs(:,c) c*ones(size(Sav_Epoachs,1),1) 2*ones(size(Sav_Epoachs,1),1)];
end 
stata_all_=[Inter_Results;Sav_Results];

save(['Stata' param{1} ],'stata_all_','stata_5_')
% 






%%
%Tau Analysis 

PV_INT=(Tau_INT(:,2)-Tau_INT(:,1))./Tau_INT(:,1);
PV_SAV=(Tau_SAV(:,2)-Tau_SAV(:,1))./Tau_SAV(:,1);
PV_INT(any(isnan(PV_INT), 2), :) = [];
PV_SAV(any(isnan(PV_SAV), 2), :) = [];
[h,p,ci,stats] = ttest2(PV_INT,PV_SAV);


figure
hold on

bar(1,mean(PV_INT),'facecolor','r')
bar(2,mean(PV_SAV),'facecolor','b')
for sub=1:length(PV_INT)
     plot(1,PV_INT(sub,1),'*','MarkerFaceColor',colorOrder(sub,:))
 end
 for sub=1:length(PV_SAV)
     plot(2,PV_SAV(sub,1),'*','MarkerFaceColor',colorOrder(sub,:))
 end
 
errorbar(1,mean(PV_INT),std(PV_INT)/sqrt(length(PV_INT)),'.','LineWidth',2,'Color','k')
errorbar(2,mean(PV_SAV),std(PV_SAV)/sqrt(length(PV_SAV)),'.','LineWidth',2,'Color','k')
legend('Interference','Savings')

ylabel('PV=tau_{A2}-tau_{A1}/tau_{A1}')
title([param; 't-test p-value=' num2str(p)])



% set(gca,'Xtick',[1.5 4.5],'XTickLabel',{'A1' 'A2'},'fontSize',12)

%% 
%Washout without outliers

Inter_Epoachs(any(isnan(Inter_Epoachs), 2), :) = [];
Sav_Epoachs(any(isnan(Sav_Epoachs), 2), :) = [];
cond=6;
g=0;

figure
hold on
g=[1 3.5 6 8.5 11 13.5 16];

for c=1:cond
    
    bar(g(c),nanmean(Inter_Epoachs(:,c)),'facecolor','r')
    bar(g(c)+1,nanmean(Sav_Epoachs(:,c)),'facecolor','b')
    errorbar(g(c),mean(Inter_Epoachs(:,c)),std(Inter_Epoachs(:,c))/sqrt(length(Inter_Epoachs(:,c))),'.','LineWidth',2,'Color','k')
    errorbar(g(c)+1,mean(Sav_Epoachs(:,c)),std(Sav_Epoachs(:,c))/sqrt(length(Sav_Epoachs(:,c))),'.','LineWidth',2,'Color','k')
    
    for sub=1:size(Inter_Epoachs,1)
        plot(g(c),Inter_Epoachs(sub,c),'*','MarkerFaceColor',colorOrder(sub,:))
    end
    
    for sub=1:size(Sav_Epoachs,1)
        plot(g(c)+1,Sav_Epoachs(sub,c),'*','MarkerFaceColor',colorOrder(sub,:))
        
    end
    
end

set(gca,'Xtick',g+.5,'XTickLabel',{'Early A1' 'SS A1','Early A2', 'SS A2', 'Washout 5', 'Washout All'},'fontSize',12)
legend('Interference','Savings')
ylabel(param{1})




%%
%Number of strides to mid perturation 
%Line Plots

figure
title('Interference')
hold on
SubjInter=[setdiff(1:9, 2)];
for i=1:8
    Inames=['I' num2str(SubjInter(i))];
    plot(1,NumbStridesToMidPert(i,1),'ob',2,NumbStridesToMidPert(i,2),'or')
    if NumbStridesToMidPert(i,1)>NumbStridesToMidPert(i,2)
        line([1 2],[NumbStridesToMidPert(i,1) NumbStridesToMidPert(i,2)],'Color','green','LineStyle','-')
    else
        line([1 2],[NumbStridesToMidPert(i,1) NumbStridesToMidPert(i,2)],'Color','black','LineStyle','-')
    end
    text(1.001,NumbStridesToMidPert(i,1),['\leftarrow' Inames])
    
end
ylabel(param)
legend('A1','A2')


figure
title('Savings')
hold on
SubjSav=[1:9];
for i=1:9
    Snames=['S' num2str(SubjSav(i))];
    plot(1,NumbStridesToMidPert(i,3),'ob',2,NumbStridesToMidPert(i,4),'or')
    if NumbStridesToMidPert(i,3)>NumbStridesToMidPert(i,4)
        line([1 2],[NumbStridesToMidPert(i,3) NumbStridesToMidPert(i,4)],'Color','green','LineStyle','-')
    else
        line([1 2],[NumbStridesToMidPert(i,3) NumbStridesToMidPert(i,4)],'Color','black','LineStyle','-')
        
    end
    text(1.001,NumbStridesToMidPert(i,3),['\leftarrow' Snames])
    
end
legend('A1','A2')
ylabel(param)

%%
%Number of strides to mid perturation 
% Bar plots
AvgNOSToMidPert=nanmean(NumbStridesToMidPert,1);
SENOSToMidPer=nanstd(NumbStridesToMidPert,1)./sqrt(length(NumbStridesToMidPert));

figure
hold on
g=0;
for gr=1:3:6
 g=g+1;
        
 bar(gr,AvgNOSToMidPert(g),'facecolor','r')
 bar(gr+1,AvgNOSToMidPert(g+2),'facecolor','b')
 for sub=1:length(NumbStridesToMidPert)
     plot(gr,NumbStridesToMidPert(sub,g),'*','MarkerFaceColor',colorOrder(sub,:))
 end
 for sub=1:length(NumbStridesToMidPert)
     plot(gr+1,NumbStridesToMidPert(sub,g+2),'*','MarkerFaceColor',colorOrder(sub,:))
 end
 
 errorbar(gr,AvgNOSToMidPert(g),SENOSToMidPer(g),'.','LineWidth',2,'Color','k')
 errorbar(gr+1,AvgNOSToMidPert(g+2),SENOSToMidPer(g+2),'.','LineWidth',2,'Color','k')
 
end
legend('Interference','Savings')
ylabel('NumbStridesToMidPert')
title(param)
set(gca,'Xtick',[1.5 4.5],'XTickLabel',{'A1' 'A2'},'fontSize',12)






%%
%Rates plots

% pathToData= 'C:\Users\dum5\OneDrive\_Shared Drive - Interference Project - Alessandro - Dulce_\Params Files\subjectsData';
% load('StridesToRemove.mat')
% load('SlaAvgIdv.mat')
% load('NumbStridesToMidPert.mat')

INCR_EXP=1;
SINGLE_EXP=1;
FIT_MODEL=SINGLE_EXP;
npar_fit=3;
incrORdecr={'Increasing'};
subs = {9 , 8};
cond_Inter=[2,5];
cond_Sav=[2,4];

figure()
SubjInter=[setdiff(1:9, 2)];
for sub=1:subs{2}
    Inames=['I' num2str(SubjInter(sub))];
    
    A1min=nanmin( data_A1_Inter);
    A1max=nanmax(data_A1_Inter);
    A2min=nanmin( data_A2_Inter);
    A2max=nanmax(data_A2_Inter);
    absmin=min([A1min A2min]);
    absmax=nanmax([A1max A2max]);
    
    data_A1_Inter=slaI{sub,cond_Inter(1)};
    data_A1_Inter=data_A1_Inter(minNOSToRemove(sub,1):end,1);
    stride_vec_iA1=1:length(data_A1_Inter);
    subplot(2,4,sub)
    hold on
    title(['Interference  ', Inames ])
    plot(1:length(data_A1_Inter),data_A1_Inter,'or')
    
    data_A2_Inter=slaI{sub,cond_Inter(2)};
    data_A2_Inter=data_A2_Inter(minNOSToRemove(sub,1):end,1);
    stride_vec_iA2=1:length(data_A2_Inter);
    subplot(2,4,sub)
    plot(1:length(data_A2_Inter),data_A2_Inter,'ob')
    
    A1min=nanmin( data_A1_Inter);
    A1max=nanmax(data_A1_Inter);
    A2min=nanmin( data_A2_Inter);
    A2max=nanmax(data_A2_Inter);
    absmin=min([A1min A2min]);
    absmax=nanmax([A1max A2max]);
    line([NumbStridesToMidPert(sub,1) NumbStridesToMidPert(sub,1)],[absmin absmax],'Color','r' )
    line([NumbStridesToMidPert(sub,2) NumbStridesToMidPert(sub,2)],[absmin absmax],'Color','b')
    ylabel(param)
    axis tight
end
legend('A1','Mid perturabtion A1','A2','Mid perturabtion A2')



figure()
SubjSav=[1:9];
for sub=1:subs{1}
    Snames=['S' num2str(SubjSav(sub))];
    
    A1min=nanmin( data_A1_Sav);
    A1max=nanmax(data_A1_Sav);
    A2min=nanmin( data_A2_Sav);
    A2max=nanmax(data_A2_Sav);
    absmin=min([A1min A2min]);
    absmax=nanmax([A1max A2max]);
    
    subplot(3,3,sub)
    hold on
    title(['Savings  ', Snames ])
    data_A1_Sav=slaS{sub,cond_Sav(1)};
    data_A1_Sav=data_A1_Sav(minNOSToRemove(sub,2):end,1);
    plot(1:length(data_A1_Sav),data_A1_Sav,'or')
    
    
    
    
    data_A2_Sav=slaS{sub,cond_Sav(2)};
    data_A2_Sav=data_A2_Sav(minNOSToRemove(sub,2):end,1);
    subplot(3,3,sub)
    plot(1:length(data_A2_Sav),data_A2_Sav,'ob')
    
    A1min=nanmin( data_A1_Sav);
    A1max=nanmax(data_A1_Sav);
    A2min=nanmin( data_A2_Sav);
    A2max=nanmax(data_A2_Sav);
    absmin=min([A1min A2min]);
    absmax=nanmax([A1max A2max]);
    line([NumbStridesToMidPert(sub,3) NumbStridesToMidPert(sub,3)],[absmin absmax],'Color','r' )
    line([NumbStridesToMidPert(sub,4) NumbStridesToMidPert(sub,4)],[absmin absmax],'Color','b' )
    ylabel(param)
    axis tight
    
end
legend('A1','Mid perturabtion A1','A2','Mid perturabtion A2')






