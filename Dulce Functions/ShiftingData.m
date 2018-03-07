%To shift the data 
cd('C:\Users\dum5\OneDrive\_Shared Drive - Interference Project - Alessandro - Dulce_\Params Files\subjectsData')

% clc
% close all 
clear all
load('SlaAvgIdv.mat')
% load('minNOSToRemoveAle.mat')
% load('StridesToRemove.mat')
% gait_par='StepLengthAsym';
% 
%for ExtAdaptaiton= SLA-SV
% load('ExtAdaptAllData.mat') 
% % % load('StridesToRemoveExtAdpat.mat')
% gait_par='ExtAdaptWithHip';

%For step length asymmetry
% load('StridesToRemove.mat')
load('NoRemovingAnyStrides.mat')
% load('stepLengthAsymAllData.mat')
param={'stepLengthAsym'};
% eval(['slaI=[' param 'INT];'])
% eval(['slaS=[' param 'SAV];'])

%for ExtAdaptaiton= SLA-SV
% load('ExtAdaptNorm2AllData.mat') 
% % load('StridesToRemoveExtAdpat.mat')
% param={'ExtAdaptNorm2'};
% slaI=ExtAdaptINT;
% slaS=ExtAdaptSAV;

%No hip 
% load('ExtAdaptPNormAllData.mat') 
% load('StridesToRemoveExtAdpatPNomr')
% param={'ExtAdaptPNorm'};
% slaI=ExtAdaptINT;
% slaS=ExtAdaptSAV;


cond_Inter=1:6;
cond_Sav=1:5;
charGroups = {'S','I'};
indSubs = {setdiff(1:9, []), setdiff(1:8, [])};
% indSubs = {setdiff(1:9, [2 5]), setdiff(1:8, [1 7])};
subs = {length(indSubs{1}) , length(indSubs{2})};

for i=1:subs{2}
    s=indSubs{2}(i);
    %     s=i;
    data_A1_Inter=slaI{s,cond_Inter(2)};
    INTERFERENCE_DATA{1,i}=data_A1_Inter(minNOSToRemove(s,1):end,1);
    
    data_A2_Inter=slaI{i,cond_Inter(5)};
    INTERFERENCE_DATA{2,i}=data_A2_Inter(minNOSToRemove(s,1):end,1);
    
    
end



for i=1:subs{1}
    s=indSubs{1}(i);
    %     s=i;
    data_A1_Sav=slaS{s,cond_Sav(2)};
    SAVINGS_DATA{1,i}=data_A1_Sav(minNOSToRemove(s,2):end,1);
    
    data_A2_Sav=slaS{s,cond_Sav(4)};
    SAVINGS_DATA{2,i}=data_A2_Sav(minNOSToRemove(s,2):end,1);
    
end


INCR_EXP=1;
SINGLE_EXP=1;
FIT_MODEL=SINGLE_EXP;
npar_fit=3;
initial_conditions_fit
INTERFERENCE=1; SAVINGS=2;
nos=[8 8];
% Adaptation 1-------------------------------------------------------------
[iA1m, iA1SE] = average_curves(INTERFERENCE_DATA(1,:)); %avg_adapt_curve_int [iA1m,iA2std]
[sA1m, sA1std] = average_curves(SAVINGS_DATA(1,:)); %avg_adapt_curve_sav [sA1m,sA2std]
[iA2m, iA2std] = average_curves(INTERFERENCE_DATA(2,:)); %avg_adapt_curve_int [iA1m,iA2std]
[sA2m, sA2se] = average_curves(SAVINGS_DATA(2,:)); %avg_adapt_curve_sav [sA1m,sA2std]




INT_min=nanmin(iA1m);
SAV_min=nanmin(sA1m);
MinValAdapt=min([INT_min SAV_min]);


INT_max=nanmax(iA2m);
SAV_max=nanmax(sA2m);
MaxValAdapt=max([INT_max SAV_max])-MinValAdapt;

for cond=1:length(cond_Inter)
for i=1:subs{2} %8
    s=indSubs{2}(i);
    %     s=i;
   data_A1_Inter=slaI{s,cond_Inter(1)}; 
   INTERFERENCE_DATA{1,i}=(data_A1_Inter(minNOSToRemove(s,1):end,1));
    
   data_A2_Inter=slaI{s,cond_Inter(2)}; 
   INTERFERENCE_DATA{2,i}=(data_A2_Inter(minNOSToRemove(s,1):end,1)-MinValAdapt)./MaxValAdapt;
   
 
   data_A3_Inter=slaI{s,cond_Inter(3)};
   INTERFERENCE_DATA{3,i}=(data_A3_Inter(minNOSToRemove(s,1):end,1)+MinValAdapt)./MaxValAdapt; 
   
   data_A4_Inter=slaI{s,cond_Inter(4)};
   newp=linspace(MinValAdapt,0,600);
   newPert=[newp zeros(1,length(data_A4_Inter(minNOSToRemove(s,1):end,1))-600)]';
   INTERFERENCE_DATA{4,i}=((data_A4_Inter(minNOSToRemove(s,1):end,1)+newPert))./MaxValAdapt; 
%    newp=nan(1,600);
%    newPert=0;
%    INTERFERENCE_DATA{4,i}=((data_A4_Inter(minNOSToRemove(i,1):end,1)+MinValAdapt))./MaxValAdapt; 
%    INTERFERENCE_DATA{4,i}=((data_A4_Inter(minNOSToRemove(s,1):end,1)+newPert))./MaxValAdapt; 
   
    data_A5_Inter=slaI{s,cond_Inter(5)};
   INTERFERENCE_DATA{5,i}=(data_A5_Inter(minNOSToRemove(s,1):end,1)-MinValAdapt)./MaxValAdapt; 
   
   data_A6_Inter=slaI{s,cond_Inter(6)};
   INTERFERENCE_DATA{6,i}=(data_A6_Inter(minNOSToRemove(s,1):end,1)); 
  
end 
end



for cond=1:length(cond_Sav)
for i=1:subs{1}%9
    s=indSubs{1}(i);
   
   data_A1_Sav=slaS{s,cond_Sav(1)}; 
   SAVINGS_DATA{1,i}=(data_A1_Sav(minNOSToRemove(s,2):end,1));
   
   data_A2_Sav=slaS{s,cond_Sav(2)}; 
   SAVINGS_DATA{2,i}=(data_A2_Sav(minNOSToRemove(s,2):end,1)-MinValAdapt)/MaxValAdapt;
   
   data_A3_Sav=slaS{s,cond_Sav(3)}; 
   SAVINGS_DATA{3,i}=(data_A3_Sav(minNOSToRemove(s,2):end,1));
   
   data_A4_Sav=slaS{s,cond_Sav(4)};
   SAVINGS_DATA{4,i}=(data_A4_Sav(minNOSToRemove(s,2):end,1)-MinValAdapt)/MaxValAdapt;
   
   
   if s==5 
    SAVINGS_DATA{5,i}=nan(1000,1);
   else
   data_A5_Sav=slaS{s,cond_Sav(5)}; 
   SAVINGS_DATA{5,i}=(data_A5_Sav(minNOSToRemove(s,2):end,1));
   end
    
end
end

cond=[6 5];
SLA_INT_shifted=[];
SLA_INT_SSEshifted=[];


for c=1:cond(1)
  [iA1m, iA1SE] = average_curves(INTERFERENCE_DATA(c,:));  
  SLA_INT_shifted=[SLA_INT_shifted iA1m  nan(1,5)];
  SLA_INT_SSEshifted=[SLA_INT_SSEshifted iA1SE nan(1,5)];
    
end

SLA_SAV_shifted=[];
SLA_SAV_SEshifted=[];

for c=1:cond(2)
  [sA1m, sA2se] = average_curves(SAVINGS_DATA(c,:));  
  SLA_SAV_shifted=[SLA_SAV_shifted sA1m nan(1,5)];
  SLA_SAV_SEshifted=[SLA_SAV_SEshifted sA2se nan(1,5)];
    
end


slaS=SAVINGS_DATA';
slaI=INTERFERENCE_DATA';

figure
hold on 
scatter(1:length(SLA_INT_shifted),SLA_INT_shifted,'r')
scatter(1:length(SLA_SAV_shifted),SLA_SAV_shifted,'b')
% h1=boundedline(1:length(SLA_INT_shifted),SLA_INT_shifted,SLA_INT_SSEshifted,'r',1:length(SLA_SAV_shifted),SLA_SAV_shifted,SLA_SAV_SEshifted,'b');
legend('INT','SAV')
axis tight
ylabel([param 'Shifted'])
xlabel('strides')



% save([param{1} 'Shifted'], 'slaS','slaI','SLA_SAV_SEshifted','SLA_INT_SSEshifted')

% h1=boundedline(1:length(SLA_INT_shifted),SLA_INT_shifted,SLA_INT_STDshifted,'r')
% h2=boundedline(1:length(SLA_SAV_shifted),SLA_SAV_shifted,SLA_SAV_STDshifted,'b')

%,x,eval([gaitPar 'Avg(SAV,:)']),slamse(SAV,:),'b');
% plot(1:length(SLA_INT_shifted),SLA_INT_shifted,1:length(SLA_INT_shifted),SLA_INT_STDshifted)
% hold on 
% plot(1:length(SLA_SAV_shifted),SLA_SAV_shifted,1:length(SLA_SAV_shifted),SLA_SAV_STDshifted)





