%Extend of Adaptation 
clear all;close all;clc
% load('velocityContributionNorm2AllData.mat')
% load('netContributionNorm2AllData.mat')
load('velocityContributionPNormAllData.mat')
load('netContributionPNormAllData.mat')
load('NoRemovingAnyStrides.mat')

groups = 2;
conditions = 2;
% subs = [8 9];
indSubs = {setdiff(1:8, [1 7]), setdiff(1:9, 2)};
subs = [length(indSubs{1}) , length(indSubs{2})];
charGroups = {'INT','SAV'};
nconds = [6 5];
INT=1; SAV=2;

ExtAdaptSAV=cell(subs(2), nconds(SAV));
ExtAdaptINT=cell(subs(1), nconds(INT));



for gr=1:groups
for s=1:subs(gr)
    sub=indSubs{gr}(s);
%     sub=s;
for cond=1:nconds(gr) 

   eval(['ExtAdapt' charGroups{gr} '{s,cond}=netContributionPNorm' charGroups{gr} '{sub,cond}-velocityContributionPNorm' charGroups{gr} '{sub,cond};']) 
   
end
end
end

ExtAdaptAvg=netContributionPNormAvg-velocityContributionPNormAvg;

params={'netContributionPNorm','velocityContributionPNorm'};
cond_Inter=1:6;
cond_Sav=1:5;


for p=1:length(params)
    
    
eval(['slaI=[' params{p} 'INT];'])
eval(['slaS=[' params{p} 'SAV];'])


for i=1:8
    
    load([params{p} 'AllData.mat'])

    
   data_A1_Inter=slaI{i,cond_Inter(1)}; 
   INTERFERENCE_DATA{1,i,p}=(data_A1_Inter(minNOSToRemove(i,1):end,1));
    
   data_A2_Inter=slaI{i,cond_Inter(2)}; 
   INTERFERENCE_DATA{2,i,p}=(data_A2_Inter(minNOSToRemove(i,1):end));
   
 
   data_A3_Inter=slaI{i,cond_Inter(3)};
   INTERFERENCE_DATA{3,i,p}=(data_A3_Inter(minNOSToRemove(i,1):end,1)); 
   
   data_A4_Inter=slaI{i,cond_Inter(4)};
   INTERFERENCE_DATA{4,i,p}=(data_A4_Inter(minNOSToRemove(i,1):end,1)); 
   
    data_A5_Inter=slaI{i,cond_Inter(5)};
   INTERFERENCE_DATA{5,i,p}=(data_A5_Inter(minNOSToRemove(i,1):end,1)); 
   
   data_A6_Inter=slaI{i,cond_Inter(6)};
   INTERFERENCE_DATA{6,i,p}=(data_A6_Inter(minNOSToRemove(i,1):end,1)); 
  

end



for s=1:9
   
   data_A1_Sav=slaS{s,cond_Sav(1)}; 
   SAVINGS_DATA{1,s,p}=(data_A1_Sav(minNOSToRemove(s,2):end,1));
   
   data_A2_Sav=slaS{s,cond_Sav(2)}; 
   SAVINGS_DATA{2,s,p}=(data_A2_Sav(minNOSToRemove(s,2):end,1));
   
   data_A3_Sav=slaS{s,cond_Sav(3)}; 
   SAVINGS_DATA{3,s,p}=(data_A3_Sav(minNOSToRemove(s,2):end,1));
   
   data_A4_Sav=slaS{s,cond_Sav(4)};
   SAVINGS_DATA{4,s,p}=(data_A4_Sav(minNOSToRemove(s,2):end,1));
   
   
   if s==5 
    SAVINGS_DATA{5,s,p}=nan(1000,1);
   else
   data_A5_Sav=slaS{s,cond_Sav(5)}; 
   SAVINGS_DATA{5,s,p}=(data_A5_Sav(minNOSToRemove(s,2):end,1));
   end
    
end

end


cond=[6 5];
SLA_INT_shifted=[];
SLA_INT_SEshifted=[];
SLA_SAV_shifted=[];
SLA_SAV_SEshifted=[];


for p=1:length(params) 
    
SLA_INT_shifted=[];
SLA_INT_SEshifted=[];
SLA_SAV_shifted=[];
SLA_SAV_SEshifted=[];
for c=1:cond(1)
    
  [iA1m, iA1SE] = average_curves(INTERFERENCE_DATA(c,:,p));  
  SLA_INT_shifted=[SLA_INT_shifted iA1m];
  SLA_INT_SEshifted=[SLA_INT_SEshifted iA1SE];
    

end


for c=1:cond(2)

    [sA1m, sA2se] = average_curves(SAVINGS_DATA(c,:,p));  
  SLA_SAV_shifted=[SLA_SAV_shifted sA1m];
  SLA_SAV_SEshifted=[SLA_SAV_SEshifted sA2se];
    
  
end
if p==1
    VelocityAVG_INT= SLA_INT_shifted;
    VelocitySE_INT= SLA_INT_SEshifted;
    
    VelocityAVG_SAV= SLA_SAV_shifted;
    VelocitySE_SAV= SLA_SAV_SEshifted;
else
    
    NetAVG_INT=SLA_INT_shifted;
    NetSE_INT=SLA_INT_SEshifted;
    
     NetAVG_SAV= SLA_SAV_shifted;
    NetSE_SAV= SLA_SAV_SEshifted;
end

end

figure
subplot(2,1,1)
h1=boundedline(1:length( VelocityAVG_INT), VelocityAVG_INT,VelocitySE_INT,'r',1:length(NetAVG_INT), NetAVG_INT,NetSE_INT,'b');
axis tight
title('Interference')
% ylabel(params)
xlabel('strides')
subplot(2,1,2)
boundedline(1:length( VelocityAVG_SAV), VelocityAVG_SAV,VelocitySE_SAV,'r',1:length( NetAVG_SAV), NetAVG_SAV,NetSE_SAV,'b');
title('Savings')
% ylabel(params)
legend(h1,params{1},params{2})








% save('ExtAdaptPNormAllData','ExtAdaptSAV','ExtAdaptINT','ExtAdaptAvg')