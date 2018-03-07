clc
close all
clear all

%% Parameters
pathToData = 'C:\Users\salat\OneDrive\Documents\MATLAB\Research\CoreL\Interference\Interf_study_AS\subjectsData';
groups = 2;
conditions = 2;
% subs = [8 9];

indSubs = {setdiff(1:9, [1 2 8]), setdiff(1:9, 2)};
subs = [length(indSubs{1})  length(indSubs{2})];
charGroups = {'I','S'};
% params={'spatialContributionNorm2','stepTimeContributionNorm2','velocityContributionNorm2','stepLengthAsym'};
% params={'spatialContributionNorm2','stepTimeContributionNorm2','velocityContributionNorm2','netContributionNorm2'};
% params={'spatialContributionPNorm','stepTimeContributionPNorm','velocityContributionPNorm','netContributionPNorm'};
% gaitPar=params{3};
gaitPar = 'stepLengthAsym';
wSize = 5;
% woSizes = zeros(groups,subs);
REMOVE_BASELINE = 1;
colors = distinguishable_colors(15);
nconds = [6 5];
INT=1; SAV=2;


% slaS = cell(subs, nconds(SAV)); % 5 conditions in total
% slaI = cell(subs, nconds(INT)); % 6 conditions in total

eval([gaitPar 'SAV =cell(subs(2), nconds(SAV));']); %DUL
eval([gaitPar 'INT =cell(subs(1), nconds(INT));']);

% gaitPar.Sav= cell(subs(1), nconds(INT));
stridesPerCond = {[150 600 600 750 600 150], [150 600 1350 600 150]};
slam = nan(groups, subs(2), sum(stridesPerCond{1}) ) ; %Each element is a matrix

indMatS = zeros(nconds(SAV),2); %Start and stop of each condition
indMatI = zeros(nconds(INT),2); %Start and stop of each condition

%% Main loop
% figure
for gr = 1:groups
    for sub = 1:subs(gr)
        %% Load the data
        cname = [charGroups{gr} '00' num2str(indSubs{gr}(sub)) 'params.mat'];
        load(cname);
        
        %% Remove baseline
%         if REMOVE_BASELINE
%             if (indSubs{gr}(sub) <= 5)
%                 baseline_cond='Baseline';
%             else
%                 baseline_cond='TM base';
%             end
%             [adaptData]=adaptData.removeBias(baseline_cond);
              
            [adaptData]=adaptData.removeBadStrides.removeBias;
%         end
        
        %% Get all conditions' names
        conditions = adaptData.metaData.conditionName;
        
        %% Extract all epochs
        for cond=1:nconds(gr)
            %             if gr==2 && sub==5 && cond==5
            %                 adaptData.getParamInCond(gaitPar, conditions{cond})
            %             end
            %% Extract current condition
            if (strcmp(adaptData.subData.ID,'I001') && strcmp( gaitPar ,'stepLengthAsym') || strcmp(adaptData.subData.ID,'I008')) && cond==2 && strcmp( gaitPar ,'stepLengthAsym')
                load(['abrupt_A1_' adaptData.subData.ID '.mat'])
                cc = a1;
                warning(['Group = ' num2str(gr) ', Sub = ' num2str(sub) ', Cond = ' conditions{cond} ...
                    ' These data required further processing ']);
            else
                cc = adaptData.getParamInCond(gaitPar, conditions{cond});
            end
            cl = length(cc);
            cflag = (sum(isnan(cc)) >= 0.5*cl) || cl < 100; %Bad subjects
            if cflag == 1
                cc = [];
                warning(['Group = ' num2str(gr) ', Sub = ' num2str(sub) ', Cond = ' conditions{cond} ...
                    ' has been discarded because it contains more than 50% of nans ']);
            else
                %% Filter
                cc = medfilt1(cc,wSize);
                firstSamples = cc(1:wSize-1);
                cc = tsmovavg(cc,'s',wSize,1);
                cc(1:wSize-1) = firstSamples;
            end
            
            %% Store current condition
            if cond == 1
                indStart = 1;
            else
                indStart = indStop + 1;
                %indStart = stridesPerCond{gr}(cond-1) + 1;
            end
            indStop  =  indStart + stridesPerCond{gr}(cond) - 1;
            
            if gr == INT
                eval([gaitPar 'INT{sub,cond} =cc;']);
%                 [gaitPar, 'INT'] (sub,cond) =cc;
%                 gaitPar.Inter{sub,cond} = cc;
%                 slaI{tsub,cond} = cc;                %Store in ca
                indMatI(cond,:) = [indStart indStop];
            else
                
                 eval([gaitPar 'SAV{sub,cond} =cc;']);
%                 gaitPar.Sav{sub,cond} = cc;
%                 slaS{sub,cond} = cc;
                indMatS(cond,:) = [indStart indStop];
            end
            slam(gr,sub,indStart:indStop) = mynanPad(cc, stridesPerCond{gr}(cond)); %Store in a matrix
            %             woSizes(gr,sub) = length(cc);
            
            %             hold on
            %             if gr == 1
            %                 plot(cc,'r','Linewidth',1)
            %             else
            %                 plot(cc,'b','Linewidth',1)
            %             end
        end
    end
end

% save([gaitPar 'AllData'],[gaitPar 'SAV'],[gaitPar 'INT'] )

%%
eval([gaitPar 'Avg' '= squeeze(nanmean(slam,2));']) %Average across subjects
% indA1 = [151:151+599];
% [deltaSla] = findDeltaSla(slamavg(:,indA1));

slamse  = nanse(slam,2);
x=1:2850;

%% Plot averages
% figure
% plot(mS)
% hold on
% plot(mI)
% hold on
figure
% h = boundedline(x,slamavg(INT,:),slamse(INT,:),'r',x,slamavg(SAV,:),slamse(SAV,:),'b');
h = boundedline(x,eval([gaitPar 'Avg(INT,:)']),slamse(INT,:),'r',x,eval([gaitPar 'Avg(SAV,:)']),slamse(SAV,:),'b');
title('Whole Exp Adaptation')
legend(h, 'Interference','Savings')
ylabel(gaitPar)
axis tight
% save([gaitPar 'AllData'],[gaitPar 'SAV'],[gaitPar 'INT'], [gaitPar 'Avg'] )
%% Contrast normalization VS no normalization 
% Since interference people seem to be more perturbed during A1,
% I am assuming that they are generally more perturbed than savings people
% eqIntPert = 1;
% pertS = [ zeros(1,150) ones(1,600) zeros(1,1350) ones(1,600) zeros(1,150)];
% pertI = [ zeros(1,150) eqIntPert*ones(1,600) (-eqIntPert)*ones(1,600),...
%     linspace(-eqIntPert,0,600) zeros(1,150) eqIntPert*ones(1,600) zeros(1,150)];
% pBas = [pertI; pertS];
% 
% eqIntPert = 1 + deltaSla;
% pertS = [ zeros(1,150) ones(1,600) zeros(1,1350) ones(1,600) zeros(1,150)];
% pertI = [ zeros(1,150) eqIntPert*ones(1,600) (-eqIntPert)*ones(1,600),...
%     linspace(-eqIntPert,0,600) zeros(1,150) eqIntPert*ones(1,600) zeros(1,150)];
% pNorm = [pertI; pertS];
% 
% f1=figure;
% M = [1 2; 3 4]; % Group X Normalization
% grCols = {'r','b'};
% for norm=1:2
%     if norm==1
%         p2 = pBas;
%         p1 = pBas;
%     else
%         p2 = pBas;
%         p1 = pNorm;
%     end
% 
%     adaptations = p1 + slamavg;
%     
%     for gr=1:2
%         subplot(3,2,M(gr,norm))
%         plot(adaptations(gr,:),grCols{gr}), hold on;
%         plot(p2(gr,:),'k')
%         ylim([-1 1]*1.02)
%         
%             
%         subplot(3,2,5+norm-1)
%         plot(adaptations(gr,:),grCols{gr}); hold on
%     end
%     ylim([-1 1])
% 
%     
% end
% 
% perturbations = pBas;
%% Save data for fitting

%Cocatenated vectors
% adaptation = adaptations(:);
% perturbation = perturbations(:);

% 
% save('Data/DataForFittingNormPert.mat','adaptations','slamavg','slamse','perturbations','INT','SAV',...
%     'adaptation','perturbation')

%Note: Plural forms are matrices. Singular forms are vectors

