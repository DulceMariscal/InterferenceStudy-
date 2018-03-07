%% Fit different models to the interference experiment data and compare them

clc
close all
clear all

NORMPERT = 0;

my_set_default(20,3,5);
pathToData = 'C:\Users\salat\OneDrive\Research\Code\!!! Data Prep For Fitting !!!';
if NORMPERT==0
    load([pathToData '\Data\DataForFitting.mat'])
else
    load([pathToData '\Data\DataForFittingNormPert.mat'])
end
load([pathToData '\defaultColors.mat'])

%% General Initializations
gnames = {'Interference','Savings'};
modnames = {'SM2006','SSIM','ETM'};
grcols = [defaultColors(2,:); defaultColors(1,:)];
%Colors fro grCond: Groups X Conditions
grCond{1,1} = defaultColors(3,:); grCond{2,1} = defaultColors(3,:);
grCond{1,2} = defaultColors(2,:); grCond{2,2} = defaultColors(1,:);
nsig = 3; %Number of significant digits (aka precision)
[ng, nt] = size(perturbations);
a1ind = [151:750];
a2ind = [2101:2700];
aeind = [2701:2850];
ainds = [a1ind; a2ind];
lstyles = {'o', '--'; 'o', '-'};  DATA=1; FIT=2; %Data, Fit
cols2   = [.5 .5 .5; 0 0 0]; %A1, A2
lwidths = [6, 3];
myLineSpecs = {{':','--'};{':','--','-.','-.'};{'-.','-.'}};
etCols = [defaultColors(5,:); defaultColors(4,:)];
nmod = length(modnames); %N models
cX=1; cY=2; A1=1; A2=2;
nconds=2;
indok = ~isnan(adaptations);
modParams = initializeModParams(perturbations,nt,ng,nconds);
transitions = [151 750 2101 2700];

%% -------------------------------------------------------------------------
%% Select 3 models you want to compare
% SM2006 = 1, %SSIM = 2, %SSIM - CWSS = 7, %SSIM - ISE2 = 13
% modelsToCompare = [2 7 13];
% modelsToCompare = [2 13 15];
% modelsToCompare = [1 2 3];
modelsToCompare = [16];

%% -------------------------------------------------------------------------

%% Other initializations
paramsE6 = [.87 .88 .74 .17 0]; %Parameters experiment 6 ETM
%params = [alpha0 alpha180 beta0 beta180 c180];
nmods=length(modelsToCompare);
BICs=[];
TAUdata = zeros(ng,nconds);
PVdata = zeros(ng,1);
TAUfit = zeros(ng,nconds);
PVfit = zeros(ng,1);

%% Main Fitting loop
% model = 3;
% gr = 1;
% [z, e, X] = my_model_evolution(model, modParams(model).kmin, modParams(model).input(gr,:,:), modParams(model).x0);

for imod=1:nmods
    model = modelsToCompare(imod);
    
    %Fit current model
    residFcn =  @(k)two_state_lsq2(k,adaptations,modParams(model).x0,modParams(model).input,model,basStd);
    [k_min] = fmincon(residFcn, modParams(model).k0, modParams(model).Acond,...
        modParams(model).bcond , [], [],modParams(model).klb,...
        modParams(model).kub);
    %     k_min = paramsE6;
    
    %Retrieve best-fitting motor output for each group
    for gr = 1:ng
        [z, e, X, Xadd]= my_model_evolution(model,k_min, modParams(model).input(gr,:,:), modParams(model).x0, basStd);
        modParams(model).z(gr,:,:) = z';
        modParams(model).e(gr,indok(gr,:),1) = (adaptations(gr,indok(gr,:)) - z(1,indok(gr,:)))'; % This here is the fitting error
        modParams(model).X(gr,:,:) = X;
        modParams(model).Xadd(gr,:,:) = Xadd;
    end
    
    %Compute and store paramaters goodness of fit
    modParams(model).BIC = computeBIC(adaptations, modParams(model).z, modParams(model).npars);
    modParams(model).R2 = computeR2(adaptations, modParams(model).z);
    modParams(model).kmin = k_min;
    BICs = [BICs modParams(model).BIC];
end

%% Create strings parameters=value
nsig = 4;
for imod=1:nmods
    model = modelsToCompare(imod);
    modParams(model).summaryString = createString(modParams(model).parNames,modParams(model).kmin, nsig);
end

%% Plot whole eperiment results for each group
f1 = figure('Name','Data + Fits','NumberTitle','off');
M1=[1 2; 3 4; 5 6];
for imod=1:nmods
    model = modelsToCompare(imod);
    
    Coord = modParams(model).zcoord;
    for gr=1:ng
        subplot(3,2,M1(imod,gr))
        % Data
        plot(adaptations(gr,:),'ko'); hold on;
        
        % Fit
        plot(modParams(model).z(gr,:,Coord),'Color',grcols(gr,:)); hold on;
        
        % Error
        %         stem(modParams(model).e(gr,:,:),'Color',[.7 .7 .7]); hold on;
        if gr==1
            ylabel(modParams(model).name);
        end
        if imod==1
            title(gnames{gr});
        end
        xlim([0 nt]);
        ylim([-1 1]);
    end
end

[deltaBics, Ranking] = myComputeBICdiff(BICs);

%% Plot whole experiment results for each experiment (overlapping groups)
f2 = figure('Name','Fits - Groups contrasts','NumberTitle','off');
M1=[1 2; 3 4; 5 6];
for imod=1:nmods
    model = modelsToCompare(imod);
    Coord = modParams(model).zcoord;
    subplot(3,1,imod)
    for gr=1:ng
        % Fit
        hold on
        plot(modParams(model).z(gr,:,Coord),'Color',grcols(gr,:));
        if gr==1
            lab = sprintf('%s \n Rank: %s  %sBIC= %s',num2str(modParams(model).name,nsig), ...
                num2str(Ranking(imod),nsig),...
                '\Delta',...
                num2str(deltaBics(imod),nsig));
            ylabel(lab);
        end
        xlim([0 nt]);
        ylim([-1 1]);
    end
end

%% Plot A1 VS A2 - Data VS Fits
stridesToSkip = 2; %These strides will be removed to compute percent variations
f3 = figure('Name',['A1 vs A2 contrasts - Strides Removed For PV : ' num2str(stridesToSkip)  ],...
    'NumberTitle','off');
M2 = [1 2 3; 4 5 6; 7 8 9] ;
htop = zeros(1,2);
for imod=1:nmods
    model = modelsToCompare(imod);
    % Columns 1 and 2------------------------------------------------------
    Coord = modParams(model).zcoord;
    for gr=1:ng
        for cond = 1:nconds
            ci = ainds(cond,:); %Current indices
            subplot(nmods, nconds + 1, M2(imod,gr))
            
            %Data----------------------------------------------------------
            cdata = adaptations(gr,ci);
            plot(cdata,lstyles{cond,DATA},'Color',cols2(cond,:),'LineWidth',1), hold on;
            if imod==1
                coeffs = my_exp_fit2(cdata(stridesToSkip + 1:end),'Increasing','Single');
                TAUdata(gr,cond) = coeffs(3);
            end
            
            %Fit-----------------------------------------------------------
            cfit = modParams(model).z(gr,ci,Coord);
            htop(cond) = plot(cfit,...
                lstyles{cond,FIT},'Color',grCond{gr,cond},'LineWidth',lwidths(cond)); hold on;
            coeffs = my_exp_fit2(cfit(stridesToSkip + 1 : end),'Increasing','Single');
            TAUfit(model,gr,cond) = coeffs(3);
        end
        ylim([0.7 1 ]); xlim([0 600]);
        
        %Compute PV
        PVdata = percent_variation(TAUdata);
        text(0.5, 0.3,['PV_{data} = ' num2str(PVdata(gr)*100,nsig) '%'],'Units','normalized','FOntSIze',18)
        
        PVfit(model,gr,cond) = percent_variation(squeeze(TAUfit(model,gr,:))');
        text(0.5, 0.15,['PV_{fit} = ' num2str(PVfit(model,gr,cond)*100,nsig) '%'],'Units','normalized','FOntSIze',18)
        
        if imod==1
            title(gnames{gr})
        end
        %            'Data^S','x^S','x_f^S','x_s^S')
        if gr==1
            ylabel(modParams(model).name);
        end
        if imod==1
            legend('Data_{A_1}','Fit_{A_1}','Data_{A_2}','Fit_{A_2}')
        end
        uistack(htop(1), 'top')
        
    end
    
    
    
    % Column 3-------------------------------------------------------------
    subplot(nmods, nconds + 1, M2(imod,gr + 1))
    %A1(Naive)
    htopp = plot(modParams(model).z(SAV,ainds(A1,:),Coord),lstyles{A1,FIT},...
        'Color',grCond{gr,A1},'LineWidth',lwidths(A1)); hold on;
    %A2 Int and Sav
    for gr=1:2
        plot(modParams(model).z(gr,ainds(A2,:),Coord),'Color',grCond{gr,A2}); hold on;
    end
    ylim([0.7 1 ]); xlim([0 600]);
    if imod==1
        title('Comparison A1 - A2')
        legend('Fit_{A_1}','Fit_{A_2}^{INT}','Fit_{A_2}^{SAV}')
    end
    uistack(htopp, 'top')
    
end



%% Plot whole eperiment results for each group with individual states
f4 = figure('Name','Data + Fits + States','NumberTitle','off');
M1=[1 2; 3 4; 5 6];
for imod=1:nmods
    model = modelsToCompare(imod);
    Coord = modParams(model).zcoord;
    for gr=1:ng
        subplot(3,2,M1(imod,gr))
        
        % Add parameters values
        if gr==1
            text(-.35, 0.5, modParams(model).summaryString ,'Units','normalized','FontSIze',14,...
                'VerticalAlignment','middle','HorizontalAlignment','left',...
                'EdgeCOlor','b'); hold on;
        end
        
        % Data
        plot(adaptations(gr,:),'ko'); hold on;
        
        % Fit
        plot(modParams(model).z(gr,:,Coord),'Color',grcols(gr,:)); hold on;
        
        % States
        for cstate=1:modParams(model).nstates
            %Plot
            p = plot(squeeze(modParams(model).X(gr,cstate,:))'); hold on;
            
            %Retrieve plot properties
            ccol = modParams(model).colors{cstate};
            if isnan(ccol)
                ccol = grcols(gr,:);
            end
            clw = modParams(model).lw(cstate);
            if isnan(clw)
                clw = 2; %Default line width
            end
            
            %Set plot poperties
            p.LineStyle = modParams(model).lineSpecs{cstate};
            p.Color = ccol;
            p.LineWidth = clw;
        end
        
        if gr==1
            ylabel(modParams(model).name);
            legend(['Adapt_{data}','Adapt_{fit}', modParams(model).stNames])
        end
        if imod==1
            title(gnames{gr});
        end
        xlim([0 nt]);
        ylim([-1 1]);
    end
end

%% COMPARE MODELS ERRORS
f6 = figure('Name','Models errors','NumberTitle','off');
M = [1 2; 3 4];
legends1={}; legends2={};
% modCols = defaultColors(3:5,:);
CE = cell(3,2);
CSE = cell(3,2);
modCols = brewermap(3,'Set1');
for imod=1:nmods
    model = modelsToCompare(imod);
    Coord = modParams(model).zcoord;
    
    for gr=1:ng
        ce = []; cse=[];
        ce  = modParams(model).e(gr,:);
        indnan = isnan(ce);
        ce(1,indnan) = 0;
        cse = cumsum(abs(ce));
        
        % Errors----------------------------------------------------------
        subplot(2,2,M(1,gr))
        area(ce,'facecolor',modCols(imod,:),'facealpha',.5); hold on
        ylim([-.05 .05])
        if imod==1
            title(gnames(gr))
        end
        xlim([0 2850])
        
        % Cumulative errors-----------------------------------------------
        subplot(2,2,M(2,gr))
        plot(cse,'Color',modCols(imod,:)); hold on
        xlim([0 2850])
        CE{imod,gr} = ce;
        CSE{imod,gr} = cse;
    end
    legends1 = [legends1; modParams(model).name];
    legends2 = [legends2; [modParams(model).name ' - Rank: ' num2str(Ranking(imod)) ' - \DeltaBIC: ' num2str(deltaBics(imod))]];
    
end

subplot(2,2,1)
legend(legends1)
ylabel('Errors')
subplot(2,2,3)
legend('Data',legends2,'Location','northwest')
ylabel('Cum Errors')
ylim([0 60])
YL = ylim();
hold on, addTransitions(transitions,YL(1),YL(2))
subplot(2,2,4)
ylim([0 60])
YL = ylim();
hold on, addTransitions(transitions,YL(1),YL(2))


%% Plot A1 VS A2 with individual states and percent variations
%% Model 2
mspecs = {} ;
statesMarkers = {'*','o','s','d'};
model = modelsToCompare(end);
nstates = modParams(model).nstates;
Coord = modParams(model).zcoord;
t = 2:30:600;
ms=10;
f7 = figure('Name','States','NumberTitle','off');

%First two columns
for gr = 1 : ng
    subplot(1,3,gr)
    for cond = 1:nconds
        ci = ainds(cond,:); %Current indices
        
        % Extract data
        z = modParams(model).z(gr,ci,Coord);
        X = modParams(model).X(gr,:,ci);
        
        p = plot(z);hold on;
        if cond==1
            p.Color = [.5 .5 .5];
        else
            p.Color = grcols(gr,:);
        end
        
        for cstate=1:nstates
            x = squeeze(X(:,cstate,:));
            p = plot(t,x(t));
            p.Marker = statesMarkers{cstate};
            if cond==1
                p.Color = [.5 .5 .5];
            else
                p.Color = grcols(gr,:);
            end
            p.MarkerSize = ms;
        end
        
    end
    title(gnames{gr})
    
    
    %Third column
    subplot(1,3,3)
    for cond = 1:nconds
        ci = ainds(cond,:); %Current indices
        
        % Extract data
        z = modParams(model).z(gr,ci,Coord);
        X = modParams(model).X(gr,:,ci);
        
        p = plot(z);hold on;
        if cond==1
            p.Color = [.5 .5 .5];
        else
            p.Color = grcols(gr,:);
        end
        
        for cstate=1:nstates
            x = squeeze(X(:,cstate,:));
            p = plot(t,x(t));
            p.Marker = statesMarkers{cstate};
            if cond==1
                p.Color = [.5 .5 .5];
            else
                p.Color = grcols(gr,:);
            end
            p.MarkerSize = ms;
        end
        xlim([0 600])
        ylim([-0.4 1])
    end
    title('Comparison A1 VS A2')
    
end


subplot(1,3,1)
legend ('z^{A1}', 'x_{f}^{A1}','x_{s}^{A1}','x_{pe}^{A1}','x_{ne^{A1}}' ,...
    'z^{A2}', 'x_{f}^{A2}','x_{s}^{A2}','x_{pe}^{A2}','x_{ne^{A2}}')

subplot(1,3,2)
legend ('z^{A1}', 'x_{f}^{A1}','x_{s}^{A1}','x_{pe}^{A1}','x_{ne^{A1}}' ,...
    'z^{A2}', 'x_{f}^{A2}','x_{s}^{A2}','x_{pe}^{A2}','x_{ne^{A2}}')



%% Analyze LISU TV BEHAVIOR
model = 14;
if ~all(isnan(modParams(model).kmin))
    
    % Plot whole eperiment results for each group with individual states
    f8 = figure('Name','Data + Fits + States + Xadd','NumberTitle','off');
    M1=[1 2; 3 4];
    
    Coord = modParams(model).zcoord;
    for gr=1:ng
        subplot(2,2,M1(1,gr))
        
        % Add parameters values
        if gr==1
            text(-.35, 0.5, modParams(model).summaryString ,'Units','normalized','FontSIze',14,...
                'VerticalAlignment','middle','HorizontalAlignment','left',...
                'EdgeCOlor','b'); hold on;
        end
        
        % Data
        plot(adaptations(gr,:),'ko'); hold on;
        
        % Fit
        plot(modParams(model).z(gr,:,Coord),'Color',grcols(gr,:)); hold on;
        
        % States
        for cstate=1:modParams(model).nstates
            %Plot
            p = plot(squeeze(modParams(model).X(gr,cstate,:))'); hold on;
            
            %Retrieve plot properties
            ccol = modParams(model).colors{cstate};
            if isnan(ccol)
                ccol = grcols(gr,:);
            end
            clw = modParams(model).lw(cstate);
            if isnan(clw)
                clw = 2; %Default line width
            end
            
            %Set plot poperties
            p.LineStyle = modParams(model).lineSpecs{cstate};
            p.Color = ccol;
            p.LineWidth = clw;
        end
        
        % Additional variables
        subplot(2,2,M1(2,gr))
        Xadds = modParams(model).Xadd;
        nas=size(Xadds,2);
        for jj=1:nas
            plot(squeeze(Xadds(gr,jj,:))); hold on;
        end
        
        if gr==1
%             ylabel(modParams(model).name);
            legend( modParams(model).avDescr)
        end
        if imod==1
            title(gnames{gr});
        end
        xlim([0 nt]);
        ylim([-1 1]);
    end
    
    
    
end










