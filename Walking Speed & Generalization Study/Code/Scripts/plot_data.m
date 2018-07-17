%% Explore data to check consistency

clc
close all
clear all

%% What to do with the script
PLOT_DATA = 1;
CREATE_REPORT  = 0;

%% How to do it
%Parameters to extract and plot
parameters = {'StepLengthAsym','singleStanceSpeedFastAbs','singleStanceSpeedSlowAbs'};
ylimits = [-0.3 0.3; -1.5 1.5; -1.5 1.5];

%% Define folder
pathToMainFolder = 'C:\Users\salat\OneDrive\Documents\MATLAB\Walking Speed & Generalization Study\';
folderToData = [pathToMainFolder 'Data\'];
groups={'OldSlow','OldFast','YoungSlow','YoungFast'};
ngroups=length(groups);

%% Define common conditions (those that eveybody did)
CCond{1} = {'OG base', 'TM base', 'adaptation', 'catch' ,'re-adaptation' ,   'OG post' ,   'TM post'} ;

%% Initializations
nsubs=zeros(1,ngroups);
CONDS_UNION=cell(1,ngroups);
CONDS_INTERSECT=cell(1,ngroups);
CONDS_ELECTIVE=cell(1,ngroups);

%Open file for report
if(CREATE_REPORT)
    completePathToReport=[pathToMainFolder 'Reports\'];
    fileID = fopen(completePathToReport,'w');
end
%% Main cycle
for group=1:ngroups
    
    cgr=groups{group};
    cfold=[folderToData cgr '\'];
    allSubs=dir([cfold '*params.mat']);
    nsubs(group) = length(allSubs);
    
    if PLOT_DATA
        [CONDS_UNION(group), CONDS_INTERSECT(group), CONDS_ELECTIVE(group)] = find_common_conditions(allSubs, cfold);
        plot_all_subjects1(allSubs, cfold, parameters, cgr, ylimits);
        
        %         [subjects_struct_union, subjects_struct_intersect, conds_union, conds_intersect] =...
        %             fill_subject_struct2(allSubs, parameters, cfold);
        
        %         plot_all_conditions(allSubjs, parameters, cfolds);
        %         [subjects_struct, miMaMatrix, conds] = fill_subject_struct(allSubs,parameters,nsubs(group),cfold); %This makes sense only if all the subjects did the same thing
        %         plot_all_subjects(subjects_struct, miMaMatrix, parameters, nsubs(group),cfold);
    end
    
    if CREATE_REPORT
        fprintf(fileID,['Group: ' cgr '\n' 'Common conditions: '   CONDS_INTERSECT(group) '\n' ...
            'Elective conditions: ' CONDS_ELECTIVE(group) '\n\n']);
    end
    
    %     figure
    %     subplot()
    %     for sub=1:nsubs(exp)
    %
    %     end
    
    
    
    
end

OverallCommon   = my_intersect_strings(CONDS_INTERSECT);
% OverallElective = my_intersect_strings(CONDS_INTERSECT);

if CREATE_REPORT
    fprintf(fileID,['Overall: ' cgr '\n' 'Common conditions: '    '\n' ...
        'Elective conditions: \n\n']);
    
    fclose(fileID);
end
%%
% v = 1:10;
% X = [v' 2*v' (v.*v)'];   %10  X 3
% Y = [v'   v'    v'];     % 10 X 3
% figure
% gplotmatrix(X,Y)