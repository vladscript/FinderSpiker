% Function to accumulate Active Neurons in Venn Diagrams
% of each Cell For each Experiment of the same
% experimental setup, i.e., same experimental conditions
% --------------------------------------------------------------------
% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   aremerged:    if merged cells have been identified
% 
% Output varbales @ mat Files
%    Set of Cells:       Activa in PAIR of Conditions

% If merged cells
%     buiding 
% 
%% Setup
% Initial:
clear; clc;
runs=1;             % Runs Counter
auxc=1;             % Auxiliar Conditions
EXPS={};            % List Of Experiments

% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'Processed Data'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.mat',['Select .mat file, ONE by ONE'],...
    'MultiSelect', 'off',CurrentPathOK);
% Table Active Cells: TOTAL/ACTIVE/POSITIVE/NEGATIVE
Trow=zeros(1,4);
TableActive=[];
%% Setup
Nconditions=[];
NcellsOffset=0;
% SETS
GlobalCells.Ninh=[];
GlobalCells.Ndep=[];
GlobalCells.Nfac=[];
GlobalCells.Nact=[];
GlobalCells.Nunc=[];
% ACtivity
GlobalRoA.inh=[];
GlobalRoA.dep=[];
GlobalRoA.fac=[];
GlobalRoA.act=[];
GlobalRoA.unc=[];
% N cells
NcellsAct.inh=[];
NcellsAct.dep=[];
NcellsAct.fac=[];
NcellsAct.act=[];
NcellsAct.unc=[];

% Active Cells in Each Condition
GlobalActiveCells.CellsIndx=[];
GlobalActiveCells.AC_1=[];
GlobalActiveCells.AC_2=[];
%%  Loop to Accummulate Data
while MoreFiles
    load([PathName,FileName]);
    % CHECK IF MARKED CELLS (A2a+, D1+, ChAT+)
    if exist('MetaDataColocaliation','var')
        aremerged=true;         % Are there already colocated cells
    else
        aremerged=false;    	% Are there already colocated cells
    end
    [~,NC]=size(R_Condition);   % N Conditions
    Nconditions(runs)=NC;       % Nconditions per Experiment
    
    % Start as empty
    if runs==1
        CellIndexes=cell(1,NC);
        if aremerged
            MarkedCellndexes=cell(1,NC);
        end
    else
        % If Number of Conditiones Changes
        if numel(unique(Nconditions))>numel(R_Condition)
            auxc=auxc+1;
            CellIndexes{auxc}=[];
            if aremerged
                MarkedCellndexes{auxc}=[];
            end
        end
    end
    % ACCUMULATE INTEL
    % Select PAIR of CONDITIONS 
    [ActivityGroups, ActiveCells, TwoNames, RoAG]= vennRASTER(R_Condition,Names_Conditions);
    % gradient-active cell groups
    GlobalCells.Ninh=[GlobalCells.Ninh; ActivityGroups.Ninh+NcellsOffset];
    GlobalCells.Ndep=[GlobalCells.Ndep; ActivityGroups.Ndep+NcellsOffset];
    GlobalCells.Nfac=[GlobalCells.Nfac; ActivityGroups.Nfac+NcellsOffset];
    GlobalCells.Nact=[GlobalCells.Nact; ActivityGroups.Nact+NcellsOffset];
    GlobalCells.Nunc=[GlobalCells.Nunc; ActivityGroups.Nunc+NcellsOffset];
    % gradient-active cell activity 
    GlobalRoA.inh=[GlobalRoA.inh; RoAG.inh];
    GlobalRoA.dep=[GlobalRoA.dep; RoAG.dep];
    GlobalRoA.fac=[GlobalRoA.fac; RoAG.fac];
    GlobalRoA.act=[GlobalRoA.act; RoAG.act];
    GlobalRoA.unc=[GlobalRoA.unc; RoAG.unc];
    
    % N cells activity per Group
    NcellsAct.inh(runs)=numel(ActivityGroups.Ninh);
    NcellsAct.dep(runs)=numel(ActivityGroups.Ndep);
    NcellsAct.fac(runs)=numel(ActivityGroups.Nfac);
    NcellsAct.act(runs)=numel(ActivityGroups.Nact);
    NcellsAct.unc(runs)=numel(ActivityGroups.Nunc);
    
    % active cell groups
    GlobalActiveCells.CellsIndx=[GlobalActiveCells.CellsIndx, ActiveCells.CellsIndx+NcellsOffset];
    GlobalActiveCells.AC_1=[GlobalActiveCells.AC_1, ActiveCells.AC_1+NcellsOffset];
    GlobalActiveCells.AC_2=[GlobalActiveCells.AC_2, ActiveCells.AC_2+NcellsOffset];
    % Increase INDEX CELLS
    NcellsOffset=NcellsOffset+max(ActiveCells.CellsIndx);
    
    if max(ActiveCells.CellsIndx)==size(R_Condition{1},1)
        fprintf('>>Index: [OK]\n')
    else
        fprintf('>>[WARNING]\n')
    end
    
    % Displas LIST of  Selected Experiments :
    EXPS{runs,1}=Experiment
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.mat',['Select .mat file, ONE by ONE'],...
    'MultiSelect', 'off',CurrentPathOK);
end
disp('>>end.')
%% PLOT VENNS
% By Changed-Activity Groups
plotvennActivity(GlobalCells)
%% PLOT ALL CELLs ACTIVITY
plot_boxplots_cellgroups(GlobalRoA);

%% Bar graph
figure
Axis=subplot(1,1,1);
Axis.XLim=[0,4];
WidthBar=1;
Ncells=max(GlobalActiveCells.CellsIndx);
bar(1,100*numel(GlobalCells.Nfac)/Ncells,WidthBar,'g');  hold on;
bar(1,100*numel(GlobalCells.Nact)/Ncells,WidthBar,'g');

bar(2,100*numel(GlobalCells.Nunc)/Ncells,WidthBar,'b');

bar(3,100*numel(GlobalCells.Ndep)/Ncells,WidthBar,'r');
bar(3,100*numel(GlobalCells.Ninh)/Ncells,WidthBar,'r'); hold off;

Axis.YLim=[0,100];
Axis.XTick=[1,2,3];
Axis.XTickLabel={'Facilitated/Activated';'Unchanged';'Depressed/Inactivated'};
Axis.YLabel.String='% Recorded Cells';
title(sprintf('Change of activity between %s and %s',TwoNames{1},TwoNames{2}))


%% By Changed-Activity Groups
plotvennCellActive(GlobalActiveCells,TwoNames);
title(sprintf('Sets of Active Neurons in %s and %s',TwoNames{1},TwoNames{2}))
%% Save Data
okbutton = questdlg('Save data in .mat file?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    disp('>>Saving .mat data...')
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Cell_Activity_Groups_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    if ~aremerged
        save([PathSave,SaveFile],'EXPS','TwoNames','GlobalCells',...
            'GlobalRoA','NcellsAct','GlobalActiveCells');
    else
%         save([PathSave,SaveFile],'EXPS','Names_Conditions','RoA_ALL',...
%             'ISIs_ALL','TranLengths_ALL','RoT_ALL',...
%             'RoA_POS','ISIs_POS','TranLengths_POS','RoT_POS',...
%             'RoA_NEG','ISIs_NEG','TranLengths_NEG','RoT_NEG');
    end
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
% clear
fprintf('done\n')
%% END: clear all