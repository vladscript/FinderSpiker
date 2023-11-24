% Script to accumulate Rate of Ensembles(RoE)
% of each Ensmble For each Experiment of the same
% experimental setup, i.e., same experimental conditions
% 
% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   fs: Sampling Frequency
% Output
%    RoE_ALL:       Rate of Ensembles
%    IEnI:          Inter Ensembles Interval
%    EnD:           Ensembles Duration
%    Number of Ensembles per Cycle
%% Setup
clear; clc;
% Initial:
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments
% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'Processed Data'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
%% Loop to keep loading files
while MoreFiles
    load([PathName,FileName])
    Nensembles=Features_Condition.Nenscond; % N ensembles in each condition
    NC=numel(Nensembles); % N Conditions
    % Start as empty
    if runs==1
        RoE_ALL=cell(1,NC);
        IEnIs_ALL=cell(1,NC);
        EnD_ALL=cell(1,NC);
        NensCyc_ALL=cell(3,NC);
        LatencyCyc_ALL=cell(3,NC);
    end
    %%  Loop to Accummulate Data
    for c=1:NC
        % Read Ensmebles Data
        RateEnsembles=[];
        IEnI=[];
        EnD=[];
        for e=1:Nensembles(c)
            fprintf('Getting Data from %s Ensemble %i\n',Names_Conditions{c},e)
            RateEnsembles=[RateEnsembles;Features_Ensemble.Rate(c,e)];
            IEnI=[IEnI;makerowvector(Features_Ensemble.IEIsExp{c,e})'];
            EnD=[EnD;makerowvector(Features_Ensemble.EDsExp{c,e})'];
        end
        RoE_ALL{c}=[RoE_ALL{c};RateEnsembles];
        IEnIs_ALL{c}=[IEnIs_ALL{c};IEnI];
        EnD_ALL{c}=[EnD_ALL{c};EnD];
        % Sort of Cycles
        if ~isempty(Features_Condition.CyclesTable{c})
            CycleTable=Features_Condition.CyclesTable{c};
            Ncycles=size(Features_Condition.CyclesTable{c},1);
            for n=1:Ncycles
                CycleType=CycleTable{n,1};
                switch CycleType
                case 'Simple'
                    NensCyc_ALL{1,c}=[NensCyc_ALL{1,c},numel(CycleTable{n,2})];
                    LatencyCyc_ALL{1,c}=[LatencyCyc_ALL{1,c},CycleTable{n,4}-CycleTable{n,3}];
                case 'Closed'
                    NensCyc_ALL{2,c}=[NensCyc_ALL{2,c},numel(CycleTable{n,2})];
                    LatencyCyc_ALL{2,c}=[LatencyCyc_ALL{2,c},CycleTable{n,4}-CycleTable{n,3}];
                case 'Open'
                    NensCyc_ALL{3,c}=[NensCyc_ALL{3,c},numel(CycleTable{n,2})];
                    LatencyCyc_ALL{3,c}=[LatencyCyc_ALL{3,c},CycleTable{n,4}-CycleTable{n,3}];
                end
            end
        end
        % + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    end
    % Disp Experiments Selected:
    EXPS{runs,1}=Experiment
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
end
disp('>>end.')
%% Save Stuff to .mat File at Processed Data Folder
%  Color Selector:
[CM,ColorIndx]=Color_Selector(Names_Conditions);
% All the plots
plot_cdf_cell(RoE_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('Rate of Ensembles [E/min]')
plot_cdf_cell(IEnIs_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('Inter Ensembles Interval [s]')
plot_cdf_cell(EnD_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('Ensembles Duration [s]')
plot_cdf_cell(NensCyc_ALL(1,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Number of Ensmebles per Simple Cycle [Ensembles]')
plot_cdf_cell(NensCyc_ALL(2,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Number of Ensmebles per Closed Cycle [Ensembles]')
plot_cdf_cell(NensCyc_ALL(3,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Number of Ensmebles per Open Cycle [Ensembles]')
plot_cdf_cell(LatencyCyc_ALL(1,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Latency of Simple Cycle Ensmbles [s]')
plot_cdf_cell(LatencyCyc_ALL(2,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Latency of Closed Cycle Ensmbles [s]')
plot_cdf_cell(LatencyCyc_ALL(3,:),Names_Conditions,CM(ColorIndx,:));
xlabel('Latency of Open Cycle Ensmbles  [s]')
okbutton = questdlg('Save data?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Rate_Ensembles_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    disp('>>saving data...')
    save([PathSave,SaveFile],'EXPS','Names_Conditions','RoE_ALL',...
        'IEnIs_ALL','EnD_ALL','NensCyc_ALL','LatencyCyc_ALL');
    %fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
    fprintf('<a href="matlab:dos(''explorer.exe /e, %s, &'')">> Data saved here</a>\n',PathSave);
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')