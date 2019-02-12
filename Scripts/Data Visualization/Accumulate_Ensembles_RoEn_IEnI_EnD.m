% Script to accumulate Rate of Ensmbles(RoE)
% of each Ensmble For each Experiment of the same
% experimental setup, i.e., same experimental conditions

% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   fs: Sampling Frequency
% Output
%    RoE_ALL:       Rate of Ensembles
%    IEnI:          Inter Ensembles Interval
%    EnD:           Ensembles Duration

%% Setup
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
            IEnI=[IEnI;Features_Ensemble.IEIsExp{c,e}'];
            EnD=[EnD;Features_Ensemble.EDsExp{c,e}'];
        end
        RoE_ALL{c}=[RoE_ALL{c};RateEnsembles];
        IEnIs_ALL{c}=[IEnIs_ALL{c};IEnI];
        EnD_ALL{c}=[EnD_ALL{c};EnD];
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
plot_cdf_cell(RoE_ALL,Names_Conditions);
plot_cdf_cell(IEnIs_ALL,Names_Conditions);
plot_cdf_cell(EnD_ALL,Names_Conditions);
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
        'IEnIs_ALL','EnD_ALL');
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')