% Function to accumulate Rate of Activity (RoA)
% of each Cell For each Experiment of the same
% experimental setup, i.e., same experimental conditions

% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   aremerged:    if merged cells have been identified
%   fs:           Sampling Frequency
% Output varbales @ mat Files
%    RoA_ALL:       Rate of Activity
%    ISI_ALL:       Inter Calcium Transient Interval
%    TransLength:   Calcium Transient Duration
%    RoT_ALL:       Rate Of (calcium) Transients
% If merged cells
%    RoA_POS, ISI, Length & RoT: colocated cells
%    RoA_NEG, ISI, Length & RoT: uncolocated cells
%% Setup
% Initial:
clear; close all; clc;
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
    load([PathName,FileName]);
    if exist('MetaDataColocaliation')
        aremerged=true;    % Are there already colocated cells
    else
        aremerged=false;    % Are there already colocated cells
    end
    [~,NC]=size(R_Condition); % N Conditions
    % Start as empty
    if runs==1
        RoA_ALL=cell(1,NC);
        ISIs_ALL=cell(1,NC);
        TranLengths_ALL=cell(1,NC);
        RoT_ALL=cell(1,NC);
        if aremerged
            RoA_POS=cell(1,NC);
            ISIs_POS=cell(1,NC);
            TranLengths_POS=cell(1,NC);
            RoT_POS=cell(1,NC);
            RoA_NEG=cell(1,NC);
            ISIs_NEG=cell(1,NC);
            TranLengths_NEG=cell(1,NC);
            RoT_NEG=cell(1,NC);
        end
    end
    %%  Loop to Accummulate Data
    for c=1:NC
        fprintf('\nGetting Data from %s ',Names_Conditions{c})
        % Read Raster
        R_ALL=R_Condition{c};
        if aremerged
            R_POS=R_merged{c};
            R_NEG=R_nomerged{c};
        end
        % Read Sizes
        [N_ALL,F_ALL]=size(R_ALL);
        
        if aremerged
            [N_POS,F_POS]=size(R_POS);
            [N_NEG,F_NEG]=size(R_NEG);
            AllActive=[MetaDataColocaliation.PositiveCells;MetaDataColocaliation.NegativeCells];
            PosActive=MetaDataColocaliation.PositiveCells;
            NegActive=MetaDataColocaliation.NegativeCells;
        end
        % Cummulate Stuff *************************************************
        % Rate of Activity: CaTransients frames / Total Frames
        AllActive=find(sum(R_ALL,2));
        fprintf('with %i of %i Active Neurons \n',numel(AllActive),numel(sum(R_ALL,2)));
        % ONLY ACTIVE
        % RoA_ALL{c}=[RoA_ALL{c}();sum(R_ALL(AllActive,:),2)./F_ALL];
        % ALL NEURONS: (better off!) Inclueded [Zero-Rates]
        RoA_ALL{c}=[RoA_ALL{c}();sum(R_ALL,2)./F_ALL];
        if aremerged
            RoA_POS{c}=[RoA_POS{c};sum(R_merged{c},2)./F_POS];
            RoA_NEG{c}=[RoA_NEG{c};sum(R_nomerged{c},2)./F_NEG];        
        end
        % Inter Calcium Transient Interval & Calcium Transient Duration
        ISIs=[];
        TranLengths=[];
        RoTs=[];
        for i=1:N_ALL
            r=R_ALL(i,:);
            [ISIcell,TranLengthscell]=interval_duration_events(r);
            % row vectors
            ISIs=[ISIs,ISIcell]; 
            TranLengths=[TranLengths,TranLengthscell];
            RoT(i)=numel(TranLengthscell)*fs/numel(r)*60;
            % RoTs=[RoTs,RoT];
        end
        ISIs_ALL{c}=[ISIs_ALL{c};ISIs'/fs];
        TranLengths_ALL{c}=[TranLengths_ALL{c};TranLengths'/fs];
        RoT_ALL{c}=[RoT_ALL{c};RoT'];
        % + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
        % If there are merged colocated cells + + + + + + + + + + + + + 
        if aremerged
            % Positive Colocated
            ISIsPOS=[];
            TranLengthsPOS=[];
            RoTPOS=[];
            for i=1:N_POS
                r=R_POS(i,:);
                [ISIcell,TranLengthscell]=interval_duration_events(r);
                % row vectors
                ISIsPOS=[ISIsPOS,ISIcell]; 
                TranLengthsPOS=[TranLengthsPOS,TranLengthscell];
                RoTPOS(i)=numel(TranLengthscell)*fs/numel(r)*60;
            end
            ISIs_POS{c}=[ISIs_POS{c};ISIsPOS'/fs];
            TranLengths_POS{c}=[TranLengths_POS{c};TranLengthsPOS'/fs];
            RoT_POS{c}=[RoT_POS{c};RoTPOS'];
            % Negative
            ISIsNEG=[];
            TranLengthsNEG=[];
            RoTNEG=[];
            for i=1:N_NEG
                r=R_NEG(i,:);
                [ISIcell,TranLengthscell]=interval_duration_events(r);
                % row vectors
                ISIsNEG=[ISIsNEG,ISIcell]; 
                TranLengthsNEG=[TranLengthsNEG,TranLengthscell];
                RoTNEG(i)=numel(TranLengthscell)*fs/numel(r)*60;
            end
            ISIs_NEG{c}=[ISIs_NEG{c};ISIsNEG'/fs];
            TranLengths_NEG{c}=[TranLengths_NEG{c};TranLengthsNEG'/fs];
            RoT_NEG{c}=[RoT_NEG{c};RoTNEG'];
        end
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
Plot_Accumulate_CDF;
okbutton = questdlg('Save data?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\RoA_RoT_ISI_SD_ACC_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    if ~aremerged
        disp('>>saving data...')
        save([PathSave,SaveFile],'EXPS','Names_Conditions','RoA_ALL',...
            'ISIs_ALL','TranLengths_ALL','RoT_ALL');
    else
        disp('>>saving data...')
        save([PathSave,SaveFile],'EXPS','Names_Conditions','RoA_ALL',...
            'ISIs_ALL','TranLengths_ALL','RoT_ALL',...
            'RoA_POS','ISIs_POS','TranLengths_POS','RoT_POS',...
            'RoA_NEG','ISIs_NEG','TranLengths_NEG','RoT_NEG');
    end
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')