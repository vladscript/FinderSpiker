% Function to accumulate Rate of Activity (RoA)
% of each Cell For each Experiment of the same
% experimental setup, i.e., same experimental conditions

% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   aremerged:    if merged cells have been identified
%   fs: Sampling Frequencu
% Output
%    RoA_ALL:   Rate of Activity
%    ICaTI:     Inter Calcium Transient Interval
%    CaTD:      Calcium Transient Duration
% If merged cells
%    RoA_POS: colocated cells
%    RoA_NEG: uncolocated cells
%% Setup
% Initial:
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments
aremerged=false;    % Are there already colocated cells
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
    [~,NC]=size(R_Condition); % N Conditions
    % Start as empty
    if runs==1
        RoA_ALL=cell(1,NC);
        ISIs_ALL=cell(1,NC);
        TranLengths_ALL=cell(1,NC);
        if aremerged
            RoA_POS=cell(1,NC);
            ISIs_POS=cell(1,NC);
            TranLengths_POS=cell(1,NC);
            RoA_NEG=cell(1,NC);
            ISIs_NEG=cell(1,NC);
            TranLengths_NEG=cell(1,NC);
        end
    end
    %%  Loop to Accummulate Data
    for c=1:NC
        fprintf('Getting Data from %s ',Names_Conditions{c})
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
        fprintf('with %i cells\n',numel(AllActive))
        RoA_ALL{c}=[RoA_ALL{c}();sum(R_ALL(AllActive,:),2)./F_ALL];
        if aremerged
            RoA_POS{c}=[RoA_POS{c};sum(R_merged{c},2)./F_POS];
            RoA_NEG{c}=[RoA_NEG{c};sum(R_nomerged{c},2)./F_NEG];        
        end
        % Inter Calcium Transient Interval & Calcium Transient Duration
        ISIs=[];
        TranLengths=[];
        for i=1:N_ALL
            r=R_ALL(i,:);
            [ISIcell,TranLengthscell]=interval_duration_events(r);
            % row vectors
            ISIs=[ISIs,ISIcell]; 
            TranLengths=[TranLengths,TranLengthscell];
        end
        ISIs_ALL{c}=[ISIs_ALL{c};ISIs'/fs];
        TranLengths_ALL{c}=[TranLengths_ALL{c};TranLengths'/fs];
        % + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
        % If there are merged colocated cells + + + + + + + + + + + + + 
        if aremerged
            % Positive Colocated
            ISIsPOS=[];
            TranLengthsPOS=[];
            for i=1:N_POS
                r=R_POS(i,:);
                [ISIcell,TranLengthscell]=interval_duration_events(r);
                % row vectors
                ISIsPOS=[ISIsPOS,ISIcell]; 
                TranLengthsPOS=[TranLengthsPOS,TranLengthscell];
            end
            ISIs_POS{c}=[ISIs_POS{c};ISIsPOS'/fs];
            TranLengths_POS{c}=[TranLengths_POS{c};TranLengthsPOS'/fs];
            % Negative
            ISIsNEG=[];
            TranLengthsNEG=[];
            for i=1:N_NEG
                r=R_NEG(i,:);
                [ISIcell,TranLengthscell]=interval_duration_events(r);
                % row vectors
                ISIsNEG=[ISIsNEG,ISIcell]; 
                TranLengthsNEG=[TranLengthsNEG,TranLengthscell];
            end
            ISIs_NEG{c}=[ISIs_NEG{c};ISIsNEG'/fs];
            TranLengths_NEG{c}=[TranLengths_NEG{c};TranLengthsNEG'/fs];
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
plot_cdf_cell(RoA_ALL,Names_Conditions);
plot_cdf_cell(ISIs_ALL,Names_Conditions);
plot_cdf_cell(TranLengths_ALL,Names_Conditions);
okbutton = questdlg('Save data?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\RoA_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    if ~aremerged
        disp('>>saving data...')
        save([PathSave,SaveFile],'EXPS','Names_Conditions','RoA_ALL',...
            'ISIs_ALL','TranLengths_ALL');
    else
        disp('>>saving data...')
        save([PathSave,SaveFile],'EXPS','Names_Conditions','RoA_ALL',...
            'ISIs_ALL','TranLengths_ALL',...
            'RoA_POS','RoA_NEG');
    end
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')