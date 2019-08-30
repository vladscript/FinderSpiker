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
clear; clc;
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments

% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'Processed Data'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
% Table Active Cells: TOTAL/ACTIVE/POSITIVE/NEGATIVE
Trow=zeros(1,4);
TableActive=[];
%% Loop to keep loading files
while MoreFiles
    load([PathName,FileName]);
    if exist('MetaDataColocaliation','var')
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
    RasterDurations=get_raster_durations(Onsets,R_Condition,fs);
    for c=1:NC
        fprintf('\nGetting Data from %s ',Names_Conditions{c})
        % Read Raster
        R_ALL=R_Condition{c};
        % Read Denoised Signal:
        [CleanSignals_ALL,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,ESTSIGNALS,XY_selected,XY);
        
        if aremerged
            % Raster of Positive Cells
            R_POS=R_merged{c};
            % Get Coordinates
            XY_merged=XY_selected(MetaDataColocaliation.PositiveCells,:);
            % Indexes of the Original Denoised Signals and Vectors
            [CleanSignals_POS,~]=Retrieve_Selected_Signal(Onsets,R_merged,ESTSIGNALS,XY_merged,XY);
            % Raster of Positive Cells
            R_NEG=R_nomerged{c};
            XY_nomerged=XY_selected(MetaDataColocaliation.NegativeCells,:);
            % Indexes of the Original Denoised Signals and Vectors
            [CleanSignals_NEG,~]=Retrieve_Selected_Signal(Onsets,R_nomerged,ESTSIGNALS,XY_nomerged,XY);
        end
        % Read Sizes
        [N_ALL,F_ALL]=size(R_ALL);
        
        if aremerged
            [N_POS,F_POS]=size(R_POS);
            [N_NEG,F_NEG]=size(R_NEG);
            AllActive=[MetaDataColocaliation.PositiveCells;MetaDataColocaliation.NegativeCells];
            PosActive=MetaDataColocaliation.PositiveCells;
            NegActive=MetaDataColocaliation.NegativeCells;
        else
            AllActive=find(sum(R_ALL,2));
            PosActive=[];
            NegActive=[];
        end
        % Cummulate Stuff *************************************************
        % Rate of Activity: CaTransients frames / Total Frames
        Trow=[numel(sum(R_ALL,2)),numel(AllActive),numel(PosActive),numel(PosActive)];
        TableActive=[TableActive;Trow];
        fprintf('with %i of %i Active Neurons \n',numel(AllActive),numel(sum(R_ALL,2)));
        % ONLY ACTIVE
        % ALL NEURONS: (better off!) Inclueded [Zero-Rates]
        % Rate of Activity ******************************************
        RoA_ALL{c}=[RoA_ALL{c}();sum(R_ALL,2)./F_ALL];
        if aremerged
            RoA_POS{c}=[RoA_POS{c};sum(R_merged{c},2)./F_POS];
            RoA_NEG{c}=[RoA_NEG{c};sum(R_nomerged{c},2)./F_NEG];        
        end
        % Inter Calcium Transient Interval & Calcium Transient Duration
        [NoT,ITI,LT]=get_RoT(R_ALL,CleanSignals_ALL{c});
        ISIs=ITI/fs;                            % [seconds]
        TranLengths=LT/fs;                      % [seconds]
        RoTs=NoT/(size(R_Condition{c},2)/fs/60);  % [minutes]
        % Accumulate
        ISIs_ALL{c}=[ISIs_ALL{c};ISIs];
        TranLengths_ALL{c}=[TranLengths_ALL{c};TranLengths];
        RoT_ALL{c}=[RoT_ALL{c};RoTs];
        % + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
        % If there are merged colocated cells + + + + + + + + + + + + + 
        if aremerged
            % Positive Colocated
            [NoT,ITI,LT]=get_RoT(R_merged{c},CleanSignals_POS{c});
            ISIsPOS=ITI/fs;                         % [seconds]
            TranLengthsPOS=LT/fs;                   % [seconds]
            RoTPOS=NoT/(size(R_merged{c},2)/fs/60);   % [minutes]
            % Accumulate Postitive
            ISIs_POS{c}=[ISIs_POS{c};ISIsPOS];
            TranLengths_POS{c}=[TranLengths_POS{c};TranLengthsPOS];
            RoT_POS{c}=[RoT_POS{c};RoTPOS];
            % Negative
            [NoT,ITI,LT]=get_RoT(R_nomerged{c},CleanSignals_NEG{c});
            ISIsNEG=ITI/fs;                         % [seconds]
            TranLengthsNEG=LT/fs;                   % [seconds]
            RoTNEG=NoT/(size(R_nomerged{c},2)/fs/60); % [minutes]
            % Accumulate
            ISIs_NEG{c}=[ISIs_NEG{c};ISIsNEG];
            TranLengths_NEG{c}=[TranLengths_NEG{c};TranLengthsNEG];
            RoT_NEG{c}=[RoT_NEG{c};RoTNEG];
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
%% Make table of Active Neurons
expID=cell(numel(EXPS)*numel(Names_Conditions),1);
aux=1;
for n=1:numel(EXPS)
    for c=1:numel(Names_Conditions)
        expname=EXPS{n}(EXPS{n}~='\');
        condname=Names_Conditions{c}(isstrprop(Names_Conditions{c},'alphanum'));
        namerow=[expname,'-',condname];
        expID{aux}=namerow;
        aux=aux+1;
    end
end
Total=TableActive(:,1);
Active=TableActive(:,2);
Positive=TableActive(:,3);
Negative=TableActive(:,4);
TblAN=table(expID,Total,Active,Positive,Negative)

%% Save Stuff to .mat File at Processed Data Folder

okbutton = questdlg('Save data @ CSV File?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Active_Cells',TS,'.csv'];
    % Select Destiny
    % Directory (default)
    CurrentPath=pwd;
    Slshes=find(CurrentPath=='\');
    % [CurrentPath(1:Slshes(end)),'Raster Features']
    CurrentPathOKac=[CurrentPath(1:Slshes(end)),'Raster Features'];
    
    PathSave=uigetdir(CurrentPathOKac);
    writetable(TblAN,[PathSave,SaveFile],...
    'Delimiter',',','QuoteStrings',true);
    
    
    fprintf('>> Data saved @: %s\n',PathSave)
else
    fprintf('>>Unsaved data.\n')
end
fprintf('Active Neurons Insights: done\n')

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