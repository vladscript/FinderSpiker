% Function to accumulate Rate of Activity (RoA)
% of each Cell For each Experiment of the same
% experimental setup, i.e., same experimental conditions
% --------------------------------------------------------------------
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
%% Loop to keep loading files
Nconditions=[];
while MoreFiles
    load([PathName,FileName]);
    if exist('MetaDataColocaliation','var')
        aremerged=true;         % Are there already colocated cells
    else
        aremerged=false;    	% Are there already colocated cells
    end
    [~,NC]=size(R_Condition);   % N Conditions
    Nconditions(runs)=NC;       % Nconditions per Experiment
    
    % Start as empty
    if runs==1
        Corrs_ALL=cell(1,NC);
        
    else
        % If Number of Conditoines Changes
        if numel(unique(Nconditions))>numel(Corrs_ALL)
            auxc=auxc+1;
            Corrs_ALL{auxc}=[];
        end
    end
    %%  Loop to Accummulate Data
    
    
    for c=1:NC
        fprintf('\nGetting Data from %s ',Names_Conditions{c})
        AllNamesCond{runs}=Names_Conditions;
        % Read Raster
        R_ALL=R_Condition{c};
        % Correlation distribution
        X=-1*(pdist(R_ALL,'correlation')-1);

        % Read Sizes
        [N_ALL,F_ALL]=size(R_ALL);
        
        % Cummulate Stuff *************************************************
        % Rate of Activity: CaTransients frames / Total Frames
        Trow=[mean(X),median(X),mode(X),var(X),skewness(X),kurtosis(X)];
        TableActive=[TableActive;Trow];
        fprintf('with %i Active Neurons \n',N_ALL);
        % ONLY ACTIVE
        % ALL NEURONS: (better off!) Inclueded [Zero-Rates]
        % Rate of Activity ******************************************
        % Only Actual And Previous
        Corrs_ALL{c}=[Corrs_ALL{c}();X(:)];
        
    end
    % Disp Experiments Selected:
    EXPS{runs,1}=Experiment
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.mat',['Select .mat file, ONE by ONE'],...
    'MultiSelect', 'off',CurrentPathOK);
end
disp('>>end.')
%% Make table of Active Neurons
% expID=cell(numel(EXPS)*numel(Names_Conditions),1);
% AllNamesCond
aux=1;
for n=1:numel(EXPS)
    Names_Conditions=AllNamesCond{n};
    for c=1:numel(Names_Conditions)
        expname=EXPS{n}(EXPS{n}~='\');
        condname=Names_Conditions{c}(isstrprop(Names_Conditions{c},'alphanum'));
        namerow=[expname,'-',condname];
        expID{aux,1}=namerow;
        aux=aux+1;
    end
end
Mean=TableActive(:,1);
Median=TableActive(:,2);
Mode=TableActive(:,3);
Variance=TableActive(:,4);
Skewness=TableActive(:,5);
Kurtosis=TableActive(:,6);
TblAN=table(expID,Mean,Median,Mode,Variance,Skewness, Kurtosis);


%% Save Stuff to .mat File at Processed Data Folder

okbutton = questdlg('Save data @ CSV File?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Correlations_Stats',TS,'.csv'];
    % Select Destiny
    % Directory (default)
    CurrentPath=pwd;
    Slshes=find(CurrentPath=='\');
    % [CurrentPath(1:Slshes(end)),'Raster Features']
    Load_Default_Directories;
    CurrentPathOKac=[CurrentPath(1:Slshes(end)),FolderNameRaster];
    
    PathSave=uigetdir(CurrentPathOKac);
    writetable(TblAN,[PathSave,SaveFile],...
    'Delimiter',',','QuoteStrings',true);
    
    
    fprintf('>> Data saved @: %s\n',PathSave)
else
    fprintf('>>Unsaved data.\n')
end
fprintf('Correlation among neurons Insights: done\n')

%% Save Stuff to .mat File at Processed Data Folder
[CM,ColorIndx]=Color_Selector(Names_Conditions);
plot_cdf_cell(Corrs_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('Correlations')

okbutton = questdlg('Save data?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    disp('>>Saving .mat data...')
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Correlations_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    save([PathSave,SaveFile],'EXPS','Names_Conditions','Corrs_ALL');
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')