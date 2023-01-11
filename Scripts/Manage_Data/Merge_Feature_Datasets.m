% Script to MERGE Datasets of the Following Caterogry of Features:
%   'Raster Activity':      Nerual Activity
%   'General Ensembles':    Neural Ensembles Features
%   'Detailed Ensembles':   Single Neural Ensembles Features
%% Intro
fprintf('\n\nGet CSV Tables from Analsysis Features.\n')
fprintf('Table_Raster_Features_\n Ensemble_Features_\n Detailed_Ensemble_Features_\n Node_StatsFeature_Dataset_\n')
%% Select Kind Of Features to Merge
Load_Default_Directories;
KindFeatures={'Raster Activity';'General Ensembles';...
    'Detailed Ensembles';'Network Features'};
[index_var,index_CHECK] = listdlg('PromptString','Select Sort of Features:',...
            'SelectionMode','single',...
            'ListString',KindFeatures);
SelectedCatergory=KindFeatures{index_var};
%% Setup
% Initial:
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments
% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
% Select Folder
if index_var<2
    FolderDefault=FolderNameRaster;
    HeadersFeatures={'Condition','EXP_ID','RateNeurons','ActivityTimeFraction',...
        'ActiveRatioCAG','EffectiveActivity',...
        'ISImean','ISImode','ISImedian','ISIvar','ISIskew','ISIkurt',...
        'Lengthmean','Lengthmode','Lengthmedian','Lengthvar','Lengthskew','Lengthkurt',...
        'CAGmean','CAGmode','CAGmedian','CAGvar','CAGskew','CAGkurt',...
        'RoAmean','RoAmode','RoAmedian','RoAvar','RoAskew','RoAkurt',...
        'RoTmean','RoTmode','RoTmedian','RoTvar','RoTskew','RoTkurt'};
else
    FolderDefault=FolderNameEnsembles;
    if index_var<3
        HeadersFeatures={'Condition','EXP_ID','Nensembles','Threshold','Dunns','MaxIntraVec','ClassError',...
         'RateTrans','RateCycles','SimpleCycles','ClosedCycles','OpenedCycles',...
         'CAGauc','CoreSize','MaxSynLinks','MaxConn_A','MaxConn_B',...
         'EDMean','EDMode','EDMedian','EDVar','EDSkew','EDKurt',...
         'IEIMean','IEIMode','IEIMedian','IEIVar','IEISkew','IEIKurt',...
         'SynWeigthMean','SynWeigthMode','SynWeightMedian','SynWeigthVar','SynWeigthSkew','SynWeigthKurt',...
         'AlternativeIndx','RecurrentIndx'};
    elseif index_var<4
        HeadersFeatures={'Condition','EXP_ID','NeuronsRation','Time','Dominance','Rate',...
     'ensCAGauc','ensCAGmean','ensCAGmode','ensCAGmedian','ensCAGvar','ensCAGskew','ensCAGkurt'...
     'IEImean','IEImode','IEImedian','IEIvar','IEIskew','IEIkurt'...
     'EDmean','EDmode','EDmedian','EDvar','EDskew','EDkurt'};
    else
        FolderDefault=[];
        % HeadersFeatures to be READ
    end
end
CurrentPathOK=[Dirpwd(1:slashesindx(end)),FolderDefault]; 
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.csv',[SelectedCatergory,' Feature Table CSV file ONE:by:ONE'],...
    'MultiSelect', 'off',CurrentPathOK);
%% Loop to keep loading files
FeaturesSingle=table;
while MoreFiles

    % Numerical Data
    rowFeatures=readtable([PathName,FileName]);
    FeaturesSingle=[FeaturesSingle;rowFeatures]; 
    
    % Disp Experiments Selected:
    EXPS{runs,1}=FileName
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.csv',[SelectedCatergory,' Feature Table CSV file ONE:by:ONE'],...
    'MultiSelect', 'off',CurrentPathOK);
end
if index_var==4
    HeadersFeatures=FeaturesSingle.Properties.VariableNames;
end
FeaturesSingle.Properties.VariableNames=HeadersFeatures;
disp('>>end.')
%% Make and Save Table
okbutton = questdlg(['Make CSV for ',SelectedCatergory,' Features Table?']);
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SelectedCatergory(SelectedCatergory==' ')='_';
    SaveFile=[SelectedCatergory,'_Dataset_',TS,'.csv'];
    %
    CurrentPathOK=[Dirpwd(1:slashesindx(end)),FolderNameDataset];
    if exist(CurrentPathOK,'dir')==0
        % create
        fprintf('>Creating database folder:')
        mkdir(CurrentPathOK);
        fprintf('done\n')
    end
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK,['Set Directory to save Dataset of: ',SelectedCatergory]);
    disp('>>Making CSV table...')
    writetable(FeaturesSingle,[PathSave,'\',SaveFile],...
                    'Delimiter',',','QuoteStrings',true);
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')
%% END