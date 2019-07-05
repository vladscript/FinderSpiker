% Script to MERGE Datasets of the Following Caterogry of Features:
%   'Raster Activity':      Nerual Activity
%   'General Ensembles':    Neural Ensembles Features
%   'Detailed Ensembles':   Single Neural Ensembles Features
%% Select Kind Of Features to Merge
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
    FolderDefault='Raster Features';
    HeadersFeatures={'Condition','EXP_ID','RateNeurons','ActivityTimeFraction',...
        'ActiveRatioCAG','EffectiveActivity',...
        'ISImean','ISImode','ISImedian','ISIvar','ISIskew','ISIkurt',...
        'Lengthmean','Lengthmode','Lengthmedian','Lengthvar','Lengthskew','Lengthkurt',...
        'CAGmean','CAGmode','CAGmedian','CAGvar','CAGskew','CAGkurt',...
        'RoAmean','RoAmode','RoAmedian','RoAvar','RoAskew','RoAkurt',...
        'RoTmean','RoTmode','RoTmedian','RoTvar','RoTskew','RoTkurt'};
else
    FolderDefault='Ensemble Features';
    if index_var<3
        HeadersFeatures={'Condition','EXP_ID','Nensembles','Threshold','Dunns','MaxIntraVec','ClassError',...
         'RateTrans','RateCycles','SimpleCycles','ClosedCycles','OpenedCycles',...
         'CAGauc','CoreSize','MaxSynLinks','MaxConn_A','MaxConn_B',...
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
[FileName,PathName,MoreFiles] = uigetfile('*.csv',[SelectedCatergory,' Feature Database file'],...
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
    [FileName,PathName,MoreFiles] = uigetfile('*.csv',[SelectedCatergory,' Feature Database file'],...
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
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK,['Set Directory to save Dataset of :',SelectedCatergory]);
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