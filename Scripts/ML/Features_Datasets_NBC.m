% Script To Evaluate Which Category of Features
% Classify the Best the Measured Conditions
%% Setup
%% Load Dataset
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % Finder Spiker Main Folder
% Load File 
[FileName,PathName,MoreFiles] = uigetfile({'*.csv';'*.xlsx'},' Dataset file',...
    'MultiSelect', 'off',CurrentPathOK);
if strcmp(FileName(end-3:end),'xlsx')
    disp('>>Network Features')
    Xraw=readtable([PathName,FileName]);        % Table
    % Check table
    LimitRow=find(isnan(table2array (Xraw(:,3) ) ),1)-1;
    if isempty(LimitRow)
        LimitRow=numel(Xraw(:,3));
    end
    LimitCol=12;
    Xraw=Xraw(1:LimitRow,1:LimitCol);
    Y=categorical(table2array( Xraw(:,1)) );    % Labels
    EXPIDs=table2array( Xraw(:,2));             % Cell of Strings
    X=table2array( Xraw(:,3:end) );             % Dataset
else
    Xraw=readtable([PathName,FileName]);        % Table
    Y=categorical(table2array( Xraw(:,1)) );    % Labels
    EXPIDs=table2array( Xraw(:,2));             % Cell of Strings
    X=table2array( Xraw(:,3:end) );             % Dataset
end
%% Get Intel of the Feature Category
% Category
SpacesIndx=find(FileName=='_');
CategoryFeature=FileName(1:SpacesIndx(2)-1);
% Number of Features
[DatasetSize,NFeatures]=size(X);
NameFeatures=Xraw.Properties.VariableNames(3:end);
% Ignore Some Column Non-Features at Ensemble Features
if  strcmp(CategoryFeature,'General_Ensembles')
    % Ignore: 15,16 & 17 -th Feature
    ignftr=[13,14,15]; % Links/Neuron A->euron B Most Connected
    acceptcols=setdiff(1:numel(NameFeatures),ignftr);
    X=X(:,acceptcols);
    NameFeatures=NameFeatures(acceptcols);
end
%% GET BEST SUBSET from ALL FEATURES ######################################
OKFeatures = Best_Subset_Feature(Y,X,NameFeatures,EXPIDs);

