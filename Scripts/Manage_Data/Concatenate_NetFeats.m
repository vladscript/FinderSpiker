%% Script to concatenate Network Features: NetFeat Files
% Makes a TABLE of Different and MULTIPLE Network Features
%% Setup
% Initial:
clear; close all; clc;
runs=1;             % Runs Counter
FEATS={};            % List Of Features

% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'NetWorks-CSV'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.csv',' Network Feature Dataset .csv file',...
    'MultiSelect', 'off',CurrentPathOK);
%% COLLECT DATA #####################################################

deliminter=',';
auxc=1;
Ytotal=table; % All Network Features of ALl Experiments

while MoreFiles
    % GET TABLE:
    X=readtable([PathName,FileName]);
    Xnames=X.Properties.VariableNames;
    Xdata=X(:,3:end); % 
    Conds=X(:,1);
    ListExps=X(:,2);
    if auxc==1
        Ytotal=[Ytotal,Conds,ListExps,Xdata];
    else
        Ytotal=[Ytotal,Xdata];
    end
    % Table Intel:
    SpaceIndx=find(FileName=='_');
    FEATPID=FileName(SpaceIndx(1)+1:SpaceIndx(2)-1)
    
    % Disp Datasets Selected:
    FEATS{runs,1}=FileName
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.csv',' Network Feature Dataset .csv file',...
    'MultiSelect', 'off',CurrentPathOK);
    auxc=auxc+1;
end

disp('>>Loading Features Network Stats: Done.')
%% SAVE DATASET
disp('>>Saving Network Features ...')

okbutton = questdlg('Make CSV Table?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\NetFEATS_Dataset_',TS,'.csv'];
    CurrentPathOK=[Dirpwd(1:slashesindx(end))];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    writetable(Ytotal,[PathSave,SaveFile],...
                    'Delimiter',',','QuoteStrings',true);
    fprintf('>> Dataset saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved dataset.\n')
end

%% END ~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~