%% Set Filter to Ban Experiments from Classification Stage
% Due to undesired very low activity or atifacts
% Input:
%   Features Thresholds
% Output:
%   List of Experiments IDs to ban
%% Setup Criteria
% For Initial Condition Usually
ConditionNameCheck='Dyskinesia';
% Raster Activity Minimal Constraint:
RateNeurons_Threshold=0;
EffectiveActivity_Threshold=0.25;
% Ensembles General Minimal Constraint:
RateCycles=0;
%% Load Datasetes
% Raster Activity Minimal Constraint:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % Finder Spiker Main Folder
% Load File 
[FileName,PathName] = uigetfile('*.csv','Raster Activity Dataset file',...
    'MultiSelect', 'off',CurrentPathOK);
CurrentPathOK=PathName;
Xraster=readtable([PathName,FileName]);
% Labels
Yraster=categorical(table2array( Xraster(:,1)) );    % Categorical ARRAYS
EXPIDraster=table2array( Xraster(:,2));             % EXPIDs Cell of Strings
Xraster=table2array( Xraster(:,3:end) );
% Ensembles General Minimal Constraint:
% Load File 
[FileName,PathName] = uigetfile('*.csv','Ensembles General Dataset file',...
    'MultiSelect', 'off',CurrentPathOK);

Xensemble=readtable([PathName,FileName]);
% Labels
Yensemble=categorical(table2array( Xensemble(:,1)) );    % Categorical ARRAYS
EXPIDensemble=table2array( Xensemble(:,2));             % EXPIDs Cell of Strings
Xensemble=table2array( Xensemble(:,3:end) );