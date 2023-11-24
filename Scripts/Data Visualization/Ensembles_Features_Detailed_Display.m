%% Neural Ensemble Detailed Features Display
% Script that read and boxplot by EACH Ensembles Features
% between or among Conditions of one Several Experiments
%                                                   Figs
% '%Neurons','Time Active',...                      [1]
% 'Domminance','ErrorClass',...                     [1]
% 'ensCAGauc','ensCAGmean','ensCAGvar',...          [2]
% 'ensCAGskew','ensCAGkurt',...                     [2]
% 'IEImean','IEIvar','IEIskew','IEIkurt',...        [3]
% 'EDmean','EDvar','EDskew','EDkurt',...            [4]
%% Read CSV Files
NC = inputdlg('Number of Conditions: ',...
             'Detailed Neural Ensembles Features', [1 75]);
NC = str2double(NC{:});    
% Setup Conditions
Conditions_String='Condition_';
n_conditions=[1:NC]';
Conditions=[repmat(Conditions_String,NC,1),num2str(n_conditions)];
Cond_Names=cell(NC,1);
% Names=cell(NC,1);
for i=1:NC
    Cond_Names{i}=Conditions(i,:);
    Names_default{i}=['...'];
end
% 2nd Input Dialogue
name='Names';
numlines=[1 75];
Names_Conditions=inputdlg(Cond_Names,name,numlines,Names_default);
% Directory (default)
CurrentPath=pwd;
Slshes=find(CurrentPath=='\');
% [CurrentPath(1:Slshes(end)),'Raster Features']
CurrentPathOK=[CurrentPath(1:Slshes(end)),'Ensemble Features'];
%% Condition LOOP
ENSEMBLE_FEATURES={};
for i=1:NC
    % Read Names
    [FileName,PathName] = uigetfile('*.csv',['CSV files for: ',Names_Conditions{i}],...
    'MultiSelect', 'on',CurrentPathOK);
    % Loop to Features from read csv
    if iscell(FileName)
        [~,NR]=size(FileName);
    else
        NR=1;
        % FileName=FileName
        FileName=mat2cell(FileName,1);
    end
    FeaturesSingle=[];
    Ensemble_Names=cell(NR,1);
    for r=1:NR
        fprintf('>>Reading %s\n',FileName{r})
        LowLine=find(FileName{r}=='_');
        Ensemble_Names{r}=FileName{r}(1:LowLine(1)-1);
        TableChek=readtable([PathName,FileName{r}]);
        [Nr,Nc]=size(TableChek);
        if Nr>0
            rowFeatures=csvread([PathName,FileName{r}],1,0);
            FeaturesSingle=[FeaturesSingle;rowFeatures]; 
            NensemblesExp(r)=size(rowFeatures,1);
        else
            % rowFeatures=zeros(1,Nc);
            % FeaturesSingle=[FeaturesSingle;rowFeatures]; 
            NensemblesExp(r)=0;
        end
    end
    NensemblesCondition{i}=NensemblesExp;
    ENSEMBLE_FEATURES{i}=FeaturesSingle;
    ENSEMBLE_NAMES{i}=Ensemble_Names;
    CurrentPathOK=PathName;
end
%% Plot Data
%% Active Ensmeble Features
FeaturesA=figure;
FeaturesA.Name='Neural Ensemble Detailed Features';
h1=subplot(2,2,1);  % % of Neurons Used
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,1,h1,0)
title(h1,'Neurons Portion')
h2=subplot(2,2,2);  % Time Active
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,2,h2,0)
title(h2,'Time Active')
h3=subplot(2,2,3);  % Domminance
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,3,h3,0)
title(h3,'Domminance')
h4=subplot(2,2,4);  % Ensemble Rate
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,4,h4,0)
title(h4,'Ensemble Rate')

%% Ensmeble CAG Features
FeaturesB=figure;
FeaturesB.Name='CAG Ensemble Features';
g1=subplot(2,4,1);  % Ensemble AUC CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,5,g1,0)
title(g1,'Ensemble AUC CAG')
g2=subplot(2,4,2);  % Ensemble mean CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,6,g2,0)
title(g2,'Ensemble mean CAG')
g3=subplot(2,4,3);  % Ensemble mean CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,7,g3,0)
title(g3,'Ensemble mode CAG')
g4=subplot(2,4,4);  % Ensemble mean CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,8,g4,0)
title(g4,'Ensemble median CAG')
g5=subplot(2,4,5);  % Ensemble var CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,9,g5,0)
title(g5,'Ensemble var CAG')
g6=subplot(2,4,6);  % Ensemble skew CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,10,g6,0)
title(g6,'Ensemble skew CAG')
g7=subplot(2,4,7);  % Ensemble kurt CAG
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,11,g7,0)
title(g7,'Ensemble kurt CAG')

%% Inter Ensmeble Interval Satatistics
FeaturesC=figure;
FeaturesC.Name='Inter Ensemble Interval Statistics';
j1=subplot(2,3,1);  % mean IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,12,j1,0)
title(j1,'mean IEI')
j2=subplot(2,3,2);  % mode IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,13,j2,0)
title(j2,'mode IEI')
j3=subplot(2,3,3);  % median IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,14,j3,0)
title(j3,'median IEI')
j4=subplot(2,3,4);  % var IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,15,j4,0)
title(j4,'var IEI')
j5=subplot(2,3,5);  % skew IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,16,j5,0)
title(j5,'skew IEI')
j6=subplot(2,3,6);  % kurt IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,17,j6,0)
title(j6,'kurt IEI')

%% Inter Ensmeble Interval Satatistics
FeaturesD=figure;
FeaturesD.Name='Ensemble Duration Statistics';
k1=subplot(2,3,1);  % mean ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,18,k1,0)
title(k1,'mean ED')
k2=subplot(2,3,2);  % mode ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,19,k2,0)
title(k2,'mode ED')
k3=subplot(2,3,3);  % median ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,20,k3,0)
title(k3,'median ED')
k4=subplot(2,3,4);  % var IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,21,k4,0)
title(k4,'var ED')
k5=subplot(2,3,5);  % skew IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,22,k5,0)
title(k5,'skew ED')
k6=subplot(2,3,6);  % kurt IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,23,k6,0)
title(k6,'kurt ED')
%% Make and Save Table
okbutton = questdlg('Make CSV Table?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Detailed_Ensemble_Features_',TS,'.csv'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK,'Set Directory to save Detailed Ensemble Features');
    disp('>>Making CSV table...')
    TableFeatures=maketableraster_ensembles(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,NensemblesCondition,Names_Conditions);
    writetable(TableFeatures,[PathSave,SaveFile],...
                    'Delimiter',',','QuoteStrings',true);
    % fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
    fprintf('<a href="matlab:dos(''explorer.exe /e, %s, &'')">> Data saved Here</a>\n',PathSave);
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')
