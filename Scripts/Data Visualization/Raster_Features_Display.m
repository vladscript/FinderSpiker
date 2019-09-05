%% Raster Features Display#################################################
% Script that read and plot boxplots of Raster Features: [24]
% between or among Conditions for Several Experiments.
%                                                                      Figs
% 'RateNeurons','ActivityTime','ActiveRatioCAG','EffectiveActivity',... [1]
% 'ISImean','ISImode','ISIvar','ISIskew','ISIkurt',...                  [2]
% 'Lengthmean','Lengthmode','Lengthvar','Lengthskew','Lengthkurt',...   [2]
% 'CAGmean','CAGmode','CAGvar','CAGskew','CAGkurt',...                  [3]
% 'RoAmean','RoAmode','RoAvar','RoAskew','RoAkurt'                      [3]
%% Read CSV Files
NC = inputdlg('Number of Conditions: ',...
             'Raster Features', [1 75]);
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
CurrentPathOK=[CurrentPath(1:Slshes(end)),'Raster Features'];
%% Condition LOOP
RASTER_FEATURES={};
for i=1:NC
    % Read Names
    [FileName,PathName] = uigetfile('*.csv',['Select MULTIPLE CSV files for: ',Names_Conditions{i}],...
    'MultiSelect', 'on',CurrentPathOK);
    % Loop to Features from read csv
    if iscell(FileName)
        [~,NR]=size(FileName);
    else
        NR=1;
        % FileName=FileName
        FileName=mat2cell(FileName,1);
    end
    Features=[];
    Raster_Names=cell(NR,1);
    for r=1:NR
        LowLine=find(FileName{r}=='_');
        Raster_Names{r}=FileName{r}(1:LowLine(1)-1);
        rowFeatures=csvread([PathName,FileName{r}],1,0);
        Features=[Features;rowFeatures];
    end
    RASTER_FEATURES{i}=Features;
    RASTER_NAMES{i}=Raster_Names;
    CurrentPathOK=PathName;
end

%% Plot Data############################################################
%% 'Active Neurons','Active Time','MeanActivity','EffectiveActivity'
FeaturesA=figure;
FeaturesA.Name='Activity Indexes';
h1=subplot(2,2,1); % Active Neurons
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,1,h1)
title(h1,'Active Neurons')
h2=subplot(2,2,2); % Duration
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,2,h2)
title(h2,'Active Time Fraction')
h3=subplot(2,2,3); % Mean Activity
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,3,h3)
title(h3,'Active CAG Area Ratio')
h4=subplot(2,2,4); % Effective Activity
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,4,h4)
title(h4,'Effective Activity')
%% ITI & LT Statistics Features
% 'ISImean','ISImode','ISImedian','ISIvar','ISIskew','ISIkurt',...
% 'Lengthmean','Lengthmode','Lengthmedian','Lengthvar','Lengthskew','Lengthkurt'
FeaturesB=figure;
FeaturesB.Name='IEI & ED statistics';
g1=subplot(2,6,1); % ISI mean
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,5,g1)
title(g1,'IEI mean')
h1.YLim=[0,1];grid(h1,'on');

g2=subplot(2,6,2); % ISImode
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,6,g2)
title(g2,'IEI mode')

g3=subplot(2,6,3); % ISImedian
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,7,g3)
title(g3,'IEI median')

g4=subplot(2,6,4); % ISIvar
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,8,g4)
title(g4,'IEI variance')

g5=subplot(2,6,5); % ISIskew
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,9,g5)
title(g5,'IEI skewness')

g6=subplot(2,6,6); % ISIkurt
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,10,g6)
title(g6,'IEI kurtosis')

g7=subplot(2,6,7); % Length mean
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,11,g7)
title(g7,'ED mean')

g8=subplot(2,6,8); % Length mode
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,12,g8)
title(g8,'ED mode')

g9=subplot(2,6,9); % Length median
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,13,g9)
title(g9,'ED mode')

g10=subplot(2,6,10); % Length var
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,14,g10)
title(g10,'ED variance')

g11=subplot(2,6,11); % Length skew
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,15,g11)
title(g11,'ED skewness')

g12=subplot(2,6,12); % Length kurt
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,16,g12)
title(g12,'ED kurtosis')

%% CAG & RoA Features
FeaturesC=figure;
FeaturesC.Name='Raster Features: CAG & RoA Statitics';
j1=subplot(2,6,1); % CAG mean
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,17,j1)
title(j1,'CAG mean')

j2=subplot(2,6,2); % CAG mode
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,18,j2)
title(j2,'CAG mode')

j3=subplot(2,6,3); % CAG median
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,19,j3)
title(j3,'CAG median')

j4=subplot(2,6,4); % CAG var
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,20,j4)
title(j4,'CAG variance')

j5=subplot(2,6,5); % CAG skew
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,21,j5)
title(j5,'CAG skewness')

j6=subplot(2,6,6); % CAG kurt
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,22,j6)
title(j6,'CAG kurtosis')

j7=subplot(2,6,7); % RoA mean
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,23,j7)
title(j7,'RoA mean')

j8=subplot(2,6,8); % RoA mode
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,24,j8)
title(j8,'RoA mode')

j9=subplot(2,6,9); % RoA median
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,25,j9)
title(j9,'RoA median')

j10=subplot(2,6,10); % RoA var
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,26,j10)
title(j10,'RoA variance')

j11=subplot(2,6,11); % RoA skew
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,27,j11)
title(j11,'RoA skewness')

j12=subplot(2,6,12); % RoA kurt
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,28,j12)
title(j12,'RoA kurtosis')

%% RoT Features (Rate of Calcium Transients)
FeaturesD=figure;
FeaturesD.Name='Raster Features: Rate of Calcium Transients Statitics';
k1=subplot(2,3,1); % RoT mean
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,29,k1)
title(k1,'RoT mean')

k2=subplot(2,3,2); % RoT mode
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,30,k2)
title(k2,'RoT mode')

k3=subplot(2,3,3); % RoT median
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,31,k3)
title(k3,'RoT median')

k4=subplot(2,3,4); % RoT var
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,32,k4)
title(k4,'RoT variance')

k5=subplot(2,3,5); % RoT skew
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,33,k5)
title(k5,'RoT skewness')

k6=subplot(2,3,6); % RoT kurt
plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,34,k6)
title(k6,'RoT kurtosis')

%% Make and Save Table
okbutton = questdlg('Make CSV Table?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Table_Raster_Features_',TS,'.csv'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    disp('>>Making CSV table...')
    TableFeatures=maketableraster(RASTER_NAMES,RASTER_FEATURES,Names_Conditions);
    writetable(TableFeatures,[PathSave,SaveFile],...
                    'Delimiter',',','QuoteStrings',true);
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')