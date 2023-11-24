%% Neural Ensemble General Features Display
% Script that read and boxplot GENERAL Ensembles Features
% between or among Conditions of one Several Experiments
%                                                   Figs
% 'N ensembles','Threshold','Dunn',...              [1]
% 'MaxIV','ErrorClass','CoreRatio',...              [1]
% 'CAGauc','Transitions Rate','Cycles Rate',...     [2]
% 'Simple Cycles','Closed Cycles','Open Cycles',    [2]
%  Statistics of Ensemble Duration                   [3]
%  Statistics of Ensemble Intervals                  [3]
% 'SynWeight_mean','SynWeight_var',...              [4]
% 'SynWeight_skew','SynWeight_kurt'                 [4]
%% Read CSV Files
NC = inputdlg('Number of Conditions: ',...
             'General Ensemble Features', [1 75]);
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
        LowLine=find(FileName{r}=='_');
        Ensemble_Names{r}=FileName{r}(1:LowLine(1)-1);
        rowFeatures=csvread([PathName,FileName{r}],1,0);
        FeaturesSingle=[FeaturesSingle;rowFeatures];
        
    end
    ENSEMBLE_FEATURES{i}=FeaturesSingle;
    ENSEMBLE_NAMES{i}=Ensemble_Names;
    CurrentPathOK=PathName;
end
%% Plot Data
%% Clustering related Features
FeaturesA=figure;
FeaturesA.Name='Neural Ensemble Clustering Features';
h1=subplot(2,3,1);  % N-Ensembles
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,1,h1)
title(h1,'N-ensembles')
h2=subplot(2,3,2);  % CAG threshold
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,2,h2)
title(h2,'CAG threshold')
h3=subplot(2,3,3);  % Dunns Index approx
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,3,h3)
title(h3,'Dunns')
h4=subplot(2,3,4);  % Maximum Distance Intra Vectors
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,4,h4)
title(h4,'MaxDIV')
h5=subplot(2,3,5);  % Classification Error
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,5,h5)
title(h5,'Validation Error')
h6=subplot(2,3,6);  % Core Ratio 
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,12,h6)
title(h6,'Core Size')

%% Clustering related Features
FeaturesB=figure;
FeaturesB.Name='Alternance & Reverberance';
g1=subplot(2,3,1);  % CAG autocorrelation
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,11,g1)
title(g1,'CAG AUC')
g2=subplot(2,3,2);  % Transitions Rate
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,6,g2)
title(g2,'Transitions/min')
g3=subplot(2,3,3);  % Cycles Rate
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,7,g3)
title(g3,'Cycles/min')
g4=subplot(2,3,4);  % Proportion Simple Cycles
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,8,g4)
title(g4,'Simple Cycles')
g5=subplot(2,3,5);  % Proportion Closed Cycles
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,9,g5)
title(g5,'Closed Cycles')
g6=subplot(2,3,6);  % Proportion Open Cycles
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,10,g6)
title(g6,'Open Cycles')

%% Time Ensemble related Features
FeaturesE=figure;
FeaturesE.Name='Neural Ensembles Time Features';
m1=subplot(2,6,1);  % mean ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,16,m1)
title(m1,'Mean Ensemble Duration')
m2=subplot(2,6,2);  % mode ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,17,m2)
title(m2,'Mode Ensemble Duration')
m3=subplot(2,6,3);  % median ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,18,m3)
title(m3,'Median Ensemble Duration')
m4=subplot(2,6,4);  % Variance ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,19,m4)
title(m4,'Mode Ensemble Duration')
m5=subplot(2,6,5);  % mode ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,20,m5)
title(m5,'Mode Ensemble Duration')
m6=subplot(2,6,6);  % mode ED
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,21,m6)
title(m6,'Mode Ensemble Inter Ensemble Interval')
m7=subplot(2,6,7);  % mean IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,22,m7)
title(m7,'Mean Ensemble Inter Ensemble Interval')
m8=subplot(2,6,8);  % mode IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,23,m8)
title(m8,'Mode Ensemble Inter Ensemble Interval')
m9=subplot(2,6,9);  % median IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,24,m9)
title(m9,'Median Ensemble Inter Ensemble Interval')
m10=subplot(2,6,10);  % Variance IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,25,m10)
title(m10,'Mode Ensemble Inter Ensemble Interval')
m11=subplot(2,6,11);  % mode IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,26,m11)
title(m11,'Mode Ensemble Inter Ensemble Interval')
m12=subplot(2,6,12);  % mode IEI
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,27,m12)
title(m12,'Mode Ensemble Duration')


%% Network related Features
FeaturesC=figure;
FeaturesC.Name='Functional Connectivity Weight Statistics';
k1=subplot(2,3,1);  % Mean Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,28,k1)
title(k1,'mean Links Weight')
k2=subplot(2,3,2);  % Mode Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,29,k2)
title(k2,'mode Links Weight')
k3=subplot(2,3,3);  % Median Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,30,k3)
title(k3,'medianLinks Weight')
k4=subplot(2,3,4);  % Variance Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,31,k4)
title(k4,'var Links Weight')
k5=subplot(2,3,5);  % Skewness Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,32,k5)
title(k5,'skew Links Weight')
k6=subplot(2,3,6);  % Kurtosis Links Weight
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,33,k6)
title(k6,'kurt Links Weight')

%% Network related Features
FeaturesD=figure;
FeaturesD.Name='Alternative & Reactivation Ensemble Index';
l1=subplot(2,1,1);  % Alternance Index
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,34,l1)
title(l1,'Alternative Ensemble Index')
l2=subplot(2,1,2);  % Reactivation Index
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,35,l2)
title(l2,'Reactivation Ensemble Index')
%% Make and Save Table
okbutton = questdlg('Make CSV Table?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Ensemble_Features_',TS,'.csv'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK,'Set Directory to save Ensemble Features');
    disp('>>Making CSV table...')
    TableFeatures=maketableraster(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions);
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

% %% Pie Charts of Kind of Cycles
% % 'ISImean','ISImode','ISIvar','ISIskew','ISIkurt',...
% % 'Lengthmean','Lengthmode','Lengthvar','Lengthskew','Lengthkurt'
% FeaturesB=figure;
% FeaturesB.Name='Type of Cycles or Paths';
% % [Simple|Open|Closed]
% ColorPie=[44,188,187;110,59,128;165,71,11]/255;
% for c=1:NC
%     ax=subplot(1,NC,c);
%     F=ENSEMBLE_FEATURES{c}(:,6:8);
%     [Nr,~]=size(F);
%     for r=1:Nr
%         F(r,isnan(F(r,:)))=0;
%         MINF=min(F(r,F(r,:)>0));
%         if isempty(MINF)
%             MINF=0;
%         end
%         if MINF==0
%             N=0;
%         else
%             N=1/MINF;
%         end
%         auxF(r,:)=N*F(r,:);
%     end
%     Pf=pie(sum(auxF));
%     Pf(1).FaceColor=ColorPie(1,:);
%     Pf(3).FaceColor=ColorPie(2,:);
%     Pf(5).FaceColor=ColorPie(3,:);
%     title(ax,Names_Conditions{c})
%     l=legend('Simple','Open','Closed');
%     l.Location='southoutside';
%     l.Box='off';
% end

% %% Dominant Ensemble Activity
% FeaturesC=figure;
% FeaturesC.Name='Dominant Ensemble Features';
% g1=subplot(1,2,1); % ISImode
% plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,9,g1)
% title(g1,'Rate')
% 
% g2=subplot(1,2,2); % ISImode
% plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,10,g2)
% title(g2,'Dominance Index')

