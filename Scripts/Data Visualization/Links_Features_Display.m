%% Links Stats Table Maker ################################################
% Script that makes the table to merge data from Link Statistics
%                                                                      Figs
%% Read CSV Files
NC = inputdlg('Number of Conditions: ',...
             'Links Stats', [1 75]);
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
CurrentPathOK=[CurrentPath(1:Slshes(end)),'NetWorks-CSV\'];
%% Condition LOOP
LINKS_FEATURES={};
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
    Features=[];
    Raster_Names=cell(NR,1);
    for r=1:NR
        LowLine=find(FileName{r}=='_');
        Raster_Names{r}=FileName{r}(1:LowLine(1)-1);
        rowFeatures=csvread([PathName,FileName{r}],1,0);
        if numel(rowFeatures)>7
            disp('>>Partial Network: Only Positive Cells')
            rowFeatures=rowFeatures([1,3:8]);
        else
            disp('>>Total Network')
        end
        Features=[Features;rowFeatures];
    end
    LINKS_FEATURES{i}=Features;
    RASTER_NAMES{i}=Raster_Names;
    CurrentPathOK=PathName;
end

%% Plot Data############################################################
%% Links Dsitribution Statistics
FeaturesA=figure;
FeaturesA.Name='Activity Indexes';
h1=subplot(2,2,1); % Active Neurons
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,1,h1)
title(h1,'CAG Autocorrelation')
h2=subplot(2,2,2); % Duration
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,2,h2)
title(h2,'Mean Links')
h3=subplot(2,2,3); % Mean Activity
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,3,h3)
title(h3,'Mode Links')
h4=subplot(2,2,4); % Effective Activity
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,4,h4)
title(h4,'Median Links')

FeaturesB=figure;
FeaturesB.Name='Activity Indexes';
g1=subplot(2,2,1); % Active Neurons
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,5,g1)
title(g1,'Variance Links')
g2=subplot(2,2,2); % Duration
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,6,g2)
title(g2,'Skewness Links')
g3=subplot(2,2,3); % Mean Activity
plot_box(RASTER_NAMES,LINKS_FEATURES,Names_Conditions,7,g3)
title(g3,'Kurtosis Links')


%% Make and Save Table
okbutton = questdlg('Make CSV Table?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Links_Stats_Features_',TS,'.csv'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    disp('>>Making CSV table...')
    TableFeatures=maketableraster(RASTER_NAMES,LINKS_FEATURES,Names_Conditions);
    B = (TableFeatures.Var3)
    writetable(TableFeatures,[PathSave,SaveFile],...
                    'Delimiter',',','QuoteStrings',true);
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')